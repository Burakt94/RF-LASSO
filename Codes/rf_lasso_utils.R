library(partykit)
require(glmnetUtils)
require(doMC)
library(parallel)
# registerDoMC(cores = 2)

get_features <- function(mtry_control){
  if(mtry_control$selection=="random"){
    features = list()
    for(i in 1:mtry_control$mtry_rep){
      random_temp = sample(mtry_control$feature_array,mtry_control$mtry_random)
      features[[i]] = random_temp
    }
  }else if(mtry_control$selection=="conv"){
    features = list()
    for(i in 1:mtry_control$mtry_rep){
      conv_solution = array(0,dim=dim(mtry_control$feature_array))
      dim_x = dim(mtry_control$feature_array)[1]
      dim_y = dim(mtry_control$feature_array)[2]
      dim_z = dim(mtry_control$feature_array)[3]
      x_int = sort(sample(1:(dim_x-(mtry_control$mtry_size[1]-1)),1))
      y_int = sort(sample(1:(dim_y-(mtry_control$mtry_size[2]-1)),1))
      z_int = sort(sample(1:(dim_z-(mtry_control$mtry_size[3]-1)),1))
      conv_solution[(x_int:(x_int+mtry_control$mtry_size[1]-1)),
                    (y_int:(y_int+mtry_control$mtry_size[2]-1)),
                    (z_int:(z_int+mtry_control$mtry_size[3]-1))] = 1
      features[[i]] = mtry_control$feature_array[conv_solution==1]
    }
  }else if(mtry_control$selection=="all"){
    features = list()
    features[[1]] = sample(mtry_control$feature_array)
  }else if(mtry_control$selection=="random_subset"){
    features = list()
    temp_fa = sample(mtry_control$feature_array)
    features[[1]] = sample(temp_fa,mtry_control$mtry_random) 
  }
  
  return(features)
}

findsplit_lm <- function(response, data, indices,mtry_control=NULL,current_index = 0,model_type="lm"){

  lm_features = get_features(mtry_control)
  
  small_data = copy(data[indices!=0])
  
  lm_df = copy(small_data[,response,with=F])
  
  i <- 1
  info_list = lapply(lm_features,function(x){
    
    prefix_temp = paste0("V",current_index,"_",i,"_",model_type)
    i <<- i + 1 
    
    if(model_type=="lm"){
      lm_formula = as.formula(paste0(response,"~",paste0(as.character(x),collapse="+"),
                                     paste0("+I(",paste0(as.character(x),collapse = "^2)+I("),"^2)"),
                                     paste0("+I(",paste0(as.character(x),collapse = "^3)+I("),"^3)")))
      temp_lm = lm(lm_formula,data=small_data)
      feat_name = paste0(prefix_temp,"V1")
      lm_df[,paste0(feat_name):=temp_lm$fitted]
      list(vars = x,
           prefix=prefix_temp,
           formula=stripGlmLR(temp_lm))
      
    }else if(model_type=="glmnet"){
      lm_formula = as.formula(paste0(response,"~",paste0(as.character(x),collapse="+"),
                                     paste0("+I(",paste0(as.character(x),collapse = "^2)+I("),"^2)"),
                                     paste0("+I(",paste0(as.character(x),collapse = "^3)+I("),"^3)")))
      temp_lm = glmnetUtils::cv.glmnet(lm_formula,data=small_data,parallel=F,alpha=1,nlambda=50,nfolds=5)
      feat_name = paste0(prefix_temp,"V1")
      lm_df[,paste0(feat_name):=predict(temp_lm,small_data)]
      list(vars = x,
           prefix=prefix_temp,
           formula=stripGlmLR(temp_lm))
      
    }
    
  })
  
  formula <- as.formula(paste0(response,"~."))
  fit = rpart(formula,data=lm_df,cp=-1,maxdepth = 1,maxsurrogate=0,maxcompete=0)
  best_split = tryCatch({
    head(data.table(var=rownames(fit$splits),fit$splits),1)
  },error=function(e){
    return(NULL)
  })
  if(is.null(best_split)){
    return(NULL)
  }
  
  info_best = info_list[[(which(colnames(lm_df)==best_split$var)-1)]]
  best_value = 
    as.data.table(predict(info_best$formula,data,type="response"))
  
  colnames(best_value) = paste0(info_best$prefix,colnames(best_value))
  
  ## return split as partysplit object
  return(partysplit(varid = as.integer(current_index),
                    breaks = best_split$index,
                    info = list(best=best_split$var,best_value = best_value,
                                best_info=info_best)))
  
}


growtree <- function(id = 1L, response, data, indices, depth=0, minbucket = 1,max_depth=4,mtry_control=NULL,feature_selection="raw",stat_funcs =c("rowMeans"),
                     model_type="lm"){
  
  ## if reach max_tree stop here
  if(depth==max_depth){
    return(partynode(id = id))
  }
  
  ## for less than minbucket observations stop here
  if (sum(indices) < minbucket){
    return(partynode(id = id))
  } 
  ## find best split
  if(feature_selection=="stat"){
    sp <- findsplit_stat(response, data, indices,mtry_control=mtry_control,current_index = id,stat_funcs=stat_funcs)
    ## no split found, stop here
    if (is.null(sp)){
      return(partynode(id = id))
    } 
    ## actually split the data
    kidids <- .bincode(as.vector(unlist(info_split(sp)$best_value)), breaks = unique(c(-Inf, breaks_split(sp), Inf)), right_split(sp))
    
  }else if(feature_selection=="pca"){
    sp <- findsplit_pca(response, data, indices,mtry_control=mtry_control,current_index = id)
    ## no split found, stop here
    if (is.null(sp)){
      return(partynode(id = id))
    } 
    ## actually split the data
    kidids <- .bincode(as.vector(unlist(info_split(sp)$best_value)), breaks = unique(c(-Inf, breaks_split(sp), Inf)), right_split(sp))
  }else if(feature_selection=="lm"){
    sp <- findsplit_lm(response, data, indices,mtry_control=mtry_control,current_index = id,model_type=model_type)
    ## no split found, stop here
    if (is.null(sp)){
      return(partynode(id = id))
    } 
    ## actually split the data
    kidids <- .bincode(as.vector(unlist(info_split(sp)$best_value)), breaks = unique(c(-Inf, breaks_split(sp), Inf)), right_split(sp))
  }else if(feature_selection=="raw"){
    sp <- findsplit(response, data, indices,mtry_control=mtry_control)
    ## no split found, stop here
    if (is.null(sp)){
      return(partynode(id = id))
    } 
    ## actually split the data
    kidids <- kidids_split(sp, data = data)
  }
  ## set up all daugther nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))
  for (kidid in 1:length(kids)) {
    ## select observations for current node
    w <- indices
    w[kidids != kidid] <- 0
    ## get next node id
    if (kidid > 1) {
      myid <- max(nodeids(kids[[kidid - 1]]))
    } else {
      myid <- id
    }
    ## start recursion on this daugther node
    kids[[kidid]] <- growtree(id = as.integer(myid + 1), response, data, w,depth+1,minbucket,max_depth,mtry_control,feature_selection,stat_funcs,model_type)
  }
  ## return nodes
  return(partynode(id = as.integer(id), split = sp, kids = kids,
                   info = info_split(sp)))
}

mytree <- function(formula, data, indices = NULL,minbucket=1,max_depth=4,mtry_control=NULL,feature_selection = "raw",stat_funcs =c("rowMeans"),
                   model_type="lm"){
  
  ## name of the response variable
  response <- all.vars(formula)[1]
  ## data without missing values, response comes last
  data <- data[complete.cases(data), c(all.vars(formula)[-1], response),with=F]
  
  if (is.null(indices)){
    indices <- rep(1L, nrow(data))
  } 
  
  ## grow tree
  nodes <- growtree(id = 1L, response, data, indices,minbucket=minbucket,max_depth=max_depth,mtry_control=mtry_control,
                    feature_selection=feature_selection,stat_funcs=stat_funcs,model_type=model_type)
  
  if(feature_selection=="stat"){
    node_apply_counter = 0
    dummy_df = data.table(col=rep(NA,nrow(data)))
    tree_info = nodeapply(nodes,ids = nodeids(nodes),FUN = function(x){
      node_apply_counter <<- node_apply_counter + 1 
      if(is.null(x$info)){
        colnames(dummy_df) <- paste0("V",nodeids(nodes)[node_apply_counter])
        list(col=dummy_df,col_info = data.table(variable=NA,formula=NA))
      }else{
        list(col=x$info$best_value,col_info = data.table(variable=x$info$best,formula=x$info$best_info$formula))
      }
    })
    
    stat_df = tree_info[[1]]$col
    stat_info = tree_info[[1]]$col_info
    for(i in 2:length(tree_info)){
      stat_df = cbind(stat_df,tree_info[[i]]$col)
      stat_info = rbind(stat_info,tree_info[[i]]$col_info)
    }
    stat_df[,paste0(response):=data[[response]]]
    
    ## compute terminal node number for each observation
    fitted <- fitted_node(nodes, data = stat_df)
    ## return rich constparty object
    ret <- party(nodes, data = stat_df,
                 fitted = data.frame(
                   "(fitted)" = fitted,
                   "(response)" = stat_df[[response]],
                   check.names = FALSE),
                 terms = terms(formula))
    return(list(fit=as.constparty(ret),stat_info=stat_info,stat_df=stat_df))
  }else if(feature_selection=="pca"){
    node_apply_counter = 0
    dummy_df = data.table(col=rep(NA,nrow(data)))
    tree_info = nodeapply(nodes,ids = nodeids(nodes),FUN = function(x){
      node_apply_counter <<- node_apply_counter + 1 
      if(is.null(x$info)){
        colnames(dummy_df) <- paste0("V",nodeids(nodes)[node_apply_counter])
        list(col=dummy_df,col_info = NULL)
      }else{
        list(col=x$info$best_value,col_info = x$info)
      }
    })
    
    pca_df = tree_info[[1]]$col
    pca_info = list()
    pca_info[[1]] = tree_info[[1]]$col_info
    for(i in 2:length(tree_info)){
      pca_df = cbind(pca_df,tree_info[[i]]$col)
      if(!is.null(tree_info[[i]]$col_info)){
        pca_info[[length(pca_info)+1]] <- tree_info[[i]]$col_info
      }
    }
    pca_df[,paste0(response):=data[[response]]]
    
    ## compute terminal node number for each observation
    fitted <- fitted_node(nodes, data = pca_df)
    ## return rich constparty object
    ret <- party(nodes, data = pca_df,
                 fitted = data.frame(
                   "(fitted)" = fitted,
                   "(response)" = pca_df[[response]],
                   check.names = FALSE),
                 terms = terms(formula))
    return(list(fit=as.constparty(ret),pca_info=pca_info,pca_df=pca_df))
  }else if(feature_selection=="lm"){
    node_apply_counter = 0
    dummy_df = data.table(col=rep(NA,nrow(data)))
    tree_info = nodeapply(nodes,ids = nodeids(nodes),FUN = function(x){
      node_apply_counter <<- node_apply_counter + 1 
      if(is.null(x$info)){
        colnames(dummy_df) <- paste0("V",nodeids(nodes)[node_apply_counter])
        list(col=dummy_df,col_info = NULL)
      }else{
        list(col=x$info$best_value,col_info = x$info)
      }
    })
    
    lm_df = tree_info[[1]]$col
    lm_info = list()
    lm_info[[1]] = tree_info[[1]]$col_info
    for(i in 2:length(tree_info)){
      lm_df = cbind(lm_df,tree_info[[i]]$col)
      if(!is.null(tree_info[[i]]$col_info)){
        lm_info[[length(lm_info)+1]] <- tree_info[[i]]$col_info
      }
    }
    lm_df[,paste0(response):=data[[response]]]
    
    ## compute terminal node number for each observation
    fitted <- fitted_node(nodes, data = lm_df)
    ## return rich constparty object
    ret <- party(nodes, data = lm_df,
                 fitted = data.frame(
                   "(fitted)" = fitted,
                   "(response)" = lm_df[[response]],
                   check.names = FALSE),
                 terms = terms(formula))
    return(list(fit=as.constparty(ret),lm_info=lm_info,lm_df=lm_df))
  }else if(feature_selection=="raw"){
    ## compute terminal node number for each observation
    fitted <- fitted_node(nodes, data = data)
    ## return rich constparty object
    ret <- party(nodes, data = data,
                 fitted = data.frame("(fitted)" = fitted,
                                     "(response)" = data[[response]],
                                     check.names = FALSE),
                 terms = terms(formula))
    return(list(fit=as.constparty(ret)))
  }
  
}


random_forest_lm_custom_tree <- function(data,target,ntrees,minbucket=5,max_depth=4,bagging_ratio=0.632,mtry_control=NULL,sample_weights=NULL,bagging_type="random",
                                         model_type="lm"){
  
  if(is.null(sample_weights)){
    sample_weights = rep(1,nrow(data))
  }
  
  nrow_data = nrow(data)
  colnames(data)[colnames(data) == target] = "target_var"
  formula = as.formula(paste0("target_var~",paste0(setdiff(colnames(data),c("target_var")),collapse = "+")))
  
  fits = mclapply(1:ntrees,function(n){
    print(n)
    if(bagging_type=="random"){
      in_index = sample(1:nrow_data,nrow_data*bagging_ratio,prob = sample_weights)
    }else if(bagging_type=="bulk"){
      possibleIndex = seq(nrow_data - nrow_data*(bagging_ratio) + 1)
      firstIndex = sample(possibleIndex, 1,prob = sample_weights[possibleIndex])
      in_index = firstIndex:(firstIndex + nrow_data*(bagging_ratio) -1)
    }
    
    fit = mytree(formula, data=data[in_index,], indices = NULL,minbucket=minbucket,max_depth=max_depth,mtry_control = mtry_control,
                 feature_selection = "lm",model_type=model_type)
    fit$fit$data = head(fit$fit$data,1)
    return_fit = list(fit=fit$fit,info=fit$lm_info)
    rm(fit)
    gc()
    return(return_fit)
  },mc.cores = 8)
  
  
  return(list(fits=fits))
  
}

custom_tree_prediction_lm <- function(fit_list,data,type=NULL){
  if(is.null(type)){
    preds = unlist(lapply(fit_list, function(x){
      dummy_data = copy(data[,1])
      for(r in 1:length(x$info)){
        dummy_data[,paste0(x$info[[r]]$best_info$prefix,"V1"):=
                     predict(x$info[[r]]$best_info$formula,data,type="response")]
      }
      prediction = predict(x$fit,dummy_data)
      return(prediction)
    }))
    pred_df = data.table(name=as.numeric(names(preds)),pred=preds) %>% group_by(name) %>% summarise(prediction = mean(pred)) 
  }else if(type=="mean_median"){
    preds = unlist(lapply(fit_list, function(x){
      dummy_data = copy(data[,1])
      for(r in 1:length(x$info)){
        dummy_data[,paste0(x$info[[r]]$best_info$prefix,"V1"):=
                     predict(x$info[[r]]$best_info$formula,data,type="response")]
      }
      prediction = predict(x$fit,dummy_data)
      return(prediction)
    }))
    pred_df = data.table(name=as.numeric(names(preds)),pred=preds) %>% group_by(name) %>% summarise(prediction = median(pred)) 
  }else if(type=="median_median"){
    preds = unlist(lapply(fit_list, function(x){
      dummy_data = copy(data[,1])
      for(r in 1:length(x$info)){
        dummy_data[,paste0(x$info[[r]]$best_info$prefix,"V1"):=
                     predict(x$info[[r]]$best_info$formula,data,type="response")]
      }
      
      nodes = predict(x$fit,dummy_data,type="node")
      temp_pred_df = x$fit$fitted %>% group_by(node=`(fitted)`) %>% summarise(mean=mean(`(response)`),median=median(`(response)`))
      matched = temp_pred_df$median[match(nodes,temp_pred_df$node)]
      prediction = setNames(matched, 1:length(matched))
      return(prediction)
    }))
    pred_df = data.table(name=as.numeric(names(preds)),pred=preds) %>% group_by(name) %>% summarise(prediction = median(pred)) 
  }
  return(pred_df$prediction)
}
custom_tree_prediction_glmnet <- function(fit_list,data,type=NULL){
  if(is.null(type)){
    preds = unlist(lapply(fit_list, function(x){
      dummy_data = copy(data[,1])
      for(r in 1:length(x$info)){
        dummy_data[,paste0(x$info[[r]]$best_info$prefix,"lambda.1se"):=
                     predict(x$info[[r]]$best_info$formula,data,type="response")[,1]]
      }
      prediction = predict(x$fit,dummy_data)
      return(prediction)
    }))
    pred_df = data.table(name=as.numeric(names(preds)),pred=preds) %>% group_by(name) %>% summarise(prediction = mean(pred)) 
  }else if(type=="mean_median"){
    preds = unlist(lapply(fit_list, function(x){
      dummy_data = copy(data[,1])
      for(r in 1:length(x$info)){
        dummy_data[,paste0(x$info[[r]]$best_info$prefix,"lambda.1se"):=
                     predict(x$info[[r]]$best_info$formula,data,type="response")[,1]]
      }
      prediction = predict(x$fit,dummy_data)
      return(prediction)
    }))
    pred_df = data.table(name=as.numeric(names(preds)),pred=preds) %>% group_by(name) %>% summarise(prediction = median(pred)) 
  }else if(type=="median_median"){
    preds = unlist(lapply(fit_list, function(x){
      dummy_data = copy(data[,1])
      for(r in 1:length(x$info)){
        dummy_data[,paste0(x$info[[r]]$best_info$prefix,"lambda.1se"):=
                     predict(x$info[[r]]$best_info$formula,data,type="response")[,1]]
      }
      
      nodes = predict(x$fit,dummy_data,type="node")
      temp_pred_df = x$fit$fitted %>% group_by(node=`(fitted)`) %>% summarise(mean=mean(`(response)`),median=median(`(response)`))
      matched = temp_pred_df$median[match(nodes,temp_pred_df$node)]
      prediction = setNames(matched, 1:length(matched))
      return(prediction)
    }))
    pred_df = data.table(name=as.numeric(names(preds)),pred=preds) %>% group_by(name) %>% summarise(prediction = median(pred)) 
  }
  return(pred_df$prediction)
}

stripGlmLR = function(cm) {
  cm$y = c()
  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()  
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  
  
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  cm
}
