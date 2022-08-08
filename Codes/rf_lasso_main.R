# - Analysis of random forest with glmnet split
# - Use glmnet with degree 3
# - Use rf_lasso_utils.R
# - Including lag lead (temporal) features
# - WS as input


#libraries for operations
writeLines("Loading libraries start")
suppressMessages({
  options(stringsAsFactors=FALSE)
  options(dplyr.width=Inf)
  options(scipen=999)
  Sys.setenv(TZ='UTC')
  options(digits.secs=3)
  options(pillar.subtle=FALSE)
  httr::set_config(httr::config(http_version = 0))
  library(aws.s3)
  try(detach("package:energysupport",unload=TRUE))
  library(energysupport)
  library(efcast)
  library(prodfcast)
  library(data.table)
  library(tidyverse)
  library(lubridate)
  library(rpart)
})
writeLines("Loading libraries end")

result_path = "Results/"
result_path = paste0(result_path,gsub("-","",lubridate::today("Turkey"),fixed=T))
if(!dir.exists(result_path)){
  dir.create(result_path,recursive = T)
}
data_path = "Data/"
function_path = "Functions/"
test_start_date = lubridate::as_date("2021-01-01")
source(file.path(function_path,"rf_lasso_utils.R"))

# Get bash arguments
arg <- commandArgs(trailingOnly = TRUE)
if(length(arg)==0){
  model = "glmnet"
  general_depth = 4
  mtry_size_vec = "2_2_1"
  mtry_rep= 16
  selection="random_subset"
  general_ntree = 100
  sw = 4
}else{
  model = arg[1]
  general_depth = as.numeric(arg[2])
  mtry_size_vec = arg[3]
  mtry_rep= as.numeric(arg[4])
  selection=arg[5]
  general_ntree = as.numeric(arg[6])
  sw = as.numeric(arg[7])
}

# Get Production and Weather Data
weather = read.csv(file.path(data_path,"multiple_farm_weather_80m_ws.csv"))
production = readRDS(file.path(data_path,"multiple_farm_agg.rds"))
setDT(weather)
setDT(production)

# Process and impute weather data
# Only used 80m data
weather[, `:=`(rn, NULL)]
weather = melt(weather,1:3)
weather[,variable:=ifelse(variable%like%"UGRD","UGRD","VGRD")]
weather[, `:=`(forecast_timestamp, 
               lubridate::as_datetime(as.numeric(forecast_epoch), 
                                      tz = "Turkey"))]
weather[, `:=`(forecast_date, lubridate::date(forecast_timestamp))]
weather[, `:=`(forecast_hour, lubridate::hour(forecast_timestamp))]
weather[, `:=`(c("forecast_timestamp", 
                 "forecast_epoch","level"), NULL)]
weather=rename_weather_data(weather)
weather=impute_weather_data(weather)

# Calculate WS and WDIR
casted_weather = dcast(weather,date+hour+lat+lon~variable,value.var="value")
casted_weather[,WS:=sqrt(UGRD^2+VGRD^2)]
casted_weather[,WDIR:=(180/pi) * atan2(UGRD, VGRD)]
casted_weather[WDIR<0,WDIR:=WDIR+360]
casted_weather = casted_weather[,list(date,hour,lat,lon,UGRD,VGRD,WS,WDIR)]

# Process Production
production[,date_time:=energysupport::epoch_to_date_time(epoch)]
production[,date:=lubridate::as_date(date_time)]
production[,hour:=lubridate::hour(date_time)]
production = production[,list(date,hour,production)]

# Combine Production and weather data
combined = casted_weather %>% left_join(production) %>% na.omit() %>% as.data.table() 
setorder(combined,date,hour)

# Add temporal features
temporal_length = 2
for(t in 1:temporal_length){
  combined[,paste0("WS_lag",t):=lag(WS,t),by=list(lat,lon)]
  combined[,paste0("WS_lead",t):=lead(WS,t),by=list(lat,lon)]
}


# Get only WS features
val_var = c(colnames(combined)[colnames(combined)%like%"WS"])
wide_combined = copy(dcast(combined,date+hour+production~lat+lon,value.var = val_var)) %>%
  na.omit() %>%
  select(date,hour,production,everything()) %>%
  as.data.table()
setorder(wide_combined,date,hour)

# Initial population
lats = sort(unique(combined$lat))
lons = sort(unique(combined$lon))
lat_lon_matrix = expand_grid(lats,lons) %>% mutate(value=paste0(lats,"_",lons)) %>% spread(lons,value) %>% arrange(-lats) %>% select(-lats) %>% as.matrix()
lat_lon_time_array = array(dim=c(length(lats),length(lons),2*temporal_length+1))
lat_lon_time_array[,,temporal_length+1] = paste0("WS_",lat_lon_matrix)
for(t in 1:temporal_length){
  lat_lon_time_array[,,(temporal_length+1)-t] = paste0("WS_lag",t,"_",lat_lon_matrix)
  lat_lon_time_array[,,(temporal_length+1)+t] = paste0("WS_lead",t,"_",lat_lon_matrix)
}
n_lon = length(unique(combined$lon))
n_lat = length(unique(combined$lat))
n_temp = 2*temporal_length+1


######################### Random Forest Part
# Train test split
train_data = wide_combined[date<test_start_date]
test_data = wide_combined[date>=test_start_date]

data = copy(train_data)
target = "production"
feature_array = lat_lon_time_array

# Temporal bagging weight options
weight_opt = list()
data_size = nrow(train_data)
weight_opt[[1]] = rep(1,data_size)
weight_opt[[2]] = rev((data_size - 1:data_size) / data_size)
weight_opt[[3]] = rev((data_size - 1:data_size)^1.5 / data_size^1.5) 
weight_opt[[4]] = rev((data_size - 1:data_size)^2.5 / data_size^2.5) 
weight_opt[[5]] = rev((data_size - 1:data_size)^4 / data_size^4) 
weight_opt[[6]] = rev((data_size - 1:data_size)^6 / data_size^6) 
weight_opt[[7]] = rev(1 - (1:data_size)^1.5/data_size^1.5) 
weight_opt[[8]] = rev(1 - (1:data_size)^2.5/data_size^2.5) 
weight_opt[[9]] = rev(1 - (1:data_size)^4/data_size^4) 

# Define parameters according to input
mtry_size = as.numeric(unlist(strsplit(mtry_size_vec,"_")))
mtry_random = prod(mtry_size) * mtry_rep
control = list(feature_array=feature_array,mtry_size=mtry_size,mtry_rep=mtry_random,selection=selection,mtry_random=prod(mtry_size))
all_preds = tibble()

# Random Forest Part
rf_fit = ranger::ranger(production~.,data = train_data[,-c("date","hour"),with=F],max.depth = general_depth,
                        replace = F,num.trees = general_ntree,oob.error = T,seed=1,mtry = mtry_random,
                        case.weights = weight_opt[[sw]])
all_preds = all_preds %>% bind_rows(
    test_data %>% select(date,hour,production) %>% mutate(prediction=predict(rf_fit,test_data)$predictions,type="Test")%>%
      mutate(model="RF",
             nodes = "mean_mean",
             selection = selection,
             depth=general_depth,tree=general_ntree,
             mtry_replication=mtry_rep,
             mtry_size = paste0(mtry_size,collapse="_"),
             mtry_random=mtry_random,
             weight_type = paste0("Type",sw)) 
)

model_fit = random_forest_lm_custom_tree(data=train_data[,-c("date","hour"),with=F],target="production",
                                         ntrees=general_ntree,max_depth=general_depth,
                                         mtry_control=control,sample_weights = weight_opt[[sw]],
                                         model_type=model)
if(model=="lm"){
  # prediction = custom_tree_prediction_lm(model_fit$fits,test_data)
  # prediction1 = custom_tree_prediction_lm(model_fit$fits,test_data,type="mean_median")
  prediction2 = custom_tree_prediction_lm(model_fit$fits,test_data,type="median_median")
}else if(model=="glmnet"){
  # prediction = custom_tree_prediction_glmnet(model_fit$fits,test_data)
  # prediction1 = custom_tree_prediction_glmnet(model_fit$fits,test_data,type="mean_median")
  prediction2 = custom_tree_prediction_glmnet(model_fit$fits,test_data,type="median_median")
}

# rm(model_fit)
# gc()

all_preds = all_preds %>%  bind_rows(
  test_data %>% select(date,hour,production) %>% mutate(prediction=prediction2,type="Test")%>%
    mutate(model=paste0("RF_",model),
           nodes = "median_median",
           selection = selection,
           depth=general_depth,tree=general_ntree,
           mtry_replication=mtry_rep,
           mtry_size = paste0(mtry_size,collapse="_"),
           mtry_random=mtry_random,
           weight_type = paste0("Type",sw)) 
)

res = all_preds %>% group_by(model,nodes,selection,depth,tree,mtry_replication,mtry_size,mtry_random,weight_type) %>%
  summarise(wmape = sum(abs(prediction-production))/sum(production),
            prod=sum(production),mindate=min(date),nday=length(unique(date))) %>% ungroup()

print(res)

post_path = paste0("multiple_farm_rf_results_spatialtemporalws_",
                   paste("M",model,"D",general_depth,"MS",mtry_size_vec,
                         "R",mtry_rep,"S",selection,"T",sw,"NT",general_ntree,sep = "_"))
openxlsx::write.xlsx(res,file.path(result_path,paste0(post_path,".xlsx")),overwrite = T)
saveRDS(all_preds,file.path(result_path,paste0(post_path,".rds")))
# saveRDS(model_fit,file.path(result_path,paste0(post_path,"_model_fits.rds")))
