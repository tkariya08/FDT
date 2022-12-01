
#### code to put in github ##################

require(ggplot2);require(dplyr);require(reshape2);require(purrr);require(lubridate);require(tidyr)

set.seed(7389)

sites<-c("MLO","ASK","MID","WIS","AZR","NWR","SHM","BRW","ZEP","ALT") #Station Names


get_source_file <- function(station, tracegas) {
  
  #Function: calling Flask data files from stored directory
  #Input:    station,tracegas
  #Returns:  Flask data file path
  
  return(paste0("Flask_data_", 
                tracegas,"_", tolower(station), ".csv"))
}

prepData  <- function(station, tracegas) {
  
  #Function: preparing data for analysis,
  #          selecting good flask measurements,
  #          obtaining mean of the good flask pairs 
  #Input:    station,tracegas
  #Returns:  data frame with Date and corresponding mean flask measurement (CO2)
  
  data <- readr::read_csv(get_source_file(station, tracegas), skip=69,comment="#") #skip added by TK
  data <- data %>% dplyr::filter(analysis_flag=="...") %>% ##dplyr added by tk
    mutate(Date = ymd(paste(sample_year, sample_month, sample_day))) %>%
    group_by(Date) %>% dplyr::summarise(CO2 = round(mean(analysis_value),2)) %>% ungroup() %>%
    mutate(numericDate=as.numeric(Date))
  data
  
}


firstDataFit <- function(data) {
  
  #Function: Fitting quadratic polynomial + 4 harmonic sinusoidal wave function 
  #Input:    data from prepData above
  #Returns:  Residuals on which loess fitting is performed
  
  wl <- 365.25/2./pi #length of year in days

  trend<- lm(CO2 ~ poly(numericDate, 2) + 
               sin(numericDate/wl) + cos(numericDate/wl) +
               sin(numericDate/wl*2) + cos(numericDate/wl*2)+
               sin(numericDate/wl*3) + cos(numericDate/wl*3)+
               sin(numericDate/wl*4) + cos(numericDate/wl*4),
             data)
  CO2_trends<-predict(trend)
  
  data  <- data %>% mutate(CO2_trend=CO2_trends, CO2residual=CO2 - CO2_trend) #Residuals R(t) 

  data
}

optiLoess <- function(x, y, spans, degrees) {
  
  #Function: for selecting optimal parameters for loess fitting
  #Input:    numericDate, residual, range of span, degree
  #Returns:  best span
  
  require(caret)
  if (is.null(names(x))) x <- cbind(x=x)
  fitControl <- trainControl(method = "repeatedcv",
                             number = 10, repeats = 5) #trainControl decides the resampling method
  searchgrid <- expand.grid(span = spans, degree = degrees) 
  train_loess <- train(x, y,
                       method = "gamLoess",
                       tuneGrid=searchgrid,
                       na.action = na.omit, trControl=fitControl)
  
}


findbestLoessParams <- function(data){
  
  #Function: finding the best loess parameters for the  residuals
  #Input:    data with residuals, optiLoess function
  #Returns:  best span for residuals
  
  ok  <- is.finite(data$CO2residual)
  npoints <- length(ok)
  samplingInterval <- median(diff(data$numericDate[which(is.finite(data$CO2))]))
  minSpan <- 3/npoints
  maxSpan <- 360./samplingInterval/npoints
  bestLoess  <- optiLoess(data$numericDate, data$CO2residual, degrees = c(1), 
                          spans = seq(minSpan, maxSpan, length.out = 5)) 
  print(bestLoess)
  print(ggplot(bestLoess) + ggtitle(npoints))
  bestLoess
}

fitLoess <- function(data, loessParams) {
  
  #Function: Loess fitting on residuals  
  #Input:    data with residuals, best span from findbestLoessParams
  #Returns:  smoothed time series, residuals for bootstrapping 
  
  regularTS <- tibble(Date=seq(min(data$Date), max(data$Date), by="day"))
  data <- left_join(regularTS, data, by="Date") %>%
    mutate(numericDate=as.numeric(Date))
  
  okish  <- which(is.finite(data$CO2residual))
  
# Obtaining best fit (smoothed residuals)
  bestFit <- gam.lo(data$numericDate[okish], data$CO2residual[okish], 
                    span=loessParams$bestTune$span,
                    degree = loessParams$bestTune$degree, xeval = data$numericDate )

# refit the model here, so it is surely additive 
  wl  <- 365.25 / 2. / pi
  model<- lm(CO2 ~ poly(numericDate,2, raw=T) + 
               sin(numericDate/wl) + cos(numericDate/wl) +
               sin(numericDate/wl*2) + cos(numericDate/wl*2)+
               sin(numericDate/wl*3) + cos(numericDate/wl*3)+
               sin(numericDate/wl*4) + cos(numericDate/wl*4),
             data)

    coefs <- coef(model)
  data  <- data %>% mutate(
    CO2_trend = predict(model, .),
    loessModelledAnomaly = as.vector(bestFit),
    loessFit = CO2_trend + loessModelledAnomaly,
    polyTrend = coefs[1] + coefs[2] * numericDate + coefs[3]* numericDate^2 ,
    detrended = loessFit - polyTrend, 
    residual = CO2 - loessFit) ## residual used for bootstrapping
  data
  
}


boot<-function(j, data, bestLoess, bootResToModel=TRUE){
  
  # Function: resampling residuals and generating bootstrap samples.
  # Input:    j number times to be iterated.
  # Returns:  j residual bootstrap samples.
  
  cat("\r",j)
  
  if (bootResToModel) variable <-  sym("loessModelledAnomaly") 
  else variable  <- sym("CO2residual")
  
  data  <- data %>% mutate(bootResidual = sample(residual, replace=TRUE),
                           CO2residualOrig = CO2residual,
                           CO2residual = !!variable + bootResidual, ## bootResidual added to  the first prediction 
                           loessModelAnomalyOrig=loessModelledAnomaly)  
  data  <- fitLoess(data, bestLoess) %>% mutate(ensemble=j) 
  data
}

FDT_CUP <- function(data) {
  
  # Function: FDT determined CUP.
  # Input:    fitted bootstrap samples
  # Returns:  data frame with FDT CUP and zero-crossing CUP

  ### Filter years with 365/366 days (avoiding end years)
  data$Date<-as.Date(data$Date)
  data <- data  %>% group_by(year=year(Date)) %>% mutate(nPerYear=n()) %>%
    ungroup() %>% filter(nPerYear>=365)
  goodyrs <- unique(year(data$Date))
  
  threshs_peak <-threshs_valley <-  c(0,0.05,0.10,0.15,0.20) ## threshold values, can add additional values
  names(threshs_peak) <- paste0("Thp_",as.character(threshs_peak))
  names(threshs_valley) <- paste0("Thv_",as.character(threshs_valley))
  
  ##### peak, valley and CUP corresponding to different threshold ######
  slopesEtc <- data %>% group_by(year=year(Date)) %>%
    dplyr::summarise(max_slope = c(NA,detrended) %>% diff %>% min(na.rm=T),
                     max_point = c(NA,detrended) %>% diff %>% which.min(),
                     maaxx = list(c(NA,detrended) %>% diff),
                     thpeak = map_dfc(threshs_peak, ~max(which(unlist(maaxx)[1:max_point] > max_slope*.x))),
                     thvalley=map_dfc(threshs_valley, ~min(which(unlist(maaxx)[max_point:length(unlist(maaxx))] > max_slope*.x))+(max_point-1)),
                     thp_day=Date[thpeak$Thp_0.15],
                     thv_day=Date[thvalley$Thv_0],
                     thmin_day=Date[which.min(detrended)],
                     thmax_day=Date[which.max(detrended[1:which.min(detrended)])]
    )%>% unpack(c(thpeak,thvalley))%>%
    mutate(duration=Thv_0-Thp_0.15, ##CUP from selected threshold
           max_min_duration=as.Date(thmin_day)-as.Date(thmax_day)) %>% select(-c("maaxx")) 
  
  ##### upward, downward and CUP zero-crossing #####
  sub_set<-c(slopesEtc$thmax_day,data$Date[length(data$Date)])
  
  zero<-map_dfr(1:(length(sub_set)-1),
                ~data %>% filter(Date %in% sub_set[.x]:sub_set[.x+1]) %>%
                  summarise(year=year(Date[1]),
                            zero_dwn=Date[which(diff(sign(detrended))<0)]%>%
                              .[.<Date[which.min(detrended)]]%>% max,
                            zero_up=Date[which(diff(sign(detrended))>0)]%>%
                              .[.>Date[which.min(detrended)]]%>% min )
  ) %>% mutate(zero_duration=zero_up-zero_dwn)
  
  slopesEtc<-merge(slopesEtc,zero,by="year")
}

runSite <- function(site) {
  
  # Function: combine above functions to generate bootstrap sample,
  #           calculates FDT and zero-crossing CUP of bootstrap samples.  
  # Input:    station name 
  # Returns:  FDT and zero-crossing CUP for bootstrap samples
  
  df <- prepData(site,tracegas = "CO2")
  df  <- firstDataFit(df)
  bestLoess <- findbestLoessParams(df)
  dfFilled <- fitLoess(df, bestLoess)
  df4boot  <- dfFilled %>% filter(is.finite(CO2))
  
  bootCUP <- map_dfr(1:500, ~boot(.x, data=df4boot, bestLoess = bestLoess) %>%
                       FDT_CUP() %>% mutate(ensemble=.x)) 
}


cup_metrics<-runSite(sites[1]) ## stores output of runSite in cup_metrics 




