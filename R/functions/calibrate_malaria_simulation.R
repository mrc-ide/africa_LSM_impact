# parameter_set <- setup_initial
# prevalence_target <- "0.1"
# prevalence_time <- "100"

calibrate_malaria_simulation <- function(
 parameter_set, #Customised by malaria_simulation_setup
 prevalence_target, #The prevalence target to fit to. If multiple values enter as a single object separated by ; (e.g. "0.5;0.25;0.1")
 prevalence_time, #The times (in days) at which to fit the prevalence to. If multiple values enter as a single object separated by ; (e.g. "100;400;900")
 tolerance = 0.02 #The prevalence difference above or below at which it will accept the calibration
){
  
  set.seed(1)
  
  calibrate_out <- calibrate(parameters = parameter_set,
                             target = as.numeric(if(grepl(";", prevalence_target)) unlist(strsplit(prevalence_target, ";")) else prevalence_target),
                             target_tt = as.numeric(if(grepl(";", prevalence_time)) unlist(strsplit(prevalence_time, ";")) else prevalence_time),
                             summary_function = summary_pfpr_0_5,   #This returns the prevalence in 0-5 year olds
                             tolerance = tolerance,                      
                             interval = c(0.001, ifelse(prevalence_target >= 0.4, 
                                                        25000, 750))) ##upper bound needs to be high enough so negative differences are not returned in uniroot
  
  data.frame(parameter_set$input_data,
             prevalence_target = prevalence_target,
             tolerance = tolerance,
             prevalence_time = prevalence_time,
             EIR_fit = calibrate_out$root,
             iterations = calibrate_out$iter
  )
  
}




