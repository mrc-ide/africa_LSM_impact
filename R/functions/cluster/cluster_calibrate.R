# row_MINT_combinations = run_all$row_MINT_combinations[1]
# savemod = run_all$savemod[1]

cluster_calibrate <- function(row_MINT_combinations,
                              savemod){
  
  #Front part of the file repo and then create a folder if it doesnt exist
  front_part <- "1_ModelSimulations/output/EIR_calibrate/"
  if(!dir.exists(paste0(front_part, savemod, "/"))) dir.create(paste0(front_part, savemod, "/"), recursive = T)
  
  #Read in and subset parameters for baseline 
  all_MINT_combinations <- read.csv("1_ModelSimulations/input_files/all_combination_MINT_calibrate.csv")
  
  #Split row_MINT_combinations and run through
  sapply(as.numeric(strsplit(as.character(row_MINT_combinations), ";")[[1]]), function(a){
    
    message(paste0("~~~~ row ", a, " ~~~~"))
    
    this_MINT_combination <- all_MINT_combinations[a, ]
    
    #Set up initial malariasimulation 
    setup_initial <- malaria_simulation_setup(
      
      #Set up changing aspects
      anthropophagy = this_MINT_combination$anthropophagy,
      biting_inbed_indoors = this_MINT_combination$biting_inbed_indoors,
      seasonal = this_MINT_combination$seasonal,
      resistance = this_MINT_combination$resistance,
      itn_use = this_MINT_combination$itn_use,
      irs_use = this_MINT_combination$irs_use,
      net_type_use = "standard",
      
      #These dont change because we're calibrating the baseline
      itn_future = 0,
      irs_future = 0,
      net_type_future = "standard",
      years_run_past = 3,
      years_run_future = 3
      
    )
    
    #Calibrate EIR for the chosen prevalence and time
    calibrated_EIR <- calibrate_malaria_simulation(parameter_set = setup_initial,
                                                   prevalence_target = this_MINT_combination$prevalence, #Prevalence to fit to
                                                   prevalence_time = 730, #Two year prevalence as its in the middle of ITN distribution
                                                   tolerance = 0.02)
    
    calibrated_EIR$row <- a
    
    write.csv(calibrated_EIR,
              paste0(front_part, savemod, "/", a, "_of_", nrow(all_MINT_combinations), "_calibrated.csv"),
              row.names = FALSE)
    
  })
  
}