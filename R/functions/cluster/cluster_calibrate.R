# row_MINT_combinations = run_all$row_MINT_combinations[1]
# savemod = run_all$savemod[1]

cluster_calibrate <- function(SSA_adm1,
                              savemod){
  
  #Front part of the file repo and then create a folder if it doesnt exist
  front_part <- "output/EIR_calibrate/"
  if(!dir.exists(paste0(front_part, savemod, "/"))) dir.create(paste0(front_part, savemod, "/"), recursive = T)
  
  #Read in and subset parameters for baseline 
  all_SSA_location <- read.csv("data/SSA_adm1_characteristics/MAP_2020_subset.csv")
  resistance_location <- read.csv("data/resistance/resistance_updated.csv") %>%
    dplyr::filter(prob == 0.5 & year == 2020)
  
  
  #Split row_MINT_combinations and run through
  sapply(strsplit(as.character(SSA_adm1), ";")[[1]], function(a){
    
    message(paste0("~~~~ row ", a, " ~~~~"))
    
    this_SSA_combination <- all_SSA_location[all_SSA_location$GID_1 == a, ]
    #The resistance comes from old shapefiles, need to take country average to work it
    this_resistance <- resistance_location[which(resistance_location$country_name == this_SSA_combination$NAME_0), ]
    
    #Set up initial malariasimulation 
    setup_initial <- malaria_simulation_setup(
      
      SSA_adm1 = a,
      
      #Set up changing aspects
      seasonal = T,
      resistance = median(this_resistance$resistance_qs),
      itn_use = this_SSA_combination$ITN_use,
      irs_use = this_SSA_combination$IRS_use,
      net_type_use = "standard",
      species_proportion = paste(this_SSA_combination[, grepl("abundance", colnames(this_SSA_combination))], collapse = ";"),
      
      #These dont change because we're calibrating the baseline
      itn_future = this_SSA_combination$ITN_use,
      irs_future = this_SSA_combination$IRS_use,
      net_type_future = "standard",
      years_run_past = 6,
      years_run_future = 6
      
    )
    
    #Calibrate EIR for the chosen prevalence and time
    calibrated_EIR <- calibrate_malaria_simulation(parameter_set = setup_initial,
                                                   prevalence_target = this_SSA_combination$prevalence, #Prevalence to fit to
                                                   prevalence_time = 730, #Two year prevalence as its in the middle of ITN distribution
                                                   tolerance = 0.02)
    
    calibrated_EIR$row <- a
    
    write.csv(calibrated_EIR,
              paste0(front_part, savemod, "/", a, "_of_", nrow(all_MINT_combinations), "_calibrated.csv"),
              row.names = FALSE)
    
  })
  
}