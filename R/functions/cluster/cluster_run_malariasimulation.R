# first_run = 4564
# itn_future = 0.8
# irs_future = 0.5
# net_type_future = "PBO"
# itn_future_distribution = 3
# irs_future_distribution = 1
# savemod = "first_run"
# summarise_annual = T


cluster_run_malariasimulation <- function(row_SSA_adm1,
                                          row_future_interventions,
                                          years_run_future,       #How many years post intervention do you run 
                                          savemod,
                                          itn_future_distribution = 3, #How often (in years) you want to distribute nets with future coverage
                                          irs_future_distribution = 1, #How often (in years) you want to implement irs with future coverage
                                          summarise_annual = T         #Pair down and summarise data annually for MINT?
){
  
  #Front part of the file repo and then create a folder if it doesnt exist
  front_part <- paste0("1_ModelSimulations/output/malariasimulation_runs/", savemod, "/")
  
  if(!dir.exists(paste0(front_part, "all/"))) dir.create(paste0(front_part, "all/"), recursive = T)
  if(!dir.exists(paste0(front_part, "summary/"))) dir.create(paste0(front_part, "summary/"), recursive = T)
  
  #Read in and subset parameters for baseline 
  all_MINT_combinations <- read.csv("1_ModelSimulations/input_files/all_combination_MINT_calibrate.csv")
  all_MINT_combinations_interventions <- read.csv("1_ModelSimulations/input_files/all_future_MINT_intervention.csv")
  
  sapply(as.numeric(strsplit(as.character(row_MINT_combinations), ";")[[1]]), function(a){
    
    sapply(as.numeric(strsplit(as.character(row_MINT_future_interventions), ";")[[1]]), function(b){
      
      message(paste0("~~~~ baseline row ", a, " ~~~~"))
      message(paste0("~~~~ intervention row ", b, " ~~~~"))
      
      this_MINT_combination <- all_MINT_combinations[a, ]
      this_MINT_intervention <- all_MINT_combinations_interventions[b, ]
      
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
        
        #These change depending on user input of the future
        itn_future = this_MINT_intervention$itn_future,
        irs_future = this_MINT_intervention$irs_future,
        net_type_future = this_MINT_intervention$net_type_future,
        itn_future_distribution = itn_future_distribution,
        irs_future_distribution = irs_future_distribution,
        
        uncertainty_draw = this_MINT_intervention$uncertainty_draw,
        
        #How long we run for 
        years_run_past = 6,
        years_run_future = 6
        
      )
      
      #Read in EIR calibration
      EIR_calibration_done <- read.csv(list.files(paste0("1_ModelSimulations/output/EIR_calibrate/", savemod, "/"),
                                                  pattern = paste0(a, "_of_"), full.names = T))
      
      #Set equilibrium for malariasimulation parameters
      setup_calibrated <- set_equilibrium(setup_initial, EIR_calibration_done$EIR_fit)
      
      #Run model
      model_ran <- run_simulation(setup_calibrated$simulation_length, setup_calibrated) 
      
      #Set up additional columns and subset
      model_ran$prevalence_0_5 <- model_ran$n_detect_0_1825/model_ran$n_0_1825
      model_ran$prevalence_0_100 <- model_ran$n_detect_0_36500/model_ran$n_0_36500
      model_ran$year <- floor((model_ran$timestep-1)/365) - setup_calibrated$input_data$years_run_past
      
      model_ran$savemod <- savemod
      model_ran$row_total <- paste0(a, "_of_", nrow(all_MINT_combinations))
      
      if(summarise_annual == T){
        
        model_ran_subset <- model_ran %>%
          group_by(year) %>%
          summarise(savemod = unique(savemod),
                    row_total = unique(row_total),
                    prevalence_0_100 = mean(prevalence_0_100),
                    prevalence_0_5 = mean(prevalence_0_5),
                    incidence_0_100 = sum(p_inc_clinical_0_36500),
                    incidence_0_5 = sum(p_inc_clinical_0_1825)) %>%
          mutate(population = setup_initial$human_population)
        
      } else {
        
        model_ran_subset <- model_ran
        
      }
      
      write.csv(model_ran_subset,
                paste0(front_part,
                       ifelse(summarise_annual == T, "summary", "all"),
                       "/", a, "_of_", nrow(all_MINT_combinations), "_malariasimulation_output.csv"),
                row.names = FALSE)
      
    })
    
  }, simplify = FALSE)
  
}