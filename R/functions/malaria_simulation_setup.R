# # # #Test values
# seasonal = 1
# resistance = 40
# itn_use = 0.2
# irs_use = 0
# prevalence_target = 0.1
# net_type_use = "standard"
# itn_future = 0.2
# irs_future = 0
# lsm_future = 0.05
# net_type_future = "PBO"
# years_run_future = 12
# itn_future_distribution = 3
# irs_future_distribution = 1
# species_proportion = "0.5;0.25;0.25"


malaria_simulation_setup <- function(
  
  #Set up baseline levels
  seasonal,                #1 = perennial, 2 = seasonal vector populations
  resistance,              #Pyrethroid resistance, percent surviving takes in 0/20/40/60/80/100
  
  #Set up intervention values and how long they run
  itn_use,                 #Historical itn use
  irs_use,                 #Historical irs use
  net_type_use,            #Takes in standard, PBO or IG2
  years_run_past,          #The number of years to run as "burn in"
  
  itn_future,             #Set up the future itn coverage
  irs_future,             #Set up the future irs coverage
  lsm_future,             #Set up the future LSM coverage (reduction in adult emergence)
  net_type_future,        #Set up the future itn type
  years_run_future,       #How many years post intervention do you run 
  uncertainty_draw = "median", #Which uncertainty draw to take for net/irs parameters, accepts low/median/high
  
  itn_future_distribution = 3, #How often (in years) you want to distribute nets with future coverage
  irs_future_distribution = 1, #How often (in years) you want to implement irs with future coverage
  
  #These options will never change but are included for future use
  species_proportion = "0.5;0.25;0.25" #Takes in the proportion of gamb, arab, fun (in that order) separated by ;, must add to 1
  
){
  
  #Record all entry data
  input_data <- data.frame(
    seasonal,                #1 = perennial, 2 = seasonal vector populations
    resistance,              #Pyrethroid resistance, percent surviving takes in 0/20/40/60/80/100
    
    #Set up intervention values and how long they run
    itn_use,                 #Historical itn use
    irs_use,                 #Historical irs use
    net_type_use,            #Takes in standard, PBO or IG2
    years_run_past,          #The number of years to run as "burn in"
    
    itn_future,             #Set up the future itn coverage
    irs_future,             #Set up the future irs coverage
    net_type_future,        #Set up the future itn type
    years_run_future,       #How many years post intervention do you run 
    
    itn_future_distribution = itn_future_distribution, #How often (in years) you want to distribute nets with future coverage
    irs_future_distribution = irs_future_distribution, #How often (in years) you want to implement irs with future coverage
    
    uncertainty_draw = uncertainty_draw, #Which uncertainty draw to take for net/irs parameters, accepts low/median/high
    
    #These options will never change but are included for future use
    species_proportion = species_proportion #Takes in the proportion of gamb, arab, fun (in that order) separated by ;, must add to 1
  )
  
  # Load data ---------------------------------------------------------------
  #Load in the site file for seasonal parameters and subset to what we want
  site <- rio::import("1_ModelSimulations/input_files/Sites/site_file_Global_malariasimulation.csv")
  site_use <- site[which(site$seasonality == seasonal), ]

  #Load in and subset resistance
  resistance_before <- rio::import(paste0("1_ModelSimulations/input_files/input_params_v1/Pyrethroid_resistance_", 
                                      resistance, 
                                      "_input_data", 
                                      ifelse(net_type_use == "standard", "", paste0("_", net_type_use)),
                                      ".csv"))[1, ]
  
  resistance_after <- rio::import(paste0("1_ModelSimulations/input_files/input_params_v1/Pyrethroid_resistance_", 
                                          resistance, 
                                          "_input_data", 
                                          ifelse(net_type_use == "standard", "", paste0("_", net_type_future)),
                                          ".csv"))[1, ]
  
  # Set up malariasimulation baseline inputs --------------------------------
  #Set up the basic info for malariasimulation
  year <- 365
  month <- 30                 
  sim_length <- 10 * year            #Running for 10 years
  human_population <- 10000          #Human population of 10,000 - needs to be of a sufficient size to deal with stochasticity
  
  #Set up initial malariasimulation parameters
  simparams <- get_parameters(
    list(
      human_population = human_population,
      
      prevalence_rendering_min_ages = c(0,  0) * 365, #This sets up the minimum ages of prevalence variables we want to output
      prevalence_rendering_max_ages = c(5, 100) * 365, #This sets up the maximum ages, so we will have 2 variables for the 0-5 prevalence and 0-100
      
      clinical_incidence_rendering_min_ages = c(0,  0) * 365, #This sets up the minimum ages of incidence variables we want to output
      clinical_incidence_rendering_max_ages = c(5, 100) * 365, #This sets up the maximum ages, so we will have 2 variables for the 0-5 incidence and 0-100
      
      model_seasonality = TRUE, #Seasonablity is enabled, we then feed in the seasonality values
      g0 = site_use$seasonal_a0,
      g = c(site_use$seasonal_a1, site_use$seasonal_a2, site_use$seasonal_a3),
      h = c(site_use$seasonal_b1, site_use$seasonal_b2, site_use$seasonal_b3),
      
      individual_mosquitoes = FALSE, #This uses the compartmental model rather than individual mosquitoes (faster and less stochastic)
      
      bednets = TRUE
    )
  )
  
  #This is used to optimally time the implementation of IRS
  peak <- peak_season_offset(simparams) 
  
  #Input vector parameters into simparms
  simparams <- set_species(simparams,
                           list(gamb_params, arab_params, fun_params),
                           as.numeric(strsplit(species_proportion, ";")[[1]]))  #These proportions never change
  
  # Set up interventions ----------------------------------------------------
  # Set up drugs ------------------------------------------------------------
  
  #Now specify interventions
  # In the MINT tool we simulate the 42% of people receive treatment and 24% of that treatment is ACT
  # These are currently fixed but would be good to have the option to make these specific
  # set treatment
  simparams <- set_drugs(simparams, list(AL_params,    ## whichever is ACT drug
                                         SP_AQ_params)) ## whichever is non-ACT drug
  simparams <- set_clinical_treatment(simparams, 
                                      drug = 1,
                                      time = 100,
                                      coverage = 0.242248865)    # currently these are restricted but it would be a future hope to allow them to change
  simparams <- set_clinical_treatment(simparams, 
                                      drug = 2,
                                      time = 100,
                                      coverage = 0.17369117) 
  
  years_of_past_itn <- ceiling(years_run_past/3)
  years_of_future_itn <- ceiling(years_run_future/itn_future_distribution)
  
  # Set up bednets ----------------------------------------------------------
  simparams <- set_bednets(
    simparams,
    
    timesteps = c(seq(1, years_run_past, by = 3) * 365,                       #Baseline distribution of nets is every 3 years, with 3 distributions
                  seq(years_run_past + itn_future_distribution,               #Future distribution of nets depends on how often they are distributed and how long the simulation runs for
                      years_run_past + years_run_future,
                      by = itn_future_distribution) * 365),

    coverages = c(rep(itn_use, years_of_past_itn),        #Baseline value
                  rep(itn_future, years_of_future_itn)),      #Future coverage
    
    retention = 5 * year, 
    ## This is currently fixed but would 
    ## be a parameter we would like to be able to vary (between 1.5 and 15 perhaps)
    ## and potentially update once we have the option to use Amelia Bertozzi-Villa's logistic decay
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito)
    ## our estimates do not distinguish species currently, a potential thing we would update but depends on data
    dn0 = matrix(c(rep(resistance_before[, grepl("kill_gamb_ss_0", colnames(resistance_before))], years_of_past_itn), 
                   rep(resistance_after[, grepl("kill_gamb_ss_1", colnames(resistance_after))], years_of_future_itn),
                   rep(resistance_before[, grepl("kill_arab_0", colnames(resistance_before))], years_of_past_itn), 
                   rep(resistance_after[, grepl("kill_arab_1", colnames(resistance_after))], years_of_future_itn),
                   rep(resistance_before[, grepl("kill_fun_0", colnames(resistance_before))], years_of_past_itn), 
                   rep(resistance_after[, grepl("kill_fun_1", colnames(resistance_after))], years_of_future_itn)), 
                 nrow = years_of_past_itn + years_of_future_itn,  #The number of timepoints of ITN implementation
                 ncol = 3), #The number of species
    ## These change with the level of resistance
    
    rn = matrix(c(rep(resistance_before[, grepl("repel_gamb_ss_0", colnames(resistance_before))], years_of_past_itn), 
                  rep(resistance_after[, grepl("repel_gamb_ss_1", colnames(resistance_after))], years_of_future_itn),
                  rep(resistance_before[, grepl("repel_arab_0", colnames(resistance_before))], years_of_past_itn), 
                  rep(resistance_after[, grepl("repel_arab_1", colnames(resistance_after))], years_of_future_itn),
                  rep(resistance_before[, grepl("repel_fun_0", colnames(resistance_before))], years_of_past_itn), 
                  rep(resistance_after[, grepl("repel_fun_1", colnames(resistance_after))], years_of_future_itn)), 
                nrow = years_of_past_itn + years_of_future_itn,  #The number of timepoints of ITN implementation
                ncol = 3), #The number of species    ## These change with the level of resistance
    
    rnm = matrix(c(.24, .24, .24), nrow = years_of_past_itn + years_of_future_itn, ncol = 3),
    ## This is fixed ... something we prob need to investigate statistically at some point!
    
    gamman = c(rep(resistance_before$itn_halflife_0, years_of_past_itn),
               rep(resistance_after$itn_halflife_1, years_of_future_itn))
    ## These change with the level of resistance
    
  )
  
  # Setup IRS --------------------------------------------------------------
  years_of_future_irs <- ceiling(years_run_future/irs_future_distribution)
  years_of_past_irs <- ceiling(years_run_past/1)
  
  simparams <- set_spraying(
    simparams, ## details same as pyrethroid only nets
    
    timesteps = c(seq(1, years_run_past, by = 1) * 365 + peak,                       #Baseline distribution of nets is every 3 years, with 3 distributions
                  seq(years_run_past + irs_future_distribution,                      #Future distribution of nets depends on how often they are distributed and how long the simulation runs for
                      years_run_past + years_run_future,
                      by = irs_future_distribution) * 365 + peak),
    
    coverages = c(rep(irs_use, years_of_past_irs),        #Baseline value
                  rep(irs_future, years_of_future_irs)),      #Future coverage
    
    ##These will change for upper and lower impacts (in min and max files) - 
    ls_theta = matrix(rep(resistance_before$irs_decay_mort1_1, (years_of_future_irs + years_run_past) * 3),  nrow = years_of_future_irs + years_run_past, ncol = 3),
    ls_gamma = matrix(rep(resistance_before$irs_decay_mort2_1, (years_of_future_irs + years_run_past) * 3), nrow = years_of_future_irs + years_run_past, ncol = 3),
    ks_theta = matrix(rep(resistance_before$irs_decay_succ1_1, (years_of_future_irs + years_run_past) * 3), nrow = years_of_future_irs + years_run_past, ncol = 3),
    ks_gamma = matrix(rep(resistance_before$irs_decay_succ2_1, (years_of_future_irs + years_run_past) * 3),  nrow = years_of_future_irs + years_run_past, ncol = 3),
    ms_theta = matrix(rep(resistance_before$irs_decay_det1_1, (years_of_future_irs + years_run_past) * 3), nrow = years_of_future_irs + years_run_past, ncol = 3),
    ms_gamma = matrix(rep(resistance_before$irs_decay_det2_1, (years_of_future_irs + years_run_past) * 3), nrow = years_of_future_irs + years_run_past, ncol = 3) 
    ## we could change product year to year (which is not possible in MINT currently)
    ## so these rows could be distinct estimates in the future
  )
  
  simparams$simulation_length <- (years_run_future + years_run_past) * 365
  simparams$input_data <- input_data
  
  simparams
  
}
