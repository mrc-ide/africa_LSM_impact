# rm(list = ls(all = TRUE))
options(
  didehpc.username = "ah1114")

drat:::add("mrc-ide")

## information on my computer
malaria <- didehpc::path_mapping("Malaria", "U:", "//fi--didef3.dide.ic.ac.uk/Malaria", "U:")

config <- didehpc::didehpc_config(shares = malaria, use_rrq =  FALSE, 
                                  cluster = "fi--dideclusthn", cores = 1,
                                  parallel = FALSE)

packages <- c("malariasimulation", "malariaEquilibrium", "reshape2", "rio", "cali", "individual")

sources = c("1_ModelSimulations/cluster/source_files.R")

src <- conan::conan_sources(packages = packages, c("mrc-ide/cali", "mrc-ide/malariasimulation", "mrc-ide/malariaEquilibrium"))

ctx <- context::context_save(path = "new_contexts", 
                             packages = packages, 
                             sources = sources,
                             package_sources = src)

obj = didehpc::queue_didehpc(ctx, config = config)

#Set up runs
all_MINT_combinations <- read.csv("1_ModelSimulations/input_files/all_combination_MINT_calibrate.csv")
chunk_combo <- split(1:nrow(all_MINT_combinations), ceiling(seq_along(1:nrow(all_MINT_combinations))/20))

run_all <- data.frame(row_MINT_combinations = sapply(chunk_combo, function(x) paste(x, collapse = ";")),
                      savemod = "first_run",
                      stringsAsFactors = FALSE)

grp <- obj$enqueue_bulk(run_all,
                        cluster_calibrate, do_call = TRUE, progress = TRUE)

#Now relaunch things that have failed
all_file_ran <- list.files(paste0("1_ModelSimulations/output/EIR_calibrate/", run_all$savemod[1], "/"),
                           pattern = ".csv")
these_ran <- as.numeric(gsub("_of_5040_calibrated.csv", "", all_file_ran))
run_these_now <- (1:nrow(all_MINT_combinations))[-these_ran]
chunk_combo_redo <- split(run_these_now, ceiling(seq_along(run_these_now)/5))

run_all_redo <- data.frame(row_MINT_combinations = sapply(chunk_combo_redo, function(x) paste(x, collapse = ";")),
                      savemod = "first_run",
                      stringsAsFactors = FALSE)

nrow(run_all_redo)

grp <- obj$enqueue_bulk(run_all_redo,
                        cluster_calibrate, do_call = TRUE, progress = TRUE)

