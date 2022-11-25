######## LOAD PACKAGES ########

library(NetBID2) # (v2.0.3)
library(tidyverse) # (v1.3.2)

######## SET PARAMETERS ########

rm(list = ls())

#which day
Day_Number <- 0

setwd(sprintf("C:/Users/vasileioubill95/Desktop/Work_environment/D%s_Critical-NonCritical/", Day_Number))
comp_name <- 'Critical.Vs.NonCritical' # Comparison of 2 groups/conditions
sit_0 <- 'NonCritical' # controls
sit_1 <- 'Critical' # cases
Condition <- sprintf("Acuity_%s", Day_Number) # pheno-column that phenotype is specified

current_date <- format(Sys.time(), "%Y-%m-%d")
project_name <- sprintf('driver_%s',current_date)
network.project.name <- sprintf('project_%s',current_date)
network.dir <- paste(getwd(), network.project.name , sep = "/")
project_main_dir <- "Driver_output"

setwd(paste(network.dir, "SJAR", network.project.name, sep = "/"))

file.rename("output_sig", sprintf("output_sig_sjaracne_%s_out_.final", network.project.name))
file.rename("output_tf", sprintf("output_tf_sjaracne_%s_out_.final", network.project.name))

setwd(paste(network.dir, "..", sep = "/"))

# This list object (analysis.par) is an ESSENTIAL variable in driver estimation pipeline
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

# RData saved after QC in the network construction step
load(sprintf('%s/DATA/network.par.Step.exp-QC.RData',network.dir))
analysis.par$cal.eset <- network.par$net.eset

# Get network information
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)

## For backup save file ##
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')

# Merge network first
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)

# Get activity matrix
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')

# Create eset using activity matrix
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

# Create empty list to store comparison result
analysis.par$DE <- list()
analysis.par$DA <- list()

# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info[[Condition]]==sit_1)] # Experiment group
G0  <- rownames(phe_info)[which(phe_info[[Condition]]==sit_0)] # Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name=sit_1,G0_name=sit_0)
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name=sit_1,G0_name=sit_0)
# Save comparison result to list element in analysis.par, with comparison name
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

# Save Step 3 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-DA')

## write DE/DA results into csv format
write.csv(analysis.par$DE$Critical.Vs.NonCritical, file = paste(network.dir, "../Driver_output", project_name, "Data/DE_res.csv", sep = "/"), row.names = F, quote = F)
write.csv(analysis.par$DA$Critical.Vs.NonCritical, file = paste(network.dir, "../Driver_output", project_name, "Data/DA_res.csv", sep = "/"), row.names = F, quote = F)

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='transcript',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(analysis.par$DE) # Users can use index or name to get target ones
# Prepare the conversion table (OPTIONAL)
use_genes <- unique(c(analysis.par$merge.network$network_dat$source.symbol,analysis.par$merge.network$network_dat$target.symbol))
transfer_tab <- get_IDtransfer2symbol2type(from_type = 'hgnc_symbol',use_genes=use_genes)
analysis.par$transfer_tab <- transfer_tab
# Creat the final master table
analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp,DE=analysis.par$DE,DA=analysis.par$DA,
                                                  target_list=analysis.par$merge.network$target_list,
                                                  tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                                  main_id_type='hgnc_symbol')

# Path and file name of the output EXCEL file
out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)

# Save the final master table as EXCEL file
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file)

# Save Step 4 analysis.par as RData, ESSENTIAL
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')
