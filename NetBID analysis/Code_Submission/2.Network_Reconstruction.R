######## LOAD PACKAGES ########

library(NetBID2) # (v2.0.3)
library(tidyverse) # (v1.3.2)
library(icesTAF) # (v4.0.0)

rm(list = ls())

######## SET PARAMETERS ########

Day_Number <- 0

setwd("C:/Users/vasileioubill95/Desktop/") # set work directory / environment
mkdir(paste(getwd(), sprintf('Work_environment/D%s_Critical-NonCritical', Day_Number), sep = "/"))
setwd("Work_environment/D0_Critical-NonCritical")
load(sprintf("C:/Users/vasileioubill95/Desktop/Projects/Covid19Proteins/ProteomicsExpressionSet_D%s.RData", Day_Number)) # load ExpressionSet format file of COVID proteomics with absolute directory

# Define main working directory and project name
project_main_dir <- './' # user defined main directory for the project, one main directory could have multiple project folders, distinguished by project name.
current_date <- format(Sys.time(), "%Y-%m-%d") # optional, if user like to add current date to name the project folder.
project_name <- sprintf('project_%s',current_date) # project name for the project folders under main directory.

# This list object (network.par) is an ESSENTIAL variable in network construction pipeline
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

# Add the variable into network.par
network.par$net.eset <- net_eset

## For backup save file ##
NetBID.saveRData(network.par = network.par,step='exp-load')

######## PERFORM QUALITY CONTROL ########

# Get the expression matrix from ExpressionSet object
mat <- exprs(network.par$net.eset)

# Filter out genes with very low expression values (bottom 5%) in most samples (more than 90%).
choose1 <- apply(mat<= quantile(mat, probs = 0.05), 1, sum)<= ncol(mat) * 0.90
print(table(choose1))
mat <- mat[choose1,]

# Update eset with normalized expression matrix
net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                          feature_info=fData(network.par$net.eset)[rownames(mat),],
                          annotation_info=annotation(network.par$net.eset))

# Update network.par with new eset
network.par$net.eset <- net_eset

## For backup save file ##
NetBID.saveRData(network.par = network.par,step='exp-QC')

######## PREPARE FILE TO RUN SJARACNe ########

# Load database
db.preload(use_level='transcript',use_spe='human',update=FALSE)

# Converts gene ID into the corresponding TF/SIG list
use_gene_type <- 'hgnc_symbol' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)

# Select samples for analysis
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) # here is using all samples, users can modify
prj.name <- network.par$project.name # if use different samples, need to change the project name
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0.5,IQR.loose_thre = 0.5,
                 SJAR.project_name=prj.name,SJAR.main_dir=network.par$out.dir.SJAR)
