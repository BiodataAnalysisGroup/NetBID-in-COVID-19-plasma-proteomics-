######## LOAD PACKAGES ########

library(NetBID2) # (v2.0.3)
library(tidyverse) # (v1.3.2)

######## SET PARAMETERS ########

rm(list=ls())

Day_Number <- 0

setwd("C:/Users/vasileioubill95/Desktop/Work_environment/D0_Critical-NonCritical/")

comp_name <- 'Critical.Vs.NonCritical'
anno_control <- 'high in NonCritical'
anno_cases <- 'high in Critical'

# Give file path to reload `ms-tab` RData from driver inference step
analysis.par <- list()

current_date <- format(Sys.time(), "%Y-%m-%d")
project_name <- sprintf('driver_%s',current_date)
network.project.name <- sprintf('project_%s',current_date)
network.dir <- paste(getwd(), network.project.name , sep = "/")

analysis.par$out.dir.DATA <- paste(network.dir, "../Driver_output", project_name, "Data", sep = "/")
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')

############### Part I: More details about the top drivers  ###############

### QI.1: How to get the top drivers with significant differential activity (DA) in the comparison between G4 vs. other subtypes ?

ms_tab <- analysis.par$final_ms_tab ## get the master table data frame

ms_tab$originalID_label <- str_replace_all(ms_tab$originalID_label, "_SIG", "")
ms_tab$originalID_label <- str_replace_all(ms_tab$originalID_label, "_TF", "")
ms_tab <- ms_tab[!duplicated(ms_tab$originalID_label),]

sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='originalID_label',logFC_col=sprintf('logFC.%s_DA',comp_name),
                               Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0,Pv_thre=0.05,
                               main='DA Critical - NonCritical',show_label=TRUE,
                               pdf_file=sprintf('%s/vocalno_label_DA.pdf',analysis.par$out.dir.PLOT),label_cex = 1)

# Get the DE data frame of target genes
DE <- analysis.par$DE[[comp_name]]
driver_list <- rownames(sig_driver) # The rownames is the originalID_label

negative_drivers <- filter(sig_driver, logFC.Critical.Vs.NonCritical_DA < 0)
driver_list_negative <- rownames(negative_drivers)

draw.GSEA.NetBID(DE=DE,profile_col='logFC',profile_trend='pos2neg',name_col='ID',
                 driver_list = driver_list_negative,
                 show_label=ms_tab[driver_list_negative,'originalID_label'],
                 driver_DA_Z=ms_tab[driver_list_negative,sprintf('Z.%s_DA',comp_name)],
                 driver_DE_Z=ms_tab[driver_list_negative,sprintf('Z.%s_DE',comp_name)],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = anno_control,right_annotation = anno_cases,
                 main= comp_name,target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_negative.pdf',analysis.par$out.dir.PLOT))

positive_drivers <- filter(sig_driver, logFC.Critical.Vs.NonCritical_DA > 0)
driver_list_positive <- rownames(positive_drivers)

draw.GSEA.NetBID(DE=DE,profile_col='logFC',profile_trend='pos2neg',name_col='ID',
                 driver_list = driver_list_positive,
                 show_label=ms_tab[driver_list_positive,'originalID_label'],
                 driver_DA_Z=ms_tab[driver_list_positive,sprintf('Z.%s_DA',comp_name)],
                 driver_DE_Z=ms_tab[driver_list_positive,sprintf('Z.%s_DE',comp_name)],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = anno_control,right_annotation = anno_cases,
                 main= comp_name,target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_positive.pdf',analysis.par$out.dir.PLOT))

# Download gene sets from MSigDB and save as RData, creat a global variable all_gs2gene
# Keep only top 30 positive and top 30 negative drivers
gs.preload(use_spe='Homo sapiens',update=FALSE)

sig_driver_up <- filter(sig_driver, logFC.Critical.Vs.NonCritical_DA > 0)
sig_driver_up <- sig_driver_up[order(sig_driver_up$P.Value.Critical.Vs.NonCritical_DA),]
sig_driver_up <- sig_driver_up[1:30,]
driver_list_up <- rownames(sig_driver_up)

sig_driver_down <- filter(sig_driver, logFC.Critical.Vs.NonCritical_DA < 0)
sig_driver_down <- sig_driver_down[order(sig_driver_down$P.Value.Critical.Vs.NonCritical_DA),]
sig_driver_down <- sig_driver_down[1:30,]
driver_list_down <- rownames(sig_driver_down)

driver_list_all <- as.vector(rbind(driver_list_up, driver_list_down))

#Get ID conversion table
transfer_tab <- analysis.par$transfer_tab
gs.preload(use_spe='Homo sapiens',update=FALSE)

# Bubble Plot to show target genes enriched biological functions

draw.bubblePlot(driver_list= driver_list_all,show_label=ms_tab[driver_list_all,'originalID_label'],
                Z_val=ms_tab[driver_list_all,sprintf('Z.%s_DA',comp_name)],
                driver_type=ms_tab[driver_list_all,'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=10,max_gs_size=5000,use_gs=c('CP:KEGG','CP:BIOCARTA','H','CP:REACTOME'),
                top_geneset_number=30,top_driver_number=60,
                pdf_file = sprintf('%s/bubblePlot_paths.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')

dev.off()
