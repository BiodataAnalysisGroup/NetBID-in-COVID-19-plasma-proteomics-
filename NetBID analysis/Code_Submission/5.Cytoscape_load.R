######## Launch CYTOSCAPE ########

######## LOAD PACKAGES FOR CYTOSCAPE ########

library(RCy3) # (v2.14.2)
library(igraph) # (v1.3.5)
library(plyr) # (v1.8.8)
library(RColorBrewer) # (v1.1.3)
library(tidyverse) # (v1.3.2)

Day_Number <- 0
nodes_path <- "C:/Users/vasileioubill95/Desktop/Work_environment/D0_Critical-NonCritical/Driver_output/driver_2022-11-24/DATA/DA_res.csv"
sj_sig <- "C:/Users/vasileioubill95/Desktop/Work_environment/D0_Critical-NonCritical/project_2022-11-24/SJAR/project_2022-11-24/output_sig_sjaracne_project_2022-11-24_out_.final/consensus_network_ncol_.txt"

## load information for Nodes
nodes <- read.csv(nodes_path, header = T)
colnames(nodes)[1] <- "id"

nodes$id <- str_replace_all(nodes$id, "_SIG", "")
nodes$id <- str_replace_all(nodes$id, "_TF", "")
nodes <- unique(nodes)

## load information for Edges
network <- read.table(sj_sig, header = T)

network <- filter(network, p.value <= 0.05)
network <- filter(network, MI >= 0.6)

title <- sprintf("Network_D%s", Day_Number)

createNetworkFromDataFrames(nodes=nodes,edges = network, title=title, collection="Networks")

## Bofere analysis it is necessary to remove duplicated edges by Edit -> Remove Duplicate Edges
## and Select Network_Dx , click Ignore direction , press OK (Cytoscape V.3.9.1)
