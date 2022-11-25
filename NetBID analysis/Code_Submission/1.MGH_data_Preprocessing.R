######## LOAD PACKAGES ########

library(tidyverse) # (v1.3.2)
library(biomaRt) # (v2.50.3)
library(Biobase) # (v2.54.0)

rm(list=ls())

######## DATA DOWNLOAD ########

# https://olink.com/application/mgh-covid-19-study/

######## SET PARAMETERS ########

work_pathway <- "C:/Users/vasileioubill95/Desktop/Projects/Covid19Proteins/MGH_Olink_COVID_Apr_27_2021/" # Absolute path where data exist
info_input <- "C:/Users/vasileioubill95/Desktop/Projects/Covid19Proteins/MGH_Olink_COVID_Apr_27_2021/MGH_COVID_Clinical_Info.txt"
expression_input <- "C:/Users/vasileioubill95/Desktop/Projects/Covid19Proteins/MGH_Olink_COVID_Apr_27_2021/MGH_COVID_OLINK_NPX.txt"
ExpressionSet_output <- "C:/Users/vasileioubill95/Desktop/Projects/Covid19Proteins/ProteomicsExpressionSet_D0.RData" # For each Day change the name (*D0.RData, *D3.RData, *D7.RData)

Acuity <- "Acuity_0" # select the day of clinical situation [Acuity_0, Acuity_3, Aquity_7]
Day <- "D0" # select the day of clinical situation [D0=Day 0, D3=Day 3, D7=Day 7]

## set work pathway / folder

setwd(work_pathway)

## Load MGH_COVID_Clinical_Info for phenotypes

pData <- read.table(info_input, sep = ";", header = T)

pData <- pData[c('subject_id','COVID', Acuity)] # select columns of interest

# convert numbers in matrix into situations

# Covid - NonCovid
pData$COVID <- str_replace_all(pData$COVID, "0", "NonCovid")
pData$COVID <- str_replace_all(pData$COVID, "1", "Covid")

# Critical - NonCritical According to the Filbin paper
pData[[Acuity]] <- str_replace_all(pData[[Acuity]], "1", "Critical")
pData[[Acuity]] <- str_replace_all(pData[[Acuity]], "2", "Critical")
pData[[Acuity]] <- str_replace_all(pData[[Acuity]], "3", "NonCritical")
pData[[Acuity]] <- str_replace_all(pData[[Acuity]], "4", "NonCritical")
pData[[Acuity]] <- str_replace_all(pData[[Acuity]], "5", "NonCritical")

# seperate into Covid-NonCovid
pData_noCovid <- filter(pData, COVID == "NonCovid")
pData_Covid <- filter(pData, COVID == "Covid")

# keep subject ID
pData_noCovid <- pData_noCovid$subject_id
pData_Covid <- pData_Covid$subject_id


## Load MGH_COVID_OLINK_NPX data for expression information
count_proteomics <- read.table(expression_input, sep = ";", header = T)

# keep day of interest
count_proteomics <- filter(count_proteomics, Timepoint == Day)

# keep only Covid condition
count_proteomics <- count_proteomics[count_proteomics$subject_id %in% pData_Covid,]

## Make count table 
subject_vector <- sort(unique(count_proteomics$subject_id))
count_process <- count_proteomics[,c("subject_id", "Assay", "NPX")]

for (i in 1:length(subject_vector)) {
  
  if (i == 1) {
    
    count_table <- filter(count_process, subject_id == subject_vector[i] )
    names(count_table)[names(count_table) == "NPX"] <- subject_vector[i]
    rownames(count_table) <- make.names(count_table$Assay, unique = T)
    count_table <- subset(count_table, select=-c(1,2))
  } else {
    
    x <- filter(count_process, subject_id == subject_vector[i] )
    names(x)[names(x) == "NPX"] <- subject_vector[i]
    rownames(x) <- make.names(x$Assay, unique = T)
    x <- subset(x, select=-c(1,2))
    count_table <- merge(count_table,x, by = 0)
    rownames(count_table) <- count_table$Row.names
    count_table <- count_table[,-1]
    print(dim(count_table))
    
    
  }
}

count_table$genes <- rownames(count_table)

## Exclude duplicated proteins

count_table <- filter(count_table, genes != "IL6.1")
count_table <- filter(count_table, genes != "IL6.2")
count_table <- filter(count_table, genes != "IL6.3")

count_table <- filter(count_table, genes != "CXCL8.1")
count_table <- filter(count_table, genes != "CXCL8.2")
count_table <- filter(count_table, genes != "CXCL8.3")

count_table <- filter(count_table, genes != "TNF.1")
count_table <- filter(count_table, genes != "TNF.2")
count_table <- filter(count_table, genes != "TNF.3")

count_table <- subset(count_table, select = -c(genes)) # exclude column genes

expression <- as.matrix(count_table)

pData <- pData[pData$subject_id %in% colnames(expression),]
rownames(pData) <- pData$subject_id

all(colnames(expression)==rownames(pData)) # Check if match perfect, SHOULD BE TRUE
pData <- pData[c('subject_id', Acuity)]

pData_annot <- new("AnnotatedDataFrame", data = pData)
net_eset <- ExpressionSet(assayData = expression, phenoData = pData_annot)

save(net_eset, file = ExpressionSet_output)