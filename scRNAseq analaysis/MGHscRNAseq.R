#MGH scRNA-seq dataset from the Villani group
#https://www.covid19cellatlas.org/index.patient.html
library(Seurat)
library(SingleCellExperiment)
library(dior)
library(ggplot2)
library(tidyseurat)
library(Seurat)
library(tidyseurat)
library(biomaRt)
library(org.Hs.eg.db)
library(msigdbr)
library(VAM)
library(dplyr)
library(patchwork)
library(msigdbr)
library(org.Hs.eg.db)
library(Nebulosa) 


#Convert .h5ad scRNA-seq data to a Seurat object
#Import from scDIOR pipeline (h5ad -----> seurat object)
adata <- read_h5(file = './MGHdata.h5', #input from Scanpy pipeline
                 assay.name = 'RNA', 
                 target.object = 'seurat')
pbmc <- adata
rm(adata)

#select COVID-19
pbmc <- pbmc %>% filter(covid == "1")

#select Day 0
pbmc <- pbmc %>% filter(time_point == "D0") #D3 for Day 3, D7 for Day 7
pbmc_s <- pbmc %>% filter(who_0 %in% c("0", "1", "2")) #severe cases
pbmc_ns <- pbmc %>% filter(!(who_0 %in% c("0", "1", "2"))) #non-severe cases

#Basic preprocessing
D0_list <- list(pbmc_s, pbmc_ns)
for (i in 1:length(D0_list)) {
  D0_list[[i]] <- subset(D0_list[[i]], subset = n_features > 200 & n_features < 2500 & mito_pct < 10)
  D0_list[[i]] <- SCTransform( D0_list[[i]], verbose=T, vars.to.regress = "mito_pct")
  D0_list[[i]] <- RunPCA(D0_list[[i]], features = VariableFeatures(object = D0_list[[i]]))
  D0_list[[i]] <-FindNeighbors(D0_list[[i]], dims = 1:15)
  D0_list[[i]] <-FindClusters(D0_list[[i]], resolution = 0.5)
  D0_list[[i]]  <- RunUMAP(D0_list[[i]] , dims = 1:15)
}

for (i in 1:length(D0_list)) {                                          #sanity check
  print(DimPlot(D0_list[[i]], reduction = "umap", group.by = "covid"))
  print(DimPlot(D0_list[[i]], reduction = "umap", group.by = "time_point"))
  print(DimPlot(D0_list[[i]], reduction = "umap", group.by = "who_0"))
  print(DimPlot(D0_list[[i]], reduction = "umap", group.by = "patient_id")) 
}

#Split list
severe <- D0_list[[1]]
non_sev <- D0_list[[2]]
# saveRDS(severe, "D0_MGH_severe.rds")
# saveRDS(non_sev, "D0_MGH_NONsevere.rds")

#Change idents to cell annotations
Idents(severe) <- 'Annotation'
Idents(non_sev) <- 'Annotation'

#Trace NetBID drivers
#Day 0
drivers_D0 <- DA_D0 %>% filter(DA_D0$adj.P.Val < 0.05) #DA_D3 for Day 3, DA_D7 for Day 7
drivers_D0_up <- drivers_D0 %>% filter(drivers_D0$logFC > 0)
drivers_D0_down <- drivers_D0 %>% filter(drivers_D0$logFC < 0)

#positive drivers
df <- data.frame(drivers_D0_up$ID)
dfup <- df %>% filter(!(drivers_D0_up.ID %in% c('AGR2', 'KRT19', 'SFTPA2', 
                                                  'DDAH1', 'CALCA', 'SRP14'))) #not present in the scRNA-seq dataset
# For Day 3: dfup <- df %>% filter(!(drivers_D3_up.ID %in% c('THBS2', 'SPP1', 'PLA2G2A', 'CALCA', 'GDF15', 
#                                                 'IGSF3', 'TNFRSF11B', 'IL1RL1', 'STC1', 'ANGPTL4', 'CHI3L1', 'CXCL13', 'FUT3_FUT5', 'SDC1', 'CCL20')))
# For Day 7: dfup <- c("NID1", "DAG1", "HS6ST1", 
# "IL18R1", "B4GALT1", "LRPAP1", "GALNT2", "HGF", "MFGE8", "CD14", "PVR",
# "CTSB", "SIGLEC10", "CLEC11A", "CXCL16")

dfup <- as.character(dfup$drivers_D0_up.ID)
dfup[1:15] #adjust to your liking
p1 <- DotPlot(severe, features = dfup[1:15], scale.max = 60, col.min = 0, cluster.idents = T)
p2 <- DotPlot(non_sev, features = dfup[1:15], col.min = 0, cluster.idents = T)
p1 / p2

#negative drivers
df_n <- data.frame(drivers_D0_down$ID)
dfdown <- df_n %>% filter(!(drivers_D0_down.ID %in% c('PON3','KIFBP','CRH', 'FETUB', 'GH1', 'BGN', 'OBP2B')))
# For Day 3: dfdown <- df %>% filter(!(drivers_D3_down.ID %in% c('KLK8', 'FABP2', 'SPINK6', 'FETUB', 
#                                                     'OBP2B', 'CTSV', 'WFIKKN2', 'CRH', 'AFP', 
#                                                     'KLK14', 'ARG1', 'AHSP', 'PKLR', 'RBP2', 'MEP1B', 'CA1', 'DDC')))
# For Day 7: dfdown <- c("GPA33", "BRK1", "CA2", "HAGH", "LXN", "HMBS", "CASP2", "GLO1", "PPME1",
#             "IGFBP3", "AARSD1", "NSFL1C", "APRT", "XPNPEP2", "GMPR")
dfdown <- as.character(dfdown$drivers_D0_down.ID)
dfdown[1:23] #adjust to your liking
p3 <- DotPlot(severe, features = dfdown[1:15],scale.max = 60, col.min = 0, cluster.idents = T)
p4 <- DotPlot(non_sev, features = dfdown[1:15], col.min = 0, cluster.idents = T)
p3 / p4


#Pathway enrichment-severe COVID-19
#Using VAM package
#We need Ensembl IDs
Srt_s <- severe
genes_Seurat <- data.frame(rownames(x = Srt_s))
colnames(genes_Seurat) <- "external_gene_name"
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
annot<-getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
             filters = "external_gene_name",
             values = genes_Seurat$external_gene_name,
             mart=ensembl)
Biomart_input <- annot %>% distinct(external_gene_name, .keep_all = TRUE) #remove some duplicate genes
non_annot <- as.data.frame(genes_Seurat[!genes_Seurat$external_gene_name %in% annot$external_gene_name,])

annot <- Biomart_input %>% dplyr::select(external_gene_name, ensembl_gene_id) #2 columns

#Prepare "feature.data" by Load Ensembl IDs
feature.data = annot[!duplicated(annot$external_gene_name), ]
feature.data <- data.frame(feature.data)
feature.data <- feature.data[, c("ensembl_gene_id", "external_gene_name")]
ensembl.ids = feature.data[,1]
gene.names = feature.data[,2]
genes.after.QC = rownames(Srt_s@assays$RNA@counts)
indices.to.keep = unlist(sapply(genes.after.QC, function(x) {which(gene.names == x)[1]}))
ensembl.ids = ensembl.ids[indices.to.keep]
ensembl.ids <- ensembl.ids[!is.na(ensembl.ids)]##removing NA values
gene.names = gene.names[indices.to.keep]


# Load the MSigDB Hallmark collection using the msigdbr package
H.collection = msigdbr(category="H")###Change here from MSigDB other pathway enrichment databases
?msigdbr()
# Get the entrez gene IDs that are mapped to an Ensembl ID
entrez2ensembl = mappedkeys(org.Hs.egENSEMBL)
# Convert to a list
entrez2ensembl = as.list(org.Hs.egENSEMBL[entrez2ensembl])
# Convert Entrez IDs to Ensembl IDs using the org.Hs.eg.db package
msigdb.entrez.ids = H.collection$entrez_gene
num.ids = length(msigdb.entrez.ids)
msigdb.ensembl.ids = rep(NA, num.ids)
for (i in 1:num.ids) {
  entrez.id = msigdb.entrez.ids[i]
  id.index = which(names(entrez2ensembl) == entrez.id)
  if (length(id.index > 0)) {
    # only use the first mapped ensembl id
    msigdb.ensembl.ids[i] = entrez2ensembl[[id.index]][1]
  }
}
# Save the ensembl IDs in the data frame
H.collection$ensembl_gene = msigdb.ensembl.ids
# Create a gene.set.collection list of Ensembl IDs
gene.set.names = unique(H.collection$gs_name)
num.sets = length(gene.set.names)
gene.set.collection = list()
for (i in 1:num.sets) {
  gene.set.name = gene.set.names[i]
  gene.set.rows = which(H.collection$gs_name == gene.set.name)
  gene.set.ensembl.ids = H.collection$ensembl_gene[gene.set.rows]
  gene.set.collection[[i]] = gene.set.ensembl.ids
}
names(gene.set.collection) = gene.set.names
# Create the collection list required by vamForSeurat()
gene.set.collection = createGeneSetCollection(gene.ids=ensembl.ids, gene.set.collection=gene.set.collection)
length(gene.set.collection)

Srt_s = vamForSeurat(seurat.data=Srt_s,gene.set.collection=gene.set.collection,center=F, gamma=T, sample.cov=F, return.dist=T)
Srt_s@assays$VAMdist[1:5,1:5]
Srt_s@assays$VAMcdf[1:5,1:5]
Srt_s.markers = FindAllMarkers(Srt_s, assay="VAMcdf", only.pos = TRUE, logfc.threshold = 0.01)
Srt_s.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DefaultAssay(object = Srt_s) = "VAMcdf"
top.pathways <- Srt_s.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
p1 <- DoHeatmap(Srt_s, slot="data", features = top.pathways$gene, size=3, label=T) + NoLegend()
#SCF/c-Kit signaling probing##############
Srt <- severe
genes_Seurat <- data.frame(rownames(x = Srt))
colnames(genes_Seurat) <- "external_gene_name"
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
annot<-getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
             filters = "external_gene_name",
             values = genes_Seurat$external_gene_name,
             mart=ensembl)
Biomart_input <- annot %>% distinct(external_gene_name, .keep_all = TRUE) #remove some duplicate genes
non_annot <- as.data.frame(genes_Seurat[!genes_Seurat$external_gene_name %in% annot$external_gene_name,])

options(repr.plot.width=20, repr.plot.height=10)
H.collection = msigdbr(species = "Homo sapiens", category="C2", subcategory = "CP:REACTOME")###to change species, please define the 1st argument of msigdbr "species" to "Mus musculus"
H.collection_b <- H.collection %>% filter(gs_name == "REACTOME_SIGNALING_BY_SCF_KIT") %>% dplyr::select(ensembl_gene)
H.col <- data.frame(H.collection_b)
blymphocyte.gene.ids = as.character(H.col[,1])
blymphocyte.gene.ids
gene.set.name = "REACTOME_SIGNALING_BY_SCF_KIT_gene"
gene.set.id.list = list()
gene.set.id.list[[1]] = blymphocyte.gene.ids
names(gene.set.id.list)[1] = gene.set.name
gene.set.id.list 

# Create the list of gene indices required by vamForSeurat()
(gene.set.collection = createGeneSetCollection(gene.ids=ensembl.ids, gene.set.collection=gene.set.id.list))

gene.indices = gene.set.collection[[1]]
(gene.names = gene.names[gene.indices])

#Run VAM
Srt = vamForSeurat(seurat.data=Srt, gene.set.collection=gene.set.collection, center=F, gamma=T, sample.cov=F, return.dist=T)
# combined@assays$VAMdist[1:5,1:5]

Srt.markers = FindAllMarkers(Srt, assay="VAMcdf", only.pos = TRUE, logfc.threshold = 0.01)
Srt.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

DefaultAssay(object = Srt) = "VAMcdf"
top.pathways <- Srt.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DoHeatmap(Srt, slot="data", features = top.pathways$gene, size=3, label=T)
DimPlot(severe, label = T)
DimPlot(severe, group.by = "Annotation")
###manual list Signaling by SCF-KIT (Homo sapiens)
Srt <- severe
genes_Seurat <- data.frame(rownames(x = Srt))
colnames(genes_Seurat) <- "external_gene_name"
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
annot<-getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
             filters = "external_gene_name",
             values = genes_Seurat$external_gene_name,
             mart=ensembl)
Biomart_input <- annot %>% distinct(external_gene_name, .keep_all = TRUE) #remove some duplicate genes
non_annot <- as.data.frame(genes_Seurat[!genes_Seurat$external_gene_name %in% annot$external_gene_name,])

annot <- Biomart_input %>% dplyr::select(external_gene_name, ensembl_gene_id) #2 columns

#Prepare "feature.data" by Load Ensembl IDs
feature.data = annot[!duplicated(annot$external_gene_name), ]
feature.data <- data.frame(feature.data)
feature.data <- feature.data[, c("ensembl_gene_id", "external_gene_name")]
ensembl.ids = feature.data[,1]
gene.names = feature.data[,2]
genes.after.QC = rownames(Srt@assays$RNA@counts)
indices.to.keep = unlist(sapply(genes.after.QC, function(x) {which(gene.names == x)[1]}))
ensembl.ids = ensembl.ids[indices.to.keep]
ensembl.ids <- ensembl.ids[!is.na(ensembl.ids)]##removing NA values
gene.names = gene.names[indices.to.keep]

gene.set.name = "KIT RECEPTOR SIGNALING WIKI"

data1 <- read.csv("KIT signaling Wikipathways.csv", header=TRUE, stringsAsFactors=FALSE)

data1

SCF.gene.ids = data1$converted_alias #whatever is on the column name of the Ensembl csv file

# Create a collection list for this gene set based on the Ensembl IDs
gene.set.id.list = list()
gene.set.id.list[[1]] = SCF.gene.ids
names(gene.set.id.list)[1] = gene.set.name
gene.set.id.list

# Create the list of gene indices required by vamForSeurat()
(gene.set.collection = createGeneSetCollection(gene.ids=ensembl.ids, gene.set.collection=gene.set.id.list))

severe = vamForSeurat(seurat.data=severe,gene.set.collection=gene.set.collection,center=F, gamma=T, sample.cov=F, return.dist=T)
severe@assays$VAMdist[1,1:10]
severe@assays$VAMcdf[1,1:10]

#Visualize VAM scores
# pdf(paste("COVID19/", "scfLOW_Covid_d7_VAM_SCF_KIT", ".pdf",sep=""))
DefaultAssay(object = severe) = "VAMcdf"
FeaturePlot(severe, reduction="umap", features=gene.set.name)
p2 <- Nebulosa::plot_density(severe, reduction="umap", features = gene.set.name, joint = TRUE, adjust = 1.5, pal = "viridis")
p3 <- DimPlot(severe, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(severe, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'Annotation')
psev <- p1 | (p3/p2)
psev

#pathway enrichment for nonsevere
Srt <- non_sev
genes_Seurat <- data.frame(rownames(x = Srt))
colnames(genes_Seurat) <- "external_gene_name"
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
annot<-getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
             filters = "external_gene_name",
             values = genes_Seurat$external_gene_name,
             mart=ensembl)
Biomart_input <- annot %>% distinct(external_gene_name, .keep_all = TRUE) #remove some duplicate genes
non_annot <- as.data.frame(genes_Seurat[!genes_Seurat$external_gene_name %in% annot$external_gene_name,])

annot <- Biomart_input %>% dplyr::select(external_gene_name, ensembl_gene_id) #2 columns

#Prepare "feature.data" by Load Ensembl IDs
feature.data = annot[!duplicated(annot$external_gene_name), ]
feature.data <- data.frame(feature.data)
feature.data <- feature.data[, c("ensembl_gene_id", "external_gene_name")]
ensembl.ids = feature.data[,1]
gene.names = feature.data[,2]
genes.after.QC = rownames(Srt@assays$RNA@counts)
indices.to.keep = unlist(sapply(genes.after.QC, function(x) {which(gene.names == x)[1]}))
ensembl.ids = ensembl.ids[indices.to.keep]
ensembl.ids <- ensembl.ids[!is.na(ensembl.ids)]##removing NA values
gene.names = gene.names[indices.to.keep]


# Load the MSigDB Hallmark collection using the msigdbr package
H.collection = msigdbr(category="H")###this must be the way to change species
?msigdbr()
# Get the entrez gene IDs that are mapped to an Ensembl ID
entrez2ensembl = mappedkeys(org.Hs.egENSEMBL)
# Convert to a list
entrez2ensembl = as.list(org.Hs.egENSEMBL[entrez2ensembl])
# Convert Entrez IDs to Ensembl IDs using the org.Hs.eg.db package

msigdb.entrez.ids = H.collection$entrez_gene
num.ids = length(msigdb.entrez.ids)
msigdb.ensembl.ids = rep(NA, num.ids)
for (i in 1:num.ids) {
  entrez.id = msigdb.entrez.ids[i]
  id.index = which(names(entrez2ensembl) == entrez.id)
  if (length(id.index > 0)) {
    # only use the first mapped ensembl id
    msigdb.ensembl.ids[i] = entrez2ensembl[[id.index]][1]
  }
}
# Save the ensembl IDs in the data frame
H.collection$ensembl_gene = msigdb.ensembl.ids
# Create a gene.set.collection list of Ensembl IDs
gene.set.names = unique(H.collection$gs_name)
num.sets = length(gene.set.names)
gene.set.collection = list()
for (i in 1:num.sets) {
  gene.set.name = gene.set.names[i]
  gene.set.rows = which(H.collection$gs_name == gene.set.name)
  gene.set.ensembl.ids = H.collection$ensembl_gene[gene.set.rows]
  gene.set.collection[[i]] = gene.set.ensembl.ids
}
names(gene.set.collection) = gene.set.names
# Create the collection list required by vamForSeurat()
gene.set.collection = createGeneSetCollection(gene.ids=ensembl.ids, gene.set.collection=gene.set.collection)
length(gene.set.collection)

Srt = vamForSeurat(seurat.data=Srt,gene.set.collection=gene.set.collection,center=F, gamma=T, sample.cov=F, return.dist=T)
Srt@assays$VAMdist[1:5,1:5]
Srt@assays$VAMcdf[1:5,1:5]
Srt.markers = FindAllMarkers(Srt, assay="VAMcdf", only.pos = TRUE, logfc.threshold = 0.01)
Srt.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DefaultAssay(object = Srt) = "VAMcdf"
top.pathways <- Srt.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
p4 <- DoHeatmap(Srt, slot="data", features = top.pathways$gene, size=3, label=T) + NoLegend()
#SCF/c-Kit signaling probing##############
Srt <- non_sev
genes_Seurat <- data.frame(rownames(x = Srt))
colnames(genes_Seurat) <- "external_gene_name"
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
annot<-getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
             filters = "external_gene_name",
             values = genes_Seurat$external_gene_name,
             mart=ensembl)
Biomart_input <- annot %>% distinct(external_gene_name, .keep_all = TRUE) #remove some duplicate genes
non_annot <- as.data.frame(genes_Seurat[!genes_Seurat$external_gene_name %in% annot$external_gene_name,])

options(repr.plot.width=20, repr.plot.height=10)
H.collection = msigdbr(species = "Homo sapiens", category="C2", subcategory = "CP:REACTOME")###to change species, please define the 1st argument of msigdbr "species" to "Mus musculus"
H.collection_b <- H.collection %>% filter(gs_name == "REACTOME_SIGNALING_BY_SCF_KIT") %>% dplyr::select(ensembl_gene)
H.col <- data.frame(H.collection_b)
blymphocyte.gene.ids = as.character(H.col[,1])
blymphocyte.gene.ids
gene.set.name = "REACTOME_SIGNALING_BY_SCF_KIT_gene"
gene.set.id.list = list()
gene.set.id.list[[1]] = blymphocyte.gene.ids
names(gene.set.id.list)[1] = gene.set.name
gene.set.id.list 

# Create the list of gene indices required by vamForSeurat()
(gene.set.collection = createGeneSetCollection(gene.ids=ensembl.ids, gene.set.collection=gene.set.id.list))

gene.indices = gene.set.collection[[1]]
(gene.names = gene.names[gene.indices])

#Run VAM
Srt = vamForSeurat(seurat.data=Srt, gene.set.collection=gene.set.collection, center=F, gamma=T, sample.cov=F, return.dist=T)
# combined@assays$VAMdist[1:5,1:5]

Srt.markers = FindAllMarkers(Srt, assay="VAMcdf", only.pos = TRUE, logfc.threshold = 0.01)
Srt.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

DefaultAssay(object = Srt) = "VAMcdf"
top.pathways <- Srt.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DoHeatmap(Srt, slot="data", features = top.pathways$gene, size=3, label=T)
DimPlot(non_sev, label = T)
DimPlot(non_sev, group.by = "Annotation")
###manual list Signaling by SCF-KIT (Homo sapiens)
Srt <- non_sev
genes_Seurat <- data.frame(rownames(x = Srt))
colnames(genes_Seurat) <- "external_gene_name"
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
annot<-getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
             filters = "external_gene_name",
             values = genes_Seurat$external_gene_name,
             mart=ensembl)
Biomart_input <- annot %>% distinct(external_gene_name, .keep_all = TRUE) #remove some duplicate genes
non_annot <- as.data.frame(genes_Seurat[!genes_Seurat$external_gene_name %in% annot$external_gene_name,])

annot <- Biomart_input %>% dplyr::select(external_gene_name, ensembl_gene_id) #2 columns

#Prepare "feature.data" by Load Ensembl IDs
feature.data = annot[!duplicated(annot$external_gene_name), ]
feature.data <- data.frame(feature.data)
feature.data <- feature.data[, c("ensembl_gene_id", "external_gene_name")]
ensembl.ids = feature.data[,1]
gene.names = feature.data[,2]
genes.after.QC = rownames(Srt@assays$RNA@counts)
indices.to.keep = unlist(sapply(genes.after.QC, function(x) {which(gene.names == x)[1]}))
ensembl.ids = ensembl.ids[indices.to.keep]
ensembl.ids <- ensembl.ids[!is.na(ensembl.ids)]##removing NA values
gene.names = gene.names[indices.to.keep]

gene.set.name = "KIT RECEPTOR SIGNALING WIKI"

data1 <- read.csv("KIT signaling Wikipathways.csv", header=TRUE, stringsAsFactors=FALSE)

data1

SCF.gene.ids = data1$converted_alias #whatever is on the column name of the Ensembl csv file

# Create a collection list for this gene set based on the Ensembl IDs
gene.set.id.list = list()
gene.set.id.list[[1]] = SCF.gene.ids
names(gene.set.id.list)[1] = gene.set.name
gene.set.id.list

# Create the list of gene indices required by vamForSeurat()
(gene.set.collection = createGeneSetCollection(gene.ids=ensembl.ids, gene.set.collection=gene.set.id.list))

non_sev = vamForSeurat(seurat.data=non_sev,gene.set.collection=gene.set.collection,center=F, gamma=T, sample.cov=F, return.dist=T)
non_sev@assays$VAMdist[1,1:10]
non_sev@assays$VAMcdf[1,1:10]

#Visualize VAM scores
# pdf(paste("COVID19/", "scfLOW_Covid_d7_VAM_SCF_KIT", ".pdf",sep=""))
DefaultAssay(object = non_sev) = "VAMcdf"
FeaturePlot(non_sev, reduction="umap", features=gene.set.name)
p5 <- Nebulosa::plot_density(non_sev, reduction="umap", features = gene.set.name, joint = TRUE, adjust = 1.5, pal = "viridis")
p6 <- DimPlot(non_sev, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(non_sev, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'Annotation')


#Plotting
patchnosev <- p4 | (p6/p5)
