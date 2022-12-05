# NetBID-in-COVID-19-plasma-proteomics

A) Create an appropriate format for the dataset
1.MGH_data_Preprocessing.R

B) Filtering and data construction for SJARACNe algorithm.
2.Network_Reconstruction.R

C) Perform SJARACNe algorithm, extracting non-linear correlations between proteins
3.SJARACNe.sh

D) Perform Bayesian Inference algorithm for feature selection
4.Driver_Inference.R

E) Load dataset into Cytoscape for better visualization of results for annotation and Enrichment analysis.
5.Cytoscape_load.R

F) Advanced analysis for plotting and interpretation
6.Advanced_Analysis.R

For the biological informed neural network (b-iNN) is followed the <a href="https://github.com/DataX-JieHao/PASNet#interpretable-neural-network-on-the-biological-pathway-level">PasNet</a> architecture. The changes are committed to the b-iNN are described in the following steps: </p>
1. Pathway Layer
## <p> At this step it is created the REACTOM_pathway_mask.csv and GO_pathway_mask.csv documents for Reactom and GO pathway libraries, respectively. These documents are used for the two different b-iNN that are created for this study.</p>
2. Explainable AI feature
<p> Next, the Train.py file of <a href="https://github.com/DataX-JieHao/PASNet#interpretable-neural-network-on-the-biological-pathway-level">PasNet</a> modified as follows:</p>

3. Receiver Operating Characteristic plot
<p> The EvalFunc.py file of <a href="https://github.com/DataX-JieHao/PASNet#interpretable-neural-network-on-the-biological-pathway-level">PasNet</a> modified as follows:</p>

4. Train and Test data
<p>For the training and test of the model has used the MGH dataset. The MGH dataset is normalized with mean = 0 and standard deviation = 1. The produced datasets are the following: </p>
 - Train: std_train_0_0.csv
 - Test: std_test_0_0.csv

5. Execution of b-iNN
<p> For the execution of the experiment first it is executed the Run_EmpiricalSearch.py. After optimal L2 and LR attributes are detected, it is executed the Run.py.</p>

## scRNAseq from MGH
The .h5ad scRNA-seq data for the analysis are accessible in the following link: https://www.covid19cellatlas.org/index.patient.html (“00803_MGH_Broad00803_MGH_Broad”). The scRNA-seq analysis can then be recreated using the “scDIOR Python MGH.ipynB” script in Python to make an adata object and then the “MGHscRNAseq.R” script in R to convert the adata object into Seurat object and then perform downstream analysis.
## drug repurposing
D-SJARACNe was assembled using the Day 7 SJARACNe network (“Cytoscape_network_D7.cys”) with data from Zhou et al. (https://doi.org/10.1038/s41587-022-01474-0). The drug-drug interaction network was assembled by using the drugs from D-SJARACNe as input in the STITCH database. Visualizations were performed in the Gephi v0.9 platform.
