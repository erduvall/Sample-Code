# Eliza Duvall
# January 2024

## Pathway benchmarking project requires perturbation of networks to compare analysis
## This way we know the ground truth

#### Notes
## Expression data (pasc.h5seurat) is large and must be opened on QUEST (northwetsern's super-computer) :) 
## Enrichment browser uses the gene sets then reconstructs the networks before applying to various methods
## We want to perturb genes with high centrality and their neighbors 

#### Steps
## get list of gene sets and decide "middle size" boundaries
## call in expression data (counts)
## look at average expression across all genes in a pathway to filter pathways again
## look at distribution of expression in each gene per pathway
## visualize final pathway candidates


###########
# in quest
###########
# > module load R/4.2.3
# > module load hdf5/1.8.19-serial
# > srun -A b1167 --partition=b1167 -N 1 -n 1 --mem=64G --time=4:00:00 --pty bash -i 
# > R


###########
# in R (on quest)
###########
library(EnrichmentBrowser)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(Seurat)
library(SeuratDisk)

## create list of all gene sets
gs <- EnrichmentBrowser::getGenesets(org = "hsa", db = "kegg", return.type = c("list", "GeneSetCollection"))
# hist(lengths(gs), breaks = 200)
## filter gene sets to have between 50-150 genes in them 
gs <- gs[which(lengths(gs) < 150 & lengths(gs) > 50)]
# hist(lengths(gs), breaks = 50)

# convert gene names to align with the expression data (counts)
# gsIndex is a numerical value to grab the correct index of the list gs
gs_symb <- lapply(1:length(gs), function(gsIndex){
  mapIds(org.Hs.eg.db, keys = gs[[gsIndex]], keytype = "ENTREZID", column = "SYMBOL")
})
names(gs_symb) <- names(gs)
# get single list of all the gene symbols used across all gene sets
allgenes <- unlist(gs_symb) %>% unique()

# call in expression data 
counts <- LoadH5Seurat("pasc.h5seurat")
# subset to the genes in gs_symb
counts_subset <- counts[which(row.names(counts) %in% allgenes),]

# subset for the healthy people!
healthy <- which(counts@meta.data[,c("Status")] == "Healthy")
counts_subset <- counts_subset[,healthy]

# subset for the cell type with the most amount of cells
counts_subset@meta.data[,c("cell_type")] %>% table()
tram1Cells <- which(counts_subset@meta.data[,c("cell_type")] == "TRAM-1")
counts_subset <- counts_subset[,tram1Cells]
dim(counts_subset)
## 4825 6120

# save this expression data as a subset seurat object
healthyTRAM1_counts <- counts_subset
SaveH5Seurat(healthyTRAM1_counts, name = "preprocessing/healthyTRAM1_counts.h5Seurat")

# function that returns the average gene expression for each gene in a pathway
getMEDGeneExpr <- function(gsIndex){
  subCounts <- counts_subset[which(row.names(counts_subset) %in% gs_symb[[gsIndex]]),]
  subCounts <- GetAssayData(subCounts, assay = "RNA", slot = "counts")
  geneMedians <- apply(subCounts, 1, median, na.rm=T)
  return(geneMedians)
}
#### NOTE: I originally used mean but there were a few genes highly expressed and
#### most genes were not expressed.... so I switched to MEDIAN <3 

# apply function to return a list of the average gene expression
geneMedians <- lapply(1:length(gs_symb), function(gsIndex){
  getMEDGeneExpr(gsIndex)
}) 
names(geneMedians) <- names(gs_symb)

## save the geneMedians as an RDA file! 
save(geneMedians, file = "preprocessing/gsExprMedians_perGene.RData")
save(gs_symb, file = "preprocessing/gs_bySymb.RData")


###########
# local
###########

# setwd("MYPATH")
# pull those saved files to local computer and load :) 
load("gsExprMeans_perGene.RData")
load("gs_bySymb.RData")

# get means for each gene set
gsMedians <- sapply(1:length(geneMedians), function(gsIndex){
  mean(geneMedians[[gsIndex]])
})

hist(gsMedians, breaks = 25)
potenPways <- gs_symb[which(gsMedians > 0)]

# look at expression distribution for potential pathways
res <- geneMedians[(names(potenPways)[1])]
hist(unlist(res), breaks = 10)

#### NOTE: I explored means and medians
#### PATHWAY TO PERTURB (potentially): HSA00190 Oxidative Phosphorylation








