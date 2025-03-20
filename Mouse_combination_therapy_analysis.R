
library(Seurat)
library(patchwork)
library(dittoSeq)
library(ggplot2)
library(Seurat)
library(dplyr)
library(sctransform)
library(cowplot)
library(BRETIGEA)
library(knitr)
library(stringr)
library(patchwork)

## Input snRNA-seq data and add metadata
# input males 
YL01_M100 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_M100/")
YL01_M100 <- CreateSeuratObject(counts = YL01_M100, project = "YL01_M100")
YL01_M220 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_M220/")
YL01_M220 <- CreateSeuratObject(counts = YL01_M220, project = "YL01_M220")
YL01_M37 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_M37/")
YL01_M37 <- CreateSeuratObject(counts = YL01_M37, project = "YL01_M37")
YL01_M265 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_M265/")
YL01_M265 <- CreateSeuratObject(counts = YL01_M265, project = "YL01_M265")
male <- merge(YL01_M100, y = c(YL01_M220, YL01_M37, YL01_M265), add.cell.ids = c("M100", "M220", "M37", "M265"), project = "male")
male@meta.data$group <- "vehicle"
YL01_M413 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_M413/")
YL01_M413 <- CreateSeuratObject(counts = YL01_M413, project = "YL01_M413")
YL01_M6 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_M6/")
YL01_M6 <- CreateSeuratObject(counts = YL01_M6, project = "YL01_M6")
YL01_M254 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_M254/")
YL01_M254 <- CreateSeuratObject(counts = YL01_M254, project = "YL01_M254")
YL01_M48 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_M48/")
YL01_M48 <- CreateSeuratObject(counts = YL01_M48, project = "YL01_M48")
male1 <- merge(YL01_M413, y = c(YL01_M6, YL01_M254, YL01_M48), add.cell.ids = c("M413", "M6", "M254", "M48"), project = "male")
male1@meta.data$group <- "treatment"
male <- merge(male, y = c(male1), project = "male")
male@meta.data$sex <- "male"
rm(list=setdiff(ls(), "male"))

# input females  
YL01_F199 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_F199/")
YL01_F199 <- CreateSeuratObject(counts = YL01_F199, project = "YL01_F199")
YL01_F145 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_F145/")
YL01_F145 <- CreateSeuratObject(counts = YL01_F145, project = "YL01_F145")
YL01_F352 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_F352/")
YL01_F352 <- CreateSeuratObject(counts = YL01_F352, project = "YL01_F352")
YL01_F161 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_F161/")
YL01_F161 <- CreateSeuratObject(counts = YL01_F161, project = "YL01_F161")
female <- merge(YL01_F199, y = c(YL01_F145, YL01_F352, YL01_F161), add.cell.ids = c("F199", "F145", "F352", "F161"), project = "female")
female@meta.data$group <- "vehicle"
YL01_F308 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_F308/")
YL01_F308 <- CreateSeuratObject(counts = YL01_F308, project = "YL01_F308")
YL01_F312 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_F312/")
YL01_F312 <- CreateSeuratObject(counts = YL01_F312, project = "YL01_F312")
YL01_F361 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_F361/")
YL01_F361 <- CreateSeuratObject(counts = YL01_F361, project = "YL01_F361")
YL01_F44 <- Read10X(data.dir = "~/data/DRsinglecell/oct3/YL01_F44/")
YL01_F44 <- CreateSeuratObject(counts = YL01_F44, project = "YL01_F44")
female1 <- merge(YL01_F308, y = c(YL01_F312, YL01_F361, YL01_F44), add.cell.ids = c("F308", "F312", "F361", "F44"), project = "female")
female1@meta.data$group <- "treatment"
female <- merge(female, y = c(female1), project = "female")
female@meta.data$sex <- "female"

# merge data into one seurat object
DRsn <- merge(male, y = c(female), project = "DRsn-seq")



## Standard processing: remove poor quality cells and genes expressed in less than 10 cells
DRsn <- PercentageFeatureSet(DRsn, pattern = "^mt-", col.name = "percent.mt")
# Add number of genes per UMI for each cell to metadata
DRsn$log10GenesPerUMI <- log10(DRsn$nFeature_RNA) / log10(DRsn$nCount_RNA)
VlnPlot(DRsn1, features = c("log10GenesPerUMI", "percent.mt"), ncol = 3, group.by = "orig.ident")
VlnPlot(DRsn1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3,  group.by = "orig.ident")
# subset by RNA count and gene feature counts
DRsn1 <- subset(DRsn, subset = nFeature_RNA > 200 & log10GenesPerUMI > 0.85 & nCount_RNA > 500)

#Gene-level filtering
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = DRsn1, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
DRsn1 <- CreateSeuratObject(filtered_counts, meta.data = DRsn1@meta.data)
rm(counts, filtered_counts, nonzero, keep_genes)


## Standard processing: normalization and scale.data
DRsn2 <- SCTransform(DRsn1, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:10, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:10, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)




## Cell type annotations
#input cell type marker genes - including hippocampal neuronal subtype markers
Dentate.Gyrus<- list(c("Ntng1",  "Plk5",  "Tanc1", "Ahcyl2", "Slc4a4", "Rfx3"))
Oligodendrocytes<- list(c("Mbp",  "Plp1", "Mag"))
CA1.Pyramidal.Neurons<- list(c("Ccdc88c", "Galntl6", "Man1a"))
CA3.Pyramidal.Neurons<- list(c("Trhde", "Rnf182", "Nectin3"))
Interneurons<- list(c("Nxph1", "Grik1", "Sst", "Cnr1", "Vip", "Npas1"))
Subiculum<- list(c("Tshz2", "Sgcz", "Meis2"))
Oligodendrocyte.Precursor<- list(c("Vcan", "Pdgfra", "Arhgap31"))
Microglia<- list(c("Arhgap45", "Ly86",  "Cd37"))
Astrocytes<-list(c("Slc1a2",  "Atp1a2",  "Phkg1",   "Slc7a11", "Atp1a2",  "Ranbp3l"))
Fibroblast.like<- list(c("Slc38a2", "Mgp", "Sned1"))
Choroid.Plexus<- list(c("Col9a3", "Ttr", "Htr2c"))
Ex.Neurons <- list(c("Slc17a6", "Slc17a7", "Neurod6", "Nrgn", "Camk2a", "Satb2", "Col5a1", "Sdk2", "Nefm"))
In.Neurons <- list(c("Slc32a1", "Gad1", "Gad2", "Tac1", "Penk", "Sst", "Npy", "Mybpc1", "Pvalb", "Gabbr2"))
# for general cell types 
#Celltype <- c("Ex.Neurons", "In.Neurons", "Astrocytes", "Microglia", "Fibroblast.like", "Oligodendrocytes", "Oligodendrocyte.Precursor", "Choroid.Plexus")


DefaultAssay(DRsn2) <- "SCT"
dataset2 <- DRsn2
# for neuronal subtype identification in hippocampal samples
Celltype <- c( "Dentate.Gyrus", "CA1.Pyramidal.Neurons", "CA3.Pyramidal.Neurons", "Interneurons", "Subiculum", "Astrocytes", "Microglia", "Fibroblast.like", "Oligodendrocytes", "Oligodendrocyte.Precursor", "Choroid.Plexus")
library(stringr)
for(key in Celltype){
  dataset2<-AddModuleScore(object=dataset2, features = get(key), pool = NULL,  nbin = 24,  ctrl = 30,  k = FALSE,  assay = NULL, name = key,  seed = 1,  search = FALSE)
}
#Make cell type column
new_meta<- DRsn2@meta.data
colnames(new_meta)
new_meta[is.na(new_meta)] <- 0
new_meta$cell_type<- names(new_meta[, 12:22])[apply(new_meta[, 12:22],1,which.max)]
table(new_meta$cell_type)
#Max and second max
new_meta$x1<- apply(new_meta[, 12:22], 1, max)
new_meta$x2<- apply(new_meta[, 12:22], 1, function(x) x[order(x)[10]]) #in decreasing order
new_meta$x1_x2<- (new_meta$x1 - new_meta$x2)/new_meta$x1
new_meta$x1_x2<- ifelse(new_meta$x1_x2 < 0.1, 1,0) #if the difference between top two scores is less than 10%, 1, if not 0
table(new_meta$x1_x2) #11385 0, 319 1
new_meta$x3<- apply(new_meta[, 12:16], 1, max)
new_meta$x4<- apply(new_meta[, 12:16], 1, function(x) x[order(x)[4]]) #in decreasing order
new_meta$x5<- (new_meta$x4 - new_meta$x3)/new_meta$x4
new_meta$x5<- ifelse(new_meta$x5 < 0.1, 1,0) #if the difference between top two scores is less than 10%, 1, if not 0
#Make second column including hybrid
new_meta$cell_type2<- new_meta$cell_type
new_meta$cell_type2<- ifelse(new_meta$x1_x2 == 1, "hybrid", new_meta$cell_type2 )
new_meta$cell_type2<- ifelse(new_meta$x1_x2 == 1 & new_meta$x5 == 1, "neuronal.hybrid", new_meta$cell_type2 )
new_meta$cell_type2<- ifelse(new_meta$x1 < 0 , "negative", new_meta$cell_type2 )
table(new_meta$cell_type2)
new_meta$cell_type2<- gsub("1", "", new_meta$cell_type2)
#add cell type annotation as metadata
DRsn2 <-AddMetaData(DRsn2, new_meta$cell_type2, col.name = 'neuron_sub_type')
DRsn2 <-AddMetaData(DRsn2, new_meta$cell_type, col.name = 'cell_type')



## Identify differential expressed genes in combination-treated versus vehicle-treated


Celltype <- c("Ex.Neurons", "In.Neurons", "Astrocytes", "Microglia", "Oligodendrocytes", "Oligodendrocyte.Precursor")
library(stringr)
DEG.list.FDR2 <- NULL
DEG.list.FDR<- NULL
for (key in Celltype)  {
  Idents(DRsn2) <- "cell_type_ident"
  subcelltype <- subset(DRsn2, idents= key)
  Idents(subcelltype) <- "Integrated_case_ctl"
  DEG.list <- FindMarkers(subcelltype, ident.1= "treatment", ident.2 ="vehicle", verbose = FALSE, test.use = "MAST", logfc.threshold = 0.1)
  DEG.list$gene <- row.names(DEG.list)
  DEG.list.FDR <- subset(DEG.list, DEG.list$p_val_adj< 0.05)
  DEG.list.FDR$dir<- ifelse(DEG.list.FDR$avg_log2 < 0, "neg","pos")
  DEG.list.FDR$celltype<- key
  DEG.list.FDR2<- rbind(DEG.list.FDR2, DEG.list.FDR)
}









### Mapping between genes names in human and mouse
#input human AD signatures
mzl<- read.csv(("~/Downloads/10.27DEGs/2021.10.29.Integrated.0.1.diag.DEG.MAST.csv"))
#input human and mouse gene reference
mart<- read.csv(("~/Downloads/mart_export.csv"))
mart<-subset(mart,!duplicated(mart$Gene.name))

mzl$Gene.name <- mzl$gene
library(dplyr, lib.loc = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library")
mzl <- inner_join(mzl, mart, by = 'Gene.name')


