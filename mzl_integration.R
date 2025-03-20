# Preparing AD snRNA-seq datasets for integration

library(DittoSeq)
library(ggplot2)
library(Seurat) 
library(dplyr)
library(sctransform)
library(cowplot)
library(BRETIGEA)
library(knitr)
library(stringr)
library(patchwork)

### Dataset 1 - Lau
## Input Lau et al dataset from GEO GSE157827_RAW
# create seurat object for Lau et al dataset
AD1 <- Read10X(data.dir = "../Lau/GSE157827_RAW/AD1/")
AD1 <- CreateSeuratObject(counts = AD1, project = "AD1")
AD2 <- Read10X(data.dir = "../Lau/GSE157827_RAW/AD2/")
AD2 <- CreateSeuratObject(counts = AD2, project = "AD2")
AD4 <- Read10X(data.dir = "../Lau/GSE157827_RAW/AD4/")
AD4 <- CreateSeuratObject(counts = AD4, project = "AD4")
AD5 <- Read10X(data.dir = "../Lau/GSE157827_RAW/AD5/")
AD5 <- CreateSeuratObject(counts = AD5, project = "AD5")
AD8 <- Read10X(data.dir = "../Lau/GSE157827_RAW/AD8/")
AD8 <- CreateSeuratObject(counts = AD8, project = "AD8")
AD9 <- Read10X(data.dir = "../Lau/GSE157827_RAW/AD9/")
AD9 <- CreateSeuratObject(counts = AD9, project = "AD9")
AD10 <- Read10X(data.dir = "../Lau/GSE157827_RAW/AD10/")
AD10 <- CreateSeuratObject(counts = AD10, project = "AD10")
AD13 <- Read10X(data.dir = "../Lau/GSE157827_RAW/AD13/")
AD13 <- CreateSeuratObject(counts = AD13, project = "AD13")
AD19 <- Read10X(data.dir = "../Lau/GSE157827_RAW/AD19/")
AD19 <- CreateSeuratObject(counts = AD19, project = "AD19")
AD20 <- Read10X(data.dir = "../Lau/GSE157827_RAW/AD20/")
AD20 <- CreateSeuratObject(counts = AD20, project = "AD20")
AD21 <- Read10X(data.dir = "../Lau/GSE157827_RAW/AD21/")
AD21 <- CreateSeuratObject(counts = AD21, project = "AD21")

C3 <- Read10X(data.dir = "../Lau/GSE157827_RAW/C3/")
C3 <- CreateSeuratObject(counts = C3, project = "C3")
C7 <- Read10X(data.dir = "../Lau/GSE157827_RAW/C7/")
C7 <- CreateSeuratObject(counts = C7, project = "C7")
C11 <- Read10X(data.dir = "../Lau/GSE157827_RAW/C11/")
C11 <- CreateSeuratObject(counts = C11, project = "C11")
C12 <- Read10X(data.dir = "../Lau/GSE157827_RAW/C12/")
C12 <- CreateSeuratObject(counts = C12, project = "C12")
C14 <- Read10X(data.dir = "../Lau/GSE157827_RAW/C14/")
C14 <- CreateSeuratObject(counts = C14, project = "C14")
C15 <- Read10X(data.dir = "../Lau/GSE157827_RAW/C15/")
C15 <- CreateSeuratObject(counts = C15, project = "C15")
C16 <- Read10X(data.dir = "../Lau/GSE157827_RAW/C16/")
C16 <- CreateSeuratObject(counts = C16, project = "C16")
C17 <- Read10X(data.dir = "../Lau/GSE157827_RAW/C17/")
C17 <- CreateSeuratObject(counts = C17, project = "C17")
C18 <- Read10X(data.dir = "../Lau/GSE157827_RAW/C18/")
C18 <- CreateSeuratObject(counts = C18, project = "C18")

# add metadata about case/ctl from original study- "diagnosis" lables original AD classification by Lau et al.
AD.final <- merge(AD1, y = c(AD2, AD4, AD5, AD8, AD9, AD10, AD13, AD19, AD20, AD21), add.cell.ids = c("AD1", "AD2", "AD4", "AD5", "AD8", "AD9", "AD10", "AD13", "AD19", "AD20", "AD21"), project = "AD")
AD.final@meta.data$diagnosis <- "AD"
ctl.final<- merge(C3, y = c(C7, C11, C12, C14, C15, C16, C17, C18), add.cell.ids = c("C3", "C7", "C11", "C12", "C14", "C15", "C16", "C17", "C18"), project = "ctl")
ctl.final@meta.data$diagnosis <- "ctl"
lau<-merge(AD.final, y=ctl.final, add.cell.ids = c("AD", "ctl"), project = "Lau")

# add other metadata
meta.lau <- read.csv("~/data/Lau/Lau_metadata.csv")
meta.lau$sex <- factor(meta.lau$SEX, levels = c("F","M"), labels = c("female","male"))
meta.lau$apoE <- factor(meta.lau$APOE, levels = c("E2/E3", "E2/E4", "E3/E3", "E3/E4", "E4/E4"), labels = c("23", "24", "33", "34", "44"))
meta.lau$cerad <- factor(meta.lau$AGE.RELATED.PLAQUE.SCORE, levels = c("0","A", "B", "C"), labels = c("4","3", "2", "1"))
individualID<- as.data.frame(lau$orig.ident, colnames="ID")
colnames(individualID) <- "ID"
meta.lau<- inner_join(individualID, meta.lau,by = "ID", copy = FALSE)
lau<- AddMetaData(lau, meta.lau$sex, col.name = 'sex')
lau<- AddMetaData(lau, meta.lau$apoE, col.name = 'apoE')
lau<- AddMetaData(lau, meta.lau$Braak.tangle.stage, col.name = 'braak')
lau<- AddMetaData(lau, meta.lau$cerad, col.name = 'cerad')


### Dataset 2 - Mathys
## Input Mathys el al 2019 dataset from AD knowledge Portal - ROSMAP
# prepare metadata
meta.mathys<- read.table("~/2021June/mathy/filtered_column_metadata.txt",header= TRUE, sep = "\t")
M_meta<- read.csv("~/2021June/mathy/raw /Mathy_pathology.csv" )
meta.mathys<- inner_join(meta.mathys,M_meta,by = "projid", copy = FALSE)
meta.mathys$sex <- factor(meta.mathys$msex, levels = c(0,1), labels = c("female","male"))
meta.mathys$orig.ident<- paste("M", meta.mathys$Subject, sep = "-")
rosmap <-read.csv("~/2021June/Zhou/ROSMAP_clinical.csv")
meta.mathys<- inner_join(meta.mathys,rosmap,by = "projid", copy = FALSE)

# Prepare count data and create seurat object for Mathys dataset
genes<- read.table("~/data/mathy/raw /features.tsv.gz",header= FALSE, sep = "\t")
barcodes<- as.data.frame(math_meta$TAG)
write.table(barcodes, file='barcodes.tsv', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
features<- read.table("filtered_gene_row_names.txt",header= FALSE, sep = "\t")
write.table(genes, file='features.tsv', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
math <- Read10X(data.dir = "~/2021June/mathy/raw /",  gene.column=1)
math <- CreateSeuratObject(counts = math, project = "Mathys")
# Add metadata
math <- AddMetaData(math, meta.mathys$orig.ident, col.name = 'orig.ident')
math <- AddMetaData(math, meta.mathys$broad.cell.type, col.name = 'M_cell_type')
math <- AddMetaData(math, meta.mathys$sex, col.name = 'sex')
# "diagnosis" labels original case/control classification by Mathys et al.
math<- AddMetaData(math, meta.mathys$pathologic.diagnosis.of.AD, col.name = 'diagnosis')
math<- AddMetaData(math, meta.mathys$apoe_genotype, col.name = 'apoE')
math<- AddMetaData(math, meta.mathys$braaksc, col.name = 'braak')
math<- AddMetaData(math, meta.mathys$ceradsc, col.name = 'cerad')




### Dataset 3 - Zhou 
## Input Zhou el al 2020 dataset from AD knowledge Portal - ROSMAP
# input metadata
zhou.meta<- read.csv("~/data/Zhou/snRNAseqAD_TREM2_biospecimen_metadata.csv")
rosmap <-read.csv("~/data/Zhou/ROSMAP_clinical.csv")
zhou.meta<- inner_join(zhou.meta,rosmap,by = "individualID", copy = FALSE)

# create Seurat object
data_dir<-'data/Zhou/gene expression/'
AD1 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/AD1/")
AD1 <- CreateSeuratObject(counts = AD1, project = "AD1")
AD2 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/AD2/")
AD2 <- CreateSeuratObject(counts = AD2, project = "AD2")
AD3 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/AD3/")
AD3 <- CreateSeuratObject(counts = AD3, project = "AD3")
AD5 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/AD5/")
AD5 <- CreateSeuratObject(counts = AD5, project = "AD5")
AD7 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/AD7/")
AD7 <- CreateSeuratObject(counts = AD7, project = "AD7")
AD8 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/AD8/")
AD8 <- CreateSeuratObject(counts = AD8, project = "AD8")
AD9 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/AD9/")
AD9 <- CreateSeuratObject(counts = AD9, project = "AD9")
AD10 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/AD10/")
AD10 <- CreateSeuratObject(counts = AD10, project = "AD10")
AD11 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/AD11/")
AD11 <- CreateSeuratObject(counts = AD11, project = "AD11")
AD12 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/AD12/")
AD12 <- CreateSeuratObject(counts = AD12, project = "AD12")
AD13 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/AD13/")
AD13 <- CreateSeuratObject(counts = AD13, project = "AD13")

C1 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/C1/")
C1 <- CreateSeuratObject(counts = C1, project = "C1")
C2 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/C2/")
C2 <- CreateSeuratObject(counts = C2, project = "C2")
C3 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/C3/")
C3 <- CreateSeuratObject(counts = C3, project = "C3")
C4 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/C4/")
C4 <- CreateSeuratObject(counts = C4, project = "C4")
C5 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/C5/")
C5 <- CreateSeuratObject(counts = C5, project = "C5")
C6 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/C6/")
C6 <- CreateSeuratObject(counts = C6, project = "C6")
C7 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/C7/")
C7 <- CreateSeuratObject(counts = C7, project = "C7")
C8 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/C8/")
C8 <- CreateSeuratObject(counts = C8, project = "C8")
C9 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/C9/")
C9 <- CreateSeuratObject(counts = C9, project = "C9")
C11 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/C11/")
C11 <- CreateSeuratObject(counts = C11, project = "C11")
C12 <- Read10X(data.dir = "../ubuntu/data/Zhou/gene expression/C12/")
C12 <- CreateSeuratObject(counts = C12, project = "C12")

# "diagnosis" labels original case/control classification by Zhou et al
AD.final <- merge(AD1, y = c(AD2, AD3, AD5, AD7, AD8, AD9, AD10, AD11, AD12, AD13), add.cell.ids = c("AD1", "AD2", "AD3", "AD5", "AD7", "AD8", "AD9", "AD10", "AD11", "AD12", "AD13"), project = "AD")
AD.final@meta.data$diagnosis <- "AD"
ctl.final<- merge(C1, y = c(C2, C3, C4, C5, C6, C7, C8, C9, C11, C12), add.cell.ids = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C11", "C12"), project = "ctl")
ctl.final@meta.data$diagnosis <- "ctl"
Zhou<-merge(AD.final, y=ctl.final, add.cell.ids = c("AD", "ctl"), project = "Zhou")

# add metadata to object
individualID<- as.data.frame(Zhou$orig.ident, colnames="specimenID")
colnames(individualID) <-c("specimenID")
zhou.meta<- inner_join(individualID, zhou.meta, by = "specimenID")
zhou.meta$orig.ident<- paste("Z", zhou.meta$specimenID, sep = "-")
zhou<- AddMetaData(zhou, zhou.meta$orig.ident, col.name = 'orig.ident')
sex <- factor(zhou.meta$msex, levels = c(0,1), labels = c("female","male"))
zhou<- AddMetaData(zhou, sex, col.name = 'sex')
apoe<- zhou.meta$apoe_genotype
zhou<- AddMetaData(zhou, apoe, col.name = 'apoE')
zhou<- AddMetaData(zhou, zhou.meta$braaksc, col.name = 'braak')
zhou<- AddMetaData(zhou, zhou.meta$ceradsc, col.name = 'cerad')





### Integration by canonical correlation analysis from Seurat package

## Prepare datasets for integration
zhou@meta.data$dataset <-"Zhou"
math@meta.data$dataset <-"Mathys"
lau@meta.data$dataset <-"Lau"
# merge data
mzl<-merge(zhou, y= c(math, lau), add.cell.ids = c("zhou","math", "lau"), project = "dataset")

#An object of class Seurat 
# 45119 features across 328950 samples within 1 assay 
#Active assay: RNA (45119 features, 0 variable features)

# remove poor quality cells and low expressed features
mzl$log10GenesPerUMI <- log10(mzl$nFeature_RNA) / log10(mzl$nCount_RNA)
# Compute percent mito ratio
mzl$mitoRatio <- PercentageFeatureSet(object = mzl, pattern = "^MT-")
mzl$mitoRatio <- mzl@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- mzl@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

mzl<- subset(x = mzl, subset= (nCount_RNA >= 500) & 
               (nFeature_RNA >= 250) & 
               (log10GenesPerUMI > 0.90) & 
               (mitoRatio < 0.1))


counts <- GetAssayData(object = mzl, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 100

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
mzl <- CreateSeuratObject(filtered_counts, meta.data = mzl@meta.data)

## run integration algorithm - CCA
mzl <- SplitObject(mzl, split.by = "dataset")
mzl <- lapply(X = mzl, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})
features <- SelectIntegrationFeatures(object.list = mzl)

anchors <- FindIntegrationAnchors(object.list = mzl, anchor.features = features)
save.image("~/data/2021July/2021.09.22.mzlall.cca.RData")

mzl <- IntegrateData(anchorset = anchors)
DefaultAssay(mzl) <- "integrated"

# Run the standard workflow for visualization and clustering
mzl <- ScaleData(mzl, verbose = FALSE)
mzl <- RunPCA(mzl, npcs = 30, verbose = FALSE)
mzl <- RunUMAP(mzl, reduction = "pca", dims = 1:30)
mzl <- FindNeighbors(mzl, reduction = "pca", dims = 1:30)
mzl <- FindClusters(mzl, resolution = 0.5)


# Add metadata - case/control standardization across datasets
mzl$Integrated_case_ctl<- ifelse(mzl$braak<= 3 & mzl$cerad >= 3, "ctl", "NA" )
mzl$Integrated_case_ctl<- ifelse(mzl$braak >= 4 & mzl$cerad <= 2 , "AD", mzl$Integrated_case_ctl )


