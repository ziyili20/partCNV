# partCNV
An R package to improve the detection of locally aneuploid cells by incorporating cytogenetics information.

# Introduction

Single-cell RNA sequencing is becoming an increasingly common tool to investigate the cellular population and patients' outcome in cancer research. However, due to the sparse data and the complex tumor microenvironment, it is challenging to identify neoplastic cells that play important roles in tumor growth and disease progression. This challenge is exaggerated in the research of blood cancer patients, from whom the neoplastic cells can be highly similar to the normal cells.

In this package we present partCNV/partCNVH, a statistical framework for rapid and accurate detection of aneuploid cells with local copy number deletion or amplification. Our method uses an EM algorithm with mixtures of Poisson distributions while incorporating **cytogenetics information** to guide the classification (partCNV). When applicable, we further improve the accuracy by integrating a Hidden Markov Model for feature selection (partCNVH).


# Understand your cytogenetics data

First of all, our method can only be used to infer cell status if the patient has a part of the cells with cytogenetics alterations. For example, if a patient has cytogenetics data as ``46,XY,del(20)(q11.1q13.1)[5]/46,XY[15]``, it means that out of 20 metaphases, 5 (i.e., 25%) of the cells has deletion on chromosome 20 in region q11.1 to q13.1. If the patient has a normal cytogenetics, e.g., 46,XY[20], or all altered cells, e.g., 46,XY,del(20)(q11.1q13.1)[25], there won't be any need to apply the proposed method. 

Second, when you have a complicated cytogenetics feature, use them one by one to identify the desired cell group. For example, ``47,XY,+8[5]/46,idem,del(8)(q21.2q24.3)/46,XY[7]``, we start with the chromosome amplification on chromosome 8 excluding the region (q21.2q24.3) to identify the cells with this alteration using the proposed method. After the cells with chromosome 8 amplification are identified, we apply the proposed method to identify cells with del(8)(q21.2q24.3). 

Lastly, our tool is primarily designed for analyzing the scRNA-seq data from patients with hematologic malignancies, such as MDS and AML. For solid tumors such as lung cancer and breast cancer, it is better to use CNV based method such as [copyKAT](https://github.com/navinlabcode/copykat) and inferCNV `r Biocpkg("inferCNV")`.


# Load example data

The analysis can start from whole-genome scRNA-seq data or a subset of scRNA-seq data based on the location of interest. In this package, we prepared a whole-genome scRNA-seq data:


<!-- Assign captions to figures in the code chunk option `fig.cap` to automatically number them, and to be able to reference them, see Figure \@ref(fig:plot). The figure label is generated from the code chunk label by prefixing it with `fig:`. -->

```{r loadData, echo=TRUE}
library(partCNV)
data(SimData)
dim(SimData)
SimData[1:5,1:5]
```

This simple example dataset contains 24519 genes and 400 cells. If you start with Seurat object (after quality control), run:

```{r runseurat, echo=TRUE, eval=FALSE}
library(Seurat)
Seurat_obj <- NormalizeData(Your_SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
Counts = Seurat_obj@assays$RNA@counts
```

# Run partCNV

Since this is simulation data, prior knowledge (e.g., cytogenetics data in real studies) shows that this example data has 40% of cells with deletion on chromosome 20 q11.1 to q13.1. Let's get started with locating this region.

```{r s1, echo=TRUE, eval=TRUE}
res <- GetCytoLocation(cyto_feature = "chr20(q11.1-q13.1)")
```

Then let's subset the data and normalize the data:
```{r s2, echo=TRUE, eval=TRUE}
GEout <- GetExprCountCyto(cytoloc_output = res, Counts = SimData, normalization = TRUE, qt_cutoff = 0.99)
```
For this function, the qt_cutoff is to filter out the cells with very low expressions. 0.99 here means that we filter out cells that only express in 1% (1-0.99) cells. Make this `qt_cutoff` larger if your total gene number within the region is small. 

Now we apply partCNV:
```{r s3, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, results='hide'}
pcout <- partCNV(int_counts = GEout$ProcessedCount,
                 cyto_type = "del",
                 cyto_p = 0.40)
```

Understand the results:
```{r s4, echo=TRUE, eval=TRUE}
table(pcout)
sum(pcout==1)/length(pcout)
p1 <- sum(pcout==1)/length(pcout)
```
39.5% of cells are labeled as locally aneuploid (1) and others are diploid (0).

Let's visualize it:
```{r s5, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, results='hide'}
library(Seurat)
sim_seurat <- CreateSeuratObject(counts = SimData)
sim_seurat <- NormalizeData(sim_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
sim_seurat <- FindVariableFeatures(sim_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sim_seurat)
sim_seurat <- ScaleData(sim_seurat, features = all.genes)
sim_seurat <- RunPCA(sim_seurat, features = VariableFeatures(object = sim_seurat))
sim_seurat <- RunUMAP(sim_seurat, dims = 1:10)
sim_seurat <- AddMetaData(
    object = sim_seurat,
    metadata = pcout,
    col.name = "partCNV_label"
)
sim_seurat$partCNV_label <- factor(sim_seurat$partCNV_label, levels = c(1,0), labels = c(
    "aneuploid", "diploid"
))
```

```{r s6, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE}
library(ggplot2)
DimPlot(sim_seurat, reduction = "umap", group = "partCNV_label") + ggtitle(paste0("partCNV (", signif(p1,2)*100, "%)"))+ theme(legend.position = "none")
```

# Run partCNVH

Compared with partCNV, partCNVH added an additional step of feature selection. This is especially helpful if your cytogenetics provide a very broad region and part of it does not have chromosomal alterations. The first two steps of using partCNVH are exactly the same as using partCNV.


```{r s11, echo=TRUE, eval=TRUE}
res <- GetCytoLocation(cyto_feature = "chr20(q11.1-q13.1)")
GEout <- GetExprCountCyto(cytoloc_output = res, Counts = SimData, normalization = TRUE, qt_cutoff = 0.99)
```

For this function, the qt_cutoff is to filter out the cells with very low expressions. 0.99 here means that we filter out cells that only express in 1% (1-0.99) cells. Make this `qt_cutoff` larger if your total gene number within the region is small. 

Now we apply partCNVH:
```{r s13, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, results='hide'}
pcHout <- partCNVH(int_counts = GEout$ProcessedCount,
                 cyto_type = "del",
                 cyto_p = 0.40)
```

Understand the results (in pcHout, EMlabel is the partCNV label and EMHMMlabel is the partCNVH label).
```{r s14, echo=TRUE, eval=TRUE}
table(pcHout$EMHMMlabel)
sum(pcHout$EMHMMlabel==1)/length(pcHout$EMHMMlabel)
p2 <- sum(pcHout$EMHMMlabel==1)/length(pcHout$EMHMMlabel)
```
41.25% of cells are labeled as locally aneuploid (1) and others are diploid (0).

Let's visualize it:
```{r s15, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, results='hide'}
# I commented these steps because they are exactly the same as partCNV run. 
# library(Seurat)
# sim_seurat <- CreateSeuratObject(counts = SimData)
# sim_seurat <- NormalizeData(sim_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
# sim_seurat <- FindVariableFeatures(sim_seurat, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(sim_seurat)
# sim_seurat <- ScaleData(sim_seurat, features = all.genes)
# sim_seurat <- RunPCA(sim_seurat, features = VariableFeatures(object = sim_seurat))
# sim_seurat <- RunUMAP(sim_seurat, dims = 1:10)
sim_seurat <- AddMetaData(
    object = sim_seurat,
    metadata = pcHout$EMHMMlabel,
    col.name = "partCNVH_label"
)
sim_seurat$partCNVH_label <- factor(sim_seurat$partCNVH_label, levels = c(1,0), labels = c(
    "aneuploid", "diploid"
))
```

```{r s16, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE}
library(ggplot2)
DimPlot(sim_seurat, reduction = "umap", group = "partCNVH_label") + ggtitle(paste0("partCNVH (", signif(p2,2)*100, "%)"))+ theme(legend.position = "none")
```
