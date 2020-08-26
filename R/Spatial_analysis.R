library(Seurat)
library(MAESTRO)
library(ggplot2)
library(cowplot)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(patchwork)
library(optparse)

option_list = list(
  make_option(c("--data"), type = "character", default = "",
              action = "store", help = "The 10X output path."
  ),
  make_option(c("--prefix"), type = "character", default = "MAESTRO",
              action = "store", help = "The prefix of the output files."
  ),
  make_option(c("--outdir"), type = "character", default = "Result/Analysis",
              action = "store", help = "The directory where the output files are stored."
  ),
  make_option(c("--signature"), type = "character", default = "",
              action = "store", help = "The cell signature file for celltype annotation. Default is built-in CIBERSORT immune cell signature."
  )
)

SpatialRunSeurat<- function(object, project, signatures = "human.immune.CIBERSORT"){
  object@project.name = project
  
  ## ---------------- QC ---------------- ##
  p1 <- VlnPlot(object, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  p2 <- SpatialFeaturePlot(object, features = "nCount_Spatial") + theme(legend.position = "right")
  p <- wrap_plots(p1, p2)
  ggsave(paste0(object@project.name, "_QC_nCount.png"), p,  width=6, height=4.5)

  p1 <- VlnPlot(object, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
  p2 <- SpatialFeaturePlot(object, features = "nFeature_Spatial") + theme(legend.position = "right")
  p <- wrap_plots(p1, p2)
  ggsave(paste0(object@project.name, "_QC_nFeature.png"), p,  width=6, height=4.5)
  
  
  ## ---------------- Clustering -------------- ##
  object <- SCTransform(object, assay = "Spatial", verbose = FALSE)
  object <- RunPCA(object, assay = "SCT", verbose = FALSE)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object, verbose = FALSE)
  object <- RunUMAP(object, reduction = "pca", dims = 1:30)
  p3 <- DimPlot(object, reduction = "umap", label = TRUE)
  p4 <- SpatialDimPlot(object, label = TRUE, label.size = 3)
  p <- wrap_plots(p3, p4)
  ggsave(paste0(object@project.name, "_clustering.png"), p,  width = 10, height=4.5)
  
  cluster.genes <- FindAllMarkersMAESTRO(object)
  cluster.genes <- cluster.genes[cluster.genes$p_val_adj<1E-5, ]
  write.table(cluster.genes, paste0(object@project.name, "_DiffGenes.tsv"), quote = F, sep = "\t")
  
  object = RNAAnnotateCelltype(object, genes = cluster.genes, signatures = signatures)
  p5 <- DimPlot(object, group.by = "assign.ident", reduction = "umap", label = TRUE)
  p6 <- SpatialDimPlot(object, group.by = "assign.ident",label = TRUE, label.size = 3)
  p <- wrap_plots(p5, p6)
  ggsave(paste0(object@project.name, "_annotation.png"), p,  width = 11, height=4.5)  
  
  p1 <- VlnPlot(object, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  p2 <- VlnPlot(object, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
  p <- plot_grid(p1, p2, ncol = 1)
  ggsave(paste0(object@project.name, "_QC_Vlnplot_cluster.png"), p,  width=6, height=6)
  
  return(list(object=object, genes = cluster.genes))
}

argue = parse_args(OptionParser(option_list = option_list, usage = "Run spatial transcriptomic analysis pipeline"))

data_dir = argue$data
SeuratObj = Load10X_Spatial(data_dir)
setwd(argue$outdir)
prefix = argue$prefix
sigfile = argue$signature

Result = SpatialRunSeurat(SeuratObj, project = "10X_brain_coronal", signatures = sigfile)
saveRDS(Result, paste0(Result$object@project.name, "_Spatial_Object.rds"))
