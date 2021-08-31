# This code will generate figures related to single nuclei RNA sequencing data obtained from Organoids.
# We begin with loading pre-processed RDS file that contains merged Seurat object for all samples.

suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
suppressMessages(require(gridExtra))
suppressMessages(require(ggplot2))
suppressMessages(require(cowplot))
suppressMessages(require(pheatmap))
suppressMessages(require(stringr))
suppressMessages(require(knitr))
suppressMessages(require(data.table))
suppressMessages(require(reticulate))
suppressMessages(require(xlsx))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(dplyr))
source("vignette_plotting_functions.R")
suppressMessages(require(DoubletFinder))

gDriveDir <- '/Users/yo4230sh/GDrive/'
dataDir <- paste0(gDriveDir, 'Seq110_ZNF558_Organoids/')

inputDir <- ./Data/
inputRDSfile <- file.path(inputDir, "results-merged.rds")
outDir <- file.path('/tmp))
dir.create(outDir, recursive = T,showWarnings = F)

### Reading RDS file

rdsData <- readRDS(inputRDSfile)
cds <- rdsData[[1]]
markers <- rdsData[[2]]

metaDf <- cds@meta.data
metaDf$treatment <- sapply(str_split(metaDf$sample, '_'), function(x) { x[2] })
metaDf$timepoint <-  sapply(str_split(metaDf$sample, '_'), function(x) { x[3] })
cds@meta.data <- metaDf

### Adding condition
condition <- rep("gCTR", length(cds@active.ident))
condition[cds@meta.data$treatment == 'g2'] <- "KD g2_g3"
condition[cds@meta.data$treatment == 'g3'] <- "KD g2_g3"
table(condition)
cds@meta.data$condition <- factor(condition)

### save annotated clusters
inferredTax2 <- c("PG_1", "Neurons", "PG_2", "New neurons", "NPCs",
                  "Unknown_1", "Unknown_2", "Unknown_3")
elems <- unique(inferredTax2)
clusterIdsForTaxonomy <- lapply(elems, function(x) { which(x == inferredTax2) - 1 })
names(clusterIdsForTaxonomy) <- elems

### setup the identity
cds@meta.data$inferredTaxonomy2 <- inferredTax2[cds@meta.data$umap_0.1]

#### Figure 6D,6E, supplementary figure 6D  

p2 <- DimPlot(cds, reduction="umap", group.by = "inferredTaxonomy2", label=T, repel=T) + NoLegend() + theme_classic()
p2 <- p2 +   ggtitle("Celltypes characterisation") + theme(plot.title = element_text(hjust = 0.5))
ggsave(p2, filename = paste0(outDir, '/celltypes_znf558.pdf'), width = 9, height = 6)

p6 <- DimPlot(cds, reduction = "umap", group.by = "timepoint") + NoLegend() + theme_classic()
p6 <- p6 +   ggtitle("Timepoints distribution") + theme(plot.title = element_text(hjust = 0.5))
ggsave(p6, filename = paste0(outDir, '/timepoint_harmony.pdf'), width = 9, height = 6)


p4 <- DimPlot(cds, reduction = "umap", group.by = "condition") + NoLegend() + theme_classic()
p4 <- p4 +   ggtitle("Merged Guides and control distribution") + theme(plot.title = element_text(hjust = 0.5))
ggsave(p4, filename = paste0(outDir, '/Mergedguide_control_znf558.pdf'), width = 9, height = 6)


####################################################################################
##### Identify markers that are differentially expressed in LacZ versus combined guides
####################################################################################

p_adj <- 0.01
clusterIdsForTaxonomy_subset <- clusterIdsForTaxonomy[- c(7,8)]

#names(data.integrated_m@meta.data$inferredTaxonomy2)

markers_tbi <- lapply(names(clusterIdsForTaxonomy_subset), function(curr_cell_type, p_adj=0.01) {
  flog.info("Finding organoids markers for cell_type: %s", curr_cell_type )
  curr_cluster_ids <- clusterIdsForTaxonomy_subset[[curr_cell_type]]
  flog.info("ClusterIDs: %s", paste(curr_cluster_ids, collapse=','))

  curr_markers <- FindMarkers(cds, subset.ident = curr_cluster_ids, group.by = "treatment", ident.1 = "LacZ")
  flog.info("cell_type: %s found markers: %d", curr_cell_type, dim(curr_markers)[1])
  curr_markers <- curr_markers[curr_markers$p_val_adj < p_adj, ]
  flog.info("filtered markers: %d at p_val_adj: %f", dim(curr_markers)[1], p_adj)
  curr_markers$cell_type <- curr_cell_type
  curr_markers$gene <- rownames(curr_markers)
  curr_markers
})
names(markers_tbi) <- names(clusterIdsForTaxonomy_subset)
markers_tbi_table_guide_combined <- Reduce(rbind, markers_tbi)
write.table(markers_tbi_table_guide_combined, file = file.path(outDir, "CombinedGuides_Markers.tsv"), row.names = F, quote = F, sep = '\t')

#### adding cell cycles scores

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cds <- CellCycleScoring(cds, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cyclingplot <- FeaturePlot(cds, features = c("S.Score","G2M.Score"), reduction = "umap", cols = c("lightblue","red"))
ggsave(cyclingplot, filename = paste0(outDir, '/CellcyclingScores.png'), width = 9, height = 6)

cycling <- rep("no", length(cds@active.ident))
cycling[cds@meta.data$Phase == 'G2M'] <- "yes"
cycling[cds@meta.data$Phase == 'S'] <- "yes"
table(cycling)

cds@meta.data$cycling <- factor(cycling)

cycling_status <- DimPlot(cds,group.by = "cycling", split.by = "condition") + theme_gray()
cycling_status <- cycling_status + ggtitle("Cycling Non cycling cells") + theme(plot.title = element_text(hjust = 0.5))
ggsave(cycling_status, filename = paste0(outDir, '/Yes_No_cellcycles.pdf'), width = 14, height = 6)

### Supplementary figure 6F

test <- subset(x = cds, subset = SPATA18 > 0)

spata18_treatment <- VlnPlot(test, features = c("SPATA18"), group.by = "inferredTaxonomy2",  split.by = "condition") + stat_summary(fun = median, geom='point', size = 5, shape = 95, colour = "blue",
                                                                                                   position = position_dodge(1)) + scale_fill_brewer(palette="Blues") + theme_classic()
ggsave(spata18_treatment, filename = paste0(outDir, '/SPATA18_treatment.pdf'), width = 14, height = 6)

### Figure 6E, supplementary figure 6C
clusterId <- c("1","3","4")
dotplotgenelist <- c("MKI67","PAX6","SOX2","GLI3","MAP2","CUX2","NCAM1","RBFOX3","GFAP")
dotplot_selectedcluster <- DotPlot(cds, features = dotplotgenelist, idents = clusterId,
              group.by = "inferredTaxonomy2") + RotatedAxis() + theme_classic() + coord_flip()
dotplot_selectedcluster <- dotplot_selectedcluster +  xlab("Genes") +  ylab("Celltypes") +
  ggtitle("Candidate marker expression") + theme(plot.title = element_text(hjust = 0.5))

ggsave(dotplot_selectedcluster, filename = paste0(outDir, '/dotplot_selectedmarkers.pdf'), width = 8, height = 7)

dotplot_allcluster_genes <- c("MKI67","SOX2","GLI3","NES","OTX2","MAP2","CUX2","RBFOX3","GFAP", "AQP4","SLC1A3","OLIG2","VCAN","MYLPF","MYL1","CXCR4","SOX17","POU5F1","NANOG","AIF1")

dotplot_allcluster <- DotPlot(cds, features = dotplot_allcluster_genes,
              group.by = "inferredTaxonomy2") + RotatedAxis() + theme_classic() + coord_flip()
dotplot_allcluster <- dotplot_allcluster +  xlab("Genes") +  ylab("Celltypes") +
  ggtitle("Candidate marker expression") + theme(plot.title = element_text(hjust = 0.5))

ggsave(dotplot_allcluster, filename = paste0(outDir, '/dotplot_allcluster.pdf'), width = 12, height = 9)

### Figure 6F

clusterId <- c(3)
Newbornneurons_glist <- c("DCDC2","CEP126","CDC14A","DNAI1","CUX2","MAP2")

Newborn_violin <- VlnPlot(cds, features = Newbornneurons_glist, group.by = "condition", idents = clusterId, pt.size = 0.01, combine = FALSE)
Newborn_violin <- lapply(Newborn_violin, function(x) {x  +  stat_summary(fun = median, geom='point', size = 5, shape = 95, colour = "blue",
                                                                         position = position_dodge(1))})
Newborn_violin <-  lapply(Newborn_violin, function(x) {x  +  scale_fill_brewer(palette="Blues") + theme_classic()})
  Newborn_violin <- lapply(Newborn_violin, function(x) {x  + theme(plot.title = element_text(hjust = 0.5))})
testplot_newborn <- plot_grid(plotlist=Newborn_violin, ncol = 3, nrow = 3)
  ggsave(testplot_newborn, filename = paste0(outDir, '/NewbornViolin.pdf'), width = 12, height = 8)

clusterId <- c(1)
Neurons_glist <- c("NCAM1","MYT1L","MAP2","RBFOX3","GRIA2","CUX2")
Neurons_violin <- VlnPlot(cds, features = Neurons_glist, group.by = "condition", idents = clusterId, pt.size = 0.01, combine = FALSE)
Neurons_violin <- lapply(Neurons_violin, function(x) {x  +  stat_summary(fun = median, geom='point', size = 5, shape = 95, colour = "blue",
                                                                         position = position_dodge(1))})
Neurons_violin <-  lapply(Neurons_violin, function(x) {x  +  scale_fill_brewer(palette="Blues") + theme_classic()})
  Neurons_violin <- lapply(Neurons_violin, function(x) {x  + theme(plot.title = element_text(hjust = 0.5))})
testplot_neurons <- plot_grid(plotlist=Neurons_violin, ncol = 3, nrow = 3)
  ggsave(testplot_neurons, filename = paste0(outDir, '/NeuronsViolin.pdf'), width = 12, height = 8)

clusterId <- c(4)
NPCs_glist <- c("BCAS3","SMARCA2","FLRT2","CRB1","FOXP2","TOX")
NPCs_violin <- VlnPlot(cds, features = NPCs_glist, group.by = "condition", idents = clusterId, pt.size = 0.01, combine = FALSE)
NPCs_violin <- lapply(NPCs_violin, function(x) {x  +  stat_summary(fun = median, geom='point', size = 5, shape = 95, colour = "blue",
                                                                         position = position_dodge(1))})
NPCs_violin <-  lapply(NPCs_violin, function(x) {x  +  scale_fill_brewer(palette="Blues") + theme_classic()})
  NPCs_violin <- lapply(NPCs_violin, function(x) {x  + theme(plot.title = element_text(hjust = 0.5))})
testplot_npcs <- plot_grid(plotlist=NPCs_violin, ncol = 3, nrow = 3)
  ggsave(testplot_npcs, filename = paste0(outDir, '/NPCsViolin.pdf'), width = 12, height = 8)

### Figure 6E
znf558 <- FeaturePlot(cds, features = "ZNF558", reduction = "umap", split.by = "condition", pt.size = 0.5)
ggsave(znf558, filename = paste0(outDir, '/ZNF558_expression.pdf'), width = 11, height = 5)
