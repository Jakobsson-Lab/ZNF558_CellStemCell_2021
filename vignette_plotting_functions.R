# Title     : Functions needed to reproduce plots in Downstreamplotting
# Objective : TODO
# Created by: yo4230sh
# Created on: 24/03/2021

########## fix me: move this to a separate funtion dircetory.

createTable <- function(dataFrameName, columnsToShow) {
  tab1 = as.table(table(dataFrameName[, columnsToShow]))
  tab1
}
showTableAsHeatmap <- function(inputTable, normalize_row=F) {
  require(pheatmap)
  if(normalize_row) {
    z <- rowSums(inputTable)
    w <- round( inputTable/rowSums( inputTable)*1000)/10.0
    new_rr <- rownames(w)
    for(i in 1:dim(w)[1]) {
      new_rr[i] = paste0(rownames(w)[i], " (", as.character(z[i]), ")")
    }
    rownames(w) <- new_rr
  } else {
    w <- inputTable
  }
  pheatmap(w, scale='row', cluster_rows = F, cluster_cols = F, display_numbers=w, legend=F, fontsize_number=10, cellheight = 20)
}


make_vln_plots <- function(cds, outDir) {
  tabInferredTaxonomy2 = createTable(cds@meta.data, c("sample", "inferredTaxonomy2"))
  zz <- round( tabInferredTaxonomy2/rowSums(tabInferredTaxonomy2)*1000)/10.0
  p5 <- showTableAsHeatmap( zz, normalize_row = T) # getting percentage till one decimal point
  ggsave(p5, filename = paste0(outDir, '/cellFractionPerSample.pdf'), width = 10, height = 7)

  tab1 <- createTable(cds@meta.data, c("treatment", "inferredTaxonomy2"))
  p6 <- showTableAsHeatmap(tab1, normalize_row = T)
  ggsave(p6, filename = paste0(outDir, '/cellFractionPD_Control.pdf'), width = 6, height = 3)

  #### Plot some QC for each sample
  v1 <- VlnPlot(cds, features = "nFeature_RNA", pt.size = 0.05, group.by = "inferredTaxonomy2") +
    theme_gray() + xlab("Samples") + ggtitle("Number of Genes Detected") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(v1,filename = paste0(outDir, '/NumberOfGenessPerCellType.pdf'), width = 18, height = 6)

  ## + theme(axis.text.x = element_text(angle = 45))
  ## +  theme (axis.title.x = element_blank())
  v2 <- VlnPlot(cds, features = "nCount_RNA", pt.size = 0.05, group.by = "sample") +
    theme_gray() + xlab("Samples") + ggtitle("Number of UMIs Detected") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(v2,filename = paste0(outDir, '/NumberOfUMIsPerSample.pdf'), width = 18, height = 6)

  v3 <- VlnPlot(cds, features = "nFeature_RNA", pt.size = 0.05, group.by = "sample") +
    theme_gray() + xlab("Samples") + ggtitle("Number of Genes Detected") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(v3,filename = paste0(outDir, '/NumberOfGenesPerSample.pdf'), width = 18, height = 6)

  tabInferredTaxonomy2
}


####################################################################################
#### Markers Per Cell type
####################################################################################

known_markers <- list()
known_markers[["Oligodendrocytes"]] <- c("PLP1","MOG","MBP","PMP2")
known_markers[["Interneurons"]] <- c("GAD1","GAD2","CALB2","CNR1")
known_markers[["Excitatory"]] <- c("RBFOX3","GRIN1","MAP2","HS3ST2")
known_markers[["Astrocytes"]] <- c("AQP4","GJA1","SLC1A3","EDNRB")
known_markers[["Microglia"]] <- c("CTSS","C1QB","P2RY12","ITGAX")
known_markers[["OPCs"]] <- c("COL9A1","VCAN","COL6A1","PDGFRA")
known_markers[["Endothelial"]] <- c("CDH5","CLDN5")
known_markers[["VLMCs"]] <- c("FBLN1","DCN","COL1A1","FBLN1")
known_markers[["EntericNeurons"]] <- c("MCTP2","KLHL1","ITGA11","THEMIS")
known_markers[["Pvalb_Interneurons"]] <- c("PVALB","ERBB4","TAC1","PTHLH")

make_feature_plot_for_known_markers <- function(cell_type, markers_for_this_cell_type, outDir) {
  flog.info("cell_type: %s markers: %s",
            cell_type, paste(markers_for_this_cell_type, sep=' ', collapse = ','))
  p <- FeaturePlot(cds, features = markers_for_this_cell_type,
                   reduction = "umap", cols = c("grey","red"), combine = FALSE)
  q <- lapply(p, function(x) {x  + theme_gray()  + theme(plot.title = element_text(hjust = 0.5))})
  Markersplot1 <- plot_grid(plotlist=q, ncol = 2, nrow = 2)
  file_name <- paste0( cell_type, "Markers.pdf")
  flog.info("saving to file: %s directory: %s", file_name, outDir)
  ggsave(Markersplot1, filename = file.path(outDir, file_name), width = 14, height = 7)
  Markersplot1
}

# Method 1
# markersPlotList <- list()
# for(curr_cell_type in names(known_markers)) {
#   z <- make_feature_plot_for_known_markers(curr_cell_type, known_markers[[curr_cell_type]], outDir)
#   markersPlotList[[curr_cell_type]] <- z
# }

# Method 2
#markersPlotList <- lapply(names(known_markers), function(x) {
#  make_feature_plot_for_known_markers(x, known_markers[[x]], outDir)
#})


getMarkersAndPlot <- function(cds, clusterId) {
  currMarkers = FindMarkers(cds, subset.ident = clusterId,  group.by = "treatment", ident.1 = "TBI")
  currMarkers <- currMarkers[currMarkers$p_val_adj < padj, ]
  if(dim(currMarkers)[1] > 0) {
    kable(currMarkers)
    p <- VlnPlot(cds, split.by = "treatment", idents = clusterId, features = rownames(currMarkers),pt.size = F, same.y.lims = T, combine = T)
    p <- p + labs(x=NULL, y=NULL) + NoLegend()
    #p <- p + labs(x=NULL, y=NULL)
    p
  } else {
    NULL
  }
}
MedianVioSeu <- function() stat_summary(fun.y = mean, geom='point', size = 15, colour = "darkred", shape = 95)

get_matrix_of_overlaps <- function(markers_X, markers_Y) {
  #' Writes out two plots (one showing intersections between each marker list)
  result = matrix(data=0, nrow = length(markers_X), ncol = length(markers_Y))
  for(i in 1:length(markers_X)) {
    for(j in 1:length(markers_Y)) {
      result[i, j] = length(intersect(unique(markers_X[[i]]$gene), unique(markers_Y[[j]]$gene)))
    }
  }
  name_updater <- function(x, n) {paste0(n, " (", dim(x)[1], ")") }
  rownames(result) <- mapply(name_updater, markers_X, names(markers_X))
  colnames(result) <- mapply(name_updater, markers_Y, names(markers_Y))
  result
}


plot_matrix <- function(m, xlab=NULL, ylab=NULL, title=NULL) {
  p <- ggplot(melt(m), aes(Var1,Var2, fill=value, label=round(value, digits=3) )) + geom_tile()
  p <- p + theme(axis.text.x = element_text(angle = 90))
  if(!is.null(xlab)) { p <- p + xlab(xlab)}
  if(!is.null(ylab)) { p <- p + ylab(ylab)}
  if(!is.null(title)) { p <- p + ggtitle(title)}
  p <- p + geom_text(col = "black")
  p <- p + scale_fill_gradient2(low = "gray", mid = "white", high = "red")
  p <- p + scale_color_gradient2(low = "gray", mid = "white", high = "red")

  p
}
