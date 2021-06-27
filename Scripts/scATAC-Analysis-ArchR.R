library(ArchR)

set.seed(1)

#Get Input Files Valid Barcodes see 10x Info
validBC <- getValidBarcodes(
	csvFiles = list.files("/Volumes/JG_SSD_2/Analysis/LKB1_Mouse", pattern = ".csv", full.names = TRUE),
	sampleNames = gsub(".singlecell.csv", "", list.files("/Volumes/JG_SSD_2/Analysis/LKB1_Mouse", pattern = ".csv", full.names = FALSE))
)
inputFiles <- getInputFiles("/Volumes/JG_SSD_2/Analysis/LKB1_Mouse/")

#Set Defaults
addArchRGenome("mm10")
addArchRThreads(10)

#Create Arrows
ArrowFiles <- createArrowFiles(
	inputFiles = inputFiles,
	validBarcodes = validBC
)

#Identify Doublets
doubletScores <- addDoubletScores(list.files(pattern = ".arrow"))

#Create ArchR Project
proj <- ArchRProject(list.files(pattern = ".arrow"))

#Filter Doublets
proj <- filterDoublets(proj)

#Filter Cells
proj <- proj[proj$TSSEnrichment >= 6]
proj <- proj[proj$nFrags <= 50000]

#QC Plots
p1 <- plotGroups(proj, colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.35, pal = paletteDiscrete(proj$Sample, set="grove")) + 
	scale_color_manual(values = paletteDiscrete(proj$Sample, set = "grove")) + scale_fill_manual(values = paletteDiscrete(proj$Sample, set = "grove"))
p2 <- plotGroups(proj, colorBy = "cellColData", name = "nFrags", plotAs = "violin", alpha = 0.35, pal = paletteDiscrete(proj$Sample, set="grove")) + 
	scale_color_manual(values = paletteDiscrete(proj$Sample, set = "grove")) + scale_fill_manual(values = paletteDiscrete(proj$Sample, set = "grove"))
plotPDF(p1, p2, name = "Plot-QC-Stats.pdf", width = 4, height = 4, ArchRProj = proj)

#Reduce Dims
proj <- addIterativeLSI(proj,
	iterations = 2,
	clusterParams = list(
		resolution = c(0.1), 
        sampleCells = 10000, 
        n.start = 10
    ),
    varFeatures = 25000,
    excludeChr = c("chrX", "chrY"),
    force = TRUE
)

#Add Clusters
proj <- addClusters(proj, resolution = 0.4, force = TRUE, reducedDims = "IterativeLSI")

#Add UMAP
proj <- addUMAP(proj, force = TRUE, reducedDims = "IterativeLSI", minDist = 0.6, name = "UMAP4")

#Add TSNE
proj <- addTSNE(proj, force = TRUE, reducedDims = "IterativeLSI")

#Add Genotype Cluster Info
proj$Genotype <- stringr::str_split(proj$Sample, pattern = "\\_", simplify = TRUE)[,2]
proj$GenotypeCluster <- paste0(proj$Clusters, ".", proj$Genotype)

#Rename Clusters
labels <- paste0("Cluster", 1:12)
names(labels) <- paste0("Cluster", c(7, 6, 8, 9, 11, 5, 4, 10, 3, 1, 2, 12))
proj$Clusters2 <- mapLabels(proj$Clusters, labels, names(labels))

#Plot UMAP
p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "grove")
p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "bear")
p3 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "calm")
p4 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "kelly")
p5 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "paired")
p6 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", size = 0.5, discreteSet = "circus")
plotPDF(p1, p2, p3, p4, p5, p6, name = "Plot-UMAP", width = 4, height = 4, ArchRProj = proj)



p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters2", size = 0.5, pal = palLKB1, labelMeans = FALSE, embedding = "UMAP2")
p3 <- plotEmbedding(proj, colorBy = "cellColData", name = "Genotype", pal = c("KPL" = "firebrick3", "KP" = "dodgerblue3"), size = 0.5)
p4 <- plotEmbedding(proj, colorBy = "cellColData", name = "log10(nFrags)")
p5 <- plotEmbedding(proj, colorBy = "cellColData", name = "TSSEnrichment")
p6 <- plotEmbedding(proj, colorBy = "cellColData", name = "DoubletEnrichment")
plotPDF(p1, p2, p3, p4, p5, p6, name = "Plot-UMAP", width = 4, height = 4, ArchRProj = proj)

#Impute Weights
proj <- addImputeWeights(proj)

#######################################################################################
# Analysis
#######################################################################################

#Plot Genes
markerGenes <- c("Sox9", "Sox17", "Sox2", "Nr4a1", "Nr4a2", "Nr4a3", "Il33", "Sik1", "Nkx2-1", "Hmga2", "Sox21", "Sox18")
p1 <- plotEmbedding(
	proj, 
	bins = 150, 
	continuousSet = "comet",
	colorBy = "GeneScoreMatrix", 
	name = markerGenes, 
	imputeWeights = getImputeWeights(proj),
	quantCut = c(0, 1)
)
plotPDF(p1, name = "Plot-UMAP-Markers-Comet", width = 4, height = 4, ArchRProj = proj)

markerGenes <- c("Sox9", "Sox17", "Sox2", "Nr4a1", "Nr4a2", "Nr4a3", "Il33", "Sik1", "Nkx2-1", "Hmga2", "Sox21", "Sox18")
p1 <- plotEmbedding(
	proj, 
	bins = 150, 
	continuousSet = "horizonExtra",
	colorBy = "GeneScoreMatrix", 
	name = markerGenes, 
	imputeWeights = getImputeWeights(proj)
)
plotPDF(p1, name = "Plot-UMAP-Markers-Horizon", width = 4, height = 4, ArchRProj = proj)

#######################################################################################
# Save Results
#######################################################################################
#saveRDS(proj, "Save-Project.rds")
#save.image("Save-ArchR-Feb21-807.Rdata")


#######################################################################################
# Call Peaks
#######################################################################################
proj <- addGroupCoverages(proj, groupBy = "Clusters2")
proj <- addReproduciblePeakSet(proj, groupBy = "Clusters2", pathToMacs2 = "/Users/jeffreygranja/Library/Python/2.7/bin/macs2")
proj <- addPeakMatrix(proj)


#######################################################################################
# Add Peak Annotation
#######################################################################################
proj <- addMotifAnnotations(proj, name = "Motif1", version = 1, force = TRUE)


#######################################################################################
# Add Motif Deviations
#######################################################################################
proj <- addDeviationsMatrix(proj, peakAnnotation = "Motif1", matrixName = "Motif1Matrix", force = TRUE)

#######################################################################################
# Motif Deviations
#######################################################################################

motifs <- grep("z:",getFeatures(proj, select = "Sox17|Nkx2|Runx|Tead|Smarcc|Jun|Fos", useMatrix = "Motif1Matrix"),value=TRUE)
motifs <- sort(motifs)

markerGenes <- c("Sox9", "Sox17", "Sox2", "Nr4a1", "Nr4a2", "Nr4a3", "Il33", "Sik1", "Nkx2-1", "Hmga2", "Sox21", "Sox18")
p1 <- plotEmbedding(
	proj, 
	bins = 150, 
	continuousSet = "solarExtra",
	colorBy = "Motif1Matrix", 
	name = motifs, 
	imputeWeights = getImputeWeights(proj),
	quantCut = c(0, 1)	
)
plotPDF(p1, name = "Plot-UMAP-MotifMarkers", width = 4, height = 4, ArchRProj = proj)

#######################################################################################
# Motif Footprint
#######################################################################################

pal <- paletteDiscrete(gtools::mixedsort(unique(proj$Clusters2)))
pal[1] <- pal["Cluster12"]
pal["Cluster12"] <- "#3D3D3D"
palLKB1 <- pal

seSubtract <- plotFootprints(
  proj, 
  inputSE = seNone,
  pal = palLKB1,
  normMethod = "subtract",
  groupBy = "Clusters2",
  positions = getPositions(proj, "Motif1")[c("Sox17_843", "Onecut2_1054", "Nkx21_446")]
)