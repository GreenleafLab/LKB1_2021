library(ArchR)

set.seed(1)

groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowMeans(mat[, which(groups == x), drop = F], 
                na.rm = na.rm)
        }
        else {
            rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

#Input Summarized Experiment
se <- readRDS("output/ATAC_Analysis/Main/ATAC-Summarized-Experiment.rds")

#Get Variable Peaks
gM <- groupMeans(assay(se), colData(se)$Group)
rV <- matrixStats::rowVars(gM)
varPeaks <- head(order(rV,decreasing=TRUE), 10000)

#PCA of Log2 Normalized Data
pcaVar <- prcomp(t(assays(se[varPeaks,])[["log2norm"]]))

#PCA Info
pVar <- round(100*pcaVar$sdev^2/sum(pcaVar$sdev^2),3)
pVar <- paste0(paste0("PC", seq_along(pVar))," ",pVar,"%")
plotPCA <- data.frame(pcaVar$x)

#Palette
pal <- paletteDiscrete(values = plotPCA$Group)

#Plot To PDF
pdf("~/Desktop/PCA-10k-Variable-Peaks.pdf", width = 4, height = 4, useDingbats = FALSE)

p <- ggplot(plotPCA, aes(PC1,PC2,color=Group)) + 
	geom_point(size=2.5) + theme_ArchR() + 
	xlab(pVar[1]) + ylab(pVar[2]) +
	scale_color_manual(values = pal)
ArchR:::.fixPlotSize(p, plotWidth = 4, plotHeight = 4, height=1)
grid::grid.newpage()

p <- ggplot(plotPCA, aes(PC1,PC3,color=Group)) + 
	geom_point(size=2.5) + theme_ArchR() + 
	xlab(pVar[1]) + ylab(pVar[3]) +
	scale_color_manual(values = pal)
ArchR:::.fixPlotSize(p, plotWidth = 4, plotHeight = 4, height=1)
grid::grid.newpage()

p <- ggplot(plotPCA, aes(PC2,PC3,color=Group)) + 
	geom_point(size=2.5) + theme_ArchR() + 
	xlab(pVar[2]) + ylab(pVar[3]) +
	scale_color_manual(values = pal)
ArchR:::.fixPlotSize(p, plotWidth = 4, plotHeight = 4, height=1)

dev.off()
