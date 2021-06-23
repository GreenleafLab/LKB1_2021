library(ArchR)
library(ArchRx)

setwd("~/Documents/LKB1/Human_OE_Exp2/")

outDiri <- "~/Desktop/Human"
dir.create(outDiri)

se <- readRDS("output/ATAC_Analysis/Main/ArchRSE_Log2Norm.rds")
seCV <- readRDS("output/ATAC_Analysis/Main/CisBP-chromVAR.rds")
tfMatches <- readRDS("output/ATAC_Analysis/Main/CisBP-Motif-Matches.rds")
geneMatches <- TxDbGeneMatches(rowRanges(se), "hg38")


######################
# PCA
######################

set.seed(1)
se <- se[,grep("GFP",colnames(se))]
gM <- groupMeans(assay(se), colData(se)$Group)
rV <- matrixStats::rowVars(gM)
varPeaks <- head(order(rV,decreasing=TRUE), 10000)

#gM <- groupMeans(assay(se), colData(se)$Group)[varPeaks, ]
#pcaVar <- prcomp(t(gM))
pcaVar <- prcomp(t(assays(se[varPeaks,])[["log2norm"]]))

pVar <- round(100*pcaVar$sdev^2/sum(pcaVar$sdev^2),3)
pVar <- paste0(paste0("PC", seq_along(pVar))," ",pVar,"%")
plotPCA <- data.frame(pcaVar$x)
plotPCA$Group <- gsub("_Rep1|_Rep2", "",rownames(plotPCA))

pal <- paletteDiscrete(values = plotPCA$Group)

pdf("~/Desktop/Human/PCA-10k-Variable-Peaks-2.pdf", width = 4, height = 4, useDingbats = FALSE)
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


pdf("~/Desktop/Human/Dist-PCA-ATAC.pdf", width = 5.5, height = 5)

cm1 <- as.dist(philentropy::distance(pcaVar$x[,1:3], method = "euclidean"))
cm <- as.matrix(cm1)
diag(cm) <- NA
rownames(cm) <- colnames(se)
colnames(cm) <- colnames(se)

pheatmap::pheatmap(cm, 
         clustering_distance_rows = cm1, 
         clustering_distance_cols = cm1,
					color = paletteContinuous(), 
					border_color = "black", 
         cutree_cols = 2,
         cutree_rows = 2,
         show_colnames = T,
         main = "EuclideanDistance")

ct <- cutree(hclust(cm1), 2)
names(ct) <- colnames(se)

cm1 <- as.dist(philentropy::distance(pcaVar$x[,1:3], method = "cosine"))
cm <- as.matrix(cm1)
diag(cm) <- NA
rownames(cm) <- colnames(se)
colnames(cm) <- colnames(se)

pheatmap::pheatmap(cm, 
         clustering_distance_rows = 1-cm1, 
         clustering_distance_cols = 1-cm1,
         cutree_cols = 2,
         cutree_rows = 2,
					color = paletteContinuous(), 
					border_color = "black", 
         show_colnames = T,
         main = "CosineSimilarity")

dev.off()



cm1 <- dist(t(groupMeans(assay(se),colData(se)$Group)))
cm <- as.matrix(cm1)
diag(cm) <- NA
pdf("~/Desktop/Human/Correlation-ATAC-2.5.pdf", width = 5.5, height = 3.75)

pheatmap::pheatmap(cm, 
         clustering_distance_rows = cm1, 
         clustering_distance_cols = cm1,
					color = paletteContinuous(), 
					border_color = "black", 
         show_colnames = F)

dev.off()


cm <- as.dist(philentropy::distance(t(gM), method = "kulczynski_d"))

cm1 <- as.dist(cor(groupMeans(assay(se),colData(se)$Group)))
cm <- as.matrix(cm1)
diag(cm) <- NA
pdf("~/Desktop/Human/Correlation-ATAC-3.pdf", width = 5.5, height = 3.75)

pheatmap::pheatmap(cm, 
         clustering_distance_rows = 1-cm1, 
         clustering_distance_cols = 1-cm1,
					color = paletteContinuous(), 
					border_color = "black", 
         show_colnames = F)

dev.off()


 [1] "euclidean"         "manhattan"         "minkowski"         "chebyshev"        
 [5] "sorensen"          "gower"             "soergel"           "kulczynski_d"     
 [9] "canberra"          "lorentzian"        "intersection"      "non-intersection" 
[13] "wavehedges"        "czekanowski"       "motyka"            "kulczynski_s"     
[17] "tanimoto"          "ruzicka"           "inner_product"     "harmonic_mean"    
[21] "cosine"            "hassebrook"        "jaccard"           "dice"             
[25] "fidelity"          "bhattacharyya"     "hellinger"         "matusita"         
[29] "squared_chord"     "squared_euclidean" "pearson"           "neyman"           
[33] "squared_chi"       "prob_symm"         "divergence"        "clark"            
[37] "additive_symm"     "kullback-leibler"  "jeffreys"          "k_divergence"     
[41] "topsoe"            "jensen-shannon"    "jensen_difference" "taneja"           
[45] "kumar-johnson"     "avg"


pdf("~/Desktop/Human/Correlation-ATAC-4.pdf", width = 5.5, height = 3.75)
cm <- 1-cor(assay(se))
diag(cm) <- NA

pheatmap::pheatmap(
	cm, 
	color = paletteContinuous(), 
	border_color = "black", 
	show_colnames = F
)
dev.off()


diffATAC <- readRDS("output/ATAC_Analysis/Main/EdgeR-Pairwise-Summary.rds")

diffMat <- diffATAC %>% 
	{rbind(
		(assays(.)[["log2FoldChange"]] >= 0.5 & assays(.)[["padj"]] <= 0.05) %>% colMeans,
		(assays(.)[["log2FoldChange"]] <= -0.5 & assays(.)[["padj"]] <= 0.05) %>% colMeans
		)
	}

rownames(diffMat) <- c("Up", "Down")
diffMat[2,] <- -1 * diffMat[2,]
df <- rbind(
	reshape2::melt(diffMat[,grep("LKB1", colnames(diffMat))])
)

df$Var2 <- factor(paste0(df$Var2), 
	levels = c(
		c("H1650_LKB1_GFP", "H1975_LKB1_GFP", "H358_LKB1_GFP", "H2009_LKB1_GFP", "H1437_LKB1_GFP", "A549_LKB1_GFP", "H460_LKB1_GFP", "H1355_LKB1_GFP")
	  )
	)

pdf("~/Desktop/Human/Relative-ATAC-Differences-LKB1.pdf", width = 4, height = 5)
ggplot(df, aes(Var2, value, fill=paste0(Var1))) + 
	geom_bar(stat = "identity") + 
	ArchR:::theme_ArchR(xText90 = TRUE) + 
	scale_fill_manual(values = paletteDiscrete(values=df$Var1)) +
	ylab("Proportion of Differential Peaks") + 
	ylim(c(-0.06, 0.12))
dev.off()
 



diffATAC <- readRDS("output/ATAC_Analysis/Main/EdgeR-Pairwise-Summary.rds")

diffMat <- diffATAC %>% 
	{rbind(
		(assays(.)[["log2FoldChange"]] >= 0.5 & assays(.)[["padj"]] <= 0.05) %>% colMeans,
		(assays(.)[["log2FoldChange"]] <= -0.5 & assays(.)[["padj"]] <= 0.05) %>% colMeans
		)
	}

rownames(diffMat) <- c("Up", "Down")
diffMat[2,] <- -1 * diffMat[2,]
df <- rbind(
	reshape2::melt(diffMat[,grep("KEAP1", colnames(diffMat))])
)

df$Var2 <- factor(paste0(df$Var2), 
	levels = c(
		gsub("LKB1", "KEAP1",c("H1650_LKB1_GFP", "H1975_LKB1_GFP", "H358_LKB1_GFP", "H2009_LKB1_GFP", "H1437_LKB1_GFP", "A549_LKB1_GFP", "H460_LKB1_GFP", "H1355_LKB1_GFP"))
	  )
	)

pdf("~/Desktop/Human/Relative-ATAC-Differences-KEAP1.pdf", width = 4, height = 5)
ggplot(df, aes(Var2, value, fill=paste0(Var1))) + 
	geom_bar(stat = "identity") + 
	ArchR:::theme_ArchR(xText90 = TRUE) + 
	scale_fill_manual(values = paletteDiscrete(values=df$Var1)) +
	ylab("Proportion of Differential Peaks") + 
	ylim(c(-0.06, 0.12))
dev.off()
 
















diffMat <- diffATAC %>% 
	{rbind(
		(assays(.)[["log2FoldChange"]] >= 0.5 & assays(.)[["padj"]] <= 0.05) %>% colMeans,
		(assays(.)[["log2FoldChange"]] <= -0.5 & assays(.)[["padj"]] <= 0.05) %>% colMeans
		)
	}

rownames(diffMat) <- c("Up", "Down")
diffMat[2,] <- -1 * diffMat[2,]
df <- rbind(
	reshape2::melt(diffMat[,grep("KEAP1", colnames(diffMat))])
)

df$Var2 <- factor(paste0(df$Var2), 
	levels = c(
		gsub("LKB1", "KEAP1", c("H1650_LKB1_GFP", "H1975_LKB1_GFP", "H358_LKB1_GFP", "H2009_LKB1_GFP", "H1437_LKB1_GFP", "A549_LKB1_GFP", "H460_LKB1_GFP", "H1355_LKB1_GFP"))
	  )
	)

for(i in seq_along(colnames(diffATAC))){

	print(i)

	pdf(paste0("~/Desktop/Human/",colnames(diffATAC)[i],".pdf"), width = 4, height = 4)
	ggPlotDifferential(
		log2Mean = assays(diffATAC[,colnames(diffATAC)[i]])[["log2Mean"]][,1],
		log2FoldChange = assays(diffATAC[,colnames(diffATAC)[i]])[["log2FoldChange"]][,1],
		padj = assays(diffATAC[,colnames(diffATAC)[i]])[["padj"]][,1],
		padjCutOff = 0.05,
		log2FCCutOff = 0.5,
		ratioYX = 3.5/4,
		ylim = c(-6, 6),
		pal = c("Increasing"="dodgerblue3","Decreasing"="firebrick3","Constitutive"="lightgrey")) %>% 
		ArchR:::.fixPlotSize(plotWidth = 3.5, plotHeight=4) %>% print
	dev.off()

}


H1650, H1975, H358, H2009, H1437, A549, H460, H1355



mUp <- assays(diffATAC)[["log2FoldChange"]] >= 0.5 & assays(diffATAC)[["padj"]] <= 0.05
mDo <- assays(diffATAC)[["log2FoldChange"]] <= -0.5 & assays(diffATAC)[["padj"]] <= 0.05

tabUp <- table(rowSums(mUp[,c("H1437_LKB1_GFP", "A549_LKB1_GFP", "H460_LKB1_GFP", "H1355_LKB1_GFP")]))
tabDo <- table(rowSums(mDo[,c("H1437_LKB1_GFP", "A549_LKB1_GFP", "H460_LKB1_GFP", "H1355_LKB1_GFP")]))

sum(tabUp[c(3,4,5)]) / sum(tabUp[c(2,3,4,5)])
sum(tabDo[c(3,4,5)]) / sum(tabDo[c(2,3,4,5)])

v1 <- (sum(tabUp[c(3,4,5)]) + sum(tabDo[c(3,4,5)])) / (sum(tabUp[c(2,3,4,5)]) + sum(tabDo[c(2,3,4,5)]))
v2 <- (sum(tabUp[c(4,5)]) + sum(tabDo[c(4,5)])) / (sum(tabUp[c(2,3,4,5)]) + sum(tabDo[c(2,3,4,5)]))
v3 <- (sum(tabUp[c(5)]) + sum(tabDo[c(5)])) / (sum(tabUp[c(2,3,4,5)]) + sum(tabDo[c(2,3,4,5)]))




tUp <- rev(cumsum(rev(tabUp))) / nrow(mUp)
tDo <- rev(cumsum(rev(tabDo))) / nrow(mDo)

tUp <- tUp[-1] / tUp[2]
tDo <- tDo[-1] / tDo[2]


cm <- cor(assays(diffATAC)[["log2FoldChange"]])

pdf("~/Desktop/Human/diff-cor.pdf", width = 7, height = 3.75)
cm <- cor(assays(diffATAC)[["log2FoldChange"]][,grep("LKB1", colnames(diffATAC))])
diag(cm) <- NA
cm[cm < 0] <- 0
cm[cm > 0.2] <- 0.2
pheatmap::pheatmap(cm, color = paletteContinuous(), border_color = "black", 
	show_colnames = F)
dev.off()


df2 <- data.frame(x=c(tUp[-1], -tDo[-1]))
df2$N <- paste0("N",rep(seq_len(4), 2))
df2$Var1 <- "Up"
df2$Var1[df2$x < 0] <- "Do"

pdf("~/Desktop/Human/Relative-ATAC-Differences-Consistent.pdf", width = 4, height = 5)
ggplot(df2, aes(N, x, fill=Var1)) + 
	geom_bar(stat = "identity") + 
	ArchR:::theme_ArchR(xText90 = TRUE) + 
	scale_fill_manual(values = paletteDiscrete(values=df2$Var1)) +
	ylab("Proportion of Differential Peaks") + 
	scale_y_continuous(limits = c(-0.125, 0.175), breaks = seq(-10,10,0.05))
dev.off()






mUp <- assays(diffATAC)[["log2FoldChange"]] >= 0.5 & assays(diffATAC)[["padj"]] <= 0.05
mDo <- assays(diffATAC)[["log2FoldChange"]] <= -0.5 & assays(diffATAC)[["padj"]] <= 0.05

tabUp <- table(rowSums(mUp[,c("H1437_KEAP1_GFP", "A549_KEAP1_GFP", "H460_KEAP1_GFP", "H1355_KEAP1_GFP")]))
tabDo <- table(rowSums(mDo[,c("H1437_KEAP1_GFP", "A549_KEAP1_GFP", "H460_KEAP1_GFP", "H1355_KEAP1_GFP")]))

tUp <- rev(cumsum(rev(tabUp))) / nrow(mUp)
tDo <- rev(cumsum(rev(tabDo))) / nrow(mDo)



m <- assay(seCV)[,grep("LKB1", colnames(seCV), value = TRUE)]






motifs <- c("Sox17", "Nkx21", "Tead3", "Mef2", "Foxa", "Zeb1", "Runx", "Nr4")
motifs2 <- lapply(seq_along(motifs), function(x){
	grep(motifs[x], rownames(seCV), value = TRUE, ignore.case = TRUE)
}) %>% unlist


for(i in motifs2){

	print(i)

	motif <- i

	df <- reshape2::melt(assays(seCV[motif, ])[[2]])
	df <- df[!grepl("KEAP1", df$Var2),]
	df$Group <- gsub("_Rep2","",gsub("_Rep1","",paste0(df$Var2)))
	df$Group2 <- gsub("GFP|LKB1","",paste0(df$Group))
	df$Type <- "LKB1"
	df$Type[grep("GFP",df$Group)] <- "GFP"

	df2 <- lapply(split(df, df$Group2), function(x){
		x$value <- x$value - mean(x$value[grep("GFP", x$Var2)])
		x
	}) %>% Reduce("rbind",.) %>% data.frame

	df3 <- lapply(split(df2, df2$Group), function(x){
		x$value <- mean(x$value)
		x[1,]
	}) %>% Reduce("rbind",.) %>% data.frame

	pdf(paste0("~/Desktop/Human/ChromVAR-",motif,"-Human.pdf"), width = 5, height = 5, useDingbats = FALSE)
	p <- ggplot(df[grep("GFP", df$Var2),], aes(Group, value, color = Type)) + 
			ArchR:::theme_ArchR(xText90 = TRUE) + 
			geom_hline(yintercept = 0, size = 0.25) + 
			#geom_bar(data = df3, aes(x = Group, y = value, fill = Type), stat = "identity") + 
			geom_point(color="black", size = 1) +
			scale_color_manual(values=c("GFP" = "firebrick2", "LKB1" = "dodgerblue3")) +
			scale_fill_manual(values=c("GFP" = "firebrick2", "LKB1" = "dodgerblue3")) +
			guides(color = FALSE) + 
			ylab(paste0("chromVAR Deviations ", motif))
	print(ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5, height=0.5))
	grid::grid.newpage()
	p <- ggplot(df2, aes(Group, value, color = Type)) + 
			ArchR:::theme_ArchR(xText90 = TRUE) + 
			geom_hline(yintercept = 0, size = 0.25) + 
			geom_bar(data = df3, aes(x = Group, y = value, fill = Type), stat = "identity") + 
			geom_point(color="black", size = 1) +
			scale_color_manual(values=c("GFP" = "firebrick2", "LKB1" = "dodgerblue3")) +
			scale_fill_manual(values=c("GFP" = "firebrick2", "LKB1" = "dodgerblue3")) +
			guides(color = FALSE) + 
			ylab(paste0("chromVAR Deviations ", motif))
	print(ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5, height=0.5))
	dev.off()


}






motifs <- c("Sox17", "Nkx21", "Tead3", "Mef2", "Foxa", "Zeb1", "Runx", "Nr4")
motifs2 <- lapply(seq_along(motifs), function(x){
	grep(motifs[x], rownames(seCV), value = TRUE, ignore.case = TRUE)
}) %>% unlist


for(i in motifs2){

	print(i)

	motif <- i

	df <- reshape2::melt(assays(seCV[motif, ])[[2]])
	df <- df[!grepl("LKB1", df$Var2),]
	df$Group <- gsub("_Rep2","",gsub("_Rep1","",paste0(df$Var2)))
	df$Group2 <- gsub("GFP|KEAP1","",paste0(df$Group))
	df$Type <- "KEAP1"
	df$Type[grep("GFP",df$Group)] <- "GFP"

	df2 <- lapply(split(df, df$Group2), function(x){
		x$value <- x$value - mean(x$value[grep("GFP", x$Var2)])
		x
	}) %>% Reduce("rbind",.) %>% data.frame

	df3 <- lapply(split(df2, df2$Group), function(x){
		x$value <- mean(x$value)
		x[1,]
	}) %>% Reduce("rbind",.) %>% data.frame

	pdf(paste0("~/Desktop/Human/ChromVAR-",motif,"-Human-KEAP1.pdf"), width = 5, height = 5, useDingbats = FALSE)
	p <- ggplot(df[grep("GFP", df$Var2),], aes(Group, value, color = Type)) + 
			ArchR:::theme_ArchR(xText90 = TRUE) + 
			geom_hline(yintercept = 0, size = 0.25) + 
			#geom_bar(data = df3, aes(x = Group, y = value, fill = Type), stat = "identity") + 
			geom_point(color="black", size = 1) +
			scale_color_manual(values=c("GFP" = "firebrick2", "LKB1" = "dodgerblue3")) +
			scale_fill_manual(values=c("GFP" = "firebrick2", "LKB1" = "dodgerblue3")) +
			guides(color = FALSE) + 
			ylab(paste0("chromVAR Deviations ", motif))
	print(ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5, height=0.5))
	grid::grid.newpage()
	p <- ggplot(df2, aes(Group, value, color = Type)) + 
			ArchR:::theme_ArchR(xText90 = TRUE) + 
			geom_hline(yintercept = 0, size = 0.25) + 
			geom_bar(data = df3, aes(x = Group, y = value, fill = Type), stat = "identity") + 
			geom_point(color="black", size = 1) +
			scale_color_manual(values=c("GFP" = "firebrick2", "KEAP1" = "dodgerblue3")) +
			scale_fill_manual(values=c("GFP" = "firebrick2", "KEAP1" = "dodgerblue3")) +
			guides(color = FALSE) + 
			ylab(paste0("chromVAR Deviations ", motif))
	print(ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5, height=0.5))
	dev.off()


}


o <- c("H1650_LKB1_GFP", "H1975_LKB1_GFP", "H358_LKB1_GFP", "H2009_LKB1_GFP", "H1437_LKB1_GFP", "A549_LKB1_GFP", "H460_LKB1_GFP", "H1355_LKB1_GFP")
o <- gsub("_LKB1_GFP", "", o)


unique(unlist(lapply(o, function(x) grep(x, paste0(df2$Group), value = TRUE))))

paste0(df2$Group)


motifs <- c("Sry", "Sox17", "Nkx21", "Tead3", "Mef2", "Foxa", "Zeb1", "Runx", "Nr4")
motifs2 <- lapply(seq_along(motifs), function(x){
	grep(motifs[x], rownames(seCV), value = TRUE, ignore.case = TRUE)
}) %>% unlist


for(i in motifs2){

	print(i)

	motif <- i

	df <- reshape2::melt(assays(seCV[motif, ])[[2]])
	df$Group <- gsub("_Rep2","",gsub("_Rep1","",paste0(df$Var2)))
	df$Group2 <- gsub("GFP|KEAP1|LKB1","",paste0(df$Group))
	df$Type <- "KEAP1"
	df$Type[grep("GFP",df$Group)] <- "GFP"
	df$Type[grep("LKB1",df$Group)] <- "LKB1"
	df$Type3 <- "2Deficient"
	df$Type3[grep("1650|1975|358|2009",df$Group)] <- "1Proficient"

	df2 <- lapply(split(df, df$Group2), function(x){
		x$value <- x$value - mean(x$value[grep("GFP", x$Var2)])
		x
	}) %>% Reduce("rbind",.) %>% data.frame

	df2$Group <- factor(paste0(df2$Group), levels = unique(unlist(lapply(o, function(x) grep(x, paste0(df2$Group), value = TRUE)))))

	df3 <- lapply(split(df2, df2$Group), function(x){
		x$value <- mean(x$value)
		x[1,]
	}) %>% Reduce("rbind",.) %>% data.frame

	plotDF <- df[grep("GFP", df$Var2),]

	pdf(paste0("~/Desktop/Human/ChromVAR-",motif,"-Human-All.pdf"), width = 5, height = 5, useDingbats = FALSE)
	p <- ggplot(data = plotDF, aes(Type3, value, color = Type3, fill = Type3)) +
		geom_violin(alpha=0.5) + 
    		theme_ArchR() +
    		geom_point() +
    		geom_boxplot(size = 0.5, outlier.size = 0, 
              outlier.stroke = 0, fill = NA) + 
	    ylab("ChromVAR Deviation Scores") +
	    scale_fill_manual(values = c("1Proficient" = "dodgerblue4", "2Deficient" = "firebrick3")) +
	    scale_color_manual(values = c("1Proficient" = "dodgerblue4", "2Deficient" = "firebrick3")) + 
	    ggtitle(paste0(motif,
	      "\nPvalueKS = ", ks.test(plotDF$value[plotDF$Type3=="1Proficient"], plotDF$value[plotDF$Type3!="1Proficient"])$p.value,
	      "\nPvalueTT = ", t.test(plotDF$value[plotDF$Type3=="1Proficient"], plotDF$value[plotDF$Type3!="1Proficient"])$p.value
	      )
	    )

	print(ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5, height=0.5))
	grid::grid.newpage()
	p <- ggplot(df2, aes(Group, value, color = Type)) + 
			ArchR:::theme_ArchR(xText90 = TRUE) + 
			geom_hline(yintercept = 0, size = 0.25) + 
			geom_bar(data = df3, aes(x = Group, y = value, fill = Type), stat = "identity") + 
			geom_point(color="black", size = 1) +
			scale_color_manual(values=c("GFP" = "firebrick2", "LKB1" = "dodgerblue3", "KEAP1" = "forestgreen")) +
			scale_fill_manual(values=c("GFP" = "firebrick2", "LKB1" = "dodgerblue3", "KEAP1" = "forestgreen")) +
			guides(color = FALSE) + 
			ylab(paste0("chromVAR Deviations ", motif))
	print(ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5, height=0.5))
	grid::grid.newpage()
	p <- ggplot(df2[!grepl("KEAP1",df2$Group),], aes(Group, value, color = Type)) + 
			ArchR:::theme_ArchR(xText90 = TRUE) + 
			geom_hline(yintercept = 0, size = 0.25) + 
			geom_bar(data = df3[!grepl("KEAP1",df3$Group),], aes(x = Group, y = value, fill = Type), stat = "identity") + 
			geom_point(color="black", size = 1) +
			scale_color_manual(values=c("GFP" = "firebrick2", "LKB1" = "dodgerblue3", "KEAP1" = "forestgreen")) +
			scale_fill_manual(values=c("GFP" = "firebrick2", "LKB1" = "dodgerblue3", "KEAP1" = "forestgreen")) +
			guides(color = FALSE) + 
			ylab(paste0("chromVAR Deviations ", motif))
	print(ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5, height=0.5))

	dev.off()


}



###################
# Restore MA Plots
###################

ggPlotDifferential <- function(pse = NULL, log2Mean = NULL, 
  log2FoldChange = NULL, padj = NULL, plot = "MA", log2FCCutOff = 1, xlim = NULL, ylim = NULL,
  padjCutOff = 0.1, limLog2mean = c(0.015, 0.99999), limLog2FC = 0.99999, limPadj = 0.99999, n = NULL, size = 0.5, alpha = 0.75, rastr = TRUE, 
  pal = c("Increasing"="firebrick3","Decreasing"="dodgerblue3","Constitutive"="lightgrey"), ...){

  suppressPackageStartupMessages(require(ggplot2))

  #Create Plot DF
  if(is.null(pse)){
    if(is.null(log2Mean) | is.null(log2FoldChange)){
      stop("need both log2Mean and log2FoldChange vectors!")
    }
    if(is.null(padj)){
      padj <- rep(1,length(log2FoldChange))
    }
  df <- data.frame(log2Mean = log2Mean, log2FoldChange = log2FoldChange, padj = padj)   
  }else{
    df <- data.frame(assay(pse))
    if(!all(c("log2Mean", "log2FoldChange", "padj") %in% colnames(df))){
      stop("Not pairwiseTest Summarized Experiment! Need log2Mean log2FC and padj column names in assay(pse)!")
    }
    df <- data.frame(df)[,c("log2Mean", "log2FoldChange", "padj")]
  }

  #Filter out NAs
  df <- df[complete.cases(df),]
  df$padj[df$padj < 10^-100] <- 10^-100

  #Define Differential Peaks
  df$class <- "Constitutive"
  df$class[which( df[,"log2FoldChange"] >= log2FCCutOff & df[,"padj"] <= padjCutOff ) ] <- "Increasing"
  df$class[which( df[,"log2FoldChange"] <= -log2FCCutOff & df[,"padj"] <= padjCutOff ) ] <- "Decreasing"

  #Subset so less overlapping
  if(!is.null(n)){
  df <- nSample(df,n,"r")
  }

  #Order Points
  df <- df[order(df$padj,decreasing=TRUE),]

  #Calc limits
  if(!is.null(limLog2mean)){
    limLog2mean2 <- quantile(df$log2Mean, limLog2mean)
  }

  if(!is.null(limLog2FC)){
    limLog2FC2 <- c(-quantile(abs(df$log2FoldChange), limLog2FC),quantile(abs(df$log2FoldChange), limLog2FC))
  }

  if(!is.null(limPadj)){
    limPadj2 <- c(0, quantile(-log10(df$padj), limPadj))
  }

  #Color
  df$color <- df$class
  colors <- pal

  if(!is.null(log2FCCutOff) & log2FCCutOff != 0){
    lfcint <- c(-log2FCCutOff, 0, log2FCCutOff)
  }else{
    lfcint <- 0
  }

  up <- length(which(paste0(df$class)=="Increasing"))
  do <- length(which(paste0(df$class)=="Decreasing"))

  title <- sprintf("%s total features tested\n%s (%s) significant features up\n%s (%s) significant features down", 
  	nrow(df), up, round(up/nrow(df), 3), do, round(do/nrow(df), 3))
  if(tolower(plot) == "ma"){
    
  	if(!is.null(xlim)){
  		limLog2mean2 <- xlim
  	}

   	if(!is.null(ylim)){
  		limLog2FC2 <- ylim
  	}

    out <- ArchR::ggPoint(df$log2Mean, df$log2FoldChange, df$color, pal = colors, xlabel = "Log2 Mean", ylabel = "Log2 Fold Change",
      size = size, alpha = alpha, title = title, xlim = limLog2mean2, ylim = limLog2FC2, rastr = rastr,
      labelMeans = FALSE,  ...) +
      geom_hline(yintercept = lfcint, lty = "dashed", size = 0.25) + scale_y_continuous(limits = limLog2FC2, breaks = seq(-50,50,2), expand = FALSE)
  
  }else if(tolower(plot) == "volcano"){
   
    out <- ggPoint(df$log2FoldChange, -log10(df$padj), df$color, pal = colors, xlabel = "Log2 Fold Change", ylabel = "-Log10 Pvalue",
      size = size, alpha = alpha, title = title, xlim = limLog2FC2, ylim = limPadj2, rastr = rastr, ...) +
      geom_vline(xintercept = lfcint, lty = "dashed") + scale_x_continuous(limits = limLog2FC2, breaks = seq(-50,50,2), expand = FALSE)
  
  }else if(tolower(plot) == "onetoone"){
  }

  out <- out + theme(legend.position = "bottom")

  return(out)
}



devMouse <- readRDS("/Volumes/JG_SSD_1/Data/LKB1/Initial_2695_2804_2381_2565/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")

mouseT2 <- grep("EtOH", grep("2695|2804", colnames(devMouse),value=TRUE), value = TRUE)
mouseT1 <- grep("TAM", grep("2695|2804", colnames(devMouse),value=TRUE), value = TRUE)

type1 <- grep("1650|1975|358|2009", colnames(seCV), value = TRUE)
type2 <- colnames(seCV)[!grepl("1650|1975|358|2009", colnames(seCV))]

diffMouse <- rowMeans(assays(devMouse)[[2]][,mouseT2]) - rowMeans(assays(devMouse)[[2]][,mouseT1])
diffHuman <- rowMeans(assays(seCV)[[2]][,type2]) - rowMeans(assays(seCV)[[2]][,type1])

TFHuman <- tolower(stringr::str_split(rownames(seCV), pattern = "\\_", simplify=TRUE)[,1])
TFMouse <- tolower(stringr::str_split(rownames(devMouse), pattern = "\\_", simplify=TRUE)[,1])

common <- intersect(tolower(TFHuman), tolower(TFMouse))

matchTF <- lapply(seq_along(common), function(x){

  TFHumanx <- rownames(seCV)[which(TFHuman %in% common[x])]
  TFMousex <- rownames(devMouse)[which(TFMouse %in% common[x])]

  df <- expand.grid(Human = TFHumanx, Mouse = TFMousex)
  df$Common <- common[x]
  df$N <- nrow(df)
  df

}) %>% Reduce("rbind", .)


colorList <- list(
  "NKX" = c("NKX"),
  "NFI" = c("NFI"),
  "MEF2" = c("MEF2"),
  "TEAD" = c("TEAD"),
  "RUNX" = c("RUNX"),
  "SOX" = c("SOX", "SRY"),
  "MESP" = c("MESP"),
  "ZEB" = c("ZEB", "SNAI"),
  "FOXA" = c("FOXA"),
  "NR4A" = c("NR4A"),
  "BZIP" = c("JUN","FOS","BATF")
)

plotDF <- data.frame(
  x = -diffHuman[paste0(matchTF[,1])],
  y = -diffMouse[paste0(matchTF[,2])],
  common = matchTF[,3],
  human = matchTF[,1],
  mouse = matchTF[,2],
  N = matchTF$N
)
plotDF$Sum <- abs(plotDF$x) + abs(plotDF$y)
plotDF <- plotDF[order(plotDF$Sum, decreasing= TRUE), ]
plotDF <- plotDF[!duplicated(plotDF$common),]
plotDF <- plotDF[plotDF$N < 10, ]


pal <- paletteDiscrete(values = names(colorList))
pal["NA"] <- "grey"

plotDF$Family <- "NA"

for(i in seq_along(colorList)){
  idx <- unlist(lapply(colorList[[i]], function(x)grep(x,plotDF$common, value = TRUE, ignore.case = TRUE)))
  plotDF$Family[which(plotDF$common %in% idx)] <- names(colorList)[i]
}

plotDF[-which(plotDF$Sum > 10 & abs(plotDF$x) > 1 & abs(plotDF$y) > 1), "Family"] <- "NA"
plotDF <- rbind(plotDF[plotDF$Family=="NA",], plotDF[plotDF$Family!="NA",])

pdf("~/Desktop/Plot-HumanCellLines-Mouse-Human.pdf", width = 5, height = 5, useDingbats = FALSE)

ggPoint(
  x = plotDF$x,
  y = plotDF$y,
  plotDF$Family,
  pal =pal,
  labelAsFactors = FALSE,
  labelMeans = FALSE,
  size = 2
) + geom_hline(yintercept = 0, lty = "dashed") +
geom_vline(xintercept = 0, lty = "dashed") +
#ggrepel::geom_label_repel(data = labelDF, aes(label = common), color = "black") +
xlab("ChromVAR Human Type1 - Type2") +
ylab("ChromVAR Mouse Type1 - Type2")

dev.off()













































































































































































plotPCA$Group <- gsub("_Rep2","",gsub("_Rep1","",gsub("_2","",gsub("_1","",gsub("S2695_TKO_","",gsub("_Tam","",gsub("_EtOH","",rownames(plotPCA))))))))
plotPCA$Type <- "EtOH"
plotPCA$Type[grep("Tam", rownames(plotPCA))] <- "Tam"

pal <- c(
	"sgSafe" = "#FDE400",
	"sgLkb1" = "#F27D2A",
	"sgPrkaa" = "#208A42",
	"sgSik" = "#BF6CAB",
	"sgNuak" = "#D8A767",
	"sgMark" = "#8A9FD1"
)

pal <- paletteDiscrete(values = plotPCA$Group)

pVar <- head(round(100*pcaVar$sdev^2/sum(pcaVar$sdev^2),3),15)
df <- data.frame(pc=paste0("PC",seq_along(pVar)), var=pVar)
df$pc <- factor(paste0("PC",seq_along(pVar)), levels = paste0("PC",seq_along(pVar)))
pdf("~/Desktop/sgSik2/PCA-PercentVar-Set1.pdf", width = 4, height = 4, useDingbats = FALSE)
p <- ggplot(df, aes(pc,var)) + 
	geom_bar(stat="identity", fill="dodgerblue4") +
	ArchR::theme_ArchR(xText90=TRUE) + 
	ylab("Percent Variance Explained") +
	scale_y_continuous(limits = c(0,max(pVar)*1.05),expand = c(0,0))
ArchR:::.fixPlotSize(p, plotWidth = 4, plotHeight = 4, height=0.5)
dev.off()

pdf("~/Desktop/sgSik2/PCA-10k-Variable-Peaks-Set1.pdf", width = 4, height = 4, useDingbats = FALSE)
p <- ggplot(plotPCA, aes(PC1,PC2,color=Group)) + 
	geom_point(size=1.5) + theme_ArchR() + 
	xlab(pVar[1]) + ylab(pVar[2]) +
	scale_color_manual(values = pal)
ArchR:::.fixPlotSize(p, plotWidth = 4, plotHeight = 4, height=1)
dev.off()

PC1Safe <- mean(plotPCA[plotPCA$Group=="S2695_sgSafe" & plotPCA$Type=="EtOH",1])

 [1] "S2695_sgLkb1.EtOH"   "S2695_sgLkb1.Tam"    "S2695_sgSafe.EtOH"  
 [4] "S2695_sgSafe.Tam"    "S2695_sgSik1.EtOH"   "S2695_sgSik1.Tam"   
 [7] "S2695_sgSik1_3.EtOH" "S2695_sgSik1_3.Tam"  "S2695_sgSik2.EtOH"  
[10] "S2695_sgSik2.Tam"    "S2695_sgSik3.EtOH"   "S2695_sgSik3.Tam" 

o <-  c("S2695_sgSafe.EtOH", "S2695_sgLkb1.EtOH", "S2695_sgSik1_3.EtOH", "S2695_sgSik1.EtOH", "S2695_sgSik2.EtOH", "S2695_sgSik3.EtOH",
	"S2695_sgSafe.Tam", "S2695_sgLkb1.Tam", "S2695_sgSik1_3.Tam", "S2695_sgSik1.Tam", "S2695_sgSik2.Tam", "S2695_sgSik3.Tam")

plotPCA$GroupType <- factor(paste0(plotPCA$Group,".",plotPCA$Type), levels = o)

plotPCA2 <- lapply(split(plotPCA, plotPCA$GroupType), function(x){
	x$PC1[1] <- mean(x$PC1)
	x[1,]
}) %>% Reduce("rbind",.) %>% data.frame

pdf("~/Desktop/sgSik2/PCA-Bar-10k-Variable-Peaks-Set1.pdf", width = 4, height = 4, useDingbats = FALSE)
p <- ggplot(plotPCA, aes(x = GroupType, y = PC1 - PC1Safe, fill = Type)) + 
	geom_hline(yintercept = 0, size = 0.25) + 
	geom_bar(data = plotPCA2, aes(x = GroupType, y = PC1 - PC1Safe, fill = Type), stat = "identity") + 
	geom_point(size = 0.5) + 
	scale_fill_manual(values = c("EtOH" = "firebrick2", "Tam" = "dodgerblue3")) +
	ylab(paste0("Delta PC1 ", pVar[1])) + 
	ArchR::theme_ArchR(xText90=TRUE)
ArchR:::.fixPlotSize(p, plotWidth = 4, plotHeight = 4, height=0.75)
dev.off()



set <- 1
set.seed(1)
gM <- groupMeans(assay(se), colData(se)$Group)
rV <- matrixStats::rowVars(gM)
varPeaks <- head(order(rV,decreasing=TRUE), 10000)

ceiling <- 2
set.seed(1)

mat1 <- assays(se[varPeaks,])[["log2norm"]]
mat1 <- groupMeans(mat1, colData(se)[,"Group"])
mat1 <- rowZscores(mat1) * 2/1.5
mat1[mat1 > ceiling] <- ceiling
mat1[mat1 < -ceiling] <- -ceiling

j <- 5
km <- kmeans(mat1, j, 100, 10)

colData <- DataFrame(
	row.names = colnames(mat1),  
	Group = gsub("_2","",gsub("_1","",gsub("S2695_","",gsub("_Tam","",gsub("_EtOH","",colnames(mat1)))))),
	Treat = "EtOH"
	)
colData$Treat[grep("Tam", colnames(mat1))] <- "Tam"


cm <- colorMapAnno(colData)
cm[[1]][gsub("S2695_","", names(pal))] <- pal


g <- grep("EtOH",rownames(colData),value=TRUE)
o <- c("S2695_sgSafe_1_EtOH", "S2695_sgLkb1_1_EtOH", "S2695_sgSik1_3_1_EtOH", "S2695_sgSik1_1_EtOH", "S2695_sgSik2_1_EtOH", "S2695_sgSik3_1_EtOH")


o <- c(
	paste0("S2695_TKO_sgSafe_",set,"_EtOH"),
	paste0("S2695_TKO_sgLkb1_",set,"_EtOH"),
	paste0("S2695_TKO_sgSik_",set,"_EtOH"),
	paste0("S2695_TKO_sgPrkaa_",set,"_EtOH"),
	paste0("S2695_TKO_sgMark_",set,"_EtOH"),
	paste0("S2695_TKO_sgNuak_",set,"_EtOH")
)

#Heatmap
pdf(file.path(outDiri,paste0("Heatmap-Var-Peaks-Set",set,".pdf")), width = 8, height = 12, useDingbats=FALSE)

ht1 <- simpleHeatmap(
	mat1[,o],
	scale = FALSE,
	split = km$cluster,
	colData = colData[o,],
	clusterRows = TRUE,
	clusterCols = FALSE,
	colorMap = cm,
	labelCols = TRUE,
	showColDendrogram = TRUE,
	draw = FALSE
)

o2 <- gsub("EtOH","Tam",o)

ht2 <- simpleHeatmap(
	mat1[,o2],
	scale = FALSE,
	split = km$cluster,
	colData = colData[o2,],
	clusterRows = TRUE,
	clusterCols = FALSE,
	colorMap = cm,
	labelCols = TRUE,
	showColDendrogram = TRUE,
	draw = FALSE
)

draw(ht1 + ht2, padding = unit(c(45, 45, 45, 
            45), "mm"), heatmap_legend_side = "bot", annotation_legend_side = "bot")

dev.off()


#Compute Motifs
featureList <- lapply(seq_along(unique(km$cluster)), function(x){
	featureX <- featureEnrichment(tfMatches,varPeaks[which(km$cluster==unique(km$cluster)[x])], seq_len(nrow(tfMatches)))
	readr::write_tsv(featureX, file.path(outDiri,paste0("K",j,"-Kmeans_",unique(km$cluster)[x],paste0("-Motifs-Enriched-Set",set,".tsv"))))
})

#Compute Gene Matches
featureList <- lapply(seq_along(unique(km$cluster)), function(x){
	featureX <- featureEnrichment(geneMatches,varPeaks[which(km$cluster==unique(km$cluster)[x])], seq_len(nrow(geneMatches)))
	readr::write_tsv(featureX, file.path(outDiri,paste0("K",j,"-Kmeans_",unique(km$cluster)[x],paste0("-Genes-Enriched-Set",set,".tsv"))))
})

peaks <- rownames(mat1) %>% stringr::str_split(pattern="_",simplify=TRUE) %>% {GRanges(.[,1],IRanges(as.integer(.[,2]),as.integer(.[,3])))}
peaks$Kmeans <- km$cluster
saveRDS(peaks, paste0("~/Desktop/sgSik/K5-Kmeans-PeakRanges_",set,".rds"))



#######################################################
#######################################################

motifs <- c("Sox17_843", "Nkx21_446", "Tead3_1304", "Mef2a_717", "Foxa1_348", "Zeb1_162")
o <-  c("S2695_sgSafe_EtOH", "S2695_sgLkb1_EtOH", "S2695_sgSik1_3_EtOH", "S2695_sgSik1_EtOH", "S2695_sgSik2_EtOH", "S2695_sgSik3_EtOH",
	"S2695_sgSafe_Tam", "S2695_sgLkb1_Tam", "S2695_sgSik1_3_Tam", "S2695_sgSik1_Tam", "S2695_sgSik2_Tam", "S2695_sgSik3_Tam")

for(i in motifs){

	print(i)

	motif <- i

	df <- reshape2::melt(assays(seCV[motif, ])[[2]])
	df$Group <- gsub("_2_","_",gsub("_Rep2","",gsub("_Rep1","",paste0(df$Var2))))
	df$Type <- "EtOH"
	df$Type[grep("Tam",df$Group)] <- "Tam"
	#df$Type <- apply(stringr::str_split(df$Group,pattern="\\_",simplify=TRUE)[,2:3],1,function(x)paste(x,collapse="_"))
	#df$Type <- factor(df$Type, c("KP_Primary","KPL_Primary","KP_Met","KPL_Met"))

	set <- 2

	#df <- df[grep("Safe|Lkb1|Sik",df$Group,ignore.case = T),]
	df$Group <- factor(paste0(df$Group), levels = o)

	val1 <- df[intersect(grep("Safe",df$Group), grep("EtOH",df$Group)),"value"]
	val2 <- df[intersect(grep("Safe",df$Group), grep("Tam",df$Group)),"value"]

	df2 <- lapply(split(df, df$Group), function(x){
		x$value[1] <- mean(x$value)
		x[1,]
	}) %>% Reduce("rbind",.) %>% data.frame

	valueET <- df2[df2$Group=="S2695_sgSafe_EtOH","value"]

	pdf(paste0("~/Desktop/sgSik2/ChromVAR-",motif,"-TKO.pdf"), width = 5, height = 5, useDingbats = FALSE)
	p <- ggplot(df, aes(Group, value - valueET, color = Type)) + 
			ArchR:::theme_ArchR(xText90 = TRUE) + 
			geom_hline(yintercept = 0, size = 0.25) + 
			geom_bar(data = df2, aes(x = Group, y = value - valueET, fill = Type), stat = "identity") + 
			geom_point(color="black", size = 1) +
			scale_color_manual(values=c("EtOH" = "firebrick2", "Tam" = "dodgerblue3")) +
			scale_fill_manual(values=c("EtOH" = "firebrick2", "Tam" = "dodgerblue3")) +
			guides(color = FALSE) + 
			ylab(paste0("chromVAR Deviations ", motif))
	print(ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5, height=0.5))
	dev.off()


}

#######################################################
#######################################################


pdf("~/Desktop/Correlation-ATAC-sgSik-Set2.pdf", width = 7, height = 3.75)
cm <- cor(assay(se))
diag(cm) <- NA
pheatmap::pheatmap(cm, color = paletteContinuous(), border_color = "black", 
	show_colnames = F)
dev.off()





seDiff <- readRDS("output/ATAC_Analysis/Main/EdgeR-Pairwise-Summary.rds")

#########################################
#########################################
dfChange <- data.frame(
	peaksUp = colSums(assays(seDiff)[["log2FoldChange"]] >= 0.5 & assays(seDiff)[["padj"]] <= 0.05),
	peaksDo = colSums(assays(seDiff)[["log2FoldChange"]] <= -0.5 & assays(seDiff)[["padj"]] <= 0.05)
)
dfChange$peaksUpP <- round(100 * dfChange$peaksUp / nrow(seDiff), 2)
dfChange$peaksDoP <- round(100 * dfChange$peaksDo / nrow(seDiff), 2)
dfChange$name <- rownames(dfChange)

idx <- grep("Restore", rownames(dfChange), value = TRUE)

dfPlot <- dfChange[idx, ]
dfPlot[, 4] <- -dfPlot[,4]
dfPlot <- dfPlot[-nrow(dfPlot), ]

df2 <- data.frame(p = c(dfPlot[,3], dfPlot[,4]), name = c(dfPlot$name, dfPlot$name))
df2$type <- ifelse(df2$p > 0, "Increasing", "Decreasing")
df2$name <- factor(df2$name, levels = c("Restore_sgSafe", "Restore_sgLkb1", "Restore_sgSik1_3","Restore_Sik1","Restore_sgSik2","Restore_sgSik3"))

pdf("~/Desktop/sgSik/Test-Change-Sik.pdf", width = 3, height = 4)
ggplot(df2, aes(name, p, fill = type)) +
	geom_bar(stat="identity", position="identity") +
	ArchR::theme_ArchR(xText90 = TRUE) +
	scale_y_continuous(limits=c(-12,8),breaks=seq(-12,12,4),expand=c(0,0)) + 
	ylab("Percent Differential Peaks") +
	scale_fill_manual(values=c("Increasing"="firebrick3", "Decreasing" = "dodgerblue4"))
dev.off()





#########################################
#########################################

pdf("~/Desktop/sgSik/Differential-sgSik-1-3.pdf", width = 4, height = 4)
ggPlotDifferential(
	log2Mean = assays(seDiff[,"Restore_sgSik1_3"])[["log2Mean"]][,1],
	log2FoldChange = assays(seDiff[,"Restore_sgSik1_3"])[["log2FoldChange"]][,1],
	padj = assays(seDiff[,"Restore_sgSik1_3"])[["padj"]][,1],
	padjCutOff = 0.05,
	log2FCCutOff = 0.5,
	ratioYX = 3.5/4,
	ylim = c(-4, 4),
	pal = c("Increasing"="dodgerblue3","Decreasing"="firebrick3","Constitutive"="lightgrey")) %>% ArchR:::.fixPlotSize(plotWidth = 3.5, plotHeight=4)
dev.off()

pdf("~/Desktop/sgSik/Differential-sgSik-1.pdf", width = 4, height = 4)
ggPlotDifferential(
	log2Mean = assays(seDiff[,"Restore_Sik1"])[["log2Mean"]][,1],
	log2FoldChange = assays(seDiff[,"Restore_Sik1"])[["log2FoldChange"]][,1],
	padj = assays(seDiff[,"Restore_Sik1"])[["padj"]][,1],
	padjCutOff = 0.05,
	log2FCCutOff = 0.5,
	ratioYX = 3.5/4,
	ylim = c(-4, 4),
	pal = c("Increasing"="dodgerblue3","Decreasing"="firebrick3","Constitutive"="lightgrey")) %>% ArchR:::.fixPlotSize(plotWidth = 3.5, plotHeight=4)
dev.off()


###################
# Restore MA Plots
###################

ggPlotDifferential <- function(pse = NULL, log2Mean = NULL, 
  log2FoldChange = NULL, padj = NULL, plot = "MA", log2FCCutOff = 1, xlim = NULL, ylim = NULL,
  padjCutOff = 0.1, limLog2mean = c(0.015, 0.99999), limLog2FC = 0.99999, limPadj = 0.99999, n = NULL, size = 0.5, alpha = 0.75, rastr = TRUE, 
  pal = c("Increasing"="firebrick3","Decreasing"="dodgerblue3","Constitutive"="lightgrey"), ...){

  suppressPackageStartupMessages(require(ggplot2))

  #Create Plot DF
  if(is.null(pse)){
    if(is.null(log2Mean) | is.null(log2FoldChange)){
      stop("need both log2Mean and log2FoldChange vectors!")
    }
    if(is.null(padj)){
      padj <- rep(1,length(log2FoldChange))
    }
  df <- data.frame(log2Mean = log2Mean, log2FoldChange = log2FoldChange, padj = padj)   
  }else{
    df <- data.frame(assay(pse))
    if(!all(c("log2Mean", "log2FoldChange", "padj") %in% colnames(df))){
      stop("Not pairwiseTest Summarized Experiment! Need log2Mean log2FC and padj column names in assay(pse)!")
    }
    df <- data.frame(df)[,c("log2Mean", "log2FoldChange", "padj")]
  }

  #Filter out NAs
  df <- df[complete.cases(df),]
  df$padj[df$padj < 10^-100] <- 10^-100

  #Define Differential Peaks
  df$class <- "Constitutive"
  df$class[which( df[,"log2FoldChange"] >= log2FCCutOff & df[,"padj"] <= padjCutOff ) ] <- "Increasing"
  df$class[which( df[,"log2FoldChange"] <= -log2FCCutOff & df[,"padj"] <= padjCutOff ) ] <- "Decreasing"

  #Subset so less overlapping
  if(!is.null(n)){
  df <- nSample(df,n,"r")
  }

  #Order Points
  df <- df[order(df$padj,decreasing=TRUE),]

  #Calc limits
  if(!is.null(limLog2mean)){
    limLog2mean2 <- quantile(df$log2Mean, limLog2mean)
  }

  if(!is.null(limLog2FC)){
    limLog2FC2 <- c(-quantile(abs(df$log2FoldChange), limLog2FC),quantile(abs(df$log2FoldChange), limLog2FC))
  }

  if(!is.null(limPadj)){
    limPadj2 <- c(0, quantile(-log10(df$padj), limPadj))
  }

  #Color
  df$color <- df$class
  colors <- pal

  if(!is.null(log2FCCutOff) & log2FCCutOff != 0){
    lfcint <- c(-log2FCCutOff, 0, log2FCCutOff)
  }else{
    lfcint <- 0
  }

  up <- length(which(paste0(df$class)=="Increasing"))
  do <- length(which(paste0(df$class)=="Decreasing"))

  title <- sprintf("%s total features tested\n%s (%s) significant features up\n%s (%s) significant features down", 
  	nrow(df), up, round(up/nrow(df), 3), do, round(do/nrow(df), 3))
  if(tolower(plot) == "ma"){
    
  	if(!is.null(xlim)){
  		limLog2mean2 <- xlim
  	}

   	if(!is.null(ylim)){
  		limLog2FC2 <- ylim
  	}

    out <- ArchR::ggPoint(df$log2Mean, df$log2FoldChange, df$color, pal = colors, xlabel = "Log2 Mean", ylabel = "Log2 Fold Change",
      size = size, alpha = alpha, title = title, xlim = limLog2mean2, ylim = limLog2FC2, rastr = rastr, ...) +
      geom_hline(yintercept = lfcint, lty = "dashed", size = 0.25) + scale_y_continuous(limits = limLog2FC2, breaks = seq(-50,50,2), expand = FALSE)
  
  }else if(tolower(plot) == "volcano"){
   
    out <- ggPoint(df$log2FoldChange, -log10(df$padj), df$color, pal = colors, xlabel = "Log2 Fold Change", ylabel = "-Log10 Pvalue",
      size = size, alpha = alpha, title = title, xlim = limLog2FC2, ylim = limPadj2, rastr = rastr, ...) +
      geom_vline(xintercept = lfcint, lty = "dashed") + scale_x_continuous(limits = limLog2FC2, breaks = seq(-50,50,2), expand = FALSE)
  
  }else if(tolower(plot) == "onetoone"){
    df$log2MeanTop <- df$log2Mean - df$log2FoldChange * 0.5
    df$log2MeanBottom <- df$log2Mean + df$log2FoldChange * 0.5
    limits <- quantile(c(df$log2MeanTop,df$log2MeanBottom), c(0.0001,max(limLog2mean)))
    out <- ggPoint(df$log2MeanTop, df$log2MeanBottom, df$color, pal = colors, xlabel = "Bottom Log2 Mean", ylabel = "Top Log2 Mean",
      size = size, alpha = alpha, title = title, xlim = limits, ylim = limits, rastr = rastr) +
      geom_abline(slope = 1, intercept = lfcint, lty = "dashed")
  }

  out <- out + theme(legend.position = "bottom")

  return(out)
}













#token <- "a0cc5ca693a5853677443fb9d9028121323c9f28"
#devtools::install_github("https://github.com/jgranja24/ArchRBrowser.git", auth_token = token)

library(ArchR)
library(ArchRx)
library(ggplot2)
library(ArchRBrowser)
set.seed(1)

rleToDF <- function(rle){
    runL <- c(0, cumsum(runLength(rle)))
    runLAll <- c(runL, runL - 1, runL + 1) #Hackish / Lazy attempt
    runLAll <- runLAll[runLAll > 0 & runLAll <= sum(runLength(rle))]
    runLAll <- runLAll[order(runLAll, decreasing = FALSE)]
    runValue <- as.numeric(rle[runLAll])
    df <- data.frame(x = runLAll, y = runValue)     
    return(df)
}

bwToDF <- function (bwList, gr, name){
    gr <- subsetSeqnames(gr[1], as.character(seqnames(gr[1])))
    rleList <- lapply(bwList, function(x) rtracklayer::import.bw(x, 
        selection = gr, as = ("Rle"))[[as.character(seqnames(gr))]])
    rleMean <- Reduce("+", rleList)/length(rleList)
    df <- rleToDF(rleMean)
    df$chr <- as.character(seqnames(gr))
    if (missing(name)) {
        df$name <- bwList[1]
    }
    else {
        df$name = name
    }
    return(df)
}

ArchRHeatmapTrack2 <- function(df, gr, resolution = 100){
    bwTile <- function(bwList, gr, tile, int = TRUE){
      rleList <- lapply(bwList, function(x) rtracklayer::import.bw(x, selection = gr, as = ("Rle"))[[as.character(seqnames(gr)[1])]])
      rleMean <- Reduce("+",rleList) / length(rleList)
      rleVal <- rleMean[tile]
      if(int){
        rleVal <- as.integer(rleVal)
      }
      return(rleVal)
    }
    df <- data.frame(df)
    #only take first GRange seqname
    gr <- subsetSeqnames(gr,as.character(seqnames(gr[1])))
    tile <- slidingWindows(gr, resolution, resolution)
    names(tile) <- paste0("p", seq_len(length(gr)))
    tile <- unlist(tile)
    s <- start(tile)
    id <- names(tile)
    mat <- lapply(unique(df$name), function(x){ 
        bwTile(paste0(df$bw[which(df$name == x)]), gr, s)
    }) %>% Reduce("rbind", .)
    rownames(mat) <- unique(df$name)
    matList <- lapply(unique(id),function(x){
        mat[,which(id==x)]
    })
    return(matList)
}

ArchRPlotBWTracks <- function (df, gr, color_map = NULL, ...){
    df <- data.frame(df)
    stopifnot(any(colnames(df) %in% "name"))
    stopifnot(any(colnames(df) %in% "bw"))
    gr <- subsetSeqnames(gr[1], as.character(seqnames(gr[1])))
    dfList <- lapply(unique(df$name), function(x) {
        dfx <- bwToDF(paste0(df$bw[which(df$name == x)]), gr, name = paste0(x))
        dfaddS <- data.frame(x=start(gr), y=0.00001,chr=seqnames(gr)[1], name = paste0(x))
        dfaddE <- data.frame(x=end(gr), y=0.00001,chr=seqnames(gr)[1], name = paste0(x))
        rbind(dfx, dfaddS, dfaddE)
    })
    dfAll <- Reduce("rbind", dfList)
    dfAll$name <- factor(dfAll$name, levels = unique(df$name))
    if (any(colnames(df) %in% "color")) {
        color_map <- paste0(df$color)
        names(color_map) <- paste0(df$name)
        p <- ArchRPlotTracks(dfAll, gr = gr, color_map = color_map, ...)
    }
    else {
        p <- ArchRPlotTracks(dfAll, gr = gr, ...)
    }
    return(p)
}

ArchRPlotTracks <- function (df, gr, color_map = NULL, type = "P", title = "", font_size = 7, 
    border_width = 0.4, tick_width = 0.4, axis_tick_length_cm = 0.25, 
    facet_font_size = 7, legend = TRUE, xlim = NULL, ylim = NULL, 
    hide_title = FALSE, hide_x = FALSE, hide_y = FALSE){

    gr <- subsetSeqnames(gr[1], as.character(seqnames(gr[1])))
    stopifnot(any(colnames(df) %in% "name"))
    stopifnot(any(colnames(df) %in% "x"))
    stopifnot(any(colnames(df) %in% "y"))
    stopifnot(any(colnames(df) %in% "chr"))
    stopifnot(any(df$chr %in% as.character(seqnames(gr))))
    if (is.null(color_map)) {
        color_map <- paletteDiscrete(set = "stallion", df$name)
    }
    if (!is.null(ylim)) {
        ylim <- quantile(df$y, ylim)
        df$y[df$y < ylim[1]] <- ylim[1]
        df$y[df$y > ylim[2]] <- ylim[2]
    }
    p <- ggplot(df, aes_string("x", "y", color = "name", fill = "name")) + 
        #geom_area(stat = "identity") +
        facet_grid(name ~ .) + ylab("Normalized Reads") + 
        scale_color_manual(values = color_map) + scale_fill_manual(values = color_map) + 
        scale_x_continuous(limits = c(start(gr), end(gr)), expand = c(0, 0)) + 
        theme_ArchR(base_size = font_size, base_rect_size = border_width, 
        base_line_size = tick_width, legend_position = "right", 
        axis_tick_length_cm = axis_tick_length_cm) + theme(panel.spacing = unit(0, 
        "lines"), axis.title.x = element_blank(), strip.text = element_text(size = facet_font_size, 
        color = "black", margin = margin(0, 0.35, 0, 0.35, "cm")), 
        strip.background = element_rect(color = "black")) + ggtitle(paste0(paste0(as.character(seqnames(gr)), 
        ":", start(gr) - 1, "-", end(gr), " min=0 max=",round(max(df$y),3)," ",title)))
    if (!legend) {
        p <- p + theme(legend.position = "none")
    }
    if (type == "P") {
        p <- p + geom_ribbon(aes(x = x, ymax = y, ymin = 0, fill = name))
    }
    if (!is.null(ylim)) {
        p <- p + scale_y_continuous(limits = ylim, expand = c(0, 0))
    }
    if (!is.null(xlim)) {
        p <- p + scale_x_continuous(limits = xlim, expand = c(0, 0))
    }
    if (hide_x) {
        p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
    }
    if (hide_y) {
        p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
    }
    if (hide_title) {
        p <- p + theme(title = element_blank())
    }
    return(p)
}

sampleNonOutliers <- function(values, n, outliers = 0.05){
    tabValues <- table(values)
    idx <- which(tabValues >= quantile(tabValues, outliers) & tabValues <= quantile(tabValues, 1 - outliers))
    valuesFilter <- names(tabValues)[idx]
    nSample(valuesFilter, n)
}

ArchRSCHeatmapTracks <- function(fragDF, gr, by = "RG", n = 100, pal = c("red","purple"), windowSize = 100, dedup = TRUE, font_size = 9, border_width = 0.4, tile_size = 0.5, ylab = "", xlab = "", facet_font_size = 6, hide_x = FALSE, hide_y = TRUE){
  scHeatDF <- function(fragDF, gr, by = "RG", fragMax = 1, windowSize = 100, dedup = TRUE, n = 250){
    hdf <- lapply(seq_along(fragDF$fragments), function(x){
      message(sprintf("%s of %s",x,length(fragDF$fragments)))
      fragx <- readRDS(paste0(fragDF$fragments[x]))
      message(paste0(length(fragx)," fragments..."))
      if(dedup){
        message("Dedup...")
        fragx <- unique(fragx)
      }
      if(!is.null(n)){
        id <- sampleNonOutliers(mcols(fragx)[,by], n) #nSample(unique(mcols(fragx)[,by]),n)
        fragx <- fragx[mcols(fragx)[,by] %in% id]
      }
      #make windows
      message("Windows...")
      windows <- slidingWindows(gr, windowSize, windowSize)[[1]]
      overlapDF <- suppressWarnings(DataFrame(findOverlaps(windows,fragx,ignore.strand=TRUE)))
      overlapDF$name <- mcols(fragx)[overlapDF[, 2], "RG"]
      overlapDF$id <- match(overlapDF$name, id)
      #overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
      m <- Matrix::sparseMatrix(i = overlapDF[, 1], j = overlapDF[, 4], x = rep(1, nrow(overlapDF)), dims = c(length(windows), n))
      m <- m[,sample(seq_col(m),ncol(m))]
      message("Heatmap...")
      heatDF <- data.frame(x=sparseRows(m),y=sparseCols(m),z=sparseVals(m))
      heatDF$z[heatDF$z > 1] <- 1
      heatDF$bp <- start(resize(windows,1,"center"))[heatDF$x]
      heatDF$name <- fragDF$name[x]
      heatDF$color <- fragDF$color[x]
      heatDF$nCells <- length(unique(mcols(fragx)[,by]))
      heatDF
    }) %>% Reduce("rbind",.)
    return(hdf)
  }
  hdf <- scHeatDF(fragDF=fragDF, gr=gr, by=by, windowSize=windowSize, dedup=dedup, n=n)
  hdf$name <- factor(hdf$name,levels=paste0(fragDF$name))
  colpal <- unique(hdf[,c("name","color")])[,2]
  names(colpal) <- unique(hdf[,c("name","color")])[,1]
  gg <- ggplot(hdf, aes(x=bp, y=y, width = windowSize, fill = name, color = name)) + 
    geom_tile(size = tile_size) + facet_grid(name ~ ., scales="free_y") + 
    theme_ArchR(base_line_size = border_width, base_rect_size = border_width, base_size = font_size,) + 
    scale_fill_manual(values=colpal) +
    scale_color_manual(values=colpal) +
    coord_cartesian(ylim=c(0,n),xlim=c(start(gr),end(gr)),expand=F) + 
    theme(legend.position = "right") +
    ylab(ylab) + xlab(xlab) +
    theme(panel.spacing = unit(0, "lines"), 
          axis.title.x = element_blank(), 
          strip.text = element_text(
            size = facet_font_size, 
            color = "black", 
            margin = margin(0, 0.35, 0, 0.35, "cm")
          ), strip.background = element_rect(color = "black")) 
  if (hide_x) {
    gg <- gg + theme(axis.title.x = element_blank(), axis.text.x  = element_blank(), axis.ticks.x = element_blank())
   }

  if (hide_y) {
    gg <- gg + theme(axis.title.y = element_blank(), axis.text.y  = element_blank(), axis.ticks.y= element_blank())
  }
  return(gg)
}

ArchRGeneTracks <- function (obj, gr, colors = c("red","blue"),font_size = 9, border_width = 0.4, segment_size = 2, box_size = 4, label_size = 2){
    gr <- subsetSeqnames(gr[1], as.character(seqnames(gr[1])))
    genes <- sort(sortSeqlevels(obj$genes), ignore.strand = T)
    exons <- sort(sortSeqlevels(obj$exons), ignore.strand = T)
    genesO <- data.frame(subsetByOverlaps(genes, gr, ignore.strand = T))
    if (nrow(genesO) > 0) {
        exonsO <- data.frame(subsetByOverlaps(exons, gr, ignore.strand = T))
        exonsO <- exonsO[which(exonsO$symbol %in% genesO$symbol), 
            ]
        genesO$facet = "TxDb Genes"
        genesO$start <- matrixStats::rowMaxs(cbind(genesO$start, 
            start(gr)))
        genesO$end <- matrixStats::rowMins(cbind(genesO$end, 
            end(gr)))
        genesO$cluster <- 0
        for (i in seq_len(nrow(genesO))) {
            if (i == 1) {
                genesO$cluster[i] <- 1
            }
            else {
                for (j in seq_len(max(genesO$cluster))) {
                  j_end <- rev(genesO$end)[match(rev(seq_len(max(genesO$cluster)))[j], 
                    rev(genesO$cluster))]
                  if (genesO$start[i] > j_end + median(genesO$width)) {
                    genesO$cluster[i] <- rev(genesO$cluster)[match(rev(seq_len(max(genesO$cluster)))[j], 
                      rev(genesO$cluster))]
                  }
                }
                if (genesO$cluster[i] == 0) {
                  genesO$cluster[i] <- genesO$cluster[i - 1] + 
                    1
                }
            }
        }
        exonsO$cluster <- genesO$cluster[match(exonsO$symbol, 
            genesO$symbol)]
        col_values <- c(`-` = colors[2], `+` = colors[1], `*` = colors[1])
        p <- ggplot(data = genesO, aes(color = strand, fill = strand)) + 
            facet_grid(facet ~ .) + ylim(c(0.5, max(genesO$cluster) + 
            0.5)) + geom_segment(data = genesO[which(as.character(genesO$strand) != 
            "-"), ], aes(x = start, xend = end, y = cluster, 
            yend = cluster, color = strand), size = segment_size) + 
            geom_segment(data = genesO[which(as.character(genesO$strand) == 
                "-"), ], aes(x = end, xend = start, y = cluster, 
                yend = cluster, color = strand), size = segment_size) + 
            geom_segment(data = exonsO, aes(x = start, xend = end, 
                y = cluster, yend = cluster, color = strand), 
                size = box_size) + scale_x_continuous(limits = c(start(gr), 
            end(gr)), expand = c(0, 0)) + theme_ArchR(base_size = font_size, 
            base_line_size = border_width, base_rect_size = border_width) + 
            scale_color_manual(values = col_values, guide = F) + 
            scale_fill_manual(values = col_values) + theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
            theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
                axis.ticks.y = element_blank()) + guides(fill = guide_legend(override.aes = list(colour = NA, 
            shape = "c", size = 4))) + theme(legend.text = element_text(size = font_size))
        if (length(which(genesO$strand != "-")) > 0) {
            p <- p + ggrepel::geom_label_repel(data = genesO[which(genesO$strand != 
                "-"), ], aes(x = start, y = cluster, label = symbol, 
                color = strand, fill = NA), segment.color = "grey", 
                nudge_x = -0.01 * (end(gr) - start(gr)), nudge_y = -0.25, 
                size = label_size, direction = "x")
        }
        if (length(which(genesO$strand == "-")) > 0) {
            p <- p + ggrepel::geom_label_repel(data = genesO[which(genesO$strand == 
                "-"), ], aes(x = end, y = cluster, label = symbol, 
                color = strand, fill = NA), segment.color = "grey", 
                nudge_x = +0.01 * (end(gr) - start(gr)), nudge_y = 0.25, 
                size = label_size, direction = "x")
        }
    }else {
        df <- data.frame(facet = "TxDb Genes", start = 0, end = 0, 
            strand = "*", symbol = "none")
        col_values <- c(`*` = colors[1])
        p <- ggplot(data = df, aes(start, end, fill = strand)) + 
            geom_point() + facet_grid(facet ~ .) + theme_ArchR(base_size = font_size, 
            base_line_size = border_width, base_rect_size = border_width) + 
            scale_color_manual(values = col_values) + scale_x_continuous(limits = c(start(gr), 
            end(gr)), expand = c(0, 0)) + theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
            theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
                axis.ticks.y = element_blank())
    }
    return(p)
}

ArchRPlotDiffTracks <- function (df, gr, color_map = NULL, ...){
    df <- data.frame(df)
    stopifnot(any(colnames(df) %in% "name"))
    stopifnot(any(colnames(df) %in% "bw"))
    gr <- subsetSeqnames(gr[1], as.character(seqnames(gr[1])))
    dfList <- lapply(unique(df$name), function(x) {
        dfx <- bwToDF(paste0(df$bw[which(df$name == x)]), gr, name = paste0(x))
        dfaddS <- data.frame(x=start(gr), y=0.00001,chr=seqnames(gr)[1], name = paste0(x))
        dfaddE <- data.frame(x=end(gr), y=0.00001,chr=seqnames(gr)[1], name = paste0(x))
        rbind(dfx, dfaddS, dfaddE)
    })
    dfAll <- Reduce("rbind", dfList)
    dfAll$name <- factor(dfAll$name, levels = unique(df$name))
    dfAll$diff <- df$diff[match(dfAll$name,df$name)]
    print(head(dfAll))
    if (any(colnames(df) %in% "color")) {
        color_map <- paste0(df$color)
        names(color_map) <- paste0(df$name)
        p <- ArchRPlotDTracks(dfAll, gr = gr, color_map = color_map, ...)
    }
    else {
        p <- ArchRPlotDTracks(dfAll, gr = gr, ...)
    }
    return(p)
}


ArchRPlotDTracks <- function (df, gr, color_map = NULL, type = "P", title = "", font_size = 7, 
    border_width = 0.4, tick_width = 0.4, axis_tick_length_cm = 0.25, 
    facet_font_size = 7, legend = TRUE, xlim = NULL, ylim = NULL, 
    hide_title = FALSE, hide_x = FALSE, hide_y = FALSE){

    gr <- subsetSeqnames(gr[1], as.character(seqnames(gr[1])))
    stopifnot(any(colnames(df) %in% "name"))
    stopifnot(any(colnames(df) %in% "x"))
    stopifnot(any(colnames(df) %in% "y"))
    stopifnot(any(colnames(df) %in% "chr"))
    stopifnot(any(df$chr %in% as.character(seqnames(gr))))
    if (is.null(color_map)) {
        color_map <- paletteDiscrete(set = "stallion", df$name)
    }
    if (!is.null(ylim)) {
        ylim <- quantile(df$y, ylim)
        df$y[df$y < ylim[1]] <- ylim[1]
        df$y[df$y > ylim[2]] <- ylim[2]
    }
    p <- ggplot(df, aes_string("x", "y", color = "name", fill = "name")) + 
        geom_line() + facet_grid(diff ~ .) + ylab("Normalized Reads") + 
        scale_color_manual(values = color_map) + scale_fill_manual(values = color_map) + 
        scale_x_continuous(limits = c(start(gr), end(gr)), expand = c(0, 0)) + 
        theme_ArchR(base_size = font_size, base_rect_size = border_width, 
        base_line_size = tick_width, legend_position = "right", 
        axis_tick_length_cm = axis_tick_length_cm) + theme(panel.spacing = unit(0, 
        "lines"), axis.title.x = element_blank(), strip.text = element_text(size = facet_font_size, 
        color = "black", margin = margin(0, 0.35, 0, 0.35, "cm")), 
        strip.background = element_rect(color = "black")) + ggtitle(paste0(paste0(as.character(seqnames(gr)), 
        ":", start(gr) - 1, "-", end(gr), " min=0 max=",round(max(df$y),3)," ",title)))
    if (!legend) {
        p <- p + theme(legend.position = "none")
    }
    if (type == "P") {
        p <- p + geom_ribbon(aes(x = x, ymax = y, ymin = 0, fill = name))
    }
    if (!is.null(ylim)) {
        p <- p + scale_y_continuous(limits = ylim, expand = c(0, 0))
    }
    if (!is.null(xlim)) {
        p <- p + scale_x_continuous(limits = xlim, expand = c(0, 0))
    }
    if (hide_x) {
        p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
    }
    if (hide_y) {
        p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
    }
    if (hide_title) {
        p <- p + theme(title = element_blank())
    }
    return(p)
}

theme_ArchR <- function(...){
    ArchRx::theme_ArchR(...)
}

########################
# Start Here
########################

plotDir <- "~/Desktop"

#Get Gene Info
genes <- getTxDbGenes("mm10")
geo <- getTxDbGeneExons(txdb = getTxDb("mm10"), orgdb=getOrgDb("mm10"))

#Get BigWigs"/Volumes/JG_SSD_1/Data/LKB1/TKO_Experiment/"
bwPath <- "/Volumes/JG_SSD_1/Data/LKB1/Sik_Knockout/output/ATAC_Analysis/BigWigs/Samples"
bws <- list.files(bwPath)

bwDF1 <- data.frame(
	name = gsub("_Rep1", "", gsub("_Rep2", "", gsub("-10M-50bp.insertion.bw", "", bws))),
	bw = file.path(bwPath, bws),
	color = ifelse(grepl("TAM", ignore.case=TRUE, bws) & grepl("2695|2804",ignore.case=TRUE, bws) , "dodgerblue3", "firebrick3")
)


#o2 <- gsub("EtOH", "Tam", o)


bwDF <- bwDF1#[grep("_1_", bwDF1$name), ]

bwDF$name <- factor(paste0(bwDF$name), levels = c(o, o2))

bwDF <- bwDF[order(match(paste0(bwDF$name), c(o, o2))),]
#Determine Plot Regions
geneToPlot <- "Sox17"

nUp <- 260000
nDo <- 75000

region <- genes[which(mcols(genes)$symbol == geneToPlot)] %>% 
    resize(1, "start") %>% 
    extendGR(upstream = nUp, downstream = nDo)
strand(region) <- "*"

ggGT <- ArchRGeneTracks(geo, region, colors = c("black","darkgrey"), label_size = 3) + 
    guides(fill = FALSE, colour = FALSE) +
    theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))

ggATAC <- ArchRPlotBWTracks(bwDF, region, hide_y = TRUE, hide_x = TRUE, ylim = c(0,0.998)) +
    theme(strip.text.y = element_text(angle = 0)) + guides(fill = FALSE, colour = FALSE) +
    theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))

pdf(file.path(plotDir, paste0(geneToPlot,"-",nUp,"-",nDo,"-SgSik.pdf")), height = 4, width = 5, onefile = FALSE)
ggAlignPlots(ggATAC, ggGT, sizes = c(3, 1, 1),type = "v")
dev.off()











