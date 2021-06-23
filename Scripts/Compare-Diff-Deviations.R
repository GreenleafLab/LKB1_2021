library(ArchR)
library(chromVAR)

computeDevDiff <- function(
  dev1 = NULL,
  up1 = NULL,
  do1 = NULL,
  name1 = NULL,
  dev2 = NULL,
  up2 = NULL,
  do2 = NULL,
  name2 = NULL,
  nPerm = 10000,
  mouseHuman = FALSE,
  assay = "z",
  threads = 10
  ){

  set.seed(1)

  #Suset deviations
  dev1s <- assays(dev1[, colnames(dev1) %in% c(up1,do1)])[[assay]]
  dev2s <- assays(dev2[, colnames(dev2) %in% c(up2,do2)])[[assay]]

  if(mouseHuman){

    TF1 <- tolower(stringr::str_split(rownames(dev1s), pattern = "\\_", simplify=TRUE)[,1])
    TF2 <- tolower(stringr::str_split(rownames(dev2s), pattern = "\\_", simplify=TRUE)[,1])
    common <- intersect(tolower(TF1), tolower(TF2))

    matchTF <- lapply(seq_along(common), function(x){

      TF1x <- rownames(dev1)[which(TF1 %in% common[x])]
      TF2x <- rownames(dev2)[which(TF2 %in% common[x])]

      df <- expand.grid(TF1 = TF1x, TF2 = TF2x)
      df$Common <- common[x]
      df$N <- nrow(df)
      df

    }) %>% Reduce("rbind", .)


  }else{

    common <- intersect(rownames(dev1), rownames(dev2))
    common2 <- tolower(stringr::str_split(common, pattern = "\\_", simplify=TRUE)[,1])
    tabCommon2 <- table(common2)

    matchTF <- data.frame(Name1 = name1, TF1 = common, Name2 = name2, TF2 = common, Common = common, N = tabCommon2[common2])

  }

  message("Testing ", nrow(matchTF), " Deviation Differences")

  diffDF <- ArchR:::.safelapply(seq_len(nrow(matchTF)), function(x){

    if(x %% 50 == 0){
      message(x, " of ", nrow(matchTF))
    }
  
    d1 <- dev1s[paste0(matchTF$TF1[x]), ]
    t1 <- rep(2, length(d1))
    t1[names(d1) %in% do1] <- 1

    d2 <- dev2s[paste0(matchTF$TF2[x]), ]
    t2 <- rep(2, length(d2))
    t2[names(d2) %in% do2] <- 1

    diff1 <- mean(d1[t1==2]) - mean(d1[t1==1])
    pval1 <- t.test(d1[t1==2], d1[t1==1])$p.value
    diff2 <- mean(d2[t2==2]) - mean(d2[t2==1])
    pval2 <- t.test(d2[t2==2], d2[t2==1])$p.value
    distDev <- (diff1)^2 + (diff2)^2

    nullDist <- lapply(seq_len(nPerm), function(y){

      s1 <- sample(t1, length(t1))
      s2 <- sample(t2, length(t2))

      (mean(d1[s1==2]) - mean(d1[s1==1]))^2 + (mean(d2[s2==2]) - mean(d2[s2==1]))^2

    }) %>% unlist

    p <- min((sum(nullDist > distDev) + 1)/ nPerm , 1)

    data.frame(
      pval = p, 
      dist = distDev, 
      bgdDistMean = mean(nullDist), 
      bgdDistSd = sd(nullDist),
      zDist = (distDev - mean(nullDist)) / sd(nullDist),
      zPval = 2 * pnorm(-abs((distDev - mean(nullDist)) / sd(nullDist))),
      diff1 = diff1,
      pval1 = pval1,
      diff2 = diff2,
      pval2 = pval2
    )

  }, threads = threads) %>% Reduce("rbind", .)


  allDF <- cbind(matchTF, diffDF)
  
  qvalObj <- qvalue::qvalue(pmin(allDF$pval, 1))
  allDF$nPerm <- nPerm
  allDF$Pvalue <- qvalObj$pvalues
  allDF$Qvalue <- qvalObj$qvalues
  allDF$LocalFDR <- qvalObj$lfdr

  qvalObj <- qvalue::qvalue(pmin(allDF$zPval, 1))
  allDF$PvalueZ <- qvalObj$pvalues
  allDF$QvalueZ <- qvalObj$qvalues
  allDF$LocalFDRZ <- qvalObj$lfdr

  allFilered <- allDF[order(allDF$dist, decreasing = TRUE), ]
  allFilered <- allFilered[!duplicated(allFilered$Common), ]

  out <- list(filtered = allFilered, all = allDF)

}

setwd("/Volumes/JG_SSD_1/Data/LKB1/DiffCV")

colorList <- list(
  "NKX" = c("NKX"),
  "NFI" = c("NFI"),
  "MEF2" = c("MEF2"),
  "TEAD" = c("TEAD"),
  "RUNX" = c("RUNX", "RUNT"),
  "SOX" = c("SOX", "SRY"),
  #"MESP" = c("MESP"),
  "ZEB" = c("ZEB", "SNAI", "MESP"),
  "FOXOther" = c("FOX"),
  "FOXA" = c("FOXA1", "FOXA2"), 
  "NR4A" = c("NR4A"),
  "BZIP" = c("JUN","FOS","BATF","FOSL", "ATF3", "JDP", "BACH"),
  "ONECUT" = c("ONECUT")
)

pal <- paletteDiscrete(values = names(colorList))
pal["SigOther"] <- "gray47"
pal["NA"] <- "lightgrey"
pal["NKX"] <- "#00ADEE"
pal["NFI"] <- "#C06CAB"
pal["MEF2"] <- "#FFDE17"
pal["TEAD"] <- "#F15A29"
pal["RUNX"] <- "#8B5E3C"
pal["SOX"] <- "#00A551"
pal["ZEB"] <- "#92278F"
pal["FOXOther"] <- "bisque"
pal["FOXA"] <- "#272E6A"
pal["NR4A"] <- "#EC008C"
pal["BZIP"] <- "#D51F26"
pal["ONECUT"] <- "slategray1"


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

#1. Compare Mouse Lines to TCGA
dev1 <- readRDS("../Initial_2695_2804_2381_2565/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")
up1 <- grep("TAM", grep("2695|2804", colnames(dev1),value=TRUE), value = TRUE)
do1 <- grep("EtOH", grep("2695|2804", colnames(dev1),value=TRUE), value = TRUE)

dev2 <- readRDS("../TCGA_LUAD_Analysis/Save-ChromVAR-Human-Deviatons.rds")
up2 <- colnames(dev2)[colData(dev2)$Group %in% paste0("LUAD_", c(8,21,11,17,3,22,5,18,12,20,14,9,4))]
do2 <- colnames(dev2)[colData(dev2)$Group %in% paste0("LUAD_", c(1,16,10,15,19,6,13,2))]

diffDev_MouseLine_HumanTCGA <- computeDevDiff(
  dev1 = dev1,
  up1 = up1,
  do1 = do1,
  name1 = "MouseLines",
  dev2 = dev2,
  up2 = up2,
  do2 = do2,
  name2 = "HumanTCGA",
  nPerm = 1000,
  mouseHuman = TRUE,
  assay = "z",
  threads = 14
)

saveRDS(diffDev_MouseLine_HumanTCGA, "Save_diffDev_MouseLine_HumanTCGA.rds")

#2. Compare Human Lines to TCGA
dev1 <- readRDS("../Human_OE_Exp2/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")
dev1 <- dev1[,grep("GFP", colnames(dev1))]
up1 <- grep("1650|1975|358|2009", colnames(dev1), value = TRUE)
do1 <- colnames(dev1)[!grepl("1650|1975|358|2009", colnames(dev1))]

dev2 <- readRDS("../TCGA_LUAD_Analysis/Save-ChromVAR-Human-Deviatons.rds")
up2 <- colnames(dev2)[colData(dev2)$Group %in% paste0("LUAD_", c(8,21,11,17,3,22,5,18,12,20,14,9,4))]
do2 <- colnames(dev2)[colData(dev2)$Group %in% paste0("LUAD_", c(1,16,10,15,19,6,13,2))]

diffDev_HumanLine_HumanTCGA <- computeDevDiff(
  dev1 = dev1,
  up1 = up1,
  do1 = do1,
  name1 = "HumanLines",
  dev2 = dev2,
  up2 = up2,
  do2 = do2,
  name2 = "HumanTCGA",
  nPerm = 1000,
  mouseHuman = FALSE,
  assay = "z",
  threads = 14
)

saveRDS(diffDev_HumanLine_HumanTCGA, "Save_diffDev_HumanLine_HumanTCGA.rds")

#3. Compare KP vs KPL Mets - KP vs KPL Tumors
dev1 <- readRDS("../KP_KPL_Tumors/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")
up1 <- grep("KPL_Met", colnames(dev1), value = TRUE)
do1 <- colnames(dev1)[!grepl("KPL_Met", colnames(dev1))]

dev2 <- readRDS("../KP_KPL_Tumors/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")
up2 <- grep("KP_Met", colnames(dev2), value = TRUE)
do2 <- colnames(dev2)[!grepl("KP_Met", colnames(dev2))]

diffDev_KP_KPL_TM <- computeDevDiff(
  dev1 = dev1,
  up1 = up1,
  do1 = do1,
  name1 = "KP_MT",
  dev2 = dev2,
  up2 = up2,
  do2 = do2,
  name2 = "KPL_MT",
  nPerm = 1000,
  mouseHuman = FALSE,
  assay = "z",
  threads = 14
)

saveRDS(diffDev_KP_KPL_TM, "Save_diffDev_KP_KPL_TM.rds")

#4. Compare KP - KPL Tumors vs TCGA 
dev1 <- readRDS("../KP_KPL_Tumors/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")
up1 <- grep("KP_Primary", colnames(dev1), value = TRUE)
do1 <- grep("KPL_Primary", colnames(dev1), value = TRUE)

dev2 <- readRDS("../TCGA_LUAD_Analysis/Save-ChromVAR-Human-Deviatons.rds")
up2 <- colnames(dev2)[colData(dev2)$Group %in% paste0("LUAD_", c(8,21,11,17,3,22,5,18,12,20,14,9,4))]
do2 <- colnames(dev2)[colData(dev2)$Group %in% paste0("LUAD_", c(1,16,10,15,19,6,13,2))]

diffDev_MT_HT <- computeDevDiff(
  dev1 = dev1,
  up1 = up1,
  do1 = do1,
  name1 = "MouseTumors",
  dev2 = dev2,
  up2 = up2,
  do2 = do2,
  name2 = "HumanTCGA",
  nPerm = 1000,
  mouseHuman = TRUE,
  assay = "z",
  threads = 14
)

saveRDS(diffDev_MT_HT, "Save_diffDev_MT_HT.rds")

#5. Compare KP - KPL Mets vs TCGA 
dev1 <- readRDS("../KP_KPL_Tumors/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")
up1 <- grep("KP_Met", colnames(dev1), value = TRUE)
do1 <- grep("KPL_Met", colnames(dev1), value = TRUE)

dev2 <- readRDS("../TCGA_LUAD_Analysis/Save-ChromVAR-Human-Deviatons.rds")
up2 <- colnames(dev2)[colData(dev2)$Group %in% paste0("LUAD_", c(8,21,11,17,3,22,5,18,12,20,14,9,4))]
do2 <- colnames(dev2)[colData(dev2)$Group %in% paste0("LUAD_", c(1,16,10,15,19,6,13,2))]

diffDev_MM_HT <- computeDevDiff(
  dev1 = dev1,
  up1 = up1,
  do1 = do1,
  name1 = "MouseMets",
  dev2 = dev2,
  up2 = up2,
  do2 = do2,
  name2 = "HumanTCGA",
  nPerm = 1000,
  mouseHuman = TRUE,
  assay = "z",
  threads = 14
)

saveRDS(diffDev_MM_HT, "Save_diffDev_MM_HT.rds")

#6. Compare KP - KPL Tumors vs In Vitro 
dev1 <- readRDS("../KP_KPL_Tumors/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")
up1 <- grep("KP_Primary", colnames(dev1), value = TRUE)
do1 <- grep("KPL_Primary", colnames(dev1), value = TRUE)

dev2 <- readRDS("../Initial_2695_2804_2381_2565/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")
up2 <- grep("TAM", grep("2695|2804", colnames(dev2),value=TRUE), value = TRUE)
do2 <- grep("EtOH", grep("2695|2804", colnames(dev2),value=TRUE), value = TRUE)

diffDev_MT_ML <- computeDevDiff(
  dev1 = dev1,
  up1 = up1,
  do1 = do1,
  name1 = "MouseTumors",
  dev2 = dev2,
  up2 = up2,
  do2 = do2,
  name2 = "MouseLines",
  nPerm = 1000,
  mouseHuman = FALSE,
  assay = "z",
  threads = 14
)

saveRDS(diffDev_MT_ML, "Save_diffDev_MT_ML.rds")

#7. Compare KP - KPL Mets vs In Vitro 
dev1 <- readRDS("../KP_KPL_Tumors/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")
up1 <- grep("KP_Met", colnames(dev1), value = TRUE)
do1 <- grep("KPL_Met", colnames(dev1), value = TRUE)

dev2 <- readRDS("../Initial_2695_2804_2381_2565/output/ATAC_Analysis/Main/CisBP-chromVAR.rds")
up2 <- grep("TAM", grep("2695|2804", colnames(dev2),value=TRUE), value = TRUE)
do2 <- grep("EtOH", grep("2695|2804", colnames(dev2),value=TRUE), value = TRUE)

diffDev_MM_ML <- computeDevDiff(
  dev1 = dev1,
  up1 = up1,
  do1 = do1,
  name1 = "MouseMets",
  dev2 = dev2,
  up2 = up2,
  do2 = do2,
  name2 = "MouseLines",
  nPerm = 1000,
  mouseHuman = FALSE,
  assay = "z",
  threads = 14
)

saveRDS(diffDev_MM_ML, "Save_diffDev_MM_ML.rds")

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

#1. Compare Mouse Lines to TCGA
plotDF <- readRDS("Save_diffDev_MouseLine_HumanTCGA.rds")[[1]]

plotDF$Family <- "NA"
plotDF$Sig <- FALSE
minCut <- 1
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Sig"] <- TRUE
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Family"] <- "SigOther"

for(i in seq_along(colorList)){
  idx <- unlist(lapply(colorList[[i]], function(x)grep(x,plotDF$Common, value = TRUE, ignore.case = TRUE)))
  plotDF$Family[which(plotDF$Common %in% idx)] <- names(colorList)[i]
}

plotDF <- plotDF[which(plotDF$Family != "FOXOther"), ]
plotDF[which(!plotDF$Sig), "Family"] <- "NA"

idxCenter <- which(abs(plotDF$diff1) < 5 & abs(plotDF$diff2) < 5)
plotCenter <- plotDF[idxCenter, ]
plotCenter <- plotCenter[order(plotCenter$dist, decreasing = TRUE), ]
plotCenter <- plotCenter[sample(seq_len(nrow(plotCenter)), nrow(plotCenter)), ]
plotDF2 <- plotDF[-idxCenter,]
plotDF2 <- plotDF2[order(plotDF2$dist, decreasing = FALSE), ]
plotDF <- rbind(plotDF2 , plotCenter)

max <- 1000
plotDF$diff1[plotDF$diff1 > max] <- max
plotDF$diff1[plotDF$diff1 < -max] <- -max
plotDF$diff2[plotDF$diff2 > max] <- max
plotDF$diff2[plotDF$diff2 < -max] <- -max

pdf("Plot-diffDev_MouseLine_HumanTCGA.pdf", width = 5, height = 5, useDingbats = FALSE)

xlim <- extendrange(plotDF$diff2, f = 0.05)
ylim <- extendrange(plotDF$diff1, f = 0.01)
ylim[2] <- ylim[2] * 1.1

p <- ggplot(data = plotDF, aes(diff2, diff1, color = Family, fill = Family)) +
  geom_point(size = 2, shape = 21, stroke = 0.2, color = "black") + 
  scale_fill_manual(values = pal) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = 0, lty = "dashed") +
  xlab("ChromVAR Human TCGA Type1 - Type2") +
  ylab("ChromVAR Mouse Lines Type1 - Type2") + 
  coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim) + 
  theme_ArchR() + 
  ggtitle(max)

ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5)

dev.off()


#2. Compare Human Lines to TCGA
plotDF <- readRDS("Save_diffDev_HumanLine_HumanTCGA.rds")[[1]]
plotDF$Family <- "NA"
plotDF$Sig <- FALSE
minCut <- 1
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Sig"] <- TRUE
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Family"] <- "SigOther"

for(i in seq_along(colorList)){
  idx <- unlist(lapply(colorList[[i]], function(x)grep(x,plotDF$Common, value = TRUE, ignore.case = TRUE)))
  plotDF$Family[which(plotDF$Common %in% idx)] <- names(colorList)[i]
}

plotDF <- plotDF[which(plotDF$Family != "FOXOther"), ]
plotDF[which(!plotDF$Sig), "Family"] <- "NA"

idxCenter <- which(abs(plotDF$diff1) < 5 & abs(plotDF$diff2) < 5)
plotCenter <- plotDF[idxCenter, ]
plotCenter <- plotCenter[order(plotCenter$dist, decreasing = TRUE), ]
plotCenter <- plotCenter[sample(seq_len(nrow(plotCenter)), nrow(plotCenter)), ]
plotDF2 <- plotDF[-idxCenter,]
plotDF2 <- plotDF2[order(plotDF2$dist, decreasing = FALSE), ]
plotDF <- rbind(plotDF2 , plotCenter)
plotDF <- rbind(plotDF[plotDF$Family %in% c("NA", "SigOther"),], plotDF[plotDF$Family %ni% c("NA", "SigOther"),])

max <- 30
plotDF$diff1[plotDF$diff1 > max] <- max
plotDF$diff1[plotDF$diff1 < -max] <- -max
plotDF$diff2[plotDF$diff2 > max] <- max
plotDF$diff2[plotDF$diff2 < -max] <- -max

pdf("Plot-diffDev_HumanLine_HumanTCGA.pdf", width = 5, height = 5, useDingbats = FALSE)

xlim <- extendrange(plotDF$diff2, f = 0)
ylim <- extendrange(plotDF$diff1, f = 0)

p <- ggplot(data = plotDF, aes(diff2, diff1, color = Family, fill = Family)) +
  geom_point(size = 2, shape = 21, stroke = 0.2, color = "black") + 
  scale_fill_manual(values = pal) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = 0, lty = "dashed") +
  xlab("ChromVAR Human TCGA Type1 - Type2") +
  ylab("ChromVAR Human Cell Line Type1 - Type2") + 
  coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim) + 
  theme_ArchR() + 
  ggtitle(max)

ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5)

dev.off()

#3. Compare KP Mets - Tumors vs KP LMets - Tumors
plotDF <- readRDS("Save_diffDev_KP_KPL_TM.rds")[[1]]

plotDF$Family <- "NA"
plotDF$Sig <- FALSE
plotDF[plotDF$Qvalue < 0.05 & sqrt(plotDF$dist) > 5, "Sig"] <- TRUE
plotDF[plotDF$Qvalue < 0.05 & sqrt(plotDF$dist) > 5, "Family"] <- "SigOther"

for(i in seq_along(colorList)){
  idx <- unlist(lapply(colorList[[i]], function(x)grep(x,plotDF$Common, value = TRUE, ignore.case = TRUE)))
  plotDF$Family[which(plotDF$Common %in% idx)] <- names(colorList)[i]
}

plotDF <- plotDF[which(plotDF$Family != "FOXOther"), ]
plotDF[which(!plotDF$Sig), "Family"] <- "NA"

idxCenter <- which(abs(plotDF$diff1) < 3 & abs(plotDF$diff2) < 3)
plotCenter <- plotDF[idxCenter, ]
plotCenter <- plotCenter[order(plotCenter$dist, decreasing = TRUE), ]
plotCenter <- plotCenter[sample(seq_len(nrow(plotCenter)), nrow(plotCenter)), ]
plotDF2 <- plotDF[-idxCenter,]
plotDF2 <- plotDF2[order(plotDF2$dist, decreasing = FALSE), ]
plotDF <- rbind(plotDF2 , plotCenter)
plotDF <- rbind(plotDF[plotDF$Family %in% c("NA", "SigOther"),], plotDF[plotDF$Family %ni% c("NA", "SigOther"),])


max <- 40
plotDF$diff1[plotDF$diff1 > max] <- max
plotDF$diff1[plotDF$diff1 < -max] <- -max
plotDF$diff2[plotDF$diff2 > max] <- max
plotDF$diff2[plotDF$diff2 < -max] <- -max

pdf("Plot-diffDev_KP_KPL_TM.pdf", width = 5, height = 5, useDingbats = FALSE)

xlim <- extendrange(c(-max,max), f = 0)
ylim <- extendrange(c(-max,max), f = 0)

p <- ggplot(data = plotDF, aes(diff2, diff1, color = Family, fill = Family)) +
  geom_point(size = 3, shape = 21, stroke = 0.2, color = "black") + 
  scale_fill_manual(values = pal) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = 0, lty = "dashed") +
  xlab("ChromVAR KP Met - Tumor") +
  ylab("ChromVAR KPL Met - Tumor") + 
  coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim) + 
  theme_ArchR() + 
  geom_abline(slope = 1, intercept = 0, lty = "dashed") +
  ggtitle(max)

ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5)

dev.off()

#4. Compare KP - KPL Tumors vs TCGA 
plotDF <- readRDS("Save_diffDev_MT_HT.rds")[[1]]

plotDF$Family <- "NA"
plotDF$Sig <- FALSE
minCut <- 1
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Sig"] <- TRUE
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Family"] <- "SigOther"

for(i in seq_along(colorList)){
  idx <- unlist(lapply(colorList[[i]], function(x)grep(x,plotDF$Common, value = TRUE, ignore.case = TRUE)))
  plotDF$Family[which(plotDF$Common %in% idx)] <- names(colorList)[i]
}

plotDF <- plotDF[which(plotDF$Family != "FOXOther"), ]
plotDF[which(!plotDF$Sig), "Family"] <- "NA"

idxCenter <- which(abs(plotDF$diff1) < 3 & abs(plotDF$diff2) < 3)
plotCenter <- plotDF[idxCenter, ]
plotCenter <- plotCenter[order(plotCenter$dist, decreasing = TRUE), ]
plotCenter <- plotCenter[sample(seq_len(nrow(plotCenter)), nrow(plotCenter)), ]
plotDF2 <- plotDF[-idxCenter,]
plotDF2 <- plotDF2[order(plotDF2$dist, decreasing = FALSE), ]
plotDF <- rbind(plotDF2 , plotCenter)
plotDF <- rbind(plotDF[plotDF$Family %in% c("NA", "SigOther"),], plotDF[plotDF$Family %ni% c("NA", "SigOther"),])


max <- 20
plotDF$diff1[plotDF$diff1 > max] <- max
plotDF$diff1[plotDF$diff1 < -max] <- -max
plotDF$diff2[plotDF$diff2 > max] <- max
plotDF$diff2[plotDF$diff2 < -max] <- -max

pdf("Plot-diffDev_MT_HT.pdf", width = 5, height = 5, useDingbats = FALSE)

xlim <- extendrange(plotDF$diff2, f = 0)
ylim <- extendrange(plotDF$diff1, f = 0)

p <- ggplot(data = plotDF, aes(diff2, diff1, color = Family, fill = Family)) +
  geom_point(size = 2, shape = 21, stroke = 0.2, color = "black") + 
  scale_fill_manual(values = pal) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = 0, lty = "dashed") +
  xlab("ChromVAR Human TCGA Type1 - Type2") +
  ylab("ChromVAR KP - KPL Tumor") + 
  coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim) + 
  theme_ArchR() +
  ggtitle(max)

ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5)

dev.off()

#5. Compare KP - KPL Mets vs TCGA 
plotDF <- readRDS("Save_diffDev_MM_HT.rds")[[1]]

plotDF$Family <- "NA"
plotDF$Sig <- FALSE
minCut <- 1
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Sig"] <- TRUE
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Family"] <- "SigOther"

for(i in seq_along(colorList)){
  idx <- unlist(lapply(colorList[[i]], function(x)grep(x,plotDF$Common, value = TRUE, ignore.case = TRUE)))
  plotDF$Family[which(plotDF$Common %in% idx)] <- names(colorList)[i]
}

plotDF <- plotDF[which(plotDF$Family != "FOXOther"), ]
plotDF[which(!plotDF$Sig), "Family"] <- "NA"

idxCenter <- which(abs(plotDF$diff1) < 3 & abs(plotDF$diff2) < 3)
plotCenter <- plotDF[idxCenter, ]
plotCenter <- plotCenter[order(plotCenter$dist, decreasing = TRUE), ]
plotCenter <- plotCenter[sample(seq_len(nrow(plotCenter)), nrow(plotCenter)), ]
plotDF2 <- plotDF[-idxCenter,]
plotDF2 <- plotDF2[order(plotDF2$dist, decreasing = FALSE), ]
plotDF <- rbind(plotDF2 , plotCenter)
plotDF <- rbind(plotDF[plotDF$Family %in% c("NA", "SigOther"),], plotDF[plotDF$Family %ni% c("NA", "SigOther"),])


max <- 50
plotDF$diff1[plotDF$diff1 > max] <- max
plotDF$diff1[plotDF$diff1 < -max] <- -max
plotDF$diff2[plotDF$diff2 > max] <- max
plotDF$diff2[plotDF$diff2 < -max] <- -max

pdf("Plot-diffDev_MM_HT.pdf", width = 5, height = 5, useDingbats = FALSE)

xlim <- extendrange(plotDF$diff2, f = 0)
ylim <- extendrange(plotDF$diff1, f = 0)

p <- ggplot(data = plotDF, aes(diff2, diff1, color = Family, fill = Family)) +
  geom_point(size = 2, shape = 21, stroke = 0.2, color = "black") + 
  scale_fill_manual(values = pal) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = 0, lty = "dashed") +
  xlab("ChromVAR Human TCGA Type1 - Type2") +
  ylab("ChromVAR KP - KPL Mets") + 
  coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim) + 
  theme_ArchR() +
  ggtitle(max)

ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5)

dev.off()

#6. Compare KP - KPL Mets vs In Vitro 
plotDF <- readRDS("Save_diffDev_MM_ML.rds")[[1]]

plotDF$Family <- "NA"
plotDF$Sig <- FALSE
minCut <- 1
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Sig"] <- TRUE
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Family"] <- "SigOther"

for(i in seq_along(colorList)){
  idx <- unlist(lapply(colorList[[i]], function(x)grep(x,plotDF$Common, value = TRUE, ignore.case = TRUE)))
  plotDF$Family[which(plotDF$Common %in% idx)] <- names(colorList)[i]
}

plotDF <- plotDF[which(plotDF$Family != "FOXOther"), ]
plotDF[which(!plotDF$Sig), "Family"] <- "NA"

idxCenter <- which(abs(plotDF$diff1) < 3 & abs(plotDF$diff2) < 3)
plotCenter <- plotDF[idxCenter, ]
plotCenter <- plotCenter[order(plotCenter$dist, decreasing = TRUE), ]
plotCenter <- plotCenter[sample(seq_len(nrow(plotCenter)), nrow(plotCenter)), ]
plotDF2 <- plotDF[-idxCenter,]
plotDF2 <- plotDF2[order(plotDF2$dist, decreasing = FALSE), ]
plotDF <- rbind(plotDF2 , plotCenter)
plotDF <- rbind(plotDF[plotDF$Family %in% c("NA", "SigOther"),], plotDF[plotDF$Family %ni% c("NA", "SigOther"),])


max <- 50
plotDF$diff1[plotDF$diff1 > max] <- max
plotDF$diff1[plotDF$diff1 < -max] <- -max
plotDF$diff2[plotDF$diff2 > max] <- max
plotDF$diff2[plotDF$diff2 < -max] <- -max

pdf("diffDev_MM_ML.pdf", width = 5, height = 5, useDingbats = FALSE)

xlim <- extendrange(plotDF$diff2, f = 0)
ylim <- extendrange(plotDF$diff1, f = 0)

p <- ggplot(data = plotDF, aes(diff2, diff1, color = Family, fill = Family)) +
  geom_point(size = 2, shape = 21, stroke = 0.2, color = "black") + 
  scale_fill_manual(values = pal) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = 0, lty = "dashed") +
  xlab("ChromVAR MouseLine LKB1 Proficient - LKB1 Deficient") +
  ylab("ChromVAR KP - KPL Mets") + 
  coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim) + 
  theme_ArchR() +
  ggtitle(max)

ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5)

dev.off()

#7. Compare KP - KPL Tumors vs In Vitro 
plotDF <- readRDS("Save_diffDev_MT_ML.rds")[[1]]

plotDF$Family <- "NA"
plotDF$Sig <- FALSE
minCut <- 1
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Sig"] <- TRUE
plotDF[plotDF$Qvalue < 0.05 & abs(plotDF$diff1) > minCut & abs(plotDF$diff2) > minCut & sqrt(plotDF$dist) > 5, "Family"] <- "SigOther"

for(i in seq_along(colorList)){
  idx <- unlist(lapply(colorList[[i]], function(x)grep(x,plotDF$Common, value = TRUE, ignore.case = TRUE)))
  plotDF$Family[which(plotDF$Common %in% idx)] <- names(colorList)[i]
}

plotDF <- plotDF[which(plotDF$Family != "FOXOther"), ]
plotDF[which(!plotDF$Sig), "Family"] <- "NA"

idxCenter <- which(abs(plotDF$diff1) < 3 & abs(plotDF$diff2) < 3)
plotCenter <- plotDF[idxCenter, ]
plotCenter <- plotCenter[order(plotCenter$dist, decreasing = TRUE), ]
plotCenter <- plotCenter[sample(seq_len(nrow(plotCenter)), nrow(plotCenter)), ]
plotDF2 <- plotDF[-idxCenter,]
plotDF2 <- plotDF2[order(plotDF2$dist, decreasing = FALSE), ]
plotDF <- rbind(plotDF2 , plotCenter)
plotDF <- rbind(plotDF[plotDF$Family %in% c("NA", "SigOther"),], plotDF[plotDF$Family %ni% c("NA", "SigOther"),])


max <- 20
plotDF$diff1[plotDF$diff1 > max] <- max
plotDF$diff1[plotDF$diff1 < -max] <- -max
plotDF$diff2[plotDF$diff2 > max] <- max
plotDF$diff2[plotDF$diff2 < -max] <- -max

pdf("diffDev_MT_ML.pdf", width = 5, height = 5, useDingbats = FALSE)

xlim <- extendrange(plotDF$diff2, f = 0)
ylim <- extendrange(plotDF$diff1, f = 0)

p <- ggplot(data = plotDF, aes(diff2, diff1, color = Family, fill = Family)) +
  geom_point(size = 2, shape = 21, stroke = 0.2, color = "black") + 
  scale_fill_manual(values = pal) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = 0, lty = "dashed") +
  xlab("ChromVAR MouseLine LKB1 Proficient - LKB1 Deficient") +
  ylab("ChromVAR KP - KPL Tumors") + 
  coord_equal(ratio = diff(xlim)/diff(ylim), xlim = xlim, ylim = ylim) + 
  theme_ArchR() +
  ggtitle(max)

ArchR:::.fixPlotSize(p, plotWidth = 5, plotHeight = 5)

dev.off()








































