library(ArchR)
library(chromVAR)
library(qvalue)

computeDevDiff <- function(
  dev1 = NULL, #ChromVAR Deviatons 1
  up1 = NULL, #Samples in 'up' category 1
  do1 = NULL, #Samples in 'down' category 1
  name1 = NULL, #Name 1
  dev2 = NULL, #ChromVAR Deviatons 2
  up2 = NULL, #Samples in 'up' category 2
  do2 = NULL, #Samples in 'down' category 2
  name2 = NULL, #Name 2
  nPerm = 10000, #Number of Permutations
  mouseHuman = FALSE, #Mouse vs Human
  assay = "z", #Deviation Assay Z scores
  threads = 10 #Threads
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

#Example Analysis

#chromVAR Results1
dev1 <- readRDS("../Initial_2695_2804_2381_2565/output/ATAC_Analysis/Main/CisBP-chromVAR.rds") #chromVAR SE
up1 <- grep("TAM", grep("2695|2804", colnames(dev1),value=TRUE), value = TRUE)
do1 <- grep("EtOH", grep("2695|2804", colnames(dev1),value=TRUE), value = TRUE)

#chromVAR Results1
dev2 <- readRDS("../TCGA_LUAD_Analysis/Save-ChromVAR-Human-Deviatons.rds") #chromVAR SE
up2 <- colnames(dev2)[colData(dev2)$Group %in% paste0("LUAD_", c(8,21,11,17,3,22,5,18,12,20,14,9,4))]
do2 <- colnames(dev2)[colData(dev2)$Group %in% paste0("LUAD_", c(1,16,10,15,19,6,13,2))]

#Compute Differential Devations Differences
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

#Save
saveRDS(diffDev_MouseLine_HumanTCGA, "Save_diffDev_MouseLine_HumanTCGA.rds")
















