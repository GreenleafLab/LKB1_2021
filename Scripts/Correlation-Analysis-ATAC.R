library(ArchR)

set.seed(1)

#Input Summarized Experiment
se <- readRDS("output/ATAC_Analysis/Main/ATAC-Summarized-Experiment.rds")

#Correlation Heatmap
library(pheatmap)
cm <- cor(assay(se))
diag(cm) <- NA
pdf("~/Desktop/Correlation-Heatmap.pdf", width = 6, height = 6)
pheatmap(cm, color = paletteContinuous(), border_color = "black")
dev.off()
