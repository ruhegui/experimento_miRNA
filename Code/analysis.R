#Verificar BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#Paquetes
list.of.packages = c("tximeta", "tximport", "limma", "edgeR", "tidyverse", "org.Mm.eg.db", "statmod", "pheatmap", "ggplotify", "ggrepel", "SummarizedExperiment", "patchwork", "xlsx", "ragg", "OmnipathR")
#Instalación por CRAN o Bioconductor
new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
if(length(new.packages)> 0) {
  for(i in new.packages) {
    if(i %in% available.packages()[,1]){ #Chequea si el paquete está en un repositorio CRAN y lo instala
      install.packages(i,dependencies=TRUE)
    }else {BiocManager::install(i, update = TRUE, ask = FALSE, version = BiocManager::version()) #Instala por BiocManager si no está en un repositorio CRAN
    }}
}
invisible(lapply(list.of.packages, FUN=library, character.only=TRUE))
rm(list.of.packages, new.packages)
setwd("/home/guille/RNA/RNA_2024/exp_mirna/exp_miRNA") #main wd
files = data.frame(read.table(file="Metadata/metadata.txt", header = TRUE, stringsAsFactors = F))
suppressPackageStartupMessages(library(SummarizedExperiment))
se = tximeta(files)                   
se <- addExons(se)
gse <- summarizeToGene(se, assignRanges="abundant") #Asigna rangos según la isoforma más abundante de los transcritos en vez de desde la isoforma más en 5' a la isoforma más en 3'. Los creadores lo recomiendan.
gse <- addIds(gse, "GENENAME", gene=TRUE)
gse <- addIds(gse, "SYMBOL", gene=TRUE)
gse <- addIds(gse, "ENTREZID", gene = TRUE)
y <- makeDGEList(gse)
rm(se)
sampleinfo <- read.delim(file ="Metadata/sampleinfo.txt", sep = " ")
sampleinfo[,2]
group = sampleinfo$group
group <- factor(group)
y$samples$group = group
y$samples

design = model.matrix(~0 + group)
colnames(design) = levels(group)
print(design)
keep <- filterByExpr(y, design)
print(summary(keep))
y <- y[keep, ]
points <- c(1:5) # Formas
colors <- c(1:5) #Colores
mds = plotMDS(y, col=colors[group], pch=points[group])

mds_data <- data.frame(
  Dim1 = mds$x,
  Dim2 = mds$y,
  Sample = colnames(mds$distance.matrix.squared),
  Group = as.factor(y$samples$group)
)
axis_label = round(mds$var.explained[1:2]*100)
library(ggConvexHull)
mds_data %>% ggplot(aes(x = Dim1, y = Dim2, col = Group, shape = Group)) +
  geom_point(size = 2) +
  labs(title = "MDS Analysis", x = paste0("Dimension 1: ", axis_label[1], "%"), y = paste0("Dimension 2: ", axis_label[2], "%"), ) +
  theme_minimal() +
  scale_shape_manual(values = rep(15:17, len = 7)) + 
  geom_convexhull(aes(fill=Group, 
                      colour= Group),
                  alpha = 0.2) + 
  theme(plot.title = element_text(hjust = 0.5))

y <- estimateDisp(y, design, method = "ROBUST")
plotBCV(y)
fit <- glmQLFit(y, design, method = "ROBUST")
plotQLDisp(fit)
logcpm = cpm(y, log=TRUE)
rownames(logcpm) = y$genes$SYMBOL #Nombre de genes
colnames(logcpm) =  paste(y$samples$group, 1:3, sep = "-")
group
contrast = makeContrasts(
  ATvsPreT = "((AterT0 + AterT15) - (preT0 + preT15))", #1
  T15vsT0 = "(AterT15 + preT15) - (AterT0 + preT0)", #2
  AT15vs0 = "AterT15 - AterT0", #3
  PreT15vs0 = "preT15 - preT0", #4
  AvsPreT15 = "AterT15 - preT15", #5
  AvsPreT0 = "AterT0 - preT0", #6
  "preT15 - control", #7
  levels=design)
colnames(contrast)[1]
res <- glmQLFTest(fit, contrast = contrast[,6])
res_corrected = topTags(res, n = Inf)
res$comparison
head(res_corrected, 20)
is.de = decideTests(res, adjust.method = "BH", p.value = 0.05, lfc = 1)
plotMD(res, status = is.de)

data = res_corrected$table
data$DE <- "NO"
data$DE[data$logFC > 1 & data$PValue < 0.05] <- "UP"
data$DE[data$logFC < -1 & data$PValue < 0.05] <- "DOWN"
ggplot(data = data, aes(x=logFC, y =-log10(PValue), col=DE, label = ifelse(abs(logFC) >= 2 & PValue <0.05, as.character(SYMBOL),  ''))) + geom_point() +
  geom_text_repel(hjust = 0, nudge_x = 0.1, color = "black") +
  theme_minimal() +
  labs(title = as.character(res$comparison)) +
  geom_hline(yintercept=-log10(0.05)) +
  geom_vline(xintercept=1) +
  geom_vline(xintercept=-1)
colnames(res_corrected)
res_corrected$table[1,c(2,5)]
setwd("Results/")
experimento = "exp_miRNA"
dir.create(experimento, showWarnings = FALSE)
if (!file.exists(paste0(experimento, "/",make.names(ProfessR::fix.names(res$comparison)), ".xlsx"))) {
  write.xlsx2(x = res_corrected[!is.na(res_corrected$table$gene_name) & res_corrected$table$PValue <= 0.05 ,c(2,5,16,18,19)], file = paste0(experimento, ".xlsx"), sheetName = colnames(contrast)[6], col.names = T, row.names = F, append = TRUE, )
  } else {
  print("Ya existe")
}
setwd("..")
?write.xlsx2
