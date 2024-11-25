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
files = data.frame(read.table(file="Metadata/metadata.txt", header = TRUE, stringsAsFactors = F))
suppressPackageStartupMessages(library(SummarizedExperiment))
se = tximeta(files)                   
se <- addExons(se)
gse <- summarizeToGene(se, assignRanges="abundant") #Asigna rangos según la isoforma más abundante de los transcritos en vez de desde la isoforma más en 5' a la isoforma más en 3'. Los creadores lo recomiendan.
y <- makeDGEList(gse)
rm(se)
sampleinfo <- read.delim(file ="Metadata/sampleinfo.txt", sep = ",")
group = paste0(sampleinfo[,1],sampleinfo[,2])
group[1:3] = rep("control",3)
group = factor(group) %>% relevel(., ref = "control")
term = factor(sampleinfo[,1]) %>% relevel(., ref = "control")
time = factor(sampleinfo[,2]) %>% relevel(., ref = "control")
y$samples$group = group
y$samples$term = term
y$samples$time = time
y$samples

design = model.matrix(~0 + group)
colnames(design) = levels(group)
print(design)
keep <- filterByExpr(y, design)
print(summary(keep))
y <- y[keep, ]
points <- c(1:3) # Formas
colors <- c(1:3) #Colores
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
rownames(logcpm) = y$genes$gene_name #Nombre de genes
colnames(logcpm) =  paste(y$samples$group, 1:3, sep = "-")
pairwisecomb = combn(levels(group), 2, function(x) paste(x[2], "-", x[1], sep = "")); pairwisecomb
contrast = makeContrasts(contrasts = pairwisecomb, levels=design)
RESs = vector(length = length(pairwisecomb))
names(RESs) = colnames(contrast)
for (i in c(1:length(RESs))){
  RESs[i] = topTags(glmQLFTest(fit, contrast = contrast[,i]), n = Inf)
}
DEGs = lapply(RESs, function(i) i %>% dplyr::filter(FDR <=0.1 ) %>% dplyr::select(description, gene_name, entrezid, logFC, PValue, FDR))
DEGs = DEGs[lapply(DEGs,nrow)>0]
data = lapply(RESs, function(i) i %>%
                dplyr::mutate(DE = case_when(logFC > 0 & FDR <0.1 ~ "UP",
                                             logFC < 0 & FDR <0.1 ~ "DOWN",
                                             .default = "NO")))
volcanoplot = function(data, name){
  ggplot(data = data, aes(x=logFC, y =-log10(PValue), col=DE)) + geom_point() +
    scale_color_manual(values = c("DOWN" = "firebrick", "UP" = "dodgerblue", "NO" = "grey")) +
    labs(title = name) +
    theme_minimal()
}
plots <- lapply(names(data), function(name) {
  volcanoplot(data[[name]], name)
})
combined_plot <- wrap_plots(plots)
print(combined_plot)

DEGs

DEG_selection = lapply(DEGs, function(i)
  logcpm[na.omit(i$gene_name),])
heatmap = function(data, samples, rows, title){
  as.ggplot(pheatmap(data[c(1:min(rows, 100)),samples], scale = "row", 
                     clustering_method = "complete",
                     display_numbers = F,
                     border_color = NA, cluster_cols = T, cutree_cols = 2, cutree_rows = 2, show_rownames = T,
                     #annotation_col = annotation, show_rownames = F, annotation_names_col = F,
                     #annotation_row = setNames(data.frame(Cluster = as.factor(cutree(q$tree_row, k=2))), "Cluster"), 
                     annotation_names_row = F, main = title,
                     legend_labels = F,))
}

heatplots = lapply(names(DEG_selection[unlist(lapply(DEG_selection, is.matrix))]), function(name) {
  string = str_split_1(name, pattern = "-")
  cols = grep(paste0("\\b", string, "\\b", collapse = "|"), y$samples$group)
  heatmap(DEG_selection[[name]], cols, nrow(DEG_selection[[name]]), name)
  
})

combined_heatplot <- wrap_plots(heatplots)
print(combined_heatplot)
