setwd ("~/Desktop/excerpt")

unicount <- read.csv("excerpt_reads.csv", row.names=1)

BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("EnhancedVolcano")
BiocManager::install("ggtree")

install.packages("remotes")
library("remotes")
remotes::install_github('YuLab-SMU/ggtree')


#######################
# CHAMAR BIBLIOTECAS #
#######################
library("DESeq2")
library ("ggplot2")
library("ggfortify")
#library('EnhancedVolcano')
library("ggrepel")
library("dplyr")
library(forcats)

#################################
# PREPARACAO DAS TABELAS INPUT #
#################################
# Retirada da coluna NAuni
uni_uni <- select (unicount, -"NAuni")

#Preparacao da tabela de readcounts (countData)
mRC_uni <- as.matrix(unicount)
mRC_uni <- as.data.frame (mRC_uni)

keep_miRNA <- rowSums(counts(count.data.set) >= 10) >= 5

filtered_df <- mRC_uni[rowSums(mRC_uni >= 10) >= 5, ]

meta_data <- read.csv("cat_Tratamento.csv")


colnames(mRC_uni) <- gsub("TotalCount", "", colnames(mRC_uni))
colnames(mRC_uni)



#################################
# condition biperiden control #
###################################

count.data.set <- DESeqDataSetFromMatrix(countData=filtered_df, 
                                         colData=meta_data, design= ~ Treatment)

keep_miRNA <- rowSums(counts(count.data.set) >= 10) >= 5
count.data.set <- count.data.set[keep_miRNA,]


count.data.set.object <- DESeq(count.data.set)

## estimating size factors

vsd <- varianceStabilizingTransformation(count.data.set.object)
norm.data = assay(vsd)
normalized <- as.data.frame(norm.data)


head(norm.data)

write.table(norm.data, sep="\t",file="Norm_data_all_genes_NO_counts_cut_off.txt", row.names=TRUE,col.names=NA,quote=FALSE)
sampleDists <- dist(t(norm.data),  method = "euclidean")

reversed_rows_columns = (t(norm.data))
reversed_rows_columns[1:5,1:5]
sampleDists
clusters=hclust(sampleDists)
plot(clusters)
plotPCA(vsd, intgroup=c("Treatment")) 
#plotPCA(vsd, intgroup=c("sex")) 
#plotPCA(vsd, intgroup=c("batch"))
#plotPCA(vsd, intgroup=c("PTE")) 


res <- results(count.data.set.object, contrast=c("Treatment", "Biperiden", "Placebo"))


res
summary(res)
res = na.omit(res)
head(res)
resfiltered = res[res$pvalue <= 0.05,]
head(resfiltered)
summary(resfiltered)
sum(resfiltered$padj > 0.05)

res.filtered.ordered = resfiltered[order(resfiltered$pvalue),]


head(res.filtered.ordered)
write.table(res.filtered.ordered, sep="\t",file="Results_filtered_ordered.txt", row.names=TRUE,col.names=NA,quote=FALSE)


require(ggplot2)

resOrdFilt.data.frame = as.data.frame(res.filtered.ordered)
head(resOrdFilt.data.frame)
ggplot(resOrdFilt.data.frame, aes(x=log2FoldChange)) +
  geom_histogram()
ggplot(resOrdFilt.data.frame, aes(x=log2FoldChange)) +
  geom_histogram(bins=100) +
  scale_x_continuous( name ="log2FoldChange",breaks=c(-5,-1,-0.5,0, 0.5,1,5), limits = c(-10,10))
resOrdFiltLFC = subset(resOrdFilt.data.frame,log2FoldChange >0.5 | log2FoldChange < -0.5)
dim(resOrdFiltLFC)
ggplot(resOrdFiltLFC, aes(x=log2FoldChange)) +
  geom_histogram(bins=100) +
  scale_x_continuous( name ="log2FoldChange",breaks=c(-5,-1,-0.5,0, 0.5,1,5), limits = c(-10,10))
write.table(resOrdFiltLFC, sep="\t",file="Results_filtered_ordered_LFC.txt", row.names=TRUE,col.names=NA,quote=FALSE)
summary(res.filtered.ordered)
Sign_genes_up = subset(resOrdFilt.data.frame, log2FoldChange >0.5)
Sign_genes_down = subset(resOrdFilt.data.frame, log2FoldChange <(-0.5))













#PCAtools: Outras informações sobre as componentes principais

a <- as.data.frame(norm.data)

a <- a %>% select(order(colnames(a)))

library (tidyverse)

b <- a %>%
  add_column(add_column = "constant_value")

row.names(a) <- unique(row.names(a))

row.names(a) <- gsub("\\:.*","", row.names(a))


meta.data <-  read.csv("meta_data_3.csv", row.names = 1)


p <- pca(a, metadata = meta.data, removeVar = 0.1)

screeplot(p, 1:8,axisLabSize = 18, titleLabSize = 22)


## Nota: O resultado mostra que o PCA (PC1xPC2) não serve muito pra analisar essas amostras, por ter muita variação explicadas pelas PCs após PC2
## Além das baixas explained variation de PC1


### Métodos para ver quantas PCs manter

horn <- parallelPCA(norm.data)
horn$n
# resultado: 6

elbow <- findElbowPoint(p$variance)
elbow
# resultado: 7


screeplot(p,  components = getComponents(p, 1:15),
          vline = c(horn$n, elbow)) +
  
  geom_label(aes(x = horn$n, y = 50,
                 label = 'Horn\'s', vjust = -4, size = 8)) +
  geom_label(aes(x = elbow, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))

## biplot - ver os principais genes e que determinam as coordenadas
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)


## biplot - ver os principais genes e que determinam as coordenadas - Edição das componentes
biplot(p, showLoadings = TRUE, x = 'PC2', y = 'PC3',
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5) 

#### Nota: O PC1 ficou menor que o mostrado pelo PlotPCA do DESeq2.


pairsplot(p)
plotloadings(p,1:6, labSize = 3)








require("EnhancedVolcano")

res <- results(count.data.set.object, contrast=c("condition", "Biperiden", "Control"))


res
summary(res)
res = na.omit(res)
head(res)
resfiltered = res[res$pvalue <= 0.05,]
head(resfiltered)
summary(resfiltered)
sum(resfiltered$pvalue > 0.05)

res.filtered.ordered = resfiltered[order(resfiltered$pvalue),] 
head(res.filtered.ordered)
write.table(res.filtered.ordered, sep="\t",file="Results_filtered_ordered.txt", row.names=TRUE,col.names=NA,quote=FALSE)


require(ggplot2)

resOrdFilt.data.frame = as.data.frame(res.filtered.ordered)

head(resOrdFilt.data.frame)
EnhancedVolcano(resOrdFilt.data.frame,
                
                lab = rownames(resOrdFilt.data.frame),
                
                x = "log2FoldChange",
                y = "padj",
                FCcutoff = 1,
                pCutoff = 0.05)

# write.table(Sign_genes_up, sep="\t",file="UpBiperiden.csv", row.names=TRUE,col.names=NA,quote=FALSE)
# write.table(Sign_genes_down, sep="\t",file="DownBiperiden.csv", row.names=TRUE,col.names=NA,quote=FALSE)





####################################
# ANALISE  DE EXPRESSAO DIFERENCIAL #
####################################
#DESeq2 quick start -> design da analise
dds_uni <- DESeqDataSetFromMatrix(countData = mRC_uni,
                            colData = meta_data,
                            design = ~ condition)


#Selecionar apenas miRNAs com pelo menos 5 amostras com no minimo 10 reads 
keep_uni <- rowSums(counts(dds_uni) >= 10) >= 5
dds_uni <- dds_uni[keep_uni,]

#Extrair a matriz de reads normalizada
dds_uni <- estimateSizeFactors(dds_uni)
norm_mRC_uni <- counts(dds_uni,normalized = T)

#Renomear as linhas e colunas pra ficar mais bonito
#rownames(norm_mRC_uni) <- gsub("[:MIMAT].*", "", rownames(norm_mRC_uni))

ntd_uni <- normTransform(dds_uni)

#Definir "Control" como ref para comparacao
dds_uni$condition <- relevel(dds_uni$condition, ref = "Placebo")

# lists the coefficients
dds_uni <- DESeq(dds_uni)
resultsNames(dds_uni)

a <- as.data.frame(norm_mRC_uni)
write.csv (a, "normalized_counts_DESeq2.csv")

#############################
### COMPARACAO Biperiden x Control #####
#############################
Bip_Ctrl_uni <- results(dds_uni, name="condition_Biperiden_vs_Placebo")
Bip_Ctrl_uni <- as.data.frame(Bip_Ctrl_uni)
summary(Bip_Ctrl_uni)

##VOLCANO PLOT##
#Criar coluna com nome de miRNAs
Bip_Ctrl_uni$miRNA <- rownames(Bip_Ctrl_uni)
Bip_Ctrl_uni$miRNA <- sub("[:MIMAT].*", "",Bip_Ctrl_uni$miRNA)


#Retirada de NAs
sBip_Ctrl_uni <- Bip_Ctrl_uni[complete.cases(Bip_Ctrl_uni), ]

#Cria??o de coluna de up ou down regulated
sBip_Ctrl_uni$diffexpressed <- "No significant\ndifferences"
sBip_Ctrl_uni$diffexpressed[sBip_Ctrl_uni$log2FoldChange > 1 & sBip_Ctrl_uni$padj < 0.05] <- "Up regulated"
sBip_Ctrl_uni$diffexpressed[sBip_Ctrl_uni$log2FoldChange < -1 & sBip_Ctrl_uni$padj < 0.05] <- "Downregulated \nin biperiden-treated\npatients"

sBip_Ctrl_uni$delabel <- NA
sBip_Ctrl_uni$delabel[sBip_Ctrl_uni$diffexpressed != "No difference"] <- sBip_Ctrl_uni$miRNA[sBip_Ctrl_uni$diffexpressed != "Sem diferença significativa"]


#Grafico
ggplot(sBip_Ctrl_uni, aes(x=log2FoldChange, y=-log10(pvalue), col=str_wrap(diffexpressed,20), label=delabel)) + 
  geom_point(aes(color = diffexpressed),size=2.6) + theme(legend.position = "right") + 
  theme_minimal() + 
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.001), col="black") + xlim (-5,5) + 
  scale_color_manual(values= c("purple", "black", "coral1")) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        axis.text=element_text(size=20),axis.title=element_text(size=17)) +
  labs(color='') + guides(colour = guide_legend(override.aes = list(size=3)))

ggplot(sBip_Ctrl_uni, aes(x=log2FoldChange, y=-log10(padj), col=str_wrap(diffexpressed,20))) + 
  geom_point(aes(color = diffexpressed),size=2.6) + theme(legend.position = "right") + 
  theme_minimal() + 
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black") + xlim (-5,5) + 
  scale_color_manual(values= c("purple", "black", "coral1")) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
        axis.text=element_text(size=20),axis.title=element_text(size=17)) +
  labs(color='') + guides(colour = guide_legend(override.aes = list(size=3)))

fig2a <- ggplot(sBip_Ctrl_uni, aes(x=log2FoldChange, y=-log10(pvalue), col=str_wrap(diffexpressed,20), label=delabel)) + 
  geom_point(aes(color = diffexpressed),size=2.6) + theme(legend.position = "right") + 
  theme_minimal() + 
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.001), col="black") + xlim (-5,5) + 
  scale_color_manual(values= c("purple", "black", "coral1")) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        axis.text=element_text(size=20),axis.title=element_text(size=17)) +
  labs(color='') + guides(colour = guide_legend(override.aes = list(size=3)))

dev.size("in")

ggsave("Volcano_Bip.jpeg", fig2a, width = 7.96, height = 6.8, dpi = 300)

write.csv(sBip_Ctrl_uni, "DE_Bip.csv")



## SALVAR OS RESULTADOS EM CSV ##
write.csv(Bip_Ctrl_uni, "BiperidenXCtrlDEG.csv", col.names=TRUE, row.names = TRUE)



#####################################
##Violin Plot - miRNAs Dif. expr.##



#miRNA 9-5p#
norm_mRCt_uni <- as.data.frame (t(norm_mRC_uni))
colnames(norm_mRCt_uni) <- gsub("[-]", "_", colnames(norm_mRCt_uni))
colnames(norm_mRCt_uni) <- gsub ("\\:.*","",colnames(norm_mRCt_uni))

miR <- select(norm_mRCt_uni, hsa_miR_9_5p)
miR <- arrange(miR , row.names.data.frame(miR))

Group = meta_data$condition

vp <- cbind( "miRNA" = miR$hsa_miR_9_5p, "group" = Group)
vp <- as.data.frame(vp)
vp$mimi <- as.numeric (vp$miRNA)



write.csv(vp, "mir-9-5p.csv")



ggplot(vp, aes(x=(vp$group), y=vp$mimi, color=group)) + 
  geom_violin(scale="width") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("miRNA counts") + xlab ("Grupo") + ggtitle("miRNA mir-9-5p") +
  theme(plot.title = element_text(hjust = 0.5))



ggplot(vp, aes(x=fct_rev(group), y=mimi, color=group)) + 
  geom_violin(scale='width')+
  theme_bw() +
  geom_dotplot(binaxis='y', stackdir='center') +
  scale_color_manual(values=c("grey1","purple")) +
  ylab("Normalized miR-9-5p\nread counts") + xlab ("") + labs(color='')+ ylim(0,220)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=17)) +
  guides(fill='none',color='none')

dev.size("in")

fig2b <- ggplot(vp, aes(x=fct_rev(group), y=mimi, color=group)) + 
  geom_violin(scale='width')+
  theme_bw() +
  geom_dotplot(binaxis='y', stackdir='center') +
  scale_color_manual(values=c("grey1","purple")) +
  ylab("Normalized miR-9-5p\nread counts") + xlab ("") + labs(color='')+ ylim(0,220)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=17)) +
  guides(fill='none',color='none')



ggsave("violin_plot.jpeg", fig2b, width = 7.96, height = 6.8, dpi = 300)


##################


ggplot(vp, aes(x=(vp$group), y=vp$mimi, color=group)) + 
  geom_violin() +
  geom_violin(trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("miRNA counts") + xlab ("Grupo") + ggtitle("miRNA mir-9-5p") +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(vp, aes(x=fct_rev(group), y=mimi, color=group)) + 
  geom_violin(scale='width')+
  theme_bw() +
  geom_dotplot(binaxis='y', stackdir='center') +
  scale_color_brewer(palette="PuRd") +
  ylab("Normalized miR-9-5p\nread counts") + xlab ("") + labs(color='')+ ylim(0,220)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=17),) +
  guides(fill='none',color='none')

dev.size("in")

fig2b <- ggplot(vp, aes(x=fct_rev(group), y=mimi, color=group)) + 
  geom_violin(scale='width')+
  theme_bw() +
  geom_dotplot(binaxis='y', stackdir='center') +
  scale_color_brewer(palette="PuRd") +
  ylab("Normalized miR-9-5p\nread counts") + xlab ("") + labs(color='')+ ylim(0,220)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=17)) +
  guides(fill='none',color='none')


ggsave("violin_plot.jpeg", fig2b, width = 7.96, height = 6.8, dpi = 300)


#miRNA 17-5p#
Group = meta_data$condition
mir <- norm_mRCt_uni$hsa_miR_17_5p
vp <- cbind( "miRNA" = mir, "group" = Group)
vp <- as.data.frame(vp)
vp$mimi <- as.numeric (vp$miRNA)

write.csv(vp, "mir-17-5p.csv")

ggplot(vp, aes(x=vp$group, y=vp$mimi, color=group)) + 
  geom_violin() +
  geom_violin(trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("miRNA counts") + xlab ("Grupo") + ggtitle("miRNA mir-17-5p") +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(vp, aes(x=vp$group, y=vp$mimi, color=group)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("counts normalizados") + xlab ("Grupo") + ggtitle("miRNA mir-17-5p") +
  theme(plot.title = element_text(hjust = 0.5))

p <- ggplot(vp, aes(x=vp$group, y=vp$mimi, color=group)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("counts normalizados") + xlab ("Grupo") + ggtitle("miRNA mir-17-5p") +
  theme(plot.title = element_text(hjust = 0.5))


library(forcats)
p + aes(x =forcats::fct_rev(reorder(group,group)))


#miRNA 3615#
norm_mRCt_uni <- as.data.frame (t(norm_mRC_uni))
colnames(norm_mRCt_uni) <- gsub("[-]", "_", colnames(norm_mRCt_uni))
Group = meta_data$condition
mir <- norm_mRCt_uni$hsa_miR_3615
vp <- cbind( "miRNA" = mir, "group" = Group)
vp <- as.data.frame(vp)
vp$mimi <- as.numeric (vp$miRNA)

ggplot(vp, aes(x=vp$group, y=vp$mimi, color=group)) + 
  geom_violin() +
  geom_violin(trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("miRNA counts") + xlab ("Grupo") + ggtitle("miRNA mir-3615") +
  theme(plot.title = element_text(hjust = 0.5))



ggplot(vp, aes(x=vp$group, y=vp$mimi, color=group)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("miRNA counts") + xlab ("Grupo") + ggtitle("miRNA mir-3615") +
  theme(plot.title = element_text(hjust = 0.5))




##Heatmap##
library("pheatmap")
select_uni <- order(rowMeans(counts(dds_uni,normalized=TRUE)),
                    decreasing=TRUE)[1:50]

df_uni <- as.data.frame(colData(dds_uni)[,c("sex","condition")])
df_uni2 <- select(df_uni, condition)


#rownames(df_uni) <- gsub("[:MIMAT].*", "", rownames(df_uni))
#rownames(ntd_uni) <- gsub("[:MIMAT].*", "", rownames(ntd_uni))
#colnames(ntd_uni) <- gsub("[:Total].*","", colnames(ntd_uni))

pheatmap(assay(ntd_uni)[select_uni,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df_uni2, fontsize_row=8, fontsize_col
         =8, main="Heatmap - top 50 expressed miRNA")



############################################################
############################################################
############################################################

## GRAFICOS
#PCA -> usando o valor dos reads absolutos
rownames(sdesign_uni) <- gsub("[:Total].*","", rownames(sdesign_uni))

pca_condition_uni <- counts(dds_uni, normalized=T)
pca_cond_uni <- prcomp(t(pca_condition_uni))

## para calcular a % que o PC representa:
pca_cond.var_uni <- pca_cond_uni$sdev^2
pca_cond.var.per_uni <- round(pca_cond.var_uni/sum(pca_cond.var_uni)*100,1)

### PCA - Batch effect ###
pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,2])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$batch)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Batch") +
  labs(title = "PCA - Batch Effect",
       subtitle = "Plot of PC1 x PC2, coloured by Sequencing Batch") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,2], Y=pca_cond_uni$x[,3])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$batch)) +
  geom_text() +
  xlab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Batch") +
  labs(title = "PCA - Batch Effect",
       subtitle = "Plot of PC2 x PC3, coloured by Sequencing Batch") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,3], Y=pca_cond_uni$x[,4])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$batch)) +
  geom_text() +
  xlab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) +
  ylab(paste("PC4 - ", pca_cond.var.per_uni[4], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Batch") +
  labs(title = "PCA - Batch Effect",
       subtitle = "Plot of PC3 x PC4, coloured by Sequencing Batch") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,3])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$batch)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Batch") +
  labs(title = "PCA - Batch Effect",
       subtitle = "Plot of PC1 x PC3, coloured by Sequencing Batch") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

### PCA - SEXO ###
pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,2])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$sex)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Sex") +
  labs(title = "PCA - Sex",
       subtitle = "Plot of PC1 x PC2, coloured by Sex") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,2], Y=pca_cond_uni$x[,3])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$sex)) +
  geom_text() +
  xlab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Sex") +
  labs(title = "PCA - Sex",
       subtitle = "Plot of PC2 x PC3, coloured by Sex") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,3], Y=pca_cond_uni$x[,4])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$sex)) +
  geom_text() +
  xlab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) +
  ylab(paste("PC4 - ", pca_cond.var.per_uni[4], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Sex") +
  labs(title = "PCA - Sex",
       subtitle = "Plot of PC3 x PC4, coloured by Sex") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,3])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$sex)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Sex") +
  labs(title = "PCA - Sex",
       subtitle = "Plot of PC1 x PC3, coloured by Sex") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

### PCA - CONDITION ###

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,2])
pca_cond <- select (pca_cond.data_uni, X, Y)


ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$condition)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Condition") +
  labs(title = "PCA - Condition",
       subtitle = "Plot of PC1 x PC2, coloured by Condition") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))




pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,2], Y=pca_cond_uni$x[,3])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$condition)) +
  geom_text() +
  xlab(paste("PC2 - ", pca_cond.var.per_uni[2], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(title = "PCA - Condition",
       subtitle = "Plot of PC2 x PC3, coloured by Condition") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,3], Y=pca_cond_uni$x[,4])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$condition)) +
  geom_text() +
  xlab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) +
  ylab(paste("PC4 - ", pca_cond.var.per_uni[4], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Condition") +
  labs(title = "PCA - Condition",
       subtitle = "Plot of PC3 x PC4, coloured by Condition") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

pca_cond.data_uni <- data.frame(Sample = rownames(pca_cond_uni$x), X=pca_cond_uni$x[,1], Y=pca_cond_uni$x[,3])
ggplot(data = pca_cond.data_uni, aes(x = X, y = Y, label = Sample, color = meta_data$condition)) +
  geom_text() +
  xlab(paste("PC1 - ", pca_cond.var.per_uni[1], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_cond.var.per_uni[3], "%", sep = "")) + 
  theme_bw() +
  labs(colour = "Condition") +
  labs(title = "PCA - Condition",
       subtitle = "Plot of PC1 x PC3, coloured by Condition") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))


#############################################
# BATCH #
#########################################


count.data.set <- DESeqDataSetFromMatrix(countData=mRC_uni, 
                                         colData=meta_data, design= ~ batch)

count.data.set.object <- DESeq(count.data.set)

## estimating size factors

vsd <- varianceStabilizingTransformation(count.data.set.object)
norm.data = assay(vsd)

head(norm.data)

sampleDists <- dist(t(norm.data),  method = "euclidean")

reversed_rows_columns = (t(norm.data))
reversed_rows_columns[1:5,1:5]
sampleDists
plotPCA(vsd, intgroup=c("batch")) 


res <- results(count.data.set.object, contrast=c("batch", "1", "2"))


res
summary(res)
res = na.omit(res)
head(res)
resfiltered = res[res$padj <= 0.05,]
head(resfiltered)
summary(resfiltered)
sum(resfiltered$padj > 0.05)


#############################################
# SEX #
#########################################


count.data.set <- DESeqDataSetFromMatrix(countData=mRC_uni, 
                                         colData=meta_data, design= ~ sex)

count.data.set.object <- DESeq(count.data.set)

## estimating size factors

vsd <- varianceStabilizingTransformation(count.data.set.object)
norm.data = assay(vsd)

head(norm.data)

sampleDists <- dist(t(norm.data),  method = "euclidean")

reversed_rows_columns = (t(norm.data))
reversed_rows_columns[1:5,1:5]
sampleDists
clusters=hclust(sampleDists)
plot(clusters)
plotPCA(vsd, intgroup=c("sex")) 


res <- results(count.data.set.object, contrast=c("sex", "M", "F"))


res
summary(res)
res = na.omit(res)
head(res)
resfiltered = res[res$padj <= 0.05,]
head(resfiltered)
summary(resfiltered)
sum(resfiltered$padj > 0.05)


#############################################
# Alcohol #
#########################################


count.data.set <- DESeqDataSetFromMatrix(countData=mRC_uni, 
                                         colData=meta_data, design= ~ Alcohol)

count.data.set.object <- DESeq(count.data.set)

## estimating size factors

vsd <- varianceStabilizingTransformation(count.data.set.object)
norm.data = assay(vsd)

head(norm.data)

#write.table(norm.data, sep="\t",file="Norm_data_all_genes_NO_counts_cut_off.txt", row.names=TRUE,col.names=NA,quote=FALSE)
sampleDists <- dist(t(norm.data),  method = "euclidean")

reversed_rows_columns = (t(norm.data))
reversed_rows_columns[1:5,1:5]
sampleDists
clusters=hclust(sampleDists)
plot(clusters)
plotPCA(vsd, intgroup=c("Alcohol")) 


res <- results(count.data.set.object, contrast=c("Alcohol", "Yes", "No"))


res
summary(res)
res = na.omit(res)
head(res)
resfiltered = res[res$padj <= 0.05,]
head(resfiltered)
summary(resfiltered)
sum(resfiltered$padj > 0.05)



#############################################
# DEATH #
#########################################


count.data.set <- DESeqDataSetFromMatrix(countData=mRC_uni, 
                                         colData=meta_data, design= ~ Death)

count.data.set.object <- DESeq(count.data.set)

## estimating size factors

vsd <- varianceStabilizingTransformation(count.data.set.object)
norm.data = assay(vsd)

head(norm.data)

write.table(norm.data, sep="\t",file="Norm_data_all_genes_NO_counts_cut_off.txt", row.names=TRUE,col.names=NA,quote=FALSE)
sampleDists <- dist(t(norm.data),  method = "euclidean")

reversed_rows_columns = (t(norm.data))
reversed_rows_columns[1:5,1:5]
sampleDists
clusters=hclust(sampleDists)
plot(clusters)
plotPCA(vsd, intgroup=c("Death")) 


res <- results(count.data.set.object, contrast=c("Death", "Dead", "Alive"))


res
summary(res)
res = na.omit(res)
head(res)
resfiltered = res[res$padj <= 0.05,]
head(resfiltered)
summary(resfiltered)
sum(resfiltered$padj > 0.05)





#############################################
# PTE #
#########################################

count.data.set <- DESeqDataSetFromMatrix(countData=mRC_uni, 
                                         colData=meta_data, design= ~ PTE)

count.data.set.object <- DESeq(count.data.set)

## estimating size factors

vsd <- varianceStabilizingTransformation(count.data.set.object)
norm.data = assay(vsd)

head(norm.data)

write.table(norm.data, sep="\t",file="Norm_data_all_genes_NO_counts_cut_off.txt", row.names=TRUE,col.names=NA,quote=FALSE)
sampleDists <- dist(t(norm.data),  method = "euclidean")

reversed_rows_columns = (t(norm.data))
reversed_rows_columns[1:5,1:5]
sampleDists
clusters=hclust(sampleDists)
plot(clusters)
plotPCA(vsd, intgroup=c("PTE")) 


nd <- as.data.frame(norm.data)


res <- results(count.data.set.object, contrast=c("PTE", "With_PTE", "Without_PTE"))


res
summary(res)
res = na.omit(res)
head(res)
resfiltered = res[res$pvalue <= 0.05,]
head(resfiltered)
summary(resfiltered)
sum(resfiltered$padj > 0.05)

res.filtered.ordered = resfiltered[order(resfiltered$padj),] 
head(res.filtered.ordered)
resOrdFilt.data.frame = as.data.frame(res.filtered.ordered)

#############################################
# Age_group #
#########################################

count.data.set <- DESeqDataSetFromMatrix(countData=mRC_uni, 
                                         colData=meta_data, design= ~ Age_group)

count.data.set.object <- DESeq(count.data.set)

## estimating size factors

vsd <- varianceStabilizingTransformation(count.data.set.object)
norm.data = assay(vsd)

head(norm.data)

#write.table(norm.data, sep="\t",file="Norm_data_all_genes_NO_counts_cut_off.txt", row.names=TRUE,col.names=NA,quote=FALSE)
sampleDists <- dist(t(norm.data),  method = "euclidean")

reversed_rows_columns = (t(norm.data))
reversed_rows_columns[1:5,1:5]
sampleDists
clusters=hclust(sampleDists)
plot(clusters)
plotPCA(vsd, intgroup=c("Age_group")) 

plotPCA(vsd, intgroup=c("Age")) 

pca <- plotPCA(vsd, intgroup=c("Age")) 

