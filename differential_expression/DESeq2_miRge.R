setwd ("~/miRge/")

unicount <- read.csv("miR.Counts.csv", row.names = 1)

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

#Preparacao da tabela de readcounts (countData)
mRC_uni <- as.matrix(unicount)
mRC_uni <- as.data.frame (mRC_uni)

#write.csv(mRC_uni, "4Runs_Reads.csv", row.names = T, col.names = T)

meta_data <- read.csv("meta_data_sort.csv")
# Ver se não é melhor o meta_data_3

colnames(mRC_uni) <- gsub("TotalCount", "", colnames(mRC_uni))

colnames(mRC_uni)

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
b <- as.data.frame(norm_mRC_uni)

# write.csv(b, "normalized_counts.csv")

#Renomear as linhas e colunas pra ficar mais bonito
#rownames(norm_mRC_uni) <- gsub("[:MIMAT].*", "", rownames(norm_mRC_uni))

ntd_uni <- normTransform(dds_uni)

#Definir "Control" como ref para comparacao
dds_uni$condition <- relevel(dds_uni$condition, ref = "Placebo")

# lists the coefficients
dds_uni <- DESeq(dds_uni)
resultsNames(dds_uni)

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

library(stringr)
#Grafico
ggplot(sBip_Ctrl_uni, aes(x=log2FoldChange, y=-log10(pvalue), col=str_wrap(diffexpressed,20))) + 
  geom_point(aes(color = diffexpressed),size=2.6) + theme(legend.position = "right") + 
  theme_minimal() + 
  geom_vline(xintercept=c(-1, 1), col="black") +
  geom_hline(yintercept=-log10(0.001), col="black") + xlim (-5,5) + 
  scale_color_manual(values= c("purple", "black", "coral1")) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="none",
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


fig2a <- ggplot(sBip_Ctrl_uni, aes(x=log2FoldChange, y=-log10(pvalue), col=str_wrap(diffexpressed,20))) + 
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
miR <- select(norm_mRCt_uni, hsa_miR_9_5p)

Group = meta_data$condition

g<-as.data.frame(Group)

vp <- cbind( "miRNA" = miR$hsa_miR_9_5p, "group" = Group)
vp <- as.data.frame(vp)
vp$mimi <- as.numeric (vp$miRNA)



#write.csv(vp, "mir-9-5p.csv")

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
  ylab("Normalized miR-9-5p\nread counts") + xlab ("") + labs(color='')+ ylim(0,1700)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=17)) +
  guides(fill='none',color='none')

dev.size("in")

fig2b <- ggplot(vp, aes(x=fct_rev(group), y=mimi, color=group)) + 
  geom_violin(scale='width')+
  theme_bw() +
  geom_dotplot(binaxis='y', stackdir='center') +
  scale_color_manual(values=c("grey1","purple")) +
  ylab("Normalized miR-9-5p\nread counts") + xlab ("") + labs(color='')+ ylim(0,1700)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=17)) +
  guides(fill='none',color='none')



ggsave("violin_plot.jpeg", fig2b, width = 7.96, height = 6.8, dpi = 300)


#hsa-miR-3615-3p#
norm_mRCt_uni <- as.data.frame (t(norm_mRC_uni))
colnames(norm_mRCt_uni) <- gsub("[-]", "_", colnames(norm_mRCt_uni))
Group = meta_data$condition
mir <- norm_mRCt_uni$hsa_miR_3615_3p

vp <- cbind( "miRNA" = mir, "group" = Group)
vp <- as.data.frame(vp)
vp$mimi <- as.numeric (vp$miRNA)

write.csv(vp, "mir-3615-3p.csv")

ggplot(vp, aes(x=vp$group, y=vp$mimi, color=group)) + 
  geom_violin() +
  geom_violin(trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("miRNA counts") + xlab ("Grupo") + ggtitle("miRNA mir-425-3p") +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(vp, aes(x=vp$group, y=vp$mimi, color=group)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("counts normalizados") + xlab ("Grupo") + ggtitle("miRNA mir-3615-3p") +
  theme(plot.title = element_text(hjust = 0.5))

p <- ggplot(vp, aes(x=vp$group, y=vp$mimi, color=group)) + 
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_color_brewer(palette="Dark2") +
  ylab("counts normalizados") + xlab ("Grupo") + ggtitle("miRNA mir-3615-3p") +
  theme(plot.title = element_text(hjust = 0.5))


library(forcats)
p + aes(x =forcats::fct_rev(reorder(group,group)))


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


meta_data$batch <- as.character(meta_data$batch)

count.data.set <- DESeqDataSetFromMatrix(countData=mRC_uni, 
                                         colData=meta_data, design= ~ condition)



count.data.set.object <- DESeq(count.data.set)

vsd <- varianceStabilizingTransformation(count.data.set.object)
norm.data = assay(vsd)

## GRAFICOS
##########################################################################################

norm.data = assay(vsd)

head(norm.data)

sampleDists <- dist(t(norm.data),  method = "euclidean")

reversed_rows_columns = (t(norm.data))
reversed_rows_columns[1:5,1:5]
sampleDists
plotPCA(vsd, intgroup=c("batch")) 
plotPCA(vsd, intgroup=c("sex")) 
plotPCA(vsd, intgroup=c("PTE")) 
plotPCA(vsd, intgroup=c("condition")) 



## Para maior análise dos PCs, usando o pacote PCAtools

library ("PCAtools")

library(genefilter)
library(ggplot2)
library(ggrepel)
library(ggfortify)


####################################

#PCAtools: Outras informações sobre as componentes principais

p <- pca(norm.data, metadata = meta_data, removeVar = 0.1)

screeplot(p, 1:8,axisLabSize = 18, titleLabSize = 22)


## Nota: O resultado mostra que o PCA (PC1xPC2) não serve muito pra analisar essas amostras, por ter muita variação explicadas pelas PCs após PC2
## Além das baixas explained variation de PC1


### Métodos para ver quantas PCs manter

horn <- parallelPCA(norm.data)
horn$n
# resultado: 8

elbow <- findElbowPoint(p$variance)
elbow
# resultado: 9


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

#######################################
##########################################################################################


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
plotPCA(vsd, intgroup=c("sex")) 
plotPCA(vsd, intgroup=c("PTE")) 
plotPCA(vsd, intgroup=c("condition")) 

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

continuous_var_categories <- cut(continuous_var, breaks = 3) # Adjust breaks as needed

# Extract PCA coordinates
pca_data <- data.frame(pca$rotation, continuous_var_categories)

# Plot PCA using ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = continuous_var_categories)) +
  geom_point() +
  scale_color_manual(values = c("red", "green", "blue")) +  # Specify color scale here
  labs(title = "PCA Plot", x = "PC1", y = "PC2")


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



#################################
# condition biperiden control #
###################################

count.data.set <- DESeqDataSetFromMatrix(countData=mRC_uni, 
                                         colData=meta_data, design= ~ condition)

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
plotPCA(vsd, intgroup=c("condition")) 



res <- results(count.data.set.object, contrast=c("condition", "Biperiden", "Control"))


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



