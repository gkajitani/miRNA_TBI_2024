setwd ("~/DE_miRNAs/miRWalk")

library (dplyr)
library (ReactomePA)
require("enrichplot")
library(org.Hs.eg.db)
require("clusterProfiler")
library(DOSE)
library (ggplot2)
library (remotes)
library (multienrichjam)
library(igraph)

options(enrichplot.colours = c("red","blue"))

### mirWalk mir-9-5p

df1 <- read.csv ("miRWalk_miRNA_Targets_miR-9-5p.csv")

df1_binding95 <- filter(df1, bindingp > 0.95)


df2 <- filter (df1_binding95, !validated == "")

df2 <- df2 %>% distinct(genesymbol)

df3 <- filter(df1, !validated=="")

df3 <- df3 %>% distinct(genesymbol)

write.csv(df2,"list_miR_9_5p.csv")

################


GOBP <- enrichGO(gene         = df2$genesymbol,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

cluster_summary_GOBP_miR_9_5p <- data.frame(GOBP)
df <- as.data.frame(GOBP)
GOBP_f <- gsfilter(GOBP, by = 'Count', min = 4)
dfGOBP <- as.data.frame(GOBP_f)



# write.csv (dfGOBP, "GOBP_4.csv")

#filtered_df <- df[!grepl("3/71", df$GeneRatio), ]
#filtered_df <- filtered_df[!grepl("2/71", filtered_df$GeneRatio), ]
#df2 <- enrichDF2enrichResult(filtered_df)

g <- cnetplot(GOBP, showCategory = 10, layout = "kk" , categorySize="geneNum", cex_label_category=1, cex_label_gene=0.8, categoryColor="pvalue")

g  

c <- g + scale_color_manual(values = c("black","purple")) 
g2 <- g + scale_colour_gradientn(name = "p.value", colours = c("blue","white","red"))
g2
g

c

edo <- pairwise_termsim(GOBP)
edo2 <- pairwise_termsim(GOBP_f)

emapplot(edo, ShowCategory=2, cluster.params = list(n = 3),layout="nicely", cex_label_category=0.7)
emapplot_cluster(edo, cluster.params = list(n = 3, color = NULL, label_words_n = 4, label_format = 30))

treeplot(edo, cluster.params = list(n = 4), nWords=4)


emapplot(edo2, ShowCategory=1, layout="nicely", cex_label_category=0.7, cluster.params = list(n = 4))
emapplot_cluster(edo2, cluster.params = list(n = 3), nWords=4)

selected_pathways_cnet <- c("response to decreased oxygen levels", "cellular response to chemical stress",
                            "negative regulation of neurogenesis", "negative regulation of nervous system development", 
                         "cellular response to environmental stimulus", "positive regulation of neuron death")

fold_change_geneList <- read.csv("adjusted_p.csv")
head(fold_change_geneList)

cplot <- cnetplot (GOBP, showCategory = selected_pathways_cnet, font.size=22, foldChange = fold_change_geneList, categorySize="pvalue")
cplot <- cnetplot (GOBP, showCategory = selected_pathways_cnet, font.size=22)

cplot
cpl <- cplot+scale_color_manual(values = c("black","purple"))

cpl

colorPalette <- function(p_adjusted) {
  # Replace 'p_adjusted' with the vector of your p-adjusted values
  colors <- colorRampPalette(c("blue", "red"))(length(p_adjusted))
  return(colors)
}
category_colors <- colorPalette(p_adjusted)

# Extract the ggplot object from the cnetplot output
plot_object <- cnet$plot

# Modify the plot aesthetics to use the category colors
plot_object <- plot_object +
  scale_fill_manual(values = category_colors)  # Modify the fill colors based on category_colors

# Print the modified plot
print(plot_object)


cpl

ggsave(cpl, filename = "cnetplot.png", bg = 'white',width=16, height=9.8)


selected_pathways_tree <- c(
  "response to hypoxia", "response to decreased oxygen levels",
"cellular response to chemical stress", 
  "response to amyloid-beta", 
  "regulation of cysteine-type endopeptidase activity involved in apoptotic process", 
  "negative regulation of neurogenesis",
  "negative regulation of nervous system development", 
  "cellular response to environmental stimulus", "miRNA metabolic process", "positive regulation of neuron death", 
  "cellular response to peptide", "positive regulation of endopeptidase activity", "regulation of neurogenesis",
  "negative regulation of cell development", "cellular response to oxidative stress",
  "regulation of endopeptidase activity"
  )
treeplot(edo2, cluster.params = list(n = 4), nWords=0, showCategory=selected_pathways_tree)
t <- treeplot(edo2, cluster.params = list(n = 4), nWords=0, showCategory=selected_pathways_tree)
ggsave (t, filename = "treeplot.png", bg= 'white', width = 10.6, height=9)


selected_pathways <- c("response to amyloid-beta", "myelination in peripheral nervous system", "negative regulation of neurogenesis", 
                       "positive regulation of neuron death", "negative regulation of gliogenesis",
                       "peripheral nervous system axon regeneration","Schwann cell proliferation","anterograde axonal transport")



barplot(GOBP,  showCategory = selected_pathways, font.size=22) 
b <- barplot(GOBP,  showCategory = selected_pathways, font.size=22) 
ggsave(b, filename = "barplot.png", bg = 'white',width=14.5, height=9)


write.table(cluster_summary_GOBP_miR_9_5p , file = "GOBP_miR_9_5p_miRWalk.csv", row.names=FALSE)

##Get the Entrez gene IDs associated with those symbols
EG_ID = mget(df2$genesymbol, revmap(org.Hs.egSYMBOL),ifnotfound=NA)

KEGG <- enrichKEGG(gene         = EG_ID,
                   organism = "hsa",
                   keyType       = 'kegg',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
)

KEGG_cluster_summary_miR_9_5p <- data.frame(KEGG)
write.table(KEGG_cluster_summary_miR_9_5p , file = "KEGG_miR_9_5p_miRWalk.csv", row.names=FALSE)


reactome <- enrichPathway(EG_ID,
                          organism = "human",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
)

reactome_summary_miR_9_5p <- data.frame (reactome)
write.table(reactome_summary_miR_9_5p , file = "Reactome_miR_9_5p_miRWalk.csv", row.names=FALSE)



############
# miR-3615 #
############

### mirWalk mir-9-5p

df1 <- read.csv ("miRWalk_miRNA_Targets_miR_3615.csv")

df1_binding95 <- filter(df1, bindingp > 0.95)


df2 <- filter (df1_binding95, !validated == "")

df2 <- df2 %>% distinct(genesymbol)


column_name <- "genesymbol"

# Keep only the first occurrence of each value in the specified column
df_no_duplicates <- df2[!duplicated(df2[[column_name]]), ]

write.csv (df_no_duplicates,"list_miR_3615_complete.csv")

df3 <- filter(df1, !validated=="")

df3 <- df3 %>% distinct(genesymbol)



write.csv(df2,"list_miR_3615.csv")



GOCC <- enrichGO(gene         = df2$genesymbol,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

cluster_summary_GOBP_miR_3615 <- data.frame(GOCC)
df <- as.data.frame(GOCC)
GOCC_f <- gsfilter(GOCC, by = 'Count', min = 4)
dfGOCC <- as.data.frame(GOCC_f)




GOBP <- enrichGO(gene         = df3$genesymbol,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

cluster_summary_GOBP_miR_3615 <- data.frame(GOBP)
df <- as.data.frame(GOBP)
GOBP_f <- gsfilter(GOBP, by = 'Count', min = 4)
dfGOBP <- as.data.frame(GOBP_f)

##Get the Entrez gene IDs associated with those symbols
EG_ID = mget(df2$genesymbol, revmap(org.Hs.egSYMBOL),ifnotfound=NA)

KEGG <- enrichKEGG(gene         = EG_ID,
                   organism = "hsa",
                   keyType       = 'kegg',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
)

KEGG_cluster_summary_miR_3615 <- data.frame(KEGG)


reactome <- enrichPathway(EG_ID,
                          organism = "human",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
)

reactome_summary_miR_3615 <- data.frame (reactome)

