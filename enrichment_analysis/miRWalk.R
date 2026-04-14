setwd ("~/DE_miRNAs/miRWalk")


library (dplyr)
library (ReactomePA)
require("enrichplot")
library(org.Hs.eg.db)
require("clusterProfiler")

### mirWalk miR-9-5p

df1 <- read.csv ("miRWalk_miRNA_Targets_miR-9-5p.csv")

df1_binding95 <- filter(df1, bindingp > 0.95)

df2 <- filter (df1_binding95, !validated == "")

df2 <- df2 %>% distinct(genesymbol)

df3 <- read.csv ("miRWalk_miRNA_Targets_miR-9-5p.csv")
df3_1 <- filter(df3, bindingp == 1, miRDB ==1)
df3_2 <- filter(df3, bindingp == 1, TargetScan ==1)
df3_3 <- filter (df3, !validated == "")
combined_df3 <- rbind(df3_1, df3_2, df3_3)
df3f <- combined_df3 %>% distinct(genesymbol)
df3 <- df3 %>% distinct(genesymbol)

df4 <- read.csv ("miRWalk_miRNA_Targets_miR-9-5p.csv")
df4 <- filter (df4, !validated == "")
df4 <- df4 %>% distinct(genesymbol)
