library (tidyverse)
library(stringr)

setwd ("~/miRNA_counts/")

average1 <- read.csv("average_counts.csv")

average2 <- select (average1, ReferenceID, Average)

average3 <- arrange(average2, desc(Average))

average3$ReferenceID <- gsub("\\:.*","",average3$ReferenceID)

average4 <- slice(average3, (1:30))

write.csv(average4, "30reads.csv")

library(forcats)

average4 %>%
  ggplot(aes(x = fct_rev(fct_reorder(ReferenceID, Average)), y = Average)) +
  geom_col() + ylab("Read counts average") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab ("")



