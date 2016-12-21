library(dplyr)
library(ggplot2)
library(ggrepel)

frame_del <- read.table("omega.inframe_del.count.txt", sep = "\t", header = TRUE) %>%
  mutate(total_num = X00 + X01 + X10 + X11) %>% filter(total_num >= 8)


pca <- prcomp(t(frame_del[,2:5] / rowSums(frame_del[,2:5])))

pca2 <- data.frame(Gene = frame_del$Gene_Symbol, 
                   PC1 = pca$rotation[,1],
                   PC2 = pca$rotation[,2])


ggplot(pca2, aes(x = PC1, y = PC2)) + 
  geom_point(colour = "red") +
  geom_text_repel(data = pca2, 
                  aes(x = PC1, y = PC2, label = Gene)) +
  labs(x = "PC1", y = "PC2") +
  theme_bw() +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)))

