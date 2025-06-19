library(ggplot2)

# Leer tabla de PCA
pca_data <- read.table("PCA_SNPs_under_sel_pop.txt", header=TRUE, stringsAsFactors = FALSE)

# Graficar PC1 vs PC2 por población
ggplot(pca_data, aes(x=PC1, y=PC2, color=POPULATION)) +
  geom_point(size=3) +
  theme_minimal() +
  labs(title="PCA all data", x="PC1", y="PC2") +
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=10))
