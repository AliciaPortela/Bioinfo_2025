library(ggplot2)

setwd("")
# Leer tabla de PCA
pca_data <- read.table("PCA_data.txt", header=TRUE, stringsAsFactors = FALSE)

# Leer tabla de poblaciones
pop_data <- read.table("popfile.txt", sep = ",", header=TRUE, stringsAsFactors = FALSE)

# Unir tablas por la columna (IID)
pca_pop <- left_join(pca_data, pop_data, by = "IID")

# Graficar PC1 vs PC2 por continente
plot <- ggplot(pca_pop, aes(x=PC1, y=PC2, color=CONTINENT)) +
  geom_point(size=3) +
  theme_minimal() +
  labs(title="PCA_all", x="PC1", y="PC2") +
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=10)) + stat_ellipse()

ggsave("PCA_all.png", plot = plot, width = 12, height = 10, dpi = 300)
