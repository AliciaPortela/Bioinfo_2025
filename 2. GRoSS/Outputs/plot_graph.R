setwd("C:/Users/Alicia/Desktop/AntropoGeo/Inputs_GRoSS")
library(igraph)

# Cargar como lista de aristas
g <- read_graph("human.graph", format = "ncol")  # O "edgelist"

# Plot
plot(g, vertex.label.cex = 0.8, edge.arrow.size = 0.4)

library(ape)

tree <- read.tree("human.graph")
plot(tree)

lines <- readLines("human.graph")

# Extrae las aristas (líneas que comienzan con 'edge')
edges <- grep("^edge", lines, value = TRUE)
edge_df <- do.call(rbind, strsplit(edges, "\t"))
edge_df <- edge_df[, 3:4]  # Solo columnas de origen y destino

# Crea el grafo dirigido
g <- graph_from_edgelist(as.matrix(edge_df), directed = TRUE)

plot(g,
     layout = layout_as_tree(g),
     vertex.label.cex = 0.8,
     vertex.size = 15,
     vertex.color = "lightblue",
     edge.arrow.size = 0.4,
     main = "Relaciones Genealógicas de Unidades Evolutivas")