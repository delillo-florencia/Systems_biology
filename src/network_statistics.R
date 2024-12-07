library(igraph)
library(ggraph)

hemoglobin <- data.frame(
    from = c("ALPHA_GLOBIN", "ALPHA_GLOBIN", "BETA_GLOBIN"),
    to = c("ALPHA_GLOBIN", "BETA_GLOBIN", "BETA_GLOBIN")
    # and any other column 
    )

g <- graph_from_data_frame(hemoglobin,directed=FALSE)

# --------------DEGREE------------
# Average degree
mean(degree(g)) 

# Plot Degree distribution
df <- data.frame(k = degree(g))
ggplot(df, aes(x = k)) +
  geom_density() + 
  geom_vline(xintercept = length(neighbors(graph = g, v = "D")), color = "red")

#Average Clustering coefficient
transitivity(g)

# Diameter
diameter(g)

#Shortest path
shortest_paths_all <- distances(g)

#Betweeness

node_betweenness <- betweenness(g)

#neighbours
C_neighbors<-neighbors(g, "C", mode = "all")  # Index of neighbours by index

#------------------Plots for detecting network type---------------

#--------------------k vs P (k)----------------
df <- data.frame(table(degree(graph = g))/sum(degree(graph = g)))
df$Var1 <- as.numeric(df$Var1)
colnames(df) <- c("k", "Pk")
ggplot(data.frame(df), aes(x = log(k), y = log(Pk))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)
#--------------------k vs C(k)----------------
df <- data.frame(k = degree(graph = g), c = transitivity(graph = g, type = "local"))
df2 <- NULL
for(i in unique(df$k)) {
  df2 <- rbind(df2, c(i, mean(df[df$k==i,]$c, na.rm = TRUE)))
}
colnames(df2) <- c("k", "Ck")
ggplot(data.frame(df2), aes(x = log(k), y = log(Ck))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)