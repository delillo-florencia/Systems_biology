library(igraph)
library(ggraph)
library(DiscoNet)
#---------------Decompose network-------------
connected_graphs <- decompose(g)  
connected_graphs_weak <- decompose(g,mode="weak")  
connected_graphs_strong <- decompose(g,mode="strong")  

#-----------------------VIRTUAL PULL DOWN----------
#Seeds: proteins of interest
seeds <- c("ALDH1A2", "BMP2", "CXADR", "GATA4", "HAS2", "NF1", "NKX2-5", "PITX2", "PKD2", "RXRA", "TBX1", "TBX2", "ZFPM1", "ZFPM2")

#load(file='/home/projects/22140/inweb_reduced.Rdata') (this loads the database)


network_ex2 <- virtual_pulldown(seed_nodes = seeds,
                 database = db_reduced, 
                 id_type = "hgnc", 
                 zs_confidence_score = 0.156)

interactions <- data.frame(network_ex2$network)
node_attributes <- data.frame(network_ex2$node_attributes)

g <- graph_from_data_frame(interactions, directed = FALSE, vertices = node_attributes)
g1 <- relevance_filtering(g, 0)
g2 <- relevance_filtering(g, 0.5)
g3 <- relevance_filtering(g, 1)


# Community detection
communities <- community_detection(g1, algorithm = "mcode")

# remember generally protein complexes ~ 20 nodes

#----------and after this... FGSEA over relevant communities