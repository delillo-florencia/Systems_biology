library(igraph)
library(ggraph)

#--------------------- Network to graph ------------------------
hemoglobin <- data.frame(
    from = c("ALPHA_GLOBIN", "ALPHA_GLOBIN", "BETA_GLOBIN"),
    to = c("ALPHA_GLOBIN", "BETA_GLOBIN", "BETA_GLOBIN"))
    # and any other column 

g <- graph_from_data_frame(hemoglobin,directed=FALSE)

plot(g) #--> This plots also self interactions - ggplot no

#-----------Visualizing the graph using ggraph------------------
#Label:  Modify it to show the GeneIDs (you will have to find out how yourself - HINT: try to search google for "ggraph node labels")
#Fill Color:Color the nodes based on the "Catalytic" variable in the
#Node size:  Make the size of the nodes correspond to the length of the amino acid sequence

network <- data.frame(
    from = c("A" , "A", "A", "B"," B" ,"B" ,"C"  ,"D ","D", "E"),
    to = c("A" ,"C" ,"E" ,"B" ,"C", "E", "C", "C", "D","E"),
    score = c(1,0.5,1,1,NA,NA,1,1,1,1))


g <- graph_from_data_frame(network,directed=FALSE)


plot(g) + # add inside plot if needed , 
  geom_edge_link() + 
  geom_node_point(size = 4, color = "blue") +  # Adjust size and color of nodes
  geom_node_text(aes(label = name), vjust = 1, hjust = 1)
                            
#layout = layout_with_fr
#layout_with_kk,               
#layout = layout_in_circle,
#layout = layout_as_tree,
#layout = layout_with_dh,   
#layout = layout_with_lgl,  
#layout = layout_on_grid,   
#layout = layout_randomly,    


g2 <- delete_edges(g, which(E(g)$score < 1))
ggraph(g2, layout = "linear", circular = TRUE) + 
  geom_edge_link() + 
  geom_node_point() +
  geom_node_text(aes(label = names(V(g))), repel = TRUE)


#----------------Deleting nodes from graph----------
updated_g <- delete_vertices(graph = g, v = !names(V(g)) %in% names(neighbors(graph = g, v = "ALPHA_GLOBIN")))

