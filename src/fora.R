
library(msigdbr)
library(fgsea)

## This requires:
# target --> list of genes we want to check
# background --> universe of genes

#First,  prepare the biological process gene sets for over representation analysis using "fora"

BP_df = msigdbr(species = "S. cerevisiae", category = "C5", subcategory = "BP")
BP_list = split(x = BP_df$ensembl_gene, f = BP_df$gs_name)

MF_df = msigdbr(species = "S. cerevisiae", category = "C5", subcategory = "MF")
MF_list = split(x = MF_df$ensembl_gene, f = MF_df$gs_name)

CC_df = msigdbr(species = "S. cerevisiae", category = "C5", subcategory = "CC")
CC_list = split(x = CC_df$ensembl_gene, f = CC_df$gs_name)


BP_gsea <- fora(BP_list, target, background)
MF_gsea <- fora(MF_list, target, background)
CC_gsea <- fora(CC_list, target, background)

##Take a look at the results. As mentioned, "fora" produces a table with 
#p-values, adjusted p-values, gene set size, and the size of the overlap 
#between the gene set and the target list, but not the enrichment.


### Calculate enrichment

compute_enrichment <- function(ora, n_genes, n_universe) {
  ora$relative_risk <-
    (ora$overlap / ora$size) /
    (n_genes / n_universe)   
  return(ora)
}

BP_gsea <- compute_enrichment(BP_gsea, length(target), length(background))
MF_gsea <- compute_enrichment(MF_gsea, length(target), length(background))
CC_gsea <- compute_enrichment(CC_gsea, length(target), length(background))


### Plot results

# Function to get top N enriched terms based on adjusted p-values
get_top_enriched_terms <- function(gsea_results, n = 10) {
  # Filter for significant terms (adjusted p-value < 0.05)
  significant_terms <-  gsea_results[gsea_results$padj < 0.05, ]
  
  # Sort by adjusted p-value
  sorted_terms <- significant_terms[order(significant_terms$padj), ]
  
  # Return top N terms
  return(head(sorted_terms, n))
}

top_BP_terms <- get_top_enriched_terms(BP_gsea, 10)
top_MF_terms <- get_top_enriched_terms(MF_gsea, 10)
top_CC_terms <- get_top_enriched_terms(CC_gsea, 10)

# Display results
print("Top 10 Biological Process GO Terms:")
print(top_BP_terms)

print("Top 10 Molecular Function GO Terms:")
print(top_MF_terms)

print("Top 10 Cellular Component GO Terms:")
print(top_CC_terms)