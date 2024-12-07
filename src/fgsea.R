
library(msigdbr)
library(fgsea)


#NOTE:You can sort as whatever you want - logFC, expression, etc. 
#Just need a ranked list.
#First,  prepare the biological process gene sets for over representation analysis using "fora"

BP_df = msigdbr(species = "S. cerevisiae", category = "C5", subcategory = "BP")
BP_list = split(x = BP_df$ensembl_gene, f = BP_df$gs_name)

MF_df = msigdbr(species = "S. cerevisiae", category = "C5", subcategory = "MF")
MF_list = split(x = MF_df$ensembl_gene, f = MF_df$gs_name)

CC_df = msigdbr(species = "S. cerevisiae", category = "C5", subcategory = "CC")
CC_list = split(x = CC_df$ensembl_gene, f = CC_df$gs_name)


# Assuming 'log2FC' is the log2(Fold Change) column and 'SysName' contains gene names
log2FC_values <- expr$log2fc  # Replace with your column name for log2FC
gene_names <- expr$SysName     # Replace with your gene names column

# Sort the genes by log2(Fold Change) in descending order (so highly upregulated genes come first)
ranked_log2FC <- log2FC_values[order(-log2FC_values)]  # Negative to sort in descending order
ranked_gene_names <- gene_names[order(-log2FC_values)]

# Create a named vector with log2FC values and the gene names
named_log2FC <- setNames(ranked_log2FC, ranked_gene_names)


gsea_log2FC <- fgsea(BP_list, named_log2FC)

# Filter and sort by significance (padj), and select the top results
top_10_log2FC_upregulated <- gsea_log2FC %>%
  filter(NES > 0) %>%
  arrange(padj) %>%
  head(10)

# Print the top 10 most significantly over-represented biological processes
print(top_10_log2FC_upregulated)