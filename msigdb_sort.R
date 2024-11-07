# If condition `T` is met, execute the following code block
if(T){
  
  # Retrieve and display available MSigDB collections
  msigdbr::msigdbr_collections()
  
  # Load gene sets for Homo sapiens from MSigDB database
  all_gene_sets_hs = msigdbr::msigdbr(species = "Homo sapiens")
  dim(all_gene_sets_hs)  # Display dimensions of the loaded dataset
  saveRDS(all_gene_sets_hs, file="./datasets/all_gene_sets_hs_msigdb.rds")  # Save dataset to an RDS file
  
  # Load gene sets for Mus musculus (mouse) from MSigDB database
  all_gene_sets_hs = msigdbr::msigdbr(species = "Mus musculus")
  dim(all_gene_sets_hs)  # Display dimensions of the loaded dataset
  
  # Display counts of gene set subcategories and categories
  table(all_gene_sets_hs$gs_subcat)
  table(all_gene_sets_hs$gs_cat)
  
  # Filtering: Identify pathways related to "APOPTOSIS"
  # (Adjust `pattern` parameter to search for different keywords if desired)
  h2 <- all_gene_sets_hs
  
  # Check the number of unique pathways
  length(unique(h2$gs_name))
  
  # Check the number of unique genes within all pathways (both human and mouse symbols)
  length(unique(h2$human_gene_symbol))
  length(unique(h2$gene_symbol))
  
  # Display gene set subcategory counts
  table(h2$gs_subcat)
  
  # Further filtering for specific categories and subcategories (C3, TFT:GTRD, TFT:TFT_Legacy)
  library(dplyr)
  h2[h2$gs_cat=='C3', ] %>% .$gs_subcat %>% table()
  
  h2[h2$gs_cat=='C3' & 
       h2$gs_subcat=='TFT:GTRD' | h2$gs_subcat=='TFT:TFT_Legacy', ] %>%
    .$gs_subcat %>% table()
  
  # Check working directory
  getwd()
  
  # Optional: Export data to an Excel file or save as an RDS file
  # openxlsx::write.xlsx(h2, "allPathway.xlsx")
  # save(h2, file = "allPathway.rds")
  
  # Convert the data into a nested list structure for organized output
  nested_list <- h2 %>%
    group_by(gs_cat, gs_name) %>%
    summarise(gene_symbol = list(gene_symbol), .groups = 'drop') %>%
    group_by(gs_cat) %>%
    summarise(gs_name = list(set_names(gene_symbol, gs_name)), .groups = 'drop') %>%
    deframe()
  
  # Save the nested list to a .rdata file
  save(nested_list, file = "./datasets/msigdbr-allPathway.rdata")
  
  # Display sample results from the nested list
  nested_list$C1[[1]]
  head(names(nested_list$C1))
  
}