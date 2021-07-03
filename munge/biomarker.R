process_biomarker_output <- function(df, num_cols =44) {
  df %>% select(1:num_cols) %>% 
    pivot_longer(names_to = c("Feature_Number", ".value"),
                 names_pattern = "feature(.*)_(.*)",
                 cols = 5:num_cols) %>% 
    mutate(Feature_Number = as.numeric(Feature_Number)+ 1,
           Type = word(name,start = -1, sep = "_"),
           Feature_Name = word(name,end = -2, sep = "_")) %>% 
    select(-name)
}

#Adjust plotting code

gene_models <- read_tsv("./data/interim/biomarker/latent-representation_v4-achilles-ensemble.tsv") %>% process_biomarker_output()
fn_models <- read_tsv("./data/interim/biomarker/latent-representation_v4-latent-var-ensemble.tsv") %>% process_biomarker_output()
