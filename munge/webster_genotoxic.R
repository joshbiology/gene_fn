#munge reference webster run for genotoxic dataset
#requires ./munge/helpers.R to be loaded

load("./cache/genotoxic_input.RData")
mat_path <- file.path(".", "data", "interim", "webster_genotoxic_freeze_2021-6-08-durocher-K=10-T=2.mat")
webster_genotoxic <- graphdl_to_factorized(import_graphdl(R.matlab::readMat(mat_path)),
                                        colnames(genotoxic_input),
                                        rownames(genotoxic_input))


#Manual ordering and annotation
function_names <- c("Polyamine",
                    "Homologous recombination",
                    "ATRi vulnerability",
                    "Nucleotide excision repair",
                    "Nedd. resistance",
                    "Common TSG",
                    "End joining",
                    "Fanconi anemia",
                    "Fork quality control",
                    "Common Essential")

dna_damage_functions <- c(4, 8, 2, 7, 9)

resistance_functions <- c(3, 5, 1)

common_functions <- c(10, 6)

function_order <- c(dna_damage_functions, resistance_functions, common_functions)
