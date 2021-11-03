#create folders
library(ProjectTemplate)
load.project()


all_folders <- list(file.path(".", "output"),
                    file.path(".", "cache"),
                    file.path(".", "output", "01-synthetic_example"),
                    file.path(".", "output", "02-genotoxic"),
                    file.path(".", "output", "03-depmap_preprocess"),
                    file.path(".", "output", "04-depmap_postprocess"),
                    file.path(".", "output", "05-function_annot"),
                    file.path(".", "output", "06-sparse_codes"),
                    file.path(".", "output", "07-complexes"),
                    file.path(".", "output", "08-landscape"),
                    file.path(".", "output", "09-projection"),
                    file.path(".", "output", "supplement"),
                    file.path(".", "output", "portal_assets"))
walk(all_folders, create_output_folder)