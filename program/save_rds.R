library(getopt)
library(GENOVA)
library(strawr)

spec = matrix(c(
                'rds', 'd', 1, "character",
                'resolution', 'r', 1, "integer",
                'hic_file', 'h', 1, "character",
                'sample_name', 'n', 1, "character",
                ), byrow=TRUE, ncol=4)
opt = getopt(spec)




matrix_test <- load_contacts(
    signal_path = opt$hic_file,
    resolution = opt$resolution,
    balancing = TRUE,
    sample_name = opt$sample_name
)

saveRDS(matrix_JAK2_DddA11_rep1_10k, file = opt$rds)