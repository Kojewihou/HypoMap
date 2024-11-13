INPUT <- "/path/to/solution/matrix/example.csv"
OUTPUT <- "/path/to/output/directory/output.csv"


# Get the file path from command-line arguments
leiden_matrix <- as.matrix(read.csv(INPUT, row.names=1))

# Source MRtree.R function from image
source('clustering_utils/_mrtree.R')

# Run MRtree
out <- mrtree(
  leiden_matrix,
  consensus = FALSE,
  sample.weighted = TRUE,
  augment.path = FALSE,
  verbose = TRUE,
  n.cores = 40
)

# Return reconciled MRtree
write.csv(out$labelmat.recon, OUTPUT)
