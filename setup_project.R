## First github project
## Creating the project directory structure

folders <- c(
  "data/raw",
  "data/processed",
  "data/synthetic",
  "R",
  "analysis",
  "figures",
  "documentation"
)

for (folder in folders) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}
