packages <- c("rmarkdown", "readr", "dplyr", "ggplot2", "scales", "reticulate")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
