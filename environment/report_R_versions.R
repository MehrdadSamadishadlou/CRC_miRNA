packages <- c(
  "GEOquery",
  "limma",
  "dplyr",
  "glmnet",
  "clusterProfiler",
  "org.Hs.eg.db"
)

print(R.version.string)

for (package_name in packages) {
  cat(package_name, ": ")
  print(packageVersion(package_name))
}

sessionInfo()
