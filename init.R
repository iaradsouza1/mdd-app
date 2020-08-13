packages <- c("shiny", "shinythemes", "dplyr", "ggplot2", "DT", "ggrepel")

for (pack in packages) {
  if(!require(pack)) {
    install.packages(pack)
  }
  require(pack)
}