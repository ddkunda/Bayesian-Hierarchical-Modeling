# init.R
#
# Example R code to install packages if not already installed
#
my_packages = c("shiny","rstan","ggplot2","parallel","StanHeaders","posterior")
install_if_missing = function(p) {
  if (p %in% rownames(installed.packages()) == FALSE) {
    install.packages(p)
  }
}
invisible(sapply(my_packages, install_if_missing))


Sys.setenv(MAKEFLAGS = "-j2")  # Speed up compilation
Sys.setenv(R_MAKEVARS_USER = "~/.R/Makevars")

# Create .R/Makevars if it doesn’t exist
makevars_path <- file.path(Sys.getenv("HOME"), ".R", "Makevars")
if (!file.exists(makevars_path)) {
  writeLines("CXX17 = g++", makevars_path)
}
