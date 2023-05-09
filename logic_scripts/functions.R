### Define functions

# To format P value columns
p_formatter <- function(x) {
  format(x, digits = 2, scientific = TRUE)
}

# To replicate default ggplot colour assignment
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# To remove isolated nodes
remove_isolated <- function(net) {
  require(igraph)
  Isolated <- which(degree(net) == 0)
  net <- delete.vertices(net, Isolated)
  return(net)
}
