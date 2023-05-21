### Define functions

# To replicate default ggplot colour assignment
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# To remove isolated nodes
remove_isolated <- function(net) {
  require(igraph)
  isolated <- which(degree(net) == 0)
  net <- delete.vertices(net, isolated)
  return(net)
}

# Tooltip over checkboxes
radioTooltip <- function(id, choice, title,
  placement = "bottom", trigger = "hover", options = NULL) {
  require(shinyBS)

  options <- shinyBS:::buildTooltipOrPopoverOptionsList(
    title, placement, trigger, options)
  options <- paste0("{'", paste(
    names(options), options, sep = "': '", collapse = "', '"), "'}")

  bs_tag <- shiny::tags$script(shiny::HTML(paste0("
    $(document).ready(function() {
      setTimeout(function() {
        $('input', $('#", id, "')).each(function(){
          if(this.getAttribute('value') == '", choice, "') {
            opts = $.extend(", options, ", {html: true});
            $(this.parentElement).tooltip('destroy');
            $(this.parentElement).tooltip(opts);
          }
        })
      }, 500)
    });
  ")))
  htmltools::attachDependencies(bs_tag, shinyBS:::shinyBSDep)
}
