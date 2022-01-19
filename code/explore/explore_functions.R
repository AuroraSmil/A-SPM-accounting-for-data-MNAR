# Explore functions

histogram_compare_missing = function(data, variable, title, ylim_l, ylim_u){
  # returns a histogram comparing present and drop-out participants
  # data: contains the hunt cohort of interest
  # variable: variable we want to explore
  # title: plot tiltle
  # ylim_l: lower limit y axis
  # ylim_u: upper limit y axis
  q_1 <- ggplot(data = data, aes_string(x = variable, group = as.factor("missing")))
  q_1 <- q_1 + geom_histogram(data = data[missing == 1], 
                              aes(y = ..density.., colour = as.factor(missing), fill = as.factor(missing)),alpha = 0.5, bins = 100)
  q_1 <- q_1 + geom_histogram(data = data[missing == 0], 
                              aes(y = ..density.., colour = as.factor(missing), fill = as.factor(missing)),alpha = 0.5, bins = 100)
  q_1 <- q_1 + theme(legend.title = element_blank(), legend.position = "bottom")#, axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  q_1 <- q_1 + scale_fill_discrete(labels = c( "Pressent", "Missing")) 
  q_1 <- q_1 + scale_colour_discrete(guide = 'none')
  q_1 <- q_1 + xlab(unname(unname(TeX(c(title))))) + ylab("Density") + ylim(ylim_l, ylim_u)
  return(q_1)
}

smooth_density_compare_missing = function(data, variable, title, ylim_l, ylim_u){
  # returns a smoothed density plot comparing present and drop-out participants
  # data: contains the hunt cohort of interest
  # variable: variable we want to explore
  # title: plot tiltle
  # ylim_l: lower limit y axis
  # ylim_u: upper limit y axis
  q_1 <- ggplot(data = data, aes_string(x = variable, colour = as.factor("missing")))
  q_1 <- q_1 + geom_density(data = data[missing == 1],aes(y = ..density.., colour = as.factor(missing), linetype = as.factor(missing)))
  q_1 <- q_1 + geom_density(data = data[missing == 0],aes(y = ..density.., colour = as.factor(missing), linetype = as.factor(missing)))
  q_1 <- q_1 + theme(legend.title = element_blank(), legend.position = "bottom")
  q_1 <- q_1 + scale_colour_discrete(labels = c( "Pressent", "Missing")) 
  q_1 <- q_1 + labs(colour  = "Guide name", linetype = "Guide name")
  q_1 <- q_1 + scale_linetype_manual(labels = c( "Pressent", "Missing"), values = c("solid","dashed"))
  q_1 <- q_1 + xlab(unname(unname(TeX(c(title))))) + ylab(element_blank()) + ylim(ylim_l, ylim_u)
  return(q_1)
}



g_legend<-function(a.gplot){
  # returns a legend to be used by arrangegrop
  # copied from https://stackoverflow.com/questions/61232092/ggplot-custom-legend-instead-of-default
  # a.gplot, plot from which one wants to take the legend
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}