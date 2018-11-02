#Plots Toy Examples
plot.heat = function(X){
  X = as.matrix(X)
  x = melt(X)
  p_X_heat = ggplot(data = x, aes(x=X2, y=-X1, fill=value)) +
    theme_bw() +
    geom_tile(show.legend = F) +
    xlab("q") +
    ylab("p")
  return(p_X_heat)
}

plot.scat = function(X,X_recons){
  X_scat = data.frame(true=c(X),reconstruction=c(X_recons))
  p_X = ggplot(X_scat,aes(x=true,y=reconstruction))+geom_point(color="darkgrey",alpha=0.075) + #geom_point(color="cadetblue3",alpha=0.075) +
    theme_bw() +
    geom_abline(slope = 1,linetype=2)+
    theme(axis.title=element_text(size=15,face="bold"))
  return(p_X)
}

plot.hist = function(X,X_recons,normal=FALSE){
  X_dif = data.frame(difference=c(X-X_recons))
  p1 = ggplot(X_dif, aes(x=difference))+
    geom_histogram(aes(y = ..density..,fill=..count..),binwidth = 0.25)+
    ggtitle("Histogram") +
    theme_bw()+
    theme(axis.title=element_text(size=15,face="bold"))
  if (normal==TRUE){
    p1 = p1+stat_function(fun=function(x) dnorm(x), linetype=1, colour="dark blue", size =1.25)
  }
  return(p1)
}

plot.trace = function(X,X_recons){
  XPlot = X_recons
  for (i in 1:ncol(X_recons)){
    XPlot[,i] = X_recons[,i]/X[i]
  }
  meltXplot <- melt(XPlot)
  p_X = ggplot(meltXplot,aes(x=Var1,y=value,colour=as.factor(Var2),group=Var2)) + geom_line(size =1.5)+
    theme_bw() +
    xlab("Iteration") +  theme(legend.position="none") +ylab("Parameter estimation / True Parameter")
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
