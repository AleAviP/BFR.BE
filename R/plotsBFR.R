#Plots Toy Examples
plot.heat <- function(X,Xlab="",Ylab="",limit=c(-3,3),rotation=FALSE){
  X = as.matrix(X)
  colnames(X)<-NULL
  rownames(X)<-NULL
  p<-nrow(X)
  q<-ncol(X)
  index<-apply(X,2,sum)!=0
  if(sum(index!=0)!=q){
    X<-cbind(X[,c(1:q)[apply(X,2,sum)!=0]],matrix(ncol=(q-sum(apply(X,2,sum)!=0)),nrow=p,0))
  }
  if(rotation==TRUE){
    X<-varimax(X)$loadings%*%diag(q)
  }
  x = melt(X)
  p_X_heat = ggplot(data = x, aes(x=X2, y=-X1, fill=value)) +
    theme_bw() +
    geom_tile(show.legend = F) +
    xlab(Xlab) +
    ylab(Ylab) +
    theme(axis.title=element_text(size=14,face="bold"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    scale_fill_gradient2(limits=limit)
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

plot.heat.cov <- function(M,Sigma,Xlab="",Ylab="",limit=c(-8,8)){
  M = as.matrix(M)
  colnames(M)<-NULL
  rownames(M)<-NULL
  Sigma = as.matrix(Sigma)
  colnames(Sigma)<-NULL
  rownames(Sigma)<-NULL
  X = M%*%t(M)+Sigma
  x = melt(X)
  p_X_heat = ggplot(data = x, aes(x=X2, y=-X1, fill=value)) +
    theme_bw() +
    geom_tile(show.legend = F) +
    xlab(Xlab) +
    ylab(Ylab) +
    theme(axis.title=element_text(size=14,face="bold"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    scale_fill_gradient2(limits=limit)
  return(p_X_heat)
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

