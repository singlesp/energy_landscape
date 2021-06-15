library(R.matlab)
library(ggplot2)
library(ggsignif)
library(ggthemes)

numClusters = 4; #make sure this is reflected in the file name below!!!!!!!
a<-readMat("/Users/sps253/Documents/brain_states-master/results/example/ViolinData_bp22_k4.mat")


LSDfo1<-a$LSDfo
PLfo1<-a$PLfo
clus<-a$clusters
LSDdt1<-a$LSDdt
PLdt1<-a$PLdt
LSDar1<-a$LSDar
PLar1<-a$PLar

s <- 1 #starting index
e <- 15 #ending index
len<-15 #length(LSDfo1[1,s:e])

#option to rename clusters (they will be plotted alphabetically)
#clus=character()
#clus[1]<-"MS-1a"
#clus[2]<-"MS-1b"
#clus[3]<-"MS-2a"
#clus[4]<-"MS-2b"



GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity",  ..., draw_quantiles = NULL,
                               trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, 
        inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#Fraction Occupancy
y<-NULL
x<-NULL
Condition<-NULL

for(i in 1:numClusters){y=c(y, LSDfo1[i,s:e], PLfo1[i,s:e])
x=c(x, rep(clus[i], len*2))
Condition=c(Condition, rep('LSD', len), rep('PL', len))
}

my_data = data.frame(y,x,Condition)

p_value<-NULL
for(i in 1:numClusters){
  p_value<-rbind(p_value,t.test(LSDfo1[i,s:e],PLfo1[i,s:e],paired = TRUE, alternative = "two.sided")$p.value)
} #clusters 2, 4, 5 are sig before correction
p_value
p.adjust(p_value,n=3*numClusters,method="BH") #5 after correction

p <- ggplot(my_data, aes(x, y, fill = Condition)) + geom_split_violin() 
#No correction for multiple comparison in these p-values. * = p<0.05 before correction
label.df <- data.frame(x = c(clus[1],clus[2],clus[3],clus[4]),
                       y = c(0.44,0.44,0.44,0.44),
                       Condition = c('LSD','PL'))
p + geom_text(data = label.df, label = c("*","*","*","*"), size=22) + 
  labs(x = "", y = "Fractional Occupancy") + 
  theme_classic(base_size = 22,base_line_size = 1) + ggtitle("Fractional Occupancy") + ylim(0, 0.45)


#Dwell Time

y<-NULL
x<-NULL
Condition<-NULL

for(i in 1:numClusters){y=c(y, LSDdt1[i,s:e], PLdt1[i,s:e])
x=c(x, rep(clus[i], len*2))
Condition=c(Condition, rep('LSD', len), rep('PL', len))
}

my_data = data.frame(y,x,Condition)


p_value<-NULL
for(i in 1:numClusters){
  p_value<-rbind(p_value,t.test(LSDdt1[i,s:e],PLdt1[i,s:e],paired = TRUE, alternative = "two.sided")$p.value)
} #clusters 1, 3, 4  sig before correction
p_value
p.adjust(p_value,n=3*numClusters,method="BH") #none sig after correction (2 p = 0.064)

p <- ggplot(my_data, aes(x, y, fill = Condition)) + geom_split_violin() 
#These p-values survive BH correction for multiple comparisons * = p<0.05
label.df <- data.frame(x = c(clus[1],clus[2],clus[3],clus[4]),
                       y = c(14,14,14,14),
                       Condition = c('LSD','PL'))
p + geom_text(data = label.df, label = c("*","*","","**"), size=22) + 
  labs(x="", y = "Dwell Time (s)") + 
  theme_classic(base_size = 22,base_line_size = 1) + ggtitle("Dwell Time") + ylim(0, 15)


#Appearance Rate

y<-NULL
x<-NULL
Condition<-NULL

for(i in 1:numClusters){y=c(y, LSDar1[i,s:e], PLar1[i,s:e])
x=c(x, rep(clus[i], len*2))
Condition=c(Condition, rep('LSD', len), rep('PL', len))
}

my_data = data.frame(y,x,Condition)


p_value<-NULL
for(i in 1:numClusters){
  p_value<-rbind(p_value,t.test(LSDar1[i,s:e],PLar1[i,s:e],paired = TRUE, alternative = "two.sided")$p.value)
} #clus 1, 4, 6 sig before correction
p_value
p.adjust(p_value,n=3*numClusters,method="BH") #none sig after

p <- ggplot(my_data, aes(x, y, fill = Condition)) + geom_split_violin() 
#These p-values survive BH correction for multiple comparisons* = p<0.05
label.df <- data.frame(x = c(clus[1],clus[2],clus[3],clus[4]),
                       y = c(3.9,3.9,3.9,3.9),
                       Condition = c('LSD','PL'))
p + geom_text(data = label.df, label = c("","","",""), size=22) + 
  labs(x="", y = "Appearances/min") + 
  theme_classic(base_size = 22,base_line_size = 1) + ggtitle("Appearance Rate") + ylim(0, 4)

