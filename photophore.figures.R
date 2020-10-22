library(tidyverse)
library(lemon)


setwd("/Users/biology/Dropbox/Documents/Photophore Size Data/manuscript/figures")

spab <- function(tax=NULL,new.col=NULL ){
  nms <-  tax
  gen <- gsub("(\\w+) (\\w+)","\\1",nms)
  sp <- gsub("(\\w+) (\\w+)","\\2",nms)
  dt <- data.table(name=nms)
  for(i in 1:max(str_length(gen))){
    srch <- paste0("(\\w{",i,"}).+")
    gen2 <- gsub(srch,"\\1.",gen)
    dt[,paste0("ab",i):=gen2]
  }
  dt.m <- melt(dt,c("name"))
  dt.m[,min.n:=gsub("ab","",variable)]
  dt.m[,un:=!value %in% value[duplicated(value)]]
  dt.m <- dt.m[grepl("\\.",value)]
  ab <- dt.m[,.(ab=value[min(which(un=="TRUE"))]),by=.(name)][,sp.ab:=paste0(ab," ",sp)]
  return(ab)
}

d <- list(
  aic.weights.ply, #for model plot
  freq.ply, # for histograms of photophore distributions
  plot.dat, # for flux radial plot
  dist.melt, # for model explainations 
  dens.m, #  density data for plots
  flux.dens.m, # for flex density barplots of each sp.
  dat, #size data
  dat.dist.m, #for size dist along body depth
  size.m, # for fize of each species
  phot.dat, #for density
  phot.dat.m, #for density
  ang.dist.m,
  dat.dist
)


saveRDS(d,"photophore.data.Oct20.RDS")



d.n <- .(
  aic.weights.ply, #for model plot
  freq.ply, # for histograms of photophore distributions
  plot.dat, # for flux radial plot
  dist.melt, # for model explainations 
  dens.m, #  density data for plots
  flux.dens.m, # for flex density barplots of each sp.
  dat, #size data
  dat.dist.m,
  size.m,
  phot.dat,
  phot.dat.m, #for density
  ang.dist.m,#for angle plots
  dat.dist
)


names(d) <- d.n

#stomias fix 
d$dat.dist.m[grep("boa",sp2)]
d$ang.dist.m[grep("boa",sp)]
#stomias fix
d$dat.dist.m %>% 
  mutate_if(is.factor,funs(str_replace_all(., "ferox", ""))) %>% 
  data.table::as.data.table() -> d$dat.dist.m

d$dat.dist.m %>% 
  mutate_if(is.character,funs(str_replace_all(., "ferox", ""))) %>% 
  data.table::as.data.table() -> d$dat.dist.m


d$size.m %>% 
  mutate_if(is.factor,funs(str_replace_all(., "ferox", ""))) %>% 
  data.table::as.data.table() -> d$size.m


#stomias fix angle
d$ang.dist.m %>% 
  mutate_if(is.character,funs(str_replace_all(., "ferox", ""))) %>% 
  data.table::as.data.table() -> d$ang.dist.m
 


#species table
sp.tab <- d$plot.dat[,.(sp,sp2,sp3)][!duplicated(sp2),][,sp.ab:=spab(sp2)$sp.ab]

sp.inc <- sp.tab[!sp.ab%in% c("Ic. ovatus","Ar. aculeatus","Si. elongatus","V. nimbaria","Po. thaeocoryla")]$sp.ab


#### bar plot scale 

x.scale <- scale_x_discrete(limits=c(rev(levels(freq.ply$Var1)),"0"),breaks = c(0,20,40,60,80,100))

#### upper right theme

ur.theme <- theme(strip.background = element_blank(),strip.text.x = element_text(size=16,lineheight=10.0,vjust=0.4,hjust=0.2),strip.text.y = element_text(size=16, angle=-0,lineheight=10.0,vjust=0.1),axis.text.x=element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank())

### main theme

m.theme <- theme(strip.background = element_blank(),strip.text.x = element_text(size=16,lineheight=10.0,vjust=0.4,hjust=0.2),strip.text.y = element_text(size=16, angle=-0,lineheight=10.0,vjust=0.1))

 #### model plot ####

d$aic.weights.ply <- merge(d$aic.weights.ply,sp.tab[,.(sp3,sp.ab)],by.x="sp2",by.y="sp3")

d$aic.weights.ply[sp.ab=="F. boureei",se:=m*0.001]

model.plot <- ggplot(d$aic.weights.ply, aes(x = variable, y=m))+geom_bar(stat = "identity",fill="gray70")+geom_errorbar(aes(ymin = m-se, ymax = m+se))+facet_rep_wrap(sp.ab~., ncol=4)+theme_classic(15)+coord_flip()+xlab("model")+scale_x_discrete()+theme(strip.background = element_blank(),strip.text.x = element_text(size=18,lineheight=10.0,vjust=0.3,hjust=0.2),strip.text.y = element_text(size=24, angle=-0,lineheight=10.0,vjust=0.05))+ylab("mean AICc")


pdf("model.plot.pdf",w=12,h=7)
print(model.plot)

dev.off()

#### distribution w body depth ####
### make panel w density mean

d$phot.dat<- merge(d$phot.dat,sp.tab[,.(sp3,sp.ab)],by.x="sp3",by.y="sp3")


phot.dist.m <- d$phot.dat[,.(n=length(Y2)),by=.(sp.ab,y.bin,spec)][,.(m=mean(n),se=se(n)),by=.(sp.ab,y.bin)]

samp <- d$phot.dat[,.(samp.max=1-(5.25/bd),samp.min=0),by=.(sp.ab,spec)][,.(m.max=mean(samp.max)*100,m.min=0),by=(sp.ab)]
samp[m.max<0,m.max:=0.01]
samp[,c("m.min","m.max"):=list(as.factor(round_any(m.min,10)),as.factor(round_any(m.max,10)))]


error.hist <- ggplot()+geom_segment(data=samp,aes(x=m.max,xend=m.max,y=0,yend=max(phot.dist.m$m)),linetype=4)+geom_bar(data=phot.dist.m, aes(x = y.bin, y=m),stat = "identity",fill=viridis(100,option="D")[65],inherit.aes = F)+geom_errorbar(data=phot.dist.m,aes(x=y.bin,ymin = m-se, ymax = m+se),inherit.aes=F)+facet_rep_wrap(sp.ab~.,ncol=3)+theme_classic(16)+coord_flip()+xlab("depth (%)")+ylab("mean count (+/- s.e)")+scale_x_discrete(breaks = c(0,20,40,60,80,100),limits = rev(levels(phot.dist.m$y.bin)))+theme(strip.background = element_blank(),strip.text.x = element_text(size=16,lineheight=10.0,vjust=0.4,hjust=0.2),strip.text.y = element_text(size=16, angle=-0,lineheight=10.0,vjust=0.1))


pdf("dist.plot.pdf",h=6,w=5)
print(error.hist)
dev.off()
                                                                                      
#### model examples #####
p.dist <-  ggplot(data=d$dist.melt)+geom_line(aes(x=x/10,y=value,color=variable),size=2,col=v.col)+theme_classic(10)+theme(legend.position="none")+ylim(c(0,.4))+facet_wrap(variable~., ncol=1)+xlab("vertical position")+ylab("density")+coord_flip()+scale_x_reverse()+theme(strip.background = element_blank(),strip.text.x = element_text(size=16,lineheight=10.0,vjust=0.4),strip.text.y = element_text(size=16, angle=-0,lineheight=10.0))

pdf("dist.models.pdf",h=6,w=3)
print(p.dist)
dev.off()






########### SIZE ############
#############################

#### size body depth bar plot ######## 

d$dat.dist.m<- merge(d$dat.dist.m,sp.tab[,.(sp3,sp.ab)],by.x="sp2",by.y="sp3")

blank.d <- tail(d$dat.dist.m,5)
blank.d[,sp.ab:=strrep(" ",1:5)]
blank.d[,c("m","se"):=list(NA,NA)]

size.dat2 <- rbind(d$dat.dist.m[sp.ab %in% sp.inc],blank.d)



#d$dat.dist.m[sp.ab %in% sp.inc]
size.error.dist <- ggplot()+geom_segment(data=samp,aes(x=m.max,xend=m.max,y=0,yend=max(d$dat.dist.m$m)),linetype=4)+geom_bar(data=size.dat2, aes(x = factor(y.bin), y=m),stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(data=size.dat2,aes(x = factor(y.bin),ymin = m-se, ymax = m+se))+facet_rep_wrap(sp.ab~.,ncol=3,scales = "free")+theme_classic(16)+coord_flip()+xlab("depth (%)")+ylab("diam. (+/- s.e)")+scale_x_discrete(limits=rev(levels(freq.ply$Var1)),breaks = c(20,40,60,80,100))+theme(strip.background = element_blank(),strip.text.x = element_text(size=16,lineheight=10.0,vjust=0.4,hjust=0.2),strip.text.y = element_text(size=16, angle=-0,lineheight=10.0,vjust=0.1)) 

#### Size barplot for each species###

d$size.m<- merge(d$size.m,sp.tab[,.(sp3,sp.ab)],by.x="sp4",by.y="sp3")
d$size.m$sp.ab <- factor(d$size.m$sp.ab,levels=d$size.m$sp.ab[order(d$size.m$size.m)])


size.plot <- ggplot(d$size.m[sp.ab %in% sp.inc], aes(x = sp.ab, y=size.m))+geom_bar(stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(aes(ymin = size.m-se, ymax = size.m+se))+ylab("mean diam. (um, +/- s.e.)")+theme_classic(15)+theme(axis.text.x = element_text(angle = 75, hjust = .95,vjust=1))+xlab("")


pdf("size.plot.pdf",w=8,h=6)
print(size.plot)
dev.off()



######## size panel figure mean and distribution ##

size.panel2 <- ggplot()+geom_segment(data=samp[ sp.ab %in% unique(samp[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x=m.max,xend=m.max,y=0,yend=max(d$dat.dist.m[sp.ab %in% sp.inc]$m)),linetype=4)+geom_bar(data=d$dat.dist.m[ sp.ab %in% unique(d$dat.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]], aes(x = factor(y.bin), y=m),stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(data=d$dat.dist.m[ sp.ab %in% unique(d$dat.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x = factor(y.bin),ymin = m-se, ymax = m+se))+facet_rep_wrap(sp.ab~.,ncol=1)+theme_classic(16)+coord_flip()+xlab("depth (%)")+ylab("mean diam. (+/- s.e)")+scale_x_discrete(limits=rev(levels(freq.ply$Var1)),breaks = c(20,40,60,80,100))+ur.theme


size.panel3 <- ggplot()+geom_segment(data=samp[sp.ab %in% sp.inc & !sp.ab %in% unique(samp[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x=m.max,xend=m.max,y=0,yend=max(d$dat.dist.m[sp.ab %in% sp.inc]$m)),linetype=4)+geom_bar(data=d$dat.dist.m[sp.ab %in% sp.inc & !sp.ab %in% unique(d$dat.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]], aes(x = factor(y.bin), y=m),stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(data=d$dat.dist.m[sp.ab %in% sp.inc & !sp.ab %in% unique(d$dat.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x = factor(y.bin),ymin = m-se, ymax = m+se))+facet_rep_wrap(sp.ab ~ .,ncol=3) + coord_capped_cart(bottom='both', left='both')+theme_bw()+theme_classic(16)+coord_flip()+xlab("body depth (%)")+ylab("mean diam. (+/- s.e)")+scale_x_discrete(limits=rev(levels(freq.ply$Var1)),breaks = c(20,40,60,80,100))+m.theme



t_row <- plot_grid(size.plot,size.panel2,labels = c('A', 'B'), label_size = 12,rel_widths = c(.65,0.37),rel_heights = c(1.2,1),greedy = T,align = 'v',axis = "t")

size.panel.plot <- plot_grid(t_row, size.panel3, labels = c(''), label_size = 12,rel_heights = c(.33,0.67),ncol=1,greedy = T,align = 'v',axis = "r")

pdf("size.panel.plot.pdf",w=8,h=11)
print(size.panel.plot)
dev.off()



################################################### 
######## density plots                           ##
###################################################
#### density barplot #####

### need to remove sp w no photophores and add sample lines
dens.exp <- expression(paste('mean density (',mm^-2,', +/- s.e.)',sep=''))

phot.dat.m2<- phot.dat2[!is.infinite(dens),.(dens.m=mean(dens),se=se(dens)),by=.(sp,sp3)]

dens.dat2 <- phot.dat[,.(dens=length(X)/(diff(range(X))*diff(range(Y2)))),by=.(sp,sp2,spec,region)]

dens.dat3 <- dat[,.(dens=length(X)/(diff(range(X))*diff(range(Y)))),by=.(sp,sp2,spec,region)]

phot.dat.m3<- dens.dat3[!is.infinite(dens),.(dens.m=mean(dens),se=se(dens)),by=.(sp,sp2)]

d$dens.m <- phot.dat.m3
d$dens.m<- merge(d$dens.m,sp.tab[,.(sp2,sp.ab)],by.x="sp2",by.y="sp2")
d$dens.m$sp.ab <- factor(d$dens.m$sp.ab,levels=d$dens.m$sp.ab[order(d$dens.m$dens.m)])



####### do we want density or number . . . same thing really.

dens.plot <- ggplot(d$dens.m[sp.ab%in%sp.inc], aes(x = sp.ab, y=dens.m))+geom_bar(stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(aes(ymin = dens.m-se, ymax = dens.m+se))+xlab("species")+ylab(dens.exp)+theme_classic(15)+theme(axis.text.x = element_text(angle = 75, hjust = .95,vjust=1))

pdf("dens.plot.pdf",w=8,h=6)
print(dens.plot)
dev.off()


### modify size error plot

count.exp <- 'mean count (+/- s.e.)'

dens.panel2 <- ggplot()+geom_segment(data=samp[sp.ab %in% unique(samp[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x=m.max,xend=m.max,y=0,yend=max(phot.dist.m[sp.ab %in% sp.inc]$m)),linetype=4)+geom_bar(data=phot.dist.m[ sp.ab %in% unique(phot.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]], aes(x = factor(y.bin), y=m),stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(data=phot.dist.m[ sp.ab %in% unique(phot.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x = factor(y.bin),ymin = m-se, ymax = m+se))+facet_rep_wrap(sp.ab~.,ncol=1)+theme_classic(16)+coord_flip()+xlab("depth (%)")+ylab("")+x.scale+ur.theme


dens.panel3 <-  ggplot()+geom_segment(data=samp[sp.ab %in% sp.inc & !sp.ab %in% unique(samp[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x=m.max,xend=m.max,y=0,yend=max(phot.dist.m[sp.ab %in% sp.inc]$m)),linetype=4)+geom_bar(data=phot.dist.m[sp.ab %in% sp.inc & !sp.ab %in% unique(phot.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]], aes(x = factor(y.bin), y=m),stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(data=phot.dist.m[sp.ab %in% sp.inc & !sp.ab %in% unique(phot.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x = factor(y.bin),ymin = m-se, ymax = m+se))+facet_rep_wrap(sp.ab~.,ncol=3)+theme_classic(16)+coord_flip()+xlab("depth (%)")+ylab(count.exp)+x.scale+m.theme

dens.t_row <- plot_grid(dens.plot,dens.panel2,labels = c('A', 'B'), label_size = 12,rel_widths = c(.6,0.35),greedy = T)

dens.panel.plot <- plot_grid(dens.t_row, dens.panel3, labels = c(''), label_size = 12,rel_heights = c(.3,0.67),ncol=1,greedy = T)

pdf("dens.panel.plot.pdf",w=8,h=11)
print(dens.panel.plot)
dev.off()


###########
########### Angle ############
#############################

#### size body depth bar plot ######## 

if(!"sp.ab"%in% colnames(d$ang.dist.m)) d$ang.dist.m<- merge(d$ang.dist.m,sp.tab[,.(sp3,sp.ab)],by.x="sp2",by.y="sp3")

blank.d <- tail(d$ang.dist.m,5)
blank.d[,sp.ab:=strrep(" ",1:5)]
blank.d[,c("m","se"):=list(NA,NA)]

ang.dat2 <- rbind(d$ang.dist.m[sp.ab %in% sp.inc],blank.d)


deg.exp <- expression(paste('mean angle (',degree,', +/- s.e.)',sep=''))

#d$dat.dist.m[sp.ab %in% sp.inc]
ang.error.dist <- ggplot(ang.dat2, aes(x = factor(y.bin), y=m))+geom_bar(stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(aes(ymin = m-se, ymax = m+se))+facet_rep_wrap(sp.ab~.,ncol=3)+theme_classic(16)+coord_flip()+xlab("depth (%)")+ylab(deg.exp)+scale_x_discrete(limits=rev(levels(freq.ply$Var1)),breaks = c(20,40,60,80,100))+theme(strip.background = element_blank(),strip.text.x = element_text(size=16,lineheight=10.0,vjust=0.4,hjust=0.2),strip.text.y = element_text(size=16, angle=-0,lineheight=10.0,vjust=0.1)) 

#### Size barplot for each species###
ang.dat <- dat[,.(m=mean(ang2),se=se(ang2)),by=.(sp2)]
ang.dat <- merge(ang.dat,sp.tab,by="sp2")

add.d <- d$size.m[sizeL.m==0,.(sp2,sp,sp3,sp.ab),][,c("m","se"):=list(0,0)]

ang.dat <- ang.dat[sp.ab %in% sp.inc]

#sort order
ang.dat$sp.ab <- factor(ang.dat$sp.ab,levels=ang.dat$sp.ab[order(ang.dat$m)])

#reorder with null spec first
lv <- c(levels(ang.dat$sp.ab)[!levels(ang.dat$sp.ab)%in%sp.inc],levels(ang.dat$sp.ab)[levels(ang.dat$sp.ab)%in%sp.inc])
ang.dat$sp.ab <- factor(ang.dat$sp.ab,levels=lv)


ang.plot <- ggplot(ang.dat, aes(x = sp.ab, y=m))+geom_bar(stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(aes(ymin = m-se, ymax = m+se))+ylab(deg.exp)+theme_classic(15)+theme(axis.text.x = element_text(angle = 75, hjust = .95,vjust=1))+xlab("")


pdf("ang.plot.pdf",w=8,h=6)
print(ang.plot)
dev.off()



######## ang panel figure mean and distribution ##

ang.panel2 <- ggplot()+geom_segment(data=samp[sp.ab %in% unique(samp[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x=m.max,xend=m.max,y=min(d$ang.dist.m[sp.ab %in% sp.inc]$m),yend=max(d$ang.dist.m[sp.ab %in% sp.inc]$m)),linetype=4)+geom_bar(data=d$ang.dist.m[ sp.ab %in% unique(d$ang.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]], aes(x = factor(y.bin), y=m),stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(data=d$ang.dist.m[ sp.ab %in% unique(phot.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x = factor(y.bin),ymin = m-se, ymax = m+se))+facet_rep_wrap(sp.ab~.,ncol=1)+theme_classic(16)+coord_flip()+xlab("depth (%)")+ylab("")+x.scale+ur.theme

ang.panel3 <-  ggplot()+geom_segment(data=samp[sp.ab %in% sp.inc & !sp.ab %in% unique(samp[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x=m.max,xend=m.max,y=min(d$ang.dist.m[sp.ab %in% sp.inc]$m),yend=max(d$ang.dist.m[sp.ab %in% sp.inc]$m)),linetype=4)+geom_bar(data=d$ang.dist.m[sp.ab %in% sp.inc & !sp.ab %in% unique(d$ang.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]], aes(x = factor(y.bin), y=m),stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(data=d$ang.dist.m[sp.ab %in% sp.inc & !sp.ab %in% unique(d$ang.dist.m[sp.ab %in% sp.inc]$sp.ab)[1:2]],aes(x = factor(y.bin),ymin = m-se, ymax = m+se))+facet_rep_wrap(sp.ab~.,ncol=3)+theme_classic(16)+coord_flip()+xlab("depth (%)")+ylab(deg.exp)+x.scale+m.theme


ang.t_row <- plot_grid(ang.plot,ang.panel2,labels = c('A', 'B'), label_size = 12,rel_widths = c(.6,0.35),greedy = T)

ang.panel.plot <- plot_grid(ang.t_row, ang.panel3, labels = c(''), label_size = 12,rel_heights = c(.3,0.67),ncol=1,greedy = T)

pdf("ang.panel.plot.pdf",w=8,h=11)
print(ang.panel.plot)
dev.off()



###########
########### Flux ############
#############################

#### radial flux plot#####

d$plot.dat<- merge(d$plot.dat,sp.tab[,.(sp3,sp.ab)],by.x="sp3",by.y="sp3")


p2 <- ggplot(data=d$plot.dat[!is.na(bins.ang)&flux.s!=0], aes(x=bins.ang, y=flux.s)) +# Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", aes(fill=log.flux)) +
  theme_minimal(10) +
  #scale_x_discrete(limits = rev(levels(plot.dat$bins.ang)))+
  coord_polar(start = -rad(12.25),direction = -1) +facet_rep_wrap(sp.ab~.,ncol=4)
p.flux <- p2+ scale_fill_viridis(begin = 0.2,end = 0.8,guide = guide_colourbar(direction = "horizontal"))+xlab("")+ylab("normalized log flux")+theme(legend.position =c(0.2,-0.075),plot.margin=unit(c(1, 1, 2.5, 0.5), units="line"))+ scale_x_discrete(breaks=seq(-180,180,45))

pdf("angle.flux.pdf",h=5,w=10)
print(p.flux)
dev.off()


#### mean flux density for each species  ######## 

d$flux.dens.m<- merge(d$flux.dens.m,sp.tab[,.(sp3,sp.ab)],by.x="sp4",by.y="sp3")
d$flux.dens.m$sp.ab <- factor(d$flux.dens.m$sp.ab,levels=d$flux.dens.m$sp.ab[order(d$flux.dens.m$dens.m)])


fdens.exp <- expression(paste('flux density (log(photons) ',mm^-2,')',sep=''))

### are these data log transformed?


f.dens.plot <- ggplot(d$flux.dens.m[sp.ab %in% sp.inc], aes(x = sp.ab, y=dens.m))+geom_bar(stat = "identity",fill=viridis(100,option="D")[65])+geom_errorbar(aes(ymin = dens.m-se, ymax = dens.m+se))+xlab("species")+ylab(fdens.exp)+theme_classic(15)+theme(axis.text.x = element_text(angle = 75, hjust = .95,vjust=1))

pdf("flux.dens.plot.pdf",w=8,h=6)
print(f.dens.plot)
dev.off()

### flux along body

flux.dist <- d$dat.dist[,.(flux.sum=sum(log(flux))),by=.(sp,sp3,spec,y.bin)]

flux.sa <- d$dat.dist[,.(sa=diff(range(X))*diff(range(Y))),by=.(sp,sp3)]
flux.dist2 <- merge(flux.dist,flux.sa)[,flux.dens:=flux.sum/(sa/10)]

flux.dist.m<- flux.dist2[!is.infinite(flux.dens),.(m=mean(flux.dens),se=se(flux.dens)),by=.(sp,sp3,y.bin)]

flux.dist.m<-merge(flux.dist.m,sp.tab[,.(sp,sp.ab)],by.x="sp3",by.y="sp")


#### flux dist is same as size so just angle and flux is fine, so make plot of just flux angle

##### inset for panel plot ##########
in.dat <- d$plot.dat[,.(bins.ang=seq(-180,180,22.5), flux.s=1)]
inset <- ggplot(data=in.dat[!is.na(bins.ang)&flux.s!=0], aes(x=bins.ang, y=flux.s)) +theme_minimal(10) +coord_polar(start = -rad(0),direction = -1)+ylab("")+xlab("")+ scale_x_continuous(breaks=seq(-180,157.5,45))+theme(legend.position = "none")


flux.comb <- ggdraw(p.flux) +draw_plot(inset, x = .4, y=-.045, width = .42,height =.42)


panel.plot <- flux.comb+draw_image(
  "~/Dropbox/Documents/Photophore Size Data/astrophylopic.png"
  ,width = 0.16,height = 0.16,x=0.545,y=0.11)

pdf("flux.panel.plot.pdf",w=10,h=6)
print(panel.plot)
dev.off()



######################################## 
#####  angle radial plot on phy ######
########################################

###ggtree with flux direction
# Make the plot
p <- ggplot(data=d$plot.dat, aes(x=as.factor(bins.ang), y=flux.s)) +# Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", aes(fill=log.flux)) +
  #  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #plot.margin = unit(rep(-.15,4), "cm") ,
    panel.spacing = unit(-2, "lines"),
    strip.text = element_text(size = rel(1.2), vjust = -3.0),
    #strip.text = element_blank()
  ) +
  #scale_x_discrete(limits = rev(levels(plot.dat$bins.ang)))+
  coord_polar(start = 0,direction = -1) 

p<- p+facet_wrap(sp~.,ncol=3)+ scale_fill_viridis(begin = 0.2,end = 0.8)
p

# extract Legend 
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend <- g_legend(p)

# plot trees with continuous and discrete tip data
psub <- ggtree(tree3,ladderize=T)+scale_y_reverse()+geom_tiplab(offset=2.5, size=4)


p <- ggplot(data=d$plot.dat, aes(x=as.factor(bins.ang), y=flux.s)) +# Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", aes(fill=log.flux)) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #plot.margin = unit(rep(-5,10), "cm") ,
    panel.spacing = unit(-1.5, "lines"),
    #strip.text = element_text(size = rel(0.85), vjust = 0),
    strip.text = element_blank(),
    legend.position="none"
  ) +
  #scale_x_discrete(limits = rev(levels(plot.dat$bins.ang)))+
  coord_polar(start = -rad(12.25),direction = -1) 

p <- p+facet_wrap(sp~.,ncol=1,as.table = T)+ scale_fill_viridis(begin = 0.2,end = 0.8)
p

wh <- 1-1/(nrow(tip.dat))
vp <- viewport(width =1, height = 1.125, y =0.5, x = 0.93)

pdf("angle.flux.tree.pdf")
psub
print(p, vp = vp)
vp2 <- viewport(width =2, height = 2.0, y =0.5, x = 0.3)
print(legend,vp=vp2)
pushViewport(vp2)
grid.draw(legend) 
dev.off()


#### combine data for table
dens.m
size.m
flux.dens.m
ang.dat

setkeyv(d$dens.m,"sp2")
setkeyv(size.m,"sp2")
setkeyv(flux.dens.m,"sp2")
setkeyv(ang.dat,"sp2")

#ashley needs Avg SL, Avg. Number of Photophores, Avg. Size of Photophores, Avg. Photophore density, and Avg. Photophore orientation angle 


sl.m <- sl.dat[,.(sl.m=mean(sl),se.sl=se(sl)),by=species]
n.dat <- phot.dat[,.(n=length(unique(spec))),by=sp2]

setkeyv(sl.m,"species")
setkeyv(n.dat,"sp2")

dat.table <- sl.m[n.dat,]
dat.table <- dat.table[size.m[d$dens.m[flux.dens.m[ang.dat,],],],]

write_csv(dat.table,"data.summary.table.csv")
