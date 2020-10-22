library(scales)        # working with ggplot2 for label formatting
library(gridExtra)     # working with ggplots for arranging plots
library(ggthemes)      # clean theme for ggplot2
library(viridis)       # color palette
library(ggplot2)
library(MASS)
library(grid)
library(qpcR)#for dAIC
require(RSAGA) #for z from matrix
require(reshape2)
require(AICcmodavg)
require(MuMIn)
library(VGAMdata)
library(data.table)
library(bbmle)
require(phangorn)
require(ape)
require(geiger)
require(phytools)
library(plyr)
library(geiger)
library(phylobase)
library(picante)
library(rfishbase)
library(rgbif)
library(AICcmodavg)
library(nlme)
library(ggimage)
library(lmerTest)
library(sjstats)
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")
library(ggtree)
library(cowplot)
library(easynls)
library(lme4)
library(broom)


library("fitdistrplus") #for distribution fitting


#suprisingly, R doesn't have a built in std erro function
se <- function(x)(sd(x)/sqrt(length(x)))

#AICc operations for fitdist objects
my.AICc <- function(x){-2*x$loglik+2*length(coef(x))*(x$n/(x$n-length(coef(x))-1))}#AIC for fitdist objects

dAIC <- function(x){ sapply(x,function(y) y-min(x,na.rm = T))}

#pass dAICc to AICw
AICw <- function(x){
  sum.k <-sum(exp(-0.5*x),na.rm=T)
  sapply(x,function(y) exp(y*-0.5)/sum.k)#here
}


rad <- function(x){
  if(abs(x)>180){x <- x-180}
  return(x*(pi/180))}

#read in latest tree-plotting data to synchronize species w size and count data
plot.dat2 <- readRDS("~/Dropbox/Documents/Photophore Size Data/Data/plot.dat.RDS")

setwd("~/Dropbox/Documents/Photophore Size Data/Data/CPK")
dat.files <- list.files(pattern="_")
dat.files <- dat.files[grepl(".txt|.xls",dat.files)]
### SL data
sl.dat <- data.table(read.csv("specimen.data.csv")) #prolly don't need
sl.dat[,sp2:=gsub("(\\w)\\w+ (\\w+)","\\1\\2",species)]

#body depth
bd <- readRDS("~/Dropbox/Documents/Photophore Size Data/Data/body.depth.RDS")

bd[,sp2:=gsub("(\\w)\\w+ (\\w+)","\\1\\2",species)]

######################################## 
######## COMPUTE COUNT AND FREQ ########
########################################

#store data in a list
dat <- list()
freq <- list()
count<- list()
for(i in dat.files){
  print(i)
  met.dat <- unlist(strsplit(i,"_")) 
  spec <- paste0(met.dat[1],met.dat[2]) #retrieve fish number
  spec2 <- as.numeric(gsub("sp\\.","",met.dat[2])) #retrieve fish number
  sp <- met.dat[1]
  region <- as.numeric(gsub("R(\\d)\\.\\w+","\\1",met.dat[3]))#retieve section
  if(grepl("xls",i)){
    print(met.dat[3])
    region <- as.numeric(gsub("R(\\d)","\\1",met.dat[3]))
  } #retieve section
  #sl <- as.numeric(gsub("mmSL","",met.dat[4]))
  sl<- sl.dat[sp2==sp & spec==spec2]$sl
  sp2 <- sl.dat[sp2==sp & spec==spec2]$species
  if(sp2%in%bd$species){
  dat.i <- read.csv(i,sep="\t")[,c("X","Y")] #read in data
  
  if(max(dat.i$X)>100){
    dat.i$Y <- dat.i$Y/1040*5.25
    dat.i$X <- dat.i$X/1388*7.01
  }
  spc <- spec
 
  if(nrow(bd[sp2==sp&spec==spec2])!=0){
    
 dat.i$Y2 <-   bd[sp2==sp&spec==spec2]$bd-max(dat.i$Y)+dat.i$Y
 dat.i$bd <- bd[sp2==sp&spec==spec2]$bd
 
 #if field is longer than BD, use Y data
 if(diff(range((dat.i$Y)))>bd[sp2==sp&spec==spec2]$bd)  dat.i$Y2 <- dat.i$Y
 }else{
   dat.i$Y2 <- NA
   dat.i$bd <- NA
   }
  
  
  dat.i <- dat.i[dat.i$X!=0,] #some dat files w every other row zeros
  dat[[i]] <- data.frame(dat.i,sp=sp,sp2=sp2,spec=spec,region=region,sp.region=paste0(sp,region),sp.region.spec=paste0(sp,region,spec))#combine data and take only first two cols of raw data
  freq.i <- table(cut(dat.i$Y2,seq(0,3.5,0.35),labels=seq(10,100,10)))#cut Y values into 10 bins by 0.35 units, 
  freq[[i]] <- data.frame(freq.i,sp=sp,spec=spec,region=region,sp.region=paste0(sp,region))
  count[[i]] <- data.frame(count=nrow(dat.i),sp=sp,spec=spec,region=region,sl) #count data
  }
}

#bind data from list to one data frame
phot.dat <- data.table(do.call(rbind,dat))#x,y position
freq.dat <- data.table(do.call(rbind,freq)) #frequency
count.dat <- data.table(do.call(rbind,count)) #count

#fix metallicus/vallaintii
phot.dat[sp2=="Bathophilus vaillanti",c("sp","sp2"):=list("Bmetallicus","Bathophilus metallicus")]

phot.dat[,sp3:=gsub("(\\w).+ (\\w+)","\\1\\. \\2",sp2)]



######### no midlateral section
no2 <- c("Pguernei","Ifasciola","Ebarbatum")
#only mid lateral
phot3 <- phot.dat[sp%in%no2 & region==3] 
phot3[,region:=2]
phot.dat <- phot.dat[!sp%in%no2] 
phot.dat <- rbind(phot.dat[region==2],phot3)

phot.dat[,scale.y:=(Y2/bd),by=.(sp,spec,region)]

phot.dat[,y.bin:=cut(scale.y,breaks=seq(0,1,0.1),labels = seq(10,100,10),right = F),by=.(sp,spec)]

phot.dat2 <- phot.dat[,.(dens=length(X)),by=.(sp,sp3,spec,y.bin)]

sa <- phot.dat[,.(sa=diff(range(X))*diff(range(Y))),by=.(sp,sp3)]
phot.dat2 <- merge(phot.dat2,sa)[,dens:=dens/(sa/10)]

phot.dat.m<- phot.dat2[!is.infinite(dens),.(dens.m=mean(dens),se=se(dens)),by=.(sp,sp3,y.bin)]

freq3 <- freq.dat[sp%in%no2 & region==3] 
freq3[,region:=2]
freq.dat <- freq.dat[!sp%in%no2] 
freq.dat <- rbind(freq.dat[region==2],freq3)

dens.dat <- phot.dat[,.(dens=length(X)/(diff(range(X))*diff(range(Y2)))),by=.(sp,sp2,spec,region)]

dens.m <- dens.dat[,.(dens.m=mean(dens),se=se(dens)),by=.(sp,sp2)]

sd.m <- mean(dens.dat[,.(sd=sd(dens)),by=.(sp,sp2)]$sd,na.rm = T)

nulls <- data.table(sp=c("Selongatum","Iovatus","Pthaeocoryla","Aaculeatus"),sp2=c("Sigmops elongatus","Ichthyococcus ovatus","Polymetme thaeocoryla","Argyropelecus aculeatus"),dens.m=0,se=NA)

dens.m <- rbind(dens.m,nulls)
dens.m[,sp3:=gsub("\\s","_",sp2)]

######################################## 
##### Data for freq distributions ######
########################################

phot.dat <- phot.dat[sp2%in%plot.dat$sp2]

#histogram of count according to bins
histo <- ggplot(phot.dat, aes(x = scale.y))+geom_histogram(bins=10,aes(fill = ..count..))+scale_fill_viridis(name="count", end=0.9)+facet_wrap(sp~.,ncol=2 )+xlim(c(1,0))+theme_classic()+coord_flip()


#cast data with mean and se
freq.ply <- ddply(freq.dat,.(Var1,sp,region,sp.region),summarise,mean.freq=mean(Freq),se.freq=se(Freq)) #mean freq for specie and region by 10 bins, rounded to nearest integer

#cast data with mean and se
count.ply <- ddply(count.dat,.(sp,region,count),summarise,mean.freq=mean(count),se.freq=se(count)) #mean freq for specie and region by 10 bins, rounded to nearest integer

#let's use AIC to find which distributions fit best
aic.weights <- list()
dAICs <- list()
ks.stats <- list()
fit.2P <- "gamma" #choice of 2 parameter models "gamma","lognormal","Weibull" [gamma fits really well]

disc=F

pdf("dist.plots.pdf")
for(i in unique(phot.dat$sp.region.spec)){
  fit.dat <- subset(phot.dat,sp.region.spec==i)$Y2
  met.dat <- subset(phot.dat,sp.region.spec==i)[,c("sp","region","spec")][2,]#retrieve metadat
  fit.dat.rev <- abs(fit.dat-max(fit.dat)) #reversed to fit distribution
  fit.dat.unif <- abs(fit.dat-mean(fit.dat))
  
  if(length(fit.dat)>=10){
    fitU <- fitdist(fit.dat,"unif",discrete = disc)
    fitU$loglik <-  -length(fit.dat)*log(max(fit.dat)-min(fit.dat)) #from: https://www.rdocumentation.org/packages/ExtDist/versions/0.6-3/topics/Uniform
    

    
  
    scale(fit.dat,center = F)
    i.titt <- gsub("(\\w)(\\w+)sp\\.(\\d)","\\1\\. \\2 specimen \\3",met.dat$spec)
    hist(fit.dat, main=i.titt)
    #choose one of the following 2 parameter models (shape included)
    #fit2P <- fitdist(fit.dat.rev, fit.2P)#a two parameter dist model
    fitEg <- fitdist(fit.dat.rev, "exp",discrete = disc) #growth, dec exp fit to data reversed (increases ventrally)
    fitEd <- fitdist(fit.dat, "exp",discrete = disc) #decay w/ data as is (decr ventrally)
    fitN <- fitdist(fit.dat, "norm",discrete = disc)
    fitG <- fitdist(fit.dat, "gamma",discrete = disc) #don't need to rev data for gamma
    fitLN <- fitdist(fit.dat.rev+0.00001, "lnorm") #don't need to rev data for gamma

    
   fits <- list(fitU,fitN,fitEg,fitEd,fitLN)
    aic <- unlist(lapply(fits,my.AICc))
    #unlist(lapply(fits,function(x) x$aic))
    fit.names <- c("unif","norm","exp.g","exp.d","rev.ln")
    #ks <- gofstat(fits,fitnames=fit.names) #retrieve KS stats
  ##aic <- ks$aic
    names(aic) <- fit.names
    lapply(fits,gofstat)
    #aic <- gofstat(fits,fitnames=fit.names) #retieve AIC stats
    aic.weights[[i]] <- data.frame(t(AICw(dAIC(aic))),met.dat)
    dAICs[[i]] <- data.frame(t(dAIC(aic)),met.dat)
    #ks <- gofstat(fits,fitnames=fit.names) #retrieve KS stats
    #ks.stats[[i]] <- data.frame(t(ks),met.dat)
    
  }
  else{aic.weights[[i]] <- NA;dAICs[[i]] <- NA;
  }
}
dev.off()

aic.weights <- do.call(rbind,aic.weights)
aic.weights <- aic.weights[complete.cases(aic.weights),]
dAICs <- do.call(rbind,dAICs)

phot.dat <- data.table(phot.dat)

#ks.test(fit.dat,fit.g)
#ks.test(fit.dat,"gamma")

#denscomp(fits, legendtext=fit.names) #could use to compare density distributions, might be nice to show

#melt the AIC data
aic.weights.melt <- data.table(melt(aic.weights,id.vars=c("sp","spec","region")))

phot.dat[,.(.N),by=.(sp,spec,region)]

aic.weights.ply <- ddply(aic.weights.melt,.(sp,variable,region),summarise,m=mean(value,na.rm = T),se=se(value))
head(aic.weights.melt)#cast data for mean +/se for AICw

aic.ply2 <- subset(aic.weights.ply,sp %in% c("Aniger","Cpliopterus","Csloani","Tdentex"))
levels(aic.ply2$region) <- levels(aic.weights.ply$region) <- c("anterior","mid lateral","posterior")

#levels(aic.ply2$sp)[c(1,3,4,17)] <- c("A. niger", "C. pliopterus", "C. sloani", "T. dentex")


aic.melt <- melt(aic.weights.ply[,-5],id.vars=.(sp,region,variable))


aic.cast <- acast(aic.melt,region+sp~variable,function(x) round(mean(x),2))


aic.table <- xtable(aic.cast,caption="AIC weight comparisions between models of photophore vertical distribution in 17 species of dragonfishes")

print.xtable(aic.table,file="aic.table.tex")

aic.weights.ply <- data.table(aic.weights.ply)

aic.weights.ply$sp2 <- gsub("(\\w)(\\w+)","\\1\\. \\2",aic.weights.ply$sp)

# ##### data for distribution curve examples
v.col <- viridis(100,option="D")[65]

dist.dat <- data.frame(x=seq(0.1,10,0.01),unif=dunif(seq(0.1,10,0.01),min = 0.1,max=10),rev.ln=dlnorm(seq(10,0.1,-0.01),meanlog = 0.6),exp.d=dexp(seq(0.1,10,0.01),1.1),exp.g=dexp(seq(10,0.1,-0.01),1.1),norm=dnorm(seq(0.1,10,0.01),5,2))

dist.melt <- melt(dist.dat,id.vars="x")

dist.col <- data.frame(dist=levels(dist.melt$variable),
                       col1=c(v.col,"black","black","black","black"),
                       col2=c("black",v.col,"black","black","black"),
                       col3=c("black","black",v.col,"black","black"),
                       col4=c("black","black","black",v.col,"black"),
                       col5=c("black","black","black","black",v.col))

dist.alph <- data.frame(dist=levels(dist.melt$variable),a1=c(1,0.5,0.5,0.5,0.5),a2=c(0.5,1,0.5,0.5,0.5),a3=c(0.5,0.5,1,0.5,0.5),a4=c(0.5,0.5,0.5,1,0.5),a4=c(0.5,0.5,0.5,0.5,1))

dist.col.l <- list()
dist.alph.l <- list()
for(i in levels(dist.melt$variable)){
  dist.col.l[[i]] <- as.matrix(dist.col[dist.col$dist==i,2:6])
  dist.alph.l[[i]] <- as.matrix(dist.alph[dist.alph$dist==i,2:6])
}

  
  
#freq.ply2 <- subset(freq.ply,sp %in% c("Aniger","Cpliopterus","Csloani","Tdentex"))
levels(freq.ply$region) <- c("anterior","mid lateral","posterior")

levels(freq.ply$sp) <- gsub("(\\w)(\\w+)","\\1\\. \\2",levels(freq.ply$sp),ignore.case = F)
levels(count.ply$sp) <- gsub("(\\w)(\\w+)","\\1\\. \\2",levels(count.ply$sp),ignore.case = F)
levels(aic.weights.ply$sp) <- gsub("(\\w)(\\w+)","\\1\\. \\2",levels(aic.weights.ply$sp),ignore.case = F)

#count.ply2 <- subset(count.ply,sp %in% c("Aniger","Cpliopterus","Csloani","Tdentex"))
levels(count.ply$region) <- c("anterior","mid lateral","posterior")

#histogram date with plotted error
plot.dat2$sp4 <- gsub("(\\w)\\w+ (\\w+)","\\1\\. \\2",plot.dat2$sp2)
freq.ply <- data.table(freq.ply)



######################################## 
##### Density data     ################
########################################
dens.m$sp4 <- gsub("(\\w)(\\w+)","\\1\\. \\2",dens.m$sp) 

dens.m$sp4 <- factor(dens.m$sp4,levels=dens.m$sp4[order(dens.m$dens.m)])


######################################## 
####### photophore size data   #####
######################################## 

setwd("~/Dropbox/Documents/Photophore Size Data/Data")

sl.dat <- data.table(read.csv("~/Dropbox/Documents/Photophore Size Data/Data/CPK/specimen.data.csv")) #prolly don't need
sl.dat[,sp2:=gsub("(\\w)\\w+ (\\w+)","\\1\\2",species)]

sl.dat[,sp2:=gsub(" ","",sp2)]

#number of lines at end that contain sample meta data
last.dat <- c(rep(3,2),3,2,2,3,2,2,2,2,2,2)

size.met <- read.csv("size.met.dat.csv")
## function to run through table to retrieve XY and area

#body depth
bd <- readRDS("~/Dropbox/Documents/Photophore Size Data/Data/body.depth.RDS")

bd[,sp2:=gsub("(\\w)\\w+ (\\w+)","\\1\\2",species)]

comb.dat <- function(x,refs){
  refs$refs <- as.numeric(as.character(refs$refs))
  refs$angle <- as.numeric(as.character(refs$angle))
  if(refs$refs[1]>0){x <- x[-c(refs$refs),]}#remove ref angles and dimension measurements
  area <- x[nrow(x),]$Length*x[nrow(x)-1,]$Length
  x <- x[-c(nrow(x),nrow(x)-1),]
  if(length(refs$refs)==1){x$Angle <- x$Angle-refs$angle} #if only one ref angle, pass through df
  if(length(refs$refs)>1){
    for(n in 1:length(refs$refs)){
      if(n==length(refs$refs)){x$Angle[refs$refs[n]:nrow(x)]<- x$Angle[refs$refs[n]:nrow(x)]-refs$angle[n]}else{x$Angle[refs$refs[n]:(refs$refs[n+1]-1)] <- x$Angle[refs$refs[n]:(refs$refs[n+1]-1)]-refs$angle[n]}
    }
  }	
  #print(x)
  l.ang <- x[seq(2,nrow(x),2),]
  xy <-x[seq(1,nrow(x),2),] 	
  return(data.frame(xy[,c("X","Y")],l.ang[,c("Angle","Length")]))
}

dat <- list()
freq <- list()
size<- list()
dat.files <- list.files(pattern=".xls")
dat.files <- dat.files[grep("SIZE", dat.files)]

#dat.files <- dat.files[grep("Lglad",dat.files)]
for(i in dat.files){
  meta.dat <- unlist(strsplit(i,"_")) 
  sp <- meta.dat[1] #retrieve fish number

  spec <- spec2 <- meta.dat[2] #retrieve specimen
  

  region <- as.numeric(unlist(strsplit(meta.dat[5],"R")))[2] #retieve section
  sl <- as.numeric(unlist(strsplit(meta.dat[4],"mmSL")))
  dat.i <- read.csv(i,sep="\t")#read in data
  
  size.met.dat <- subset(size.met,dat.files==i)
  sp2 <- size.met.dat$sp
  sp3 <- gsub("(\\w+) .+","\\1",size.met.dat$sp)
  ref.ang <-if(size.met.dat$vert==T){
    data.frame(refs=0,angle=0)
    
  }else{refs <- as.numeric(unlist(strsplit(as.character(size.met.dat$ref.angle),split=",")))
  data.frame(refs,angle=dat.i[refs,"Angle"])}
  
  if(size.met.dat$vert==T) dat.i$Angle <- dat.i$Angle+90
  print(spec)
  print(ref.ang)
  
  dat.i <- data.table(dat.i)
  dat[[i]] <- data.frame(comb.dat(dat.i,ref.ang),sp=sp,spec=spec,region=region,sp2=sp2,sp3=sp3,sl=sl)
  #freq.i <- table(cut(dat.i$Y,seq(0,3.5,0.35),labels=seq(10,100,10)))#cut Y values into 10 bins by 0.35 units, 
  #freq[[i]] <- data.frame(freq.i,sp=sp,spec=spec,region=region,sp.region=paste0(sp,region))
  #count[[i]] <- data.frame(count=nrow(dat.i),sp=sp,spec=spec,region=region,sl) #count data
}

dat <- data.table(do.call(rbind,dat))
dat <- merge(dat,bd[,.(species,spec,bd)],by.x=c("sp2","spec"),by.y=c("species","spec"))

dat.melt.w <- melt(dat,id.vars=c("sp","spec","sp2","region","Y","X","Length"))

b <- seq(0,max(dat$Length),length.out=20)
dat$bins <- with(dat, cut(Length, breaks=b))
levels(dat$bins) <- as.factor(round(b,2))


dat.ply <- ddply(dat,.(sp,sp2,bins),summarise,mean=mean(length(Length)),se=se(length(na.exclude(Length))))


######################################## 
####### flux data   #####
######################################## 

flux.dat <- read.csv("mes.case.flux.txt",sep="\t")[c("X","Y")]
colnames(flux.dat) <- c("Length","flux")
flux.dat$flux <- 10^flux.dat$flux

flux.lm <- glm(log(flux)~Length,flux.dat,family=gaussian)
summary(flux.lm)


dat$log.L <- log(dat$Length)
dat$flux <- exp(predict(flux.lm,newdata = data.frame(Length=dat$Length)))
ang.flux.lm.1 <- glm(log(flux)~Angle,dat,family=gaussian)
ang.flux.lm.2 <- glm(log(flux)~Angle+Y,dat,family=gaussian)
ang.flux.lm.3 <- glm(log(flux)~Y,dat,family=gaussian)
summary(ang.flux.lm.1)
summary(ang.flux.lm.2)
summary(ang.flux.lm.3)

dat[,ang2:=Angle-90]
dat[,pos:=ifelse(ang2>-45 & ang2<45,"bottom","other")]

qplot(data=dat,x=pos,y=log(flux))+geom_boxplot()
summary(aov(ang2~pos,dat))

qplot(data=dat,x=ang2,log(flux))+geom_smooth(method="lm") #like overall this plot
qplot(data=dat,x=Y,y=log(flux),col=spec)+geom_smooth(method="lm") 
b.ang <- seq(-180,180,22.5) #22.5 deg bins
dat$bins.ang <- with(dat, cut(ang2, breaks=b.ang))
levels(dat$bins.ang) <- as.factor(round(b.ang,2))

dat$plot <- "data"
blank.dat <- tail(dat,17) #change if not 22.5 deg bins
blank.dat$plot <- "blank"
blank.dat$sp <- "Tdentex"

blank.dat$bins.ang <- as.factor(seq(-180,180,22.5))
blank.dat$log.flux <- 0


#plot.dat <- dat[complete.cases(dat),]
plot.dat <- copy(dat) 


#plot.dat <- data.table(rbind(plot.dat,blank.dat))
plot.dat <- plot.dat[,.(f.sum=sum(flux),f.den=sum(flux)/(diff(range(X))*diff(range(Y)))),by=.(bins.ang,sp,sp2,sp3,spec)] #density per mm
plot.dat <- plot.dat[,.(f.mean=mean(f.sum),f.mean.dens=mean(f.den)),by=.(bins.ang,sp,sp2,sp3)]
plot.dat <- plot.dat[f.mean!=0]
plot.dat$log.flux <- log(plot.dat$f.mean)

sum.dat <- plot.dat[,.(f.sum=sum(log.flux)),by=.(sp)]

plot.dat <- merge(plot.dat,sum.dat,all.x=T,by="sp")

plot.dat <- plot.dat[,flux.s:=(log.flux/f.sum),by=.(bins.ang,sp)]

#polymetme has large photophores counted; remove
plot.dat[sp=="Polymetme thaeocoryla",c("f.mean","f.mean.dens","log.flux","f.sum","flux.s"):=list(0,0,0,0,0)]

#plot.dat <- merge(plot.dat,dens.m[,c("sp","sp2","sp3","dens.m")],by="sp") #species names

#add sp to plot.dat that don't have acc photophores and remove those from tree that we don't have data
null.dat <- tail(plot.dat,1)
null.dat[,c("f.mean","f.mean.dens","log.flux","f.sum","flux.s"):=list(0,0,0,0,0)]
plot.dat <-rbind(plot.dat,null.dat) 
plot.dat[.N,c("sp3","sp2","sp"):=list("Sigmops","Sigmops elongatus","Selongatus")]
plot.dat <-rbind(plot.dat,null.dat) 
plot.dat[.N,c("sp3","sp2","sp"):=list("Argyropelecus" ,"Argyropelecus aculeatus","Aaculeatus")]
plot.dat <-rbind(plot.dat,null.dat) 
plot.dat[.N,c("sp3","sp2","sp"):=list("Ichthyococcus" ,"Ichtyococcus ovatus","Iovatus")]


plot.dat$gen <-  plot.dat$sp3

none <- setdiff(tree2$tip.label,plot.dat$gen)


tree3 <- drop.tip(tree2,none)
setdiff(tree3$tip.label,plot.dat$gen)

#plot.dat <-plot.dat[sp!="Iovatus"&sp!="Vnimbaria"&sp!="Pthaeocoryla"] 

#### depth, downward flux PIC

dat[,down:=ifelse(ang2>-45&ang2<45,"down","other")]

dat$gen <-  gsub("(\\w+) .+","\\1",dat$sp2)

null.dat <- tail(dat,1)
null.dat[,c("log.L","flux"):=list(0,0)]
dat <-rbind(dat ,null.dat) 
dat[.N,c("sp2","sp","gen"):=list("Sigmops elongatus","Selongatus","Sigmops")]

dat <-rbind(dat,null.dat) 
dat[.N,c("sp2","sp","gen"):=list("Argyropelecus aculeatus","Aaculeatus","Argyropelecus")]

dat <-rbind(dat,null.dat)
dat[.N,c("sp2","sp","gen"):=list("Ichthyococcus ovatus","Iovatus","Ichthyococcus")]


#polymetme only large photophores
dat[gen=="Polymetme",c("log.L","flux"):=list(0,0)]


none <- setdiff(tree2$tip.label,dat$gen)

setdiff(dat$gen,tree2$tip.label)
tree3 <- drop.tip(tree2,none)
plot(tree3)

dat <- dat[gen%in%tree3$tip.label]

#flux in down position (ventral flux)

flux.pos <- dat[,.(v.flux=sum(flux[pos=="bottom"]/sum(flux))),by=.(sp,sp2,gen,spec)]
flux.pos[,d.flux:=1-v.flux]

flux.m <- flux.pos[,.(m.vflux=mean(v.flux),se.vflux=se(v.flux),m.dflux=mean(d.flux),se.dflux=se(d.flux)),by=gen]


flux.m[is.na(m.vflux),c("m.vflux","se.vflux","m.dflux","se.dflux"):=list(0,0,0,0)]

#qplot(data=flux.m,x=gen,y=m.dflux)

vflux.plot <- flux.m$m.vflux
dflux.plot <- flux.m$m.dflux
names(vflux.plot) <- flux.m$gen
names(dflux.plot) <- flux.m$gen


######################################## 
##### data for angle radial plot on phy ######
########################################



###ggtree with flux direction
plot.dat <-plot.dat[sp=="Iovatus"|sp=="Vnimbaria"|sp=="Pthaeocoryla",flux.s:=0]


#### tip data

tip.dat <- data.table(psub$data[psub$data$isTip==T,])
tip.dat <- tip.dat[order(tip.dat$parent),]

tip.dat$y.pos  <- 1-tip.dat$y/(max(tip.dat$y)+1)
tip.dat$x.pos  <- 0.91
tip.dat$sp <- tip.dat$label


plot.dat$sp <- factor(plot.dat$gen,levels = tip.dat$sp[order(tip.dat$y.pos,decreasing = T)])
saveRDS(plot.dat,"plot.dat.RDS")

#plot.dat <- readRDS("plot.dat.RDS") #read plot dat from 9Mar20

plot.dat[,bins.ang:=as.numeric(as.character(bins.ang))]
plot.dat <- plot.dat[!duplicated(plot.dat)]
plot.dat[,sp3:=gsub("(\\w).+ (\\w+)","\\1\\. \\2",sp2)]
plot.dat[log.flux==0,log.flux:=NA]


######################################## 
######## size distributon data ########
########################################

bd <- readRDS("body.depth.RDS")
dat.dist <- copy( dat)

dat.dist[,Y2:=bd-max(Y)+Y,.(sp,spec,region)]
#dat.dist[,scale.y:=(Y-min(Y))/max(Y-min(Y)),by=.(sp,spec,region)]
dat.dist[,scale.y:=(Y2/bd),by=.(sp,spec,region)]
dat.dist[,y.bin:=cut(scale.y,breaks=seq(0,1,0.1),labels = seq(10,100,10)),by=.(sp,spec)]
dat.dist[,logL:=log(Length*1000)]

#remove large photophores

dat.dist[, q:=between(Length,quantile(Length,probs =c(0.01,0.8))[1],quantile(Length,probs =c(0.01,0.8))[2]),,by=.(sp,spec)]


dat.dist.m <- dat.dist[q==T,.(m=mean(Length*1000),se=se(Length*1000)),by=.(sp,y.bin)]
dat.dist.m[]

dat.dist.m$sp2 <- gsub("(\\w)(\\w+)","\\1\\. \\2",dat.dist.m$sp)


#### work body depth into scale.y calculation
dat.dist <- merge(dat.dist,bd,by.x=c("sp2","spec"),by.y=c("species","spec"))
dat.dist[,scale.y2:=(diff(range(Y))/bd),by=.(sp,spec,region)]

dat.dist[,.(min=min(scale.y2),max=max(scale.y2)),by=.(sp,spec,region)]

### prepare angle data for plot

ang.dist.m <- dat.dist[q==T,.(m=mean(ang2),se=se(ang2)),by=.(sp,y.bin)]
ang.dist.m$sp2 <- gsub("(\\w)(\\w+)","\\1\\. \\2",ang.dist.m$sp)



#AIC, p, and r^ values of lm and log lm models


mod.comb <- function(dts){
  
  dt2 <- setDT(unlist(dts[-1], recursive = FALSE), check.names = F)[]
  
  dt3 <- dt2[, .SD, .SDcols =! names(dt2) %like% ".gen|.sp"]
  dt4 <- cbind(dts$AICd$gen,dt3)
  return(dt4)
}

anova.dt <- function(m1,m2){
  anova.m <- anova(m1,m2)
  return(anova.m$`Pr(>F)`)
}

#use "z" to ..... account for specimen, lm and exp p values include z as factor
size.aic <- mod.fit2(dt=dat.dist,x="scale.y",y = "Length",z="spec",by=c("sp")) 
size.aic2 <- mod.comb(size.aic)

size.aic3 <- cbind(size.aic$AICd$sp,size.aic2)
saveRDS(size.aic3,"size.AIC.models.RDS")


dat.dist[,q:=between(Length,quantile(Length,probs =c(0.01,0.8))[1],quantile(Length,probs =c(0.01,0.8))[2]),by=.(sp) ]


AICc.dt <- function(x){
  unlist(MuMIn::AICc(x)[,2])
}

###lmer p values (fraught), but see: https://stats.stackexchange.com/questions/95054/how-to-get-an-overall-p-value-and-effect-size-for-a-categorical-factor-in-a-mi
###uses lmertest package

lmer.p <- function(x){
  t <- anova(x)
  return(t$`Pr(>F)`)
}




dat.dist2 <- copy(dat.dist[q==T &!sp%in% c("Selongatus","Aaculeatus","Iovatus","Pthaeocoryla")])

size.aic <- dat.dist2[q==T,.(lm.aic=MuMIn::AICc(lmer(Length~scale.y+(1|spec))),log.aic=MuMIn::AICc(lmer(log(Length+1e-4)~scale.y+(1|spec)))), by=.(sp)]

size.aic[,c("lm.w","log.w"):=round(exp(.SD*-0.5)/sum(exp(-0.5*.SD)),3),by=sp]
p.vals <- dat.dist2[q==T,.(lm.pval=lmer.p(lmer(Length~scale.y+(1|spec))),log.pval=lmer.p(lmer(log(Length+1e-4)~scale.y+(1|spec)))), by=.(sp)]

#from ?r2_nakagawa {performance}
#Marginal and conditional r-squared values for mixed models are calculated based on Nakagawa et al. 2017. For more details on the computation of the variances, see get_variance. 

#The marginal r-squared considers only the variance of the fixed effects, while the conditional r-squared takes both the fixed and random effects into account. The random effect variances are actually the mean random effect variances, thus the r-squared value is also appropriate for mixed models with random slopes or nested random effects (see Johnson 2014)

r2 <- dat.dist2[q==T,.(
  lm.r2=r2(lmer(Length~scale.y+(1|spec)))[2][1],
  log.r2=r2(lmer(log(Length+1e-4)~scale.y+(1|spec)))[2][1]
), by=.(sp)]
size.aic2 <- merge(size.aic,r2)

size.aic2 <- merge(size.aic2,p.vals)



dat.dist2[,ang3:=abs(ang2)]


ang.aic <- dat.dist2[q==T,.(lm.aic=MuMIn::AICc(lmer(ang3~scale.y+(1|spec))),log.aic=MuMIn::AICc(lmer(log(ang3+1e-4)~scale.y+(1|spec)))), by=.(sp)]

ang.aic[,c("lm.w","log.w"):=round(exp(.SD*-0.5)/sum(exp(-0.5*.SD)),3),by=sp]
p.vals <- dat.dist2[q==T,.(lm.pval=lmer.p(lmer(ang3~scale.y+(1|spec))),log.pval=lmer.p(lmer(log(ang3+1e-4)~scale.y+(1|spec)))), by=.(sp)]


r2 <- dat.dist2[q==T,.(
  lm.r2=r2(lmer(ang3~scale.y+(1|spec)))[2][1],
  log.r2=r2(lmer(log(ang3+1e-4)~scale.y+(1|spec)))[2][1]
), by=.(sp)]
ang.aic2 <- merge(ang.aic,r2)

ang.aic2 <- merge(ang.aic2,p.vals)


pdf("ang.plots.pdf")
for(i in unique(dat.dist2$sp) ) {
  p1 <- qplot(data=dat.dist2[sp==i],x=scale.y,y=ang3,col=spec)+geom_smooth(method="lm")+ggtitle(paste0(i," lm p=",round(ang.aic2[sp==i]$lm.pval,4)))
  p2 <- qplot(data=dat.dist2[sp==i],x=scale.y,y=log(ang3),col=spec)+geom_smooth(method="lm")+ggtitle(paste0(i," log p=", round(ang.aic2[sp==i]$log.pval,4)))
  print(p1)
  print(p2)
}
dev.off()


pdf("size.plots.pdf")
for(i in unique(dat.dist2$sp) ) {
  p1 <- qplot(data=dat.dist2[sp==i],x=scale.y,y=Length,col=spec)+geom_smooth(method="lm")+ggtitle(paste0(i," lm p=",round(ang.aic2[sp==i]$lm.pval,4)))
  p2 <- qplot(data=dat.dist2[sp==i],x=scale.y,y=log(Length),col=spec)+geom_smooth(method="lm")+ggtitle(paste0(i," log p=", round(ang.aic2[sp==i]$log.pval,4)))
  print(p1)
  print(p2)
}
dev.off()

phot.no <- dat.dist[sp%in%dat.dist2$sp,.(n=.N),by=.(sp,spec)][,.(no.spec=length(unique(spec)),mean.n=round(mean(n),1),sd.n=round(sd(n),1)),by=sp]

size.aic2 <- merge(phot.no,size.aic2,by="sp")
saveRDS(size.aic2,"size.AIC.models.RDS")
ang.aic2 <- merge(phot.no,ang.aic2,by="sp")
saveRDS(ang.aic2,"ang.AIC.models.RDS")
t <- sapply(size.aic2,function(x) paste0(unlist(x)))
write.csv(t,"size.AIC.models.csv",row.names=F)
t2 <- sapply(ang.aic2,function(x) paste0(unlist(x)))
write.csv(t2,"ang.AIC.models.csv",row.names=F)

######################################## 
######## Flux density data ########
########################################


flux.dens <- dat[,.(f.den=sum(log(flux))/(diff(range(X))*diff(range(Y)))),by=.(sp,sp2,spec)] #log flux

flux.dens.m <- flux.dens[,.(dens.m=mean(f.den),se=se(f.den)),by=.(sp,sp2)]
flux.dens.m$sp4 <- gsub("(\\w)(\\w+)","\\1\\. \\2",flux.dens.m$sp)
flux.dens.m[dens.m==-Inf,c("dens.m","se"):=list(0,0)]

flux.dens.m$sp4 <- factor(flux.dens.m$sp4,levels=flux.dens.m$sp4[order(flux.dens.m$dens.m)])


######################################## 
######## prepare size data ########
########################################

dat[, q:=between(Length,quantile(Length,probs =c(0.01,0.8))[1],quantile(Length,probs =c(0.01,0.8))[2]),,by=.(sp,sp2,spec)]

size.m <- dat[q==T& log.L!=-Inf,.(sizeL.m=mean(log(Length*1000),na.rm=T),seL=se(log(Length*1000)),size.m=mean(Length*1000,na.rm=T),se=se(Length*1000),sl.m=mean(sl,na.rm=T)),by=.(sp,sp2)]

size.m$sp4 <- gsub("(\\w)(\\w+)","\\1\\. \\2",size.m$sp)
size.m[sp%in%c("Pthaeocoryla","Selongatus","Aaculeatus","Iovatus"),c("size.m","se","sizeL.m","seL"):=list(0,0,0,0)]
size.m$sp4 <- factor(size.m$sp4,levels=size.m$sp4[order(size.m$size.m)])


spec.sum <- phot.dat[,.(n=length(X)),by=.(sp2,spec)][,.(n=length(spec),mean.n=mean(n),sd.n=sd(n)),by=sp2]

nulls <- data.table(sp2=flux.dens.m[dens.m==0]$sp2, n=3,mean.n=0,sd.n=0)

spec.sum <- rbind(spec.sum,nulls)
sl <- sl.dat[,.(sl=paste(range(round(sl,1)),collapse= "—")),by=.(species)]
spec.sum <- merge(spec.sum,sl,by.x="sp2",by.y="species")



write.csv(spec.sum,"specimen.summary.csv")


bd <- data.table(read.csv("body.depth.csv"))
bd3 <- merge(sl.dat,bd,by.x="species",by.y="sp2")

bd2<- bd3[,.(bd=paste(range(round(sl*BD,1)),collapse= "—")),by=.(species)]

bd3 <- bd3[,.(bd=round(sl*BD,1)),by=.(species,spec)]

saveRDS(bd3,"body.depth.RDS")

#save data
#saveRDS(phot.dat,"phot.dat.RDS")
#saveRDS(dat,"size.dat.RDS")
