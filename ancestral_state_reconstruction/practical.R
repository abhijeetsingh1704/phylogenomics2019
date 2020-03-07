#prepare desktop

library(phytools)
library(car)
library(mnormt)
library(caper)
library(geiger)
library(nlme)
library(phylolm)
libray(phangorn)
libray(castor)
library(diversitree)
library(mvSLOUCH)

#########################
####CONTINIUOS TRAITS####
#########################

# for futher simulations see:
# https://lukejharmon.github.io/ilhabela/instruction/2015/07/02/simulating-Brownian-motion/

#check how covariance and variance is generated along the tree and at the tip
tree<-pbtree(n=6,tip.label=LETTERS[1:6])
X<-fastBM(tree,nsim=100)
plot(rotateNodes(tree,"all"),edge.width=2,no.margin=TRUE)
scatterplotMatrix(t(X))
hist(t(X)[,"A"]) 
hist(X[,1])

#MODELLING TRAIT EVOLUTION
nspec=200
tree<-pbtree(n=nspec,scale=100)
x1<-fastBM(tree,a=2,sig2=0.01, internal=TRUE)
x2<-fastBM(tree,a=2,sig2=0.01, mu=0.05, internal=TRUE)
x3<-fastBM(tree,a=2,theta=2,alpha=0.8,sig2=0.1, internal=TRUE, model = "OU")

par(mfrow=c(1,3))
phenogram(tree,x1, spread.labels=F, ylim=c(0, max(as.numeric(x3),as.numeric(x2),as.numeric(x1))), main="BM")
phenogram(tree,x2, spread.labels=F, ylim=c(0, max(as.numeric(x3),as.numeric(x2),as.numeric(x1))), main="BM with trend")
phenogram(tree,x3, spread.labels=F, ylim=c(0, max(as.numeric(x3),as.numeric(x2),as.numeric(x1))), main="OU")

#check out how different models fit with the simulated data
fit1=fitContinuous(tree,x1[1:nspec])
fit2=fitContinuous(tree,x1[1:nspec],model="trend")
fit3=fitContinuous(tree,x1[1:nspec],model="OU")
setNames(c(fit1$opt$aicc,fit2$opt$aicc,fit3$opt$aicc),c("BM","trend","OU"))


###PAGEL'S PARAMETERS (try to play with different n. of species)
nspec=100
tree<-pbtree(n=nspec,scale=100)

##transform lambda
tree1=rescale(tree, "lambda", 0)
tree2=rescale(tree, "lambda", 0.5)

#plot
layout(matrix(1:3,1,3))
plot(tree1, edge.width=1.5, cex=1.2); mtext("Lambda = 0")
plot(tree2, edge.width=1.5, cex=1.2); mtext("Lambda = 0.5")
plot(tree, edge.width=1.5, cex=1.2); mtext("Lambda = 1 (default)")

#simulate data
x<-fastBM(tree,a=0,sig2=0.01, internal=F)
x1<-fastBM(tree1,a=0,sig2=0.01, internal=F)
x2<-fastBM(tree2,a=0,sig2=0.01, internal=F)

#estimate lambda for the simulated data
fitContinuous(tree,x, model="lambda")$opt$lambda
fitContinuous(tree,x1, model="lambda")$opt$lambda
fitContinuous(tree,x2, model="lambda")$opt$lambda

#check likelihood surface
xdata=data.frame(x)
xdata$species=rownames(xdata)
tree.x=comparative.data(phy = tree, data = xdata, names.col =species)
res<-pgls(x~ 1, data = tree.x, lambda='ML')

xdata1=data.frame(x1)
xdata1$species=rownames(xdata1)
tree.x1=comparative.data(phy = tree, data = xdata1, names.col =species)
res1<-pgls(x1~ 1, data = tree.x, lambda='ML')

xdata2=data.frame(x2)
xdata2$species=rownames(xdata2)
tree.x2=comparative.data(phy = tree, data = xdata2, names.col =species)
res2<-pgls(x2~ 1, data = tree.x2, lambda='ML')

layout(matrix(1:3,1,3))
plot(pgls.profile(res, 'lambda'))
plot(pgls.profile(res1, 'lambda'))
plot(pgls.profile(res2, 'lambda'))

#estimate phylogenetic signal
phylosig(tree,x,method="lambda")$lambda
phylosig(tree,x1,method="lambda")$lambda
phylosig(tree,x2,method="lambda")$lambda

phylosig(tree,x,method="K", test=T)
phylosig(tree,x1,method="K", test=T)
phylosig(tree,x2,method="K", test=T)


##transform delta & kappa and repeat the above (try to simulate data and fit models as above)
tree1=rescale(tree, "delta", 0.1)
tree2=rescale(tree, "delta", 2)

layout(matrix(1:3,1,3))
plot(tree1, edge.width=1.5, cex=1.2); mtext("delta = 0.1")
plot(tree, edge.width=1.5, cex=1.2); mtext("delta = 1 (default)")
plot(tree2, edge.width=1.5, cex=1.2); mtext("delta = 2")

tree1=rescale(tree, "kappa", 0)
tree2=rescale(tree, "kappa", 2)

layout(matrix(1:3,1,3))
plot(tree1, edge.width=1.5, cex=1.2); mtext("kappa = 0")
plot(tree, edge.width=1.5, cex=1.2); mtext("kappa = 1 (default)")
plot(tree2, edge.width=1.5, cex=1.2); mtext("kappa = 2")




#########################
####DISCRETE   TRAITS####
#########################



nspec=100
tree<-pbtree(n=nspec,scale=100)
x1<-rTraitDisc(tree, model = "ER")
x2<-rTraitDisc(tree, model = "SYM")
x3<-rTraitDisc(tree, model = "ARD", rate=c(0.1,0.4))

fit1<-fitDiscrete(tree,x1,model="ER")
fit2<-fitDiscrete(tree,x1,model="SYM")
fit3<-fitDiscrete(tree,x1,model="ARD")
fit4<-fitDiscrete(tree,x1,model=matrix(c(0,0,1,0),2,2))

layout(matrix(1:4,1,4))
plot(fit1, signif=5, main=paste("ER model, LL =", round(fit1$opt$lnL,3)))
plot(fit2, signif=5, main=paste("SYM model, LL =", round(fit2$opt$lnL,3)))
plot(fit3, signif=5, main=paste("ARD model, LL =", round(fit3$opt$lnL,3)))
plot(fit4, signif=5, main=paste("IRR model, LL =", round(fit4$opt$lnL,3)))

aicc<-setNames(
    c(fit1$opt$aicc,fit2$opt$aicc,fit3$opt$aicc,fit4$opt$aicc),
    c("ER","SYM","ARD","IRR"))
aic.w(aicc)

LRT=2*(fit1$opt$lnL-fit3$opt$lnL)
pchisq(abs(LRT), 1, lower.tail=FALSE)

LRT=2*(fit1$opt$lnL-fit4$opt$lnL)
pchisq(abs(LRT), 1, lower.tail=FALSE)


#more than 2 states (you can also use fitMk from phytools)
k=3
tree<-pbtree(n=nspec,scale=100)
x1<-rTraitDisc(tree, k=k,model = "ER")
x2<-rTraitDisc(tree, k=k, model = "SYM", rate=sample(seq(0.1, 1, by=0.1), k*(k - 1)/2))
x3<-rTraitDisc(tree, k=k, model = "ARD", rate=sample(seq(0.1, 1, by=0.1), k*(k - 1)))

Q = castor::get_random_mk_transition_matrix(k, rate_model="ER")
tip_states = LETTERS[simulate_mk_model(tree, Q)$tip_states]
names(tip_states)=tree$tip.label

fit1<-fitDiscrete(tree,x1,model="ER")
fit2<-fitDiscrete(tree,x2,model="SYM")
fit3<-fitDiscrete(tree,x3,model="ARD")

layout(matrix(1:3,1,3))
plot(fit1, signif=5, main="ER")
plot(fit2, signif=5, main="SYM")
plot(fit3, signif=5, main="ARD")


#try larger k and define an ordered matrix
model<-matrix(c(0,1,0,0,
    2,0,3,0,
    0,4,0,5,
    0,0,6,0),4,4)
x1<-rTraitDisc(tree, k=4,model = "ER")    
fit1<-fitDiscrete(tree,x1,model="ER")
fit4<-fitDiscrete(tree,x1,model=model)

layout(matrix(1:2,1,2))
plot(fit1, signif=5, main="ER")
plot(fit2, signif=5, main="oredered")

LRT=2*(fit1$opt$lnL-fit4$opt$lnL)
pchisq(abs(LRT), 1, lower.tail=FALSE)



######################################
####ANCESTRAL TRAIT RECONSTRUCTION####
######################################

#continiuous characters
#unknown ancestral states
nspec=100
tree<-pbtree(n=nspec,scale=100)
x1<-fastBM(tree,a=2,sig2=0.01, internal=F)

fit<-fastAnc(tree,x1,vars=TRUE,CI=TRUE)
print(fit,printlen=10)

range(x1)
fit$CI95[1,]

## projection of the reconstruction onto the edges of the tree
to_plot<-contMap(tree,x1,plot=F)
plot(to_plot,type="fan",legend=0.7*max(nodeHeights(tree)),sig=2,fsize=c(0.7,0.9))

#known ancestral states
nspec=100
tree<-pbtree(n=nspec,scale=1)
x1<-fastBM(tree,internal=TRUE)
a<-x1[as.character(1:tree$Nnode+Ntip(tree))]
x1<-x1[tree$tip.label]

fit=ace(x1, tree)

plot(a,fit$ace,xlab="True states",ylab="Estimated states",bg="black",cex=1.4,pch=21)
lines(range(c(x,a)),range(c(x,a)),lty="dashed",col="red") ## 1:1 line
segments(x0=a, y0=fit$CI95[,1], x1=a, y1=fit$CI95[,2])


#primate data
tree=read.tree("~/Desktop/data/primate/primate_tree.phy")
xdata=read.table(file="~/Desktop/data/primate/primate_spec.txt",header=T,sep="\t",fill=T, strip.white=T, dec = ".")

rownames(xdata)=xdata$species
xdata=xdata[tree$tip.label,]
x=xdata$lg.body
names(x)=rownames(xdata)

to_plot<-contMap(tree,x,plot=F)
plot(to_plot,type="fan",legend=0.7*max(nodeHeights(tree)),sig=2,fsize=c(0.7,0.9))



#discrete characters
trees=read.nexus("~/Desktop/data/Estrildid/23755.nex")
xdata=read.table(file="~/Desktop/data/Estrildid/dt(85species).txt",header=T,sep="\t",fill=T, strip.white=T, dec = ".")

tree=trees[[1]]

rownames(xdata)=xdata$name
xdata=xdata[tree$tip.label,]
x=xdata$Fsongbi
names(x)=rownames(xdata)

plotTree(tree, type = "fan", fsize = 0.7, ftype="i")
cols <- setNames(palette()[1:length(unique(x))], sort(unique(x)))
tiplabels(pie = x, piecol = cols, cex = 0.2)

fit=ace(x, tree, model = "ER")
#fit=rerootingMethod(tree, x, model = "ER") marginal likelihoods
nodelabels(node = as.numeric(names(fit$ace)), pie = fit$ace, piecol = cols, cex = 0.5)
add.simmap.legend(leg= c("song", "no song"),colors = cols, prompt = FALSE, x = 0.9 * par()$usr[1], y = -max(nodeHeights(tree)))

mtree <- make.simmap(tree, x, model = "ER")
plotSimmap(mtree, cols, type = "fan", fsize = 0.9, ftype = "i")

mtrees <- make.simmap(tree, x, model = "ER", nsim = 25)
par(mfrow = c(5, 5))
plotSimmap(mtrees, cols, lwd = 1, ftype = "off")
dev.off()

pd <- describe.simmap(mtrees, plot = F)
plotSimmap(mtree, cols, type = "fan", fsize = 0.9, ftype = "i")
nodelabels(pie = pd$ace, piecol = cols, cex = 0.5)
add.simmap.legend(leg= c("no song", "song"),colors = cols, prompt = FALSE, x = 0.9 * par()$usr[1], y = -max(nodeHeights(tree)))



#import Estrildid data $ choose a character and create a vector
trees=read.nexus("~/Desktop/data/Estrildid/23755.nex")
xdata=read.table(file="~/Desktop/data/Estrildid/dt(85species).txt",header=T,sep="\t",fill=T, strip.white=T, dec = ".")
tree=trees[[1]]

rownames(xdata)=xdata$name
xdata=xdata[tree$tip.label,]
x=xdata$Fsongbi
names(x)=rownames(xdata)

#import primate data & choose a character and create a vector
tree=read.tree("~/Desktop/data/primate/primate_tree.phy")
xdata=read.table(file="~/Desktop/data/primate/primate_spec.txt",header=T,sep="\t",fill=T, strip.white=T, dec = ".")

rownames(xdata)=xdata$species
xdata=xdata[tree$tip.label,]
x=xdata$lg.body
names(x)=rownames(xdata)

#import Stenus beetle data & choose a character and create a vector
tree=read.nexus("~/Desktop/data/Stenus_data/14species_tree.nex")
xdata<-read.table(pipe("pbpaste"), header=TRUE, sep="\t",fill=T, strip.white=T) #copy from Excel
xdata <-read.table("clipboard", header=T, sep="\t",fill=T, strip.white=T) #for Windows

rownames(xdata)=xdata$species
xdata=xdata[tree$tip.label,]
x=xdata$BM
names(x)=rownames(xdata)

library(knitr)
setwd("~/Desktop/data")
knit("parctical_garamszegi.html" , "practical_garamszegi_out.html" )


