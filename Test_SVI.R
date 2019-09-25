library(dplyr)
library(ggplot2)
library(reshape2)
library(ggplot2)
source('svi_h.R')

# read one sample ivs at tau = 0.12
tau = 0.12
c <- read.csv('iv.csv')
# measure estimation time
ptm <- proc.time()
param <- fit_ga_svi(c$k,c$iv,tau)
#show how long it takes for estimate the parameters
print(proc.time()-ptm)

ptm <- proc.time()
param <- fit_svih(c$k,c$iv,tau, 0.95)
#show how long it takes for estimate the parameters
print(proc.time()-ptm)

print(param)

# SVI parameterized IV estimation
c$svi <- svi_fun(c$k,as.numeric(param[1,2:7]))/sqrt(tau)

# plot fitting results
ggplot(c, aes(x = k, y = iv)) + geom_point()  + 
  geom_line(aes(y = svi), color = 'red')
c$err = c$iv - c$svi

ggplot(c, aes(x = k, y = err)) + geom_point()
print(sd(c$err))

# read a chain
library(parallel)
c2 <- read.csv('chain.csv')
# suppose now H = 0.65. This can be done by using genetic algo
ptm <- proc.time()
para <- par_fit_svi(c2, 0.65)
print(proc.time()-ptm)

print(para)



vol_svi<-matrix(0,dim(c2)[1],dim(para)[1])
Expiry<-unique(c2$tau)
strikes<-seq(-.2,.1,.01)

cbind(para,Expiry)
vol_surface<-svi_fun(strikes,as.numeric(para[1,2:7]))/sqrt(Expiry[1])
  
for(i in 2 : length(Expiry))
{
  vol_surface<-cbind(vol_surface,
        svi_fun(strikes,as.numeric(para[i,2:7]))/sqrt(Expiry[i]))
}
colnames(vol_surface)<-as.Date(as.Date(unique(c2[1,2]))+round(Expiry*360))
rownames(vol_surface)<-strikes

xymelt <- melt(vol_surface, id.vars=rownames(vol_surface))

xymelt$Var2<-as.Date(xymelt$Var2,origin="1970-01-01")


ggplot(xymelt, aes(x = Var1, y = value,  color=as.factor(Var2))) +geom_line()


### Specific Tenor

tenor<-2
c2_<-filter(c2,tau==Expiry[tenor])
param2 <- fit_svih(c2_$k,c2_$iv,c2_$tau, 0.95)
c2_$svi <- svi_fun(c2_$k,as.numeric(param2[1,2:7]))/sqrt(c2_$tau)
ggplot(c2_, aes(x = k, y = iv)) + geom_point()  + 
  geom_line(aes(y = svi), color = 'red')

c2_$err <- c2_$iv - c2_$svi

param3<-fit_ga_svi(c2_$k, c2_$iv, c2_$tau)
c2_$svi2 <- svi_fun(c2_$k,as.numeric(param3[1,2:7]))/sqrt(c2_$tau)
ggplot(c2_, aes(x = k, y = iv)) + geom_point()  + 
  geom_line(aes(y = svi2), color = 'red')

c2_$err_gen <- c2_$iv - c2_$svi2


ggplot(c2_, aes(x = k, y = iv)) + geom_point()  + 
  geom_line(aes(y = svi2,x=k), color = 'red')+
  geom_line(aes(y = svi,x=k), color = 'blue')+
  scale_color_discrete(name = "Fit", labels = c("Nelder-Mead", "Genetic"))


