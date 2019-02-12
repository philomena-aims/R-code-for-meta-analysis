#R version 3.4.4 (2018-03-15) -- "Someone to Lean On"
#Copyright (C) 2018 The R Foundation for Statistical Computing
#Platform: x86_64-w64-mingw32/x64 (64-bit)

#load the packages "MBESS","meta", "effsize","xtable", "orddom"
##################################################################

#setting the seed so that we can replicate the same estimates
set.seed(1000) 
#we define the values of ES, mu1,mu2 sd1, and sd2 
#such that it satisfy the equation ES=(mu1-mu2)/sd1 or sd2)
ES=0.8 #pre specified True effect size
mx=80  # population mean of sample x
sdx=7  #population standard deviation of sample x
my=74.4 #population mean of sample y
sdy=7 #population standard deviation of sample y
####################################################################

nsims<-10000 #number of simulated experiments

####################################################################

#setting up empty containers for the simulations

#set up empty container for the simulated Effect Size (cohen's d)
cd.es<-numeric(nsims) 
#set up empty container for the simulated Effect Size (Hedge's g)
hg.es<-numeric(nsims)
#set up empty container for the simulated Effect size (Glass's d)
G.es<-numeric(nsims)
#set up empty container for the simulted Effect size (Cliff d)
C.es<-numeric(nsims)
# set up empty container for simulated Effect size 
#(Probability of Superiority PS)
E.es<-numeric(nsims)
#set up empty container for simulated variance of 
#Cliff d effect size
C.es.v<-numeric(nsims)
#set up empty container for simulated Mann Whitney statistic 
f<-numeric(nsims)

# A function that take seven arguments to generate
#contaminated normal samples
rcnorm <-function(n,mu1,mu2,sig1, sig2,cpct){
  zo <- rnorm(n, mean = mu1, sd = sig1)
  z1 <- rnorm(n, mean = mu2, sd = sig2)
  flag <- rbinom(n,size=1, prob=cpct)
  z <- zo*(1-flag) + z1*flag
  return(z)
}

for(i in 1:nsims){ # for each simulated experiment
  samplesize1<-120   #c(10,15,30,40,120) #sample size of study 1
  samplesize2<-80   #c(8,20,18,50,80) #sample size of study 2
  
  #produce simulated participants from the normal distribution
  x<-rnorm(n=samplesize1, mean = mx, sd=sdx) 
  y<-rnorm(n=samplesize2, mean=my, sd=sdy)
  
  #produce simulated participants from a contaminated normal distribution
  #cpct is the percetage of contamination expressed in a fraction 
  #x<-rcnorm(n=samplesize1, mx, mx-10, sdx, sdx-0.5, cpct= 0.10)
  #y<-rcnorm(n=samplesize2, my, my-10, sdy, sdy-0.5, cpct=0.10)
  
  ########################################################################
  #Estimating Effect size
  
  #Use MBESS to calc d unbiased(Cohen's d)
  cd.es[i]<-smd(Mean.1=mean(x),Mean.2=mean(y),s.1=sd(x), 
                s.2=sd(y),n.1=samplesize1,n.2=samplesize2,Unbiased=FALSE) 
  
  #Use MBESS to calc d biased(Hedges' g)
  hg.es[i]<-smd(Mean.1=mean(x), Mean.2=mean(y), s.1=sd(x), 
                s.2=sd(y), n.1=samplesize1, n.2=samplesize2, Unbiased=TRUE) 
  
  #Use MBESS to calc Glass (Glass's g)
  G.es[i]<-smd.c(Group.T=x, Group.C = y) 
  
  #use effsize to calc Cliff delta effect size
  C.es[i]<-cliff.delta(x,y, return.dm = TRUE)$estimate
  C.es.v[i] <-cliff.delta(x,y, return.dm = TRUE)$var
  
  #calc the Mann Whitney U test
  f[i]<- wilcox.test(x,y)$statistic
  #estimating the Probability of Superiority effect size 
  E.es[i]<-f[i]/(samplesize1*samplesize2)
}

##################################################################
#Estimating Standard deviation of the various effect sizes

#standard deviation of cohen's d
cd.se<-sqrt((((samplesize1+samplesize2)/samplesize1*samplesize2))
            +(cd.es^2/(2*(samplesize1+samplesize2)))) 

#standard deviation of Hedge's g
J<-1-3/(4*(samplesize1+samplesize2-2)-1) #correction for bias
hg.se<-sqrt((((samplesize1+samplesize2)/samplesize1*samplesize2))
            +(cd.es^2/(2*(samplesize1+samplesize2)))*J^2) 

#Standard deviation of Glass effect size
G.se<-sqrt((((samplesize1+samplesize2)/samplesize1*samplesize2))
           +(G.es^2/(2*(samplesize2-1))))

#standard deviation of Mann Whitney
E.se<-rep(sqrt((samplesize1*samplesize2*(samplesize1 +
          samplesize2 + 1))/12)/(samplesize1*samplesize2),nsims)

# standard deviation of cliff
C.se <- sqrt(C.es.v)

####################################################################
# A function that take two argument to do a Meta-Analysis

metafun <-function(Efs, EFsd){
  metat<-metagen(Efs,EFsd)
  return(metat)
}  

##################################################################
#creating a dataframe for the collectors

#table of effect size
data <- data.frame(n1=samplesize1, n2=samplesize2, 
                    cohensd=metafun(cd.es,cd.se)$TE.fixed,
                    Hedgesg = metafun(hg.es, hg.se)$TE.fixed,
                    Glassd = metafun(G.es, G.se)$TE.fixed,
                    CliffsD = metafun(C.es, C.se)$TE.fixed,
                    MannU = metafun(E.es, E.se)$TE.fixed)

#Converted cliff delta and PS Effects Size

dataa<-data.frame(n1=samplesize1, n2=samplesize2, 
                   Cohens=data$cohensd, Hedges=data$Hedgesg, Glass=data$Glassd, 
                   Cliff=delta2cohd(data$CliffsD),
                   PS=qnorm(data$MannU)*sqrt(2))



##########################################################################################
#table for standard deviation
data1<-data.frame(n1=samplesize1, n2=samplesize2, 
                   cohensd=metafun(cd.es,cd.se)$seTE.fixed,
                   Hedgesg = metafun(hg.es, hg.se)$seTE.fixed,
                   Glassd = metafun(G.es, G.se)$seTE.fixed,
                   CliffsD = metafun(C.es, C.se)$seTE.fixed,
                   MannU = metafun(E.es, E.se)$seTE.fixed)


data11<-data.frame(n1=samplesize1, n2=samplesize2, 
                   Cohens=data1$cohensd, Hedges=data1$Hedgesg, Glass=data1$Glassd,
                   Cliff=delta2cohd(data1$CliffsD),
                   PS=qnorm(data1$MannU)*sqrt(2))






#############################################################################################

#confidence interval
data2 <- data.frame(n1=samplesize1, n2=samplesize2, 
                    cohensdL=metafun(cd.es,cd.se)$lower.fixed,
                    
                    cohensdU= metafun(cd.es,cd.se)$upper.fixed,
                    
                    HedgesgL = metafun(hg.es, hg.se)$lower.fixed,
                    
                    HedgesgU=metafun(hg.es, hg.se)$upper.fixed,
                    GlassdL = metafun(G.es, G.se)$lower.fixed,
                    
                    GlassdU= metafun(G.es, G.se)$upper.fixed,
                    
                    CliffsDL = metafun(C.es,C.se)$lower.fixed,
                    
                    CliffsDU= metafun(C.es,C.se)$upper.fixed,
                    
                    MannUL = metafun(E.es, E.se)$lower.fixed,
                    
                    MannUU= metafun(E.es, E.se)$upper.fixed)


data22 <- data.frame(n1=samplesize1, n2=samplesize2, 
                    CohensL=data2$cohensdL, CohensU=data2$cohensdU, hedgesL=data2$HedgesgL, 
                    hedgesU=data2$HedgesgU, GlassL=data2$GlassdL, 
                    GlassU=data2$GlassdU, CliffL=delta2cohd(data2$CliffsDL),
                    CliffU=delta2cohd(data2$CliffsDU), PSL=qnorm(data2$MannUL)*sqrt(2), 
                    PSU=qnorm(data2$MannUU)*sqrt(2))
                    
 

################################################################
# writing the dataframes to a tables in latex code.

#xdata<-xtable(data, caption = "effect size of studies " , 
              
              #digits=c(0,0,0,3,3,3,3,3))
#print(xdata)


xdataa<-xtable(dataa, caption = "Converted Effect size of Studies", 
               digits=c(0,0,0,3,3,3,3,3))
print(xdataa)


#xdata1<-xtable(data1, caption = "Standard deviation of Effects sizes " ,
              #digits=c(0,0,0,3,5,6,6,4))
#print(xdata1) 



xdata11<-xtable(data11, caption = "Standard deviation of Effects sizes " ,
               digits=c(0,0,0,3,5,6,6,4))
print(xdata11)


#xdata2<-xtable(data2, caption = "Confidence Intervals" ,
               
              # digits=c(0,0,0,3,3,3,3,3,3,3,3,3,3))
#print(xdata2)


xdata22<-xtable(data22, caption = "Confidence Intervals" ,
               
               digits=c(0,0,0,3,3,3,3,3,3,3,3,3,3))
print(xdata22)

##################################################################
