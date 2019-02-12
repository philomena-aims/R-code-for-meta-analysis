#load package orddom

################################################################

# specified Population effect size
ES=0.2 # Popu. effect size for cohen's, Hedges and Glass
ESC<-cohd2delta(ES)   # converted popu. effect size for cliff
ESP<-pnorm(ES/sqrt(2))    # converted popu. effect size for PS

####################################################################


#Meta-Analysis of the five studies for the overall Meta-Analysis, Cohens d

#M.es.cd is the Estimated mean effect size of Cohens d
M.es.cd<-c(29.830, 34.511, 28.571, 33.380, 28.873) # to be typed
n1=c(10,15,30,40,120) #sample sizes in group 1
n2<-c(8,20,18,50,80) #sample sizes in group 2

#calculating the variance
M.es.cd.v<-(((n1+n2)/n1*n2))+(M.es.cd^2/(2*(n1+n2)))

#calculating the standard error
M.cd.se<-c(sqrt(M.es.cd.v))

#performing Meta-Analysis
meta<-metagen(M.es.cd, M.cd.se, comb.fixed = T)
meta

###################################################################

#Analyses of Analysis

#Meta-Analysis of the five studies for the overall Meta-Analysis Hedges' g

#M.es.hg is the Estimated Mean effects size of Hedges' g
M.es.hg<-c(28.466, 33.739, 28.107, 33.095, 28.764)
n1=c(10,15,30,40,120)
n2<-c(8,20,18,50,80)

#Calculate Variance of the Hedges' g
J<-1-3/(4*(n1+n2-2)-1) #correction for bias
M.es.hg.v<-(((n1+n2)/n1*n2))+(M.es.cd^2/(2*(n1+n2)))*J^2

#Calculate Standard errors of Hedges' g
M.hg.se<-c(sqrt(M.es.hg.v))

meta1<-metagen(M.es.hg, M.hg.se, comb.fixed = T)
meta1

#collectors are the "statistical bias" and "standard error" 
#for the estimated effect 
#size and "confidance band width" and "confidence interval
#coverage" for the confidence
#interval

####################################################################

#Meta-Analysis of the five studies for the overall Meta-Analysis, Glass

#M.es.G is the Estimated Mean effects size of Glass
M.es.G<-c(88.347, 90.352, 89.902, 91.328, 91.270)

n1=c(10,15,30,40,120)
n2<-c(8,20,18,50,80)

#calculating the variance
M.es.G.v<-(((n1+n2)/n1*n2))+(M.es.G^2/(2*(n1-1)))

#calculating the standard error
M.G.se<-c(sqrt(M.es.G.v))

#performing Meta-Analysis
meta2<-metagen(M.es.G, M.G.se, comb.fixed = T)
meta2

###################################################################

#Meta-Analysis of the five studies for the overall Meta-Analysis, Cliff Delta

#M.es.C is the estimated Mean effect size of Cliff delta
M.es.C<-c(1.000, 1.000, 1.000, 1.000, 1.000)

#M.se.C is the estimated standard deviation of Cliff delta
M.se.C<-c(0.0001, 0.00001, 0.00001, 0.000001, 0.0000001)

#performing Meta-Analysis
meta3<-metagen(M.es.C, M.se.C, comb.fixed = T)
meta3

####################################################################

#Meta-Analysis of the five studies for the overall Meta-Analysis, Mannu

#M.es.M is the estimated mean effect size of Probability of Superiority
M.es.M<-c(1.000, 1.000, 1.000, 1.000, 1.000)
#M.se.M is the standard deviation of Probability of Superiority
M.se.M<-c(0.001, 0.001, 0.001, 0.001, 0.0004)

#performing Meta-Analysis
meta4<-metagen(M.es.M, M.se.M, comb.fixed = T)
meta4

###############################################################
#collectors

#Estimated Mean Effects size
data5<-data.frame(cohens=meta$TE.fixed, hedges=meta1$TE.fixed,
                  glass=meta2$TE.fixed, cliff=meta3$TE.fixed, mann=meta4$TE.fixed)
data5

#Standard deviation
data51<-data.frame(cohens=meta$seTE.fixed, hedges=meta1$seTE.fixed,
                   glass=meta2$seTE.fixed, cliff=meta3$seTE.fixed, mann=meta4$seTE.fixed)
data51

#confidence interval
data52<-data.frame(cohensL=meta$lower.fixed,cohensU=meta$upper.fixed,
                   
                   hedgesL=meta1$lower.fixed,hedgesU=meta1$upper.fixed,
                   GlassL=meta2$lower.fixed,GlassU=meta2$upper.fixed,
                   CliffL=meta3$lower.fixed, cliffU=meta3$upper.fixed,
                   MannL=meta4$lower.fixed,MannU=meta4$upper.fixed)
data52


##############################################################

#statistical bias (ES is the value of the population effect size in cohen's d)
data6<-data.frame(ES, Contamination=0, 
                  Cohensd=meta$TE.fixed-ES, Hedgesg=meta1$TE.fixed-ES,
                  Glassg=meta2$TE.fixed-ES, Cliff=delta2cohd(meta3$TE.fixed)-ES,
                  Mannu=(qnorm(meta4$TE.fixed)*sqrt(2))-ES )


#Standard Errors
data7<-data.frame(ES, Contamination=0, 
                  Cohensd=meta$seTE.fixed, Hedgesg=meta1$seTE.fixed,
                  Glassd=meta2$seTE.fixed, Cliif=delta2cohd(meta3$seTE.fixed), Mannu=qnorm(meta4$seTE.fixed)*sqrt(2))

#confidence interval width
data8 <- data.frame(ES, Contamination=0, 
                    Cohensd=(meta$upper.fixed-meta$lower.fixed),   
                    Hedgesg = (meta1$upper.fixed-meta1$lower.fixed),
                    Glassd = (meta2$upper.fixed-meta2$lower.fixed),
                    CliffsD = (delta2cohd(meta3$upper.fixed-meta3$lower.fixed)),
                    MannU = (qnorm(meta4$upper.fixed)*sqrt(2))-(qnorm(meta4$lower.fixed)*sqrt(2)))

#percentage error
data9 <- data.frame(ES,contamination=0,
                    Cohensd=((ES-meta$TE.fixed)/ES*100),
                    Hedgesg=((ES-meta1$TE.fixed)/ES*100),
                    Glassd=((ES-meta2$TE.fixed)/ES*100),
                    CliffsD=((delta2cohd((ESC-meta3$TE.fixed))/ES)*100),
                    MannU=(ES-(qnorm(meta4$TE.fixed)*sqrt(2))/ES)*100) 

#################################################################

#tables in latex
xdata6<-xtable(data6, caption = "Statistical bias",
               digits=c(0,2,0,3,3,3,3,3))
xdata6

xdata7<-xtable(data7, caption = "Standard errors", 
               digits=c(0,2,0,3,3,3,9,9))
xdata7

xdata8<-xtable(data8, caption = "Confidence interval width",
               digits=c(0,2,0,3,3,9,9,9))
xdata8

xdata9<-xtable(data9, caption = "Percentage error",
               digits=c(0,2,0,5,5,5,8,5))
xdata9

#######################################################################