#
# R-source for the figures for the Esk paper
#   Stephenson D.B., B. Casati, C.A.T. Ferro and C.A. Wilson, 2008: 
#   The extreme dependency score: a non-vanishing measure for forecasts of rare events.
#   Meteorol. Appl., 15, 41-50. 
#
# Author: B. Casati, Ouranos (Casati.Barbara@ouranos.ca)
#
-------------------------------------------------------------------------------------------------

#
# scan in data file with precipitation accumulations
#

inp = read.table("ESKDALEMUIR_T+06",header=T)

#
# eliminate missing values
#

esk = data.matrix(inp)
esk[esk<0] = NA
esk = na.omit(esk)

# esk has 6337 entries. after eliminaing NA, 6266 entries

#
# define Forecast and Observation
#

Obs = esk[,2]
For = esk[,3]

#
# define persistence
#
# the time series has a jump on the 31/12/2000
# and a jump at every 14th and 21st of Jan-Nov of each year
# for a total of 111 jumps
#
# elimination of jumps and pairing of Obs with Persistence
#

index =  numeric()
dates = inp[1:6337,1]
for(i in seq(1,6336)){
	if((dates[i]+6)!=dates[i+1]){
		YYYYMMDDprev = floor(dates[i]/100)
		YYYYMMDDafter = floor(dates[i+1]/100)
		hhprev = dates[i] - YYYYMMDDprev*100
		hhafter = dates[i+1] - YYYYMMDDafter*100
		if((hhprev != 18)|(hhafter != 0)|((YYYYMMDDprev+1)!=YYYYMMDDafter)){
			DDprev = YYYYMMDDprev - floor(YYYYMMDDprev/100)*100
			MMprev = floor(YYYYMMDDprev/100)-floor(YYYYMMDDprev/10000)*100
			YYYYprev = floor(YYYYMMDDprev/10000)
			DDafter = YYYYMMDDafter - floor(YYYYMMDDafter/100)*100
			MMafter = floor(YYYYMMDDafter/100)-floor(YYYYMMDDafter/10000)*100
			YYYYafter = floor(YYYYMMDDafter/10000)
		if((MMprev == 1)|(MMprev == 3)|(MMprev == 5)|(MMprev == 7)|(MMprev == 8)|(MMprev == 10)){
			if((hhprev != 18)|(hhafter != 0)|(DDprev != 31)|(DDafter != 1)|((MMprev+1) != MMafter)){
				print(c(i,dates[i],dates[i+1]))
				index = c(index,i)}}
		if((MMprev == 4)|(MMprev == 6)|(MMprev == 9)|(MMprev == 11)){
			if((hhprev != 18)|(hhafter != 0)|(DDprev != 30)|(DDafter != 1)|((MMprev+1) != MMafter)){
				print(c(i,dates[i],dates[i+1]))
				index = c(index,i)}}
		if((MMprev == 2)&(YYYYMMDDprev != 2000)){
			if((hhprev != 18)|(hhafter != 0)|(DDprev != 28)|(DDafter != 1)|((MMprev+1) != MMafter)){
				print(c(i,dates[i],dates[i+1]))
				index = c(index,i)}}
		if((MMprev == 2)&(YYYYMMDDprev == 2000)){
			if((hhprev != 18)|(hhafter != 0)|(DDprev != 29)|(DDafter != 1)|((MMprev+1) != MMafter)){
				print(c(i,dates[i],dates[i+1]))
				index = c(index,i)}}
		if(YYYYprev != YYYYafter){
			if((hhprev != 18)|(hhafter != 0)|(DDprev != 31)|(DDafter != 1)|(MMprev != 12)|(MMafter != 1)|
				((YYYYprev+1) != YYYYafter)){
				print(c(i,dates[i],dates[i+1]))
				index = c(index,i)}}
			}
		}
	}
rm(dates,i,YYYYMMDDprev,YYYYMMDDafter,hhafter,hhprev,DDprev,MMprev,YYYYprev,DDafter,MMafter,YYYYafter)
jumpdates = cbind(inp[index,1],inp[index+1,1])

# index is a vector containing the entries of inp where a date jumps occur

indices.obs.p = c(2,sort(c(index,index+2)),6337)

date.YYYYMMDDhh = c()
For.p = c()
Obs.p = c()
for(i in seq(1,length(indices.obs.p),2)){
	date.YYYYMMDDhh = c(date.YYYYMMDDhh, inp[indices.obs.p[i]:indices.obs.p[i+1],1])
	For.p = c(For.p, inp[(indices.obs.p[i]-1):(indices.obs.p[i+1]-1),2])
	Obs.p = c(Obs.p, inp[indices.obs.p[i]:indices.obs.p[i+1],2])
	}
	
	inp.p = cbind(date.YYYYMMDDhh, Obs.p, For.p)
rm(i,date.YYYYMMDDhh,For.p,Obs.p)

# inp.p has 6225 entries

#
# eliminate missing values
#

	esk.p = data.matrix(inp.p)
	esk.p[esk.p<0] = NA
	esk.p = na.omit(esk.p)

# esk.p has 6225 entries. after eliminating NA, 6113 entries

#
# define Persistence Forecast and Observation
#

Obs.p = esk.p[,2]
For.p = esk.p[,3]

#
# Dither
#

dither.esk.obs<-function(x)
	{
	res <- x - floor(x)
	noise <- runif(length(x))
	noise[x==0] <- noise[x==0] * 0.025
	noise[x==0.05] <- noise[x==0.05] * 0.075 - 0.025
	noise[x==0.2] <- noise[x==0.2] * 0.2 - 0.1
	noise[x==0.4] <- noise[x==0.4] * 0.2 - 0.1
	noise[x==0.6] <- noise[x==0.6] * 0.2 - 0.1
	noise[x==0.8] <- noise[x==0.8] * 0.2 - 0.1
	noise[x==0.1] <- 0
	noise[x==0.3] <- 0
	noise[x==0.5] <- 0
	noise[x==0.7] <- 0
	noise[x==0.9] <- 0
	noise[(x<0.95)&(x>0.05)&(x - round(x,1)!=0)]  <- 0
	noise[x==0.95]<- noise[x==0.95] * 0.1 - 0.05
	noise[x==1]<-noise[x==1] * 0.5
	noise[(x>1)&((x - round(x))==0)] <- noise[(x>1)&(x - round(x)==0)] - 0.5
	noise[(x>1)&((x - round(x))!=0)&(res < 0.95)] <- 
		noise[(x>1)&((x - round(x))!=0)&(res < 0.95)]*0.2 -0.1
	noise[(x>1)&(res > 0.9)] <- noise[(x>1)&(res > 0.9)]*0.1-0.05
	noise[(x>1)&(res == 0.1)] <- 0
	noise[(x>1)&(res == 0.5)] <- 0
	x <- x + noise
	x
	}

	dither.esk.for<-function(x)
	{
	noise = runif(length(x))
	noise[x==0] = noise[x==0] * 0.025
	noise[(x>0)&(x<1)] = noise[(x>0)&(x<1)] * 0.05 - 0.025
	noise[x==1] = noise[x==1] * 0.075 - 0.025
	noise[x>1] = noise[x>1] * 0.1 - 0.05
	noise[(x>1)&((x - round(x,1))!=0)] = 0
	x = x + noise
	x
	}

Obs.d = dither.esk.obs(Obs)
For.d = dither.esk.for(For)
Obs.p.d = dither.esk.obs(Obs.p)
For.p.d = dither.esk.obs(For.p)

pobs=(rank(Obs)-1)/(length(Obs)-1)
pfor=(rank(For)-1)/(length(For)-1)
pobs.d=(rank(Obs.d)-1)/(length(Obs.d)-1)
pfor.d=(rank(For.d)-1)/(length(For.d)-1)
pobs.p.d=(rank(Obs.p.d)-1)/(length(Obs.p.d)-1)
pfor.p.d=(rank(For.p.d)-1)/(length(For.p.d)-1)

---------------------------------------------------------------------------------------

#
# Time series                Figure 1
#

postscript("Esk.fig1.ps",horizontal=FALSE)
	par(mfrow=c(2,1))
	par(pty='m')
	par(oma=c(6,2,6,2))
	plot(Obs.d, xlab="", ylab="", type="h", axes=FALSE)
	axis(1, at=seq(1, length(Obs.d), length=6), labels=seq(1998,2003,1))
	axis(2)
	box()
	mtext("a) Eskdalemuir observations",side=3,line=2,cex=2)
	mtext("precipitation (mm)",side=2,line=3,cex=1.5)
	mtext("time (year)",side=1,line=3,cex=1.5)
	plot(For.d, xlab="", ylab="", axes=FALSE,type="h")
	axis(1, at=seq(1, length(For.d), length=6), labels=seq(1998,2003,1))
	axis(2)
	box()
	mtext("b) Eskdalemuir T+6 forecasts",side=3,line=2,cex=2)
	mtext("precipitation (mm)",side=2,line=3,cex=1.5)
	mtext("time (year)",side=1,line=3,cex=1.5)
dev.off()

#
# Scatter-plot For.d vs Obs.d                   Figure 2
#

postscript("Esk.fig2.ps",horizontal=FALSE)
	par(pty='s')
	par(mfrow=c(1,1))
	plot(Obs.d,For.d, xlab="",ylab="", 
	xlim=c(0,max(Obs.d,For.d)), ylim= c(0,max(Obs.d,For.d)),pch=16)
	mtext("forecast (mm)",side=2,line=3,cex=1.5)
	mtext("observation (mm)",side=1,line=3,cex=1.5)
dev.off()	

#
# Empirical distribution functions of For.d and Obs.d              Figure 3
#

postscript("Esk.fig3.ps",horizontal=FALSE)
	par(pty='s')
	par(mfrow=c(1,1))
	plot(sort(Obs.d),sort(pobs.d),type="l",lwd=3, lty=2,
		xlab="",ylab="")
	lines(sort(For.d),sort(pfor.d),type="l",lwd=3,lty=1)
	mtext("1 - base rate",side=2,line=3,cex=1.5)
	mtext("Precipitation amount (mm)",side=1,line=3,cex=1.5)
	abline(1,0)
dev.off()

#
# Scatter-plot probabilities          Figure 4
#

postscript("Esk.fig4.ps",horizontal=FALSE)
	par(pty='s')
	par(mfrow=c(1,1))
	plot(pobs.d,pfor.d,xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),pch=20)
	mtext("cumulative forecast probability",side=2,line=3,cex=1.5)
	mtext("cumulative observation probability",side=1,line=3,cex=1.5)
dev.off()

------------------------------------------------------------------------------------	

#
# Associate thresholds with cumulative probabilities
# 
# ths = c(seq(0, 0.9,0.1), seq(0.91, 0.99, 0.01), seq(0.991, 0.999, 0.001),1)
#

ths<-(seq(1,length(Obs.d),50)-1)/(length(Obs.d)-1)
ths = ths[quantile(Obs.d,ths)>0.025]

# qqplots to compare values of different samples:
# Obs.d has the max value of accumulated precipitation 
# For.p.d has the max quantile value (obs 
# values but sample size is smaller) 
# 
# cp.th.For.d = cbind(ths,quantile(For.d,ths))
# cp.th.Obs.d = cbind(ths,quantile(Obs.d,ths))
# cp.th.For.p.d = cbind(ths,quantile(For.p.d,ths))
# cp.th.Obs.p.d = cbind(ths,quantile(Obs.p.d,ths))


# column 1: cumulative probability
# column 2: threshold corresponding to the comulative probability

#
# Evaluation of the entries of the contingency table (matrix ct.cp)
#
# row 1: threshold (cumulative probability)
# row 2-3-4-5: H,FA,M,CR or a,b,c,d (probabilities).
#

CT.cp<-function(ths,F, O)
	{	
	ct.cp = c()
 	For.cp = (rank(F)-1)/(length(F)-1)
	Obs.cp = (rank(O)-1)/(length(O)-1)
	for(t in ths) 
		{
		a = length(For.cp[(For.cp>t)&(Obs.cp>t)])/length(F)
		b = length(For.cp[(For.cp>t)&(Obs.cp<=t)])/length(F)
		c = length(For.cp[(For.cp<=t)&(Obs.cp>t)])/length(F)
		d = length(For.cp[(For.cp<=t)&(Obs.cp<=t)])/length(F)
		ct.cp = c(ct.cp, t, a,b,c,d)
		}
	ct.cp = matrix(ct.cp, nrow = 5)
	ct.cp
	}

# ct.cp.FO.T6 = CT.cp(ths,For, Obs)
# ct.cp.FOpers.T6 =  CT.cp(ths,For.p, Obs.p)

ct.cp.FO.d.T6 = CT.cp(ths,For.d, Obs.d) 
ct.cp.FOpers.d.T6 =  CT.cp(ths,For.p.d, Obs.p.d)

ct = ct.cp.FO.d.T6 		# = CT.cp(For.d, Obs.d)
ctp = ct.cp.FOpers.d.T6 	# = CT.cp(For.d.p, Obs.d.p)

#
# Evaluation of contingency table scores
#

br = ct[2,]+ct[4,]
H = ct[2,]/(ct[2,]+ct[4,])
F = ct[3,]/(ct[3,]+ct[5,]) 
TS = ct[2,]/(ct[2,]+ct[3,]+ct[4,])
PC = ct[5,]+ct[2,]
arand = (ct[2,]+ct[4,])*(ct[2,]+ct[3,]) 
drand = (ct[5,]+ct[4,])*(ct[5,]+ct[3,]) 
ETS = (ct[2,]-arand)/(ct[2,]+ct[3,]+ct[4,]-arand)
HSS = (PC-arand-drand)/(1-arand-drand)
PSS = (PC-arand-drand)/(1-(ct[5,]+ct[4,])^2-(ct[2,]+ct[4,])^2)
OR = (ct[2,]*ct[5,])/(ct[3,]*ct[4,])
YQ = (OR-1)/(OR+1)
chi = 2 - log(ct[5,])/log(ct[3,] + ct[5,])
chibar = 2*log(ct[2,]+ct[4,])/log(ct[2,]) - 1

brp = ctp[2,]+ctp[4,]
Hp = ctp[2,]/(ctp[2,]+ctp[4,])
Fp = ctp[3,]/(ctp[3,]+ctp[5,]) 
TSp = ctp[2,]/(ctp[2,]+ctp[3,]+ctp[4,])
PCp = ctp[5,]+ctp[2,]
arandp = (ctp[2,]+ctp[4,])*(ctp[2,]+ctp[3,]) 
drandp = (ctp[5,]+ctp[4,])*(ctp[5,]+ctp[3,]) 
ETSp = (ctp[2,]-arandp)/(ctp[2,]+ctp[3,]+ctp[4,]-arandp)
HSSp = (PCp-arandp-drandp)/(1-arandp-drandp)
PSSp = (PCp-arandp-drandp)/(1-(ctp[5,]+ctp[4,])^2-(ctp[2,]+ctp[4,])^2)
ORp = (ctp[2,]*ctp[5,])/(ctp[3,]*ctp[4,])
ORp[ctp[2,]==0]<- NA
YQp = (ORp-1)/(ORp+1)
chip = 2 - log(ctp[5,])/log(ctp[3,] + ctp[5,])
chibarp = 2*log(ctp[2,]+ctp[4,])/log(ctp[2,]) - 1
chibarp[ctp[2,]==0]<- NA

---------------------------------------------------------------------------------------------------

#
# Plot of joint probabilities vs base rate in log scale             Figure 5
#

postscript("Esk.fig5.ps",horizontal=FALSE)
	par(pty='s')
	par(mfrow=c(1,1))
	plot(1,1,type="n",log="xy",xlab="",ylab="",xlim=c(min(1-ct[1,]),max(1-ct[1,])),ylim=c(min(ct[2:3,]),1))
	lines(1-ct[1,],ct[2,],lty=2,lwd=3)
	lines(1-ct[1,],ct[3,],lty=3,lwd=3)
	lines(1-ct[1,],ct[5,],lty=4,lwd=3)
	abline(0,1,lwd=1,lty=1)
	mtext("joint probabilities",side=2,line=3,cex=1.5)
	mtext("base rate",side=1,line=3,cex=1.5)
#	legend(min(1-ct[1,]),0.5,lty=c(2,3,4),legend=c("a","b=c  ","d"),lwd=3,cex=1.5)
dev.off()

#
# Plots of hit rate, false alarm rate and ROC curves                Figure 6
#


postscript("Esk.fig6.ps",horizontal=FALSE)
	par(pty='s')
	par(mfrow=c(2,2))
	par(oma=c(6,2,6,2))

# hit rate

	plot(quantile(Obs.d,ths),H,type="l",lwd=3,
		xlab="threshold (mm)",ylab="hit rate (H)",ylim=c(0,0.8))
	lines(quantile(Obs.d,ths),Hp,lty=2,lwd=3)
	lines(quantile(Obs.d,ths),1-ths,lty=3,lwd=3)
#	legend(x=8,y=0.8,legend=c("Met Office T+6", "Persistence","Random (H=p)"),lty=c(1,2,3),lwd=3,cex=0.8)
	axis(3,at=quantile(Obs.d,ths),labels=round(ths,3))
	mtext("1 - base rate",side=3,line=3,cex=0.8)
	mtext("a)                                      ",side=3,line=2,cex=2)

# false alarm rate

	plot(quantile(Obs.d,ths),F,type="l",lwd=3,
		xlab="threshold (mm)",ylab="false alarm rate (F)",ylim=c(0,0.54))
	lines(quantile(Obs.d,ths),Fp,lty=2,lwd=3)
	lines(quantile(Obs.d,ths),1-ths,lty=3,lwd=3)
#	legend(x=8,y=0.5,legend=c("Met Office T+6", "Persistence","Random (F=p)"),lty=c(1,2,3),lwd=3,cex=0.8)
	axis(3,at=quantile(Obs.d,ths),labels=round(ths,3))
	mtext("1 - base rate",side=3,line=3,cex=0.8)
	mtext("b)                                      ",side=3,line=2,cex=2)

# ROC curve

	plot(F,H,xlab="F",ylab="H",ylim=c(0,1),xlim=c(0,1),type="l",lwd=3)
	lines(Fp,Hp,lty=2,lwd=3)
	abline(0,1,lty=3,lwd=3)
	title("ROC curve")
#	legend(x=0.3,y=0.3,legend=c("Met Office T+6", "Persistence","Random (F=H)"),lty=c(1,2,3),lwd=3,cex=0.8)
	mtext("c)                                      ",side=3,line=2,cex=2)

# ROC curve on log axes

	plot(log(F/(1-F)),log(H/(1-H)),xlab="log F/(1-F)",ylab="log H/(1-H)",
		ylim=c(-7,2),xlim=c(-7,2),type="l",lwd=3)
	lines(log(Fp/(1-Fp)),log(Hp/(1-Hp)),lty=2,lwd=3)
	abline(0,1,lty=3,lwd=3)
	title("logistic ROC curve")
#	legend(x=-4,y=-4,legend=c("Met Office T+6", "Persistence","Random (F=H)"),lty=c(1,2,3),lwd=3,cex=0.8)
	mtext("d)                                     ",side=3,line=2,cex=2)

dev.off()

#
# Plots of PC, PSS, ETS, HSS vs threshold                Figure 7
#

postscript("Esk.fig7.ps",horizontal=FALSE)
	par(pty='m')
	par(mfrow=c(2,2))
	par(oma=c(6,2,6,2))

# PC
	plot(quantile(Obs.d,ths),PC,type="l",lwd=3,
		xlab="threshold (mm)",ylab="Proportion Correct (PC)",ylim=c(0.5,1.0))
	lines(quantile(Obs.d,ths),PCp,lty=2,lwd=3)
	lines(quantile(Obs.d,ths),arand+drand,lty=3,lwd=3)
	axis(3,at=quantile(Obs.d,ths),labels=round(ths,3))
	mtext("1 - base rate",side=3,line=3,cex=0.8)
	mtext("a)                                      ",side=3,line=2,cex=2)

# PSS 
	plot(quantile(Obs.d,ths),PSS,type="l",lwd=3,
		xlab="threshold (mm)",ylab="Peirce Skill Score (PSS)",ylim=c(0,max(PSS)))
	lines(quantile(Obs.d,ths),PSSp,lty=2,lwd=3)
	lines(quantile(Obs.d,ths),rep(0,length(ths)),lty=3,lwd=3)
	axis(3,at=quantile(Obs.d,ths),labels=round(ths,3))
	mtext("1 - base rate",side=3,line=3,cex=0.8)
	mtext("b)                                      ",side=3,line=2,cex=2)

# ETS
	plot(quantile(Obs.d,ths),ETS,type="l",lwd=3,
		xlab="threshold (mm)",ylab="Equitable Threat Score (ETS)",ylim=c(0,max(ETS)))
	lines(quantile(Obs.d,ths),ETSp,lty=2,lwd=3)
	lines(quantile(Obs.d,ths),rep(0,length(ths)),lty=3,lwd=3)
	axis(3,at=quantile(Obs.d,ths),labels=round(ths,3))
	mtext("1 - base rate",side=3,line=3,cex=0.8)
	mtext("c)                                      ",side=3,line=2,cex=2)

# HSS
	plot(quantile(Obs.d,ths),HSS,type="l",lwd=3,
		xlab="threshold (mm)",ylab="Heidke Skill Score (HSS)",ylim=c(0,max(HSS)))
	lines(quantile(Obs.d,ths),HSSp,lty=2,lwd=3)
	lines(quantile(Obs.d,ths),rep(0,length(ths)),lty=3,lwd=3)
	axis(3,at=quantile(Obs.d,ths),labels=round(ths,3))
	mtext("1 - base rate",side=3,line=3,cex=0.8)
	mtext("d)                                      ",side=3,line=2,cex=2)

dev.off()

#
# Plots of log odds ratio vs threshold                Figure 8
#

postscript("Esk.fig8.ps",horizontal=FALSE)
	par(pty='m')
	par(mfrow=c(1,1))
	par(oma=c(6,2,6,2))
	plot(quantile(Obs.d,ths),log(OR),type="l",lwd=3,
		xlab="",ylab="",ylim=c(0,max(log(OR))))
	lines(quantile(Obs.d,ths),log(ORp),lty=2,lwd=3)
	lines(quantile(Obs.d,ths),rep(0,length(ths)),lty=3,lwd=3)
	axis(3,at=quantile(Obs.d,ths),labels=round(ths,3))
	mtext("threshold (mm)",side=1,line=3,cex=1.5)
	mtext("log odds ratio",side=2,line=3,cex=1.5)
	mtext("1 - base rate",side=3,line=3,cex=1.5)
dev.off()

#
# Plots of EDS vs threshold                Figure 9
#

varchibar = ((H*(1-H))/(length(Obs.d)*br))*((2*log(br))/(H*(log(br*H))^2))^2
varchibar.p = ((Hp*(1-Hp))/(length(Obs.p.d)*brp))*((2*log(brp))/(H*(log(brp*Hp))^2))^2
varchibarrand = (1-br)/(4*length(Obs.d)*(br*log(br))^2)
varchibarrandp = (1-brp)/(4*length(Obs.p.d)*(brp*log(brp))^2)
					   
postscript("Esk.fig9.ps",horizontal=FALSE)
	par(pty='m')
	par(mfrow=c(1,1))
	par(oma=c(6,2,6,2))
	plot(quantile(Obs.d,ths),chibar,type="l",lwd=3,
		xlab="",ylab="",ylim=c(0,max(chibar+1.96*sqrt(varchibar))))
	lines(quantile(Obs.d,ths),(chibar+1.96*sqrt(varchibar)))
	lines(quantile(Obs.d,ths),(chibar-1.96*sqrt(varchibar)))
	lines(quantile(Obs.d,ths),chibarp,lty=2,lwd=3)
	lines(quantile(Obs.d,ths),(chibarp+1.96*sqrt(varchibarp)),lty=2)
	lines(quantile(Obs.d,ths),(chibarp-1.96*sqrt(varchibarp)),lty=2)
	lines(quantile(Obs.d,ths),rep(0,length(ths)),lty=3,lwd=3)
#	lines(quantile(Obs.d,ths),1.96*sqrt(varchibarrand),lty=3)
#	lines(quantile(Obs.d,ths),-1.96*sqrt(varchibarrand),lty=3)	
#	lines(quantile(Obs.d,ths),1.96*sqrt(varchibarrandp),lty=3)
#	lines(quantile(Obs.d,ths),-1.96*sqrt(varchibarrandp),lty=3)	
	axis(3,at=quantile(Obs.d,ths),labels=round(ths,3))
	mtext("threshold (mm)",side=1,line=3,cex=1.5)
	mtext("Extreme Dependency Score (EDS)",side=2,line=3,cex=1.5)
	mtext("1 - base rate",side=3,line=3,cex=1.5)
dev.off()

postscript("Esk.fig9b.ps",horizontal=FALSE)
	par(pty='m')
	par(mfrow=c(1,1))
	par(oma=c(6,2,6,2))
	plot(quantile(Obs.d,ths),chibar,type="l",lwd=3,
		xlab="",ylab="",ylim=c(0,max(chibar+1.96*sqrt(varchibar))))
	lines(quantile(Obs.d,ths),(chibar+1.96*sqrt(varchibar)))
	lines(quantile(Obs.d,ths),(chibar-1.96*sqrt(varchibar)))
	lines(quantile(Obs.d,ths),chibarp,lty=2,lwd=3)
	lines(quantile(Obs.d,ths),(chibarp+1.96*sqrt(varchibarp)),lty=2)
	lines(quantile(Obs.d,ths),(chibarp-1.96*sqrt(varchibarp)),lty=2)
	lines(quantile(Obs.d,ths),rep(0,length(ths)),lty=3,lwd=3)
	lines(quantile(Obs.d,ths),rep(0,length(ths)),lty=3,lwd=3)
	lines(quantile(Obs.d,ths),1.96*sqrt(varchibarrand),lty=3,col=2)
	lines(quantile(Obs.d,ths),1.96*sqrt(varchibarrandp),lty=3,col=4)	
	axis(3,at=quantile(Obs.d,ths),labels=round(ths,3))
	mtext("threshold (mm)",side=1,line=3,cex=1.5)
	mtext("Extreme Dependency Score (EDS)",side=2,line=3,cex=1.5)
	mtext("1 - base rate",side=3,line=3,cex=1.5)
dev.off()
