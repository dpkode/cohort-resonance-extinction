# Plot SALMOD simulation results for various climate scenarios ####

# written by Chris Mosser
# modified by Patrick Kilduff 
# in November 2012

#Click on the console window
#Click file, change dir
#Navigate to where this file is stored
#Highlight code and press 'ctl'-'r' to run the code
#########################
#########################

setwd("~/NEP_salmon/chap_3/data/simulation_results/")

#Load climate data into memory
data = read.csv("CW_extract.csv")
datam <- melt(data, id = "YEAR")
ggplot(data = datam, aes(x = YEAR, y = value)) + 
  geom_line() +
  facet_wrap(~variable, nrow = 4)
  
year = data[1:90, 1]
A2_CNRMCM3 = data[1:90, 2]
A2_GFDLCM21 = data[1:90, 3]
A2_MIROC32MED = data[1:90, 4]
A2_MPIECHAM5 = data[1:90, 5]
A2_NCARCCSM3 = data[1:90, 6]
A2_NCARPCM1 = data[1:90, 7]
B1_CNRMCM3 = data[1:90, 8]
B1_GFDLCM21 = data[1:90, 9]
B1_MIROC32MED = data[1:90, 10]
B1_MPIECHAM5 = data[1:90, 11]
B1_NCARCCSM3 = data[1:90, 12]
B1_NCARPCM1 = data[1:90, 13]

#########################
#########################

#Print 4x3 plot of baseline climate data
#tiff(filename="fig8.tiff", width=3000, height=1800, res=300, pointsize=8)
op <- par(mfcol = c(4,3))

plot(year, A2_CNRMCM3, type='h', ylab="Proportion", main='A2 CNRMCM3', xlab='')
plot(year, B1_CNRMCM3, type='h', ylab="Proportion", main='B1 CNRMCM3', xlab='')
plot(year, A2_GFDLCM21, type='h', main='A2 GFDLCM21', xlab='', ylab='Proportion')
plot(year, B1_GFDLCM21, type='h', main='B1 GFDLCM21', xlab='Year', ylab='Proportion')
plot(year, A2_MIROC32MED, type='h', main='A2 MIROC32MED', xlab='', ylab='')
plot(year, B1_MIROC32MED, type='h', main='B1 MIROC32MED', xlab='', ylab='')
plot(year, A2_MPIECHAM5, type='h', main='A2 MPIECHAM5', xlab='', ylab='')
plot(year, B1_MPIECHAM5, type='h', main='B1 MPIECHAM5', xlab='Year', ylab='')
plot(year, A2_NCARCCSM3, type='h', main='A2 NCARCCSM3', xlab='', ylab='')
plot(year, B1_NCARCCSM3, type='h', main='B1 NCARCCSM3', xlab='', ylab='')
plot(year, A2_NCARPCM1, type='h', main='A2 NCARPCM1', xlab='', ylab='')
plot(year, B1_NCARPCM1, type='h', main='B1 NCARPCM1', xlab='Year', ylab='')
par(op)
#dev.off()

############################
############################
############################
############################

#Load simulation results into memory

BAU = read.csv("BAU_Extract.csv")/15000
CW = read.csv("CW_Extract.csv")/15000
Forecast = read.csv("Forecast_Extract.csv")/15000
FE = read.csv("FE_Extract.csv")/15000
ND = read.csv("ND_Extract.csv")/15000
Shade = read.csv("Shade_Extract.csv")/15000
RP = read.csv("RP_Extract.csv")/15000

############################
############################
############################
############################



############################
############################
############################
############################



############################
############################
############################
############################
#Proportional spawner survival graph
#Must load 'data1', 'data2', and 'data3' above to use this code

par(mfcol = c(1,1))
#x1 = (colSums(data1[,2:13])-colSums(data[,2:13]))/colSums(data[,2:13])
#x2 = (colSums(data2[,2:13])-colSums(data[,2:13]))/colSums(data[,2:13])
#x3 = (colSums(data3[,2:13])-colSums(data[,2:13]))/colSums(data[,2:13])
#x4 = (colSums(dataBAU[,2:13])-colSums(data[,2:13]))/colSums(data[,2:13])
#x5 = (colSums(dataForcast[,2:13])-colSums(data[,2:13]))/colSums(data[,2:13])

x1 = (colSums(CW[,2:13])-colSums(BAU[,2:13]))/colSums(BAU[,2:13])
x2 = (colSums(Forecast[,2:13])-colSums(BAU[,2:13]))/colSums(BAU[,2:13])
x3 = (colSums(FE[,2:13])-colSums(BAU[,2:13]))/colSums(BAU[,2:13])
x4 = (colSums(Shade[,2:13])-colSums(BAU[,2:13]))/colSums(BAU[,2:13])
x5 = (colSums(RP[,2:13])-colSums(BAU[,2:13]))/colSums(BAU[,2:13])
x6 = (colSums(ND[,2:13])-colSums(BAU[,2:13]))/colSums(BAU[,2:13])

dat = rbind(x1, x2, x3, x4, x5, x6)
labs = c('A2 CNRMCM3', 'A2 GFDLCM21', 'A2 MIROC32MED', 'A2 MPIECHAM5', 'A2 NCARCCSM3', 'A2 NCARPCM1', 'B1 CNRMCM3', 'B1 GFDLCM21', 'B1 MIROC32MED', 'B1 MPIECHAM5', 'B1 NCARCCSM3', 'B1 NCARPCM1')
barplot(as.matrix(dat), beside=TRUE, cex.names=0.75, ylim=c(-0.2, 0.8), space=c(0,3), names.arg=labs, col=c("red", "green", "blue", "yellow", "purple", "black"), xlab="Climate Scenario", ylab="Proportional Difference in Spawner Survival")
abline(h=0)
legend("topright", c("CW", "Forecast", "FE", "Shade", "RP", "ND"), bty='n', fill=c("red", "green", "blue", "yellow", "purple", "black"))



############################
############################
############################
############################
tiff(filename="fig8a.tiff", width=3000, height=1800, res=300, pointsize=8)
#This is a modified version of the 4x3 graph where the data
# is truncated at the extinction year

#Highlight from here-down to 'STOP'

#Calculate year to extinction
extinct = function(dats, crit, cons)
{
	#dats = data set
	#crit = critical pop size
	#cons = consectutive years of reaching critical pop size

	# makes abundance < crit = 0
	dats = round(dats*15000)
	dat_crit = (dats > crit)*dats 
	

	ret = 0
	for (i in (cons-1):87)
	{
		if((dat_crit[i] < 1) && (dat_crit[i+1] < 1) && (dat_crit[i+2] < 1) && (dat_crit[i+3] < 1))
		{
			ret = i-1
			break
		}
		else { ret = 90 }
	} 
	ret

}

extincts = function(datas, crit, cons)
{
	cols = ncol(datas)
	ext = matrix(0, 1,cols)
	for (j in 1:cols)
	{
		ext[j] = extinct(datas[, j], crit, cons)
	}
	ext
}

criteria = 0

base = extincts(data[,2:13], criteria, 4)
man1 = extincts(data1[,2:13], criteria, 4)
man2 = extincts(data2[,2:13], criteria, 4)
man3 = extincts(data3[,2:13], criteria, 4)

ext = rbind(base, man1, man2, man3)
ext

#4x3 plot of climate data
par(mfcol = c(4,3))

plot(year[1:base[1]], A2_CNRMCM3[1:base[1]], type='h', ylab="Proportion", main='A2 CNRMCM3', xlab='', xlim=c(2010, 2100))
plot(year[1:base[7]], B1_CNRMCM3[1:base[7]], type='h', ylab="Proportion", main='B1 CNRMCM3', xlab='', xlim=c(2010, 2100))
plot(year[1:base[2]], A2_GFDLCM21[1:base[2]], type='h', main='A2 GFDLCM21', xlab='', ylab='Proportion', xlim=c(2010, 2100))
plot(year[1:base[8]], B1_GFDLCM21[1:base[8]], type='h', main='B1 GFDLCM21', xlab='Year', ylab='Proportion', xlim=c(2010, 2100))
plot(year[1:base[3]], A2_MIROC32MED[1:base[3]], type='h', main='A2 MIROC32MED', xlab='', ylab='', xlim=c(2010, 2100))
plot(year[1:base[9]], B1_MIROC32MED[1:base[9]], type='h', main='B1 MIROC32MED', xlab='', ylab='', xlim=c(2010, 2100))
plot(year[1:base[4]], A2_MPIECHAM5[1:base[4]], type='h', main='A2 MPIECHAM5', xlab='', ylab='', xlim=c(2010, 2100))
plot(year[1:base[10]], B1_MPIECHAM5[1:base[10]], type='h', main='B1 MPIECHAM5', xlab='Year', ylab='', xlim=c(2010, 2100))
plot(year[1:base[5]], A2_NCARCCSM3[1:base[5]], type='h', main='A2 NCARCCSM3', xlab='', ylab='', xlim=c(2010, 2100))
plot(year[1:base[11]], B1_NCARCCSM3[1:base[11]], type='h', main='B1 NCARCCSM3', xlab='', ylab='', xlim=c(2010, 2100))
plot(year[1:base[6]], A2_NCARPCM1[1:base[6]], type='h', main='A2 NCARPCM1', xlab='', ylab='', xlim=c(2010, 2100))
plot(year[1:base[12]], B1_NCARPCM1[1:base[12]], type='h', main='B1 NCARPCM1', xlab='Year', ylab='', xlim=c(2010, 2100))

# STOP
# End of 4x3 graph
dev.off()

############################
############################
############################
############################

tiff(filename="fig9_10a.tiff", width=3000, height=1500, res=300, pointsize=8)
par(mfcol = c(1,1))
x1 = (colSums(data1[,2:13])-colSums(data[,2:13]))/colSums(data[,2:13])
x2 = (colSums(data2[,2:13])-colSums(data[,2:13]))/colSums(data[,2:13])
x3 = (colSums(data3[,2:13])-colSums(data[,2:13]))/colSums(data[,2:13])
dsum = colSums(data[,2:13])
d1sum = colSums(data1[,2:13])
d2sum = colSums(data2[,2:13])
d3sum = colSums(data3[,2:13])
dat = rbind(x1, x2, x3)
#y=c(x1[1], x2[1], x3[1])

labs = c('A2 CNRMCM3', 'A2 GFDLCM21', 'A2 MIROC32MED', 'A2 MPIECHAM5', 'A2 NCARCCSM3', 'A2 NCARPCM1', 'B1 CNRMCM3', 'B1 GFDLCM21', 'B1 MIROC32MED', 'B1 MPIECHAM5', 'B1 NCARCCSM3', 'B1 NCARPCM1')
#labs = c(1:12)



criteria = 0

base = extincts(data[,2:13], criteria, 4)
man1 = extincts(data1[,2:13], criteria, 4)
man2 = extincts(data2[,2:13], criteria, 4)
man3 = extincts(data3[,2:13], criteria, 4)
#BAU = extincts(dataBAU[,2:13], criteria, 4)
#Forcast = extincts(dataForcast[,2:13], criteria, 4)

ext = rbind(base, man1, man2, man3)


barplot(as.matrix(ext), col=c("gray22", "gray47", "gray72", "gray97"), beside=TRUE, cex.names=0.75, ylim=c(0, 100), space=c(0,3), names.arg=labs, xlab="Climate Scenario", ylab="Years Until Extinction Threshold Met")
#abline(h=90)
abline(h=0)
legend(bg='white', "topleft", c("Climate", "Management 1", "Management 2", "Management 3"), bty='n', fill=c("gray22", "gray47", "gray72", "gray97"))

dev.off()
############################
############################
############################
############################

#Fig 4a

tiff(filename="fig4.tiff", width=3000, height=1500, res=300, pointsize=8)
cx = 1.5
par(mfcol = c(2,1))#use for the next 2 graphs
b0 = 115.08
b1 = -5.42

T = seq(15, 25, 0.05)

Pr = 1-1/(1+exp(-(b0+(b1*T))))
Pr

plot(T, Pr, xlab = "Temperature (C)", lwd=3, cex.lab=cx, cex.axis=cx, ylab = "Weekly Mortality Rate", type='l', xlim=c(18, 24))
mtext(side=1, line=3, adj=0, "A", cex=cx)
############################
############################
############################
############################

#Fig 4b including WEAP historic

obs = c(139, 1716, 5616, 207, 309, 122, 316, 489)
mod1 = c(682, 1349, 4953, 543, 1094, 623, 39, 132)
mod2 = c(5909, 1928, 2034, 285, 221)
yrs = 2001:2008


#spawners = rbind(obs, mod1, mod2)

plot(yrs, obs, lwd=0.75, pch=1, xlab="Year", cex.axis=cx, ylab="Pre-spawn Mortality (Females)", ylim=c(0, 6000), cex.lab=cx)
points(yrs, mod1, pch=2, lwd=0.75)
points(yrs[1:5], mod2, pch=7, lwd=0.75)
#legend("topright", pch=c(1,2,7), c("Simulated mortality using \n simulated temperature and flow", "Simulated mortality using \n observed temperature and flow", "Observed mortality data \n"), bty='n', cex=1.5, y.intersp=1.5 )
mtext(side=1, line=3, adj=0, "B", cex=cx)
dev.off()
############################
############################
############################
############################
library(survey)
mod = lm(obs~mod1)
cof = c(-242.989, 1.1533)
diff = c(0, 1)
se = c(233.5410, 0.1233)
tstat = (cof[1]-diff[1])/se[1]

pt()

x=regTermTest(mod, "mod1", method="Wald")
summary(mod)