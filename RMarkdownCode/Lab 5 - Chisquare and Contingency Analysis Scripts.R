# Lab 5 scripts
setwd("~/Dropbox/Biometrics Labs/Lab 5")

d<-read.csv("days of the week.csv")
str(d)
chisq.test(xtabs(~d$day))



d<-read.csv("CCSVI.csv")
str(d)
plot(d)
xtabs(~d$MS + d$CCSVI)

chisq.test(xtabs(~d$MS+d$CCSVI))


x<-rpois(10000, lambda=0.604)
plot(density(x))
     hist(x)


# truffles
d<-read.csv("truffles.csv")
str(d)
x<-d$Frequency
x

# pois between 1 & 2#:
ppois(0, 0.604)
ppois(1, 0.604)-ppois(0, 0.604)
ppois(2, 0.604)-ppois(1, 0.604)
ppois(3, 0.604)-ppois(2, 0.604)
ppois(4, 0.604)-ppois(3, 0.604)

p0<-ppois(0, 0.604)
p1_4<-ppois(1:4, lambda=0.604) - ppois(0:3, lambda=0.604)
p<-c(p0, p1_4)
p
# alternatively, use the dpois function to obtain identical result!
dpois(0:4, lambda=0.604)



freq<-p*sum(d$Frequency)
freq # one of these is too low (i.e. <5)
newfreq<-freq[1:4]
newfreq[4]<-newfreq[4]+freq[5]
newfreq

newx<-c(203, 39, 18, 13+15)
newx

chisq.test(newx, p=newfreq, rescale.p=T)


