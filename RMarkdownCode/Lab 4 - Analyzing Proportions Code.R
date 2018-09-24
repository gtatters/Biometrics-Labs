# Lab 4 Code
#  Lab Objectives:
#  ̈ Use SPSS to conduct binomial testing
#  ̈ Use SPSS to compute a complete binomial distribution
#  ̈ Use Excel to calculate confidence intervals for proportions using the Agresti-Coull method
#  ̈ Practice using SPSS and Excel

setwd("~/Documents/Teaching/Biometrics - BIOL 3P96/Labs/Lab 4")

d<-read.csv("Earturn.csv")

ears<-xtabs(~d$Ear)
success.fail<-xtabs(~d$Ear)

binom.test(success.fail, p = 0.5,
           alternative = c("two.sided"),
           conf.level = 0.95)
# alternatively, we can denote x and n:
  
binom.test(x=6, n=19+6, p = 0.5,
           alternative = c("two.sided"),
           conf.level = 0.95)


# Spermatogenesis genes
# A study of genes involved in spermatogenesis was carried out to test the hypothesis 
# that spermatogenesis genes should occur disproportionately often on the X chromosome,
# which contains 6.1% of the total genes. The researchers found that 40% of the
# analyzed spermatogenesis genes were on the X chromosome.

d<-read.csv("Spermatogenesis.csv")
counts<-xtabs(~d$chromosome)
barplot(counts, beside=T, xlab="Chromosome", ylab="Number of genes")

counts<-table(d$onX,d$chromosome )
counts

barplot(counts, beside=T, xlab="Chromosome", ylab="Number of genes", 
        legend=rownames(counts), ylim=c(0,12))

x=10  
n=sum(counts)
p=0.061
binom.test(x, n, p, alternative = "two.sided")$p.value

# How many specifically on Y?  Make use of the sum() command and a boolean expression
sum(d$chromosome=="Y")
x=3
n=sum(counts)
p=0.012
# if I set p=0.38, I get the answer Lori has, but I don't know how you arrive at p=0.38
binom.test(x, n, p=0.38, alternative = "two.sided")$p.value

3/(25-10)

# Lost wallets
# In a study in Scotland, researchers left a total of 240 wallets around Edinburgh,
# as though the wallets were lost. Each contained contact information including an
# address. Data on the number of wallets that were returned is found in the 
# file Wallets.csv.

d<-read.csv("Wallets.csv")
str(d)
counts<-table(d)
counts
p=0.5
result<-binom.test(counts, p, alternative = "two.sided")
result$statistic
result$p.value





#  Groundhog day stats
d<-read.csv("GroundhogDay.csv")
str(d)
counts<-table(d$Result)
# According to http://www.stormfax.com/ghogday.htm, the groundhog's accuracy is 
# 39%, even though the actual spring data were not provided.
sum(counts)
# From the table above, there are 10 records with no data and 130 records in total 
# thus n = 130-10 = 120
n=130-10
# Success?  39% is the success rate, thus:
x=.39*n
# note this yields a fraction, not a whole number, so we are not working with the 
# precise data, and the binomial test will assume we need whole numbers
x<-round(0.39*n)
x

binom.test(x, n, p=0.5, alternative="two.sided")

# Punxsutawney Phil's prediction rate fall below that of random chance



# Wine example?  I need SPSS to follow this example.  Maybe it means using
# rbinom?

ranbinom<-rbinom(100000, 15, 0.5)
h <- hist(ranbinom,  plot=FALSE)
cuts <- cut(h$breaks, c(-Inf,4,Inf))
hist(ranbinom, probability=T, col=cuts, main="Probability of >= 4 people in red",
     xlab="Number of people who prefer expensive wine")
# lines(density(ranbinom, bw=1), col="blue")

f<-hist(ranbinom, breaks=20, freq=F)
f$breaks
f$counts
lines(f)
curve(dbinom(xx, 15, 0.5))
plot(density(xx))

pbinom(4,15,0.5, lower.tail=T)

# Heterosygous crosses
p=0.25
n=20

rbinom(100000,20,prob=p)


# Toast example
x<-9821
n<-9821+6101
pp<-(x+2)/(n+4)
CR<-1.96*sqrt(pp*(1-pp)/(n+4))
x/n - CR
x/n + CR

counts<-rbinom(100000, n, p=x/n)/n
100000*0.025
100000*0.975
sort(counts)[2500]
sort(counts)[97500]



# Country data

d<-read.csv("Countries.csv")
str(d)
I<-d$Immunization_DPT
h<-hist(I, breaks=seq(0,100,10))
h$breaks

summary(I)
mean(I, na.rm=TRUE)
sd(I, na.rm=TRUE)



# trouble with binom test p values for Y chromosome
xsuccesses <- 0:25
probx <- dbinom(xsuccesses, size = 25, prob = 0.012)
data.frame(xsuccesses, probx)

xsuccesses
probx

# Use these probabilities to calculate the P-value corresponding to an observed 3 spermatogenesis genes on the Y chromosome. 
# Remember to multiply the probability of 10 or more successes by 2 for the two-tailed test result.

2 * sum(probx[xsuccesses >= 3])

(1-sum(probx[!xsuccesses==3]))*2

binom.test(3, 10, p=0.2, alternative="two.sided")$p.value
binom.test(3, 10, p=0.2, alternative="greater")$p.value*2

