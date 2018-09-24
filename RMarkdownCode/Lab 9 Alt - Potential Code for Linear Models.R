# ------------------------------------------------------------

# What I would like to propose are 1 or 2 additional/substitute/future labs that take
# the ANOVA/lm lab 8 and extend them, in order to demonstrate that most statistics
# could be replaced with glms!  Future use in R will be on the higher level stats.
# I think a good place to start would be to demonstrate that a t-test can be 
# performed using the linear model, lm() command
# otherwise, there are good examples from Whitlock and Schluter on multivariate 
# lms that introduce the broader applicability of lm and glm approaches. 
# Drawing heavily from Chapter 17 and 18 + other sources (TBD)


setwd("~/Dropbox/Biometrics Labs/Lab 9 Alt - Linear and Multivariable Models")


# Figure 17.5-4 (left). <a href="../wp-content/data/chapter17/chap17f5_4BlueTitCapColor">Blue tit cap color
# Create a residual plot to check assumptions. The data are cap color of offspring 
# and parents in a sample of the blue tit.
# Read and inspect the data.

capcolour <- read.csv("capcolour.csv")

# Basic scatter plot with regression line.

plot(offspringColorScore ~ midparentColorScore, data = capColor)
capColorRegression <- lm(offspringColorScore ~ midparentColorScore, data = capColor)
abline(capColorRegression)

# Residual plot. Residuals are calculated from the model object created from the previous lm command, using the residuals command.

plot(residuals(capColorRegression) ~ midparentColorScore, data = capColor)
abline(c(0,0))

# Commands for a fancier residual plot using more options are here.

plot(residuals(capColorRegression) ~ midparentColorScore, data = capColor, 
     pch = 16, col = "firebrick", las = 1, cex = 1.5, bty = "l", 
     xlab = "Midparent color score", ylab = "Residuals")
abline(c(0,0))

# ------------------------------------------------------------


# Figure 17.6-3. <a href="../wp-content/data/chapter17/chap17f6_3IrisPollination.csv">Iris pollination
# Use a square root transformation to improve the fit to assumptions of linear 
# regression. The data are number of pollen grains received and floral tube length of
# an iris species.
# Read and inspect data.

iris <- read.csv("iris.csv")
str(iris)

# Scatter plot with regression line.

plot(grainsDeposited ~ tubeLengthMm, data = iris)
irisRegression <- lm(grainsDeposited ~ tubeLengthMm, data = iris)
abline(irisRegression)

# Residual plot.

plot(residuals(irisRegression) ~ tubeLengthMm, data = iris)
abline(c(0,0))

# Commands for a residual plot with more options are here.

plot(residuals(irisRegression) ~ tubeLengthMm, data = iris, pch = 16,  
     col = "firebrick", las = 1, cex = 1.5, bty = "l",
     xlab = "Floral tube length (mm)", ylab = "Residuals")
abline(c(0,0))

# Square root transformation.

iris$sqRootGrains <- sqrt(iris$grainsDeposited + 1/2)

# Scatter plot using transformed data, with new regression line added.

plot(sqRootGrains ~ tubeLengthMm, data = iris)
irisRegressionSqrt <- lm(sqRootGrains ~ tubeLengthMm, data = iris)
abline(irisRegressionSqrt)

# Residual plot based on the transformed data.

plot(residuals(irisRegressionSqrt) ~ tubeLengthMm, data = iris)
abline(c(0,0))

# Instructions for a residual plot with more options are here.

plot(residuals(irisRegressionSqrt) ~ tubeLengthMm, data = iris,  
     pch= 16, col = "firebrick", las = 1, cex=1.5, bty = "l", 
     xlab = "Floral tube length (mm)", ylab = "Residuals")
abline(c(0,0))

# ------------------------------------------------------------

# Figure 17.8-1. <a href="../wp-content/data/chapter17/chap17f8_1IronAndPhytoplanktonGrowth.csv">Iron and phytoplankton growth
# Fit a nonlinear regression curve having an asymptote (Michaelis-Menten curve). 
# The data are the relationship between population growth rate of a phytoplankton
# in culture and the concentration of iron in the medium.
# Read and examine data.

phytoplankton <- read.csv("phytoplankton.csv")

# Scatter plot.

plot(phytoGrowthRate ~ ironConcentration, data = phytoplankton)

# Instructions for a scatter plot with more options are here.

plot( phytoGrowthRate ~ ironConcentration, data = phytoplankton, pch = 16, 
      col = "firebrick", las = 1, cex = 1.5, bty = "l", 
      ylab = "Growth rate (no./day)",
      xlab = expression(paste("Iron concentration (", mu, "mol)")) )

#   	
# Fit a Michaelis-Menten curve to the phytoplankton data using the nls
# (nonlinear least squares). To fit the curve, provide a formula that also 
# includes symbols (letters) for the parameters to be estimated. In the following
# function, we use "a" and "b" to indicate the parameters of the Michaelis-Menten
# curve we want to estimate. The function includes an argument where we must provide 
# an initial guess of parameter values. The value of the initial guess is not so 
# important -- here we choose a=1 and b=1 as initial guesses. The first function
# below carried out the model fit and save the results in an R object named phytoCurve.

phytoCurve <- nls(phytoGrowthRate ~ a*ironConcentration / ( b+ironConcentration), 
                  data = phytoplankton, list(a = 1, b = 1))

# Obtain the parameter estimates using the summary command to , including 
# standard errors and t-tests of null hypotheses that parameter values are zero.

summary(phytoCurve)

# Add the nonlinear regression curve to scatter plot. If necessary, redraw the scatter plot before issuing the following commands.

xpts <- seq(min(phytoplankton$ironConcentration), 
            max(phytoplankton$ironConcentration), length.out = 100)
ypts <- predict(phytoCurve, new = data.frame(ironConcentration = xpts))
lines(xpts, ypts)

# Many of the functions that can be used to extract results from a saved lm object
# work in the same way when applied to an nls object, such as predict, residuals, 
# and coef.

# ------------------------------------------------------------

# Figure 17.8-2. <a href="../wp-content/data/chapter17/chap17f8_2PondPlantsAndProductivity.csv">Pond plants and productivity
# Fit a quadratic curve to the relationship between the number of plant species present 
# in ponds and pond productivity.
# Read and examine data.

pondproductivity <- read.csv("pondproductivity.csv")

# Scatter plot.

plot(species ~ productivity, data = pondProductivity)

# Commands for a fancier scatter plot are here.

plot(species ~ productivity, data = pondProductivity, pch = 16, 
     col = "firebrick", las = 1, cex = 1.5, bty = "l", 
     ylab = "Number of species", xlab = "Productivity (g/15 days)" )

#   	
# Fit a quadratic curve to the data. Here, the single variable productivity in the data frame is included in the formula both as itself and as the squared term, productivity2. To make the squared term work, we need to wrap the term with I(). The results of the model fit are saved in an R object productivityCurve.

productivityCurve <- lm(species ~ productivity + I(productivity^2), 
                        data = pondProductivity)

# Show estimates of the parameters of the quadratic curve (regression coefficients) are obtained as follows, along with standard errors and t-tests.

summary(productivityCurve)

# Add quadratic regression curve to scatter plot. If necessary, redraw the previous scatter plot before issuing the following commands.

xpts <- seq(min(pondProductivity$productivity), max(pondProductivity$productivity), 
            length.out = 100)
ypts <- predict(productivityCurve, new = data.frame(productivity = xpts))
lines(xpts, ypts)

# ------------------------------------------------------------

# Example 17.8. <a href="../wp-content/data/chapter17/chap17e8ShrinkingSeals.csv">Incredible shrinking seal
# Fit a formula-free curve (cubic spline) to the relationship between body 
# length and age for female fur seals.
# Read and inspect data

shrink <- read.csv("shrink.csv")

# Scatter plot.

plot(length ~ ageInDays, data = shrink, pch = ".")

# Commands for a fancier scatter plot with more options are here.

plot(jitter(length, factor = 2) ~ ageInDays, data = shrink, pch = ".", 
     col = "firebrick", las = 1, bty = "l", ylab = "Body length (cm)",
     xlab = "Female age (days)")

#   	
# Fit a cubic spline. The argument df stands for effective degrees of freedom, which allows you to control how complicated the curve should be. The simplest possible curve is a straight line, which has df=2 (one for slope and another for intercept). More complex curves require more df. Here we fit a very complicated curve.

shrinkCurve <- smooth.spline(shrink$ageInDays, shrink$length, df = 40)

# Add curve to scatter plot. If needed, replot the scatter plot before issuing the commands below.

xpts <- seq(min(shrink$ageInDays), max(shrink$ageInDays), length.out = 1000)
ypts <- predict(shrinkCurve, xpts)$y
lines(xpts, ypts)

# ------------------------------------------------------------

# Figure 17.9-1. <a href="../wp-content/data/chapter17/chap17f9_1GuppyColdDeath.csv">Guppy cold death
# Fit a logistic regression to the relationship between guppy mortality and
# duration of exposure to a temperature of 5 degrees C.
# Read and inspect the data.

guppy <- read.csv("guppy.csv")

# Draw a frequency table of mortality at different exposure times.

table(guppy$mortality, guppy$exposureDurationMin, 
      dnn = c("Mortality","Exposure (min)"))

# Scatter plot of the data.

plot(mortality ~ jitter(exposureDurationMin), data = guppy)

# Commands for a fancier scatter plot, with more options, are here.

plot(jitter(mortality, factor = 0.1) ~ jitter(exposureDurationMin, factor = 1),  
     data = guppy, col = "firebrick", las = 1, bty = "l", ylab = "Mortality",
     xlab = "Duration of exposure (min)")

#   	
# Fit a logistic regression.

guppyGlm <- glm(mortality ~ exposureDurationMin, data = guppy,
                family = binomial(link = logit))

# Add logistic regression curve to scatter plot.

xpts <- seq(min(guppy$exposureDurationMin), max(guppy$exposureDurationMin), 
            length.out = 100)
ypts <- predict(guppyGlm, newdata = data.frame(exposureDurationMin = xpts), 
                type = "response")
lines(xpts, ypts)

# Table of regression coefficients, with parameter estimates (Table 17.9-2).

summary(guppyGlm)

# 95% confidence intervals for parameter estimates. It is necessary to load the MASS package first. 

library(MASS)
confint(guppyGlm)

# Predict probability of mortality (mean mortality) for a given x-value, 10 min duration, including standard error of prediction.

predict(guppyGlm, newdata = data.frame(exposureDurationMin = 10),
        type = "response", se.fit = TRUE)

# Estimate the LD50, the dose at which half the individuals are predicted to die from exposure.

library(MASS)
dose.p(guppyGlm, p = 0.50)

# Analysis of deviance table, with a of the null hypothesis of zero slope (Table 17.9-3).

anova(guppyGlm, test = "Chi")









# Multiple explanatory variables: R code for Chapter 18 examples
# Download the R code on this page as a single file <a href="../wp-content/rcode/chap18.r">here
# ------------------------------------------------------------

# Figure 18.1-1 <a href="../wp-content/data/chapter17/chap17e3PlantDiversityAndStability.csv">
# Modeling with linear regression
# Compare the fits of the null and univariate regression models to data on the 
# relationship between stability of plant biomass production and the initial number 
# of plant species assigned to plots. The data are from Example 17.3.
# Read and inspect data.

prairie <- read.csv(url("http://www.zoology.ubc.ca/~schluter/WhitlockSchluter/wp-content/data/chapter17/chap17e3PlantDiversityAndStability.csv"))
head(prairie)

# Take the log-transformation of stability.

prairie$logStability <- log(prairie$biomassStability)
head(prairie)

# Fit the null model to the data, which in simple linear regression 
# is a line of 0 slope.

prairieNullModel <- lm(logStability ~ 1, data = prairie)

# Fit the full model, which includes the treatment variable.

prairieRegression <- lm(logStability ~ nSpecies, data = prairie)

# Scatter plot to compare models visually.

plot(logStability ~ nSpecies, data = prairie, bty = "l", col="firebrick", 
     pch = 1, las = 1, cex = 1.5, xlim = c(0,16), xaxp = c(0,16,8), 
     xlab = "Species number treatment", 
     ylab = "Log-transformed ecosystem stability")
abline(prairieNullModel, lty = 2)
abline(prairieRegression)

# The F-test of improvement in fit of the full (regression) model.

anova(prairieRegression)

# ------------------------------------------------------------

# Figure 18.1-2. <a href="../wp-content/data/chapter15/chap15e1KneesWhoSayNight.csv">Generalizing linear regression
# Compare the fits of the null and single-factor ANOVA model to data on phase shift
# in the circadian rhythm of melatonin production in participants given alternative
# light treatments. The data are from Example 15.1. 
# Read and inspect data.

circadian <- read.csv(url("http://www.zoology.ubc.ca/~schluter/WhitlockSchluter/wp-content/data/chapter15/chap15e1KneesWhoSayNight.csv"))

# Set the order of groups for tables and graphs.

circadian$treatment <- factor(circadian$treatment, 
                              levels = c("control", "knee","eyes")) 

# Fit the null model to the data. In single-factor ANOVA, this fits a constant (grand mean) to all groups.

circadianNullModel <- lm(shift ~ 1, data = circadian)

# Fit the full model to the data, which includes the treatment effect.

circadianAnova <- lm(shift ~ treatment, data = circadian)

# Strip chart to compare models visually. The dashed line is the null model. The value adjustAmount is used to control the width of the horizontal line segments.

par(bty="l")
adjustAmount <- 0.15
stripchart(shift ~ treatment, data = circadian, method = "jitter",
           vertical = TRUE, las = 1, pch = 1, xlab = "Light treatment",
           ylab = "Shift in circadian rhythm (h)", col = "firebrick", 
           cex = 1.2)

xpts <- as.numeric(circadian$treatment)
ypts <- predict(circadianNullModel)
segments(xpts - adjustAmount, ypts, xpts + adjustAmount, ypts, lty = 2)

xpts <- as.numeric(circadian$treatment)
ypts <- predict(circadianAnova)
segments(xpts - adjustAmount, ypts, xpts + adjustAmount, ypts, col = "firebrick")

# F-test of improvement in fit of the full model.

anova(circadianAnova)

# ------------------------------------------------------------

# Example 18.2. <a href="../wp-content/data/chapter18/chap18e2ZooplanktonDepredation.csv">Zooplankton depredation
# Analyze data from a randomized block experiment, which measured the effects of 
# fish abundance (the factor of interest) on zooplankton diversity.
# Treatments were repeated at 5 locations in a lake (blocks). 
# Read and inspect data.

zooplankton <- read.csv(url("http://www.zoology.ubc.ca/~schluter/WhitlockSchluter/wp-content/data/chapter18/chap18e2ZooplanktonDepredation.csv"))
head(zooplankton)

# Set the order of groups in tables and plots.

zooplankton$treatment <- factor(zooplankton$treatment, levels = c("control","low","high"))

# Table of the data (Table 18.2-1). One measurement is available from each combination of treatment and block

tapply(zooplankton$diversity, list(Treatment = zooplankton$treatment, 
                                   Location = zooplankton$block), unique)

# A blocking variable is typically analyzed as a random effect. We need to load the nlme package.

library(nlme)

# Fit the null model, which includes block but leaves out the treatment variable.

zoopNullModel <- lme(diversity ~ 1, random = ~ 1| block, data = zooplankton)

# Fit the full model, which includes treatment.

zoopRBModel <- lme(diversity ~ treatment, random = ~ 1| block, data = zooplankton)

# 	
# Visualize model fits, beginning with the null model. The result here differs from that in the book, which shows block means. Lines are coincident when we analyze the data using lme, which estimates hardly any variance among the block means.

interaction.plot(zooplankton$treatment, zooplankton$block, 
                 response = predict(zoopNullModel), ylim = range(zooplankton$diversity),
                 trace.label = "Block", las = 1)
points(zooplankton$treatment, zooplankton$diversity, 
       pch = as.numeric(zooplankton$block))

# Visualize model fits, continuing with the full model. With the treatment term included, lme estimates that there is indeed variation among block means.

interaction.plot(zooplankton$treatment, zooplankton$block, 
                 response = predict(zoopRBModel), ylim = range(zooplankton$diversity),
                 trace.label = "Block", las = 1)
points(zooplankton$treatment, zooplankton$diversity, 
       pch = as.numeric(zooplankton$block))

# F--test of improvement in fit of the full model. This represents the test of treatment effect. Notice that R will not test the random (block) effect of an lme model.
anova(zoopRBModel)

# ------------------------------------------------------------

# Example 18.3 <a href="../wp-content/data/chapter18/chap18e3IntertidalAlgae.csv">Intertidal interaction zone
# Analyze data from a factorial experiment investigating the effects of herbivore
# presence, height above low tide, and the interaction between these factors,
# on abundance of a red intertidal alga using two-factor ANOVA.
# Read and inspect the data.

algae <- read.csv(url("http://www.zoology.ubc.ca/~schluter/WhitlockSchluter/wp-content/data/chapter18/chap18e3IntertidalAlgae.csv"))
head(algae)

# Fit the null model having both main effects but no interaction term. Note that we use lm because both factors are fixed effects.

algaeNoInteractModel <- lm(sqrtArea ~ height + herbivores, data = algae)

# Fit the full model, with interaction term included.

algaeFullModel <- lm(sqrtArea ~ height * herbivores, data = algae)

# Visualize the model fits, beginning with the no-interaction model.

interaction.plot(algae$herbivores, algae$height, response = predict(algaeNoInteractModel), 
                 ylim = range(algae$sqrtArea), trace.label = "Height", las = 1,
                 ylab = "Square root surface area (cm)", xlab = "Herbivore treatment")
adjustAmount = 0.05
points(sqrtArea ~ c(jitter(as.numeric(herbivores), factor = 0.2) + adjustAmount), 
       data = subset(algae, height == "low"))
points(sqrtArea ~ c(jitter(as.numeric(herbivores), factor = 0.2) - adjustAmount), 
       data = subset(algae, height == "mid"), pch = 16)

# Visualize the model fits, continuing with the full model including the interaction term.

interaction.plot(algae$herbivores, algae$height, response = predict(algaeFullModel), 
                 ylim = range(algae$sqrtArea), trace.label = "Height", las = 1,
                 ylab = "Square root surface area (cm)", xlab = "Herbivore treatment")
adjustAmount = 0.05
points(sqrtArea ~ c(jitter(as.numeric(herbivores), factor = 0.2) + adjustAmount), 
       data = subset(algae, height == "low"))
points(sqrtArea ~ c(jitter(as.numeric(herbivores), factor = 0.2) - adjustAmount), 
       data = subset(algae, height == "mid"), pch = 16)

# 	
# Test the improvement in fit of the model including the interaction term. 

anova(algaeNoInteractModel, algaeFullModel)

# Test all terms in the model in a single ANOVA table. Most commonly, this is done using either "Type III" sums of squares (see footnote 5 on p 618 of the book) or "Type I" sums of squares (which is the default in R). In the present example the two methods give the same answer, because the design is completely balanced, but this will not generally happen when the design is not balanced.
# Here is how to test all model terms using "Type III" sums of squares. We need to include a contrasts argument for the two categorical variables in the lm command. Then we need to load the car package and use its Anova command. Note that "A" is in upper case in Anova - a very subtle difference.

algaeFullModelTypeIII <- lm(sqrtArea ~ height * herbivores, data = algae, 
                            contrasts = list(height = contr.sum, herbivores = contr.sum))
library(car)
Anova(algaeFullModelTypeIII, type = "III") # note "A" in Anova is capitalized

# Here is how we test all model terms using "Type I" (sequential) sums of squares. Note that "a" is in lower case in anova.

anova(algaeFullModel)

# A residual plot from the full model.

plot( residuals(algaeFullModel) ~ fitted(algaeFullModel) )
abline(0,0)

# A normal quantile plot of the residuals.

qqnorm(residuals(algaeNoInteractModel), pch = 16, col = "firebrick", 
       las = 1, ylab = "Residuals", xlab = "Normal quantile", main = "")

# ------------------------------------------------------------

# Example 18.4 <a href="../wp-content/data/chapter18/chap18e4MoleRatLayabouts.csv">Mole-rat layabouts
# Analyze a factor while adjusting for a covariate, comparing energy expenditure of two castes of naked mole-rat while adjusting for differences in body mass using analysis of covariance ANCOVA.
# Read and inspect the data.

moleRat <- read.csv(url("http://www.zoology.ubc.ca/~schluter/WhitlockSchluter/wp-content/data/chapter18/chap18e4MoleRatLayabouts.csv"))
head(moleRat)

# We are going to sort the data according to the value of the ln of body mass. This simplifies graphing of the model fits. The graph commands below assume that the data are sorted in this way.

moleRatSorted <- moleRat[ order(moleRat$lnMass), ]

# Scatter plot of the data.

plot(lnEnergy ~ lnMass, data = moleRat, type = "n", las = 1, bty = "l")
points(lnEnergy ~ lnMass, data = subset(moleRatSorted, caste == "worker"), 
       pch = 1, col = "firebrick")
points(lnEnergy ~ lnMass, data = subset(moleRatSorted, caste == "lazy"), 
       pch = 16, col = "firebrick")

# Fit models to the data, beginning with the model lacking an interaction term. Use lm because caste and mass are fixed effects. Save the predicted values in the data frame.

moleRatNoInteractModel <- lm(lnEnergy ~ lnMass + caste, data = moleRatSorted)
moleRatSorted$predictedNoInteract <- predict(moleRatNoInteractModel)

# Fit the full model, which includes the interaction term. Again, save the predicted values in the data frame.

moleRatFullModel <- lm(lnEnergy ~ lnMass * caste, data = moleRatSorted)
moleRatSorted$predictedInteract <- predict(moleRatFullModel)

# Visualize the model fits, beginning with the fit of the no-interaction model. Redraw the scatter plot (see commands above), if necessary, before issuing the following commands.

lines(predictedNoInteract ~ lnMass, data = subset(moleRatSorted, 
                                                  caste == "worker"), lwd = 1.5)
lines(predictedNoInteract ~ lnMass, data = subset(moleRatSorted, 
                                                  caste == "lazy"), lwd = 1.5)

# Visualize the fit of the full model, including the interaction term. Redraw the scatter plot, if necessary, before issuing the following commands.

lines(predictedInteract ~ lnMass, data = subset(moleRatSorted, 
                                                caste == "worker"), lwd = 1.5, lty = 2)
lines(predictedInteract ~ lnMass, data = subset(moleRatSorted, 
                                                caste == "lazy"), lwd = 1.5, lty = 2)

# Test the improvement in fit of the full model, including the interaction term. This is a test of the interaction term only.

anova(moleRatNoInteractModel, moleRatFullModel)

# Test for differences in ln body mass between castes, assuming that no interaction term is present in the mole rat population (i.e., assuming that the two castes hve equal slopes). Most commonly this is done using either "Type III" sums of squares (see footnote 5 on p 618 of the book) or "Type I" sums of squares (the default in R). The two methods do not give identical answers here because the design is not balanced (in a balanced design, each value of the x-variable would have the same number of y-observations from each group).
# Test using "Type III" sums of squares. We need to include a contrasts argument for the categorical variable in the lm command. Then we need to load the car package and use its Anova command. Note that "A" is in upper case in Anova().

moleRatNoInteractModelTypeIII <- lm(lnEnergy ~ lnMass + caste, data = moleRat, 
                                    contrasts = list(caste = contr.sum))
library(car)
Anova(moleRatNoInteractModelTypeIII, type = "III") # note "A" in Anova is capitalized

# Test using "Type I" (sequential) sums of squares. Make sure that the covariate (lnMass) comes before the factor (caste) in the lm formula, as shown. Note that "a" is in lower case in anova.

moleRatNoInteractModel <- lm(lnEnergy ~ lnMass + caste, data = moleRat)
anova(moleRatNoInteractModel)

# Residual plot from the linear model. 

plot( residuals(moleRatNoInteractModel) ~ fitted(moleRatNoInteractModel) )
abline(0,0)

# Normal quantile plot of residuals.

qqnorm(residuals(moleRatNoInteractModel), pch = 16, col = "firebrick", 
       las = 1, ylab = "Residuals", xlab = "Normal quantile", main = "")


