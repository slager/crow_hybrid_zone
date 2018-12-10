library(magrittr)
library(nlme)

read.csv("random effects variance.csv",header=T,
         colClasses=c("character","factor",
                      "numeric","numeric","factor")) -> df

m1 <- lme(N ~ nbc_cbc, random = ~1|pop, data=df)
m2 <- lme(N ~ nbc_cbc, random = ~1|pop, weights=varIdent(form=~1|nbc_cbc), data=df)
anova(m1, m2) -> a
m1
m2 # See variance function parameters, 9x
a
summary(a)
str(a)
a$`p-value`
a$L.Ratio

#Compare 2 nested models.  One with a random effect for population.
#Another with random effect for population, and variance function allowing for different variances between nbc_cbc, and all other populations