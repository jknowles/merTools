library('lme4')
mydata <- Penicillin
set.seed(10)
mydata$income <- rnorm(144)
mydata$school <- Penicillin$plate
mydata$state <- Penicillin$sample
mydata$SAT <- Penicillin$diameter
m1 <- lm(SAT ~ income, data=mydata[mydata$state=="A",])
summary(m1)$coefficients
fm <- lmer(SAT ~ income + (income|school) + (1|state), data=mydata)
