# I'm testing merMod objects

library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
summary(fm1)# (with its own print method)

REextract(fm1)
