
set.seed(51315)
# library(lme4)
# data(grouseticks)
# grouseticks$HEIGHT <- scale(grouseticks$HEIGHT)
# grouseticks <- merge(grouseticks, grouseticks_agg[, 1:3], by = "BROOD")
# #
# # Build out models
# form <- TICKS ~ YEAR + HEIGHT +(1|BROOD) + (1|INDEX) + (1|LOCATION)
# glmer3Lev  <- glmer(form, family="poisson",data=grouseticks,
#                     control = glmerControl(optimizer="Nelder_Mead",
#                                            optCtrl=list(maxfun = 1e5)))
# # GLMER 3 level + slope
# form <- TICKS ~ YEAR + HEIGHT +(1 + HEIGHT|BROOD) + (1|LOCATION) + (1|INDEX)
# glmer3LevSlope  <- glmer(form, family="poisson",data=grouseticks,
#                          control = glmerControl(optimizer="bobyqa",
#                                                 optCtrl=list(maxfun = 1e5)))
#
# # GLMER 2 level
# # data(VerbAgg)
# # fmVA <- glmer(r2 ~ Anger + Gender + btype + situ +
# #                 (1|id) + (1|item), family = binomial, data =
# #                 VerbAgg)
#
# # Sleepstudy
# lmerSlope1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#
#
#
# predict(glmer3Lev, newdata = grouseticks[2, ])
# predict(glmer3Lev, newdata = grouseticks[2, ], re.form = NA)
# predict(lmerSlope1, newdata=sleepstudy[10,])
# predict(lmerSlope1, newdata=sleepstudy[10,], re.form=NA)
