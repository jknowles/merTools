# Test REmargins

fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

context("Test random effect marginalization works")

mfx <- REmargins(merMod = fm1, newdata = sleepstudy[1:10,])


ggplot(mfx) + aes(x = case, y = fit_Subject) +
  geom_line() +
  facet_wrap(~obs)
