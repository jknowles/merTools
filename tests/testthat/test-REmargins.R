# Test REmargins

test_that("Text marginalized effects object has the correct dimensions", {
  skip_on_travis()
  skip_on_cran()
  set.seed(51315)
  fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  #

  #
  mfx <- REmargins(merMod = fm1, newdata = sleepstudy[1:10,])
  #
  # mfx has a row for each unique combo of row in newdata, breaks, grouping_var, and term
  expect_equal(nrow(mfx), 10 * length(unique(mfx$breaks)) * length(unique(mfx$grouping_var)) *
                 length(unique(mfx$term)))
  expect_equal(ncol(mfx), 17)


})

# ggplot(out_w) + aes(x = obs, y = fit_Subject) +
#   geom_line() +
#   facet_wrap(~case)

