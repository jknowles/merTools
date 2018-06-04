## ----setup, echo = FALSE, message=FALSE, warning=FALSE, results='hide'----
knitr::opts_chunk$set(
  cache=FALSE,
  comment="#>",
  collapse=TRUE, 
  echo=TRUE, 
  fig.width = 7
)
library(knitr); library(merTools)

## ------------------------------------------------------------------------
library(lme4)
head(InstEval)
str(InstEval)

## ------------------------------------------------------------------------
m1 <- lmer(y ~ service + lectage + studage + (1|d) + (1|s), data=InstEval)

## ------------------------------------------------------------------------
library(merTools)
fastdisp(m1)

## ------------------------------------------------------------------------
feEx <- FEsim(m1, 1000)
cbind(feEx[,1] , round(feEx[, 2:4], 3))

## ------------------------------------------------------------------------
library(ggplot2)
ggplot(feEx[feEx$term!= "(Intercept)", ]) + 
  aes(x = term, ymin = median - 1.96 * sd, 
      ymax = median + 1.96 * sd, y = median) + 
  geom_pointrange() + 
  geom_hline(yintercept = 0, size = I(1.1), color = I("red")) + 
  coord_flip() + 
  theme_bw() + labs(title = "Coefficient Plot of InstEval Model", 
                    x = "Median Effect Estimate", y = "Evaluation Rating")

## ------------------------------------------------------------------------
plotFEsim(feEx) + 
  theme_bw() + labs(title = "Coefficient Plot of InstEval Model", 
                    x = "Median Effect Estimate", y = "Evaluation Rating")


## ------------------------------------------------------------------------
reEx <- REsim(m1)
head(reEx)

## ------------------------------------------------------------------------
table(reEx$term)
table(reEx$groupFctr)

## ---- eval=FALSE, echo = TRUE--------------------------------------------
#  lattice::dotplot(ranef(m1, condVar=TRUE))

## ------------------------------------------------------------------------
p1 <- plotREsim(reEx)
p1

## ------------------------------------------------------------------------
example1 <- draw(m1, type = 'random')
head(example1)

## ------------------------------------------------------------------------
# predict it
predict(m1, newdata = example1)
# change values
example1$service <- "1"
predict(m1, newdata = example1)

## ------------------------------------------------------------------------
example2 <- wiggle(example1, varlist = "lectage", 
          valueslist = list(c("1", "2", "3", "4", "5", "6")))

example2

## ------------------------------------------------------------------------
example2$yhat <- predict(m1, newdata = example2)

ggplot(example2, aes(x = lectage, y = yhat)) + geom_line(aes(group = 1)) + 
  theme_bw() + ylim(c(1, 5)) +
  geom_hline(yintercept = mean(InstEval$y), linetype = 2) + 
  geom_hline(yintercept = mean(InstEval$y) + sd(InstEval$y), linetype = 3) + 
  geom_hline(yintercept = mean(InstEval$y) - sd(InstEval$y), linetype = 3)


## ------------------------------------------------------------------------
example3 <- draw(m1, type = 'average')
example3

## ------------------------------------------------------------------------
example3 <- wiggle(example1, varlist = "service", 
          valueslist = list(c("0", "1")))
example3$yhat <- predict(m1, newdata = example3)

ggplot(example3, aes(x = service, y = yhat)) + geom_line(aes(group = 1)) + 
  theme_bw() + ylim(c(1, 5)) +
  geom_hline(yintercept = mean(InstEval$y), linetype = 2) + 
  geom_hline(yintercept = mean(InstEval$y) + sd(InstEval$y), linetype = 3) + 
  geom_hline(yintercept = mean(InstEval$y) - sd(InstEval$y), linetype = 3)


## ------------------------------------------------------------------------
REquantile(m1, quantile = 0.25, groupFctr = "s")
REquantile(m1, quantile = 0.25, groupFctr = "d")

## ------------------------------------------------------------------------
example4 <- draw(m1, type = 'average')
example4 <- wiggle(example4, varlist = "s", 
                      list(REquantile(m1, quantile = seq(0.1, 0.9, .1), 
                                 groupFctr = "s")))

example4$yhat <- predict(m1, newdata = example4)

ggplot(example4, aes(x = reorder(s, -yhat), y = yhat)) + 
  geom_line(aes(group = 1)) + 
  theme_bw() + ylim(c(1, 5)) +
  geom_hline(yintercept = mean(InstEval$y), linetype = 2) + 
  geom_hline(yintercept = mean(InstEval$y) + sd(InstEval$y), linetype = 3) + 
  geom_hline(yintercept = mean(InstEval$y) - sd(InstEval$y), linetype = 3)

## ------------------------------------------------------------------------
subExample <- list(studage = "2", lectage = "4")
example5 <- draw(m1, type = 'average', varList = subExample)
example5

## ------------------------------------------------------------------------
data(VerbAgg)
m2 <- glmer(r2 ~ Anger + Gender + btype + situ +
                (1|id) + (1 + Gender|item), family = binomial, 
              data = VerbAgg)

example6 <- draw(m2, type = 'average', varList = list("id" = "149"))
example6$btype <- "scold"
example6$situ <- "self"


tempdf <- wiggle(example6, varlist = "Gender", list(c("M", "F")))
tempdf <- wiggle(tempdf, varlist = "item", 
                    list(unique(VerbAgg$item)))
tempdf$yhat <- predict(m2, newdata = tempdf, type = "response", 
                       allow.new.levels = TRUE)

ggplot(tempdf, aes(x = item, y = yhat, group = Gender)) + 
  geom_line(aes(color = Gender))+ 
  theme_bw() + ylim(c(0, 1)) + 
  theme(axis.text.x = element_text(angle = 20, hjust=1), 
        legend.position = "bottom") + labs(x = "Item", y = "Probability")


## ------------------------------------------------------------------------
exampPreds <- predictInterval(m2, newdata = tempdf, 
                              type = "probability", level = 0.8)

tempdf <- cbind(tempdf, exampPreds)

ggplot(tempdf, aes(x = item, y = fit, ymin = lwr, ymax = upr, 
                   group = Gender)) + 
  geom_ribbon(aes(fill = Gender), alpha = I(0.2), color = I("black"))+ 
  theme_bw() + ylim(c(0, 1)) + 
  theme(axis.text.x = element_text(angle = 20), 
        legend.position = "bottom")+ labs(x = "Item", y = "Probability")

## ------------------------------------------------------------------------
exampPreds <- predictInterval(m2, newdata = tempdf, 
                              type = "probability", 
                              include.resid.var = FALSE, level = 0.8)
tempdf <- cbind(tempdf[, 1:8], exampPreds)

ggplot(tempdf, aes(x = item, y = fit, ymin = lwr, ymax = upr, 
                   group = Gender)) + 
  geom_ribbon(aes(fill = Gender), alpha = I(0.2), color = I("black"))+ 
  geom_line(aes(color = Gender)) + 
  theme_bw() + ylim(c(0, 1)) + 
  theme(axis.text.x = element_text(angle = 20), 
        legend.position = "bottom") + labs(x = "Item", y = "Probability")


