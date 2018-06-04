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
data(hsb)

# Create a function to randomly assign NA

add_NA <- function(x, prob){
  z <- rbinom(length(x), 1, prob = prob)
  x[z==1] <- NA
  return(x)
}

hsb$minority <- add_NA(hsb$minority, prob = 0.05)
table(is.na(hsb$minority))

hsb$female <- add_NA(hsb$female, prob = 0.05)
table(is.na(hsb$female))

hsb$ses <- add_NA(hsb$ses, prob = 0.05)
table(is.na(hsb$ses))

hsb$size <- add_NA(hsb$size, prob = 0.05)
table(is.na(hsb$size))


## ----impute, message=FALSE-----------------------------------------------
# Load imputation library
library(Amelia)
# Declare the variables to include in the imputation data
varIndex <- names(hsb)
# Declare ID variables to be excluded from imputation
IDS <- c("schid", "meanses")
# Imputate
impute.out <- amelia(hsb[, varIndex], idvars = IDS, 
                         noms = c("minority", "female"), 
                         m = 5)
summary(impute.out)

## ------------------------------------------------------------------------
fmla <- "mathach ~ minority + female + ses + meanses + (1 + ses|schid)"
mod <- lmer(fmla, data = hsb)
modList <- lmerModList(fmla, data = impute.out$imputations)

## ------------------------------------------------------------------------
fixef(mod) # model with dropped missing
fixef(modList)

## ------------------------------------------------------------------------
VarCorr(mod) # model with dropped missing
VarCorr(modList) # aggregate of imputed models

## ------------------------------------------------------------------------
lapply(modList, fixef)

## ------------------------------------------------------------------------
fixef(modList[[1]])
fixef(modList[[2]])

## ------------------------------------------------------------------------
print(modList)

## ------------------------------------------------------------------------
summary(modList)

## ------------------------------------------------------------------------
fastdisp(modList)

## ------------------------------------------------------------------------
modelRandEffStats(modList)
modelFixedEff(modList)
VarCorr(modList)

## ------------------------------------------------------------------------
modelInfo(mod)

## ------------------------------------------------------------------------
lapply(modList, modelInfo)

