# Import variables from function
library(ggplot2)
library(shiny)

merMod <<- .shinyMerPar$merMod
if(is.null(.shinyMerPar$simData)){

} else {
  newdata <<- .shinyMerPar$simData
}

if (!exists("newdata")) {
  df.choices <- c("Model Frame"   = "orig",
                  "Random Obs"    = "rand",
                  "Average Obs"   = "mean")
} else {
  df.choices <- c("User Supplied" = "user",
                  "Model Frame"   = "orig",
                  "Random Obs"    = "rand",
                  "Average Obs"   = "mean")
}


