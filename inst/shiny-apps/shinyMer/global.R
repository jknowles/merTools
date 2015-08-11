# Import variables from function

merMod <<- .shinyMerPar$merMod
if(is.null(.shinyMerPar$simData)){

} else {
  newdata <<- .shinyMerPar$simData
}

if (!exists("simData")) {
  df.choices <- c("Model Frame"   = "orig",
                  "Random Obs"    = "rand",
                  "Average Obs"   = "mean")
} else {
  df.choices <- c("User Supplied" = "user",
                  "Model Frame"   = "orig",
                  "Random Obs"    = "rand",
                  "Average Obs"   = "mean")
}


