#' @title Clean up variable names in data frames
#' @name sanitizeNames
#' @description Strips out transformations from variable names in data frames
#' @param data a data.frame
#' @return a data frame with variable names cleaned to remove factor() construction
sanitizeNames <- function(data){
  badFac <- grep("factor\\(", names(data))
  for(i in badFac){
    names(data)[i] <- gsub("factor\\(", "", names(data)[i])
    names(data)[i] <- gsub("\\)", "", names(data)[i])
  }
  row.names(data) <- NULL
  return(data)
}

#' @title Remove attributes from a data.frame
#' @name stripAttributes
#' @description Strips attributes off of a data frame that come with a merMod model.frame
#' @param data a data.frame
#' @return a data frame with variable names cleaned to remove all attributes except for
#' names, row.names, and class
stripAttributes <- function(data){
  attr <- names(attributes(data))
  good <- c("names", "row.names", "class")
  for(i in attr[!attr %in% good]){
    attr(data, i) <- NULL
  }
  return(data)
}


#' @title Draw a single observation out of an object matching some criteria
#' @name draw
#' @description Draw is used to select a single observation out of an R object.
#' Additional parameters allow the user to control how that observation is
#' chosen in order to manipulate that observation later. This is a generic
#' function with methods for a number of objects.
#' @param object the object to draw from
#' @param type what kind of draw to make. Options include random or average
#' @param varList a list specifying filters to subset the data by when making the
#' draw
#' @param seed numeric, optional argument to set seed for simulations, ignored if type="average"
#' @return a data.frame with a single row representing the desired observation
#' @details In cases of tie, ".", may be substituted for factors.
#' @export draw
#' @rdname draw
draw <- function(object, type = c("random", "average"),
                 varList = NULL, seed = NULL){
  UseMethod("draw")
}

#' @title Draw an observation from a merMod object
#' @rdname draw
#' @method draw merMod
#' @export
#' @import lme4
#' @examples
#' fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' # Random case
#' draw(fm1, type = "random")
#' # Average
#' draw(fm1, type = "average")
#' # Subset
#' draw(fm1, type = "average", varList = list("Subject" = "308"))
#'
draw.merMod <- function(object, type = c("random", "average"),
                        varList = NULL, seed = NULL){
  type <- match.arg(type, c("random", "average"), several.ok = FALSE)
  if(type == 'random'){
    out <- randomObs(object, varList, seed)
    return(out)
  } else if(type == 'average'){
    out <- averageObs(object, varList)
    return(out)
  }
}

#' @title Select a random observation from model data
#' @name randomObs
#' @description Select a random observation from the model frame of a merMod
#' @param merMod an object of class merMod
#' @param varList optional, a named list of conditions to subset the data on
#' @param seed numeric, optional argument to set seed for simulations
#' @return a data frame with a single row for a random observation, but with full
#' factor levels. See details for more.
#' @details Each factor variable in the data frame has all factor levels from the
#' full model.frame stored so that the new data is compatible with predict.merMod
randomObs <- function(merMod, varList, seed = NULL){
  if(!missing(varList)){
    data <- subsetList(merMod@frame, varList)
  }

  if (!is.null(seed))
    set.seed(seed)
  else if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)

  out <- data[sample(1:nrow(data), 1),]
  chars <- !sapply(out, is.numeric)
  for(i in names(out[, chars])){
    out[, i] <- superFactor(out[, i], fullLev = unique(merMod@frame[, i]))
  }
  out <- stripAttributes(out)
  return(out)
}

#' @title Collapse a dataframe to a single average row
#' @name collapseFrame
#' @description Take an entire dataframe and summarize it in one row by using the
#' mean and mode.
#' @param data a data.frame
#' @return a data frame with a single row
#' @details Each character and factor variable in the data.frame is assigned to the
#' modal category and each numeric variable is collapsed to the mean. Currently if
#' mode is a tie, returns a "."
collapseFrame <- function(data){
  chars <- !sapply(data, is.numeric)
  chars <- names(data[, chars, drop = FALSE])
  nums <- sapply(data, is.numeric)
  nums <- names(data[, nums, drop = FALSE])
  numDat <- apply(data[, nums, drop = FALSE], 2, mean)
  statmode <- function(x){
    z <- table(as.vector(x))
    m <- names(z)[z == max(z)]
    if (length(m) == 1) {
      return(m)
    }
    return(".")
  }
  charDat <- apply(data[, chars, drop = FALSE], 2, statmode)
  cfdata <- cbind(as.data.frame(t(numDat)), as.data.frame(t(charDat)))
  cfdata <- cfdata[, names(data)]
  return(cfdata)
}


#' @title Subset a data.frame using a list of conditions
#' @name subsetList
#' @description Split a data.frame by elements in a list
#' @param data a data.frame
#' @param list a named list of splitting conditions
#' @return a data frame with values that match the conditions in the list
subsetList <- function(data, list){
  for(i in names(list)){
    data <- split(data, data[, i])
    data <- data[[list[[i]]]]
    data <- as.data.frame(data)
  }
  return(data)
}

#' @title Find the average observation for a merMod object
#' @name averageObs
#' @description Extract a data frame of a single row that represents the
#' average observation in a merMod object. This function also allows the
#' user to pass a series of conditioning argument to calculate the average
#' observation conditional on other characteristics.
#' @param merMod a merMod object
#' @param varList optional, a named list of conditions to subset the data on
#' @return a data frame with a single row for the average observation, but with full
#' factor levels. See details for more.
#' @details Each character and factor variable in the data.frame is assigned to the
#' modal category and each numeric variable is collapsed to the mean. Currently if
#' mode is a tie, returns a "." Uses the collapseFrame function.
averageObs <- function(merMod, varList = NULL){
  if(!missing(varList)){
    data <- subsetList(merMod@frame, varList)
    if(nrow(data) < 20 & nrow(data) > 2){
      warning("Subset has less than 20 rows, averages may be problematic.")
    }
    if(nrow(data) <3){
      warning("Subset has fewer than 3 rows, computing global average instead.")
      data <- merMod@frame
    }
  } else{
    data <- merMod@frame
  }
  out <- collapseFrame(data)
  reTerms <- names(ngrps(merMod))
  for(i in 1:length(reTerms)){
    out[, reTerms[i]] <- REquantile(merMod = merMod,
                                        quantile = 0.5, groupFctr = reTerms[[i]])
    out[, reTerms[i]] <- as.character(out[, reTerms[i]])
  }
  chars <- !sapply(out, is.numeric)
  for(i in names(out[, chars])){
    out[, i] <- superFactor(out[, i], fullLev = unique(merMod@frame[, i]))
  }
  out <- stripAttributes(out)
  out <- out[, names(merMod@frame)]
  return(out)
}


#' @title Create a factor with unobserved levels
#' @name superFactor
#' @description Create a factor variable and include unobserved levels
#' for compatibility with model prediction functions
#' @param x a vector to be converted to a factor
#' @param fullLev a vector of factor levels to be assigned to x
#' @return a factor variable with all observed levels of x and all levels
#' of x in fullLev
#' @export
#' @examples
#' regularFactor <- c("A", "B", "C")
#' regularFactor <- factor(regularFactor)
#' levels(regularFactor)
#' # Now make it super
#' newLevs <- c("D", "E", "F")
#' regularFactor <- superFactor(regularFactor, fullLev = newLevs)
#' levels(regularFactor) # now super
superFactor <- function(x, fullLev){
  x <- as.character(x)
  if("factor" %in% class(fullLev)){
    fullLev <- unique(levels(fullLev))
  }
  unobsLev <- unique(x)[!unique(x) %in% fullLev]
  x <- factor(x, levels = c(fullLev, unobsLev),
              labels = c(fullLev, unobsLev))
  return(x)
}

#' @title Randomly reorder a dataframe
#' @name shuffle
#' @description Randomly reorder a dataframe by row
#' @param data a data frame
#' @return a data frame of the same dimensions with the rows reordered
#' randomly
shuffle <- function(data){
  return(data[sample(nrow(data)),])
}

#' @title Assign an observation to different values
#' @name wiggle
#' @description Creates a new data.frame with copies of the original observation,
#' each assigned to a different user specified value of a variable. Allows the
#' user to look at the effect of changing a variable on predicted values.
#' @param data a data frame with one or more observations to be reassigned
#' @param var a character specifying the name of the variable to adjust
#' @param values a vector with the variables to assign to var
#' @return a data frame with each row in data assigned to all values for
#' the variable chosen
#' @details If the variable specified is a factor, then wiggle will return it
#' as a character.
#' @export
#' @examples
#' data(iris)
#' wiggle(iris[3,], "Sepal.Width", values = c(1, 2, 3, 5))
#' wiggle(iris[3:5,], "Sepal.Width", values = c(1, 2, 3, 5))
wiggle <- function(data, var, values){
  tmp.data <- data
  while(nrow(data) < length(values) * nrow(tmp.data)){
    data <- rbind(data, tmp.data)
  }
  data[, var] <- rep(values, each = nrow(tmp.data))
  if(any(class(tmp.data[, var]) %in% c("factor", "ordered"))){
    data[, var] <- superFactor(data[, var],
                               fullLev = levels(tmp.data[, var]))

  }
  return(data)
}

#' @title Identify group level associated with RE quantile
#' @name REquantile
#' @description For a user specified quantile (or quantiles) of the random effect
#' terms in a merMod object. This allows the user to easily identify the obsevation
#' associated with the nth percentile effect.
#' @param merMod a merMod object with one or more random effect levels
#' @param quantile a numeric vector with values between 0 and 100 for quantiles
#' @param groupFctr a character of the name of the random effect grouping factor to extract
#' quantiles from
#' @param term a character of the random effect to extract for the grouping factor
#' specified. Default is the intercept.
#' @return a vector of the level of the random effect grouping term that corresponds
#' to each quantile
#' @export
#' @examples
#' fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' REquantile(fm1, quantile = 0.25, groupFctr = "Subject")
#' REquantile(fm1, quantile = 0.25, groupFctr = "Subject", term = "Days")
REquantile <- function(merMod, quantile, groupFctr, term = "(Intercept)"){
  if(any(quantile > 1 | quantile < 0)){
    stop("Quantiles must be specified on the range 0-1")
  }
  myRE <- ranef(merMod)[[groupFctr]]
  if(is.null(myRE)){
    stop("Random effect group name not found. Please respecify grouping factor.")
  }
  myRE.tmp <- try(myRE[order(myRE[, term]), ,drop = FALSE], silent = TRUE)
  if(class(myRE.tmp) != "data.frame"){
    term1 <- names(myRE)[1]
    myRE.tmp <- try(myRE[order(myRE[, term1]), ,drop = FALSE], silent = TRUE)
    warning(paste0(term, " not found in random effect terms. Returning first term, ",
          term1,", for grouping factor, ", groupFctr, ", instead."))
  }
  myRE <- myRE.tmp; myRE.tmp <- NULL
  nobs <- nrow(myRE)
  if(nobs < 20){
    message("Number of observations < 20, random effect quantiles may not be well-defined.")
  }
  obsnum <- floor(quantile * nobs)
  return(rownames(myRE)[obsnum])
}
