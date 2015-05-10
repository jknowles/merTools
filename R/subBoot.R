# subBoot
# Sleepstudy


#' Title
#'
#' @param model
#'
#' @return
#' @export
#'
#' @examples
# thetaExtract <- function(model){
#   return(model@theta)
# }

#' Title
#'
#' @param model
#' @param n
#' @param FUN
#' @param B
#'
#' @return
#' @export
#'
#' @examples
# subBoot <- function(model, n, FUN, B){
#   resultMat <- matrix(FUN(model), nrow = 1)
#   tmp <- matrix(data=NA, nrow=B, ncol=ncol(resultMat))
#   resultMat <- rbind(resultMat, tmp); rm(tmp)
#     for(i in 1:B){
#       rows <- as.numeric(row.names(model@frame))
#       mysamp <- as.character(sample(rows,  n, replace=TRUE))
#       # http://proteo.me.uk/2013/12/fast-subset-selection-by-row-name-in-r/
#       newdata <- model@frame[match(mysamp, rownames(model@frame)),]
#       # Only for lmerMod
#       tmpMod <- lmer(formula(model), data = newdata)
#       resultMat[i + 1, ] <- FUN(tmpMod)
#     }
#   resultMat <- as.data.frame(resultMat)
#   resultMat$sample <- c("original", 1:B)
#   return(resultMat)
# }


