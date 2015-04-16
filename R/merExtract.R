#' @title REextract
#' @name Extracts random effects
#' @description Extracts random effect terms from an lme4 model
#' @param mod a merMod object from the lme4 package
#' @importFrom plyr adply
#' @return a data frame with random effects from lme4
#' @export
REextract <- function(mod){
  out <- lme4::ranef(mod, condVar = TRUE)
  out.se <- plyr::adply(attr(out[[1]], which = "postVar"), c(3),
                        function(x) sqrt(diag(x)))
  out.pt <- out[[1]]
  names(out.se)[-1] <- paste0(names(out.pt), "_se")
  newout <- cbind(out.pt, out.se[, -1])
  return(newout)
}
#
#
# Use the Gelman sim technique to build empirical Bayes estimates
# Uses the sim function in the arm package
REsim <- function(x, nsims, OR = FALSE){
  mysim <- sim(x, n.sims = nsims)
  reDims <- length(mysim@ranef)
  tmp.out <- vector("list", reDims)
  for(i in c(1:reDims)){
    tmp.out[[i]] <- plyr::adply(mysim@ranef[[i]], c(2, 3), plyr::each(c(mean, median, sd)))
    tmp.out[[i]]$level <- paste0("Level ", i)
    tmp.out[[i]]$level <- as.character(tmp.out[[i]]$level)
    tmp.out[[i]]$X1 <- as.character(tmp.out[[i]]$X1)
    tmp.out[[i]]$X2 <- as.character(tmp.out[[i]]$X2)
  }
  dat <- do.call(rbind, tmp.out)
  if(OR == TRUE){
    dat$median <- exp(dat$median)
    dat$mean <- exp(dat$mean)
    dat$sd <- NA # don't know how to do SE of odds ratios currently
    return(dat)
  } else{
    return(dat)
  }
}
#
#
# FEsim <- function(x, nsims){
#   mysim <- sim(x, n.sims = nsims)
#   means <- apply(mysim@fixef, MARGIN = 2, mean)
#   medians <- apply(mysim@fixef, MARGIN = 2, median)
#   sds <- apply(mysim@fixef, MARGIN =2, sd)
#   dat <- data.frame(var = names(means), meanEff = means, medEff = medians,
#                     sdEff = sds, row.names=NULL)
#   return(dat)
# }
#
# # http://stats.stackexchange.com/questions/56525/standard-deviation-of-random-effect-is-0
# # http://stat.columbia.edu/~jcliu/paper/HierarchicalPrior.pdf
# # http://www.stat.columbia.edu/~radon/paper/paper.pdf (example)
#
# # Plot the effects
# plotMCMCre <- function(dat, scale, var, sd, sigmaScale = NULL, oddsRatio = FALSE,
#                        labs = NULL){
#   if(!missing(sigmaScale)){
#     dat[, sd] <- dat[, sd] / sigmaScale
#     dat[, var] <- dat[, var] / sigmaScale
#   }
#   dat[, sd] <- dat[, sd] * scale
#   dat[, "ymax"] <- dat[, var] + dat[, sd]
#   dat[, "ymin"] <- dat[, var] - dat[, sd]
#   hlineInt <- 0
#   if(oddsRatio == TRUE){
#     dat[, "ymax"] <- exp(dat[, "ymax"])
#     dat[, var] <- exp(dat[, var])
#     dat[, "ymin"] <- exp(dat[, "ymin"])
#     hlineInt <- 1
#   }
#   if(!missing(labs)){
#     xlabs.tmp <- element_text(face = "bold")
#     xvar <- labs
#   } else {
#     xlabs.tmp <- element_blank()
#     xvar <- "id"
#   }
#
#   dat[order(dat[, var]), "id"] <- c(1:nrow(dat))
#   ggplot(dat, aes_string(x = xvar, y = var, ymax = "ymax",
#                          ymin = "ymin")) +
#     geom_pointrange(alpha = I(0.4)) + theme_dpi() + geom_point() +
#     labs(x = "Group", y = "Effect Range", title = "Effect Ranges") +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           axis.text.x = xlabs.tmp, axis.ticks.x = element_blank()) +
#     geom_hline(yintercept = hlineInt, color = I("red"), size = I(1.1))
# }
