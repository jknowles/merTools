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
#
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
#
#
