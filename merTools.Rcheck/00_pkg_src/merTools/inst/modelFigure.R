# # Not working currently
#
# # Diagramming a model
# library(DiagrammeR)
#
#
# # Get the 'nycflights13' package if not already installed
# # install.packages('nycflights13')
#
# # Get the 'lubridate' package if not already installed
# # install.packages('lubridate')
#
# # Get the latest build of the 'DiagrammeR' package from GitHub
# devtools::install_github('rich-iannone/DiagrammeR')
#
# library("nycflights13")
# library("lubridate")
# library("DiagrammeR")
# library("pipeR")
#
# # Choose a day from 2013 for NYC flight data
# # (You can choose any Julian day, it's interesting to see results for different days)
# day_of_year <- 10
#
# # Get a data frame of complete cases (e.g., flights have departure and arrival times)
# nycflights13 <-
#   nycflights13::flights[which(complete.cases(nycflights13::flights) == TRUE), ]
#
# # Generate a POSIXct vector of dates using the 'ISOdatetime' function
# # Columns 1, 2, and 3 are year, month, and day columns
# # Column 4 is a 4-digit combination of hours (00-23) and minutes (00-59)
# date_time <-
#   data.frame("date_time" =
#                ISOdatetime(year = nycflights13[,1],
#                            month = nycflights13[,2],
#                            day = nycflights13[,3],
#                            hour = gsub("[0-9][0-9]$", "", nycflights13[,4]),
#                            min = gsub(".*([0-9][0-9])$", "\\1", nycflights13[,4]),
#                            sec = 0, tz = "GMT"))
#
# # Add the POSIXct vector 'date_time' to the 'nycflights13' data frame
# nycflights13 <- cbind(date_time, nycflights13)
#
# # Select flights only from the specified day of the year 2013
# nycflights13_day <-
#   subset(nycflights13,
#          date_time >= ymd('2013-01-01', tz = "GMT") + days(day_of_year - 1) &
#            date_time < ymd('2013-01-01', tz = "GMT") + days(day_of_year))
#
# # Create the 'nodes' data frame where at least one column is named "nodes" or "node_id"
# # Column 12 is the 3-letter code for the airport departing from
# # Column 13 is for the airport arriving to
# # (Option: change df to 'nycflights13_day' and only airports used for the day will be included)
# nodes_df <- create_nodes(nodes = unique(c(nycflights13[,12],
#                                           nycflights13[,13])),
#                          label = FALSE)
#
# # The 'edges' data frame must have columns named 'edge_from' and 'edge_to'
# # The color attribute is determined with an 'ifelse' statement, where
# # column 8 is the minutes early (negative values) or minutes late (positive values)
# # for the flight arrival
# edges_df <- create_edges(edge_from = nycflights13_day[,12],
#                          edge_to = nycflights13_day[,13],
#                          color = ifelse(nycflights13_day[,8] < 0,
#                                         "green", "red"))
#
# # Set the graph diagram's default attributes for...
#
# # ...nodes
# node_attrs <- c("style = filled", "fillcolor = lightblue",
#                 "color = gray", "shape = circle", "fontname = Helvetica",
#                 "width = 1")
#
# # ...edges
# edge_attrs <- c("arrowhead = dot")
#
# # ...and the graph itself
# graph_attrs <- c("layout = circo",
#                  "overlap = false",
#                  "fixedsize = true",
#                  "ranksep = 3",
#                  "outputorder = edgesfirst")
#
# # Generate the graph diagram in the RStudio Viewer.
# # The green lines show flights that weren't late (red indicates late arrivals)
# # This graph is for a single day of flights, airports that are unconnected on a
# # given day may be destinations on another day
# create_graph(nodes_df = nodes_df, edges_df = edges_df,
#              graph_attrs = graph_attrs, node_attrs = node_attrs,
#              edge_attrs = edge_attrs, directed = TRUE) %>>%
#   render_graph(width = 1200, height = 800)
