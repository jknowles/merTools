#Helpers

# Function to take only rows that form distinct levels of factors
trimModelFrame <- function(data){
  # Identify numerics
  nums <- sapply(data, is.numeric)
  vars <- names(nums[!nums == TRUE])
  dataList <- vector(mode = "list", length = length(vars))
  names(dataList) <- vars
    for(i in vars){
      dataList[[i]] <- data[!duplicated(data[, i]),]
    }
    newdat <- do.call(rbind, dataList)
    newdat <- newdat[!duplicated(newdat),]
    return(newdat)
}
