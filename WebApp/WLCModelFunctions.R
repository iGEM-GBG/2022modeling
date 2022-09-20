
# Model polymer endpoint distances and endpoint connection
onePolymer <- function(some_text) {
  text <- some_text
  
  number <- 4
  resultList <- list("text" = text, "number" = number)
  return (resultList)
}

# Model distance and connection for two polymers
twoPolymers <- function(some_text) {
  text <- some_text
  number <- 8
  resultList <- list("text" = text, "number" = number)
  return (resultList)
}

onePolymerSupport <- function(start) {
  
}

twoPolymersSupport <- function() {
  
}