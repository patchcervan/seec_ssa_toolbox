# Convinienent function to calculate Euclidean distance in a dplyr workflow

calcDist <- function(s, x, y){
    sqrt((s[1] - x)^2 + (s[2] - y)^2)
}
