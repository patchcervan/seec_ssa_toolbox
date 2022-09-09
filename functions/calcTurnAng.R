# Function to calculate turning angles. Adapted from moveHMM
# https://github.com/TheoMichelot/moveHMM/blob/master/R/turnAngle.R

calcTurnAng <- function(s1, s2, s3){

    s1 = as.numeric(s1); s2 = as.numeric(s2); s3 = as.numeric(s3)

    v <- c(s2[1] - s1[1], s2[2] - s1[2])

    if(!any(is.na(s3))){
        w <- c(s3[1] - s2[1], s3[2] - s2[2])
        ang <- atan2(w[2], w[1]) - atan2(v[2], v[1])

        if(ang <= -pi)
            ang <- ang + 2*pi
        if(ang > pi)
            ang <- ang - 2*pi

    } else {
        ang <- 0
    }

    return(ang)

}
