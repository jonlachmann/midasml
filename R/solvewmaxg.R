solvewmaxg <- function (gstart, gend, gamma, lb, rb, gw, pf, xy, nvars) {
    #IMPLICIT NONE
    #! -------- INPUT VARIABLES -------- !
    #INTEGER :: gstart, gend, nvars
    #DOUBLE PRECISION :: gamma, lb, rb, gw, pf(nvars), xy(nvars)
    #! -------- LOCAL DECLARATIONS -------- !
    #INTEGER :: stopflag, indexi, iters
    #DOUBLE PRECISION ::  tol = 1E-13, mp, fl, fm, fr, tmpl, tmpm, tmpr

    tol <- 1E-13
    stopflag <- 0
    iters <- 0
    while (TRUE) {
        mp <- 0.5 * (lb + rb)
        fl <- 0
        fm <- 0
        fr <- 0
        tmpl <- 0
        tmpm <- 0
        tmpr <- 0
        for (indexi in gstart:gend) {
            tmpl <- abs(xy[indexi]) - gamma * lb * pf[indexi]
            tmpm <- abs(xy[indexi]) - gamma * mp * pf[indexi]
            tmpr <- abs(xy[indexi]) - gamma * rb * pf[indexi]
            if (tmpl > 0) {
                fl <- fl + tmpl * tmpl
            }
            if (tmpm > 0) {
                fm <- fm + tmpm * tmpm
            }
            if (tmpr > 0) {
                fr <- fr + tmpr * tmpr
            }
        }
        fl <- fl - (1 - gamma) * (1 - gamma) * lb * lb * gw * gw
        fm <- fm - (1 - gamma) * (1 - gamma) * mp * mp * gw * gw
        fr <- fr - (1 - gamma) * (1 - gamma) * rb * rb * gw * gw
        if (fl * fm < 0) {
            if (abs(lb - mp) > tol) {
                rb <- mp
            } else {
                stopflag <- 1
            }
        } else {
            if (fm * fr < 0) {
                if (abs(mp - rb) > tol) {
                    lb <- mp
                } else {
                    stopflag <- 1
                }
            } else {
                stopflag <- 1
            }
        }
       if (stopflag == 1) break
       if (iters > 1000) break
       iters <- iters + 1
    }
    rb <- mp
  return(rb)
}