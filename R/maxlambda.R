maxlambda <- function (nvars, nobs, x, y, gamma, gindex, ngroups, pf) {
#    IMPLICIT NONE
#    ! -------- INPUT VARIABLES -------- !
#    INTEGER :: nvars, nobs, ngroups, gindex(ngroups)
#    DOUBLE PRECISION :: x(nobs,nvars), y(nobs), pf(nvars), gamma, maxlam
#    ! -------- LOCAL DECLARATIONS -------- !
#    INTEGER :: k, c, nzvars
#    INTEGER :: gstart, gend, gs, gj
#    DOUBLE PRECISION :: gw, xy(nvars), r(nobs)
#    DOUBLE PRECISION :: wmaxg(ngroups), lb, rb

#    ! -------- BEGIN PROGRAM -------- !
    c <- 0
    wmaxg <- numeric(ngroups)
    r <- y
    nzvars <- 0
    for (k in seq_len(nvars)) {
        if (pf[k] == 0) {
            nzvars <- nzvars + 1
        }
    }
    if (nzvars != 0) {
        # CALL rnz(nvars, nobs, nzvars, y, x, r, pf)
    }
    xy <- (t(x) %*% r) / nobs


    if (gamma == 1) {
        maxlam <- max(abs(xy))
    } else {
        for (k in seq_len(ngroups)) {
            gend <- gindex[k]
            if (k == 1) {
                gstart <- 1
            } else {
                gstart <- gindex[k - 1] + 1
            }
            gs <- gend - gstart + 1
            gw <- 0
            for (gj in gstart:gend) {
                gw <- gw + pf[gj]
            }
            gw <- sqrt(gw)
            if (gw == 0) {
                wmaxg[k] <- 0
            } else {
                if (gamma == 0) {
                    # !rb = NORM2(xy(gstart:gend))
                    rb <- sqrt(sum(xy[gstart:gend] * xy[gstart:gend]))
                    wmaxg[k] <- rb / gw
                } else {
                    lb <- 0
                    rb <- max(abs(xy[gstart:gend]))/gamma
                    rb <- solvewmaxg(gstart, gend, gamma, lb, rb, gw, pf, xy, nvars)
                    wmaxg[k] <- rb
                }
            }
        }
        maxlam <- max(wmaxg)
    }
    # !--- ADD SMALL NUMBER TO ENSURE b = 0 @ maxlam (DUE TO NUMERICAL IMPRESSION)
    maxlam <- maxlam + 1E-5
    return(maxlam)
}