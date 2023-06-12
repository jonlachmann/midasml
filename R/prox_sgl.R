prox_sgl <- function (gstart, gend, nvars, nobs, x, r, b, al, gamma, pf, peps, gw, step) {
    big <- 9.9E30
    bold <- numeric(nvars)
    s <- step
    #! -------- BEGIN PROGRAM -------- #!
    while (TRUE) {
        bold[gstart:gend] <- b[gstart:gend]
        #!--------- LASSO PART ----------#!
        for (g in gstart:gend) {
            u <- b[g] + s * sum(x[, g] * r) / nobs
            v <- abs(u) - s * al * gamma * pf[g]
            if (v > 0) {
                tmp <- v * sign(u)
            } else {
                tmp <- 0
            }
            b[g] <- tmp
        }
        #!--------- g-LASSO PART ----------#!
        #! L2 norm of b_g
        normg <- sqrt(sum(b[gstart:gend] * b[gstart:gend]))
        #! Initialize storage vars
        maxg <- 0
        vg <- s * gw * al * (1 - gamma) / normg
        if (normg == 0) {
            vg <- big
        }

        for (g in gstart:gend) {
            scl <- 1 - pf[g] * vg
            scl <- max(scl, 0)
            #!l_2,1 norm map
            tmp <- scl * b[g]
            d <- tmp - bold[g]
            r <- r - x[,g] * d
            maxg <- max(maxg, abs(d))
            b[g] <- tmp
        }
        #!--------- CHECK CONVERGENCE ----------#!
        if (maxg < peps) break
    }
  return(list(b = b, r = r))
}
