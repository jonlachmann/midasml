prox_sgl <- function (gstart, gend, nvars, nobs, x, r, b, al, gamma, pf, peps, gw, step) {
    big <- 9.9E30
    bold <- numeric(nvars)
    s <- step
    vg_const <- s * gw * al * (1 - gamma)
    #! -------- BEGIN PROGRAM -------- #!
    aaa <- 1
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
        if (normg != 0) {
            vg <- vg_const / normg
        } else {
            vg <- big
        }
        #! Initialize storage vars

        if (normg == 0) {
            vg <- big
        }

        scl <- pmax(1 - pf[gstart:gend] * vg, 0)
        b[gstart:gend] <- b[gstart:gend] * scl
        d <- b[gstart:gend] - bold[gstart:gend]
        r <- r - x[, gstart:gend] %*% d
        maxg <- max(abs(d))

        #!--------- CHECK CONVERGENCE ----------#!
        if (maxg < peps) break
        aaa <- aaa + 1
    }
    cat("\n", aaa)
  return(list(b = b, r = r))
}
