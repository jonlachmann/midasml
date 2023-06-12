prox_sgl <- function (gstart, gend, nobs, x, r, b, al, gamma, pf, peps, gw, step) {
    big <- 9.9E30
    s <- step
    vg_const <- s * gw * al * (1 - gamma)
    #! -------- BEGIN PROGRAM -------- #!
    aaa <- 1
    while (TRUE) {
        bold <- b
        #!--------- LASSO PART ----------#!
        u <- b + s * t(x[, gstart:gend, drop = FALSE]) %*% r / nobs
        v <- abs(u) - s * al * gamma * pf[gstart:gend]
        v <- pmax2(v, 0)
        b <- v * sign(u)

        #!--------- g-LASSO PART ----------#!
        #! L2 norm of b_g
        normg <- sqrt(sum(b * b))
        if (normg != 0) {
            vg <- vg_const / normg
        } else {
            vg <- big
        }

        scl <- pmax2(1 - pf[gstart:gend] * vg, 0)
        b <- b * scl
        d <- b - bold
        r <- r - x[, gstart:gend] %*% d
        maxg <- max(abs(d))

        #!--------- CHECK CONVERGENCE ----------#!
        if (maxg < peps) break
        aaa <- aaa + 1
    }
    cat("\n", aaa)
  return(list(b = b, r = r))
}

pmax2 <- function(k,x) (x+k + abs(x-k))/2
