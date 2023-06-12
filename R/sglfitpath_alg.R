# ! --------------------------------- sglfitpathF --------------------------------- !
sglfitpath_alg <- function(
    maj, gamma, ngroups, gindex, nobs, nvars, x, y, pf, dfmax,
    pmax, nlam, flmin, ulam, eps, peps, maxit, nalam, b0, beta, nbeta, alam, intr) {
  mnlam <- 6
  steps <- stepSize(x, gindex, ngroups, nobs)
  b <- numeric(nvars + 1)
  oldbeta <- numeric(nvars + 1)
  beta <- matrix(0, pmax, nlam)
  m <- numeric(pmax)
  mm <- numeric(nvars)
  npass <- 0
  ni <- npass
  mnl <- min(mnlam, nlam)
  ju <- numeric(nvars)
  r <- y

  # ! ----------------- LAMBDA LOOP (OUTMOST LOOP) ------------------- !
  for (l in seq_len(nlam)) {
    al <- ulam[l]
    ctr <- 0
    pln <- 0
    # ! ------------------ OUTER LOOP -------------------- !
    while (TRUE) {
      if (intr == 1) oldbeta[1] <- b[1]
      if (ni > 0) oldbeta[m[1:ni] + 1] <- b[m[1:ni] + 1]
      # ! ----------------- MIDDLE LOOP -------------------- !
      while (TRUE) {
        npass <- npass + 1
        dif <- 0
        if (intr == 1) oldbeta[1] <- b[1]
        if (ni > 0) oldbeta[m[1:ni] + 1] <- b[m[1:ni] + 1]
        pln <- pln + 1
        # ! ----------------- GROUP LOOP -------------------- !
        for (k in seq_len(ngroups)) {
          gend <- gindex[k]
          if (k == 1) {
            gstart <- 1
          } else {
            gstart <- gindex[k - 1] + 1
          }
          gw <- 0
          for (gj in gstart:gend) {
            gw <- gw + pf[gj]
          }
          gw <- sqrt(gw)
          skip <- 1
          if (pln == 1) {
            skip <- 0
          }
          if (ju[k] == 1) {
            skip <- 0
          }
          if (skip == 0) {
            # ! --- sg-LASSO PROXIMAL MAP --- !
            prox_res <- prox_sgl(gstart, gend, nvars, nobs, x, r, b[1:nvars + 1], al, gamma, pf, peps, gw, steps[k])
            r <- prox_res$r
            b[1:nvars + 1] <- prox_res$b
            # ! UPDATE REMAINING VARIABLES
            for (g in gstart:gend) {
              if (abs(b[g + 1]) > 0) {
                if (pln == 1) {
                  ju[k] <- 1
                }
                d <- oldbeta[g + 1] - b[g + 1]
                dif <- max(dif, maj[g] * d^2)
                if (mm[g] == 0) {
                  ni <- ni + 1
                  if (ni > pmax) break
                  mm[g] <- ni
                  m[ni] <- g # ! RECORD ACTIVE VARIABLES
                }
              }
            }
          } # ! ----------> END GROUP UPDATES
        } # ! ----------> END GROUP LOOP
        if (ni > pmax) break
        if (intr == 1) {
          oldb <- oldbeta[1]
          u <- 0
          while (TRUE) { # ! BEGIN GRADIENT DESCENT
            d <- sum(r) / nobs
            if (d^2 < eps) break
            b[1] <- b[1] + d
            r <- r - d
          } # ! END GRADIENT DESCENT
          d <- b[1] - oldb
          if (abs(d) > 0) dif <- max(dif, d^2)
        }
        if (dif < eps) break
      } # ! ----------> END MIDDLE LOOP
      if (ni > pmax) break
      # ! -------------- FINAL CHECK ---------------- !
      vrg <- 1
      if (intr == 1) {
        if ((b[1] - oldbeta[1])^2 >= eps) vrg <- 0
      }
      for (j in seq_len(ni)) {
        if ((b[m[j] + 1] - oldbeta[m[j] + 1])^2 >= eps) {
          vrg <- 0
          break
        }
      }
      if (vrg == 1) break
      ctr <- ctr + 1
      if (ctr > maxit) {
        return()
      }
    } # ! -------> END OUTER LOOP
    # ! ----------- FINAL UPDATE & SAVE RESULTS ------------ !
    if (ni > pmax) {
      break
    }
    if (ni > 0) beta[seq_len(ni), l] <- b[m[seq_len(ni)] + 1]
    nbeta[l] <- ni
    b0[l] <- b[1]
    alam[l] <- al
    nalam <- l
    if (l < mnl) next
    if (flmin >= 1) next
    me <- sum(abs(beta[1:ni, l]) > 0)
    if (me > dfmax) break
  } # ! -------> END LAMBDA LOOP

  return(list(
    alam = alam,
    beta = beta,
    b0 = b0,
    nalam = nalam,
    nbeta = nbeta,
    ibeta = m,
    jerr = 0,
    npass = npass
  ))
}

stepSize <- function (x, gindex, ngroups, nobs) {
  steps <- numeric(ngroups)
  for (k in seq_len(ngroups)) {
    gend <- gindex[k]
    if (k == 1) {
      gstart <- 1
    } else {
      gstart <- gindex[k - 1] + 1
    }
    xtx <- t(x[, gstart:gend]) %*% x[, gstart:gend]
    steps[k] <- 1 / sqrt(sum(xtx / nobs * xtx / nobs))
  }
  return(steps)
}
