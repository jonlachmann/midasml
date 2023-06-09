! --------------------------------- sglfitpathF --------------------------------- !
sglfitpath_alg <- function (maj, gamma, ngroups, gindex, nobs, nvars, x, y, ju, pf, dfmax,
 pmax, nlam, flmin, ulam, eps, peps, maxit, nalam, b0, beta, m, nbeta, alam,
 npass, jerr, intr) {

  #IMPLICIT NONE
  #! -------- INPUT VARIABLES -------- !
  #INTEGER :: mnl, nobs, nvars, dfmax, pmax, nlam, maxit, nalam, npass, jerr, intr, ngroups
  #INTEGER :: ju(nvars), m(pmax), nbeta(nlam), gindex(ngroups)
  #DOUBLE PRECISION :: eps, gamma, peps, steps(ngroups)
  #DOUBLE PRECISION :: x(nobs, nvars), y(nobs), maj(nvars)
  #DOUBLE PRECISION :: pf(nvars)
  #DOUBLE PRECISION :: beta(pmax, nlam), b0(nlam)
  #DOUBLE PRECISION :: ulam(nlam), alam(nlam)
  #! -------- LOCAL DECLARATIONS -------- !
  #INTEGER,  PARAMETER :: mnlam = 6
  #DOUBLE PRECISION :: d, dif, oldb, u, al, flmin
  #DOUBLE PRECISION,  DIMENSION(:), ALLOCATABLE :: b, oldbeta, r
  #DOUBLE PRECISION :: gw, tmp
  #INTEGER :: gstart, gend
  #INTEGER :: k, j, l, g, vrg, ctr, ierr, ni, me, pln
  #INTEGER :: gs, gj, skip
  #INTEGER,  DIMENSION(:),  ALLOCATABLE :: mm
  #! -------- ALLOCATE VARIABLES -------- !
  #ALLOCATE(b(0:nvars), STAT=jerr)
  #ALLOCATE(oldbeta(0:nvars), STAT=ierr)
  #jerr = jerr + ierr
  #ALLOCATE(mm(1:nvars), STAT=ierr)
  #jerr = jerr + ierr
  #ALLOCATE(r(1:nobs), STAT=ierr)
  #jerr = jerr + ierr
  #if (jerr /= 0) RETURN
  #! ---------------- INITIALIZATION ---------------- !
  b <- 0
  oldbeta <- 0
  m <- 0
  mm <- 0
  npass <- 0
  ni <- npass
  mnl <- min(mnlam, nlam)
  ju <- 0
  r <- y
  for (k in seq_len(ngroups)) {
    gend <- gindex[k]
    if (k == 1) {
        gstart <- 1
    } else {
        gstart <- gindex[k-1] + 1
    }
    tmp <- sum((t(x[, gstart:gend]) %*% x[, gstart:gend]) / nobs * (t(x[, gstart:gend]) %*% x[, gstart:gend]) / nobs)
    steps[k] <- 1 / sqrt(tmp)
  }

  # ! ----------------- LAMBDA LOOP (OUTMOST LOOP) ------------------- !
  for (l in seq_len(nlam)) {
        al <- ulam[l]
        ctr <- 0
        pln <- 0
        # ! ------------------ OUTER LOOP -------------------- !
        while (TRUE) {
            if (intr == 1) oldbeta[0] <- b[0]
            if (ni > 0) oldbeta[m[1:ni]] <- b[m[1:ni]]
            # ! ----------------- MIDDLE LOOP -------------------- !
            while (TRUE) {
                npass <- npass + 1
                dif <- 0
                if (intr == 1) oldbeta[0] <- b[0]
                if (ni > 0) oldbeta[m[1:ni]] <- b[m[1:ni]]
                pln <- pln + 1
                # ! ----------------- GROUP LOOP -------------------- !
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
                    skip <- 1
                    if (pln == 1) {
                        skip <- 0
                    }
                    if (ju[k] == 1) {
                        skip <- 0
                    }
                    if (skip == 0) {
                        # ! --- sg-LASSO PROXIMAL MAP --- !
                        ## TODO: CALL prox_sgl(gstart, gend, nvars, nobs, x, r, b(1:nvars), al, gamma, pf, peps, gw, steps(k))
                        # ! UPDATE REMAINING VARIABLES
                        for (g in gstart:gend) {
                            # !if (abs(b(g))<eps) b(g) = 0
                            if (abs(b[g]) > 0) {
                                if (pln == 1) {
                                    ju[k] <- 1
                                }
                                d <- oldbeta[g] - b[g]
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
                }# ! ----------> END GROUP LOOP
                if (ni > pmax) break
                if (intr == 1) {
                    oldb <- oldbeta[0]
                    u <- 0
                    while (TRUE) { #! BEGIN GRADIENT DESCENT
                      d <- sum(r) / nobs
                      if (d^2 < eps) break
                      b[0] <- b[0] + d
                      r <- r - d
                    }# ! END GRADIENT DESCENT
                    d <- b[0] - oldb
                    if (abs(d) > 0) dif <- max(dif, d^2)
                }
                if (dif < eps) break
            }# ! ----------> END MIDDLE LOOP
            if (ni > pmax) break
            # ! -------------- FINAL CHECK ---------------- !
            vrg <- 1
            if (intr == 1) {
            if ((b[0] - oldbeta[0])^2 >= eps) vrg <- 0
            }
            for (j in seq_len(ni)) {
                if ((b[m[j]] - oldbeta[m[j]])^2 >= eps) {
                    vrg <- 0
                    break
                }
            }
            if (vrg == 1) break
            ctr <- ctr + 1
            if (ctr > maxit) {
             jerr <- - l
             return()
            }
        }# ! -------> END OUTER LOOP
        # ! ----------- FINAL UPDATE & SAVE RESULTS ------------ !
        if (ni > pmax) {
         jerr <- - 10000 - l
         break
        }
        if (ni > 0) beta[1:ni, l] <- b[m[1:ni]]
        nbeta[l] <- ni
        b0[l] <- b(0)
        alam[l] <- al
        nalam <- l
        if (l < mnl) next
        if (flmin >= 1) next
        me <- sum(abs(beta[1:ni, l]) > 0)
        if (me > dfmax) break
  }# ! -------> END LAMBDA LOOP

  return()
}