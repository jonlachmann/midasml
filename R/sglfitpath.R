sglfitpath <- function(x, y, nlam, flmin, ulam, isd, intr, nf, eps, peps, dfmax, pmax, jd, 
                    pf, gindex, ngroups,  maxit, gamma, nobs, nvars, vnames) {
    #################################################################################
    # data setup
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    
    # gamma setup
    if (gamma < 0 || gamma > 1) stop("gamma must be in [0,1]")  
    
	# gindex should contain only the variable index at the end of each group - no overlap is allowed 
    if (any(diff(gindex)>1)) 
        stop("only adjacent group memberships are allowed")
    
    gindex <- which(c(diff(gindex),1)==1) 
    #################################################################################
    # call Fortran
    if (nf == 0) {

        ju <- chkvars2(x)

        std <- standard2(x, isd, intr)

        if (max(pf) < 0) {
            stop("ERROR")
        }

        pf <- pmax(pf, 0)

        if (ulam[1] == -1) {
            maxlam <- maxlambda(nvars, nobs, std$x, y, gamma, gindex, ngroups, pf)

            ulam <- numeric(nlam)
            ulam[1] <- maxlam
            for (j in 2:nlam) {
                tmp <- log(maxlam) + (log(maxlam * flmin) - log(maxlam)) * (j - 1) / (nlam - 1)
                ulam[j] <- exp(tmp)
            }
        }

        fit2 <- .Fortran("sglfitF", as.double(gamma), as.integer(ngroups), as.integer(gindex),
                        as.integer(nobs), as.integer(nvars), as.matrix(std$x), y, pf, dfmax, pmax, nlam, flmin, ulam,
                        eps, as.double(peps), intr, maxit, nalam = integer(1), b0 = double(nlam),
                        beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                        alam = double(nlam), npass = integer(1), jerr = integer(1), ju = as.integer(ju), maj = as.double(std$maj),
                        PACKAGE = "midasml")

        fit <- sglfitpath_alg(std$maj, gamma, ngroups, gindex, nobs, nvars, std$x, y, pf, dfmax,
                    pmax, nlam, flmin, ulam, eps, peps, maxit, nalam, nlam, beta, nlam, nlam, intr)

        xmean <- attr(std$x, "scaled:center")
        xnorm <- attr(std$x, "scaled:scale")

        beta <- matrix(fit$beta, ncol = fit$nalam)
        for (l in seq_len(fit$nalam)) {
            nk <- fit$nbeta[l]
            if (isd == 1) {
                 for (j in seq_len(nk)) {
                      beta[j, l] <- beta[j, l] / xnorm[fit$ibeta[j]]
                 }
            }
            if (intr == 1) {
                fit$b0[l] <- fit$b0[l] - sum(beta[seq_len(nk), l] * xmean[fit$ibeta[seq_len(nk)]])
            }
        }
        fit$beta <- as.numeric(beta)

    } else {
        fit <- .Fortran("panelsglfitF", as.integer(nf), as.integer(nobs/nf), as.double(gamma), as.integer(ngroups), as.integer(gindex),
                        as.integer(nobs), as.integer(nvars), as.matrix(x), y, pf, dfmax, pmax, nlam, flmin, ulam, 
                        eps, as.double(peps), isd, intr, maxit, nalam = integer(1), b0 = double(nlam), a0 = double(nf * nlam),
                        beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                        alam = double(nlam), npass = integer(1), jerr = integer(1), 
                        PACKAGE = "midasml")
        
    }
    #################################################################################
    # output
    if (nf == 0){
        nf <- intr
    }
    fit$nf <- nf
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    outlist$dimx <- c(nobs, nvars)
    class(outlist) <- c("sglpath")
    outlist
} 
