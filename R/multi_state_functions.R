#' sinlink
#' @description link transformation functions - transform beta's into real parameters
#' @param x
#' @references Miller, D. A. W., C. S. Brehem, J. E. Hines, J. D. Nichols, and R. N. Fisher. 2012. Joint estimation of habitat dynamics and species interactions: disturbance reduces co-occurrence of non-native predators with an endangered toad. Journal of Applied Ecology 81:1288-1297.

#' @export

sinlink = function(x) {
	return(sin(x)/2 + 0.5)
}

#' logitlink
#' @description link transformation functions - transform beta's into real parameters
#' @param x
#' @references Miller, D. A. W., C. S. Brehem, J. E. Hines, J. D. Nichols, and R. N. Fisher. 2012. Joint estimation of habitat dynamics and species interactions: disturbance reduces co-occurrence of non-native predators with an endangered toad. Journal of Applied Ecology 81:1288-1297.

#' @export

logitlink = function(x) {
	return(1/(1 + exp(-x)))
}

#' linkfn
#' @description link transformation functions - transform beta's into real parameters
#' @param x
#' @references Miller, D. A. W., C. S. Brehem, J. E. Hines, J. D. Nichols, and R. N. Fisher. 2012. Joint estimation of habitat dynamics and species interactions: disturbance reduces co-occurrence of non-native predators with an endangered toad. Journal of Applied Ecology 81:1288-1297.

#' @export

linkfn = function(x) {
	return(c(
		sin(x[1:76])/2 + 0.5,
		1/(1 + exp(-x[77:104])),
		sin(x[105:180])/2 + 0.5,
		1/(1 + exp(-x[181:208]))
	))
}

#' Computes initial state vector
#'
#' @param beta
#' @param site_type
#' @param np2

#' @references Miller, D. A. W., C. S. Brehem, J. E. Hines, J. D. Nichols, and R. N. Fisher. 2012. Joint estimation of habitat dynamics and species interactions: disturbance reduces co-occurrence of non-native predators with an endangered toad. Journal of Applied Ecology 81:1288-1297.

#' @export

initPsi = function(beta, site_type) {
    k = site_type * np2
    w1 = beta[k + 1] # Pr(hab: wet)
    n1 = beta[k + 2] # Pr(occ nonnatives)
    t1n = beta[k + 3] # Pr(occ toads| not occ non-natives)
    t1N = beta[k + 4] # Pr(occ non-natives| occ non-natives)
    # so, init state vector = Pr(state=1),Pr(state=2)
    return(c(
		1 - w1, # prob that its not wet?
		w1 * (1 - n1) * (1 - t1n), # wet, without natives and without toads
		w1 * (1 - n1) * t1n, # wet, without natives but with toads
		w1 * n1 * (1 - t1N), # wet with natives and without non-natives
		w1 * n1 * t1N # wet with non-natives
	))
}

#' Computes transition matrix
#'
#' @param beta
#' @param yr
#' @param site_type
#' @param np2

#' @references Miller, D. A. W., C. S. Brehem, J. E. Hines, J. D. Nichols, and R. N. Fisher. 2012. Joint estimation of habitat dynamics and species interactions: disturbance reduces co-occurrence of non-native predators with an endangered toad. Journal of Applied Ecology 81:1288-1297.

#' @export

phi = function(beta, yr, site_type, np2) {
    k = 4 + site_type * np2  # (eg., k=4, nyrsM1=6, yr=1)
    #  beta[5] = Pr(hab: dry-dry)
	w0 = beta[k + yr]
    k = k + nyrsM1  
    #  beta[11]= Pr(hab: wet-wet)
	w1 = beta[k + yr]
    k = k + nyrsM1  
    #  beta[17]= Pr(non-natives: unocc-unocc|dry in i-1)
	n0d = beta[k + yr]
    k = k + nyrsM1  
    #  beta[23]= Pr(non-natives: unocc-unocc|unocc toads in i-1)
	n0 = beta[k + yr]
    k = k + nyrsM1  
    #  beta[29]= Pr(non-natives: unocc-unocc|occ toads in i-1)
	n0a = beta[k + yr]
    k = k + nyrsM1  
    #  beta[35]= Pr(non-natives: occ-occ|unocc toads in i-1)
	n1 = beta[k + yr]
    k = k + nyrsM1  
    #  beta[41]= Pr(non-natives: occ-occ|occ toads in i-1)
	n1a = beta[k + yr]
    k = k + nyrsM1  
    #  beta[47]= Pr(toads: unocc-unocc| dry in i-1)
	t0d = beta[k + yr]
    k = k + nyrsM1  
    #  beta[53]= Pr(toads: unocc-unocc| unocc non-natives in i-1)
	t0 = beta[k + yr]
    k = k + nyrsM1  
	#  beta[59]= Pr(toads: unocc-unocc| occ non-natives in i-1)
	t0a = beta[k + yr]
    k = k + nyrsM1  
    #  beta[65]= Pr(toads: occ-occ| unocc non-natives in i-1)
	t1 = beta[k + yr]
    k = k + nyrsM1
   #  beta[71]= Pr(toads: occ-occ| occ non-natives in i-1)
   t1a = beta[k + yr]
    
    # here's the transition matrix - Pr(system goes from state i to state j) where i=row, j=col
    x1 = c(
		w0,
		(1 - w0) * t0d * n0d,
		(1 - w0) * (1 - t0d) * n0d,
		(1 - w0) * t0d * (1 - n0d),
		(1 - w0) * (1 - t0d) * (1 - n0d)
	)
    x2 = c(
		1 - w1,
		w1 * t0 * n0,
		w1 * (1 - t0) * n0,
		w1 * t0 * (1 - n0),
		w1 * (1 - t0) * (1 - n0)
	)
    x3 = c(
		1 - w1,
		w1 * (1 - t1) * n0a,
		w1 * t1 * n0a,
		w1 * (1 - t1) * (1 - n0a),
		w1 * t1 * (1 - n0a)
	)
    x4 = c(
		1 - w1,
		w1 * t0a * (1 - n1),
		w1 * (1 - t0a) * (1 - n1),
		w1 * t0a * n1,
		w1 * (1 - t0a) * n1
	)
    x5 = c(
		1 - w1,
		w1 * (1 - t1a) * (1 - n1a),
		w1 * t1a * (1 - n1a),
		w1 * (1 - t1a) * n1a,
		w1 * t1a * n1a
	)
    return(rbind(x1, x2, x3, x4, x5))
}

#' Computes detection probability matrix
#' @description Pr(row,col) = Pr(det in state=col | true state=row)
#' @param beta
#' @param yr
#' @param detState
#' @param site_type

#' @references Miller, D. A. W., C. S. Brehem, J. E. Hines, J. D. Nichols, and R. N. Fisher. 2012. Joint estimation of habitat dynamics and species interactions: disturbance reduces co-occurrence of non-native predators with an endangered toad. Journal of Applied Ecology 81:1288-1297.

#' @export
	
pq = function(beta, yr, detState, site_type) {
    # if no detection, return ...
    x = id5
    if (detState < 1) {
        # q = 1-p = Pr(not detected, given current state) = ID matrix
		return(id5)
	}
    k = (site_type + 1) * np2 - 4 * nyrs
    #  Pr(detect non-natives | not detect toads      )
	pNt = beta[k + yr]
    k = k + nyrs
    #  Pr(detect non-natives | detect toads          )
	pNT = beta[k + yr]
    k = k + nyrs  
    #  Pr(detect toads       | not detect non-natives)
	pTn = beta[k + yr]
    k = k + nyrs  
    #  Pr(detect toads       | detect non-natives    )
	pTN = beta[k + yr]  
    x[1, ] = c(1, 0, 0, 0, 0)
    x[2, ] = c(0, 1, 0, 0, 0)
    x[3, ] = c(
		0,
		1 - pTn,
		pTn,
		0,
		0
	)
    x[4, ] = c(
		0,
		1 - pNt,
		0,
		pNt,
		0
	)
    x[5, ] = c(
		0,
		(1 - pTn) * (1 - pNt),
		pTn * (1 - pNT),
		(1 - pTN) * pNt,
		pTN * pNT
	)
	#   use diagonal of above matrix
    return(diag(x[, detState]))
}


#' Computes detection probability matrix
#' @description Here's the Likelihood function for the Multi-state, integrated habitat model. It's suprisingly simple, using matrix multiplication.  It's done as a standard Multi-strata capture-recapture model. This function computes likelihood value for given betas.
#' @param beta

#' @references Miller, D. A. W., C. S. Brehem, J. E. Hines, J. D. Nichols, and R. N. Fisher. 2012. Joint estimation of habitat dynamics and species interactions: disturbance reduces co-occurrence of non-native predators with an endangered toad. Journal of Applied Ecology 81:1288-1297.

#' @export

Like = function(beta) {
    # convert vector of beta's to real parameters
	realparms = linkfn(dm %*% beta)  
    xll = 0
    for (i in 1:nsites) {
        # for each site... compute detection history prob.
        prb = initPsi(realparms, stype[i])
        # history prob = init. psi * PROD{(p[j] or q[j]) * phi[j]}
		iyr = 1  
        for (j in 1:nocc) {
            # for j=1 to number of survey occasions
            prb = prb %*% pq(realparms, iyr, h[i, j], stype[i])
            if (intrvl[j] > 0) {
                # if interval between surveys > 0
                xphi = phi(realparms, iyr, stype[i])  #   then transitions can occur
                prb = prb %*% xphi
            }
            iyr = iyr + intrvl[j]
        }
        xll = xll + log(sum(prb))  #  log-likelihood = sum(log(cell-prob)) for each history
        
    }
    niter <<- niter + 1
    cat(c(niter, " xll=", xll, "\n"))
    return(-xll)
}

#' Computes detection probability matrix
#' @description This function runs a model specified by name, where name is of the form psi(X), phi(Y), p(Z), where:
#' \itemize{
#' 	\item{X is a design-matrix file for initial occupancy (psi), and specified by by the variable, modname1}
#' 	\item{Y is a design-matrix file for state transitions (phi), and specified by by the variable, modname2}
#' 	\item{Z is a design-matrix file for detection probabilities (p), and specified by by the variable, modname3}
#' }
# The model to run is defined by a single design-matrix.  This design matrix consists of R rows and C columns, where 
#' \itemize{
	#' \item{R = number of real parameters (\psi, \phi, p) and}
	#' \item{C = number of \beta's}
#' }
#'	To make input simpler, the design matrix is split into 3 parts:
#' \itemize{
	#' \item{Part 1 contains all rows and the columns which affect initial \psi.}
	#' \item {Part 2 contains all rows and the columns which affect transitions (e's and g's).}
	#' \item{Part 3 contains all rows and the columns which affect detection (p's).}
#' }
#' The 3 parts are then horizontally concatenated to produce a single design matrix.
#' @param beta

#' @references Miller, D. A. W., C. S. Brehem, J. E. Hines, J. D. Nichols, and R. N. Fisher. 2012. Joint estimation of habitat dynamics and species interactions: disturbance reduces co-occurrence of non-native predators with an endangered toad. Journal of Applied Ecology 81:1288-1297.

#' @export

occHabDyn = function(mnum, modname1, modname2, modname3) {
    modname = sprintf("%d)%s%s%s", mnum, modname1, modname2, modname3)
    # read in design-matrix, part 1
    dmname = sprintf("DM%s.csv", modname1)
    dmx = read.table(dmname, sep = ",", na = "-", header = FALSE, colClasses = c("character", "numeric"))
    dm1 = data.matrix(dmx[, 2:dim(dmx)[2]])
    lbls = dmx[, 1]
    npar <<- length(lbls)
    np2 <<- npar/2
    
    
    # read in design-matrix, part 2
    
    dmname = sprintf("DM%s.csv", modname2)
    dmx = read.table(dmname, sep = ",", na = "-", header = FALSE, colClasses = c("character", "numeric"))
    dm2 = data.matrix(dmx[, 2:dim(dmx)[2]])
    
    # read in design-matrix, part 3
    
    dmname = sprintf("DM%s.csv", modname3)
    dmx = read.table(dmname, sep = ",", na = "-", header = FALSE, colClasses = c("character", "numeric"))
    dm3 = data.matrix(dmx[, 2:dim(dmx)[2]])
    
    # join 3 parts into one big design matrix...
    dm <<- cbind(dm1, dm2, dm3)
    nest = dim(dm)[2]  #  nest = number of estimated beta's ( = number of cols in design matrix)
    b0 = rep(0, nest)  #  b0 = initial values for beta's
    
    # create a file to containe output...
    outname = sprintf("%s.out", modname)
    unlink(outname)
    zz <- file(outname, open = "wt")
    mname <<- modname
    
    # compute log-likelihood value for the initial values of beta...
    xll = Like(b0)
    x = b0
    ans = list(value = xll, minimum = xll, par = b0, estimate = b0)
    cat(c("init values:", linkfn(dm %*% b0), "\n"))
    cat(c("log-like(init)=", xll, "\n"))
    
    # call R optimization (nlm non-linear-model) routine to get maximum-likelihood estimates
    niter <<- 0
    tm1 = proc.time()[1]  #  set initial timer to see how long the optimization takes
    # do optimization using likelihood function, 'Like' and init values, b0
    ans = nlm(Like, b0, hessian = TRUE)
    x = ans$estimate
    xll = ans$minimum
    tm2 = proc.time()[1]
    sink(zz, type = c("output", "message"), append = F, split = T)
    cat(c(hdr, "\n*****************\n      ", name, "\n"))
    cat(c("intrvl:", intrvl, "\n"))
    cat(c("\n ************** ", modname, " *********\n\ndesign matrix:\n"))
    cat(c("npar=", npar, " np2=", np2, " nest=", nest, "\n"))
    print(cbind(lbls, dm), quote = FALSE)
    cat(c("init values:", linkfn(dm %*% b0), "\n"))
    cat(c("log-like(init)=", xll, "\n"))
    cat(c("optimization time:", tm2 - tm1, "\n"))
    cat(c("log-lik=", xll, "\n"))
    cat(c("beta   =", x, "\n"))
    
    # compute variance-covariance matrix of beta's using the hessian result of optimization routine
    vcb = try(solve(ans$hessian))
    if (length(vcb) < 2) 
        vcb = diag(rep(0, nest))
    # compute std. errors of beta's
    cat("\n    beta estimates and variances:\n")
    xse = sqrt(diag(vcb))
    c1 = x - 1.96 * xse
    c2 = x + 1.96 * xse
    cat("\n  estimate  std.err  95% conf. interval")
    cat(sprintf("\n%f %f %f - %f", x, xse, c1, c2))
    cat("\n")
    
    # Next, compute var-cov matrix of real parameters using the 'Delta' method VCR = der(RP wrt BP) * VCB *
    # der(RP wrt BP)' where VCR = var-cov matrix of real params der(RP wrt BP) = derivative of real params
    # with respect to beta params VCB = var-cov matrix of beta params
    
    rp = linkfn(dm %*% x)  #  get real parms from design matrix * final beta's
    
    # compute derivatives of real params with respect to beta's...
    deri = rep(0, npar * nest)
    dim(deri) = c(npar, nest)
    for (i in 1:npar) {
        t1 = linkfn(dm %*% x)
        for (j in 1:nest) {
            x[j] = x[j] + eps
            t2 = linkfn(dm %*% x)
            deri[i, j] = (t2[i] - t1[i])/eps
            x[j] = x[j] - eps
        }
    }
    
    vcr = deri %*% vcb %*% t(deri)  #   compute var-cov of real params
    rpse = sqrt(diag(vcr))  #   compute real param std. errors, rpse
    c1 = rp - 1.96 * rpse
    c2 = rp + 1.96 * rpse  #  compute confidence interval limits, c1,c2
    cat("\n    real estimates and variances:\n")
    cat("\n  parameter  estimate  std.err  95% conf. interval")
    cat(sprintf("\n%12s %f %f %f - %f", lbls, rp, rpse, c1, c2))
    cat("\n")
    sink()
    return(ans)
}
