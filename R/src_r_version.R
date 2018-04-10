

## {{{ constants }}}

## initial values and constants
dphGLM_get_constants_xxr      <- function(family, d, dj=NULL)
{
    if(family=='gaussian') {        
        sigma_beta = 10
        fix = list(alpha      = 1,
                   mu_beta    = rep(0, d+1),
                   Sigma_beta = sigma_beta * diag(d+1),
                   s2_sigma   = 10,
                   df_sigma   = 10)
    }
    if(family=='binomial') {        
        sigma_beta = 10
        fix = list(alpha      = 1,
                   mu_beta    = rep(0, d+1),
                   Sigma_beta = sigma_beta * diag(d+1))
    }

    return(fix)
}
dpGLM_get_inits_xxr <- function(K, d, fix, family)
{
    beta  = MASS::mvrnorm(n=1, mu=fix$mu_beta, Sigma = fix$Sigma_beta)
    if(family=='gaussian') {        
        theta = matrix(nrow=0,ncol=1+(d+1)+1)  ## ! fo k, (d+1) for covars plus intercept, 1 for sigma
        for (k in 1:K){
            sigma = sqrt(LaplacesDemon::rinvchisq(1, df=fix$df_sigma, scale=fix$s2_sigma))
            theta = rbind(theta, c(k=k, beta=as.vector(beta), sigma=sigma))
        }
    }
    if(family=='binomial') {        
        theta = matrix(nrow=0,ncol=1+(d+1))  ## ! fo k, (d+1) for covars plus intercept, 1 for sigma
        for (k in 1:K){
            theta = rbind(theta, c(k=k, beta=as.vector(beta)))
        }
    }
    return(theta)
}
dphGLM_check_constants_xxr <- function(family=family, d=d, dj=NULL, fix)
{
    ## check type of 'fix'
    ## ------------------
    if (class(fix) != 'list')
    {
        stop("The parameter \'fix\" must be a list")
    }
    ## check names of the constants:
    ## -----------------------------
    if (family=='guassian')
    {
        if (any(! names(fix) %in% c('mu_beta', 'Sigma_beta', 'alpha', 's2_sigma', 'df_sigma')))
        {
            stop("\n\nThe parameter \'fix\' must be a named list with the parameters of the model. The parametters are 'mu_beta', 'Sigma_beta', 'alpha', 's2_sigma', 'df_sigma' \n\n")
        }
    }
    if (family=='binomial')
    {
        if (any(! names(fix) %in% c('mu_beta', 'Sigma_beta', 'alpha')))
        {
            stop("\n\nThe parameter \'fix\' must be a named list with the parameters of the model. The parametters are 'mu_beta', 'Sigma_beta', 'alpha'\n\n")
        }
    }

    ## check type of constants
    ## ----------------
    if (!all(is.vector(fix$mu_beta), is.vector(fix$alpha)))
    {
        stop("\n\n fix$mu_theta and fix$alpha must be vectors\n\n")
    }
    if (family=='gaussian')
    {
        if (!all(is.vector(fix$s2_sigma), is.vector(fix$df_sigma)))
        {
            stop("\n\n fix$s2_sigma and fix$df_sigma must be vectors\n\n")
        }
    } 
    if (!is.matrix(fix$Sigma_beta))
    {
        stop("\n\n fix$Sigma_beta must be a matrix\n\n")
    }

    ## check the dimensions
    ## --------------------
    if (length(fix$mu_beta)!=d+1)
    {
        stop(paste0("Dimension of fix$mu_beta must be the same as the number of covariates plus one (for the intercept), that is", d+1, sep=''))
    }
    if (!all(dim(fix$Sigma_beta)==c(d+1,d+1)) )
    {
        stop(paste0("The dimension of fix$Sigma_beta must be ", d+1,' x ', d+1,  sep=''))
    }
}

## }}}
## {{{ update theta }}}

## parameters and latent variables
dpGLM_update_theta_xxr <- function(y, X, Z, K, theta,  fix, family, epsilon, leapFrog, hmc.iter)
{
    Zstar = unique(Z)

    ## gaussian
    ## --------
    if (family=='gaussian')
    {
        for (k in Zstar)
        {
            mh.accept.info$trial <<- mh.accept.info$trial + 1
            mh.accept.info$accepted <<- mh.accept.info$accepted + 1

            Xk     = matrix(X[Z==k,], nrow=sum(Z==k), ncol=ncol(X))
            yk     = y[Z==k]
            Nk     = sum(Z == k)
            ## get index of cluster k and coluns of parameters in matrix theta to be updated
            k.idx     = which(theta[,'k']==k)
            sigma.idx = which(grepl(colnames(theta), pattern='sigma'))
            beta.idx  = which(grepl(colnames(theta), pattern='beta'))
        
            ## update beta
            ## -----------
            theta[k.idx, beta.idx]  = dphGLM_update_beta_gaussian_xxr(yk=yk ,Xk=Xk, sigmak= theta[k.idx, sigma.idx], fix)
            ## update sigma
            ## ------------
            theta[k.idx, sigma.idx] = dphGLM_update_sigma_gaussian_xxr(yk=yk ,Xk=Xk, Nk=Nk, betak= theta[k.idx, beta.idx], fix)


        }
        Zstar_complement = setdiff(1:K, Zstar)
        for (k in Zstar_complement)
        {
            k.idx     = which(theta[,'k']==k)
            beta.idx  = which(grepl(colnames(theta), pattern='beta'))
            sigma.idx = which(grepl(colnames(theta), pattern='sigma'))

            ## update beta
            ## -----------
            ## get index of cluster k and coluns of parameters in matrix theta to be updated
            betaNew    = MASS::mvrnorm(1, fix$mu_beta, fix$Sigma_beta)
            theta[k.idx, beta.idx ]  = betaNew
            
            ## update sigma
            ## ------------
            df       = fix$df_sigma
            s2       = fix$s2_sigma
            sigmaNew = LaplacesDemon::rinvchisq(1, df=df, scale=s2)
            theta[k.idx, sigma.idx]  = sigmaNew
        }
    }


    ## binomial
    ## --------
    if (family=='binomial'){
        for (k in Zstar){
            Xk     = matrix(X[Z==k,], nrow=sum(Z==k), ncol=ncol(X))
            yk     = y[Z==k]
            Nk     = sum(Z == k)
            ## get index of cluster k and coluns of parameters in matrix theta to be updated
            k.idx     = which(theta[,'k']==k)
            beta.idx  = which(grepl(colnames(theta), pattern='beta'))
        
            ## update beta
            ## -----------
            theta_k_t   = theta[k.idx, beta.idx]
            theta_k_new = dphGLM_hmc_update_binomial_xxr(yk=yk ,Xk=Xk, betak= theta[k.idx, beta.idx], fix=fix, epsilon, leapFrog, hmc.iter)
            theta[k.idx, beta.idx] = theta_k_new

            ## update acceptance rate
            mh.accept.info$trial <<- mh.accept.info$trial + 1
            if (any(theta_k_t != theta_k_new)) {mh.accept.info$accepted <<- mh.accept.info$accepted + 1}
        }
        Zstar_complement = setdiff(1:K, Zstar)
        for (k in Zstar_complement){
            k.idx     = which(theta[,'k']==k)
            beta.idx  = which(grepl(colnames(theta), pattern='beta'))

            ## update beta
            ## -----------
            ## get index of cluster k and coluns of parameters in matrix theta to be updated
            betaNew    = MASS::mvrnorm(1, fix$mu_beta, fix$Sigma_beta)
            theta[k.idx, beta.idx ]  = betaNew
        }
    }

    return(theta)
}

## jags
dphGLM_jags_xxr <- function(model=NA, n.chains=4, family)
{

    model = paste0('dphGLM_', model,'.jag', sep='')
    ## start parallel processing
    ## -------------------------
    cl.chains    = makePSOCKcluster(min(4,n.chains))
    name         = 'tmp'
    sim          = parJagsModel(cl=cl.chains, file=model, data = data, name=name, inits=mcmc$inits, n.chains=mcmc$n.chains, n.adapt=mcmc$n.adapt)
    parUpdate(cl=cl.chains, object=name, n.iter=mcmc$burn.in)
    samples      = parCodaSamples(cl=cl.chains, model=name, variable.names=parameters, n.iter=mcmc$n.iter, thin=1)
    ## stop parallelization
    ## --------------------
    stopCluster(cl.chains)

    #### make here the output to match the one in mcmc_gaussian/mcmc_logistic, etc
    return(samples)

}

## }}}
## {{{ Gibbs : betas (gaussian family) }}}

## gaussian
dphGLM_update_beta_gaussian_xxr   <- function(yk,Xk, sigmak, fix){
    Sigma_beta  = fix$Sigma_beta
    Sk          = solve(Sigma_beta + t(Xk)%*%Xk)
    mu_betak    = Sk %*% t(Xk) %*% yk
    Sigma_betak = Sk * sigmak^2
    betaNew     = MASS::mvrnorm(1, mu = mu_betak, Sigma=Sigma_betak)

    return(betaNew)
}
dphGLM_update_sigma_gaussian_xxr  <- function(yk,Xk, Nk, betak, fix){
    nu = fix$df_sigma
    s2 = fix$s2_sigma
    d = ncol(Xk) 
    if(Nk-d<=0){
        m=0
    }else{
        m=(1 / (Nk - d))
    }
    s2khat =  m*( t(yk - Xk%*%betak) %*% (yk - Xk%*%betak) )
    df = nu + Nk
    scale = min(1,(nu*s2+Nk*s2khat)/df)
    sigmaNew = LaplacesDemon::rinvchisq(1, df=df, scale=scale)
    
    return(sigmaNew)
}

## }}}
## {{{ HMC   : betas (binomial family) }}}

U_xxr      <- function(theta,fix){
    X    = fix$X
    y    = fix$y
    beta = theta
    Sigma_beta     = fix$Sigma_beta
    Sigma_beta.inv = solve(Sigma_beta)
    mu_beta       = matrix(fix$mu_beta, nrow=length(fix$mu_beta))
    d              = nrow(mu_beta) -1       ## -1 to exclude the intercept
    log.p = (-(d+1)/2)*log(2*3.141593) - (1/2)*log(det(Sigma_beta)) - (1/2)*t(beta - mu_beta) %*% Sigma_beta.inv %*% (beta - mu_beta) - sum(y*log(1+exp(- X %*% beta))) - sum((1-y)*log(1+exp( X %*% beta)))
    return( - log.p )
}
Kinectic_xxr <- function(v, theta){
    return( (t(v)%*%G_xxr(theta)%*%v)/2 )
}
grad_U_xxr <- function(theta, fix){
    beta = theta
    Sigma_beta.inv = solve(fix$Sigma_beta)
    mu_beta  = fix$mu_beta
    X        = fix$X
    y        = fix$y
    
    h1 = apply(X, 2, function(Xj) y * Xj * (1/(1+exp(X %*% beta)))) 
    h2 = apply(X, 2, function(Xj) (1-y) * Xj * (1/(1+exp(-X %*% beta)))) 
    if (nrow(X)>1)
    {
        h1 = colSums(h1) 
        h2 = colSums(h2)
    }
    grad.log.p = - t(beta - mu_beta) %*% Sigma_beta.inv + h1 - h2
    grad.log.p = matrix(grad.log.p,  nrow=nrow(theta))
    return ( - grad.log.p )
}
G_xxr      <- function(theta){
    return(diag(length(theta)))
}
q_xxr     <- function(theta, fix){
    return ( matrix(MASS::mvrnorm(1, mu=rep(0, nrow(theta)), Sigma=G_xxr(theta)), nrow=nrow(theta)) )
}
dphGLM_hmc_update_binomial_xxr <- function(yk, Xk, betak, fix, epsilon=0.01, L=40, hmc.iter)
{
    fix$X = Xk
    fix$y = yk
    thetak = matrix(betak, nrow=length(betak))
    for (i in 1:hmc.iter)
    {
        betak = hmc_update_xxr(thetak, epsilon=epsilon, L=L, U_xxr=U_xxr, grad_U_xxr=grad_U_xxr, G_xxr=G_xxr, fix=fix)
        thetak = betak
    }
    return(betak)
}

## }}}
## {{{ Gibbs : pi and Z }}}

dpGLM_update_pi    <- function(Z, K, fix){
    V  = rep(NA,K)
    pi = rep(0,K)
    alpha = fix$alpha
    N = rep(NA, times=K)

    ## computing Nk
    for (k in 1:K){N[k] = sum(Z==k)}
    
    l = 2:K
    V[1]          = rbeta(1, 1 + N[1], alpha + sum(N[l]))
    ## V[1]          = rbeta(1, 1 + N[1], alpha + N[l])
    pi[1]         = V[1]
    
    for (k in 2:(K-1)){
        l = (k+1):K
        V[k]  = rbeta(1, 1 + N[k], alpha + sum(N[l]))  
        ## V[k]  = rbeta(1, 1 + N[k], alpha + N[l])  
        pi[k] = V[k] * prod( 1 - V[ 1:(k-1) ] )
    }
    V[K]  = 1
    pi[K] = V[k] * prod( 1 - V[ 1:(K-1) ] )

    return(pi)
}
dpGLM_update_Z     <- function(y,X, pi, K, theta, family){
    n    = nrow(X)
    Z    = matrix(rep(1,nrow(X)))
    phi  = matrix(NA,ncol=K,nrow=n)

    if (family=='gaussian'){
        for (k in 1:K){
            k.idx  = theta[,'k'] == k
            betak  = theta[k.idx, grepl(colnames(theta), pattern='beta') ]
            sigmak = theta[k.idx, grepl(colnames(theta), pattern='sigma') ]

            shatk = y - X %*% betak
            phi[,k] = pi[k] * dnorm(shatk, mean = 0, sd = sigmak)
        }
    }
    if (family =='binomial'){
        for (k in 1:K){
            k.idx  = theta[,'k'] == k
            betak  = theta[k.idx, grepl(colnames(theta), pattern='beta') ]

            nuk = X %*% betak
            pk  = (y==1)*(1 / (1+exp(-nuk))) + (y==0)*(1 - 1 / (1+exp(-nuk)))

            phi[,k] = pi[k] * pk
            
        }
    }

    phi=t(apply(phi,1, function(row) row/(sum(row))))

    for (i in 1:n){
        Z[i,] = sample(1:K, size=1, prob=phi[i,])
    }

    return(Z)

}

## cluster probabilities
dphGLM_update_countZik <- function(countZik, Z){
    for (i in 1:nrow(Z)){
        countZik[i,Z[i]] = countZik[i,Z[i]] + 1
    }           
    return(countZik)
}
dphGLM_get_pik         <- function(countZik){
    pik           = t( apply(countZik, 1, function(x) x/max(1,sum(x))) )
    colnames(pik) = paste0('p.Z',1:ncol(countZik), sep='')
    return(pik)
}

## }}}

H_xxr <- function(U, K){
    return( U + K )
}
hmc_update_xxr <- function(theta_t, epsilon, L, U_xxr, grad_U_xxr, G_xxr, fix)
{
    v.current  = q_xxr(theta_t, fix)

    v     = v.current
    theta = theta_t

    ## Leapfrog Method together with Modified Euller method
    v = v - (epsilon/2)*grad_U_xxr(theta, fix)
    for (l in 1:(L-1)){
        theta = theta + epsilon * solve(G_xxr(theta)) %*% v
        v     = v     - epsilon * grad_U_xxr(theta, fix)
    }
    theta = theta +   epsilon   * solve(G_xxr(theta)) %*% v
    v     = v     - (epsilon/2) * grad_U_xxr(theta, fix)
    v     = - v

    u     = runif(1,0,1)

    ## print(paste0("U new: ", U(theta, fix)))
    ## print(paste0("k new: ", Kinectic_xxr(v, theta)) )
    ## print(paste0("H new: ", H(U(theta, fix),Kinectic_xxr(v, theta)) ), sep='')
    ## print(paste0("U_t", U(theta_t, fix)), sep='') 
    ## print(paste0("k_t", Kinectic_xxr(v.current, theta_t)), sep='')
    ## print(paste0("H_t", H(U(theta_t, fix), Kinectic_xxr(v.current, theta_t))), sep='')

    alpha = exp( - H_xxr(U_xxr(theta, fix),Kinectic_xxr(v, theta)) + H_xxr(U_xxr(theta_t, fix), Kinectic_xxr(v.current, theta_t)) )
    if (u <= alpha){
        return (theta)
    }else{
        return (theta_t)
    }
}

## MCMC
dpGLM_mcmc_xxr            <- function(y, X, weights, K, fix, family, mcmc, epsilon, leapFrog, n.display, hmc.iter=1)
{
    ## Constants
    ## ---------
    d   = ncol(X) - 1                       # d is the number of covars, we subtract the intercept as X has a column with ones
    n   = nrow(X)                           # sample size
    N   = mcmc$burn.in + mcmc$n.iter          # number of interations in the MCMC algorithm
    ## fix = dphGLM_get_constants_xxr(family=family, d=d) # list with values of the parameters of the priors/hyperpriors

    ## Initialization
    ## --------------
    Z            = matrix(rep(1,nrow(X)))
    theta        = dpGLM_get_inits_xxr(K=K, d=d, family=family, fix=fix)
    countZik     = stats::setNames(data.frame(matrix(0, ncol=K,nrow=n)), nm=paste0('count.z',1:K , sep='')) 
    countZik[,1] = 1
    n.clusters   = 1                                                                                 

    ## MCMC
    ## ----
    samples            = list()
    if (family=='gaussian') {
        samples$samples    = data.table::setnames(data.table::data.table(matrix(ncol=ncol(theta),nrow=0)), c('k',paste0("beta", 1:(d+1), sep=''), 'sigma')) 
    }else{
        samples$samples    = data.table::setnames(data.table::data.table(matrix(ncol=ncol(theta),nrow=0)), c('k',paste0("beta", 1:(d+1), sep='')) ) 
    }
    colnames(theta) = names(samples$samples)
    samples$pik        = NA                                                                  
    samples$n.clusters = NA                                                                  

    ## meta
    ## ----
    if (n.display==0) n.display = mcmc$burn.in + mcmc$n.iter + 2
    Khat      = 1                ## maximum active cluster used in a iteration 
    mh.accept.info <<- list(trial = 0, accepted= 0, average.acceptance=0) 
    
    pb <- txtProgressBar(min = 0, max = mcmc$burn.in+mcmc$n.iter, style = 3, width=70)
    for (iter in 1:(mcmc$burn.in+mcmc$n.iter)){
        setTxtProgressBar(pb, iter)

        ## parameters
        ## ----------
        theta  = dpGLM_update_theta_xxr(y=y, X=X, Z=Z, K=K, theta=theta,  fix, family, epsilon, leapFrog, hmc.iter)
        pi     = dpGLM_update_pi(Z=Z, K=K, fix=fix)
        Z      = dpGLM_update_Z(y=y,X=X, pi=pi, K=K, theta=theta, family)
        ## update countZik (number of times) i was classified in cluster k
        ## ---------------
        countZik = dphGLM_update_countZik(countZik, Z)
        ## count the number of unique clusters in the iteration
        ## ----------------------------------------------------
        n.clusters = c(n.clusters, length(unique(Z)))

        ## saving samples
        ## --------------
        if(iter > mcmc$burn.in){
            thetaNew           = matrix(theta[unique(Z),], ncol=ncol(theta), nrow=length(unique(Z)))
            colnames(thetaNew) = colnames(samples$samples)
            samples$samples = rbind(samples$samples,  thetaNew)
        }


        Khat = max(Khat,length(table(Z)))
        mh.accept.info$average.acceptance <<- (mh.accept.info$average.acceptance + mh.accept.info$accepted/mh.accept.info$trial) /2
        ## Debug/Monitoring message --------------------------
        if (iter %% n.display == 0 & iter > n.display){
            msg <- paste0('\n\n-----------------------------------------------------', '\n'); cat(msg)
            msg <- paste0('MCMC in progress .... \n'); cat(msg)
            cat('\n')
            msg <- paste0('Family of the link function of the GLMM components: ', family, '\n'); cat(msg)
            msg <- paste0('Burn-in: ', mcmc$burn.in,'\n'); cat(msg)
            msg <- paste0('Number of MCMC samples per chain: ', mcmc$n.iter,'\n'); cat(msg)
            cat('\n')
            msg <- paste0('Iteration ', iter, '\n'); cat(msg)
            cat('\n')
            msg <- paste0('Acceptance rate for beta : ', round(mh.accept.info$accepted/mh.accept.info$trial,4), '\n'); cat(msg)
            msg <- paste0('Average acceptance rate  : ', round(mh.accept.info$average.acceptance,4) , '\n'); cat(msg)
            cat('\n')
            msg <- paste0('Maximum Number of cluster allowed (K): ', K,'\n', sep='');cat(msg)
            msg <- paste0('Maximum Number of cluster activated  : ', max(Khat), '\n',sep='');cat(msg)
            msg <- paste0('Current number of active clusters    : ', length(table(Z)), '\n',sep='');cat(msg)
            cat('\n')
            msg <- paste0('Distribution of total number of clusters:\n', sep='');cat(msg)
            print(formatC(table(n.clusters, dnn='')/sum(table(n.clusters)), digits=4, format='f'))
            msg <- paste0("Percentage of data classified in each clusters k at current iteraction ", sep=''); cat(msg)
            print(round(100*table(Z, dnn='')/sum(table(Z)),1))
        }
        ## ---------------------------------------------------
    }
    samples$pik        = dphGLM_get_pik(countZik)
    samples$n.clusters = n.clusters

    return(samples)
}

## {{{ doc }}}

#' Hierarchical Dirichlet Process GLM
#'
#' The function estimates a semi-parametric mixture of Generalized
#' Linear Models. It uses a (hierarchical) Dependent Dirichlet Process
#' Prior for the mixture probabilities. 
#'
#' @param formula1 a single symbolic description of the linear model of the
#'                 mixture GLM components to be fitted. The sintax is the same
#'                 as used in the \code{\link{lm}} function.
#' @param formula2 eihter NULL (default) or a single symbolic description of the
#'                 linear model of the hierarchical component of the model.
#'                 It specifies how the average parameter of the base measure
#'                 of the Dirichlet Process Prior varies linearly as a function
#'                 of group level covariates. If \code{NULL}, it will use
#'                 a single base measure to the DPP mixture model.
#' @param data a data.frame with all the variables specified in \code{formula1}
#'             and \code{formula2}
#' @param weights numeric vector with the same size as the number of rows of the data. It must contains the weights of the observations in the data set. NOTE: FEATURE NOT IMPLEMENTED YET
#' @param mcmc a list containing elements named \code{burn.in} (required, an
#'             integer greater or equal to 0 indicating the number iterations used in the
#'             burn-in period of the MCMC) and \code{n.iter} (required, an integer greater or
#'             equal to 1 indicating the number of iterations to record after the burn-in
#'             period for the MCMC).
#' @param K an integer indicating the maximum number of clusters to truncate the
#'          Dirichlet Process Prior in order to use the blocked Gibbs sampler.
#' @param fix either NULL or a list with the constants of the model. If not NULL,
#'            if must contain a vector named \code{mu_beta}, whose size must be
#'            equal to the number of covariates specified in \code{formula1}
#'            plus one for the constant term; \code{Sigma_beta}, which must be a squared
#'            matrix, and each dimension must equal to the size of the vector \code{mu_beta}; 
#'            and \code{alpha}, which must a single number. If @param family is 'gaussian',
#'            then it must also contains \code{s2_sigma} and \code{df_sigma}, both
#'            single numbers. If NULL, the defaults are \code{mu_beta=0},
#'            \code{Sigma_beta=diag(10)}, \code{alpha=1}, \code{df_sigma=10},
#'            \code{d2_sigma=10} (all with the dimension automatically set to the
#'            correct values).
#' @param family a character with either 'gaussian', 'binomial', or 'multinomial'.
#'               It indicates the family of the GLM components of the mixture model.
#' @param epsilon numeric, used when \code{family='binomial'} or \code{family='multinomial'}.
#'                It is used in the Stormer-Verlet Integrator (a.k.a leapfrog integrator)
#'                to solve the Hamiltonian Monte Carlo in the estimation of the model.
#'                Default is 0.01.
#' @param leapFrog an integer, used when \code{family='binomial'} or \code{family='multinomial'}.
#'                 It indicates the number of steps taken at each iteration Hamiltonian
#'                 Monte Carlo for the Stormer-Verlet Integrator. Default is 40.
#' @param n.display an integer indicating the iteration to display information
#'                  about the estimation process. If zero, it does not display any information.
#'                  Note: displaying informaiton at every iteration (n.display=1) may increase
#'                  the time to estimate the model slightly. 
#' @param hmc_iter an integer, used when \code{family='binomial'} or \code{family='multinomial'}.
#'                 It indicates the number of HMC interation for each Gibbs iteration.
#'                 Default is 1.
#' @return The function returns a list with elements \code{samples}, \code{pik}, \code{max_active},
#'         \code{n.iter}, \code{burn.in}, and \code{time.elapsed}. The \code{samples} element
#'         contains a MCMC object (from \pkg{coda} package) with the samples from the posterior
#'         distribution. The \code{pik} is a \code{n x K} matrix with the estimated
#'         probabilities that the observation $i$ belongs to the cluster $k$
#' 
#'
#' @details
#' The estimation is conducted using Blocked Gibbs Sampler if the output
#' variable is gaussian distributed. It uses Metropolis-Hastings inside Gibbs if
#' the output variable is binomial or multinomial distributed.
#' This is specified using the parameter \code{family}. See
#'
#' Ishwaran, H., & James, L. F., Gibbs sampling methods for stick-breaking priors,
#' Journal of the American Statistical Association, 96(453), 161–173 (2001). 
#'
#' Neal, R. M., Markov chain sampling methods for dirichlet process mixture models,
#' Journal of computational and graphical statistics, 9(2), 249–265 (2000).
#'
#' @examples
#' data    = hdpGLM_simulateData(n=5000,nCov=4, K=5, family='gaussian')
#' mcmc    = list(burn.in = 0,  n.iter = 2000)
#' samples = hdpGLM(y~., data=data$data, mcmc=mcmc, family='gaussian', n.display=30, K=100)
#'
#' summary(samples, nk=6)
#' 
#' plot(samples, K=5, true_parameters=data$parameters, plot.hist=F, title='Posterior Distribution', prop.time.active=.95)
#' 
#' plot(samples, K=5, true_parameters=data$parameters, plot.hist=F, title='Posterior Distribution', separate=T,  plot.hist.k=F)
#'  
#' \dontrun{
#' }
#' @export

## }}}
hdpGLM_xxr <- function(formula1, formula2=NULL, data, mcmc, K=50, fix=NULL, family='gaussian', epsilon=0.01, leapFrog=40, n.display=1000, hmc_iter=1)
{
    if(! family %in% c('gaussian', 'binomial', 'multinomial'))
        stop(paste0('Error: Parameter -family- must be a string with one of the following options : \"gaussian\", \"binomial\", or \"multinomial\"'))

    ## ## construct the regression matrices (data.frames) based on the formula provided
    ## ## -----------------------------------------------------------------------------
    func.call <- match.call(expand.dots = FALSE)
    mat     = .getRegMatrix(func.call, data, weights, formula_number=1)
    y       = mat$y
    X       = mat$X
    weights = ifelse(!is.null(mat$w), mat$w, rep(1,nrow(X)))
    ## Hierarchical covars
    Xj = unlist( ifelse(is.null(formula2), list(NULL), list ( getRegMatrix(func.call, data, weights, formula_number=2)$X) ) ) # list and unlist only b/c ifelse() do not allow to return NULL

    ## get constants
    ## -------------
    d   = ncol(X)  - 1                        # d is the number of covars, we subtract the intercept as X has a column with ones
    dj  = unlist( ifelse(is.null(Xj), list(NULL), list(ncol(Xj) - 1)) ) # list and unlist only b/c ifelse() do not allow to return NULL
    if (is.null(fix))
    {
        fix = dphGLM_get_constants_xxr(family=family, d=d, dj=dj) # list with values of the parameters of the priors/hyperpriors
    }else
    {
        dphGLM_check_constants_xxr(family=family, d=d, dj=dj, fix=fix)
    }
    

    ## get the samples from posterior
    ## ------------------------------
    T.mcmc  = Sys.time()
    if (is.null(Xj))
    {
        samples        =  dpGLM_mcmc_xxr( y, X,    weights, K, fix,  family, mcmc, epsilon, leapFrog, n.display, hmc_iter) #
    }else
    {
        samples        =  hdpGLM_mcmc_xxr(y, X, Xj, weights, K, fix, family, mcmc, epsilon, leapFrog, n.display, hmc_iter)
    }
    T.mcmc  = Sys.time() - T.mcmc


    ## including colnames and class for the output
    ## -------------------------------------------
    if (is.null(Xj)){
        if(family=='gaussian'){
            colnames(samples$samples)  <- c('k', paste0('beta',1:(d+1),  sep=''),'sigma')
        }else{
            colnames(samples$samples)  <- c('k', paste0('beta',1:(d+1),  sep=''))
        }
    }else{
        if(family=='gaussian'){
            colnames(samples$samples)  <- c('k', paste0('beta',1:(d+1),  sep=''), 'sigma', paste0('tau',1:(dj+1),  sep=''))
        }else{
            colnames(samples$samples)  <- c('k', paste0('beta',1:(d+1),  sep=''), paste0('tau',1:(dj+1),  sep=''))
        }
    }

    samples$samples                   = coda::as.mcmc(samples$samples)
    attr(samples$samples, 'mcpar')[2] = mcmc$n.iter - mcmc$burn.in
    class(samples)                    = 'dpGLM'

    samples$time_elapsed = T.mcmc
    return(samples)
}
