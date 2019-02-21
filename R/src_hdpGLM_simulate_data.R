
## {{{ docs }}}

#' Simulated a Data Set from hdpGLM model
#'
#' This function generates a data set from the Hierarchical
#' Dirichlet Process of Generalized Linear Model (hdpGLM)
#'
#'
#' @param n integer, the sample size of the data. If there are multiple contexts, each context will have n cases
#' @param K integer, the number of clusters. If there are multiple contexts, K is the average number of clusters across contexts, and each context gets a number of clusters sampled from a Poisson distribution, except if \code{same.K} is \code{TRUE}.
#' @param nCov integer, the number of covariates of the
#'             GLM components
#' @param nCovj an integer indicating the number of
#'              covariates determining the average parameter of the
#'              base measure of the Dirichlet process prior
#' @param J an integer representing the number of contexts 
#' @param parameters either NULL or a list with the parameters
#'                   to generate the model. If not NULL, it must
#'                   contain a sublist name beta, a vector named tau,
#'                   and a vector named pi. The sublist beta must be
#'                   a list of vectors, each one with size nCov+1
#'                   to be the coefficients of the GLM mixtures
#'                   components that will generate the data.
#'                   For the vector tau, if nCovj=0 (single-context case)
#'                   then it must be a 1x1 matrix containing 1.
#'                   If ncovj>0, it must be a (nCov+1)x(nCovj+1) matrix.
#'                   The vector pi must add up to 1 and have length K.
#' @param pi either NULL or a vector with length K that add up to 1.
#'           If not NULL, it determines the mixture probabilities
#' @param same.K boolean, used when data is sampled from more than one context. If \code{TRUE} all contexts get the same number of clusters. If \code{FALSE}, each context gets a number of clusters sampled from a Poisson distribution with expectation equals to \code{K} (current not implemented)
#' @param context.effect either \code{NULL} or a two dimensional integer vector. If it is \code{NULL}, all the coefficients (\code{beta}) of the individual level covariates are functions of context-level features (\code{tau}). If it is not \code{NULL}, the first component of the vector indicates the index of the lower level covariate (\code{X}) whose linear effect \code{beta} depends on context (\code{tau}) (0 is the intercept). The second component indicates the index context-level covariate (\code{W}) whose linear coefficient (\code{tau}) is non-zero.
#' @param same.clusters.across.contexts boolean, if \code{TRUE} all the contexts will have the same number of clusters AND each cluster will have the same coefficient \code{beta}.
#' @param context.dependent.cluster integer, indicates which cluster will be context-dependent. If \code{zero}, all clusters will be context-dependent 
#' @param seed a seed for \code{\link{set.seed}}
#' @inheritParams hdpGLM 
#'
#' @return The function returns a list with the data set, the parameters,
#'         the vector Z indicating the cluster of each observation, and
#'         the number of clusters.
#'
#' @examples
#' data = hdpGLM_simulateData(n=2000, K=2, nCov=0, family='gaussian')
#' @export
## }}}
hdpGLM_simulateData <- function(n, K, nCov, nCovj=0, J=1, parameters=NULL, pi=NULL, family, same.K=FALSE, seed=sample(1:777,1), 
                                context.effect=NULL, same.clusters.across.contexts=NULL, context.dependent.cluster=0)
{

    if(nCovj==0 | J < 2) 
        dat = dpGLM_simulateData_main(n, K, nCov, nCovj=NULL, parameters=parameters, pi=pi, family=family, seed=seed) 
    if(nCovj> 0){ 
        dat = hdpGLM_simulateData_main(n, K, nCov, nCovj=nCovj, J=J, parameters=parameters, pi=pi, family=family, same.K, seed=seed, 
                                       context.effect, same.clusters.across.contexts=same.clusters.across.contexts, context.dependent.cluster) 
    }
    return(dat)
}

## =====================================================
## dpGLM
## =====================================================
dpGLM_simulateData_main <- function(n, K, nCov, nCovj=NULL, parameters=NULL, pi=NULL, family, seed=sample(1:777,1))
{
    ## error handling
    if (n %% 2) stop("Sample size n must be an even number")
    if(! family %in% c('gaussian', 'binomial', 'multivariate')) stop(paste0('Error: family must be one element of the set : {', paste0(c('gaussian', 'binomial', 'multinomial'), collapse=','),'}'))
    if (!is.null(parameters)){
        ## check the format of parameters
        if(! all(names(parameters) %in% c('pi', 'beta') | names(parameters) %in% c('pi', 'beta', 'K')) | is.null(names(parameters)) ) stop('\n\n\"parameters\" must be a named list. The names must be \'beta\' and \'pi\'.\n ')
        if(! length(parameters$beta) == K) stop(paste0('\n\nThe element \'beta\' of the list \'parameters\' must be a list itself with size ', K, '. Each element of that list must be a vector of ', nCov+1, ' linear parameter(s).' ,  sep=''))
        if(! length(parameters$pi) == K) stop(paste0('\n\nThe element \'pi\' of the list \'parameters\' must be a vector with size ', K,'.' ,  sep=''))
        if(! sum(parameters$pi) == 1) stop ('\n\nThe element \'pi\' of the list \'parameters\' must add up to 1.')
        if(! all(parameters$beta %>% purrr::map(., ~length(.)) %>% unlist == nCov +1)) stop(paste0('\n\nAll elements of the list \'beta\' of the list \'parameters\' must have length ', nCov+1, sep=''))
    }

    if (is.null(parameters)) {
        parameters <- dpGLM_simulateParameters(K=K, nCov=nCov, pi, seed=seed) 
    }
    if(family=='gaussian'){sim_data = .dpGLM_simulateData_gaussian(n=n, K=K, nCov=nCov, parameters, seed=seed)}
    if(family=='binomial'){sim_data = .dpGLM_simulateData_binomial(n=n, K=K, nCov=nCov, parameters, seed=seed)}

    class(sim_data) = "dpGLM_data"
    return(sim_data)
}
## {{{ docs }}}
#' Simulated a Parameters used to create Data Sets from dpGLM model
#'
#' This function generates parameters that can be used to
#' simulated data sets from the Hierarchical
#' Dirichlet Process of Generalized Linear Model (dpGLM)
#'
#' @inheritParams hdpGLM_simulateData
#'
#' @return The function returns a list with the parameters
#'         used to generate data sets from the dpGLM model.
#'         This list can be used in the function \code{\link{hdpGLM_simulateData}}
#' 
#' @examples
#' parameters = dpGLM_simulateParameters(nCov=2, K=2) 
#'
#' @export
## }}}
dpGLM_simulateParameters <- function(nCov, nCovj=NULL, K=NULL, pi=NULL, seed=NULL)
{
    if(is.null(seed)) seed <- sample(1:777,1)
    parameters <- list(pi=NA ,beta = list(), K=K)

    ## pi
    ## --
    if(is.null(pi)){
        parameters$pi <- rep(1/K,K)
    }else{
        if(length(pi) != K) stop('\n\n Size of pi must be equal to K\n')
        parameters$pi <- pi
    }

    ## beta
    ## ----
    parameters$beta <- lapply(rep(nCov, K), function(nCov) stats::runif(n=nCov+1, -10,10))
    parameters$beta <- parameters$beta %>% purrr::map(., ~stats::setNames(., paste0('beta', 1:(nCov+1), sep='')))

    return(parameters)
}
.dpGLM_simulateData_gaussian <- function(n, K, nCov, parameters, seed)
{
    set.seed(seed)
    n = floor(n*parameters$pi)

    beta  = parameters$beta
    X     = list()
    y     = list()
    for (k in 1:K){
        epsilon   <- stats::rnorm(n=n[k], 0, 1)

        if(nCov!=0){
            X[[k]] <- MASS::mvrnorm(n=n[k], mu = rep(0,nCov), Sigma = diag(nCov))
            y[[k]] <- base::cbind(1,X[[k]]) %*% beta[[k]] + epsilon
        }else{
            X[[k]] <- NULL
            y[[k]] <- beta[[k]] + epsilon
        }
    }
    if(nCov!=0){
        data   = stats::setNames(data.frame(y = unlist(y), X = data.frame(do.call(base::rbind, X))), nm=c('y', paste0('X',1:nCov)))
    }else{
        data   = data.frame(y = unlist(y))
    }
    return ( list(data=data, Z = rep(1:K,n), parameters=parameters) )
}
.dpGLM_simulateData_binomial <- function(n, K, nCov, parameters,  seed)
{
    set.seed(seed)
    n = floor(n*parameters$pi)

    beta  = parameters$beta
    X     = list()
    y     = list()
    p     = list()
    nu    = list()
    for (k in 1:K){
        epsilon  <- stats::rnorm(n=n[k], 0, 1)

        if(nCov!=0){
            X[[k]]   <- MASS::mvrnorm(n=n[k], mu = rep(0,nCov), Sigma = diag(nCov))
            nu[[k]]  <- base::cbind(1,X[[k]]) %*% beta[[k]]
            p[[k]]   <- 1/(1+exp(-nu[[k]]))
        }else{
            X[[k]]   <- NULL
            nu[[k]]  <- rep(beta[[k]], n[k]) 
            p[[k]]   <- 1/(1+exp(-nu[[k]]))
        }
        y[[k]] = rep(NA, n[k])
        for (i in 1:n[k]) y[[k]][i] = stats::rbinom(n=1, size=1, prob=p[[k]][i])
    }
    if(nCov!=0){
        data   = stats::setNames(data.frame(y = unlist(y), X = data.frame(do.call(base::rbind, X))), nm=c('y', paste0('X',1:nCov)))
    }else{
        data   = data.frame(y = unlist(y))
    }
    return ( list(data=data, Z = rep(1:K,n), parameters=parameters) )
}
## to be completed
.dpGLM_simulateData_multinomial <- function(n, K=2, nCov=0, nCat, pi=pi,  family, seed){
    
    set.seed(seed)
    ## rand('state', state);
    ## randn('state', state);

    # n is the total sample size (must be an even number)
    ## J = J = nClass = 4; # dep var  : number of categories of y: y in {1,2,3,4}
    ## D = nCov = l = nVar = 5;   # covars   : dimension of X (does not include the intercept !)
    ## K = 2               # clusters : pre fixed number of components (clusters)
    ## pi = rep(1/K,K)
    ## n = 2*5000
    n = n*pi
    J = nCat

    alpha = list()
    beta  = list()
    tau   = list() 
    nu    = list() 
    mu    = list()
    sd    = list()
    X     = list()
    y     = list()
    eta   = list()
    p     = list()
    for (k in 1:K){
        ## hyperpriors on G_o (for tau and nu)
        tau[[k]]  <-  sqrt( exp( stats::rnorm(n=1, mean=0, sd=.1)) );    ## ln(tau^2) ~ N(0,.1^2)
        nu[[k]]   <-  sqrt( exp( stats::rnorm(n=1, mean=0, sd= 2)) );     ## ln(nu^2)  ~ N(0,2^2)

        ## Priors for theta defined by G_o
        mu[[k]]     <- stats::rnorm(n=nCov, mean=0, sd=1);                    ## mu          ~ N_d(0,1)   which is expect of X ~ N_d(mu, sigma^2 * I)
        sd[[k]]     <- sqrt( exp( stats::rnorm(n=nCov, mean=0, sd=2) ) );  ## ln(sigma^2) ~ N_d(0,2^2) which is var    or X ~ N_d(mu, sigma^2 * I)
        alpha[[k]]  <- stats::rnorm(n=J,mean=0, sd=tau[[k]] )
        beta[[k]]   <- MASS::mvrnorm(n=nCov,mu=rep(0,J), Sigma=nu[[k]]^2  * diag(J))

        ## generating X's
        ## --------------
        X[[k]]  <-  MASS::mvrnorm(n=n[k], mu = mu[[k]], Sigma = sd[[k]] * diag(nCov))

        ## computing y, eta = alpha + X^T beta, and p = exp(eta)/sum(exp(eta))
        ## -------------------------------------------------------------------
        eta[[k]]   <-  base::cbind(1, X[[k]]) %*% base::rbind(alpha[[k]], beta[[k]])
        eta[[k]]   <-  eta[[k]] + MASS::mvrnorm(n[k], mu=rep(0,J), Sigma = .05 * diag(J))
        p[[k]]     <-  t(apply(eta[[k]], 1, function(eta.kj) exp(eta.kj) / sum(exp(eta.kj)) ))
        y[[k]] <- rep(NA, n[k])
        for (i in 1:n[k]) y[[k]][i]  = which(stats::rmultinom(n=1, size=1, prob=p[[k]][i,]) == 1)

    }
    return ( list(y = unlist(y),
                  X = stats::setNames(data.frame(do.call(base::rbind, X)), nm=paste0('X',1:nCov) ),
                  Z = rep(1:K,n)) )
}


## =====================================================
## hdpGLM
## =====================================================
hdpGLM_simulateData_main <- function(n, K, nCov, nCovj=NULL, J, parameters=NULL, pi=NULL, family, same.K, seed=sample(1:777,1), context.effect, same.clusters.across.contexts, context.dependent.cluster)
{
    ## error handling ----------------------------------------
    if (n %% 2) stop("Sample size n must be an even number")
    if(! family %in% c('gaussian', 'binomial', 'multivariate')) stop(paste0('Error: family must be one element of the set : {', paste0(c('gaussian', 'binomial', 'multinomial'), collapse=','),'}'))
    if (!is.null(parameters)){
        ## check the format of parameters
        if(! all(names(parameters) %in% c('pi', 'beta') | names(parameters) %in% c('pi', 'beta', 'K')) | is.null(names(parameters)) ) stop('\n\n\"parameters\" must be a named list. The names must be \'beta\' and \'pi\'.\n ')
        if(! length(parameters$beta) == K) stop(paste0('\n\nThe element \'beta\' of the list \'parameters\' must be a list itself with size ', K, '. Each element of that list must be a vector of ', nCov+1, ' linear parameter(s).' ,  sep=''))
        if(! length(parameters$pi) == K) stop(paste0('\n\nThe element \'pi\' of the list \'parameters\' must be a vector with size ', K,'.' ,  sep=''))
        if(! sum(parameters$pi) == 1) stop ('\n\nThe element \'pi\' of the list \'parameters\' must add up to 1.')
        if(! all(parameters$beta %>% purrr::map(., ~length(.)) %>% unlist == nCov +1)) stop(paste0('\n\nAll elements of the list \'beta\' of the list \'parameters\' must have length ', nCov+1, sep=''))
    }
    ## -----------------------------------------------------


    if (is.null(parameters)) {parameters <- hdpGLM_simulateParameters(K=K, nCov=nCov, nCovj=nCovj, J=J, pi, same.K, seed=seed, 
                                                                      context.effect=context.effect,same.clusters.across.contexts=same.clusters.across.contexts, context.dependent.cluster=context.dependent.cluster)}
    if(family=='gaussian')   {sim_data = .hdpGLM_simulateData_gaussian(n=n, K=K, nCov=nCov, nCovj=nCovj, parameters, seed=seed)}
    if(family=='binomial')   {sim_data = .hdpGLM_simulateData_binomial(n=n, K=K, nCov=nCov, nCovj=nCovj, parameters, seed=seed)}

    class(sim_data) = "hdpGLM_data"
    return(sim_data)
}

## {{{ docs }}}

#' Simulated a Parameters used to create Data Sets from hdpGLM model
#'
#' This function generates parameters that can be used to
#' simulated data sets from the Hierarchical
#' Dirichlet Process of Generalized Linear Model (hdpGLM)
#'
#' @inheritParams hdpGLM_simulateData
#'
#' @return The function returns a list with the parameters
#'         used to generate data sets from the hdpGLM model.
#'         This list can be used in the function \code{hdpGLM_simulateData}
#' 
#' @examples
#' parameters = hdpGLM_simulateParameters(nCov=2, K=2, nCovj=3, J=20,
#'                 same.clusters.across.contexts=FALSE, context.dependent.cluster=0) 
#'
#' @export
## }}}
hdpGLM_simulateParameters <- function(nCov, K=NULL, nCovj=NULL, J=NULL, pi=NULL, same.K, seed=NULL, context.effect=NULL, same.clusters.across.contexts, context.dependent.cluster)
{
    Dw = nCovj
    Dx = nCov
    if(is.null(seed)) seed <- base::sample(1:777,1)

    ## K
    ## -
    ## if(!same.K){
    ##     Ks = stats::rpois(J, lambda=K)
    ##     Ks[Ks==0] = 1
    ## }else{
    ##     Ks = rep(K, J)
    ## }

    parameters <- list(pi=NA ,beta = list(), K=K)

    ## pi
    ## --
    if(is.null(pi)){
        parameters$pi <- rep(1/K,K)
    }else{
        if(length(pi) != K) stop('\n\n Size of pi must be equal to K\n')
        parameters$pi <- pi
    }

    ## tau
    ## ---
    if(Dw ==0) {
        parameters$tau = as.matrix(0)  
        W              = as.matrix(1)  
    }
    if(Dw  >0) {
        ## user can specify specific betas to be affected by some specific context-level covars only, as well as which context level covar will affect that beta (instead of all context-level covars)
        if(!is.null(context.effect)){
            n.of.non.zero.taus = length(context.effect$on.betas) + length(context.effect$by.Ws)
            ## all taus are zero...
            parameters$tau = matrix(0, nrow=Dw+1, ncol=Dx+1)
            ## ... but the intercept and the coefficients tau of those W's that the user wants to affect beta
            parameters$tau[c(1, context.effect$by.Ws+1), context.effect$on.betas+1] = as.matrix(MASS::mvrnorm(n=1 , mu=rep(0,n.of.non.zero.taus), Sigma=5*diag(n.of.non.zero.taus)))  
        }else{
            ## if they don't all taus are different from zero and all W's affect all betas
            parameters$tau = as.matrix(MASS::mvrnorm(n=Dw+1, mu=rep(0, Dx+1), Sigma=5*diag(Dx+1)))  
        }
        Wprime         = matrix(MASS::mvrnorm(n=J, mu=rep(0, Dw) , Sigma=5*diag(Dw)), nrow=J, ncol=Dw)
        W              = cbind(1,Wprime)
        parameters$W   = W
    }
    colnames(parameters$W)   = paste0("W", 0:Dw, sep='')
    colnames(parameters$tau) = paste0("beta", 0:Dx, sep='')
    rownames(parameters$tau) = paste0("w", 0:Dw, sep='')
    
    ## beta
    ## ----
    parameters$beta <- lapply(1:J, function(J) rep(list(rep(NA, Dx+1)),K))
    tau = parameters$tau
    for (j in 1:J) {
        for (k in 1:K) {
            parameters$beta[[j]][[k]] = MASS::mvrnorm(n=1, mu=W[j, ] %*% tau, diag(nCov+1))
            names(parameters$beta[[j]][[k]] ) = paste0('beta', 1:(nCov+1))
        }
        names(parameters$beta[[j]]) = paste0('k=', 1:K, sep='')
        
    }
    if(!is.null(same.clusters.across.contexts)){
        beta = parameters$beta[[1]]
        for (j in 1:J) {
            for (k in 1:K) {
                if(!is.null(context.effect)){
                    if (context.dependent.cluster == 0 | k != context.dependent.cluster) {
                        parameters$beta[[j]][[k]][-(context.effect[1]+1)] = beta[[k]][-(context.effect[1]+1)]
                    }else{
                        parameters$beta[[j]][[k]] = beta[[k]]
                    }
                }else{
                    parameters$beta[[j]][[k]] = beta[[k]]
                }
            }
        }
    }
    names(parameters$beta) = paste0('j=', 1:J, sep='')

    return(parameters)
}

.hdpGLM_simulateData_gaussian <- function(n, K, nCov, nCovj=NULL, parameters, seed)
{
    set.seed(seed)
    n = floor(n*parameters$pi)

    W     = parameters$W
    J     = nrow(W)
    betas = parameters$beta
    X     = matrix(NA, ncol=nCov, nrow=0) # +1 to C_i
    Z     = c()
    C     = c()
    y     = c()
    for (j in 1:J)
    {
        for (k in 1:K){
            epsilon   <- stats::rnorm(n=n[k], 0, 1)

            if(nCov!=0){
                X_tmp <- MASS::mvrnorm(n=n[k], mu = rep(0,nCov), Sigma = diag(nCov))
                y_tmp <- base::cbind(1,X_tmp) %*% betas[[j]][[k]] + epsilon
            }else{
                X_tmp <- NULL
                y_tmp <- rep(1, n[k]) * betas[[j]][[k]] + epsilon
            }
            if(!is.null(X_tmp)) X = base::rbind(X, X_tmp)
            Z = c(Z, rep(k, n[k]))
            C = c(C, rep(j, n[k]))
            y = c(y, y_tmp)
        }
    }
    if(nCov!=0){
        data = data.frame(y=y, X=X, C=C) %>%
            stats::setNames(., nm = c('y', paste0('X', 1:nCov, sep=''), 'C')) %>%
            dplyr::full_join(., parameters$W %>% tibble::as_data_frame(.) %>% dplyr::mutate(C=1:J), by='C' )
    }else{
        data  = data.frame(y = y, C=C)
    }
    return ( list(data=data, Z = Z, parameters=parameters) )
}
.hdpGLM_simulateData_binomial <- function(n, K, nCov, nCovj=NULL, parameters,  seed)
{
    set.seed(seed)
    n = floor(n*parameters$pi)
    
    W     = parameters$W
    J     = nrow(W)
    betas = parameters$beta
    X     = matrix(NA, ncol=nCov, nrow=0) # +1 to C_i
    Z     = c()
    C     = c()
    y     = c()
    for (j in 1:J)
    {
        for (k in 1:K){
            epsilon  <- stats::rnorm(n=n[k], 0, 1)
            
            if(nCov!=0){
                X_tmp  <- MASS::mvrnorm(n=n[k], mu = rep(0,nCov), Sigma = diag(nCov))
                nu     <- base::cbind(1,X_tmp) %*% betas[[j]][[k]]
                p_tmp  <- 1/(1+exp(-nu))
            }else{
                X_tmp  <- NULL
                nu     <- rep(betas[[j]][[k]], n[k]) 
                p_tmp  <- 1/(1+exp(-nu))
            }
            y_tmp = rep(NA, n[k])
            for (i in 1:n[k]) y_tmp[i] = stats::rbinom(n=1, size=1, prob=p_tmp[i])

            if(!is.null(X_tmp)) X = base::rbind(X, X_tmp)
            Z = c(Z, rep(k, n[k]))
            C = c(C, rep(j, n[k]))
            y = c(y, y_tmp)
        }
    }
    if(nCov!=0){
        data = data.frame(y=y, X=X, C=C) %>%
            stats::setNames(., nm = c('y', paste0('X', 1:nCov, sep=''), 'C')) %>%
            dplyr::full_join(., parameters$W %>% tibble::as_data_frame(.) %>% dplyr::mutate(C=1:J), by='C' )
    }else{
        data  = data.frame(y = y, C=C)
    }
    return ( list(data=data, Z = Z, parameters=parameters) )
}
## to be completed
.hdpGLM_simulateData_multinomial <- function(n, K=2, nCov=0, nCat, pi=pi,  family, seed)
{
    
    set.seed(seed)
    ## rand('state', state);
    ## randn('state', state);

    # n is the total sample size (must be an even number)
    ## J = J = nClass = 4; # dep var  : number of categories of y: y in {1,2,3,4}
    ## D = nCov = l = nVar = 5;   # covars   : dimension of X (does not include the intercept !)
    ## K = 2               # clusters : pre fixed number of components (clusters)
    ## pi = rep(1/K,K)
    ## n = 2*5000
    n = n*pi
    J = nCat

    alpha = list()
    beta  = list()
    tau   = list() 
    nu    = list() 
    mu    = list()
    sd    = list()
    X     = list()
    y     = list()
    eta   = list()
    p     = list()
    for (k in 1:K){
        ## hyperpriors on G_o (for tau and nu)
        tau[[k]]  <-  sqrt( exp( stats::rnorm(n=1, mean=0, sd=.1)) );    ## ln(tau^2) ~ N(0,.1^2)
        nu[[k]]   <-  sqrt( exp( stats::rnorm(n=1, mean=0, sd= 2)) );     ## ln(nu^2)  ~ N(0,2^2)

        ## Priors for theta defined by G_o
        mu[[k]]     <- stats::rnorm(n=nCov, mean=0, sd=1);                    ## mu          ~ N_d(0,1)   which is expect of X ~ N_d(mu, sigma^2 * I)
        sd[[k]]     <- sqrt( exp( stats::rnorm(n=nCov, mean=0, sd=2) ) );  ## ln(sigma^2) ~ N_d(0,2^2) which is var    or X ~ N_d(mu, sigma^2 * I)
        alpha[[k]]  <- stats::rnorm(n=J,mean=0, sd=tau[[k]] )
        beta[[k]]   <- MASS::mvrnorm(n=nCov,mu=rep(0,J), Sigma=nu[[k]]^2  * diag(J))

        ## generating X's
        ## --------------
        X[[k]]  <-  MASS::mvrnorm(n=n[k], mu = mu[[k]], Sigma = sd[[k]] * diag(nCov))

        ## computing y, eta = alpha + X^T beta, and p = exp(eta)/sum(exp(eta))
        ## -------------------------------------------------------------------
        eta[[k]]   <-  base::cbind(1, X[[k]]) %*% base::rbind(alpha[[k]], beta[[k]])
        eta[[k]]   <-  eta[[k]] + MASS::mvrnorm(n[k], mu=rep(0,J), Sigma = .05 * diag(J))
        p[[k]]     <-  t(apply(eta[[k]], 1, function(eta.kj) exp(eta.kj) / sum(exp(eta.kj)) ))
        y[[k]] <- rep(NA, n[k])
        for (i in 1:n[k]) y[[k]][i]  = which(stats::rmultinom(n=1, size=1, prob=p[[k]][i,]) == 1)

    }
    return ( list(y = unlist(y),
                  X = stats::setNames(data.frame(do.call(base::rbind, X)), nm=paste0('X',1:nCov) ),
                  Z = rep(1:K,n)) )
}


