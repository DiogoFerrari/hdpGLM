

.hdpGLM_get_constants      <- function(family, d, Dw)
{
    Dw = ifelse(is.null(Dw), 0, Dw)
    if(family=='gaussian') {        
        sigma_beta = 10
        sigma_tau  = 10
        fix = list(alpha      = 1,
                   mu_beta    = rep(0, d+1),
                   Sigma_beta = sigma_beta * diag(d+1),
                   s2_sigma   = 10,
                   df_sigma   = 10,
                   ## tau
                   mu_tau     = rep(0, Dw+1),
                   Sigma_tau  = sigma_tau * diag(Dw+1)
                   )
    }
    if(family=='binomial') {        
        sigma_beta = 10
        sigma_tau  = 10
        fix = list(alpha      = 1,
                   mu_beta    = rep(0, d+1),
                   Sigma_beta = sigma_beta * diag(d+1),
                   ## tau
                   mu_tau     = rep(0, Dw+1),
                   Sigma_tau  = sigma_tau * diag(Dw+1)
                   )
    }
    return(fix)
}

.dphGLM_check_constants <- function(family=family, d=d, Dw=NULL, fix)
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
        if (any(! names(fix) %in% c('mu_beta', 'Sigma_beta', 'alpha', 's2_sigma', 'df_sigma', 'mu_tau', 'Sigma_tau')))
        {
            stop("\n\nThe parameter \'fix\' must be a named list with the parameters of the model. The parametters are 'mu_beta', 'Sigma_beta', 'alpha', 's2_sigma', 'df_sigma' \n\n. If there are context-level covariates, the additional parameters are required: 'mu_tau' and 'Sigma_tau'\n\n")
        }
    }
    if (family=='binomial')
    {
        if (any(! names(fix) %in% c('mu_beta', 'Sigma_beta', 'alpha')))
        {
            stop("\n\nThe parameter \'fix\' must be a named list with the parameters of the model. The parametters are 'mu_beta', 'Sigma_beta', 'alpha'\n\n. If there are context-level covariates, the additional parameters are required: 'mu_tau' and 'Sigma_tau'\n\n")
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
    if (!is.null(Dw) & !is.vector(fix$mu_tau))
    {
        stop("\n\n fix$mu_tau must be a vector\n\n")
    }
    if (!is.null(Dw) & !is.matrix(fix$Sigma_tau))
    {
        stop("\n\n fix$Sigma_tau must be a matrix\n\n")
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
    if(!is.null(Dw))
    {
        if (length(fix$mu_tau)!=Dw+1)
        {
            stop(paste0("Dimension of fix$mu_tau must be the same as the number of context-level covariates plus one (for the intercept), that is", Dw+1, sep=''))
        }
        if (!all(dim(fix$Sigma_beta)==c(d+1,d+1)) )
        {
            stop(paste0("The dimension of fix$Sigma_tau must be ", Dw+1,' x ', Dw+1,  sep=''))
        }

    }

}

getContextIndices <- function(W)
{
    return (
        W %>%
        tibble::as_data_frame(.) %>%  
        dplyr::select(-dplyr::contains('Intercept')) %>% 
        dplyr::group_by_(., names(.)) %>%
        dplyr::group_indices(.)
    )

}
getUniqueW <- function(W, C)
{
    return(
        base::cbind(W,C) %>%
        tibble::as_data_frame(.)  %>%
        dplyr::filter(!base::duplicated(W)) %>%
        dplyr::arrange(C) %>% 
        dplyr::select(-C) %>%
        as.matrix
    )
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
#'             and \code{formula2}. Note: it is advisable to scale the variables before the estimation
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
#' @param imp.bin string, either "R" or "Cpp" indicating the language of the implementation of the binomial model.
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
#' mcmc    = list(burn.in = 0,  n.iter = 1000)
#' samples = hdpGLM(y~., data=data$data, mcmc=mcmc, family='gaussian', n.display=30, K=50)
#'
#' summary(samples)
#' summary(samples, true.beta=summary(data)$beta)
#' 
#' plot(samples)
#' 
#' plot(samples, separate=TRUE)
#'  
#' @export

## }}}
hdpGLM <- function(formula1, formula2=NULL, data, weights=NULL, mcmc, K=100, fix=NULL, family='gaussian', epsilon=0.01, leapFrog=40, n.display=1000, hmc_iter=1, imp.bin="R")
{
    if(! family %in% c('gaussian', 'binomial', 'multinomial'))
        stop(paste0('Error: Parameter -family- must be a string with one of the following options : \"gaussian\", \"binomial\", or \"multinomial\"'))

    ## Debug/Monitoring message --------------------------
    msg <- paste0('\n\n','Starting Estimation ...',  '\n\n'); cat(msg)
    ## ---------------------------------------------------

    ## ## construct the regression matrices (data.frames) based on the formula provided
    ## ## -----------------------------------------------------------------------------
    func.call <- match.call(expand.dots = FALSE)
    mat     = .getRegMatrix(func.call, data, weights, formula_number=1)
    y       = mat$y
    X       = mat$X
    weights = ifelse(!is.null(mat$w), mat$w, rep(1,nrow(X)))
    ## Hierarchical covars
    if(is.null(formula2)){
        W = NULL
        C = NULL
    }else{
        W = .getRegMatrix(func.call, data, weights, formula_number=2)$X
        C = getContextIndices(W)
        W = getUniqueW(W, C)
    } 

    ## get constants
    ## -------------
    d   = ncol(X)  - 1                        # d is the number of covars, we subtract the intercept as X has a column with ones
    Dw  = unlist( ifelse(is.null(W), list(NULL), list(ncol(W) - 1)) ) # list and unlist only b/c ifelse() do not allow to return NULL

    if (is.null(fix)) {
        fix = .hdpGLM_get_constants(family=family, d=d, Dw=Dw) # list with values of the parameters of the priors/hyperpriors
    }else {
        .dphGLM_check_constants(family=family, d=d, Dw=Dw, fix=fix)
    }
    
    ## get the samples from posterior
    ## ------------------------------
    T.mcmc  = Sys.time()
    ## Debug/Monitoring message --------------------------
    msg <- paste0('\n\n','Estimation in progress ...',  '\n\n'); cat(msg)
    ## ---------------------------------------------------
    if (is.null(W))
    {
        if (family=='binomial') {
            if (imp.bin=="R")   samples =  dpGLM_mcmc_xxr( y, X, weights, K, fix,  family, mcmc, epsilon, leapFrog, n.display, hmc_iter)
            if (imp.bin=="cpp") samples =  dpGLM_mcmc    ( y, X, weights, K, fix,  family, mcmc, epsilon, leapFrog, n.display, hmc_iter) 
        }else{
            samples        =  dpGLM_mcmc( y, X,       weights, K, fix,  family, mcmc, epsilon, leapFrog, n.display, hmc_iter) #
        }
    }else
    {
        samples        =  hdpGLM_mcmc(y, X, W, C, weights, K, fix, family, mcmc, epsilon, leapFrog, n.display, hmc_iter)
    }
    T.mcmc  = Sys.time() - T.mcmc


    ## including colnames and class for the output
    ## -------------------------------------------
    if (is.null(W)){
        if(family=='gaussian'){
            colnames(samples$samples)  <- c('k', paste0('beta',1:(ncol(samples$samples)-2),  sep=''),'sigma')
        }else{
            colnames(samples$samples)  <- c('k', paste0('beta',1:(ncol(samples$samples)-1),  sep=''))
        }
    }else{
        if(family=='gaussian'){
            colnames(samples$samples)  <- c('k', 'j', paste0('beta',1:(d+1),  sep=''), 'sigma')
            colnames(samples$tau)      <- paste0(rep(paste0('tau', 1:(Dw+1)), d+1), unlist(lapply(1:(d+1), function(d) rep(d, Dw+1))),sep='')
        }else{
            colnames(samples$samples)  <- c('k', 'j', paste0('beta',1:(d+1),  sep=''))
            colnames(samples$tau)      <- paste0(rep(paste0('tau', 1:(Dw+1)), d+1), unlist(lapply(1:(d+1), function(d) rep(d, Dw+1))),sep='')
        }
    }

    samples$samples                   = coda::as.mcmc(samples$samples)
    if (!is.null(W)){
        samples$tau                       = coda::as.mcmc(samples$tau)
        samples$context.index             = C
        samples$context.cov               = tibble::as_data_frame(cbind(C=sort(unique(C)), W))  %>% dplyr::select(-dplyr::contains("Intercept")) 
    }
    if(family=='gaussian') attr(samples$samples, 'terms') = data.frame(term = c(colnames(X), 'sigma') , Parameter=c(stringr::str_subset(colnames(samples$samples), pattern="beta"), "sigma") )
    if(family=='binomial') attr(samples$samples, 'terms') = data.frame(term = c(colnames(X)         ) , Parameter=c(stringr::str_subset(colnames(samples$samples), pattern="beta")) )
    attr(samples$samples, 'mcpar')[2] = mcmc$n.iter
    class(samples)                    = ifelse(is.null(W), 'dpGLM', 'hdpGLM')

    samples$time_elapsed = T.mcmc
    return(samples)
}

