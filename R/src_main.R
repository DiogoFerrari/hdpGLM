
get_mcmc_epsilon <- function(mcmc)
{
    if ('epsilon' %in% names(mcmc)) {
        epsilon=epsilon
    }else {
        epsilon=0.01
    } 
    return(epsilon)
}
get_mcmc_leapFrog <- function(mcmc)
{
    if ('leapFrog' %in% names(mcmc)) {
        leapFrog=leapFrog
    }else {
        leapFrog=40
    } 
    return(leapFrog)
}
get_mcmc_hmc_iter <- function(mcmc)
{
    if ('hmc_iter' %in% names(mcmc)) {
        hmc_iter=hmc_iter
    }else {
        hmc_iter=1
    } 
    return(hmc_iter)
}

.hdpGLM_get_constants      <- function(family, d, Dw)
{
    Dw = ifelse(is.null(Dw), 0, Dw)
    if(family=='gaussian') {        
        sigma_beta = 10
        sigma_tau  = 5
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
    ## if (class(fix) != 'list')
    if (!methods::is(fix, 'list'))
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
    ## options(warn=-1)
    ## on.exit(options(warn=0))
    W = W %>% tibble::as_tibble() 
    cols = names(W)
    C = W %>%  
        dplyr::distinct(.) %>%
        dplyr::mutate(C = 1:nrow(.))  %>%
        dplyr::right_join(., W, by=cols)  %>%
        dplyr::select(C)  %>%
        dplyr::pull(.)
    return (C)

}
getUniqueW <- function(W, C)
{
    return(
        base::cbind(W,C) %>%
        tibble::as_tibble(.) %>%
        dplyr::filter(!base::duplicated(W)) %>%
        dplyr::arrange(C) %>% 
        dplyr::select(-C) %>%
        as.matrix
    )
}
hdpglm_exclude_nas <- function(data, formula1, formula2, context.id)
{
    ## vars = c(all.vars(formula1), all.vars(formula2))
    if (is.null(formula2)) {
        vars = c(formula.tools::get.vars(formula1, data=data))
    }else{
        vars = c(formula.tools::get.vars(formula1, data=data), formula.tools::get.vars(formula2, data=data))
    }
    if (!is.null(context.id)) vars = c(vars, context.id)
    return(data %>% dplyr::select(vars)  %>% dplyr::filter(stats::complete.cases(.)))
}

## {{{ doc }}}

#' Fit Hierarchical Dirichlet Process GLM
#'
#' The function estimates a semi-parametric mixture of Generalized
#' Linear Models. It uses a (hierarchical) Dependent Dirichlet Process
#' Prior for the mixture probabilities. 
#'
#' @param formula1 a single symbolic description of the linear model of the
#'                 mixture GLM components to be fitted. The syntax is the same
#'                 as used in the \code{\link{lm}} function.
#' @param formula2 eihter NULL (default) or a single symbolic description of the
#'                 linear model of the hierarchical component of the model.
#'                 It specifies how the average parameter of the base measure
#'                 of the Dirichlet Process Prior varies linearly as a function
#'                 of group level covariates. If \code{NULL}, it will use
#'                 a single base measure to the DPP mixture model.
#' @param data a data.frame with all the variables specified in \code{formula1} and \code{formula2}. Note: it is advisable to scale the variables before the estimation
#' @param mcmc a named list with the following elements
#' 
#'             - \code{burn.in} (required): an integer greater or equal to 0
#'                              indicating the number iterations used in the
#'                              burn-in period of the MCMC.
#' 
#'             - \code{n.iter}  (required): an integer greater or equal to 1
#'                              indicating the number of iterations to record
#'                              after the burn-in period for the MCMC.
#' 
#'             - \code{epsilon} (optional): a positive number. Default is 0.01.
#'                               Used when \code{family='binomial'} or
#'                              \code{family='multinomial'}. It is used in the
#'                               Stormer-Verlet Integrator (a.k.a leapfrog
#'                               integrator) to solve the Hamiltonian Monte
#'                               Carlo in the estimation of the model.
#'                
#'             - \code{leapFrog} (optional) an integer. Default is 40. Used when
#'                               \code{family='binomial'} or \code{family='multinomial'}.
#'                               It indicates the number of steps taken at each
#'                               iteration of the Hamiltonian Monte Carlo for
#'                               the Stormer-Verlet Integrator.
#' 
#'             - \code{hmc_iter} (optional) an integer. Default is 1. Used when 
#'                               \code{family='binomial'} or \code{family='multinomial'}.
#'                               It indicates the number of HMC iteration(s)
#'                               for each Gibbs iteration.
#'                 
#' @param family a character with either 'gaussian', 'binomial', or 'multinomial'.
#'               It indicates the family of the GLM components of the mixture model.
#' @param K an integer indicating the maximum number of clusters to truncate the
#'          Dirichlet Process Prior in order to use the blocked Gibbs sampler.
#' @param constants either NULL or a list with the constants of the model. If not NULL,
#'            it must contain a vector named \code{mu_beta}, whose size must be
#'            equal to the number of covariates specified in \code{formula1}
#'            plus one for the constant term; \code{Sigma_beta}, which must be a squared
#'            matrix, and each dimension must be equal to the size of the vector \code{mu_beta}; 
#'            and \code{alpha}, which must be a single number. If @param family is 'gaussian',
#'            then it must also contain \code{s2_sigma} and \code{df_sigma}, both
#'            single numbers. If NULL, the defaults are \code{mu_beta=0},
#'            \code{Sigma_beta=diag(10)}, \code{alpha=1}, \code{df_sigma=10},
#'            \code{s2_sigma=10} (all with the dimension automatically set to the
#'            correct values).
#' @param context.id string with the name of the column in the data that uniquely identifies the contexts. If \code{NULL} (default) contexts will be identified by numerical indexes and unique context-level variables. The user is advised to pre-process the data to provide meaningful labels for the contexts to facilitate later visualization and analysis of the results.
#' @param weights numeric vector with the same size as the number of rows of the data. It must contain the weights of the observations in the data set. NOTE: FEATURE NOT IMPLEMENTED YET
#' @param n.display an integer indicating the number of iterations to wait before printing information
#'                  about the estimation process. If zero, it does not display any information.
#'                  Note: displaying informaiton at every iteration (n.display=1) may increase
#'                  the time to estimate the model slightly. 
#' @param imp.bin string, either "R" or "Cpp" indicating the language of the implementation of the binomial model.
#' @param na.action string with action to be taken for the \code{NA} values. (currently, only \code{exclude} is available)
#' 
#' @return The function returns a list with elements \code{samples}, \code{pik}, \code{max_active},
#'         \code{n.iter}, \code{burn.in}, and \code{time.elapsed}. The \code{samples} element
#'         contains a MCMC object (from \pkg{coda} package) with the samples from the posterior
#'         distribution. The \code{pik} is a \code{n x K} matrix with the estimated
#'         probabilities that the observation $i$ belongs to the cluster $k$
#' 
#'
#' @details
#' This function estimates a Hierarchical Dirichlet Process generalized
#' linear model, which is a semi-parametric Bayesian approach to regression
#' estimation with clustering. The estimation is conducted using Blocked Gibbs Sampler if the output
#' variable is gaussian distributed. It uses Metropolis-Hastings inside Gibbs if
#' the output variable is binomial or multinomial distributed.
#' This is specified using the parameter \code{family}. See:
#'
#' Ferrari, D. (2020). Modeling Context-Dependent Latent Effect Heterogeneity,
#' Political Analysis, 28(1), 20–46.
#' 
#' Ishwaran, H., & James, L. F., Gibbs sampling methods for stick-breaking priors,
#' Journal of the American Statistical Association, 96(453), 161–173 (2001). 
#'
#' Neal, R. M., Markov chain sampling methods for dirichlet process mixture models,
#' Journal of computational and graphical statistics, 9(2), 249–265 (2000).
#'
#' @examples
#'
#' 
#' ## Note: this example is for illustration. You can run the example
#' ## manually with increased number of iterations to see the actual
#' ## results, as well as the data size (n)
#' 
#' set.seed(10)
#' n = 300
#' data = tibble::tibble(x1 = rnorm(n, -3),
#'                       x2 = rnorm(n,  3),
#'                       z  = sample(1:3, n, replace=TRUE),
#'                       y  =I(z==1) * (3 + 4*x1 - x2 + rnorm(n)) +
#'                           I(z==2) * (3 + 2*x1 + x2 + rnorm(n)) +
#'                           I(z==3) * (3 - 4*x1 - x2 + rnorm(n)) 
#'                       ) 
#' 
#' 
#' mcmc    = list(burn.in = 0,  n.iter = 20)
#' samples = hdpGLM(y~ x1 + x2, data=data, mcmc=mcmc, family='gaussian',
#'                  n.display=30, K=50)
#' 
#' summary(samples)
#' plot(samples)
#' plot(samples, separate=TRUE)
#' 
#' ## compare with GLM
#' ## lm(y~ x1 + x2, data=data,  family='gaussian')
#' 
#'  
#' @export

## }}}
hdpGLM <- function(formula1, formula2=NULL, data, mcmc, family='gaussian', K=100,
                   context.id=NULL, constants=NULL, weights=NULL,
                   n.display=1000, na.action = "exclude", imp.bin="R")
{
    if (!is.null(weights)) {
        stop(paste0("\n\nNote: weights are not implemented ",
                    "yet. Leave \'weights=NULL.\'\n\n"))
    }
    ## other options
    ## par.default <- par(no.readonly = TRUE)
    ## on.exit(par(par.default), add=TRUE)
    ## keep all default options
    op.default <- options()
    on.exit(options(op.default), add=TRUE)
    ## keep current working folder on exit
    dir.default <- getwd()
    on.exit(setwd(dir.default), add=TRUE)
    ## no warning messages
    ## options(warn=-1)
    ## on.exit(options(warn=0))

    if(! family %in% c('gaussian', 'binomial', 'multinomial'))
        stop(paste0('Error: Parameter -family- must be a string with one of' ,
                    'the following options : \"gaussian\", \"binomial\",',
                    ' or \"multinomial\"'))

    ## MCMC tunning parameters
    ## -----------------------
    fix = constants
    epsilon  = get_mcmc_epsilon(mcmc)
    leapFrog = get_mcmc_leapFrog(mcmc)
    hmc_iter = get_mcmc_hmc_iter(mcmc)

        

    ## Debug/Monitoring message --------------------------
    msg <- paste0('\n\n','Preparing for estimation ...',  '\n\n'); cat(msg)
    ## ---------------------------------------------------

    ## make sure the data is not grouped
    ## ---------------------------------
    data = data %>% dplyr::ungroup(.) 

    ## Exclude NA values
    ## -----------------
    if (na.action == 'exclude') {
        data    = hdpglm_exclude_nas(data, formula1, formula2, context.id)
    }else{
        ## Debug/Monitoring message --------------------------
        msg <- paste0('\n','Note: only action currently available to deal ',
                      'with NA\'s is \'exclude\'',  '\n'); cat(msg)
        ## ---------------------------------------------------
    }

    ## ## construct the reg matrices (data.frames) based on the formula provided
    ## ## ----------------------------------------------------------------------
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
    J = length(unique(C))

    ## get constants
    ## -------------
    ## d is the number of covars; subtract the intercept as X has a col with 1's
    d   = ncol(X)  - 1      
    ## list and unlist only b/c ifelse() do not allow to return NULL
    Dw  = unlist( ifelse(is.null(W), list(NULL), list(ncol(W) - 1)) ) 

    if (is.null(fix)) {
        ## list with values of the parameters of the priors/hyperpriors
        fix = .hdpGLM_get_constants(family=family, d=d, Dw=Dw)
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
            if (imp.bin=="R")   samples = dpGLM_mcmc_xxr( y, X, weights, K, fix,
                                                         family, mcmc, epsilon,
                                                         leapFrog, n.display,
                                                         hmc_iter)
            if (imp.bin=="cpp") samples = dpGLM_mcmc   ( y, X, weights, K, fix,
                                                        family, mcmc, epsilon,
                                                        leapFrog, n.display,
                                                        hmc_iter) 
        }else{
            samples =  dpGLM_mcmc( y, X, weights, K, fix,  family, mcmc,
                                  epsilon, leapFrog, n.display, hmc_iter) 
        }
    }else
    {
        if (family=='binomial') {
            stop(paste0("\n\nHierarchical version of hdpGLM is not implemented",
                        " for binomial family yet. \n\nIn the current",
                        " implementation, you can use the binomial family",
                        " to estimate hdpGLM in a single context only. \n\n"))
        }else{
            samples =  hdpGLM_mcmc(y, X, W, C, weights, K, fix, family, mcmc,
                                   epsilon, leapFrog, n.display, hmc_iter)
        }
    }
    T.mcmc  = Sys.time() - T.mcmc


    ## including colnames and class for the output
    ## -------------------------------------------
    if (is.null(W)){
        if(family=='gaussian'){
            colnames(samples$samples)  <- c('k',
                                            paste0('beta',
                                                   1:(ncol(samples$samples)-2),
                                                   sep=''),'sigma')
        }else{
            colnames(samples$samples)  <- c('k',
                                            paste0('beta',
                                                   1:(ncol(samples$samples)-1),
                                                   sep=''))
        }
    }else{
        colnames(samples$tau)      <- paste0(rep(paste0('tau[', 0:(Dw)), d+1),
                                             "][",
                                             unlist(
                                                 lapply(0:(d),
                                                        function(d)
                                                            rep(d, Dw+1))),
                                             ']',sep='')
        if(family=='gaussian'){
            colnames(samples$samples)  <- c('k', 'j', paste0('beta[',0:(d), "]",
                                                             sep=''), 'sigma')
        }else{
            colnames(samples$samples)  <- c('k', 'j', paste0('beta[',0:(d), "]",
                                                             sep=''))
        }
    }

    samples$samples                   = coda::as.mcmc(samples$samples)
    if (!is.null(W)){
        samples$tau                       = coda::as.mcmc(samples$tau)
        samples$context.index             = C
        samples$context.cov               = (
          tibble::as_tibble(cbind(C=sort(unique(C)), W))
            %>% dplyr::select(-dplyr::contains("Intercept"))  
        )
        if (!is.null(context.id)) {
            context = cbind(C=C, data %>%
                                 dplyr::select(context.id) )  %>%
                dplyr::filter(!duplicated(.))
            samples$context.cov = samples$context.cov %>%
                dplyr::left_join(., context, by=c("C")) 
        }
    }
    if(family=='gaussian')
        attr(samples$samples, 'terms') = data.frame(
            term = c(colnames(X), 'sigma') ,
            Parameter=c(stringr::str_subset(colnames(samples$samples),
                                            pattern="beta"), "sigma") )
    if(family=='binomial')
        attr(samples$samples, 'terms') = data.frame(
            term = c(colnames(X)         ) ,
            Parameter=c(stringr::str_subset(colnames(samples$samples),
                                            pattern="beta")) )
    attr(samples$samples, 'mcpar')[2] = mcmc$n.iter
    class(samples)                    = ifelse(is.null(W), 'dpGLM', 'hdpGLM')
    ## include context-level term names
    if (methods::is(samples, 'hdpGLM')) {
        tau.name = samples$tau %>% colnames
        tau.idx  = tau.name  %>%
            stringr::str_extract(string=.,
                                 pattern="tau\\[.*\\]\\[") %>%
            stringr::str_replace_all(string=., pattern="tau|\\[|\\]",
                                     replacement="")  %>%
            as.integer 
        attr(samples$sample_pi_postMean, "dimensions") = c("pi", 'context')
        attr(samples$sample_pi_postMean, "dimnames") = list(1:K, 1:J)
        attr(samples$tau, "terms")  = data.frame(Parameter = tau.name,
                                                 tau.idx = tau.idx) %>%
            dplyr::mutate(
                       beta = paste0("beta",
                                     stringr::str_extract(
                                                  Parameter,
                                                  pattern="\\[[0-9]*\\]$") ) %>%
                           as.character ) %>% 
            dplyr::left_join(., attr(samples$samples, "terms") %>%
                                dplyr::rename(beta = Parameter)  %>%
                                dplyr::mutate(beta = as.character(beta))
                           , by=c('beta')) %>%
            dplyr::rename(term.beta = term) %>%
            dplyr::arrange(tau.idx)  %>% 
            dplyr::bind_cols(.,
                             data.frame(term.tau = lapply(colnames(W),
                                                          function(x)
                                                              rep(x, times=d+1) ) %>%
                                            unlist) ) 
    }

    samples$time_elapsed = T.mcmc
    attr(samples, 'formula1') = formula1
    attr(samples, 'formula2') = formula2
    attr(samples, 'family') = family
    samples$data=data
    return(samples)
}

