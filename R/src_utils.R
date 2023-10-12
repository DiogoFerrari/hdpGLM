# required (by devtools) to link the cpp code 
#' @useDynLib hdpGLM
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
#' @importFrom Hmisc wtd.quantile
#' @importFrom stats median
#' @importFrom stats sd
#' @importFrom utils head
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
NULL
 

.get_layout         <-function(n, bycol=T){
    if(bycol){
        ncols = ceiling(sqrt(n))
        nlines  = floor(sqrt(n))
        m=ifelse(n - nlines*ncols>0, 1,0)
        grid=c(nlines=nlines+m, ncols=ncols)
    }else{
        nlines = ceiling(sqrt(n))
        ncols  = floor(sqrt(n))
        m=ifelse(n - nlines*ncols>0, 1,0)
        grid=c(nlines=nlines, ncols=ncols+m)
    }
    return(grid)
}
.getRegMatrix <- function(func.call, data, weights, formula_number=1)
{
    args <- names(func.call)
    ## creating the dependent variable and the covariates matrix from the fomula 1
    f = paste0('formula', formula_number, sep='')
    idx.args  <- match(c(f,  "data", "weights"), args, 0L)
    func.call <- func.call[c(1L, idx.args)]
    names(func.call)[names(func.call)==f] = "formula"
    func.call$drop.unused.levels <- TRUE
    func.call[[1L]] <- quote(stats::model.frame)
    func.call[[3]] = quote(data)
    reg.matrix <- eval(func.call, parent.frame())
    ## response variable
    y   <- stats::model.response(reg.matrix, "numeric")
    ## weights
    w   <- as.vector(stats::model.weights(reg.matrix))
    if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector")
    offset <- as.vector(stats::model.offset(func.call))
    if (!is.null(offset)) {
        if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", length(offset), NROW(y)), domain = NA)
    }
    ## covariates
    mt1    <- attr(reg.matrix, "terms")
    if (stats::is.empty.model(mt1)) {
        x <- matrix(1, ncol=1,nrow=nrow(y))
        results <- list(coefficients = if (is.matrix(y)) matrix(, 0, 3) else numeric(), residuals = y, fitted.values = 0 * y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w != 0) else if (is.matrix(y)) nrow(y) else length(y))
        if (!is.null(offset)) {
            results$fitted.values <- offset
            results$residuals <- y - offset
        }
    } else {
        x <- stats::model.matrix(mt1, reg.matrix)
    }
    return(list(y=y, X=x, w=w))
}
.hdpGLM_select_numerical   <- function(data)
{
    idxChar <- sapply(sapply(data, class), function(x) x[[1]]) == 'factor' |
        sapply(sapply(data, class), function(x) x[[1]]) == 'character'  |
        sapply(sapply(data, class), function(x) x[[1]]) == 'ordered'
    return(data.frame(data[,!idxChar]))
}
.hdpGLM_select_categorical <- function(data)
{
    idxChar <- sapply(sapply(data, class), function(x) x[[1]]) == 'factor' |
        sapply(sapply(data, class), function(x) x[[1]]) == 'character'  |
        sapply(sapply(data, class), function(x) x[[1]]) == 'ordered'
    return(data.frame(data[,idxChar]))
}


## Startup and closing
## -------------------
.onAttach<- function(libname, pkgname)
{
    packageStartupMessage(paste0('
## ===============================================================
## Hierarchial Dirichlet Process Generalized Linear Model (hdpGLM)
## ===============================================================

Author: Diogo Ferrari
Usage : https://github.com/DiogoFerrari/hdpGLM

'))
}
.onUnload <- function (libpath)
{
    library.dynam.unload("hdpGLM", libpath)
}


## Math functions
## --------------
edist <- function(v1,v2)
{
    return(sqrt(sum( (v1 - v2)^2 )))
}


## Colors for Plot
## ---------------
addalpha <- function(colors, alpha=1.0) {
    r <- grDevices::col2rgb(colors, alpha=T)
                                        # Apply alpha
    r[4,] <- alpha*255
    r <- r/255.0
    return(grDevices::rgb(r[1,], r[2,], r[3,], r[4,]))
}
colorRampPaletteAlpha <- function(colors, n=32, interpolate='linear') {
                                        # Create the color ramp normally
    cr <- grDevices::colorRampPalette(colors, interpolate=interpolate)(n)
                                        # Find the alpha channel
    a <- grDevices::col2rgb(colors, alpha=T)[4,]
                                        # Interpolate
    if (interpolate=='linear') {
        l <- stats::approx(a, n=n)
    } else {
        l <- stats::spline(a, n=n)
    }
    l$y[l$y > 255] <- 255 # Clamp if spline is > 255
    cr <- addalpha(cr, l$y/255.0)
    return(cr)
}

## global varibles for dplyr (used only so that the check ignores it,
## it does not actually creates global variables)
## -------------------------
if(getRversion() >= "2.15.1")  utils::globalVariables(
                                          c(".",
                                            "<<-",
                                            ## "txtProgressBar",
                                            ## "setTxtProgressBar",
                                            "mh.accept.info" ,
                                            "parameter" ,
                                            "Cluster" ,
                                            "Percentage.of.Iter..Cluster.was.active" ,
                                            "HPD.l" ,
                                            "HPD.u",
                                            "HPD.lower" ,
                                            "HPD.upper",
                                            "flag" ,
                                            "Parameter",
                                            "k",
                                            "j",
                                            "Mean",
                                            "Median",
                                            "True.Cluster.match",
                                            "True",
                                            "d",
                                            "values" ,
                                            "C",
                                            "sigma",
                                            "jnext",
                                            "covars",
                                            "X2.5."  ,
                                            "X97.5.",
                                            "SD",
                                            "dw",
                                            "dx",
                                            "term",
                                            "term.tau",
                                            "term.beta",
                                            "w",
                                            "Description",
                                            "E.beta.pred",
                                            "term.label",
                                            "Parameter.facet",
                                            "Parameter.label",
                                            "beta.label",
                                            "value",
                                            "nclusters",
                                            "tau.label",
                                            "w.label",
                                            "beta.exp",
                                            "beta.exp.values",
                                            "Parameter.tau",
                                            "beta.exp.values.tau",
                                            "True.beta",
                                            "tau.idx",
                                            "pk"
                                            ))
