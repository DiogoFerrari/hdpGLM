# hdpGLM (development version)

# hdpGLM 1.0.3

- This is a minor release with typos corrected and the inclusion of information regarding the publication of the package in the Journal of Statistical Software (10.18637/jss.v107.i10)

# hdpGLM 1.0.2

- The MCMC tunning parameters of the function `hdpGLM()`, namely  `epsilon`, `leapFrog`, and `hmc_iter` were removed from that function's signature. Now, when the users want to tune those parameters, they must provide it as part of the named list that is used for the argument `mcmc`. It makes the code and the function signature cleaner.

- The function `classify()` replaced the function `hdpGLM_classify()`. The latter will be removed in future versions.

- A function `nclusters()` was added. It returns the number of clusters estimated.

- A function txtProgressBar is now loaded in the namespace. 

# hdpGLM 1.0.1

* Various improvements for efficiency and speeding up the estimation.
