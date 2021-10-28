#' Search for rigth seeds for the rigth AR simulation with arima.sin() finction using auto.arima() function
#'
#' More detailed Description
#'
#' @describeIn This ARSearch helps to Search for rigth seeds for the rigth AR simulation with arima.sin() finction using auto.arima() function
#'
#' @importFrom forecast forecast auto.arima accuracy
#'
#' @param a first seed boundary
#'
#' @param z last seed boundary
#'
#' @param n number of samples
#'
#' @param ar11 first coefficient of autoregressive
#'
#' @param ar22 second coefficient of autoregressive
#'
#' @param ar33 third coefficient of autoregressive
#'
#' @param p order of the autoregressive
#'
#' @param d degree of difference
#'
#' @param  q degree of moving average
#'
#' @param sd standard deviation of the series
#'
#' @param j1 length of character to search for in first coefficient of autoregressive
#'
#' @param j2 length of character to search for in second coefficient of autoregressive
#'
#' @param j3 length of character to search for in second coefficient of autoregressive
#'
#' @param arr1 character to search for in first coefficient of autoregressive
#'
#' @param arr2 character to search for in second coefficient of autoregressive
#'
#' @param arr3 character to search for in third coefficient of autoregressive
#'
#' @importFrom foreach `%dopar%` foreach
#'
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#'
#' @importFrom future plan multisession
#'
#' @importFrom tibble tibble
#'
#' @importFrom stats arima.sim
#'
#' @importFrom simEd set.seed
#'
#' @example
#' arsearch(a = 280000,  z = 290000, n = 10, p = 1, d = 0, q = 0, ar11 = 0.8, ar22 = 0.1, ar33 = 0.3, sd = 1, j1 = 4, j2 = 4, j3 = 4, arr1 = "0.80", arr2 = "0.30", arr3 = "0.30")
#'
#' @export
arsearch <- function(a, z, n, ar11, ar22, ar33, p, d, q, sd = 1, j1, j2, j3, arr1, arr2, arr3){

  output <- if (p == 1) {
    ar1_sim_search <- function(a, z, n, ar11, p, d, q, sd = 1, j1, arr1){
      future::plan(future::multisession)
      n_cores <- parallel::detectCores()
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cores = n_cores)

      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = c(ar11), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        } else if (all(grepl(c("ar1|intercept"), names(cf))) &
                   substr(cf["ar1"], 1, j1) %in% arr1) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]

      parallel::stopCluster(cl)
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ar1_sim_search(a = a,  z = z, n = n, p = 1, d = d, q = q, ar11 = ar11, sd = sd, j1 = j1, arr1 = arr1)

  } else if (p == 2) {

    #####################################################################################################

    ar2_sim_search <- function(a, z, n, ar11, ar22, p, d, q, sd = 1, j1, j2, arr1, arr2){
      future::plan(future::multisession)
      n_cores <- parallel::detectCores()
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cores = n_cores)

      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = c(ar11, ar22), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        }  else if (all(grepl(c("ar1|ar2|intercept"), names(cf))) &
                    substr(cf["ar1"], 1, j1) %in% arr1 & substr(cf["ar2"], 1, j2) %in% arr2) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]

      parallel::stopCluster(cl)
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ar2_sim_search(a = a,  z = z, n = n, p = 2, d = d, q = q, ar11 = ar11, ar22 = ar22, sd = sd, j1 = j1, j2 = j2, arr1 = arr1, arr2 = arr2)

  } else {

    #####################################################################################################

    ar3_sim_search <- function(a, z, n, ar11, ar22, ar33, p, d, q, sd = 1, j1, j2, j3, arr1, arr2, arr3){
      future::plan(future::multisession)
      n_cores <- parallel::detectCores()
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cores = n_cores)

      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        simEd::set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = c(ar11, ar22, ar33), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        } else if (all(grepl(c("ar1|ar2|ar3|intercept"), names(cf))) &
                   substr(cf["ar1"], 1, j1) %in% arr1 & substr(cf["ar2"], 1, j2) %in% arr2 & substr(cf["ar3"], 1, j3) %in% arr3) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]

      parallel::stopCluster(cl)
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ar3_sim_search(a = a,  z = z, n = n, p = 3, d = d, q = q, ar11 = ar11, ar22 = ar22, ar33 = ar33, sd = sd, j1 = j1, j2 = j2, j3 = j3, arr1 = arr1, arr2 = arr2, arr3 = arr3)

  }
  ##############################################################################

  return(output)
}
