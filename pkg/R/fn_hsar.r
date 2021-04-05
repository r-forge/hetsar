#-- wrapper: write the estimator as a function of the dataset, only --#
fn_hsar <- function(l_data, method = "L-BFGS-B", l_bounds = NULL, ...) {
     # -------------------------------------------------------------------------- #
     # l_data is a list containing:
     # m_y: (N,T) matrix
     # l_x: list with K-1 (N,T)-matrices of covariates (it exclude the intercept)
     # m_W: (N,N) matrix

     # -------------------------------------------------------------------------- #
     N  <- nrow(l_data[["m_y"]])
     TT <- ncol(l_data[["m_y"]])

     # add y*:=Wy
     l_data[["m_ys"]] <- l_data[["m_W"]] %*% l_data[["m_y"]] # (N,T)
    ## adding intercept to be suppressed ####
    ## (keep the original one if any) #######
     # add an intercept at the **end** of the list of covariates
     Km1 <- length(l_data[["l_x"]]) # mnemonic: "K minus 1"
     K <- Km1 + 1
     l_data[["l_x"]][[K]] <- matrix(1L, N, TT)
     #print(lobstr::ref(l_data))
     #browser()
     #stop()
    #########################################
     # -------------------------------------------------------------------------- #
     # constructs v_theta_ini
     l_theta_ini <- NULL #cc#
     v_theta_ini <- vector("numeric", N * (K + 2))
     if (is.null(l_theta_ini)) {
          v_theta_ini[1:(N + (K * N))] <- 0
          v_theta_ini[(N + (K * N) + 1):(N + (N * K) + N)] <- 1
     } else {
          v_psi_ini   <- l_theta_ini[["psi"  ]]
          m_beta_ini  <- l_theta_ini[["beta" ]]
          v_sgmsq_ini <- l_theta_ini[["sqmsq"]]

          ## reshape m_beta_ini
          m_beta_ini_tr <- t(m_beta_ini) # (K,N)
          v_beta_ini <- c(m_beta_ini_tr) # (NK,)

          ## fill v_theta_ini
          v_theta_ini[1:N]                                 <- v_psi_ini   # (N,)
          v_theta_ini[(N + 1):(N + (K * N))]               <- v_beta_ini  # (KN,)
          v_theta_ini[(N + (K * N) + 1):(N + (N * K) + N)] <- v_sgmsq_ini # (N,)
     }

     # -------------------------------------------------------------------------- #
     # bounds
     if (is.null(l_bounds)) {
          v_lb <- c(rep(-0.995, N), rep(-Inf, K * N), rep(0.001, N)) # (N(K+2),)
          v_ub <- c(rep( 0.995, N), rep( Inf, (K * N) + N))
     } else {
          v_lb_psi   <- l_bounds[["lb_psi"  ]] # (N,)
          m_lb_beta  <- l_bounds[["lb_beta" ]]
          v_lb_sgmsq <- l_bounds[["lb_sqmsq"]]

          v_ub_psi   <- l_bounds[["ub_psi"  ]]
          m_ub_beta  <- l_bounds[["ub_beta" ]]
          v_ub_sgmsq <- l_bounds[["ub_sqmsq"]]

          ## add intercept
          m_lb_beta <- cbind(m_lb_beta, -Inf) # (K,N)
          m_ub_beta <- cbind(m_ub_beta,  Inf) # (K,N)


          ## reshape m_lb_beta
          m_lb_beta_tr <- t(m_lb_beta) # (K,N)
          m_ub_beta_tr <- t(m_ub_beta) # (K,N)
          v_lb_beta <- c(m_lb_beta_tr) # (NK,)
          v_ub_beta <- c(m_ub_beta_tr) # (NK,)

          ## fill v_lb and v_ub
          v_lb <- vector("numeric", N * (K + 2))
          v_ub <- vector("numeric", N * (K + 2))

          v_lb[1:N]                                 <- v_lb_psi   # (N,)
          v_lb[(N + 1):(N + (K * N))]               <- v_lb_beta  # (KN,)
          v_lb[(N + (K * N) + 1):(N + (N * K) + N)] <- v_lb_sgmsq # (N,)

          v_ub[1:N]                                 <- v_ub_psi   # (N,)
          v_ub[(N + 1):(N + (K * N))]               <- v_ub_beta  # (KN,)
          v_ub[(N + (K * N) + 1):(N + (N * K) + N)] <- v_ub_sgmsq # (N,)
          #browser()
     }

     # -------------------------------------------------------------------------- #
     switch (method,
          "L-BFGS-B" = {
               l_out <- optim(
                    par     = v_theta_ini, 
                    fn      = fn_ofun,
                    gr      = fn_gr,
                    l_data  = l_data,
                    method  = "L-BFGS-B",
                    lower   = v_lb,
                    upper   = v_ub,
                    control = list(maxit = 1e6, ...
                         #, trace = TRUE
                         ))
               v_theta <- l_out[["par"]]
          },
          "nlminb" = {
               l_out <- nlminb(
                    start     = v_theta_ini,
                    objective = fn_ofun,
                    gradient  = fn_gr,
                    lower     = v_lb,
                    upper     = v_ub,
                    l_data    = l_data)
               v_theta <- l_out[["par"]]
          },
          "trustOptim" = {
               # I don't seem to be able to suppress the output
               #sink("/dev/null")
               l_out <- trust.optim(
                    v_theta_ini,
                    fn           = fn_ofun,
                    gr           = fn_gr,
                    hs           = fn_hessian,
                    method       = "Sparse",
                    control = list(
                        maxit = 5000,
                        report.level = 0), # do not display any information (default 2)
                    l_data       = l_data)
               #sink()
               v_theta <- l_out[["solution"]]
          }#,
          #"trust" = {
          #     l_out <- trust::trust(
          #          objfun  = fn_ofun_gr_hessian,
          #          parinit = v_theta_ini,
          #          rinit   = 1, #!! rinit??
          #          rmax    = 5, #!! rmax??
          #          l_data  = l_data)
          #     v_theta <- l_out[["argument"]]
          #}
     )
     #print(l_out)

     v_psi   <- v_theta[1:N] # (N,)
     v_beta  <- v_theta[(N + 1):(N + (K * N))]; # (KN,)
     v_sgmsq <- v_theta[(N + (K * N) + 1):(N + (N * K) + N)] # (N,)
     #print(range(v_sgmsq))

     # reshape beta
     m_beta_tr <- matrix(v_beta, K, N) # (K,N)
     m_beta <- t(m_beta_tr) # (N,K)

     # pack all estimates into a matrix with names
     m_theta <- matrix(NA_real_, N, K + 2)
     m_theta[, 1]         <- v_psi
     m_theta[, 2:(K + 1)] <- m_beta
     m_theta[, (K + 2)]   <- v_sgmsq
     #browser()
     #print(m_theta)

     # -------------------------------------------------------------------------- #
     # variance
     l_return <- fn_var_ml(m_theta, l_data) # (2,)

     # -------------------------------------------------------------------------- #
     # return
     l_return[["theta"]] <- m_theta # (3,)
     l_return
}
