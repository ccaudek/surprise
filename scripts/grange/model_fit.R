### code to fit the DSTP and SSP models 


#------------------------------------------------------------------------------
### fit the DSTP model to individual data
fit_dstp_individual <- function(data){


  #---
  # data & model parameter preparation

  # first, change the header of the data frame to match that required by
  # the flankr package
  data <- data %>%
    rename(subject = id)

  # get the list of subjects
  subjects <- unique(data$subject)

  # set random seed so user can reproduce simulation outcome
  set.seed(42)

  # during the fit, how many sets of starting parameters should be used?
  n_start_parms <- 50

  # what should the variance across starting parameters be?
  var_start_parms <- 20

  # how many trials to simulate during each iteration of the fit routine whilst
  # exploring multiple starting parameters?
  n_first_pass <- 1000

  # how many trials to simulate during the final fit routine?
  n_final_pass <- 50000

  # prepare containers for model parameters
  dstp_fit <- matrix(0, nrow = length(subjects), ncol = 10)
  colnames(dstp_fit) <- c("id", "A", "C", "muT", "muFl", "muSS", "muRS2",
                          "ter", "G2", "bBIC")
  dstp_fit <- data.frame(dstp_fit)

  #---
  # do the model fitting

  # loop over all subjects
  for(i in 1:length(subjects)){

    print(i)

    # get the current subject's data
    subject_data <- data %>%
      filter(subject == subjects[i])

    # perform the first fit
    first_fit <- fitMultipleDSTP(subject_data, var = var_start_parms,
                                 nParms = n_start_parms,
                                 nTrials = n_first_pass,
                                 multipleSubjects = FALSE)

    # perform final fit, and store the parameter values
    final_fit <- fitDSTP(subject_data, parms = first_fit$bestParameters,
                         nTrials = n_final_pass, multipleSubjects = FALSE)
    dstp_fit[i, 1] <- as.character(subjects[i])
    dstp_fit[i, 2:8] <- final_fit$bestParameters
    dstp_fit[i, 9] <- final_fit$g2
    dstp_fit[i, 10] <- final_fit$bBIC

    # update the file to csv after each subject (in case of power outage)
    write.csv(dstp_fit, "dstp_fit.csv", row.names = FALSE)
  }

  return(final_fit)

}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### fit the SSP model to individual data
fit_ssp_individual <- function(data){


  #---
  # data & model parameter preparation

  # first, change the header of the data frame to match that required by
  # the flankr package
  data <- data %>%
    rename(subject = id)

  # get the list of subjects
  subjects <- unique(data$subject)

  # set random seed so user can reproduce simulation outcome
  set.seed(42)

  # during the fit, how many sets of starting parameters should be used?
  n_start_parms <- 50

  # what should the variance across starting parameters be?
  var_start_parms <- 20

  # how many trials to simulate during each iteration of the fit routine whilst
  # exploring multiple starting parameters?
  n_first_pass <- 1000

  # how many trials to simulate during the final fit routine?
  n_final_pass <- 50000

  # prepare containers for model parameters
  ssp_fit <- matrix(0, nrow = length(subjects), ncol = 8)
  colnames(ssp_fit) <- c("id", "A", "ter", "p", "rd", "sda","G2", "bBIC")
  ssp_fit <- data.frame(ssp_fit)

  #---
  # do the model fitting

  # loop over all subjects
  for(i in 1:length(subjects)){

    print(paste("ssp_", i, sep = ""))

    # get the current subject's data
    subject_data <- data %>%
      filter(subject == subjects[i])

    # perform the first fit
    first_fit <- fitMultipleSSP(subject_data, var = var_start_parms,
                                nParms = n_start_parms,
                                nTrials = n_first_pass,
                                multipleSubjects = FALSE)

    # perform final fit, and store the parameter values
    final_fit <- fitSSP(subject_data, parms = first_fit$bestParameters,
                        nTrials = n_final_pass, multipleSubjects = FALSE)
    ssp_fit[i, 1] <- as.character(subjects[i])
    ssp_fit[i, 2:6] <- final_fit$bestParameters
    ssp_fit[i, 7] <- final_fit$g2
    ssp_fit[i, 8] <- final_fit$bBIC

    # update the file to csv after each subject (in case of power outage)
    write.csv(ssp_fit, "ssp_fit.csv", row.names = FALSE)
  }

  return(final_fit)

}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### fit the DSTP model to individual data using parallel processing 
fit_dstp_individual_par <- function(data, n_start_parms, var_start_parms, 
                                    n_first_pass, n_final_pass){
  
  
  #--- 
  # data & model parameter preparation
  
  # first, change the header of the data frame to match that required by
  # the flankr package
  data <- data %>% 
    rename(subject = id)
  
  # get the list of subjects
  subjects <- unique(data$subject)
  
  # set random seed so user can reproduce simulation outcome
  set.seed(42)
  
  # prepare containers for model parameters
  dstp_fit <- matrix(0, nrow = length(subjects), ncol = 10)
  colnames(dstp_fit) <- c("id", "A", "C", "muT", "muFl", "muSS", "muRS2", 
                          "ter", "G2", "bBIC")
  dstp_fit <- data.frame(dstp_fit)
  
  #---
  # do the model fitting
  
  # loop over all subjects
  for(i in 1:length(subjects)){
    
    print(paste("DSTP_fit_", i, sep = ""))
    
    # get the current subject's data
    subject_data <- data %>% 
      filter(subject == subjects[i])
    
    # perform the first fit 
    first_fit <- explore_parameters_dstp(subject_data, var = var_start_parms, 
                                         n_start_parms = n_start_parms, 
                                         n_trials = n_first_pass)
    
    # perform final fit, and store the parameter values
    final_fit <- fitDSTP(subject_data, parms = first_fit$best_parameters, 
                         nTrials = n_final_pass, multipleSubjects = FALSE)
    dstp_fit[i, 1] <- as.character(subjects[i])
    dstp_fit[i, 2:8] <- final_fit$bestParameters
    dstp_fit[i, 9] <- final_fit$g2
    dstp_fit[i, 10] <- final_fit$bBIC
    
    # update the file to csv after each subject (in case of power outage)
    write.csv(dstp_fit, "dstp_fit.csv", row.names = FALSE)
  }
  
  return(final_fit)
  
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### fit the DSTP model to group data using parallel processing
fit_dstp_group_par <- function(data, n_start_parms, var_start_parms,
                               n_first_pass, n_final_pass){


  #---
  # data & model parameter preparation

  # first, change the header of the data frame to match that required by
  # the flankr package
  data <- data %>%
    rename(subject = id)

  # set random seed so user can reproduce simulation outcome
  set.seed(42)

  #---
  # do the model fitting

  # perform the first fit
  first_fit <- explore_parameters_group_dstp(data,
                                             var = var_start_parms,
                                             n_start_parms = n_start_parms,
                                             n_trials = n_first_pass)


  # perform final fit, and store the parameter values
  final_fit <- fitDSTP(subject_data, parms = first_fit$best_parameters,
                       nTrials = n_final_pass, multipleSubjects = TRUE)

  return(final_fit)

}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### DSTP function to explore a wider range of starting parameters to inform
### the best parameters to use in the final fit
explore_parameters_dstp <- function(data, var_start_parms, n_start_parms, 
                                    n_trials, 
                                    start_parms = c(0.145, 0.08, 0.10, 0.07, 
                                                    0.325, 1.30, 0.240),
                                    max_parms = c(1, 1, 1, 1, 1, 2, 1)){
   
  
  
  # get a matrix of starting parameters
  parameters <- get_random_parms(start_parms, var_start_parms, max_parms,
                                 n_start_parms)

  # explore these starting parameters in parallel

  num_cores <- detectCores()
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  

  # dstp_search <- for(i in 1:n_start_parms){
  #   print(i)
  #   do_dstp_parallel(parameters[i, ], data, n_trials)
  # }

  dstp_search <- foreach(i = 1:n_start_parms,
                         .export = "do_dstp_parallel",
                         .packages = "flankr") %dopar%{
    do_dstp_parallel(parameters[i, ], data, n_trials)
  }

  stopCluster(cl)

  # collate all data back into a matrix
  dstp_search <- matrix(unlist(dstp_search), ncol = ncol(parameters) + 2,
                        byrow = TRUE)

  colnames(dstp_search) <- c("a", "c", "drift_target", "drift_flanker",
                             "drift_ss", "drift_rs2", "ter", "g2", "bBIC")
  dstp_search <- data.frame(dstp_search)
  # choose the best fitting parameters
  best <- dstp_search %>%
    top_n(-1, bBIC)

  # return them
  best <- as.numeric(best)
  best <- list(best_parameters = best[1:7],
               g2 = best[8],
               bBIC = best[9])

  return(best)
}
#------------------------------------------------------------------------------


myFun2 <- function(dfs) {
  foreach(i = seq_along(dfs)) %dopar% {
    df <- dfs[[i]]
    count(df, Species)
  }
}


#------------------------------------------------------------------------------
### DSTP function to explore a wider range of starting parameters to inform
### the best parameters to use in the final fit
explore_parameters_group_dstp <- function(data, var_start_parms, n_start_parms,
                                          n_trials,
                                          start_parms = c(0.145, 0.08, 0.10,
                                                          0.07, 0.325, 1.30,
                                                          0.240),
                                          max_parms = c(1, 1, 1, 1, 1, 2, 1)){

  # get a matrix of starting parameters
  parameters <- get_random_parms(start_parms, var_start_parms, max_parms,
                                 n_start_parms)

  # explore these starting parameters in parallel

  num_cores <- detectCores()
  registerDoParallel(num_cores)


  dstp_search <- foreach(i = 1:n_start_parms) %dopar%{
    do_dstp_parallel_group(parameters[i, ], data, n_trials)
  }

  # collate all data back into a matrix
  dstp_search <- matrix(unlist(dstp_search), ncol = ncol(parameters) + 2,
                        byrow = TRUE)

  colnames(dstp_search) <- c("a", "c", "drift_target", "drift_flanker",
                             "drift_ss", "drift_rs2", "ter", "g2", "bBIC")
  dstp_search <- data.frame(dstp_search)
  # choose the best fitting parameters
  best <- dstp_search %>%
    top_n(-1, bBIC)

  # return them
  best <- as.numeric(best)
  best <- list(best_parameters = best[1:7],
               g2 = best[8],
               bBIC = best[9])

  return(best)
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### fit the SSP model to individual data using parallel processing 
fit_ssp_individual_par <- function(data, n_start_parms, var_start_parms, 
                                   n_first_pass, n_final_pass){
  
  
  #--- 
  # data & model parameter preparation
  
  # first, change the header of the data frame to match that required by
  # the flankr package
  data <- data %>% 
    rename(subject = id)
  
  # get the list of subjects
  subjects <- unique(data$subject)
  
  # set random seed so user can reproduce simulation outcome
  set.seed(42)
  
  # prepare containers for model parameters
  ssp_fit <- matrix(0, nrow = length(subjects), ncol = 8)
  colnames(ssp_fit) <- c("id", "a", "ter", "p", "rd", "sda", "g2", "bBIC")
  ssp_fit <- data.frame(ssp_fit)
  
  #---
  # do the model fitting
  
  # loop over all subjects
  for(i in 105:length(subjects)){
    
    print(paste("SSP_fit_", i, sep = ""))
    
    # get the current subject's data
    subject_data <- data %>% 
      filter(subject == subjects[i])
    
    # perform the first fit 
    first_fit <- explore_parameters_ssp(subject_data, var = var_start_parms, 
                                        n_start_parms = n_start_parms, 
                                        n_trials = n_first_pass)
    
    
    # perform final fit, and store the parameter values
    final_fit <- fitSSP(subject_data, parms = first_fit$best_parameters, 
                        nTrials = n_final_pass, multipleSubjects = FALSE)
    ssp_fit[i, 1] <- as.character(subjects[i])
    ssp_fit[i, 2:6] <- final_fit$bestParameters
    ssp_fit[i, 7] <- final_fit$g2
    ssp_fit[i, 8] <- final_fit$bBIC
    
    # update the file to csv after each subject (in case of power outage)
    write.csv(ssp_fit, "ssp_fit.csv", row.names = FALSE)
  }
  
  return(final_fit)
  
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### SSP function to explore a wider range of starting parameters to inform
### the best parameters to use in the final fit
explore_parameters_ssp <- function(data, var_start_parms, n_start_parms, 
                                   n_trials, 
                                   start_parms = c(0.05, 0.3, 0.4, 0.05, 1.5),
                                   max_parms = c(1, 1, 1, 1, 3)){
  
  # get a matrix of starting parameters
  parameters <- get_random_parms(start_parms, var_start_parms, max_parms, 
                                 n_start_parms)
  
  
  
  num_cores <- detectCores()
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  ssp_search <- foreach(i = 1:n_start_parms,
                        .export = "do_ssp_parallel",
                        .packages = "flankr") %dopar%{
                        do_ssp_parallel(parameters[i, ], data, n_trials)
                         }
  
  stopCluster(cl)
  
  # collate all data back into a matrix
  ssp_search <- matrix(unlist(ssp_search), ncol = ncol(parameters) + 2,
                        byrow = TRUE)
  
  colnames(ssp_search) <- c("a", "ter", "p", "rd", "sda", "g2", "bBIC")
  ssp_search <- data.frame(ssp_search)
  
  # choose the best fitting parameters 
  best <- ssp_search %>% 
    top_n(-1, bBIC) 
  
  # return them
  best <- as.numeric(best)
  best <- list(best_parameters = best[1:5], 
               g2 = best[6], 
               bBIC = best[7])
  
  return(best)
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### function called during parallel fit of the dstp model
do_dstp_parallel <- function(start_parms, data, n_trials){
  
  fit <- fitDSTP(data, parms = start_parms, nTrials = n_trials,
                 multipleSubjects = FALSE)

  fit_as_vector <- c(fit$bestParameters, fit$g2, fit$bBIC)

  return(fit_as_vector)

}
#------------------------------------------------------------------------------



# #------------------------------------------------------------------------------
# ### function called during parallel fit of the dstp model
# do_dstp_parallel_group <- function(start_parms, data, n_trials){
#   
#   fit <- fitDSTP(data, parms = start_parms, nTrials = n_trials, 
#                  multipleSubjects = TRUE)
#   
#   fit_as_vector <- c(fit$bestParameters, fit$g2, fit$bBIC)
#   
#   return(fit_as_vector)
#   
# }
# #------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### function called during parallel fit of the ssp model
do_ssp_parallel <- function(start_parms, data, n_trials){
  
  fit <- fitSSP(data, parms = start_parms, nTrials = n_trials, 
                multipleSubjects = FALSE)
  
  fit_as_vector <- c(fit$bestParameters, fit$g2, fit$bBIC)
  
  return(fit_as_vector)
  
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### function to generate a matrix of random starting parameters
get_random_parms <- function(start_parms, var_start_parms, max_parms, 
                             n_start_parms){
  
  # initialise the matrix
  final_parms <- matrix(0, nrow = n_start_parms, ncol = length(start_parms))
  
  
  # get the desired SD to explore around the starting parameters
  var_start_parms <- (start_parms / 100) * var_start_parms
  
  # for each row...
  for(i in 1:n_start_parms){
    
    # generate the parameters
    curr_parms <- start_parms + rnorm(length(start_parms), 0, var_start_parms)
    
    # make sure no parameter falls outside <0 and >max_parms
    while(((min(curr_parms) < 0 | min(max_parms - curr_parms) < 0))) {
      
      curr_parms <- start_parms + rnorm(length(start_parms), 0, 
                                        var_start_parms)
      
    }
    
    # store the current vector
    final_parms[i, ] <- curr_parms
    
  }
  
  # change the first entry to match the starting parameters
  final_parms[1, ] <- start_parms
  
  # round them to 3 decimal places
  final_parms <- round(final_parms, 3)
  
  # return the matrix
  return(final_parms)
  
}
#------------------------------------------------------------------------------