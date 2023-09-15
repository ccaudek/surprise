### function file

#------------------------------------------------------------------------------
### import the relevant data
get_data <- function(){
  
  # The data have been created with the 01_, 10_, 11_ scripts in 
  # the present directory.
  mydata <- read_csv(
    here("data", "processed", "surprise_control_exps.csv")
  )
  
  return(mydata)
}

#------------------------------------------------------------------------------
### create a random alpha-numeric code of a fixed length for each subject

generate_random_code <- function(length) {
  pool <- c(0:9, letters, LETTERS)
  paste0(sample(pool, length, replace = TRUE), collapse = "")
}

#------------------------------------------------------------------------------
### tidy the flanker data

tidy_flanker <- function(data){
  
  # Be sure that sequence of trial is in the right order
  df <- data %>%
    dplyr::arrange(subj_id, date, time_of_day)
  
  # Remove blocks with less than 20 trials.
  out <- df %>% 
    group_by(subj_name, block) %>% 
    summarise(
      ntrials_per_block = n()
    )
  # hist(out$ntrials_per_block)
  df1 <- left_join(df, out)
  df2 <- df1[df1$ntrials_per_block > 20, ]
  
  # Remove subjects with less than 100 trials.
  # add number of trials for each subject
  df2 <- df2 %>% 
    group_by(subj_name) %>% 
    mutate(
      bysub_n = length(rt)
    )
  df3 <- df2 %>%
    dplyr::filter(bysub_n > 100)
  
  # Remove block 5 and trial_number == 5
  # technical problems: there should be none.
  thedat <- df3 %>%
    dplyr::filter(
      trials_after_clip != "4" & block < 5
    )
  
  # Create a random alpha-numeric code of a fixed length for each subject
  set.seed(123) # Set seed for reproducibility
  # Get the unique subject IDs
  unique_subj_ids <- unique(thedat$subj_id)
  # Generate a unique random code for each unique subject ID
  unique_codes <- 
    setNames(as.list(sapply(unique_subj_ids, function(x) generate_random_code(10))), unique_subj_ids)
  # Create a new column with the random codes, matched to the subj_id column
  thedat$new_subj_id <- unlist(unique_codes[as.character(thedat$subj_id)])
  # length(unique(thedat$new_subj_id))
  thedat$subj_id <- thedat$new_subj_id
  thedat$new_subj_id <- NULL
  
  # Correct movie_id coding.
  part_to_remove <- "/Users/lab/Documents/surprise_project/surprise/stim/"
  result <- sub(part_to_remove, "", thedat$movie_id)
  
  part_to_remove <- "/Users/lab/Documents/magic3Definitivo/stim/"
  result1 <- sub(part_to_remove, "", result)
  
  part_to_remove <- "/Users/lab/Documents/magic_control_exp/stim/"
  result2 <- sub(part_to_remove, "", result1)
  
  part_to_remove <- "/Users/corrado/Dropbox/experiments/surprise_control_exp/stim/"
  result3 <- sub(part_to_remove, "", result2)
  
  part_to_remove <- "/Users/corrado/Documents/experiments/2019/surprise/stim/"
  result4 <- sub(part_to_remove, "", result3)
  
  part_to_remove <- "C:/Users/Miranda/Desktop/matlab_script/stim/"
  result5 <- sub(part_to_remove, "", result4)
  
  part_to_remove <- "/Users/lab/Documents/experiments/2019/surprise_5/stim/"
  result6 <- sub(part_to_remove, "", result5)
  
  part_to_remove <- "/Users/lab/Desktop/Surprise/stim/"
  result7 <- sub(part_to_remove, "", result6)
  
  thedat$movie_id <- result7
  
  # data wrangling
  thedat$experiment <- factor(thedat$experiment)
  thedat$subj_id <- factor(thedat$subj_id)
  thedat$resp <- factor(thedat$resp)
  thedat$movie_id <- factor(thedat$movie_id)
  thedat$is_surprise_clip <- factor(thedat$is_surprise_clip)
  thedat$is_clip_trial <- factor(thedat$is_clip_trial)
  
  return(thedat)
}

#------------------------------------------------------------------------------
### Get the accuracy per participant and condition 

get_flanker_accuracy <- function(data, overall = TRUE){
  
  flanker_data <- data[!is.na(data$correct), ]
  
  if(overall == TRUE){
    flanker_accuracy <- flanker_data %>% 
      group_by(experiment, subj_id) %>% 
      summarise(accuracy = sum(correct) / length(correct))
  } else {
    flanker_accuracy <- flanker_data %>% 
      group_by(experiment, subj_id, is_congruent_trial) %>% 
      summarise(accuracy = sum(correct) / length(correct))
  }
  
  return(flanker_accuracy)
  
}

#------------------------------------------------------------------------------
### calculate the accuracy congruency effect

get_flanker_congruency_accuracy <- function(flanker_data){
  
  flanker_congruent_accuracy <- flanker_data %>% 
    group_by(subj_id) %>% 
    filter(is_congruent_trial == "Congruent") %>% 
    summarise(accuracy = sum(correct) / length(correct))
  
  flanker_incongruent_accuracy <- flanker_data %>% 
    group_by(subj_id) %>% 
    filter(is_congruent_trial == "Incongruent") %>% 
    summarise(accuracy = sum(correct) / length(correct))
  
  flanker_congruency_accuracy <- flanker_congruent_accuracy
  flanker_congruency_accuracy$accuracy <- 
    flanker_incongruent_accuracy$accuracy - flanker_congruent_accuracy$accuracy
  
  
  # rename the id column
  # flanker_congruency_accuracy <- flanker_congruency_accuracy %>% 
  #  rename(id = subj_id)
  
  return(flanker_congruency_accuracy)
  
}

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### calculate the rt congruency effect

get_flanker_congruency_rt <- function(flanker_data){
  
  flanker_congruent_rt <- flanker_data %>% 
    filter(accuracy == 1) %>% 
    group_by(id) %>% 
    filter(congruency == "congruent") %>% 
    summarise(rt = mean(rt))
  
  flanker_incongruent_rt <- flanker_data %>% 
    filter(accuracy == 1) %>% 
    group_by(id) %>% 
    filter(congruency == "incongruent") %>% 
    summarise(rt = mean(rt))
  
  flanker_congruency_rt <- flanker_congruent_rt
  flanker_congruency_rt$rt <- 
    flanker_incongruent_rt$rt - flanker_congruent_rt$rt
  
  flanker_congruency_rt$rt <- round(flanker_congruency_rt$rt, 0)
  
  return(flanker_congruency_rt)
  
}
#------------------------------------------------------------------------------



#-----------------------------------------------------------------------------
### function to return specific data on participant
check_participant <- function(data, id, test){
  
  if(test == "flanker"){
    temp <- data$flanker %>% 
      filter(Participant.Public.ID == id)
    return(temp)
  }
  
  if(test == "shps"){
    temp <- data$shps %>% 
      filter(Participant.Public.ID == id)
    return(temp)
  }
  
  
  if(test == "qids")
    temp <- data$qids %>% 
      filter(Participant.Public.ID == id)
  return(temp)
}
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
view_shps <- function(data){
  
  # get the shps data
  shps <- data$shps
  
  # filter out the "end of file" rows in the data set
  shps <- shps %>% filter(Experiment.ID != "NA")
  
  # get the list of unique participants
  participants <- na.exclude(unique(shps$Participant.Public.ID))
  
  # create a data frame to store data
  shps_scores <- data.frame(matrix(0, nrow = length(participants), ncol = 15))


  
    # loop over each subject 
  for(i in 1:length(participants)){
    
    # get their data
    temp <- shps %>% filter(Participant.Public.ID == participants[i])
    
    # select the relevant rows
    temp <- temp %>% filter(stringr::str_detect(Question.Key, "quantised"))
    
    # change response to numeric
    temp$Response <- as.numeric(as.character(temp$Response))
    
    shps_scores[i, 1] <- as.character(temp$Participant.Public.ID)[1]
    shps_scores[i, 2] <- temp %>% 
      filter(Question.Key == "shps_1-quantised") %>% 
      pull(Response)
    shps_scores[i, 3] <- temp %>% 
      filter(Question.Key == "shps_2-quantised") %>% 
      pull(Response)
    shps_scores[i, 4] <- temp %>% 
      filter(Question.Key == "shps_3-quantised") %>% 
      pull(Response)
    shps_scores[i, 5] <- temp %>% 
      filter(Question.Key == "shps_4-quantised") %>% 
      pull(Response)
    shps_scores[i, 6] <- temp %>% 
      filter(Question.Key == "shps_5-quantised") %>% 
      pull(Response)
    shps_scores[i, 7] <- temp %>% 
      filter(Question.Key == "shps_6-quantised") %>% 
      pull(Response)
    shps_scores[i, 8] <- temp %>% 
      filter(Question.Key == "shps_7-quantised") %>% 
      pull(Response)
    shps_scores[i, 9] <- temp %>% 
      filter(Question.Key == "shps_8-quantised") %>% 
      pull(Response)
    shps_scores[i, 10] <- temp %>% 
      filter(Question.Key == "shps_9-quantised") %>% 
      pull(Response)
    shps_scores[i, 11] <- temp %>% 
      filter(Question.Key == "shps_10-quantised") %>% 
      pull(Response)
    shps_scores[i, 12] <- temp %>% 
      filter(Question.Key == "shps_11-quantised") %>% 
      pull(Response)
    shps_scores[i, 13] <- temp %>% 
      filter(Question.Key == "shps_12-quantised") %>% 
      pull(Response)
    shps_scores[i, 14] <- temp %>% 
      filter(Question.Key == "shps_13-quantised") %>% 
      pull(Response)
    shps_scores[i, 15] <- temp %>% 
      filter(Question.Key == "shps_14-quantised") %>% 
      pull(Response)
    
  } # end of participant loop
  
  
  return(shps_scores)
}
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
view_qids <- function(data){
  
  # get the qids data
  qids <- data$qids
  
  # filter out the "end of file" rows in the data set
  qids <- qids %>% filter(Experiment.ID != "NA")
  
  # get the list of unique participants
  participants <- na.exclude(unique(qids$Participant.Public.ID))
  
  # create a data frame to store data
  qids_scores <- data.frame(matrix(0, nrow = length(participants), ncol = 17))
  
  
  
  # loop over each subject 
  for(i in 1:length(participants)){
    
    # get their data
    temp <- qids %>% filter(Participant.Public.ID == participants[i])
    
    # select the relevant rows
    temp <- temp %>% filter(stringr::str_detect(Question.Key, "quantised"))
    
    # change response to numeric
    temp$Response <- as.numeric(as.character(temp$Response))
    
    qids_scores[i, 1] <- as.character(temp$Participant.Public.ID)[1]
    qids_scores[i, 2] <- temp %>% 
      filter(Question.Key == "qids_1-quantised") %>% 
      pull(Response)
    qids_scores[i, 3] <- temp %>% 
      filter(Question.Key == "qids_2-quantised") %>% 
      pull(Response)
    qids_scores[i, 4] <- temp %>% 
      filter(Question.Key == "qids_3-quantised") %>% 
      pull(Response)
    qids_scores[i, 5] <- temp %>% 
      filter(Question.Key == "qids_4-quantised") %>% 
      pull(Response)
    qids_scores[i, 6] <- temp %>% 
      filter(Question.Key == "qids_5-quantised") %>% 
      pull(Response)
    qids_scores[i, 7] <- temp %>% 
      filter(Question.Key == "qids_6-quantised") %>% 
      pull(Response)
    qids_scores[i, 8] <- temp %>% 
      filter(Question.Key == "qids_7-quantised") %>% 
      pull(Response)
    qids_scores[i, 9] <- temp %>% 
      filter(Question.Key == "qids_8-quantised") %>% 
      pull(Response)
    qids_scores[i, 10] <- temp %>% 
      filter(Question.Key == "qids_9-quantised") %>% 
      pull(Response)
    qids_scores[i, 11] <- temp %>% 
      filter(Question.Key == "qids_10-quantised") %>% 
      pull(Response)
    qids_scores[i, 12] <- temp %>% 
      filter(Question.Key == "qids_11-quantised") %>% 
      pull(Response)
    qids_scores[i, 13] <- temp %>% 
      filter(Question.Key == "qids_12-quantised") %>% 
      pull(Response)
    qids_scores[i, 14] <- temp %>% 
      filter(Question.Key == "qids_13-quantised") %>% 
      pull(Response)
    qids_scores[i, 15] <- temp %>% 
      filter(Question.Key == "qids_14-quantised") %>% 
      pull(Response)
    qids_scores[i, 16] <- temp %>% 
      filter(Question.Key == "qids_15-quantised") %>% 
      pull(Response)
    qids_scores[i, 17] <- temp %>% 
      filter(Question.Key == "qids_16-quantised") %>% 
      pull(Response)
    
  } # end of participant loop
  
  
  return(qids_scores)
}
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
### opposite of %in% 

'%!in%' <- function(x,y)!('%in%'(x,y))
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
### assess model fit via qq plots
qq_plot <- function(human_data, dstp_model_data, ssp_model_data, 
                    n_trials = 50000){
  

  # get the congreunt quantiles from human data
  congruent_human <- flanker_data %>% 
    filter(congruency == "congruent") %>% 
    group_by(id) %>% 
    summarise(acc = mean(accuracy), 
              quant_1 = quantile(rt, probs = 0.25), 
              quant_2 = quantile(rt, probs = 0.5), 
              quant_3 = quantile(rt, probs = 0.75))
  
  
  # get the incongreunt quantiles from human data
  incongruent_human <- flanker_data %>% 
    filter(congruency == "incongruent") %>% 
    group_by(id) %>% 
    summarise(acc = mean(accuracy), 
              quant_1 = quantile(rt, probs = 0.25), 
              quant_2 = quantile(rt, probs = 0.5), 
              quant_3 = quantile(rt, probs = 0.75))
  
  
 ############################
 # do the DSTP model
 ############################
  
  # get the congruent quantiles from model data
  congruent_model_dstp <- congruent_human %>% 
    group_by(id) %>% 
    mutate(acc = 0, 
           quant_1 = 0, 
           quant_2 = 0, 
           quant_3 = 0)
  
  incongruent_model_dstp <- incongruent_human %>% 
    group_by(id) %>% 
    mutate(acc = 0, 
           quant_1 = 0, 
           quant_2 = 0, 
           quant_3 = 0)
  
  
  for(i in 1:nrow(congruent_model_dstp)){
    
    print(paste("dstp_", i, sep = ""))
    
    # get the current id 
    curr_id <- as.character(congruent_model_dstp$id[i])
    
    # get that person's best parameters
    id_parms <- dstp_model_data %>% 
      filter(as.character(id) == as.character(curr_id))
    
    # simulate their data
    id_sim <- simulateDSTP(parms = as.numeric(id_parms[2:8]), 
                           nTrials = n_trials)
    id_sim <- id_sim %>% mutate(id = curr_id)
    
    # store the congruent result
    congruent_model_dstp[i, 1:5] <- id_sim %>% 
      filter(congruency == "congruent") %>% 
      summarise(id = unique(id), 
                acc = mean(accuracy), 
                quant_1 = quantile(rt, probs = 0.25), 
                quant_2 = quantile(rt, probs = 0.5), 
                quant_3 = quantile(rt, probs = 0.75))
    
    # store the incongruent result
    incongruent_model_dstp[i, 1:5] <- id_sim %>% 
      filter(congruency == "incongruent") %>% 
      summarise(id = unique(id), 
                acc = mean(accuracy), 
                quant_1 = quantile(rt, probs = 0.25), 
                quant_2 = quantile(rt, probs = 0.5), 
                quant_3 = quantile(rt, probs = 0.75))
    
  }
  
  
  ############################
  # do the SSP model
  ############################
  
  # get the congruent quantiles from model data
  congruent_model_ssp <- congruent_human %>% 
    group_by(id) %>% 
    mutate(acc = 0, 
           quant_1 = 0, 
           quant_2 = 0, 
           quant_3 = 0)
  
  incongruent_model_ssp <- incongruent_human %>% 
    group_by(id) %>% 
    mutate(acc = 0, 
           quant_1 = 0, 
           quant_2 = 0, 
           quant_3 = 0)
  
  
  for(i in 1:nrow(congruent_model_ssp)){
    
    print(paste("ssp", i, sep = ""))
    
    # get the current id 
    curr_id <- as.character(congruent_model_ssp$id[i])
    
    # get that person's best parameters
    id_parms <- ssp_model_data %>% 
      filter(as.character(id) == as.character(curr_id))
    
    # simulate their data
    id_sim <- simulateSSP(parms = as.numeric(id_parms[2:6]), 
                          nTrials = n_trials)
    id_sim <- id_sim %>% mutate(id = curr_id)
    
    # store the congruent result
    congruent_model_ssp[i, 1:5] <- id_sim %>% 
      filter(congruency == "congruent") %>% 
      summarise(id = unique(id), 
                acc = mean(accuracy), 
                quant_1 = quantile(rt, probs = 0.25), 
                quant_2 = quantile(rt, probs = 0.5), 
                quant_3 = quantile(rt, probs = 0.75))
    
    # store the incongruent result
    incongruent_model_ssp[i, 1:5] <- id_sim %>% 
      filter(congruency == "incongruent") %>% 
      summarise(id = unique(id), 
                acc = mean(accuracy), 
                quant_1 = quantile(rt, probs = 0.25), 
                quant_2 = quantile(rt, probs = 0.5), 
                quant_3 = quantile(rt, probs = 0.75))
    
  }
  
  # translate the model data to milliseconds
  congruent_model_dstp[, 3:5] <- congruent_model_dstp[, 3:5] * 1000
  incongruent_model_dstp[, 3:5] <- incongruent_model_dstp[, 3:5] * 1000
  congruent_model_ssp[, 3:5] <- congruent_model_ssp[, 3:5] * 1000
  incongruent_model_ssp[, 3:5] <- incongruent_model_ssp[, 3:5] * 1000
  
  # do the qq_plot
  setwd(here("paper_figures"))
  
  # set what COLOUR the points should be 
  plot_col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3)
  
  png("qq_model_fit.png", width = 8.5, height = 8.5, units = "in", res = 1200)
  par(mfrow = c(4,4))
  
  # congruent_dstp
  plot(congruent_human$acc, as.numeric(congruent_model_dstp$acc), pch = 19, 
       col = plot_col, xlim = c(0.5, 1), ylim = c(0.5, 1), 
       xlab = "Observed", ylab = "Predicted", main = "Accuracy")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  text(0.55, 0.95, "A", cex = 2)
  
  plot(congruent_human$quant_1, as.numeric(congruent_model_dstp$quant_1), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "25th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  plot(congruent_human$quant_2, as.numeric(congruent_model_dstp$quant_2), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "50th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  plot(congruent_human$quant_3, as.numeric(congruent_model_dstp$quant_3), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "75th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  
  # incongruent_dstp
  plot(incongruent_human$acc, as.numeric(incongruent_model_dstp$acc), pch = 19, 
       col = plot_col, xlim = c(0.5, 1), ylim = c(0.5, 1), 
       xlab = "Observed", ylab = "Predicted", main = "Accuracy")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  text(0.55, 0.95, "B", cex = 2)
  
  plot(incongruent_human$quant_1, as.numeric(incongruent_model_dstp$quant_1), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "25th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  plot(incongruent_human$quant_2, as.numeric(incongruent_model_dstp$quant_2), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "50th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  plot(incongruent_human$quant_3, as.numeric(incongruent_model_dstp$quant_3), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "75th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  
  # congruent_ssp
  plot(congruent_human$acc, as.numeric(congruent_model_ssp$acc), pch = 19, 
       col = plot_col, xlim = c(0.5, 1), ylim = c(0.5, 1), 
       xlab = "Observed", ylab = "Predicted", main = "Accuracy")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  text(0.55, 0.95, "C", cex = 2)
  
  plot(congruent_human$quant_1, as.numeric(congruent_model_ssp$quant_1), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "25th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  plot(congruent_human$quant_2, as.numeric(congruent_model_ssp$quant_2), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "50th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  plot(congruent_human$quant_3, as.numeric(congruent_model_ssp$quant_3), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "75th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  
  # incongruent_ssp
  plot(incongruent_human$acc, as.numeric(incongruent_model_ssp$acc), pch = 19, 
       col = plot_col, xlim = c(0.5, 1), ylim = c(0.5, 1), 
       xlab = "Observed", ylab = "Predicted", main = "Accuracy")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  text(0.55, 0.95, "D", cex = 2)
  
  plot(incongruent_human$quant_1, as.numeric(incongruent_model_ssp$quant_1), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "25th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  plot(incongruent_human$quant_2, as.numeric(incongruent_model_ssp$quant_2), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "50th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  plot(incongruent_human$quant_3, as.numeric(incongruent_model_ssp$quant_3), 
       pch = 19, col = plot_col, xlim = c(0, 1200), ylim = c(0, 1200), 
       xlab = "Observed", ylab = "Predicted", main = "75th Percentile")
  abline(a = 0, b = 1, lwd = 1.5, lty = 2)
  
  dev.off()
  
}
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
do_clean_graphs <- function(){
  
  op <- par(cex.main = 1.5, mar = c(5, 5, 4, 1) + 0.1, mgp = c(3.5, 1, 0), 
            cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
  
}
#-----------------------------------------------------------------------------
