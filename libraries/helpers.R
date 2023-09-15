# The palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

correct_block <- function(df, ID) {
  
  dd <- df[df$id == ID, ]
  
  # difference with respect to the trial index
  dd$trial_index_diff <- c(-79, diff(dd$trial_within_blk))
  
  # each time trial_index_diff is negative, block is incremented by 1
  dd$new_block <- rep(0, nrow(dd))
  increment <- 0
  for (i in 1:nrow(dd)) {
    
    if (dd$trial_index_diff[i] < 0)
      increment = increment + 1
    
    dd$new_block[i] = dd$new_block[i] + increment
  }
  
  dd$block <- dd$new_block
  dd$trial_index_diff <- NULL
  dd$new_block <- NULL
  
  # remove one subject from data frame
  temp <- df[df$id != ID, ]
  # add subject with corrected block index
  new_df <- rbind(temp, dd)
  
  # return the complete data frame
  new_df
  
}



detrend_id_block <- function(df){
  
  mydata <- surprise_df
  
  # create a unique identifier for each id and block combination
  mydata[is.na(mydata)] = ''
  mydata <- unite(mydata, id_block, subject_name, block, sep = '')
  
  n_id_block <- length(unique(mydata$id_block))
  n_id_block
  
  # sub_block_index from 1 to 776
  mydata$id_block_index <- as.numeric(as.factor(as.character(mydata$id_block)))
  # sort(unique(mydata$subjIndex))
  
  out <- surprise_df %>% 
    dplyr::filter(experiment == "Control") %>% 
    group_by(id, block) %>% 
    summarise(
      n = n()
    )
  
data.frame(out)
  
  
  
temp <- surprise_df[surprise_df$id == 17 & surprise_df$block == 2, ]
tapply(temp$rtTukey, temp$block, length)
temp$time_of_day

temp$tdiff <- unlist(tapply(temp$time_of_day, INDEX = temp$block,
                          FUN = function(x) c(0, diff(as.numeric(x)))))
  
  
  mylist <- list() # create an empty list
  
  for (i in 1:n_id_block) {
    dsogg <- mydata[mydata$id_block_index == i, ]
    y <- dsogg$rtTukey 
    yd <- detrend.series(y, method = "Spline")  
    
    sigma_ratio <- sd(dsogg$rtTukey) / sd(yd)
    yy2 <- yd * sigma_ratio
    drt <- (yy2 - mean(yy2)) + mean(dsogg$rtTukey)
    
    mylist[[i]] <- drt
    rm(y, yd, dsogg)
    print(i)
  }
  
  return(mylist)
}




detrend <- function(df){

  mydata <- df
  mydata$id <- factor(mydata$id)
  nSubj <- length(unique(mydata$id))
  nSubj

  # id is not 1 ... max. subjIndex has this property
  mydata$subjIndex <- as.numeric(as.factor(mydata$id))
  # unique(mydata$subjIndex)

  mylist <- list() # create an empty list

  for (i in 1:nSubj) {
    dsogg <- mydata[mydata$subjIndex == i, ]
    y <- dsogg$rtTukey
    yd <- detrend.series(y, method = "Spline")

    sigma_ratio <- sd(dsogg$rtTukey) / sd(yd)
    yy2 <- yd * sigma_ratio
    drt <- (yy2 - mean(yy2)) + mean(dsogg$rtTukey)

    mylist[[i]] <- drt
    rm(y, yd, dsogg)
    print(i)
  }

  return(mylist)
}





to_snake_case <- function(camelcases){
  # catch some input that should be handled like underscores too
  camelcases <- stringr::str_replace_all(camelcases, "\\s+|\\.+", "_")
  # get to know, if a string starts with a small letter
  small_start <- !is.na(stringr::str_extract(camelcases, "^[a-z]"))
  # get all capital letter sequences from a string
  capitals <- stringr::str_extract_all(camelcases, "[A-Z]+")
  # Setting an underscore before capital and first letters
  starts <- purrr::pmap(list(camelcases,
                             small_start,
                             capitals),
                        function(x,y,z)
                          if (length(z) == 0) {"_"} else {
                            c("_", paste0("_", z))
                          }
  )
  # split the strings by their capital letter sequences.
  rests <- stringr::str_split(camelcases, "[A-Z]+")
  # setting all peaces together:
  # - pasting first and capital letters with the rest of the string
  # - applying tolower, remove more than one "_" and starting "_"
  corrected <- purrr::map2_chr(starts, rests, stringr::str_c, collapse = "") %>% 
    purrr::map_chr(stringr::str_to_lower) %>% 
    purrr::map_chr(~ stringr::str_replace_all(.x, "_+", "_")) %>% 
    purrr::map_chr(~ stringr::str_replace_all(.x, "^_+|_+$", ""))
  corrected
}




plot_rt_bysub <- function(d, ID) {
  
  temp <- d %>%
    dplyr::filter(
      subj_id == ID & !is.na(rtTukey)
    )
  
  publPalette <- c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
  
  p1 <- (
    ggplot(data=temp, aes(temp$rtTukey)) +
    geom_histogram(aes(y =..density..), col="red", fill="green", alpha = .2) +
    stat_function(fun=dnorm, args=list(mean=mean(temp$rtTukey), sd=sd(temp$rtTukey))) +
    labs(x="Reaction Times (ms)", y="Density") +
    # ggtitle("Raw RT") +
      scale_fill_manual(values = cbp1) + 
      scale_colour_manual(values=cbp1)
    # scale_colour_manual(values=publPalette) +
    # scale_fill_manual(values=publPalette)
  )
  
  print(p1)
  
}


# Create plot 1 -----------------------------------------------------------

create_plot_1 <- function(data) {
  
  bysub_rt <- data %>% 
    group_by(experiment, subject_name, is_congruent_trial) %>% 
    summarise(
      mlrt = mean(rtTukey, na.rm = TRUE)
    )
  
  congr_df <- bysub_rt %>% 
    dplyr::filter(is_congruent_trial == "Congruent")
  congr_df$is_congruent_trial <- NULL
  congr_df <- rename(congr_df, mlrt_con = mlrt)
  
  incon_df <- bysub_rt %>% 
    dplyr::filter(is_congruent_trial == "Incongruent")
  incon_df$is_congruent_trial <- NULL
  incon_df <- rename(incon_df, mlrt_inc = mlrt)
  
  total <- merge(
    congr_df, incon_df, by=c("experiment", "subject_name")
  )
  
  total$congr_eff <- total$mlrt_inc - total$mlrt_con
  
  p <- total %>%
    ggplot(aes(x=experiment, y=congr_eff)) +
    geom_boxplot(fill="skyblue", notch=TRUE) +
    geom_jitter(size=0.9, color="orange", width=0.1) +
    geom_hline(yintercept=0, linetype="dashed", color = "gray") +
    xlab("Experiment") +
    ylab("Congruency Effect (log ms)") +
    theme_apa()
  
  print(p)
  
}



# Create plot 2 -----------------------------------------------------------

create_plot_2 <- function(data) {
  
  bysub_acc <- data %>% 
    dplyr::filter(!is.na(rtTukey)) %>% 
    group_by(experiment, subject_name) %>% 
    summarise(
      avg_acc = mean(correct, na.rm = TRUE),
      n = n()
    )
  data.frame(bysub_acc)
  
  hist(bysub_acc$avg_acc)
  
  p <- ggplot(bysub_acc, aes(x=avg_acc)) +
    geom_histogram(
      data=subset(bysub_acc, experiment == 'Control'),
      fill = "red", alpha = 0.2) +
    geom_histogram(
      data=subset(bysub_acc, experiment == 'Surprise'),
      fill = "blue", alpha = 0.2) +
    xlab("Average Accuracy") +
    ylab("Frequency") +
    theme_apa()
  
  print(p)
  
}



# Create plot 3 -----------------------------------------------------------

create_plot_3 <- function(data) {
  
  # Is the previous trial congruent?
  data <- data %>%
    # dplyr::filter(tot_trials > 80) %>%
    group_by(experiment, ID, block) %>%
    mutate(
      is_congruent_prev_trial = lag(is_congruent_trial, 1)
    )
  
  # check it out
  if (0) {
    data.frame(
      i = final$trial_within_blk,
      cond = final$is_congruent_trial,
      prev = final$is_congruent_prev_trial
    )[1:200, ]
  }
  
  bysub_rt <- data %>%
    # dplyr::filter(trials_after_clip == 0) %>%
    group_by(
      experiment, ID, is_congruent_trial,
      is_congruent_prev_trial
    ) %>%
    summarise(
      mlrt = mean(rtTukey, na.rm = TRUE)
    )
  
  # Remove NAs
  bysub_rt <- bysub_rt[!is.na(bysub_rt$is_congruent_prev_trial), ]
  
  congr_df <- bysub_rt %>%
    dplyr::filter(is_congruent_trial == "Congruent")
  congr_df$is_congruent_trial <- NULL
  congr_df <- dplyr::rename(congr_df, mlrt_con = mlrt)
  
  incon_df <- bysub_rt %>%
    dplyr::filter(is_congruent_trial == "Incongruent")
  incon_df$is_congruent_trial <- NULL
  incon_df <- dplyr::rename(incon_df, mlrt_inc = mlrt)
  
  total <- merge(
    congr_df, incon_df,
    by = c("experiment", "ID", "is_congruent_prev_trial")
  )
  
  total$congr_eff <- total$mlrt_inc - total$mlrt_con
  
  total$experiment <- factor(total$experiment)
  
  p <- total %>%
    ggplot(aes(x=is_congruent_prev_trial, y=congr_eff)) +
    geom_boxplot(fill="skyblue", notch=FALSE) +
    geom_jitter(size=0.9, color="orange", width=0.1) +
    facet_wrap(~ experiment) +
    geom_hline(yintercept=0, linetype="dashed", color = "gray") +
    xlab("Previous Trial") +
    ylab("Congruency Effect (log ms)") +
    papaja::theme_apa() +
    theme(axis.text.x=element_text(angle=15,hjust=1))
  
  print(p)
  
}


# Create plot 4 -----------------------------------------------------------


create_plot_4 <- function(data, WHICH_DATA) {
  
  if (WHICH_DATA == "FIRST_TRIAL")
    data <- data[data$trials_after_clip == 0, ]
  if (WHICH_DATA == "SECOND_TRIAL")
    data <- data[data$trials_after_clip == 1, ]
  if (WHICH_DATA == "THIRD_TRIAL")
    data <- data[data$trials_after_clip == 2, ]
  if (WHICH_DATA == "FOURTH_TRIAL")
    data <- data[data$trials_after_clip == 3, ]
  
  quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  exp <- c("Control", "Surprise")
  df <- list()
  
  for (i_exp in 1:2) {
    
    data_exp <- data[data$experiment == exp[i_exp], ]
    
    subs_names <- unique(data_exp$subject_name)
    n_subs <- length(subs_names)
    cdf_data <- matrix(0, nrow = length(quantiles), ncol = n_subs)
    
    condition <- c("Incongruent", "Congruent")
    avg_CDF <- list()
    
    for (i_cond in 1:2) {
      data_one_cond <-
        data_exp[data_exp$is_congruent_trial == condition[i_cond], ]
      
      for (i in 1:n_subs) {
        temp_data <- subset(
          data_one_cond,
          data_one_cond$subject_name == subs_names[i]
        )
        cdf_data[, i] <- quantile(
          temp_data$rtTukey,
          quantiles,
          na.rm = TRUE
        )
      }
      avg_CDF[[i_cond]] <- cdf_data
    }
    
    # Compute the congruency effect (inc - con) by subject
    congr_eff_mat <- avg_CDF[[1]] - avg_CDF[[2]]
    
    df_plot <- data.frame(
      quantiles,
      eff = apply(congr_eff_mat, 1, mean),
      stderr = apply(congr_eff_mat, 1, sd) / sqrt(n_subs)
    )
    
    df[[i_exp]] <- df_plot
    
  }
  
  df_plot <- rbind(df[[1]], df[[2]])
  df_plot$experiment <- rep(
    c("Control", "Surprise"),
    each = length(quantiles)
  )
  
  p <- ggplot(df_plot,
              aes(x = quantiles, y = eff, color = experiment)) +
    geom_errorbar(
      aes(ymin = eff - stderr, ymax = eff + stderr), width = .05) +
    geom_line() +
    coord_cartesian(ylim = c(-35, 65)) +
    geom_point(size = 3) +
    ylab("Congruency Effect (ms)") +
    xlab("Quantiles") +
    labs(color = "Experiment") +
    theme_apa()
  
  print(p)
  
}


# Create plot 5 -----------------------------------------------------------


create_plot_5 <- function(data, WHICH_DATA) {
  
  if (WHICH_DATA == "FIRST_TRIAL")
    data <- data[data$trials_after_clip == 0, ]
  if (WHICH_DATA == "SECOND_TRIAL")
    data <- data[data$trials_after_clip == 1, ]
  if (WHICH_DATA == "THIRD_TRIAL")
    data <- data[data$trials_after_clip == 2, ]
  if (WHICH_DATA == "FOURTH_TRIAL")
    data <- data[data$trials_after_clip == 3, ]
  
  quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  
  exp_cond <- rep(NA, nrow(data))
  for (i in 1:nrow(data)) {
    if (data$experiment[i] == "control")
      exp_cond[i] <- "Control"
    else if (data$experiment[i] == "surprise" & data$is_surprise_clip[i] == "Surprise")
      exp_cond[i] <- "Surprise"
    else if (data$experiment[i] == "surprise" & data$is_surprise_clip[i] == "No Surprise")
      exp_cond[i] <- "No Surprise"
    else
      exp_cond[i] <- "ERROR"
  }
  data$exp_cond <- exp_cond
  
  exp <- c("Control", "Surprise", "No Surprise")
  df <- list()
  avg_rt_q <- list()
  
  for (i_exp in 1:3) {
    
    data_exp <- data[data$exp_cond == exp[i_exp], ]
    
    subs_names <- unique(data_exp$ID)
    n_subs <- length(subs_names)
    cdf_data <- matrix(0, nrow = length(quantiles), ncol = n_subs)
    
    condition <- c("Incongruent", "Congruent")
    avg_CDF <- list()
    
    for (i_cond in 1:2) {
      data_one_cond <- # all the raw data for the selected condition
        data_exp[data_exp$is_congruent_trial == condition[i_cond], ]
      
      for (i in 1:n_subs) {
        temp_data <- subset(
          data_one_cond,
          data_one_cond$ID == subs_names[i]
        )
        cdf_data[, i] <- quantile(
          temp_data$rtTukey,
          quantiles,
          na.rm = TRUE
        )
      }
      avg_CDF[[i_cond]] <- cdf_data # list of two matrices with the quantiles for each subject,
      # in the incongruent and congruent conditions
    }
    
    # Compute the congruency effect (incon - con) by subject
    congr_eff_mat <- avg_CDF[[1]] - avg_CDF[[2]]
    
    df_plot <- data.frame(
      quantiles,
      eff = apply(congr_eff_mat, 1, mean),
      stderr = apply(congr_eff_mat, 1, sd) / sqrt(n_subs)
    )
    
    df[[i_exp]] <- df_plot # save congruency effect for each experiment_condition
    
    # compute the average RT for each subject and each quantile
    bysub_avg_rt_mat <- (avg_CDF[[1]] + avg_CDF[[2]]) / 2
    
    avg_rt_q[[i_exp]] <- apply(bysub_avg_rt_mat, 1, mean)
    
  }
  
  df_plot <- rbind(df[[1]], df[[2]], df[[3]])
  df_plot$exp_cond <- rep(
    c("Control", "Surprise", "No Surprise"),
    each = length(quantiles)
  )
  
  df_plot$experiment <- c(
    rep("Control Experiment", length(quantiles)),
    rep("Surprise Experiment", length(quantiles)),
    rep("Surprise Experiment", length(quantiles))
  )
  
  df_plot$avg_rt_quantiles <- c(avg_rt_q[[1]], avg_rt_q[[2]], avg_rt_q[[3]])
  
  pd <- position_dodge(0.0) # move them .05 to the left and right
  
  p <- ggplot(df_plot,
              aes(x = avg_rt_quantiles, y = eff, color = exp_cond)) +
    geom_errorbar(
      aes(ymin = eff - stderr, ymax = eff + stderr), width = .05) +
    geom_line(position=pd) +
    coord_cartesian(ylim = c(-35, 65)) +
    geom_point(size = 3, position=pd) +
    facet_wrap(~ experiment) +
    ylab("Delta RT (ms): Incon - con") +
    xlab("Mean RT (msec)") +
    labs(color = "Experiment") +
    papaja::theme_apa() +
    theme(legend.position="bottom") +
    theme(legend.title=element_blank()) +
    geom_hline(yintercept=0, linetype="dashed", color = "gray", size = 1.25) +
    # scale_color_manual(values=c('#999999','#999999', '#E69F00')) 
    scale_colour_manual(values = cbp1)
    # scale_color_grey(start = 0.8, end = 0.2) 
  
  
  # p <- p + annotate("text", 
  #                   x = 600, 
  #                   y = 39, 
  #                   label = c("", "Surprise videos"))
  # 
  # p <- p + annotate("text", 
  #                   x = 600, 
  #                   y = -8, 
  #                   label = c("", "No Surprise\nvideos"))
  
  
  print(p)
  
}

