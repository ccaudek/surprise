#' Get the data for the control group run by Giulia.
#' The function returns a data.frame. This data.frame is used in the script
#' 01_import_data.R.
#' 
get_data_giulia <- function() {
  
  # Get file names and number of files of the control experiment (Giulia)
  filenames <- Sys.glob(here("data", "no_surprise_video_giulia", "*.dat"))
  n_files <- length(filenames)
  
  # Get list of data.frames.
  mylist <- list()
  for (i in 1:n_files) {
    mylist[[i]] <- read.table(filenames[i], header = TRUE, row.names=NULL)
  }
  
  new_list <- list()
  for (i in 1:length(mylist)) {
    
    foo <- mylist[[i]]
    
    subj_name_part_1 <- unique(foo$row.names)
    subj_name_part_2 <- unique(foo$subject_name)
    subj_name_string <- paste0(subj_name_part_1, subj_name_part_2)
    
    n <- length(foo$subject_name)
    
    df <- data.frame(
      subj_name = rep(subj_name_string, n),
      subject_number = 1000, # TODO
      trial = foo$iTrial,
      target_or = foo$flankersOrientation,
      flanker_or = foo$response,
      resp = foo$correct,
      correct = foo$onscreen_stim_duration,
      movie_id = foo$rt,
      rt = foo$date,
      date = foo$time_of_day,
      time_of_day = foo$deltaX,
      is_surprise_clip = foo$isClipTrial,
      is_clip_trial = foo$isSurprise,
      movie_id2 = substring(foo$rt, 57),
      subj_id = 1000, # TODO
      block = foo$block
    )
    
    new_list[[i]] <- df
    rm(df)
  }
  
  tot_df <- do.call(what = rbind, args = new_list)
  
  tot_df$subject_number <- 5000 + as.numeric(as.factor(as.character(tot_df$subj_name)))
  tot_df$subj_id <- 5000 + as.numeric(as.factor(as.character(tot_df$subj_name)))
  
  tot_df
}


add_control_subjects <- function(tot_df) {
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 1) 
  bits1 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 2) 
  bits2 <- temp[321:480, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 3) 
  bits3 <- temp[321:400, ]
  
  s68 <- bind_rows(list(bits1, bits2, bits3))
  
  s68$subj_id <- 68
  s68$subj_name <- "subj_contr_68"
  s68$subject_number <- s68$subject_number + 1000
  dim(s68)
  
  #--------------------------------------------
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 4) 
  bits4 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 5) 
  bits5 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 7) 
  bits7 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 8) 
  bits8 <- temp[321:400, ]
  
  s69 <- bind_rows(list(bits4, bits5, bits7, bits8))
  
  s69$subj_id <- 69
  s69$subj_name <- "subj_contr_69"
  s69$subject_number <- s69$subject_number + 1000
  dim(s69)
  
  #--------------------------------------------
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 9) 
  bits9 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 10) 
  bits10 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 11) 
  bits11 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 12) 
  bits12 <- temp[321:400, ]
  
  s70 <- bind_rows(list(bits9, bits10, bits11, bits12))
  
  s70$subj_id <- 70
  s70$subj_name <- "subj_contr_70"
  s70$subject_number <- s70$subject_number + 1000
  dim(s70)
  
  #--------------------------------------------
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 13) 
  bits13 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 14) 
  bits14 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 15) 
  bits15 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 16) 
  bits16 <- temp[321:400, ]
  
  s71 <- bind_rows(list(bits13, bits14, bits15, bits16))
  
  s71$subj_id <- 71
  s71$subj_name <- "subj_contr_71"
  s71$subject_number <- s71$subject_number + 1000
  dim(s71)
  
  
  #--------------------------------------------
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 17) 
  bits17 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 27) 
  bits27 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 28) 
  bits28 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 29) 
  bits29 <- temp[321:400, ]
  
  s72 <- bind_rows(list(bits17, bits27, bits28, bits29))
  
  s72$subj_id <- 72
  s72$subj_name <- "subj_contr_72"
  s72$subject_number <- s72$subject_number + 1000
  dim(s72)
  
  #--------------------------------------------
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 30) 
  bits30 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 31) 
  bits31 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 32) 
  bits32 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 33) 
  bits33 <- temp[321:400, ]
  
  s73 <- bind_rows(list(bits30, bits31, bits32, bits33))
  
  s73$subj_id <- 73
  s73$subj_name <- "subj_contr_73"
  s73$subject_number <- s73$subject_number + 1000
  dim(s73)
  
  
  #--------------------------------------------
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 35) 
  bits35 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 36) 
  bits36 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 37) 
  bits37 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 38) 
  bits38 <- temp[321:400, ]
  
  s74 <- bind_rows(list(bits35, bits36, bits37, bits38))
  
  s74$subj_id <- 74
  s74$subj_name <- "subj_contr_74"
  s74$subject_number <- s74$subject_number + 1000
  dim(s74)
  
  #--------------------------------------------
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 43) 
  bits43 <- temp[321:400, ]
  
  temp <- tot_df %>% 
    dplyr::filter(subj_id == 45) 
  bits45 <- temp[321:400, ]
  
  s75 <- bind_rows(list(bits43, bits45))
  
  s75$subj_id <- 75
  s75$subj_name <- "subj_contr_75"
  s75$subject_number <- rep(1:2, each = 80) + 1100
  dim(s75)
  
  
  temp <- bind_rows(list(s68, s69, s70, s71, s72, s73, s74, s75))
  last <- temp[complete.cases(temp), ]
  
  last
}









bayes_rsq<- function(fit) {
  
  cond_r2 <- bayes_R2(fit, re_formula = NULL, summary = TRUE)
  
  print(
    paste(
      "cond_r2 = ", round(cond_r2[1], 4), 
      ", 95% CI: [", round(cond_r2[3], 4),
      ", ", round(cond_r2[4], 4),
      "]"
    )
  )
  
  cond_minus_slopes_r2 <- 
    bayes_R2(
      fit, 
      re_formula = ~ (1 | subj_name), 
      summary = TRUE)
  
  print(
    paste(
      "cond_minus_slopes_r2 = ", round(cond_minus_slopes_r2[1], 4), 
      ", 95% CI: [", round(cond_minus_slopes_r2[3], 4),
      ", ", round(cond_minus_slopes_r2[4], 4),
      "]"
    )
  )
  
  marginal_r2 <- bayes_R2(fit, re_formula = NA, summary = TRUE)
  
  print(
    paste(
      "marginal_r2 = ", round(marginal_r2[1], 4), 
      ", 95% CI: [", round(marginal_r2[3], 4),
      ", ", round(marginal_r2[4], 4),
      "]"
    )
  )
  
}





add_block <- function(df) {
  
  df$dtime <- paste(df$date, df$time_of_day, sep = ' ')
  df$dtime <- stringr::str_replace_all(df$dtime, "March", "03")
  df$dtime <- stringr::str_replace_all(df$dtime, "April", "04")
  df$dtime <- stringr::str_replace_all(df$dtime, "May",   "05")
  df$dtime <- lubridate::dmy_hms(df$dtime)
  
  df1 <- dplyr::arrange(df, dtime)
  
  # add block
  df1$d <- c(-79, diff(df1$trial))
  
  df1$block <- rep(NA, nrow(df1))
  index_block <- 0
  for (i in 1:nrow(df1)) {
    if (df1$d[i] != 1) {
      index_block <- index_block + 1
    }
    df1$block[i] <- index_block
  }
  df1$d <- NULL
  
  df1
}



remove_eff_block_trial_one_subj <- function(d) {
  
  d$block <- factor(d$block)
  d$trials_after_clip <- factor(d$trials_after_clip)
  
  fm <- lm(
    RT ~ block * trials_after_clip,
    data = d
  )
  
  d$drt <- fm$res + mean(d$RT)
  
  return(d)
}



remove_eff_block_trial_all_subjs <- function(df){
  
  mydata <- ds_clean
  mydata$id <- factor(mydata$id)

  # id is not 1 ... max. subjIndex has this property
  mydata$subj_index <- as.numeric(as.factor(mydata$id))
  n_subj <- length(unique(mydata$subj_index))
  
  mylist <- list() # create an empty list
  
  for (i in 1:n_subj) {
    dsogg <- mydata[mydata$subj_index == i, ]
    dd <- remove_eff_block_trial_one_subj(dsogg)
    mylist[[i]] <- dd
    rm(dd, dsogg)
    print(i)
  }
  
  # for (i in 1:n_subj) {
  #   print(i)
  #   print
  #   print(dim(mylist[[i]]))
  # }

  df <- do.call(rbind.data.frame, mylist)
  
  return(df)
}



