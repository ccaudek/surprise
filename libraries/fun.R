# ------------------------------------------------------------------------------
# Custom function for the project
# 'DDM for post-error RTs in trials with different penalities'
# Written by Corrado Caudek
# This script was last modified on "Sun Oct 22 08:46:23 2017"
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' @title Select participants and trials for running the DDM in different conditions
#'
#' @param input a data.frame and the specification of the condition
#' @return output a data.frame
#'
prep_data <- function(df, 
                     EXP,
                     FAM_COMB,
                     IS_CONGR) {
  
  subdf <- df %>%
    dplyr::filter(experiment        == EXP & 
                  fam_combinations  == FAM_COMB &
                  is_congruent      == IS_CONGR)
  
  # Remove rows with NA on the rt variable
  thedat <- subdf[!is.na(subdf$rt), ]
  
  # Put data in proper format for choiceRT_ddm()
  thedat$subjID <- as.numeric(factor(thedat$sID))
  thedat$choice <- ifelse(thedat$correct == 1, 2, 1)
  thedat$RT     <- thedat$rt / 1000
  
  thedat <- as.data.frame(thedat)

  # choiceRT_ddm() requires subjID to be sorted   
  thedat <- thedat %>% 
    arrange(subjID) %>%
    dplyr::select(subjID, choice, RT) %>%
    ungroup()
  
  thedat
}


# ------------------------------------------------------------------------------
#' @title insure that each subject has at least one error
#'
#' @param input a data.frame and how to proceed
#' @return output a data.frame
#'
manage_errors_by_subj <- function(thedat, 
                                  PROCEDURE) {
  
  # find subject with all correct resonses (choice == 2)
  bysub_p <- thedat %>%
    group_by(subjID) %>%
    summarise(
      avg_p = mean(choice)
    ) %>%
    as.data.frame
  
  index <- bysub_p$avg_p == 2.00
  # list of subjID with perfect performance (no errors)
  bad_id <- bysub_p$subjID[index] 
  
  if (PROCEDURE == 1) { # code as error the first reponse of each subjects 
                        # with perfect performance
    # select only subjects with 0 errors
    no_err_df <- thedat[thedat$subjID %in% bad_id, ]
    # find number of trials for those subjects
    no_err_n <- no_err_df %>%
      group_by(subjID) %>%
      summarise(
        n = n()
      ) 
    # add vector ii with trials from 1 to n for each subject
    info <- data.frame(start = 1, len = no_err_n$n)
    no_err_df$ii <- sequence(info$len) + rep(info$start-1, info$len)
    # code as error the first trial of each subject
    index <- no_err_df$ii == 1
    no_err_df$choice[index] <- 1 # 1 = error; 2 = correct 
    # remove ii column 
    no_err_df$ii <- NULL
    # select subjects with at least one error
    with_err_df <- thedat[!thedat$subjID %in% bad_id, ]
    # combine the two dataframes
    ddm_df <- rbind(no_err_df, with_err_df)
  } else { # remove subjects with perfect performance
    ddm_df <- thedat[!thedat$subjID %in% bad_id, ]
  }
  
  ddm_df
}

