#' In the files suprise_cntr.csv and suprise_exp.csv there is a problem with the coding
#' of the block. Sometimes the same index is used for multile blocks.  This coding error
#' has been corrected in the file surprise_data.csv.



options(max.print=100000000)


for (i in 1:167) {
  surprise_df <- correct_block(surprise_df, i)
  print(i)
}


n_df <- surprise_df %>% 
  group_by(id, block) %>% 
  summarise(
    n = n()
  )

data.frame(n_df)

# add number of trials for each block

surprise_df <- surprise_df %>% 
  group_by(id, block) %>% 
  mutate(
    n_trials_by_block = n()
  )


good_blks <- surprise_df[surprise_df$n_trials_by_block > 10, ]





temp <- surprise_df[surprise_df$id == 48, ]
temp1 <- temp[!temp$block %in% c(1, 2, 9), ]
temp1$block <- temp1$block - 2

without_one <- surprise_df[surprise_df$id != 48, ]

surprise_df <- rbind(surprise_df, temp1)



write_csv(surprise_df, "surprise_data.csv")



dd <- surprise_df[surprise_df$id == 19, ]

# temporal difference between the current and the previous trial
dd$tdiff <- unlist(tapply(dd$time_of_day, INDEX = dd$block,
                            FUN = function(x) c(0, diff(as.numeric(x)))))


# difference with respect to the trial index
dd$trial_index_diff <- c(-79, diff(dd$trial_within_blk))

data.frame(
  trial_within_blk=dd$trial_within_blk, 
  trial_index_diff=dd$trial_index_diff, 
  tdiff=dd$tdiff
)

# each time trial_index_diff is negative, block is incremented by 1

dd$new_block <- rep(0, nrow(dd))
increment <- 0
for (i in 1:nrow(dd)) {
  
  if (dd$trial_index_diff[i] < 0)
    increment = increment + 1
  
    dd$new_block[i] = dd$new_block[i] + increment
}

data.frame(dd$new_block, dd$trial_within_blk)

dd$block <- dd$new_block

tapply(dd$rtTukey, dd$block, length)


# remove one subject from data frame

temp <- surprise_df[surprise_df$id != 17, ]

surprise_df <- rbind(temp, dd)
