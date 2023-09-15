

# Get the list of files
if (exp_flag == "SURPRISE") {
  filenames <- Sys.glob("../data/surprise_videos/*.dat") # both Lucrezia's and Elena's data
  n_files <- length(filenames)
} else if (exp_flag == "CONTROL") { # control experiment
  filenames <- Sys.glob("../data/no_surprise_videos/*.dat")
  n_files <- length(filenames)
} else {
  message("\nWhich experiment do you want to analyze?\n")
}


mylist <- list()
for (i in 1:n_files) {
  
  tempfile <- read.table(filenames[i], header = TRUE, row.names=NULL)
  
  if(dim(tempfile)[2] == 20)
    temp_data <- tempfile[, 2:20]
  
  mylist[[i]] <- 
    
}







file_list <- 
  list.files(path="/Users/corrado/Dropbox/papers/surprise/data/surprise_videos")

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
dataset <- data.frame()

#had to specify columns to get rid of the total column
for (i in 1:length(file_list)){
  temp_data <- read_excel(file_list[i], range = cell_cols("A:H")) #each file will be read in, specify which columns you need read in to avoid any errors
  
  temp_data <- read.table(filenames[i], header = TRUE, row.names=NULL)
  
  
  dataset <- rbind(dataset, temp_data) #for each iteration, bind the new data to the building dataset
}

