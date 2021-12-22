
# calc_rsem_fragment_length_stats.R
# Calculates fragment length mean and standard deviation 
# from RSEM model file.

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


rm(list=ls())

args <- commandArgs()
script_dir <- dirname(sub("--file=", "", args[4]))
print(script_dir)

print(args)
system(paste(c("git", "-C", script_dir, "rev-parse", "HEAD"), collapse = " "))
system(paste(c("git", "-C", script_dir, "rev-parse", "--abbrev-ref", "HEAD"), collapse = " "))

model_data <- read.csv(args[6], sep = " ", header = F)

bounds <-  model_data[3,seq(1,3)]
values <- sample(seq(bounds$V1 + 1, bounds$V2), size = 1000000, replace = T, prob =t(model_data[4,]))

print("")
print(paste("Mean:", mean(values)))
print(paste("Sd:", sd(values)))

print("Done")
