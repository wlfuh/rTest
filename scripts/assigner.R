#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("pracma"))

option_list <- list( 
  make_option(c("-i", "--iterations"), type="integer", default=500,
              help="number of probablistic assignments [default %default]"),
  make_option(c("-s", "--scale"), type="numeric", default=0.75,
              help="scale value for noise added to the data [default %default]"),
  make_option(c("-o", "--output"), type="character",default="assigned_shifts.txt",
              help="output file name [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="print header and progress information [default %default]")
)
parser <- OptionParser(usage = "%prog [options] predicted_chemical_shifts chemical_shift_peaks accuracy_file", option_list=option_list)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if(length(arguments$args) != 3) {
  cat("Incorrect number of required positional arguments\n\n")
  print_help(parser)
  stop()
} else {
  if (opt$verbose){
    cat("Project: CSAssigner -- Optimial Chemical Shift Assignment\n")
    cat("Author: William Fuh\n")
    cat("Author: Aaron T. Frank\n")
    cat(sprintf("%s\n",date()))
  }
  
  #user functions
  # Goto working directoru
  setwd("~/GitSoftware/rTest/")
  
  # load library and our custom functions
  library(pracma)
  library(plyr)
  source("scripts/library.R")
  
  # get arguments
  predicted_shifts_file <- arguments$args[1]
  refdata_file <- arguments$args[2]
  accu_file <- arguments$args[3]

  # get options
  iter <- opt$iterations
  scale <- opt$scale
  output <- opt$output
  verbose <- opt$verbose
  
  # read in reference data
  refdata <- read.table(refdata_file, header=FALSE)
  colnames(refdata) <- c("resname","resid","nucleus","cs","dummy")
  
  # larmord accuracy data
  accu <- read.table(accu_file,header=FALSE)
  colnames(accu) <- c("resname", "nucleus", "error")
  
  # assign chemical shifts
  assigned <- matchModel(filename = predicted_shifts_file, iter = iter, scale = scale)
  write.table(assigned, file = output, quote = FALSE, row.names = FALSE, col.names = FALSE)
}


