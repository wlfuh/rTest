#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GA"))
suppressPackageStartupMessages(library("bio3d"))
option_list <- list( 
  make_option(c("-a", "--atomSelection"), type="character",default="O1P",
              help="Atom selection for calculating inter-residue distances [default %default]"),
  make_option(c("-d", "--direct"), action="store_true", default=FALSE,
              help="whether to use direct, brute-force method to determine optimal charge. Default is GA based optimization [default %default]"),
  make_option(c("-i", "--iterations"), type="integer", default=500,
              help="number of iterations for GA algorithm [default %default]"),
  make_option(c("-p", "--population"), type="integer", default=20,
              help="population size for GA algorithm [default %default]"),
  make_option(c("-o", "--output"), type="character",default="charged.pdb",
              help="name of modified PDB file [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="print header and progress information [default %default]")
)
parser <- OptionParser(usage = "%prog [options] pdb_file charge_state", option_list=option_list)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if(length(arguments$args) != 2) {
  cat("Incorrect number of required positional arguments\n\n")
  print_help(parser)
  stop()
} else {
  if (opt$verbose){
    cat("Project: Eugen -- Optimial Charge Location Estimator\n")
    cat("Author: Aaron T. Frank\n")
    cat(sprintf("%s\n",date()))
  }
  
  #user functions
  source("include/library.R")

  # get arguments
  pdbfile <- arguments$args[1]
  charge_state <- abs(as.numeric(arguments$args[2]))
  
  # get options
  atomSelection <- opt$atomSelection
  direct <- opt$direct
  iterations <- opt$iterations
  population <- opt$population
  output <- opt$output
  verbose <- opt$verbose

  # read in PDB file
  r <- read.pdb(file = pdbfile)
  # make selection
  s <- atom.select(r, elety = atomSelection)    
  
  # optimize
  if (direct){
    charged <- direct_method(r, s, charge_state)
  } else {
    charged <- ga_method(r, s, charge_state, verbose)
  }

  # print out solution
  cat(sprintf("Charged Residue:\n"))
  print(charged)
  
  # write out pdb with residue modified to reflect predictions
  resids <- mutate(charged, r)
  write.pdb(pdb = r, file = output, resid = resids)
}
