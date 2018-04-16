# This script converts a constraint tree (or any tree for that matter) to a MrBayes block, to be added to the config file
library(ape)
library(paleotree)

# GET ARGUMENTS - (1) working dir, called "dir" (2) tree called "tree" in NEWICK format
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}

#require(paleotree) # load package
setwd(dir) # set working directory 
tree = read.tree(tree) # read constraint tree
block = createMrBayesConstraints(tree) # convert to MrBayes block
writeLines(block,"constraint_block.txt") # output - MrBayes block in txt format