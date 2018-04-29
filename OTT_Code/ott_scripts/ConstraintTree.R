require(diversitree)
require(paleotree)
require(VGAM) #for probit
require(diptest)
require("R.paleotree")

#-Erase # args: working_dir, id, OUT_GROUP
# args: tree_file, out_file_mb_blk

### PARSE INPUT ARGUMENTS
args <- commandArgs( TRUE )
for( i in 1:length(args) ){
  eval( parse( text = args[[i]] ) )
}

setwd(working_dir)

Tree = read.tree(tree_file) # if format other than newick we should modify this command
createMrBayesConstraints(tree)







