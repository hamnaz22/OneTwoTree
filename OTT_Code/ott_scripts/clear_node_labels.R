# Clear node labels 
# user variable is : working_dir, tree, tree_out
# read tree 
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}

require(ape)
tree = read.tree(tree)
tree$node.label = NULL
setwd(working_dir)
tree = multi2di(tree)
write.tree(tree,file=tree_out)