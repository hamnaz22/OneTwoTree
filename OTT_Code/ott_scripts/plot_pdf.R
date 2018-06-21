# Creates a PDF that shows the tree. 

# read tree, working_dir, output
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}

require(ape)
setwd(working_dir)
tree = try(read.tree(tree_file))
if (inherits(tree,"try-error")){
	tree = read.nexus(tree_file)
}

t = "phylogram"
if (length(tree$tip.label) > 150){
  t = "fan"
}

pdf(output)
plot(tree,cex=0.4,type=t)
dev.off()