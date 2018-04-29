# Creates a PDF that shows the tree. 

# read tree, working_dir, output
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}

require(ape)
setwd(working_dir)
tree = read.tree(tree)

t = "phylogram"
if (length(tree$tip.label) > 200){
  t = "fan"
}

pdf(output)
plot(tree,label.offset=1, cex=0.4,type=t)
dev.off()