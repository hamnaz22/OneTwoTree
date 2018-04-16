require(ape)

parsemb <- function(dir_path, OUT_GROUP=NULL, burnin_frac = 0.25){

  # find most likely tree
  tab1 = read.table(paste(dir_path,"mb.out.run1.p",sep=""),skip = 2,sep = "\t")
  tab2 = read.table(paste(dir_path,"mb.out.run2.p",sep=""),skip = 2,sep = "\t")

  # rm burnin
  tab1 = tab1[(round(dim(tab1)[1]*burnin_frac)+1):dim(tab1)[1],]
  tab2 = tab2[(round(dim(tab2)[1]*burnin_frac)+1):dim(tab2)[1],]

  if (max(tab1[,2]) > max(tab2[,2])){
    max_ln_index = which(max(tab1[,2])==tab1[,2])[1]
    t.file = paste(dir_path,"mb.out.run1.t",sep="")
  }else {
    max_ln_index = which(max(tab2[,2])==tab2[,2])[1]
    t.file = paste(dir_path,"mb.out.run2.t",sep="")
  }

  tree_list = read.nexus(t.file)
  tree_list = tree_list[(round(length(tree_list)*burnin_frac)+1):length(tree_list)]

  most_likely = tree_list[[max_ln_index]]
  if (!is.null(OUT_GROUP)){
    most_likely = drop.tip(most_likely,OUT_GROUP)
  }

  write.tree(most_likely, file = paste(dir_path,"parsemb_map_tree.tre",sep=""), append = FALSE)

  # write trees
  t.file1 = paste(dir_path,"mb.out.run1.t",sep="")
  t.file2 = paste(dir_path,"mb.out.run2.t",sep="")

  t1=read.nexus(t.file1)
  t2=read.nexus(t.file2)

  t1 = t1[(round(length(t1)*burnin_frac)+1):length(t1)]
  t2 = t2[(round(length(t2)*burnin_frac)+1):length(t2)]

  all_trees = c(t1,t2)
  all_trees = all_trees[sample(length(all_trees),length(all_trees))] # shuffle

  if (!is.null(OUT_GROUP)){
    dropped = rmtree(length(all_trees), length(most_likely$tip.label))
    for (i in 1:length(all_trees)){
      dropped[[i]] =  drop.tip(all_trees[[i]],OUT_GROUP)
    }
    dropped = dropped[sample(length(dropped),length(dropped))]
    write.tree(dropped, file = paste(dir_path,"parsemb_trees.tre",sep=""), append = FALSE)
  } else{
    write.tree(all_trees, file = paste(dir_path,"parsemb_trees.tre",sep=""), append = FALSE)
  }

}
