# args: working_dir, id, OUT_GROUP, scripts_dir
#require(ape) # for nexus newick conversion

### PARSE INPUT ARGUMENTS
args <- commandArgs( TRUE )
for( i in 1:length(args) ){
  eval( parse( text = args[[i]] ) )
}

setwd(working_dir)
source(paste(scripts_dir,"sa_functions.R",sep="/"))

# read input files
nc= readLines(paste(working_dir,"mb_config.nex",sep=""))
seq_file=paste(working_dir,"mb_final_seq.nex",sep="")


# post process run 1
con_file = paste(working_dir,"mb.out.con.tre",sep="")  
if (file.exists(con_file)){   
  parsemb(working_dir,OUT_GROUP=OUT_GROUP)
}else{
  cat(paste(id," run1 did not finish",sep=""), file= "parsemb.log",append = TRUE,sep = "\n")
  print(paste(id," run1 did not finish",sep=""))
}








