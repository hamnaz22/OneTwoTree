
require(diversitree)
require(hash)
require(VGAM) #for probit
require(diptest)
require("R.utils")

assign("mar_vec", c(5.1, 4.1, 4.1, 4.1), envir = .GlobalEnv)
assign("inset", c(-0.3,0), envir = .GlobalEnv)
assign("legend_cex", 0.6, envir = .GlobalEnv)

get_genus_data <- function(genera,ANALYSIS){  
  dat = matrix(NA,ncol=7,nrow=length(genera))
  load("species_per_genus.RData")
  for (t in 1:length(genera)){      
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]  
    tree_file = paste(trees_dir,"/", genus, ".tre", sep="")  
    if (file.exists(tree_file)){      
      states = get_states(genus_file,ANALYSIS) 
      
      phy_list <- read.tree(tree_file)
      species_num = species_per_genus$species_num[which(species_per_genus$Genus==genus)]
      ans = merge_tree_states(phy_list[[1]],states,species_num)      
      phy = ans[[1]]
      states = ans[[2]]
      df = ans[[3]]
      sampling.f = ans[[4]]       
      
      dat[t,] = c(df$n.sp.tree.inter, df$pct.non.na, 
                  df$pct.0, df$pct.1, df$n.0, df$n.1, df$sampling.f)                
    }
  }  
  colnames(dat)=c("n.sp", "pct.non.na", "pct.0", 
                  "pct.1", "n.0", "n.1","sampling.f")  
  dat = data.frame(genus=gsub(".csv","",genera),dat)  
  return(dat)
}

clean_data <- function(input_file, clean_file="out.csv"){
  
  f <-read.csv(input_file, header = TRUE, stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("NA","") )
  f = f[which(!is.na(f$Sexual.System) & f$ThePlantListV1.1..1.accepted.at.species.level..1.recognized.synonym.0.other.==1),]
  
  print(head(f))
  print(length(unique(f$Sexual.System[which(!is.na(f$Sexual.System))])))
  print(sort(table(f$Sexual.System)))
  
  sally.sp = paste(f$Genus,f$species)
  sally.sp = gsub(" ","_",sally.sp)
  sally.sp = gsub("\\.","",sally.sp)
  sally.sp = gsub("-","_",sally.sp)
  
  clean_df=data.frame(genus=f$Genus, species=sally.sp, Sexual.System=f$Sexual.System)
  
  write.csv(clean_df, file = clean_file, row.names = FALSE)
  
  genera = unique(f$Genus)
  for (i in 1:length(genera)){    
    genus=genera[i]  
    ind = which(f$Genus==genus)
    if (length(ind)>1){
      temp_df = f[ind,]
      write.csv(temp_df, file = paste("Genera_Cleaned/",genus,".csv",sep=""), row.names = FALSE)
    }
  }
    
}


get_states <- function(genus_file,ANALYSIS){
  # state info
  state_hash = get_hash(ANALYSIS)  
  
  dat1 <- read.csv(paste("Genera_Cleaned/", genus_file, sep=""))
  dat1 <- dat1[, c("Name", "SexSyst")]  
  dat1$Name <- gsub(" ", "_", dat1$Name)  
  i <- which(dat1$SexSyst == "")
  if (length(i) > 0)
    dat1$SexSyst[i] <- NA  
  
  sexSyst = rep(NA,length(dat1$SexSyst))
  for (i in 1:length(sexSyst)){
    if (!is.na(dat1$SexSyst[i])){
      s = clean_state(dat1$SexSyst[i])
      s1 = resolve_state(s,state_hash)
      if (is.null(s1)) s1 =NA
      sexSyst[i] = s1
      #print(c(dat1$SexSyst[i],s,s1))
    } 
  }
  
  # merge sub species and variants
  unique_sp = unique(dat1$Name)
  unique_sexSyst = rep(NA,length(unique_sp))
  for (i in 1:length(unique_sp)){
    n = unique_sp[i]
    ind = which(dat1$Name==n)
    if (length(ind) == 1){
      unique_sexSyst[i] = sexSyst[ind] 
    } else {
      #print(i)
      ambig_states = sexSyst[ind]
      #print(ambig_states)
      if (is.element(1,ambig_states) && is.element(0,ambig_states)){
        unique_sexSyst[i] = NA
      } else if (is.element(1,ambig_states)){
        unique_sexSyst[i] = 1
      } else if (is.element(0,ambig_states)){
        unique_sexSyst[i] = 0
      } else {
        unique_sexSyst[i] = NA
      }
      #print(unique_sexSyst[i]) 
    }    
  }
  states <- unique_sexSyst
  names(states) <- unique_sp
  
  return(states)
}

# 
# merge_tree_states <- function(tree_file,states){
#   df = data.frame(n.sp=NA,n.sp.tree=NA,n.sp.tree.inter=NA,pct.non.na=NA,pct.0=NA,pct.1=NA)
#   
#   #   tree_file2 = gsub("-fin-","-1st-",tree_file) # moshe needs to change this issue
#   #   if (file.exists(tree_file) || file.exists(tree_file2)){  
#   #     
#   #     if (file.exists(tree_file)) {
#   #       phy <- read.tree(tree_file)
#   #     } else{
#   #       phy <- read.tree(tree_file2)
#   #     }
#   if (file.exists(tree_file)){ 
#     
#     phy_list <- read.tree(tree_file) # make sure the input file contains more than 1 tree
#     
#     phy = phy_list[[1]]
#     
#     df$n.sp.tree = length(phy$tip.label) 
#     
#     if (length(phy$tip.label) == 0){
#       print(paste("empty tree: ",tree_file))
#       return(NULL)
#     } 
#     
#     phy <- multi2di(phy, random = TRUE)
#     
#     
#     #print(sum(phy$edge.length))
#     
#     phy$edge.length<-
#       phy$edge.length/max(branching.times(phy))
#     
#     #print(sum(phy$edge.length))
#     
#     # drop tips with subsp. or   var.
#     if (length(grep("subsp.", phy$tip.label)) > 0) phy = drop.tip(phy,phy$tip.label[grep("subsp.", phy$tip.label)])
#     if (length(grep("var.", phy$tip.label)) > 0) phy = drop.tip(phy,phy$tip.label[grep("var.", phy$tip.label)])
#     
#     df$n.sp.tree.inter = length(phy$tip.label) 
#     
#     if (length(phy$tip.label) == 0){
#       print(paste("empty tree after subsp. removal: ",tree_file))
#       return(NULL)
#     } 
#     
#     # sampling calc:
#     max.taxa = length(union(phy$tip.label,names(states)))
#     sf = length(phy$tip.label)/max.taxa
#     sampling.f = c(sf,sf)    
#     
#     # drop data not on tree
#     states <- states[which(is.element(names(states),phy$tip.label))]
#     
#     # add the taxa not in states to states
#     d1 = setdiff(phy$tip.label,names(states))
#     states2 = rep(NA,length(d1))
#     names(states2) = d1
#     states = c(states2,states)
#     
#     # reorder according to tree tips
#     states <- states[phy$tip.label]
#     
#     df$pct.non.na = round(100*length(which(!is.na(states)))/length(states)) 
#     df$pct.0 = round(100*length(which(states==0))/length(which(!is.na(states)))) 
#     df$pct.1 = round(100*length(which(states==1))/length(which(!is.na(states)))) 
#     df$n.0 = length(which(states==0))
#     df$n.1 = length(which(states==1))
#     df$sampling.f = sf
#     
#     return(list(phy,states,df,sampling.f,phy_list))
#     
#   } else {
#     print(paste("no tree: ",tree_file))
#     return(NULL)    
#   }
# }

merge_tree_states <- function(phy,states,species_num=1){
  df = data.frame(n.sp=NA,n.sp.tree=NA,n.sp.tree.inter=NA,pct.non.na=NA,pct.0=NA,pct.1=NA)
  
  df$n.sp.tree = length(phy$tip.label) 
  
  if (length(phy$tip.label) == 0){
    print(paste("empty tree: ",tree_file))
    return(NULL)
  } 
  
  # rescale
  phy <- multi2di(phy, random = TRUE)
  
  phy$edge.length<-
    phy$edge.length/max(branching.times(phy))
  
#   rescale.function = rescale(phy,"depth")
#   phy = rescale.function(1)
  
  # drop tips with subsp. or   var.
  if (length(grep("subsp.", phy$tip.label)) > 0) phy = drop.tip(phy,phy$tip.label[grep("subsp.", phy$tip.label)])
  if (length(grep("var.", phy$tip.label)) > 0) phy = drop.tip(phy,phy$tip.label[grep("var.", phy$tip.label)])
  
  df$n.sp.tree.inter = length(phy$tip.label) 
  
  if (length(phy$tip.label) == 0){
    print(paste("empty tree after subsp. removal: ",tree_file))
    return(NULL)
  } 
  
  # sampling calc:
  #max.taxa = length(union(phy$tip.label,names(states)))
  df$species_num = species_num
  sf = length(phy$tip.label)/species_num
  sampling.f = c(sf,sf)    
  
  # drop data not on tree
  states <- states[which(is.element(names(states),phy$tip.label))]
  
  # add the taxa not in states to states
  d1 = setdiff(phy$tip.label,names(states))
  states2 = rep(NA,length(d1))
  names(states2) = d1
  states = c(states2,states)
  
  # reorder according to tree tips
  states <- states[phy$tip.label]
  
  df$pct.non.na = round(100*length(which(!is.na(states)))/length(states)) 
  df$pct.0 = round(100*length(which(states==0))/length(which(!is.na(states)))) 
  df$pct.1 = round(100*length(which(states==1))/length(which(!is.na(states)))) 
  df$n.0 = length(which(states==0))
  df$n.1 = length(which(states==1))
  df$sampling.f = sf
  
  return(list(phy,states,df,sampling.f))     
}

server_commands <- function(dir_name,to_run,to_server_dir,run_file,NSTEPS=1000,server="jekyl",time_out=864000){ # 864000 = 10 days
  
  dir.create(paste(to_server_dir,"/",dir_name,sep=""))
  dir.create(paste(to_server_dir,"/",dir_name,"/temp",sep=""))
  file.copy(run_file,paste(to_server_dir,"/",dir_name,sep=""))
  
  command_file = paste(to_server_dir,"/",dir_name,"/command.txt",sep="")
  file.create(command_file)
  
  dir_path = paste( getwd(),"/",to_server_dir,"/",dir_name,"/",sep="")
  
  cat("/share/apps/R301/bin/R CMD BATCH '--args working_dir=\"",
      file= command_file,append = FALSE)
  cat(dir_path,
      file= command_file,append = TRUE)
  cat("\" NSTEPS=",
      file= command_file,append = TRUE)
  cat(NSTEPS,
      file= command_file,append = TRUE)
  cat("' ",
      file= command_file,append = TRUE)
  cat(dir_path,
      file= command_file,append = TRUE)
  cat(run_file,
      file= command_file,append = TRUE)
  cat(" ",
      file= command_file,append = TRUE)
  cat(dir_path,
      file= command_file,append = TRUE)
  cat(run_file,
      file= command_file,append = TRUE)  
  cat("out\tr",
      file= command_file,append = TRUE)
  cat(dir_name,
      file= command_file,append = TRUE)  
  
  if (server=="jekyl"){
    cat("perl /groups/itay_mayrose/nivsabath/Scripts/run_cmds_in_q_WithNameInCmd_jekyl.pl ",
        file= to_run,append = TRUE)    
  }else{ # lecs
    cat("perl /groups/itay_mayrose/nivsabath/Scripts/run_cmds_in_q_WithNameInCmd.pl ",
        file= to_run,append = TRUE)    
  }
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("command.txt ",
      file= to_run,append = TRUE)
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("temp/ itaym 1 yes ",
      file= to_run,append = TRUE)
  cat(time_out,
      file= to_run,append = TRUE)
  cat(" r",
      file= to_run,append = TRUE)
  cat(dir_name,
      file= to_run,append = TRUE)  
  cat("\n",
      file= to_run,append = TRUE)  
  
}


server_commands2 <- function(dir_name,to_run,to_server_dir,run_file,NSTEPS=1000,server="jekyl",time_out=864000){ # 864000 = 10 days
  
  file.copy(run_file,paste(to_server_dir,"/",dir_name,sep=""))  
  command_file = paste(to_server_dir,"/",dir_name,"/command.txt",sep="")
  file.create(command_file)  
  dir_path = paste( getwd(),"/",to_server_dir,"/",dir_name,"/",sep="")
  
  cat("/share/apps/R301/bin/R CMD BATCH '--args working_dir=\"",
      file= command_file,append = FALSE)
  cat(dir_path,
      file= command_file,append = TRUE)
  cat("\" NSTEPS=",
      file= command_file,append = TRUE)
  cat(NSTEPS,
      file= command_file,append = TRUE)
  cat("' ",
      file= command_file,append = TRUE)
  cat(dir_path,
      file= command_file,append = TRUE)
  cat(run_file,
      file= command_file,append = TRUE)
  cat(" ",
      file= command_file,append = TRUE)
  cat(dir_path,
      file= command_file,append = TRUE)
  cat(run_file,
      file= command_file,append = TRUE)  
  cat("out\tr",
      file= command_file,append = TRUE)
  cat(dir_name,
      file= command_file,append = TRUE)  
  
  if (server=="jekyl"){
    cat("perl /groups/itay_mayrose/nivsabath/Scripts/run_cmds_in_q_WithNameInCmd_jekyl.pl ",
        file= to_run,append = TRUE)    
  }else{ # lecs
    cat("perl /groups/itay_mayrose/nivsabath/Scripts/run_cmds_in_q_WithNameInCmd.pl ",
        file= to_run,append = TRUE)    
  }
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("command.txt ",
      file= to_run,append = TRUE)
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("temp/ itaym 1 yes ",
      file= to_run,append = TRUE)
  cat(time_out,
      file= to_run,append = TRUE)
  cat(" r",
      file= to_run,append = TRUE)
  cat(dir_name,
      file= to_run,append = TRUE)  
  cat("\n",
      file= to_run,append = TRUE)  
  
}

run_mcmc <- function(phy,states,sampling.f=1,nsteps=1000) {   
  p <- starting.point.bisse(phy)
  
  if (p[3]==0 || p[1]/p[3]>10){
    p[3:4] = p[1]/10
  }
  
  lik <- make.bisse(phy, states,sampling.f = sampling.f)  
  
  priorrate=1/(2*p);               
  prior = make.prior.exponential(priorrate)  
  w <- c(1,1,1,1,0.4,0.4) 
  
  samples <- mcmc(lik, starting.p, nsteps=NSTEPS, w=w, lower=0,
                  prior=prior,print.every=100, control=list(backend="CVODES")) 
  
  samples <- subset(samples, i > length(samples[,1])/10)
  
  samples$div0 = samples$lambda0-samples$mu0
  samples$div1 = samples$lambda1-samples$mu1
  
  return(samples)
}


run_mcmc2 <- function(lik,p,fit,nsteps,name,phy) { 
  
  tree.length <- max(branching.times(phy))
  r.est <- log(length(phy$tip.label))/tree.length
  priorpar <- r.est*2
  priorrate <- 1/priorpar
  prior = make.prior.exponential(priorrate)
  
  print(c(1 / (2 * (p[1] - p[3])),priorrate))
  
  #prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))
  tmp <- mcmc(lik, fit$par, nsteps=100, prior=prior,lower=0, w=rep(1, 6), print.every=0)
  w <- diff(sapply(tmp[2:7], range))
  samples <- mcmc(lik, fit$par, nsteps=nsteps, w=w, lower=0, prior=prior,print.every=100, control=list(backend="CVODES"))
  samples <- subset(samples, i > length(samples[,1])/10)
  samples$div0 = samples$lambda0-samples$mu0
  samples$div1 = samples$lambda1-samples$mu1
  save(samples, file = paste("mcmc2/",name,".mcmc",sep=""))
  return(samples)
}


run_mcmc3 <- function(lik,p,fit,nsteps,name,phy) { 
  
  tree.length <- max(branching.times(phy))
  r.est <- log(length(phy$tip.label))/tree.length
  priorpar <- r.est*2
  priorrate <- 1/priorpar
  
  prior = make.prior.exponential(priorrate)  
  
  tmp <- mcmc(lik, fit$par, nsteps=100, prior=prior,lower=0, w=rep(1, 6), print.every=0)
  w <- diff(sapply(tmp[2:7], range))
  samples <- mcmc(lik, fit$par, nsteps=nsteps, w=w, lower=0, prior=prior,print.every=100, control=list(backend="CVODES"))
  samples <- subset(samples, i > length(samples[,1])/10)
  samples$div0 = samples$lambda0-samples$mu0
  samples$div1 = samples$lambda1-samples$mu1
  save(samples, file = paste("mcmc/",name,".mcmc",sep=""))

  
  prior = make.prior.exponential(priorrate*2)  
  
  tmp <- mcmc(lik, fit$par, nsteps=100, prior=prior,lower=0, w=rep(1, 6), print.every=0)
  w <- diff(sapply(tmp[2:7], range))
  samples <- mcmc(lik, fit$par, nsteps=nsteps, w=w, lower=0, prior=prior,print.every=100, control=list(backend="CVODES"))
  samples <- subset(samples, i > length(samples[,1])/10)
  samples$div0 = samples$lambda0-samples$mu0
  samples$div1 = samples$lambda1-samples$mu1
  save(samples, file = paste("mcmc/",name,"2.mcmc",sep=""))

  prior = make.prior.exponential(priorrate*4)  
  
  tmp <- mcmc(lik, fit$par, nsteps=100, prior=prior,lower=0, w=rep(1, 6), print.every=0)
  w <- diff(sapply(tmp[2:7], range))
  samples <- mcmc(lik, fit$par, nsteps=nsteps, w=w, lower=0, prior=prior,print.every=100, control=list(backend="CVODES"))
  samples <- subset(samples, i > length(samples[,1])/10)
  samples$div0 = samples$lambda0-samples$mu0
  samples$div1 = samples$lambda1-samples$mu1
  save(samples, file = paste("mcmc/",name,"4.mcmc",sep=""))
  
}


plot_dens <- function(d0,d1,xlab,title="",xlim = NULL) { 
  d0 <- density(d0)
  d1 <- density(d1)  
  if (is.null(xlim)) xlim <- c(0,max(d0$x,d1$x))
  ylim <- c(0,max(d0$y,d1$y))
  plot(d0, xlim = xlim, ylim = ylim, xlab = xlab, main = title, cex.lab = 1.2)
  #put our density plots in
  polygon(d0, density = -1, col = rgb(1,0,0,0.2))
  polygon(d1, density = -1, col = rgb(0,0,1,0.2))
  ## add a legend in the corner
}

plot_dens2 <- function(d0,d1,xlab,title="",xlim = NULL) { 
  d0 <- density(d0)
  d1 <- density(d1)
  if (is.null(xlim)) xlim <- range(d0$x,d1$x)
  ylim <- range(0,d0$y, d1$y)
  plot(d0, xlim = xlim, ylim = ylim, xlab = xlab, main = title, cex.lab = 1.2)
  #put our density plots in
  polygon(d0, density = -1, col = rgb(1,0,0,0.2))
  polygon(d1, density = -1, col = rgb(0,0,1,0.2))
  ## add a legend in the corner
}

plot_mcmc_results <- function(samples,name,xlim = NULL) {   
  
  old.par <- par( no.readonly = TRUE ) 
  
  par(oma=c(1,1,4,1))
  
  par(mar=c(5,6,4,2))
  
  par(mfrow=c(2,2))
  
  plot_dens(samples$lambda0,samples$lambda1,"Speciation rate","",xlim)
  legend('topright',c(expression(lambda[0]),expression(lambda[1])), fill = c(rgb(1,0,0,0.2), rgb(0,0,1,0.2)), bty = 'n', border = NA)
  mtext("a", side = 3, line = 1, adj = 0, cex = 1)
  
  plot_dens(samples$mu0,samples$mu1,"Extinction rate","")
  legend('topright',c(expression(mu[0]),expression(mu[1])), fill = c(rgb(1,0,0,0.2), rgb(0,0,1,0.2)), bty = 'n', border = NA)
  mtext("b", side = 3, line = 1, adj = 0, cex = 1)
  
  plot_dens(samples$div0,samples$div1,"Diversification rate","")
  legend('topright',c(expression(r[0]),expression(r[1])), fill = c(rgb(1,0,0,0.2), rgb(0,0,1,0.2)), bty = 'n', border = NA)
  mtext("c", side = 3, line = 1, adj = 0, cex = 1)
  
  plot_dens(samples$q01,samples$q10,"Rate of change in character state","")
  legend('topright',c(expression(q[01]),expression(q[10])), fill = c(rgb(1,0,0,0.2), rgb(0,0,1,0.2)), bty = 'n', border = NA)
  mtext("d", side = 3, line = 1, adj = 0, cex = 1)
  
  title( name, outer = TRUE )
 
  
  par( old.par )
  
  print("MCMC Results:")
  print(paste("lambda0 > lambda1: ",100*length(which(samples$lambda0>samples$lambda1))/length(samples$lambda0),"%",sep=""))
  print(paste("mu0 > mu1: ",100*length(which(samples$mu0>samples$mu1))/length(samples$lambda0),"%",sep=""))  
  print(paste("div0 > div1: ",100*length(which(samples$div0>samples$div1))/length(samples$lambda0),"%",sep=""))
  print(paste("q01 > q10: ",100*length(which(samples$q01>samples$q10))/length(samples$lambda0),"%",sep=""))
  
  print("mean values:")
  round(c(mean(samples$lambda0),mean(samples$lambda1),mean(samples$mu0),mean(samples$mu1),mean(samples$q01),mean(samples$q10)),4)
}


my_profile_plot <- function(data,xlab,legend.txt,r.est=NULL){ 
  col.fill <- c("#FF000066", "#0000FF66")
  col.line <- c("red", "blue")
  ylab <- "Probability density"
  
  profiles.plot(data, col.line, col.fill, las=1, xlab=xlab, ylab=ylab)
  legend("topright", legend.txt, fill=col.fill,bty="n")
#   x<- seq(0,60,0.5)
#   yfunc <- 1/(2*r.est)*exp(-1/(2*r.est)*x)
#   lines(x,yfunc,type="l")
  
}


plot_mcmc_results_new <- function(samples,name="",p=NULL,xlim = NULL) {   
  
  old.par <- par( no.readonly = TRUE ) 
  
  par(oma=c(1,1,4,1))
  
  par(mar=c(5,5,2,2))
  
  par(mfrow=c(2,2))
  
  
  my_profile_plot(samples[c("lambda0", "lambda1")],"Speciation rate",c("lam0", "lam1"),p[1])
  mtext("a", side = 3, line = 1, adj = 0, cex = 1)
  
  my_profile_plot(samples[c("mu0", "mu1")],"Extinction rate",c("mu0", "mu1"),p[3])
  mtext("b", side = 3, line = 1, adj = 0, cex = 1)
  
  my_profile_plot(samples[c("div0", "div1")],"Diversification rate",c("r0", "r1"),p[1]-p[3])
  mtext("c", side = 3, line = 1, adj = 0, cex = 1)
  
  my_profile_plot(samples[c("q01", "q10")],"Transition rate",c("q01", "q10"),p[5])
  mtext("d", side = 3, line = 1, adj = 0, cex = 1)
      
 # title(name)
  title( name,line=1.2, outer = TRUE )
  par( old.par )
  
  print("MCMC Results:")
  print(paste("lambda0 > lambda1: ",round(100*length(which(samples$lambda0>samples$lambda1))/length(samples$lambda0)),"%",sep=""))
  print(paste("mu0 > mu1: ",round(100*length(which(samples$mu0>samples$mu1))/length(samples$lambda0)),"%",sep=""))
  print(paste("div0 > div1: ",round(100*length(which(samples$div0>samples$div1))/length(samples$lambda0)),"%",sep=""))
  print(paste("q01 > q10: ",round(100*length(which(samples$q01>samples$q10))/length(samples$lambda0)),"%",sep=""))
  
  print("mean values:")
  round(c(mean(samples$lambda0),mean(samples$lambda1),mean(samples$mu0),mean(samples$mu1),mean(samples$q01),mean(samples$q10)),4)
}



# ## Read in the full traits file.  Keep only the name and sex system fields.
# ## Create composite name fields, for later convenience.
# get.sexsyst <- function(genus_file)
# {
#   dat <- read.csv(paste("Genera_Cleaned/", genus_file, sep=""))
#   
#   dat <- dat[, c("Name", "SexSyst")]
#     
#   dat$Name <- gsub(" ", "_", dat$Name)
#   
#   i <- which(dat$SexSyst == "")
#   if (length(i) > 0)
#     dat$SexSyst[i] <- NA
#   
#   return(dat)
# }





plot_trait <- function(tree,trait,reg_str,title="") { 
  
  tree$tip.state = data.frame(reg=trait$Reg,x=round(trait$X,1),
                              row.names=tree$tip.label)
  
  for (i in 1:length(tree$tip.label)){
    if (!is.na(tree$tip.state$x[i])){
      tree$tip.label[i] = paste(tree$tip.state$x[i],tree$tip.label[i])
    }    
  }
  
  t = "phylogram"
  if (length(tree$tip.label) > 200){
    t = "fan"
  }
  
  col <- c("blue", "red")
  statecols <- c("lightblue", "blue")
  plot(tree, label.offset=1, cex=0.4,tip.color=col[tree$tip.state$reg+1],type=t) 
  #tiplabels(col=col[tree$tip.state$reg+1], pch=19,cex=0.6)
  #tiplabels(col=col[tree$tip.state$reg+1], pch=19,cex=0.6)
  legend("topright", reg_str,col=col,pch=19)
  mtext(title)
}



get_hash <-function(ANALYSIS){
  if (ANALYSIS==1){
    state_hash = hash( list("hermaphrodite"=0, 
                            "monoecy"=0, "andromonoecy"=0, "gynomonoecy"=0, "polygamomonoecy"=0, 
                            "gynodioecy"=0,"gynodioecious"=0, "androdioecy"=0, "polygamodioecy"=0,                            
                            
                            "dioecy"=1,                             
                            
                            "apomictic"=NA,     
                            "other"=NA) )
  }
  
  if (ANALYSIS==2){
    state_hash = hash( list("hermaphrodite"=0, 
                            "monoecy"=0, "andromonoecy"=0, "gynomonoecy"=0, "polygamomonoecy"=0, 
                            
                            "gynodioecy"=1,"gynodioecious"=1, "androdioecy"=1, "polygamodioecy"=1,                                                 
                            "dioecy"=1,                                                        
                            
                            "apomictic"=NA,     
                            "other"=NA) )    
  }
  
  if (ANALYSIS==11){
    state_hash = hash( list("hermaphrodite"=0, 
                            
                            "monoecy"=1, "andromonoecy"=1, "gynomonoecy"=1, "polygamomonoecy"=1, 
                            "gynodioecy"=1,"gynodioecious"=1, "androdioecy"=1, "polygamodioecy"=1,                                                 
                            
                            "dioecy"=2,                                                        
                            
                            "apomictic"=NA,     
                            "other"=NA) )    
  }
  if (ANALYSIS==21){
    state_hash = hash( list("hermaphrodite"=0, 
                            
                            "monoecy"=1, "andromonoecy"=1, "gynomonoecy"=1, "polygamomonoecy"=1, 
                            
                            "gynodioecy"=2,"gynodioecious"=2, "androdioecy"=2, "polygamodioecy"=2,                                                 
                            
                            "dioecy"=2,                                                        
                            
                            "apomictic"=NA,     
                            "other"=NA) )    
  }
  
  return(state_hash)
}

clean_state <-function(state){
  state = gsub(",", "|",state) 
  state = gsub("sequentially", "",state) 
  state = gsub("in part", "",state)
  state = gsub("sometimes", "",state)
  state = gsub(" ", "",state)
  sp = strsplit(state,"[|]")
  if (length(sp[[1]])>2){
    sp=sp[[1]]
    state = paste(sp[1],"|",sp[2],sep="")
  }
  return(state)
}

resolve_state <-function(state,state_hash){
  if (is.na(state)) return(NA)
  sp = strsplit(state,"[|]")
  if (length(sp[[1]])==1) return(state_hash[[state]])
  if (length(sp[[1]])==2){
    sp=sp[[1]]
    s1 = state_hash[[sp[1]]]
    s2 = state_hash[[sp[2]]]
    if (is.null(s1) || is.na(s1)) s1 = -1
    if (is.null(s2) || is.na(s2)) s2 = -1
    if (s1==0 && s2==0) return(0)
    if (s1==1 && s2==1) return(1)
    if (s1==-1 && s2==-1) return(NA)
    if ((s1==1 && s2==0) || (s1==0 && s2==1)) return(NA)
    if ((s1==1 && s2==-1) || (s1==-1 && s2==1)) return(1)
    if ((s1==0 && s2==-1) || (s1==-1 && s2==0)) return(0)
  } 
  else{
    warning("error in state: ",state)
    return(NA)
  }   
}



prepare_to_server <- function(trees_dir,genera,to_run = "to_run.txt",to_server_dir = "to_server",
                              run_file = "run_on_server.R",NSTEPS = 2000,ANALYSIS=1,factor_vec=c(1,2,4)){  
  
  load("species_per_genus.RData")
  for (t in 1:length(genera)){
    
    # file name
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    
    tree_file = paste(trees_dir,"/", genus, ".tre", sep="")  
    
    states = get_states(genus_file,ANALYSIS) 
    names(states)=gsub("-","_",names(states))
    species_num = species_per_genus$species_num[which(species_per_genus$Genus==genus)]
    ans = merge_tree_states(tree_file,states,species_num)
    
    if (!is.null(ans)){
      phy = ans[[1]]
      states = ans[[2]]
      df = ans[[3]]
      sampling.f = ans[[4]] 
      
      if (df$n.sp.tree.inter >= 10 && df$n.0 >= 2 && df$n.1 >= 2){
        for (FACTOR in factor_vec){
          name = paste(genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,sep="")
          server_commands(name,to_run,to_server_dir,run_file,NSTEPS)    
          file.create("Data.RData")
          save(phy, states, name, FACTOR,df,sampling.f,file="Data.RData")
          file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
          file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))
        }
      } 
    }
  }
}


prepare_to_server_new <- function(trees_dir,genera,to_run = "to_run.txt",to_server_dir = "to_server",
                                  run_file = "run_mcmc.R",NSTEPS = 2000,ANALYSIS=1,factor_vec=c(1,2,4),
                                  prior.type="exp.3.param.mle",start.type="best.point",seconed=FALSE,root.p=NULL){  
  load("species_per_genus.RData")
  for (t in 1:length(genera)){
    
    # file name
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    
    tree_file = paste(trees_dir,"/", genus, ".tre", sep="")  
    
    states = get_states(genus_file,ANALYSIS) 
    names(states)=gsub("-","_",names(states))
    species_num = species_per_genus$species_num[which(species_per_genus$Genus==genus)]
    ans = merge_tree_states(tree_file,states,species_num)
    
    if (!is.null(ans)){
      phy = ans[[1]]
      states = ans[[2]]
      df = ans[[3]]
      sampling.f = ans[[4]] 
      
      if (df$n.sp.tree.inter >= 10 && df$n.0 >= 2 && df$n.1 >= 2){
        for (FACTOR in factor_vec){          
          
          name = paste(genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,sep="")
          
          server_commands(name,to_run,to_server_dir,run_file,NSTEPS)    
          file.create("Data.RData")
          save(phy, states, name, FACTOR,sampling.f,root.p=root.p,
               prior.type=prior.type,start.type=start.type,file="Data.RData")
          file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
          file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))  
          
          if (start.type=="best.point"){
            points_name = paste("points/",genus,"_Analyais",ANALYSIS,"_points",sep="")            
            best_point=read_points(points_name,seconed=seconed)
            save(best_point,file="best_point.RData")
            file.copy("best_point.RData",paste(to_server_dir,"/",name,sep=""))
          }          
        }
      } 
    }
  }
}


prepare_to_server_trees <- function(trees_dir,genera,to_run = "to_run.txt",to_server_dir = "to_server",
                                    run_file = "run_mcmc.R",NSTEPS = 2000,ANALYSIS=1,factor_vec=c(1,2,4),
                                    prior.type="starting.point.bisse",start.type="starting.point.bisse",
                                    seconed=FALSE,root.p=NULL,N_TREES=20,RATIO=10,time_out=36000){# 10 hours  
  
  for (t in 1:length(genera)){
    
    # file name
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
        
    tree_file = paste(trees_dir,"/", genus, ".tre", sep="")
    
    load("species_per_genus.RData")
    
    if (file.exists(tree_file)){
      
      print(genus)    
      phy_list <- read.tree(tree_file)
      
      states = get_states(genus_file,ANALYSIS) 
      names(states)=gsub("-","_",names(states))
      
      species_num = species_per_genus$species_num[which(species_per_genus$Genus==genus)]
      for (tr in 1:N_TREES){  
        
        phy_original = phy_list[[tr]]
        
        ans = merge_tree_states(phy_original,states,species_num)
        
        phy = ans[[1]]
        states = ans[[2]]
        df = ans[[3]]
        sampling.f = ans[[4]]       
        
        if (df$n.sp.tree.inter >= 10 && df$n.0 >= 2 && df$n.1 >= 2){
          for (FACTOR in factor_vec){    
            
            name = paste(genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"_Tree",tr,sep="")
            
            server_commands(name,to_run,to_server_dir,run_file,NSTEPS,time_out=time_out)    
            file.create("Data.RData")
            save(phy, states, name, FACTOR,sampling.f,root.p=root.p,
                 prior.type=prior.type,start.type=start.type,RATIO,file="Data.RData")
            file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
            file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))  
            
            if (start.type=="best.point"){
              points_name = paste("points/",genus,"_Analyais",ANALYSIS,"_points",sep="")            
              best_point=read_points(points_name,seconed=seconed)
              save(best_point,file="best_point.RData")
              file.copy("best_point.RData",paste(to_server_dir,"/",name,sep=""))
            } 
          }
        }
      } 
    }
  }
}


prepare_to_server_trees_sim2 <- function(trees_dir,genera,to_run = "to_run.txt",to_server_dir = "to_server",
                                    run_file = "run_mcmc.R",NSTEPS = 2000,ANALYSIS=1,factor_vec=c(1,2,4),
                                    prior.type="starting.point.bisse",start.type="starting.point.bisse",
                                    seconed=FALSE,root.p=NULL,N_TREES=20,RATIO=10,time_out=36000,type=0){# 10 hours  
  
  for (t in 1:length(genera)){
    
    # file name
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    
    tree_file = paste(trees_dir,"/", genus, ".tre", sep="")
    
    load("species_per_genus.RData")
    
    if (file.exists(tree_file)){
      
      print(genus)    
      phy_list <- read.tree(tree_file)
      
      states = get_states(genus_file,ANALYSIS) 
      names(states)=gsub("-","_",names(states))
      species_num = species_per_genus$species_num[which(species_per_genus$Genus==genus)]
      
      sim.states = NULL
      while (is.null(sim.states)){
        ans = merge_tree_states(phy_list[[sample(N_TREES,1)]],states,species_num)
        res <- NULL;
        tryCatch({
          res <- evalWithTimeout({
            sim.states = sim.mk2(ans[[1]], ans[[2]], ans[[4]],type);
          }, timeout=1);
        }, TimeoutException=function(ex) {
          cat("Timeout. Skipping.\n");
        })  
      }
      states = sim.states
      
      for (tr in 1:N_TREES){  
        
        phy_original = phy_list[[tr]]
        
        ans = merge_tree_states(phy_original,states,species_num)
        
        phy = ans[[1]]
        states = ans[[2]]
        df = ans[[3]]
        sampling.f = ans[[4]]       
        
        if (df$n.sp.tree.inter >= 10 && df$n.0 >= 2 && df$n.1 >= 2){
          for (FACTOR in factor_vec){    
            
            name = paste(genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"_Tree",tr,sep="")
            
            server_commands(name,to_run,to_server_dir,run_file,NSTEPS,time_out=time_out)    
            file.create("Data.RData")
            save(phy, states, name, FACTOR,sampling.f,root.p=root.p,
                 prior.type=prior.type,start.type=start.type,RATIO,file="Data.RData")
            file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
            file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))  
            
            if (start.type=="best.point"){
              points_name = paste("points/",genus,"_Analyais",ANALYSIS,"_points",sep="")            
              best_point=read_points(points_name,seconed=seconed)
              save(best_point,file="best_point.RData")
              file.copy("best_point.RData",paste(to_server_dir,"/",name,sep=""))
            } 
          }
        }
      } 
    }
  }
}

prepare_to_server_trees_sim3 <- function(trees_dir,genera,to_run = "to_run_sim.txt",to_server_dir = "to_server_sim",
                                         run_file = "run_sim.R",ANALYSIS=1,FACTOR=2,                                         
                                         N_TREES=100,N_SIM=30,time_out=36000,type=1){# 10 hours    
  for (t in 1:length(genera)){    
    # file name
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]    
    tree_file = paste(trees_dir,"/", genus, ".tre", sep="")    
    load("species_per_genus.RData")
    
    if (file.exists(tree_file)){
      
      print(genus)    
      phy_list <- read.tree(tree_file)
      
      states = get_states(genus_file,ANALYSIS) 
      names(states)=gsub("-","_",names(states))
      species_num = species_per_genus$species_num[which(species_per_genus$Genus==genus)]  
      
      for (tr in 1:N_TREES){   
        print(c(t,tr))
        phy_original = phy_list[[tr]]        
        ans = merge_tree_states(phy_original,states,species_num)        
        phy = ans[[1]]
        states = ans[[2]]
        df = ans[[3]]
        sampling.f = ans[[4]]           
        for (sim in 1:N_SIM){             
          name = paste(genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"_Tree",tr,"_sim",sim,sep="")          
          server_commands(name,to_run,to_server_dir,run_file,time_out=time_out)    
          file.create("Data.RData")         
          save(phy, states, name,sampling.f,type,file="Data.RData")          
          file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
          file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))            
        }
      }
    }    
  }
}

check_status_sim3 <- function(genera,to_server_dir = "to_server_sim",
                              ANALYSIS=1,FACTOR=2, N_TREES=100,N_SIM=30){
  genera = gsub(".csv","",genera)
  res = matrix(0,nrow=length(genera), ncol=N_TREES)
  for (t in 1:length(genera)){    
    
    genus=genera[t]    
    print(genus)   
    
    for (tr in 1:N_TREES){ 
	  print(tr)	
      for (sim in 1:N_SIM){   
		#print(c(tr,sim))
        out_file = paste(to_server_dir,"/",genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"_Tree",tr,"_sim",
                         sim,"/states.RData",sep="")          
        if (file.exists(out_file)){      
          res[t,tr] = res[t,tr]+1
        }
      }
    }
  }    
  
  df = data.frame(genera,res)
  write.csv(df, file="status.csv") 
  return(res)
}


prepare_to_server_trees_sim4 <- function(trees_dir,genera,status,states_sim_dir = "sim_states",to_run = "to_run_sim.txt",to_server_dir = "to_server_sim",
                                         FACTOR=2,run_file = "run_mcmc.R",NSTEPS = 2000,ANALYSIS=1,
                                         prior.type="starting.point.bisse",start.type="starting.point.bisse",
                                         seconed=FALSE,root.p=NULL,RATIO=10,
                                         N_TREES=20,N_SIM=20,time_out=36000){   
  col <- c("red", "blue")
  pdf("examine.pdf")
  genera = gsub(".csv","",genera)
  
  for (t in 1:length(genera)){    
    
    genus=genera[t]
    
    tree_file = paste(trees_dir,"/", genus, ".tre", sep="")    
    load("species_per_genus.RData")
    print(genus) 
    phy_list <- read.tree(tree_file)
    
    species_num = species_per_genus$species_num[which(species_per_genus$Genus==genus)]  
    
    tree.ind=0
    for (tr in 1:N_TREES){   
      print(tr)
      sim.ind=1
      tree.ind = tree.ind+1
      
      out_file = paste(states_sim_dir,"/",genus,"_Analyais",ANALYSIS,"_Factor",
                       FACTOR,"_Tree",tree.ind,"_sim",sim.ind,"/states.RData",sep="") 
      while(!file.exists(out_file)){
        #print(c(tr,tree.ind,sim.ind))
        sim.ind=sim.ind+1          
        out_file = paste(states_sim_dir,"/",genus,"_Analyais",ANALYSIS,"_Factor",
                         FACTOR,"_Tree",tree.ind,"_sim",sim.ind,"/states.RData",sep="") 
        if (sim.ind==N_TREES){
          tree.ind = tree.ind+1
          sim.ind = 1
        }
        if (tree.ind>N_TREES)
          break;
      }
      if (file.exists(out_file)){
        load(out_file)
        sim.states = states
        rm(states)
        for (sim in 1:N_SIM){   
          
          phy_original = phy_list[[sim]]        
          ans = merge_tree_states(phy_original,sim.states,species_num)        
          phy = ans[[1]]
          states = ans[[2]]
          df = ans[[3]]
          sampling.f = ans[[4]]           
          
          name = paste(genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"_Tree",tr,"_sim",sim,sep="")          
          server_commands(name,to_run,to_server_dir,run_file,time_out=time_out)    
          file.create("Data.RData")         
          save(phy, states, name, FACTOR,sampling.f,root.p=root.p,
               prior.type=prior.type,start.type=start.type,RATIO,file="Data.RData")
          file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
          file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))   
          
          plot(phy, cex=0.6, font=1,label.offset=0.02,main=paste("tr:",tr,"  sim.ind:",sim.ind,"   sim",sim))
          tiplabels(col=col[states+1], pch=19)          
          legend("topleft", c("N","D"),col=col,pch=19) 
          
        }
      }
    }
  }
  dev.off()
}    


prepare_to_server_trees_sim5 <- function(trees_dir,genera,status,states_sim_dir = "sim_states",to_run = "to_run_sim.txt",to_server_dir = "to_server_sim",
                                         FACTOR=2,run_file = "run_mcmc.R",NSTEPS = 2000,ANALYSIS=1,
                                         prior.type="starting.point.bisse",start.type="starting.point.bisse",
                                         seconed=FALSE,root.p=NULL,RATIO=10,
                                         N_TREES=20,N_TREES2=20,N_SIM=100,time_out=36000){   
  col <- c("red", "blue")
  pdf("examine_sim5.pdf")
  genera = gsub(".csv","",genera)
  
  for (t in 1:length(genera)){    
    
    genus=genera[t]
    print(genus) 
    
    index_complete_sim = which(status[t,]>=N_SIM)
    
    if (length(index_complete_sim)<N_TREES){
      print(paste("skipping",genus,": not enough simulations done"))
      next; 
    }
    
    tree_file = paste(trees_dir,"/", genus, ".tre", sep="")    
    phy_list <- read.tree(tree_file)
    
    load("species_per_genus.RData")
    species_num = species_per_genus$species_num[which(species_per_genus$Genus==genus)]      
    
    for (tr in 1:N_TREES){   
      
      tree.ind = index_complete_sim[tr]      
      
      for (sim in 1:N_SIM){   
        print(c(t,tr,sim))
        out_file = paste(states_sim_dir,"/",genus,"_Analyais",ANALYSIS,"_Factor",
                         FACTOR,"_Tree",tree.ind,"_sim",sim,"/states.RData",sep="") 
        
        load(out_file)
        sim.states = states
        rm(states)
        
        for (tr2 in 1:N_TREES2){   
          
          phy_original = phy_list[[tr2+50]]        
          ans = merge_tree_states(phy_original,sim.states,species_num)        
          phy = ans[[1]]
          states = ans[[2]]
          df = ans[[3]]
          sampling.f = ans[[4]]           
          
          name = paste(genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"_Tree",tr,"_sim",sim,"_tr",tr2,sep="")          
          server_commands(name,to_run,to_server_dir,run_file,time_out=time_out)    
          file.create("Data.RData")         
          save(phy, states, name, FACTOR,sampling.f,root.p=root.p,
               prior.type=prior.type,start.type=start.type,RATIO,file="Data.RData")
          file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
          file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))   
          
          if(tr<3 && sim<3 && tr2<3){
            plot(phy, cex=0.6, font=1,label.offset=0.02,main=paste("tr:",tr,"  sim:",sim,"   tr2",tr2))
            tiplabels(col=col[states+1], pch=19)          
            legend("topleft", c("N","D"),col=col,pch=19) 
          }
        }
      }
    }
  }
  dev.off()
}  


sim.mk2 <- function(phy,states,sampling.f,type){
  # type: 0-shuffel, 1-unconst, 2-const tip ratio
  if (type==0){
    s.names = names(states)
    states = states[sample(length(states))]
    names(states) = s.names
    return(states)
  }
  
  n0 = length(which(states==0))
  n1 = length(which(states==1))
  
  pct.0 = 100*n0/(n0+n1)
  pct.1 = 100*n1/(n0+n1)
  
  # creat Bisse model
  lik <- make.bisse(phy, states, sampling.f = sampling.f)
  p <- starting.point.bisse(phy)
  
  if (p[3]==0 || p[1]/p[3]>10){
    p[3:4] = p[1]/10
  }
  
  print("find.mle")
  fit <- find.mle(lik, p)
  
  # reconstruction with Bisse
  print("asr.marginal")
  st <- asr.marginal(lik, coef(fit))
    
  print("sim.character")
  sim.states = NULL
  for (i in 1:10){  
    
    # set state of ancestral node
    x0=0
    if (runif(1, 0, 1) < st[2,1])
      x0=1
    
    sim.states <- sim.character(phy, fit$par[5:6], x0=x0, model="mk2")
    
    n0.sim = length(which(sim.states==0))
    n1.sim = length(which(sim.states==1))
    
    pct.0.sim = 100*n0.sim/(n0.sim+n1.sim)
    pct.1.sim = 100*n1.sim/(n0.sim+n1.sim)
    
    states.eq = 100*length(which(sim.states==states))/length(states)
    
    print(paste("dif %0:",abs(pct.0.sim-pct.0),"    similarity:",states.eq,"    n0.sim:",n0.sim,"    n1.sim:",n1.sim))
    if (type==1){
      if (n0.sim > 1 && n1.sim > 1){
        break;
      }    
    }
    if (type==2){
      if (n0.sim > 1 && n1.sim > 1 && abs(pct.0.sim-pct.0) <= 10){
        break;
      }    
    }
    sim.states = NULL
  }  
  return(sim.states)
}



make_points <- function(phy, states, sampling.f=c(1,1),name,n.points=100,
                        root.p=NULL,to_run="to_run.txt"){
  
  lik <- make.bisse(phy, states,sampling.f = sampling.f)  
  starting.p <- starting.point.bisse(phy)  
#   p <- starting.point.bisse(phy)
#   lik.0 <- make.bisse(phy, states,sampling.f = sampling.f)
#   fit.p6 <- find.mle(lik.0, p) 
#   print(fit.p6$par)
#   
#   lik=make.bd(phy)
#   fit <- find.mle(lik, c(.1, .03))
#   lik.mk2 <- make.mk2(phy, states, control=list(method="mk2"))
#   fit.mk2 <- find.mle(lik.mk2, c(.1, .1), method="subplex")
#   lik.mk1 <- constrain(lik.mk2, q10 ~ q01)
#   fit.mk1 <- find.mle(lik.mk1, 0.1)
  
#   my.starting.point = c(fit$par[1],fit$par[1],fit$par[2],fit$par[2],fit.mk1$par[1],fit.mk1$par[1])
#   print(my.starting.point)
  
  points=matrix(0,ncol=6,nrow=n.points)
  for (i in 1:6)
    points[,i]=sample(seq(0,10*max(starting.p),0.01), n.points)
  
  points[1,] = starting.p
  
  save(points,file=paste("points/",name,"_points.RData",sep=""))
  save(phy, states, sampling.f,root.p,file="data.RData")
  
  to_server_dir = paste("points/",name,"_points",sep="")
  dir.create(to_server_dir)
    
  for (i in 1:n.points){
    p=points[i,]
    
    file.create("p.RData")
    save(p,file="p.RData")
    
    p.name = paste("p",i,sep="")
    server_commands(p.name,to_run=to_run,to_server_dir,run_file="run_subplex.R",time_out=600) 
    file.copy("p.RData",paste(to_server_dir,"/",p.name,sep=""))  
    file.copy("data.RData",paste(to_server_dir,"/",p.name,sep=""))  
  }  
}

prepare_all_points <- function(trees_dir,genera,to_run = "to_run.txt",ANALYSIS=1){  
  
  load("species_per_genus.RData")
  for (t in 1:length(genera)){
    
    # file name
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    
    tree_file = paste(trees_dir,"/", genus, ".tre", sep="")  
    
    states = get_states(genus_file,ANALYSIS) 
    names(states)=gsub("-","_",names(states))
    species_num = species_per_genus$species_num[which(species_per_genus$Genus==genus)]
    ans = merge_tree_states(tree_file,states,species_num)
    
    if (!is.null(ans)){
      phy = ans[[1]]      
      states = ans[[2]]
      df = ans[[3]]
      sampling.f = ans[[4]] 
      
      if (df$n.sp.tree.inter >= 10 && df$n.0 >= 2 && df$n.1 >= 2){
        name = paste(genus,"_Analyais",ANALYSIS,sep="")
        print(genus)
        make_points(phy, states,sampling.f,name,to_run=to_run)        
        
      } 
    }
  }
}


plot_all_points <- function(genera,ANALYSIS=1){  
  
  for (t in 1:length(genera)){
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    
    name = paste("points/",genus,"_Analyais",ANALYSIS,"_points",sep="")
    
    if (file.exists(name)) read_points(name,plot=TRUE)    
  }
}

examine_all_points <- function(genera,ANALYSIS=1){  
  dd = data.frame(genus=genera,diff=NA)
  for (t in 1:length(genera)){
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    
    name = paste("points/",genus,"_Analyais",ANALYSIS,"_points",sep="")
    
    if (file.exists(name)){
      dd$diff[t] = read_points_diff(name)
    } 
  }
  print("likelihood difference between best points")
  print( dd[ order(-dd[,2]), ])
}

prepare_count <- function(trees_dir,genera,ANALYSIS=1){  
  stat = data.frame(genus=genera,n.sp=NA,n.sp.N=NA,n.sp.D=NA,n.sp.tree=NA,n.sp.tree.N=NA,n.sp.tree.D=NA)
  load("species_per_genus.RData")
  for (t in 1:length(genera)){
    
    # file name
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    
    tree_file = paste(trees_dir,"/", genus, ".tre", sep="")    
    phy_list <- read.tree(tree_file)
    
    states = get_states(genus_file,ANALYSIS) 
    names(states)=gsub("-","_",names(states))
    
    stat$n.sp[t] = length(states)
    stat$n.sp.N[t] = length(which(states==0))
    stat$n.sp.D[t] = length(which(states==1))
    species_num = species_per_genus$species_num[which(species_per_genus$Genus==genus)]
    ans = merge_tree_states(phy_list[[1]],states,species_num)
    
    if (!is.null(ans)){
      phy = ans[[1]]
      states = ans[[2]]
      df = ans[[3]]
      sampling.f = ans[[4]] 
      
      stat$n.sp.tree[t] = df$n.sp.tree.inter
      stat$n.sp.tree.N[t] = df$n.0
      stat$n.sp.tree.D[t] = df$n.1      
      
    }
  }
  return(stat)
}

prepare_to_server_sim <- function(trees_dir,genera,to_run = "to_run.txt",to_server_dir = "to_server",
                                  run_file = "run_on_server.R",NSTEPS = 2000,ANALYSIS=1,factor_vec=c(1,2,4),NSIM=100){  
  load("species_per_genus.RData")  
  for (t in 1:length(genera)){
    
    # file name
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    
    tree_file = paste(trees_dir,"/", genus, ".tre", sep="")  
    
    states = get_states(genus_file,ANALYSIS) 
    names(states)=gsub("-","_",names(states))
    species_num = species_per_genus$species_num[which(species_per_genus$Genus==genus)]
    ans = merge_tree_states(tree_file,states,species_num)
    
    if (!is.null(ans)){
      phy = ans[[1]]
      states_origin = ans[[2]]
      df = ans[[3]]
      sampling.f = ans[[4]] 
      
      if (df$n.sp.tree.inter >= 10 && df$n.0 >= 2 && df$n.1 >= 2){
        for (FACTOR in factor_vec){
          
          for (s in 1:NSIM){
            states = states_origin
            
            n = names(states)
            states = states[sample(length(states),length(states))]
            names(states) = n            
            
            name = paste(genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"_Sim",s,sep="")
            server_commands(name,to_run,to_server_dir,run_file,NSTEPS)    
            file.create("Data.RData")
            save(phy, states, name, FACTOR,df,sampling.f,file="Data.RData")
            file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
            file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))
          }         
          
        }
      } 
    }
  }
}

prepare_to_server_sim2 <- function(pars=c(0.2, 0.2, 0.03, 0.03, 0.03, 0.03),
                                   to_run = "to_run_sim2.txt",to_server_dir = "to_server_sim2",
                                   run_file = "run_on_server.R",NSTEPS = 2000,NSIM=100,max.t=30,FACTOR=2,min.taxa=100){  
  # here I simulate random trees with given pars    
  sampling.f = 1  
  
  for (s in 1:NSIM){
    phy=NULL
    while(is.null(phy)){
      phy <- tree.bisse(pars, max.t=max.t, x0=0)       
      if (!is.null(phy) && (length(phy$tip.state) < min.taxa || length(phy$tip.state) > 1.3*min.taxa || length(which(phy$tip.state==0))==length(phy$tip.state) || length(which(phy$tip.state==0))==0)){
        phy=NULL
      }
    }
    
    states <- phy$tip.state 
    name = paste("Sim2_",s,sep="")
    server_commands(name,to_run,to_server_dir,run_file,NSTEPS)    
    file.create("Data.RData")
    save(phy, states, name, FACTOR,sampling.f,file="Data.RData")
    file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
    file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))
  }  
}

prepare_to_server_sim2_factors <- function(pars=c(0.2, 0.2, 0.03, 0.03, 0.03, 0.03),
                                   to_run = "to_run_sim2.txt",to_server_dir = "to_server_sim2",
                                   run_file = "run_on_server.R",NSTEPS = 2000,NSIM=100,max.t=30,min.taxa=100){  
  # here I simulate random trees with given pars    
  sampling.f = 1  
  factor_vec=c(1,2,4)
  
  for (s in 1:NSIM){
    phy=NULL
    while(is.null(phy)){
      phy <- tree.bisse(pars, max.t=max.t, x0=0)       
      if (!is.null(phy) && (length(phy$tip.state) < min.taxa || length(phy$tip.state) > 1.3*min.taxa || length(which(phy$tip.state==0))==length(phy$tip.state) || length(which(phy$tip.state==0))==0)){
        phy=NULL
      }
    }
    for (FACTOR in factor_vec){
      states <- phy$tip.state 
      name = paste("Sim2_",s,"_Factor_",FACTOR,sep="")
      server_commands(name,to_run,to_server_dir,run_file,NSTEPS)    
      file.create("Data.RData")
      save(phy, states, name, FACTOR,sampling.f,file="Data.RData")
      file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
      file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))
    }
  }  
}

split_get <- function(a_list,split="_",get=1){
  g=strsplit(as.character(a_list),split=split,fixed=TRUE)
  x<-function(x) x[get]
  res=sapply(g, x, simplify = "vector")
  return(res)
}


prepare_to_server_sim3 <- function(phy=NULL,pars=c(0.2, 0.2, 0.03, 0.03, 0.03, 0.03),to_run = "to_run_sim3.txt",to_server_dir = "to_server_sim3",
                                   run_file = "run_on_server.R",NSTEPS = 2000,NSIM=100,min.taxa=100){  
  # here I simulate permutation trees with given pars    
  if (is.null(phy)){
    phy <- tree.bisse(pars, max.t=30, x0=0)
    col <- c("#004165", "#eaab00")
    plot(history.from.sim.discrete(phy, 0:1), phy, col=col)    
  }  
  states_origin <- phy$tip.state  
  sampling.f = 1  
  FACTOR=2
  
  for (s in 1:NSIM){
    states = states_origin
    
    n = names(states)
    states = states[sample(length(states),length(states))]
    names(states) = n            
    
    name = paste("Sim3_",s,sep="")
    server_commands(name,to_run,to_server_dir,run_file,NSTEPS)    
    file.create("Data.RData")
    save(phy, states, name, FACTOR,sampling.f,file="Data.RData")
    file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
    file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))
  }         
}




update_using_curated_data <- function(final_df1){
  f <-read.csv("Data/genus_count_manualy_currated.csv", header = TRUE,
               stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("NA",""),fill=TRUE,as.is=TRUE)
  for (i in 1:length(f$genus)){
    ind = which(final_df1$genus==f$genus[i])
    if (!is.na(f$perennial.1[i]) & f$perennial.1[i] == 1) final_df1$Life.Cycle[ind] = "perennial"
    if (!is.na(f$annual.1[i])    & f$annual.1[i] == 1)    final_df1$Life.Cycle[ind] = "annual"
    if (!is.na(f$woody.1[i])     & f$woody.1[i] == 1)     final_df1$Growth.Form[ind] = "woody"
    if (!is.na(f$herb.1[i])      & f$herb.1[i] == 1)      final_df1$Growth.Form[ind] = "herb"
  }
  
  return(final_df1)
}


update_using_all_db <- function(final_df1){
  f <-read.csv("Data/K_G_I_M_E_genusName.csv", header = TRUE,
               stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("NA",""),fill=TRUE,as.is=TRUE)
  df.sp = gsub("_"," ",final_df1$species)
  u_species = intersect(df.sp,f$Species)
  for (i in 1:length(u_species)){
    sp = u_species[i]
    ind1 = which(df.sp==sp)
    ind2 = which(f$Species==sp)
    if (!is.na(f$LIFECYCLE[ind2])  & is.na(final_df1$Life.Cycle[ind1]))  final_df1$Life.Cycle[ind1]  = f$LIFECYCLE[ind2]
    if (!is.na(f$GROWTHFORM[ind2]) & is.na(final_df1$Growth.Form[ind1])) final_df1$Growth.Form[ind1] = f$GROWTHFORM[ind2]
  }  
  return(final_df1)
}

get_ylim <- function(){
  usr <- par('usr')
  xr <- (usr[2] - usr[1]) / 27 # 27 = (100 + 2*4) / 4
  yr <- (usr[4] - usr[3]) / 27
  xlim <- c(usr[1] + xr, usr[2] - xr)
  ylim <- c(usr[3] + yr, usr[4] - yr)
  return(ylim)
}

get_xlim <- function(){
  usr <- par('usr')
  xr <- (usr[2] - usr[1]) / 27 # 27 = (100 + 2*4) / 4
  yr <- (usr[4] - usr[3]) / 27
  xlim <- c(usr[1] + xr, usr[2] - xr)
  ylim <- c(usr[3] + yr, usr[4] - yr)
  return(xlim)
}


plot_trop <- function(res,col1,col2,xlab,ylab){
  
  par(mar=mar_vec, xpd=TRUE)
  plot(res[,col1], res[,col2], 
       xlim = c(0,100),
       ylim = c(0,100),      
       xlab = xlab,
       ylab = ylab)
  c=cor.test(res[,col1], res[,col2])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""),cex=0.7)
  points(res[res$trop.dist=="trop",col1], res[res$trop.dist=="trop",col2], col = "red")
  points(res[res$trop.dist=="cosmo",col1], res[res$trop.dist=="cosmo",col2], col = "blue")
  points(res[res$trop.dist=="temp",col1], res[res$trop.dist=="temp",col2], col = "green")
  legend('topright',c("trop","cosmo","temp"), pch = 1, cex =legend_cex, inset=inset, col = c("red","blue","green"), title="Distribution")
  
}

plot_Fruit.dispersal <- function(res,col1,col2,xlab,ylab){
  
  par(mar=mar_vec, xpd=TRUE)
  plot(res[,col1], res[,col2], 
       xlim = c(0,100),
       ylim = c(0,100),      
       xlab = xlab,
       ylab = ylab)
  c=cor.test(res[,col1], res[,col2])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""),cex=0.7)
  points(res[res$Fruit.dispersal=="animals",col1], res[res$Fruit.dispersal=="animals",col2], col = "red")
  points(res[res$Fruit.dispersal=="abiotic",col1], res[res$Fruit.dispersal=="abiotic",col2], col = "blue")
  points(res[res$Fruit.dispersal=="variable",col1], res[res$Fruit.dispersal=="variable",col2], col = "green")
  legend('topright',c("animals","abiotic","variable"), pch = 1, cex =legend_cex, inset=inset, col = c("red","blue","green"), title="Dispersal")
  
}


plot_Pollination <- function(res,col1,col2,xlab,ylab){
  
  par(mar=mar_vec, xpd=TRUE)
  plot(res[,col1], res[,col2], 
       xlim = c(0,100),
       ylim = c(0,100),      
       xlab = xlab,
       ylab = ylab)
  c=cor.test(res[,col1], res[,col2])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""),cex=0.7)
  points(res[res$Pollination=="animals",col1], res[res$Pollination=="animals",col2], col = "red")
  points(res[res$Pollination=="abiotic",col1], res[res$Pollination=="abiotic",col2], col = "blue")
  points(res[res$Pollination=="variable",col1], res[res$Pollination=="variable",col2], col = "green")
  legend('topright',c("animals","abiotic","variable"), pch = 1, cex =legend_cex, inset=inset, col = c("red","blue","green"), title="Pollination")
  
}

plot_Growth.Form <- function(res,col1,col2,xlab,ylab){
  
  par(mar=mar_vec, xpd=TRUE)
  plot(res[,col1], res[,col2], 
       xlim = c(0,100),
       ylim = c(0,100),      
       xlab = xlab,
       ylab = ylab)
  c=cor.test(res[,col1], res[,col2])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""),cex=0.7)
  points(res[res$woody.c=="woody",col1], res[res$woody.c=="woody",col2], col = "red")
  points(res[res$woody.c=="herb",col1], res[res$woody.c=="herb",col2], col = "blue")
  legend('topright',c("woody","herb"), pch = 1, cex =legend_cex, inset=inset, col = c("red","blue"), title="Growth Form")
  
}



plot_monocot_dicot <- function(res,col1,col2,xlab,ylab){
  
  par(mar=mar_vec, xpd=TRUE)
  plot(res[,col1], res[,col2], 
       xlim = c(0,100),
       ylim = c(0,100),      
       xlab = xlab,
       ylab = ylab)
  c=cor.test(res[,col1], res[,col2])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
  points(res[res$monocots_dicots=="monocots",col1], res[res$monocots_dicots=="monocots",col2], col = "red")
  points(res[res$monocots_dicots=="dicots",col1], res[res$monocots_dicots=="dicots",col2], col = "blue")
  legend('topright',c("monocots","dicots"), pch = 1, cex =legend_cex, inset=inset, col = c("red","blue"), title="Group")
  
}

plot_H.pct <- function(res,col1,col2,xlab,ylab){
  
  par(mar=mar_vec, xpd=TRUE)
  plot(res[,col1], res[,col2], 
       xlim = c(0,100),
       ylim = c(0,100),      
       xlab = xlab,
       ylab = ylab)
  c=cor.test(res[,col1], res[,col2])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
  points(res[res$H.pct.c=="H",col1], res[res$H.pct.c=="H",col2], col = "red")
  points(res[res$H.pct.c=="M",col1], res[res$H.pct.c=="M",col2], col = "blue")
  legend('topright',c("Mostly H","Mostly M"), pch = 1, cex =legend_cex, inset=inset, col = c("red","blue"), title="Group")
  
}

plot_n <- function(res){
  par(mfrow=c(2,2), pty="s")
  
  plot(res[,"n.sp"], res[,"lam_f2"],main="Speciation", 
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(lam"[N] ~ "> lam"[D] ~ ")"))
  plot(res[,"n.sp"], res[,"mu_f2"], main="Extinction",
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(mu"[N] ~ "> mu"[D] ~ ")"))
  plot(res[,"n.sp"], res[,"div_f2"], main="Diversification",
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  plot(res[,"n.sp"], res[,"q_f2"], main="Transition",
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(q"[ND] ~ "> q"[DN] ~ ")"))
  
}

plot_col <- function(res,col,xlab){
  par(mfrow=c(2,2), pty="s")
  
  plot(res[,col], res[,"lam_f2"],main="Speciation", 
       ylim = c(0,100),
       xlab = xlab,
       ylab = bquote("PP(lam"[N] ~ "> lam"[D] ~ ")"))
  plot(res[,col], res[,"mu_f2"], main="Extinction",
       ylim = c(0,100),
       xlab = xlab,
       ylab = bquote("PP(mu"[N] ~ "> mu"[D] ~ ")"))
  plot(res[,col], res[,"div_f2"], main="Diversification",
       ylim = c(0,100),
       xlab = xlab,
       ylab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  plot(res[,col], res[,"q_f2"], main="Transition",
       ylim = c(0,100),
       xlab = xlab,
       ylab = bquote("PP(q"[ND] ~ "> q"[DN] ~ ")"))
  
}

plot_na <- function(res){
  par(mfrow=c(2,2), pty="s")
  
  plot(res[,"pct.non.na"], res[,"lam_f2"],main="Speciation", 
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(lam"[N] ~ "> lam"[D] ~ ")"))
  plot(res[,"pct.non.na"], res[,"mu_f2"], main="Extinction",
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(mu"[N] ~ "> mu"[D] ~ ")"))
  plot(res[,"pct.non.na"], res[,"div_f2"], main="Diversification",
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  plot(res[,"pct.non.na"], res[,"q_f2"], main="Transition",
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(q"[ND] ~ "> q"[DN] ~ ")"))
  
}

plot_box <- function(res,col="div_f2",param="r"){
  
#   res$trop.dist = factor(res$trop.dist)
#   res$Fruit.dispersal = factor(res$Fruit.dispersal)
#   res$Pollination = factor(res$Pollination)
#   res$woody.c = factor(res$woody.c)
     
  trop.dist = res[res$trop.dist!="cosmo",c(col,"trop.dist")]
  Fruit.dispersal = res[res$Fruit.dispersal!="variable",c(col,"Fruit.dispersal")]
  Pollination = res[res$Pollination!="variable",c(col,"Pollination")]
  
  ylab=paste("PP(",param,"N > ",param,"D)",sep="")
  
  par(mfrow=c(3,2), pty="m")
  
  boxplot(res[,col] ~ res$woody.c, main="Growth Form",
          names=c(paste("herb (",length(which(res$woody.c=="herb")),")",sep=""),paste("woody (",length(which(res$woody.c=="woody")),")",sep="")), 
          ylab=ylab) 
  z <- wilcox.test(res[which(res$woody.c=="woody"),col],res[which(res$woody.c=="herb"),col]) 
  mtext(paste("",round(z$p.value,4)),cex=0.8)
  
  boxplot(trop.dist[,col] ~ trop.dist$trop.dist, main="Geographical Distribution",
          names=c(paste("temp(",length(which(trop.dist$trop.dist=="temp")),")",sep=""),
                  paste("trop(",length(which(trop.dist$trop.dist=="trop")),")",sep="")), 
          ylab=ylab) 
  z <- lm(res$div_f2 ~ res$trop.dist);  az=anova(z); 
  tt <- wilcox.test(res[which(res$trop.dist=="trop"),col],res[which(res$trop.dist=="temp"),col]) 
  mtext(paste("",round(tt$p.value,4)),cex=0.8)
  #mtext(paste("anova:",round(az[1,"Pr(>F)"],4),"   t-test:",round(tt$p.value,4)),cex=0.8)
  
  boxplot(Fruit.dispersal[,col] ~ Fruit.dispersal$Fruit.dispersal, main="Fruit Dispersal",
          names=c(paste("animals(",length(which(Fruit.dispersal$Fruit.dispersal=="animals")),")",sep=""),
                  paste("abiotic(",length(which(Fruit.dispersal$Fruit.dispersal=="abiotic")),")",sep="")), 
          ylab=ylab) 
  z <- lm(res$div_f2 ~ res$Fruit.dispersal);  az=anova(z); #mtext(round(az[1,"Pr(>F)"],4),cex=0.8)
  tt <- wilcox.test(res[which(res$Fruit.dispersal=="animals"),col],res[which(res$Fruit.dispersal=="abiotic"),col]) 
  mtext(paste("",round(tt$p.value,4)),cex=0.8)
  #mtext(paste("anova:",round(az[1,"Pr(>F)"],4),"   t-test:",round(tt$p.value,4)),cex=0.8)
  
  boxplot(Pollination[,col] ~ Pollination$Pollination, main="Pollination" ,
          names=c(paste("animals(",length(which(Pollination$Pollination=="animals")),")",sep=""),
                  paste("abiotic(",length(which(Pollination$Pollination=="abiotic")),")",sep="")), 
          ylab=ylab) 
  z <- lm(res$div_f2 ~ res$Pollination);  az=anova(z); #mtext(round(az[1,"Pr(>F)"],4),cex=0.8)
  tt <- wilcox.test(res[which(res$Pollination=="animals"),col],res[which(res$Pollination=="abiotic"),col]) 
  mtext(paste("",round(tt$p.value,4)),cex=0.8)
  #mtext(paste("anova:",round(az[1,"Pr(>F)"],4),"   t-test:",round(tt$p.value,4)),cex=0.8)


  boxplot(res[,col] ~ res$H.pct.c, main="Hermaphrodite %", 
          names=c(paste("Mostly H(",length(which(res$H.pct.c=="H")),")",sep=""),
                  paste("Mostly M(",length(which(res$H.pct.c=="M")),")",sep="")),
          ylab=ylab) 
  z <- wilcox.test(res[which(res$H.pct.c=="H"),col],res[which(res$H.pct.c=="M"),col]) 
  mtext(round(z$p.value,4),cex=0.8)
  
}

plot_box2 <- function(res){
  
  res$trop.dist = factor(res$trop.dist)
  res$Fruit.dispersal = factor(res$Fruit.dispersal)
  res$Pollination = factor(res$Pollination)
  
  par(mfrow=c(2,3))
  boxplot(res[,"woody"] ~ res$trop.dist,       main="Geographical Distribution", ylab="Woody %") 
  z <- lm(res$woody ~ res$trop.dist);  az=anova(z); mtext(round(az[1,"Pr(>F)"],4),cex=0.8)
  
  boxplot(res[,"woody"] ~ res$Fruit.dispersal, main="Fruit Dispersal"          , ylab="Woody %") 
  z <- lm(res$woody ~ res$Fruit.dispersal);  az=anova(z); mtext(round(az[1,"Pr(>F)"],4),cex=0.8)
  
  boxplot(res[,"woody"] ~ res$Pollination,     main="Pollination"              , ylab="Woody %")
  z <- lm(res$woody ~ res$Pollination);  az=anova(z); mtext(round(az[1,"Pr(>F)"],4),cex=0.8)
  
  boxplot(res[,"perennial"] ~ res$trop.dist,       main="Geographical Distribution", ylab="Perennials %") 
  z <- lm(res$perennial ~ res$trop.dist);  az=anova(z); mtext(round(az[1,"Pr(>F)"],4),cex=0.8)
  
  boxplot(res[,"perennial"] ~ res$Fruit.dispersal, main="Fruit Dispersal"          , ylab="Perennials %") 
  z <- lm(res$perennial ~ res$Fruit.dispersal);  az=anova(z); mtext(round(az[1,"Pr(>F)"],4),cex=0.8)
  
  boxplot(res[,"perennial"] ~ res$Pollination,     main="Pollination"              , ylab="Perennials %")
  z <- lm(res$perennial ~ res$Pollination);  az=anova(z); mtext(round(az[1,"Pr(>F)"],4),cex=0.8)
  
}

plot_box_new <- function(res1,res2,col="div_f2"){
    
  trop.dist = res[res$trop.dist!="cosmo",c(col,"trop.dist")]
  Fruit.dispersal = res[res$Fruit.dispersal!="variable",c(col,"Fruit.dispersal")]
  Pollination = res[res$Pollination!="variable",c(col,"Pollination")]
  
  ylab=paste("PP(",col,"(N) > ",col,"(D))",sep="")
  
  par(mfrow=c(3,2), pty="m")
  
  boxplot(res[,col] ~ res$woody.c, main="Growth Form",
          names=c(paste("herb (",length(which(res$woody.c=="herb")),")",sep=""),paste("woody (",length(which(res$woody.c=="woody")),")",sep="")), 
          ylab=ylab) 
  z <- wilcox.test(res[which(res$woody.c=="woody"),col],res[which(res$woody.c=="herb"),col]) 
  mtext(paste("",round(z$p.value,4)),cex=0.8)
  
  boxplot(trop.dist[,col] ~ trop.dist$trop.dist, main="Geographical Distribution",
          names=c(paste("temp(",length(which(trop.dist$trop.dist=="temp")),")",sep=""),
                  paste("trop(",length(which(trop.dist$trop.dist=="trop")),")",sep="")), 
          ylab=ylab) 
  z <- lm(res$div_f2 ~ res$trop.dist);  az=anova(z); 
  tt <- wilcox.test(res[which(res$trop.dist=="trop"),col],res[which(res$trop.dist=="temp"),col]) 
  mtext(paste("",round(tt$p.value,4)),cex=0.8)
  #mtext(paste("anova:",round(az[1,"Pr(>F)"],4),"   t-test:",round(tt$p.value,4)),cex=0.8)
  
  boxplot(Fruit.dispersal[,col] ~ Fruit.dispersal$Fruit.dispersal, main="Fruit Dispersal",
          names=c(paste("animals(",length(which(Fruit.dispersal$Fruit.dispersal=="animals")),")",sep=""),
                  paste("abiotic(",length(which(Fruit.dispersal$Fruit.dispersal=="abiotic")),")",sep="")), 
          ylab=ylab) 
  z <- lm(res$div_f2 ~ res$Fruit.dispersal);  az=anova(z); #mtext(round(az[1,"Pr(>F)"],4),cex=0.8)
  tt <- wilcox.test(res[which(res$Fruit.dispersal=="animals"),col],res[which(res$Fruit.dispersal=="abiotic"),col]) 
  mtext(paste("",round(tt$p.value,4)),cex=0.8)
  #mtext(paste("anova:",round(az[1,"Pr(>F)"],4),"   t-test:",round(tt$p.value,4)),cex=0.8)
  
  boxplot(Pollination[,col] ~ Pollination$Pollination, main="Pollination" ,
          names=c(paste("animals(",length(which(Pollination$Pollination=="animals")),")",sep=""),
                  paste("abiotic(",length(which(Pollination$Pollination=="abiotic")),")",sep="")), 
          ylab=ylab) 
  z <- lm(res$div_f2 ~ res$Pollination);  az=anova(z); #mtext(round(az[1,"Pr(>F)"],4),cex=0.8)
  tt <- wilcox.test(res[which(res$Pollination=="animals"),col],res[which(res$Pollination=="abiotic"),col]) 
  mtext(paste("",round(tt$p.value,4)),cex=0.8)
  #mtext(paste("anova:",round(az[1,"Pr(>F)"],4),"   t-test:",round(tt$p.value,4)),cex=0.8)
  
  
  boxplot(res[,col] ~ res$H.pct.c, main="Hermaphrodite %", 
          names=c(paste("Mostly H(",length(which(res$H.pct.c=="H")),")",sep=""),
                  paste("Mostly M(",length(which(res$H.pct.c=="M")),")",sep="")),
          ylab=ylab) 
  z <- wilcox.test(res[which(res$H.pct.c=="H"),col],res[which(res$H.pct.c=="M"),col]) 
  mtext(round(z$p.value,4),cex=0.8)
  
}

# plot_box3 <- function(res){
#   
#   trop.dist = res[res$trop.dist!="cosmo",c("div_f2","trop.dist")]
#   Fruit.dispersal = res[res$Fruit.dispersal!="variable",c("div_f2","Fruit.dispersal")]
#   Pollination = res[res$Pollination!="variable",c("div_f2","Pollination")]
#   
#   res$Cat = NA
#   res$Cat[res$pct.0>=70] = "H/M"
#   res$Cat[res$pct.0<=30] = "Dioecy"
#   
#   par(mfrow=c(2,2))
#   
#   boxplot(res[,"woody"] ~ res$Cat, main="Growth Form",          
#           ylab=bquote("PP(r"[N] ~ "> r"[D] ~ ")")) 
#   z <- wilcox.test(res$div_f2[which(res$woody.c=="woody")],res$div_f2[which(res$woody.c=="herb")]) 
#   mtext(paste("",round(z$p.value,4)),cex=0.8)
#   
# }


plot_from_server <- function(trees_dir,genera,to_server_dir = "to_server",
                             ANALYSIS=1,factor_vec=c(1,2,4),to_pdf=TRUE){ 
  pdf(file=paste("res_a",ANALYSIS,to_server_dir,".pdf",sep=""))
  
  dat = matrix(NA,ncol=8,nrow=length(genera))
  res = matrix(NA,ncol=4*length(factor_vec),nrow=length(genera))
  woody_vec = matrix(NA,ncol=1,nrow=length(genera))
  load("properties.sp")
  load("properties.g")
  for (t in 1:length(genera)){
    
    # file name
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    
    ml_file = paste("from_server/",to_server_dir,"/",genus,"_Analyais",ANALYSIS,"_Factor",1,"/",
                    genus,"_Analyais",ANALYSIS,"_Factor",1,".ML",sep="")
    df.file = paste("from_server/",to_server_dir,"/",genus,"_Analyais",ANALYSIS,"_Factor",1,"/Data.RData",sep="")
    
    ind=1
    for (FACTOR in factor_vec){
      
      out_file = paste("from_server/",to_server_dir,"/",genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"/",
                       genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,".mcmc",sep="")      
      if (file.exists(out_file)){
        load(out_file)       
        lam = 100*length(which(samples$lambda0>samples$lambda1))/length(samples$lambda0)
        mu = 100*length(which(samples$mu0>samples$mu1))/length(samples$lambda0)
        div = 100*length(which(samples$div0>samples$div1))/length(samples$lambda0)
        q = 100*length(which(samples$q01>samples$q10))/length(samples$lambda0)
        
        res[t,seq(ind,ind+3)] = round(c(lam,mu,div,q))
        ind = ind+4
        if (FACTOR==2){
          
          load(df.file)         
          dat[t,] = c(df$n.sp.tree, df$n.sp.tree.inter, df$pct.non.na, 
                      df$pct.0, df$pct.1, df$n.0, df$n.1,df$sampling.f)    
          
          ###########   plotting #############################
          tt = "phylogram" 
          if (length(phy$tip.label)>120)
            tt="fan"
          statecols <- c("red","blue") 
          plot(phy, label.offset=0.02,tip.color=statecols[states+1], 
               cex=(0.2+(10/length(phy$tip.label))), main=genus,type=tt ) 
          
          tiplabels(col=statecols[states+1], pch=19, cex=(0.1+(1/length(phy$tip.label))))          
          
          legend("topleft", c("H/M","Dioecy"),col=statecols,pch=19)      
          
          plot_mcmc_results_new(samples,"")
        }
      }
    }    
  }  
  dev.off()
  
  if (to_pdf)
    pdf(file=paste("sum_a",ANALYSIS,to_server_dir,".pdf",sep=""))
  
  temp = round(cbind(dat,res),1)
  colnames(temp)=c("n.sp.tree", "n.sp.var.subsp.rm", "pct.non.na", "pct.0", 
                   "pct.1", "n.0", "n.1","sampling.f",
                   "lam_f1","mu_f1","div_f1","q_f1",
                   "lam_f2","mu_f2","div_f2","q_f2",
                   "lam_f4","mu_f4","div_f4","q_f4")
  
  all.res = data.frame(genus=gsub(".csv","",genera),temp)
   
  all.res$diff.pct = abs(all.res$div_f4-all.res$div_f1)
  
  save(all.res,file=paste("all.res_a",ANALYSIS,to_server_dir,".RData",sep=""))
  
  analyzed = which(all.res$n.sp.var.subsp.rm>0)
  #selected = which(all.res$flag==0 & all.res$flag2==0 & all.res$diff.pct<=20 & all.res$n.sp.var.subsp.rm>=10 & all.res$n.0>=2 & all.res$n.1>=2 & (all.res$n.0+all.res$n.1>=10))
  selected0 = which(all.res$n.sp.var.subsp.rm>=10 & all.res$n.0>=2 & all.res$n.1>=2 & (all.res$n.0+all.res$n.1>=10))
  selected = which(all.res$diff.pct<=20 & all.res$n.sp.var.subsp.rm>=10 & all.res$n.0>=2 & all.res$n.1>=2 & (all.res$n.0+all.res$n.1>=10))
  
  excluded = setdiff(analyzed,selected)
  
  print(paste("all:",length(analyzed),
              "  genera with n.sp >= 10 and n.sp(each state) >= 2:",length(selected0),
              "  afetr Factors filter:",length(selected)))
   
  old.par <- par( no.readonly = TRUE ) 
  
  par(mfrow=c(1,2), pty="s")
  
  h=hist(all.res[selected,"div_f2"],seq(0,100,10),
         main=paste("selected genera (",length(selected),")",sep=""),
         xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  tt=wilcox.test(all.res[selected,"div_f2"],mu=50)
  #mtext(paste("wilcox.test: p =",round(tt$p.value,4)))
  
  print(dip.test(all.res[selected,"div_f2"]))
  #print(dip.test(all.res[selected,"div_f2"], simulate=TRUE, B=5000))
  
  ylim=get_ylim()
  
  rr = par('usr')
  hist(all.res[excluded,"div_f2"],seq(0,100,10),
       main=paste("excluded genera (",length(excluded),")",sep=""),
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")")
       ,ylim=ylim)
  
  par(mfrow=c(1,2), pty="s")
  
  plot(all.res[selected,"pct.0"], all.res[selected,"q_f2"], 
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "N %",
       ylab = "q(N->D) > q(D->N) % ")
  #,main="Srict Dioecy")
  text(all.res[selected,"pct.0"], all.res[selected,"q_f2"]+2, all.res[selected,"genera"],cex=0.7)
  c=cor.test(all.res[selected,"pct.0"], all.res[selected,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
  
  plot(all.res[selected,"div_f2"], all.res[selected,"q_f2"], 
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "r(N) > r(D) %",
       ylab = "q(N->D) > q(D->N) %")
  #,main="Srict Dioecy")
  text(all.res[selected,"div_f2"], all.res[selected,"q_f2"]+2, all.res[selected,"genera"],cex=0.7)
  c=cor.test(all.res[selected,"div_f2"], all.res[selected,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
  
  
  res = merge(all.res[selected,],properties.g,by="genus")
  
  par(mfrow=c(2,2), pty="s")
  
  plot(res[,"n.sp"], res[,"lam_f2"], 
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(lam"[N] ~ "> lam"[D] ~ ")"))
  plot(res[,"n.sp"], res[,"mu_f2"], 
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(mu"[N] ~ "> mu"[D] ~ ")"))
  plot(res[,"n.sp"], res[,"div_f2"], 
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  plot(res[,"n.sp"], res[,"q_f2"], 
       ylim = c(0,100),
       xlab = "# species",
       ylab = bquote("PP(q"[ND] ~ "> q"[DN] ~ ")"))
  
   
  
  # make woody categorical
  res$woody.c = NA
  res$woody.c[res$woody>=70] = "woody"
  res$woody.c[res$woody<=30] = "herb"
  
  # make perennial categorical
  res$perennial.c = NA
  res$perennial.c[res$perennial>=70] = "perennial"
  res$perennial.c[res$perennial<=30] = "annual"
  
  res$H.pct.c = NA
  if (ANALYSIS==1){
    res$H.pct.c[res$H.pct.1>=70] = "H"
    res$H.pct.c[res$H.pct.1<=30] = "M"
  }
  if (ANALYSIS==2){
    res$H.pct.c[res$H.pct.2>=70] = "H"
    res$H.pct.c[res$H.pct.2<=30] = "M"
  }
  
  save(res,file=paste("res_a",ANALYSIS,to_server_dir,".RData",sep=""))
  
  plot_box(res)
  
  #plot_box2(res)
  
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","r(N) > r(D) % ")
  plot_Fruit.dispersal(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","r(N) > r(D) % ")
  plot_Pollination(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","r(N) > r(D) % ")
  plot_monocot_dicot(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","r(N) > r(D) % ")
  
  
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,"pct.0","div_f2","N %","r(N) > r(D) % ")
  plot_Fruit.dispersal(res,"pct.0","div_f2","N %","r(N) > r(D) % ")
  plot_Pollination(res,"pct.0","div_f2","N %","r(N) > r(D) % ")
  plot_H.pct(res,"pct.0","div_f2","N %","r(N) > r(D) % ")
 
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,"pct.0","lam_f2","N %","lam(N) > lam(D) % ")
  plot_Fruit.dispersal(res,"pct.0","lam_f2","N %","lam(N) > lam(D) % ")
  plot_Pollination(res,"pct.0","lam_f2","N %","lam(N) > lam(D) % ")
  plot_H.pct(res,"pct.0","lam_f2","N %","lam(N) > lam(D) % ")
  
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,"pct.0","mu_f2","N %","mu(N) > mu(D) % ")
  plot_Fruit.dispersal(res,"pct.0","mu_f2","N %","mu(N) > mu(D) % ")
  plot_Pollination(res,"pct.0","mu_f2","N %","mu(N) > mu(D) % ")
  plot_H.pct(res,"pct.0","mu_f2","N %","mu(N) > mu(D) % ")
  
  print(paste("number of genera with enough data:",length(which(!is.na(all.res$div_f2)))))
  print(paste("number of genera filtered out:",length(excluded)))
  print(paste("number of genera remaining:",length(selected)))
  
  print(paste("mean number of taxa per genus:", mean(res$n.sp.var.subsp.rm)))
  
  c=cor.test(res$lam_f2,res$div_f2)
  print(paste("corr(lam,div):   r  =  ",round(c$estimate,2),"          p  =  ",c$p.value,sep=""))
  c=cor.test(res$mu_f2,res$div_f2)
  print(paste("corr(mu,div):   r  =  ",round(c$estimate,2),"          p  =  ",c$p.value,sep=""))
  
  
  tt = wilcox.test(res$div_f2,mu=50)
  print(paste("wilcox.test(res$div_f2,mu=50):  p = ",tt$p.value,sep=""))
  
  
  #hist(res$n.sp.var.subsp.rm)
  
  par( old.par )
  
  if (to_pdf)
    dev.off()
}

plot_from_server_sim <- function(trees_dir,genera,to_server_dir = "to_server_sim", ANALYSIS=1,FACTOR=2,NSIM=100){ 
  
  pdf(file=paste("res_a",ANALYSIS,to_server_dir,".pdf",sep=""))
  
  res_div = matrix(NA,ncol=NSIM,nrow=length(genera))
  res_lam = matrix(NA,ncol=NSIM,nrow=length(genera))
  res_mu = matrix(NA,ncol=NSIM,nrow=length(genera))
  res_q = matrix(NA,ncol=NSIM,nrow=length(genera))
  load(paste("all.res_a",ANALYSIS,"to_server.RData",sep=""))
  
  for (t in 1:length(genera)){
    
    genus = genera[t]   
    
    for (s in 1:NSIM){
      
      out_file = paste("from_server/",to_server_dir,"/",genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"_Sim",s,"/",
                       genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"_Sim",s,".mcmc",sep="")      
      if (file.exists(out_file)){
        load(out_file)       
        lam = 100*length(which(samples$lambda0>samples$lambda1))/length(samples$lambda0)
        mu = 100*length(which(samples$mu0>samples$mu1))/length(samples$lambda0)
        div = 100*length(which(samples$div0>samples$div1))/length(samples$lambda0)
        q = 100*length(which(samples$q01>samples$q10))/length(samples$lambda0)
        
        res_div[t,s] = round(div)
        res_lam[t,s] = round(lam)
        res_mu[t,s] = round(mu)
        res_q[t,s] = round(q)
        
        #         remove(states)        
        #         
        #         df.file = paste("from_server/",to_server_dir,"/",genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"_Sim",s,"/Data.RData",sep="")
        #         load(df.file)
        #         
        #         statecols <- c("red","blue")    
        #         plot(phy, label.offset=0.008, cex=0.6,tip.color=statecols[states+1],main=paste(s,res_div[t,s])) 
        #         tiplabels(col=statecols[states+1], pch=19,cex=0.6)
        #         legend("topleft", c("H/M","Dioecy"),col=statecols,pch=19)
        
      }
    }
    
    par(mfrow=c(2,2))
    
    h=hist(res_div[t,],seq(0,100,10),
           main=paste("diversification - ",genus,sep=""),
           xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
    mtext(paste("non-na%:",all.res[which(all.res$genus==genus),"pct.non.na"],"      N%:",all.res[which(all.res$genus==genus),"pct.0"]))
    ylim=get_ylim()
    text(all.res[which(all.res$genus==genus),"div_f2"]+2,ylim[2],
         length(which(sort(res_div[t,])>all.res[which(all.res$genus==genus),"div_f2"])))
    lines(x=c(all.res[which(all.res$genus==genus),"div_f2"],all.res[which(all.res$genus==genus),"div_f2"]),
          y=c(0,100), col="red")
    lines(x=c(all.res[which(all.res$genus==genus),"div_f1"],all.res[which(all.res$genus==genus),"div_f1"]),
          y=c(0,100), col="blue")    
    lines(x=c(all.res[which(all.res$genus==genus),"div_f4"],all.res[which(all.res$genus==genus),"div_f4"]),
          y=c(0,100), col="green")
    
    h=hist(res_lam[t,],seq(0,100,10),
           main=paste("speciation - ",genus,sep=""),
           xlab = bquote("PP(lam"[N] ~ "> lam"[D] ~ ")"))
    ylim=get_ylim()
    text(all.res[which(all.res$genus==genus),"lam_f2"]+2,ylim[2],
         length(which(sort(res_lam[t,])>all.res[which(all.res$genus==genus),"lam_f2"])))
    lines(x=c(all.res[which(all.res$genus==genus),"lam_f2"],all.res[which(all.res$genus==genus),"lam_f2"]),
          y=c(0,100), col="red")
    lines(x=c(all.res[which(all.res$genus==genus),"lam_f1"],all.res[which(all.res$genus==genus),"lam_f1"]),
          y=c(0,100), col="blue")    
    lines(x=c(all.res[which(all.res$genus==genus),"lam_f4"],all.res[which(all.res$genus==genus),"lam_f4"]),
          y=c(0,100), col="green")
    
    h=hist(res_mu[t,],seq(0,100,10),
           main=paste("extinction - ",genus,sep=""),
           xlab = bquote("PP(mu"[N] ~ "> mu"[D] ~ ")"))
    ylim=get_ylim()
    text(all.res[which(all.res$genus==genus),"mu_f2"]+2,ylim[2],
         length(which(sort(res_mu[t,])>all.res[which(all.res$genus==genus),"mu_f2"])))
    lines(x=c(all.res[which(all.res$genus==genus),"mu_f2"],all.res[which(all.res$genus==genus),"mu_f2"]),
          y=c(0,100), col="red")
    lines(x=c(all.res[which(all.res$genus==genus),"mu_f1"],all.res[which(all.res$genus==genus),"mu_f1"]),
          y=c(0,100), col="blue")    
    lines(x=c(all.res[which(all.res$genus==genus),"mu_f4"],all.res[which(all.res$genus==genus),"mu_f4"]),
          y=c(0,100), col="green")
    
    h=hist(res_q[t,],seq(0,100,10),
           main=paste("transition - ",genus,sep=""),
           xlab = bquote("PP(q"[ND] ~ "> q"[DN] ~ ")"))
    ylim=get_ylim()
    text(all.res[which(all.res$genus==genus),"q_f2"]+2,ylim[2],
         length(which(sort(res_q[t,])>all.res[which(all.res$genus==genus),"q_f2"])))
    lines(x=c(all.res[which(all.res$genus==genus),"q_f2"],all.res[which(all.res$genus==genus),"q_f2"]),
          y=c(0,100), col="red")
    lines(x=c(all.res[which(all.res$genus==genus),"q_f1"],all.res[which(all.res$genus==genus),"q_f1"]),
          y=c(0,100), col="blue")    
    lines(x=c(all.res[which(all.res$genus==genus),"q_f4"],all.res[which(all.res$genus==genus),"q_f4"]),
          y=c(0,100), col="green")
    
  }
  
  dev.off()
}

plot_from_server_sim3 <- function(to_server_dir = "to_server_sim3", FACTOR=2,NSIM=100){ 
  
  pdf(file=paste("res_",to_server_dir,".pdf",sep=""))
  
  res_div = matrix(NA,ncol=NSIM,nrow=1)
  res_lam = matrix(NA,ncol=NSIM,nrow=1)
  res_mu = matrix(NA,ncol=NSIM,nrow=1)
  res_q = matrix(NA,ncol=NSIM,nrow=1)
  
  pct = matrix(NA,ncol=NSIM,nrow=1)
  for (s in 1:NSIM){
    
    out_file = paste("from_server/",to_server_dir,"/","Sim3_",s,"/",
                     "Sim3_",s,".mcmc",sep="")     
    data_file = paste("from_server/",to_server_dir,"/","Sim3_",s,"/",
                      "Data.RData",sep="")     
    if (file.exists(out_file)){
      load(out_file)       
      lam = 100*length(which(samples$lambda0>samples$lambda1))/length(samples$lambda0)
      mu = 100*length(which(samples$mu0>samples$mu1))/length(samples$lambda0)
      div = 100*length(which(samples$div0>samples$div1))/length(samples$lambda0)
      q = 100*length(which(samples$q01>samples$q10))/length(samples$lambda0)
      
      res_div[s] = round(div)
      res_lam[s] = round(lam)
      res_mu[s] = round(mu)
      res_q[s] = round(q)
      
      load(data_file)
      pct[s]=length(which(states==0))
    }else{
      print(paste("no out file:",out_file))
    }  
  }
  
  res_div = res_div[!is.na(pct)]
  res_lam = res_lam[!is.na(pct)]
  res_mu = res_mu[!is.na(pct)]
  res_q = res_q[!is.na(pct)]
  pct = pct[!is.na(pct)]
  
  par(mfrow=c(2,2))
  
  h=hist(res_div,seq(0,100,10),
         main=paste("diversification - ","Permutation of states",sep=""),
         xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))   
  lines(x=c(median(res_div),median(res_div)), y=c(0,100), col="blue")
  
  h=hist(res_lam,seq(0,100,10),
         main=paste("speciation - ","Permutation of states",sep=""),
         xlab = bquote("PP(lam"[N] ~ "> lam"[D] ~ ")")) 
  lines(x=c(median(res_lam),median(res_lam)), y=c(0,100), col="blue")
  
  h=hist(res_mu,seq(0,100,10),
         main=paste("extinction - ","Permutation of states",sep=""),
         xlab = bquote("PP(mu"[N] ~ "> mu"[D] ~ ")"))
  lines(x=c(median(res_mu),median(res_mu)), y=c(0,100), col="blue")
  
  h=hist(res_q,seq(0,100,10),
         main=paste("transition - ","Permutation of states",sep=""),
         xlab = bquote("PP(q"[ND] ~ "> q"[DN] ~ ")"))   
  lines(x=c(median(res_q),median(res_q)), y=c(0,100), col="blue")
  dev.off()
}



plot_from_server_sim2 <- function(to_server_dir = "to_server_sim2", FACTOR=2,NSIM=100,to_pdf=TRUE){ 
  
  if (to_pdf)
    pdf(file=paste("res_",to_server_dir,".pdf",sep=""))
  
  res_div = matrix(NA,ncol=NSIM,nrow=1)
  res_lam = matrix(NA,ncol=NSIM,nrow=1)
  res_mu = matrix(NA,ncol=NSIM,nrow=1)
  res_q = matrix(NA,ncol=NSIM,nrow=1)
  
  pct = matrix(NA,ncol=NSIM,nrow=1)
  changes = matrix(NA,ncol=NSIM,nrow=1)
  for (s in 1:NSIM){
    
    out_file = paste("from_server/",to_server_dir,"/","Sim2_",s,"/",
                     "Sim2_",s,".mcmc",sep="")     
    data_file = paste("from_server/",to_server_dir,"/","Sim2_",s,"/",
                      "Data.RData",sep="")     
    if (file.exists(out_file)){
      load(out_file)       
      lam = 100*length(which(samples$lambda0>samples$lambda1))/length(samples$lambda0)
      mu = 100*length(which(samples$mu0>samples$mu1))/length(samples$lambda0)
      div = 100*length(which(samples$div0>samples$div1))/length(samples$lambda0)
      q = 100*length(which(samples$q01>samples$q10))/length(samples$lambda0)
      
      res_div[s] = round(div)
      res_lam[s] = round(lam)
      res_mu[s] = round(mu)
      res_q[s] = round(q)
      
      load(data_file)
      pct[s]=100*length(which(states==0))/length(states)      
      changes[s]=count_changes(phy)
    }else{
      print(paste("no out file:",out_file))
    }  
  }
  
  res_div = res_div[!is.na(pct)]
  res_lam = res_lam[!is.na(pct)]
  res_mu = res_mu[!is.na(pct)]
  res_q = res_q[!is.na(pct)]
  changes = changes[!is.na(pct)]
  pct = pct[!is.na(pct)]
  
  par(mfrow=c(2,2))
  txt = gsub("to_server_sim2_","N=",to_server_dir)
  h=hist(res_div,seq(0,100,10),
         main=paste("diversification  ",txt,sep=""),
         xlab = bquote("PP(r"[0] ~ "> r"[1] ~ ")"))   
  lines(x=c(median(res_div),median(res_div)), y=c(0,100), col="blue")
  plot_mtext(res_div)
  
  h=hist(res_lam,seq(0,100,10),
         main=paste("speciation  ",txt,sep=""),
         xlab = bquote("PP(lam"[0] ~ "> lam"[1] ~ ")")) 
  lines(x=c(median(res_lam),median(res_lam)), y=c(0,100), col="blue")
  plot_mtext(res_lam)
  
  h=hist(res_mu,seq(0,100,10),
         main=paste("extinction  ",txt,sep=""),
         xlab = bquote("PP(mu"[0] ~ "> mu"[1] ~ ")"))
  lines(x=c(median(res_mu),median(res_mu)), y=c(0,100), col="blue")
  plot_mtext(res_mu)
  
  h=hist(res_q,seq(0,100,10),
         main=paste("transition  ",txt,sep=""),
         xlab = bquote("PP(q"[0] ~ "> q"[1] ~ ")"))   
  lines(x=c(median(res_q),median(res_q)), y=c(0,100), col="blue")
  plot_mtext(res_q)
  
  par(mfrow=c(2,2))
  
  plot_r(pct,res_div,xlab="% state 0",ylab="PP(r0 > r1)",title="diversification")
  plot_r(pct,res_lam,xlab="% state 0",ylab="PP(lam0 > lam1)",title="speciation")
  plot_r(pct,res_mu,xlab="% state 0",ylab="PP(mu0 > mu1)",title="extinction")
  plot_r(pct,res_q,xlab="% state 0",ylab="PP(q01 > q10)",title="transition")

  par(mfrow=c(2,2))
  
  plot_r(changes,res_div,xlab="number of state transitions",ylab="PP(r0 > r1)",title="diversification")
  plot_r(changes,res_lam,xlab="number of state transitions",ylab="PP(lam0 > lam1)",title="speciation")
  plot_r(changes,res_mu,xlab="number of state transitions",ylab="PP(mu0 > mu1)",title="extinction")
  plot_r(changes,res_q,xlab="number of state transitions",ylab="PP(q01 > q10)",title="transition")
  
  if (to_pdf)
    dev.off()
}

plot_r <- function(x,y,xlab,ylab,title){
  plot(x,y,xlab=xlab,ylab = ylab, main=title)  
  abline(lm(as.vector(y)~as.vector(x)), col="red")
  c=cor.test(as.vector(y),as.vector(x))
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))   
}

plot_mtext <- function(x){
  if (length(unique(x))>1){
    z <- wilcox.test(x,mu=50) 
    x[which(x==100)]=99.9
    x[which(x==0)]=0.001
    z1 <- wilcox.test(probit(x/100),mu=0) 
    mtext(paste("wilcox.test:",round(z$p.value,4),"   t-test(probit):",round(z1$p.value,4)),cex=0.7)      
  }
}

plot_mtext2 <- function(x){
  if (length(unique(x))>1){
    z <- wilcox.test(x,mu=0)   
    mtext(paste("t-test:",round(z$p.value,4)),cex=0.7)      
  }
}

count_changes <- function(phy){
  count=0
  for (i in 2:length(phy$node.state)){
    child_node = i
    edge_index = which(phy$edge[,2]==child_node+length(phy$tip.label))
    parent_node = phy$edge[edge_index,1]-length(phy$tip.label)
    if (phy$node.state[child_node] != phy$node.state[parent_node]){
      count = count+1
      #print(paste(phy$node.label[parent_node],phy$node.label[child_node]))
    }
  }
  for (i in 1:length(phy$tip.state)){
    child_node = i
    edge_index = which(phy$edge[,2]==child_node)
    parent_node = phy$edge[edge_index,1]-length(phy$tip.label)
    if (phy$tip.state[child_node] != phy$node.state[parent_node]){
      count = count+1
      #print(paste(phy$node.label[parent_node],phy$tip.label[child_node]))
    }
  }
  return(count)
}


plot_from_server_ml <- function(to_server_dir, FACTOR=2,NSIM=100,to_pdf=TRUE){ 
  
  if (to_pdf)
    pdf(file=paste("res_",to_server_dir,".pdf",sep=""))
  
  res_div = matrix(NA,ncol=NSIM,nrow=1)
  res_lam = matrix(NA,ncol=NSIM,nrow=1)
  res_mu = matrix(NA,ncol=NSIM,nrow=1)
  res_q = matrix(NA,ncol=NSIM,nrow=1)
  
  pct = matrix(NA,ncol=NSIM,nrow=1)
  changes = matrix(NA,ncol=NSIM,nrow=1)
  
  ml = matrix(NA,nrow=NSIM,ncol=6)
  ml_p = matrix(NA,nrow=NSIM,ncol=2)
  for (s in 1:NSIM){
    
    out_file = paste("from_server/",to_server_dir,"/","Sim2_",s,"/",
                     "Sim2_",s,".mcmc",sep="")     
    data_file = paste("from_server/",to_server_dir,"/","Sim2_",s,"/",
                      "Data.RData",sep="")  
    ml_file = paste("from_server/",to_server_dir,"/","Sim2_",s,"/",
                    "Sim2_",s,".ML",sep="")
    if (file.exists(out_file)){
      load(out_file)       
      lam = 100*length(which(samples$lambda0>samples$lambda1))/length(samples$lambda0)
      mu = 100*length(which(samples$mu0>samples$mu1))/length(samples$lambda0)
      div = 100*length(which(samples$div0>samples$div1))/length(samples$lambda0)
      q = 100*length(which(samples$q01>samples$q10))/length(samples$lambda0)
      
      res_div[s] = round(div)
      res_lam[s] = round(lam)
      res_mu[s] = round(mu)
      res_q[s] = round(q)
      
      load(data_file)
      pct[s]=100*length(which(states==0))/length(states)      
      changes[s]=count_changes(phy)
      
      load(ml_file)
      ml[s,1:length(fit$par)] = fit$par 
      ml_p[s,1] = anova(fit, fit.l1)[2,5]
      ml_p[s,2] = anova(fit, fit.l2)[2,5]
      
    }else{
      print(paste("no out file:",out_file))
    }  
  }
  
  res_div = res_div[!is.na(pct)]
  res_lam = res_lam[!is.na(pct)]
  res_mu = res_mu[!is.na(pct)]
  res_q = res_q[!is.na(pct)]
  changes = changes[!is.na(pct)]
  pct = pct[!is.na(pct)]
  
  #####  ML  #################################################
  
  par(mfrow=c(2,2))
  
  x=log(ml[,1]/ml[,2])
  h=hist(x,seq(-2,20,0.2),xlim=c(-2,2),
         main="speciation",
         xlab = bquote("log(lam"[0] ~ "/ lam"[1] ~ ")"))   
  lines(x=c(median(x),median(x)), y=c(0,100), col="blue")
  plot_mtext2(x)
  
  x=log(ml[,3]/ml[,4])
  h=hist(x,seq(-20,200,0.2),xlim=c(-2,2),
         main="extinction",
         xlab = bquote("log(mu"[0] ~ "/ mu"[1] ~ ")"))   
  lines(x=c(median(x),median(x)), y=c(0,100), col="blue")
  plot_mtext2(x)
  
  x=log(ml[,5]/ml[,6])
  h=hist(x,seq(-20,200,0.2),xlim=c(-2,2),
         main="transition",
         xlab = bquote("log(q"[01] ~ "/ q"[10] ~ ")"))   
  lines(x=c(median(x),median(x)), y=c(0,100), col="blue")
  plot_mtext2(x)
  
  barplot((100/length(ml_p[,1]))*c(length(which(ml_p[,1]<0.05)),
                                   length(which(ml_p[,2]<0.05))),
              names.arg=c("lambda","mu"),ylim=c(0,100),ylab="p < 0.05 (%)",
          main="significant acceptance of constrained model")
  lines(x=c(0,3), y=c(5,5), col="blue")
  
  par(mfrow=c(2,2))
  
  plot_r(pct,log((ml[,1]-ml[,3])/(ml[,2]-ml[,4])),xlab="% state 0",ylab="log(lam0-mu0/lam1-mu1)",title="diversification")
  plot_r(pct,log(ml[,1]/ml[,2]),xlab="% state 0",ylab="log(lam0/lam1)",title="speciation")
  plot_r(pct,log(ml[,3]/ml[,4]),xlab="% state 0",ylab="log(mu0/mu1)",title="extinction")
  plot_r(pct,log(ml[,5]/ml[,6]),xlab="% state 0",ylab="log(q01/q10)",title="transition")
  
  ######### MCMC ##########################################
  
  par(mfrow=c(2,2))
  
  txt = gsub("to_server_sim2_","N=",to_server_dir)
  h=hist(res_div,seq(0,100,10),
         main=paste("diversification  ",txt,sep=""),
         xlab = bquote("PP(r"[0] ~ "> r"[1] ~ ")"))   
  lines(x=c(median(res_div),median(res_div)), y=c(0,100), col="blue")
  plot_mtext(res_div)
  
  h=hist(res_lam,seq(0,100,10),
         main=paste("speciation  ",txt,sep=""),
         xlab = bquote("PP(lam"[0] ~ "> lam"[1] ~ ")")) 
  lines(x=c(median(res_lam),median(res_lam)), y=c(0,100), col="blue")
  plot_mtext(res_lam)
  
  h=hist(res_mu,seq(0,100,10),
         main=paste("extinction  ",txt,sep=""),
         xlab = bquote("PP(mu"[0] ~ "> mu"[1] ~ ")"))
  lines(x=c(median(res_mu),median(res_mu)), y=c(0,100), col="blue")
  plot_mtext(res_mu)
  
  h=hist(res_q,seq(0,100,10),
         main=paste("transition  ",txt,sep=""),
         xlab = bquote("PP(q"[01] ~ "> q"[10] ~ ")"))   
  lines(x=c(median(res_q),median(res_q)), y=c(0,100), col="blue")
  plot_mtext(res_q)
  
  
  par(mfrow=c(2,2))
  
  plot_r(pct,res_div,xlab="% state 0",ylab="PP(r0 > r1)",title="diversification")
  plot_r(pct,res_lam,xlab="% state 0",ylab="PP(lam0 > lam1)",title="speciation")
  plot_r(pct,res_mu,xlab="% state 0",ylab="PP(mu0 > mu1)",title="extinction")
  plot_r(pct,res_q,xlab="% state 0",ylab="PP(q01 > q10)",title="transition")
  
  par(mfrow=c(2,2))
  
  plot_r(changes,res_div,xlab="number of state transitions",ylab="PP(r0 > r1)",title="diversification")
  plot_r(changes,res_lam,xlab="number of state transitions",ylab="PP(lam0 > lam1)",title="speciation")
  plot_r(changes,res_mu,xlab="number of state transitions",ylab="PP(mu0 > mu1)",title="extinction")
  plot_r(changes,res_q,xlab="number of state transitions",ylab="PP(q01 > q10)",title="transition")
  
  if (to_pdf)
    dev.off()
}













plot_from_server_ml_factors <- function(to_server_dir, FACTOR=2,NSIM=100,to_pdf=TRUE){ 
  
  if (to_pdf)
    pdf(file=paste("res_",to_server_dir,".pdf",sep=""))
  
  res_div = matrix(NA,ncol=NSIM,nrow=1)
  res_lam = matrix(NA,ncol=NSIM,nrow=1)
  res_mu = matrix(NA,ncol=NSIM,nrow=1)
  res_q = matrix(NA,ncol=NSIM,nrow=1)
  
  pct = matrix(NA,ncol=NSIM,nrow=1)
  changes = matrix(NA,ncol=NSIM,nrow=1)
  
  ml = matrix(NA,nrow=NSIM,ncol=6)
  ml_p = matrix(NA,nrow=NSIM,ncol=2)
  
  all.res = matrix(NA,nrow=NSIM,ncol=10)
  
  factors_flag = matrix(0,nrow=NSIM,ncol=1)
  for (s in 1:NSIM){
    name = paste("Sim2_",s,"_Factor_",2,sep="")
    name1 = paste("Sim2_",s,"_Factor_",1,sep="")
    name4 = paste("Sim2_",s,"_Factor_",4,sep="")
    
    out_file = paste("from_server/",to_server_dir,"/",name,"/",name,".mcmc",sep="")     
    data_file = paste("from_server/",to_server_dir,"/",name,"/","Data.RData",sep="")  
    ml_file = paste("from_server/",to_server_dir,"/",name,"/",name,".ML",sep="")
    
    out_file1 = paste("from_server/",to_server_dir,"/",name1,"/",name1,".mcmc",sep="") 
    out_file4 = paste("from_server/",to_server_dir,"/",name4,"/",name4,".mcmc",sep="") 
    
    if (file.exists(out_file)){
      load(out_file)       
      lam = 100*length(which(samples$lambda0>samples$lambda1))/length(samples$lambda0)
      mu = 100*length(which(samples$mu0>samples$mu1))/length(samples$lambda0)
      div = 100*length(which(samples$div0>samples$div1))/length(samples$lambda0)
      q = 100*length(which(samples$q01>samples$q10))/length(samples$lambda0)
      
      res_div[s] = round(div)
      res_lam[s] = round(lam)
      res_mu[s] = round(mu)
      res_q[s] = round(q)
      
      load(data_file)
      pct[s]=100*length(which(states==0))/length(states)      
      changes[s]=count_changes(phy)
      
      load(ml_file)
      ml[s,] = fit$par 
      ml_p[s,1] = anova(fit, fit.l1)[2,5]
      ml_p[s,2] = anova(fit, fit.l2)[2,5]
      
      load(out_file1)       
      lam1 = 100*length(which(samples$lambda0>samples$lambda1))/length(samples$lambda0)
      mu1 = 100*length(which(samples$mu0>samples$mu1))/length(samples$lambda0)
      div1 = 100*length(which(samples$div0>samples$div1))/length(samples$lambda0)
      q1 = 100*length(which(samples$q01>samples$q10))/length(samples$lambda0)
      load(out_file4)       
      lam4 = 100*length(which(samples$lambda0>samples$lambda1))/length(samples$lambda0)
      mu4 = 100*length(which(samples$mu0>samples$mu1))/length(samples$lambda0)
      div4 = 100*length(which(samples$div0>samples$div1))/length(samples$lambda0)
      q4 = 100*length(which(samples$q01>samples$q10))/length(samples$lambda0)
      
      if (abs(lam1-lam4)>=20 || abs(mu1-mu4)>=20 || abs(div1-div4)>=20 || abs(div1-div4)>=20)
        factors_flag[s] = 1
      
      all.res[s,] = c(lam,lam1,lam4,mu,mu1,mu4,div,div1,div4,pct[s])
      
    }else{
      print(paste("no out file:",out_file))
    }  
  }
    
  #####  ML  #################################################
  
  plot_ml(ml,ml_p)
  
  par(mfrow=c(2,2))
  
  #plot_r(pct,)log((ml[,1]-ml[,3])/(ml[,2]-ml[,4])),xlab="% state 0",ylab="log(lam0-mu0/lam1-mu1)",title="diversification")
  plot_r(pct,log(ml[,1]/ml[,2]),xlab="% state 0",ylab="log(lam0/lam1)",title="speciation")
  plot_r(pct,log(ml[,3]/ml[,4]),xlab="% state 0",ylab="log(mu0/mu1)",title="extinction")
  plot_r(pct,log(ml[,5]/ml[,6]),xlab="% state 0",ylab="log(q01/q10)",title="transition")
  
  ######### MCMC ##########################################
  
  txt = gsub("to_server_sim2_","N=",to_server_dir)  
  
  plot_mcmc(res_div,res_lam,res_mu,res_q,txt) 
  
  par(mfrow=c(2,2))
  
  plot_r(pct,res_div,xlab="% state 0",ylab="PP(r0 > r1)",title="diversification")
  plot_r(pct,res_lam,xlab="% state 0",ylab="PP(lam0 > lam1)",title="speciation")
  plot_r(pct,res_mu,xlab="% state 0",ylab="PP(mu0 > mu1)",title="extinction")
  plot_r(pct,res_q,xlab="% state 0",ylab="PP(q01 > q10)",title="transition")
  
  par(mfrow=c(2,2))
  
  plot_r(changes,res_div,xlab="number of state transitions",ylab="PP(r0 > r1)",title="diversification")
  plot_r(changes,res_lam,xlab="number of state transitions",ylab="PP(lam0 > lam1)",title="speciation")
  plot_r(changes,res_mu,xlab="number of state transitions",ylab="PP(mu0 > mu1)",title="extinction")
  plot_r(changes,res_q,xlab="number of state transitions",ylab="PP(q01 > q10)",title="transition")
  
  ###############################
  
  res_div = res_div[which(factors_flag==0)]
  res_lam = res_lam[which(factors_flag==0)]
  res_mu = res_mu[which(factors_flag==0)]
  res_q = res_q[which(factors_flag==0)]
  pct = pct[which(factors_flag==0)]
  changes = changes[which(factors_flag==0)]
  
  print(paste(sum(factors_flag),"runs filtered"))
  
  txt = gsub("to_server_","after filter, N=",to_server_dir)  
  
  plot_mcmc(res_div,res_lam,res_mu,res_q,txt) 
    
  par(mfrow=c(2,2))
  
  plot_r(pct,res_div,xlab="% state 0",ylab="PP(r0 > r1)",title="diversification")
  plot_r(pct,res_lam,xlab="% state 0",ylab="PP(lam0 > lam1)",title="speciation")
  plot_r(pct,res_mu,xlab="% state 0",ylab="PP(mu0 > mu1)",title="extinction")
  plot_r(pct,res_q,xlab="% state 0",ylab="PP(q01 > q10)",title="transition")
  
  par(mfrow=c(2,2))
  
  plot_r(changes,res_div,xlab="number of state transitions",ylab="PP(r0 > r1)",title="diversification")
  plot_r(changes,res_lam,xlab="number of state transitions",ylab="PP(lam0 > lam1)",title="speciation")
  plot_r(changes,res_mu,xlab="number of state transitions",ylab="PP(mu0 > mu1)",title="extinction")
  plot_r(changes,res_q,xlab="number of state transitions",ylab="PP(q01 > q10)",title="transition")
  
  if (to_pdf)
    dev.off()
}


plot_ml <- function(ml,ml_p){
  par(mfrow=c(2,2))
  
  x=log(ml[,1]/ml[,2])
  h=hist(x,seq(-2,20,0.2),xlim=c(-2,2),
         main="speciation",
         xlab = bquote("log(lam"[0] ~ "/ lam"[1] ~ ")"))   
  lines(x=c(median(x),median(x)), y=c(0,100), col="blue")
  plot_mtext2(x)
  
  x=log(ml[,3]/ml[,4])
  h=hist(x,seq(-20,200,0.2),xlim=c(-2,2),
         main="extinction",
         xlab = bquote("log(mu"[0] ~ "/ mu"[1] ~ ")"))   
  lines(x=c(median(x),median(x)), y=c(0,100), col="blue")
  plot_mtext2(x)
  
  x=log(ml[,5]/ml[,6])
  h=hist(x,seq(-20,200,0.2),xlim=c(-2,2),
         main="transition",
         xlab = bquote("log(q"[0] ~ "/ q"[1] ~ ")"))   
  lines(x=c(median(x),median(x)), y=c(0,100), col="blue")
  plot_mtext2(x)
  
  barplot((100/length(ml_p[,1]))*c(length(which(ml_p[,1]<0.05)),
                                   length(which(ml_p[,2]<0.05))),
          names.arg=c("lambda","mu"),ylim=c(0,100),ylab="p < 0.05 (%)",
          main="significant acceptance of constrained model")
  lines(x=c(0,3), y=c(5,5), col="blue")
}

plot_mcmc <- function(res_div,res_lam,res_mu,res_q,txt){
  
  par(mfrow=c(2,2))  
  
  h=hist(res_div,seq(0,100,10),
         main=paste("diversification  ",txt,sep=""),
         xlab = bquote("PP(r"[0] ~ "> r"[1] ~ ")"))   
  lines(x=c(median(res_div),median(res_div)), y=c(0,100), col="blue")
  plot_mtext(res_div)
  
  h=hist(res_lam,seq(0,100,10),
         main=paste("speciation  ",txt,sep=""),
         xlab = bquote("PP(lam"[0] ~ "> lam"[1] ~ ")")) 
  lines(x=c(median(res_lam),median(res_lam)), y=c(0,100), col="blue")
  plot_mtext(res_lam)
  
  h=hist(res_mu,seq(0,100,10),
         main=paste("extinction  ",txt,sep=""),
         xlab = bquote("PP(mu"[0] ~ "> mu"[1] ~ ")"))
  lines(x=c(median(res_mu),median(res_mu)), y=c(0,100), col="blue")
  plot_mtext(res_mu)
  
  h=hist(res_q,seq(0,100,10),
         main=paste("transition  ",txt,sep=""),
         xlab = bquote("PP(q"[0] ~ "> q"[1] ~ ")"))   
  lines(x=c(median(res_q),median(res_q)), y=c(0,100), col="blue")
  plot_mtext(res_q)
}



plot_t_vs_size <- function(pars,t.MIN,t.MAX){
  
  tt=seq(t.MIN,t.MAX)
  vec=rep(0,length(tt))
  vec_std=rep(0,length(tt))
  N=100
  mat = matrix(0,ncol=length(tt),nrow=N)
  
  ind=1
  for (i in tt){
    v=rep(0,N)  
    for (j in 1:N){
      phy=NULL
      while(is.null(phy)){
        phy <- tree.bisse(pars, max.t=i, x0=0)  
      }
      v[j] = length(phy$tip.label) 
      mat[j,ind] = length(phy$tip.label) 
    }
    vec[ind] = mean(v)
    vec_std[ind] = sd(v)
    ind = ind + 1
    print(c(i,mean(v)))
  }
  
  m = data.frame(mat)
  names(m)=as.character(tt)
  boxplot(m,xlab="max.t",ylab="tree size",ylim=c(0,600))
  abline(h=100)
  mtext(paste("pars:",pars[1],pars[2],pars[3],pars[4],pars[5],pars[6],sep="   "))
  #errbar(tt,vec,vec+vec_std,vec-vec_std,xlab="max.t",ylab="tree size")
}

read_points <- function(name,seconed=FALSE,n.points=100,plot=FALSE){
  load(paste(name,".RData",sep=""))
  print(name)
  res = matrix(-10^100,ncol=7,nrow=n.points)
  missing_points=0
  for (i in 1:n.points){ 
    p=points[i,]  
    name1 = paste("p",i,sep="")
    
    out.file = paste(name,"/",name1,"/res.RData",sep="")  
    if (file.exists(out.file)){      
      load(out.file)
      res[i,1] = l_min
      res[i,2:7] = p_min
    } else{
      #print(paste("No file:",out.file))
      missing_points = missing_points + 1
    }  
  }
  print(paste("Number of missing points:",missing_points))
  
  res[,1] = round(res[,1])
  max_points=sort(unique(res[,1]),decreasing = TRUE)
  p1=res[res[,1]==max_points[1],2:7]
  p1=p1[1,]
  p2=res[res[,1]==max_points[2],2:7]
  
  if (!is.null(dim(p2)))
    p2=p2[1,]
  
  print(paste("p1:",paste(p1,collapse=" "),"  lik:",max_points[1]))
  print(paste("p2:",paste(p2,collapse=" "),"  lik:",max_points[2]))
  
  #   aic1=12-2*max_points[1]
  #   aic2=12-2*max_points[2]
  #   print(paste("AIC - p2 is as probable as p1:",exp((aic1-aic2)/2)))
  
  if (plot){
    plot_it <- function(points,res,ind1,ind2,xlab,ylab,main,POINTS_LIMIT=100){
      par(pty="s")
      plot(points[,ind1],points[,ind2],xlim=c(0,POINTS_LIMIT),ylim=c(0,POINTS_LIMIT),pch=1,xlab=xlab,ylab=ylab,main=main,pty="s",cex=0.3)
      points(points[res[,1]==max_points[1],ind1],points[res[,1]==max_points[1],ind2],col="red",pch=1,cex=0.3)
      points(points[res[,1]==max_points[2],ind1],points[res[,1]==max_points[2],ind2],col="blue",pch=1,cex=0.3)
      points(p1[ind1],p1[ind2],col="red",pch=8,cex=1,lw=1)
      points(p2[ind1],p2[ind2],col="blue",pch=8,cex=1,lw=1)
      text(p1[ind1]+POINTS_LIMIT/12,p1[ind2],max_points[1])
      text(p2[ind1]+POINTS_LIMIT/12,p2[ind2],max_points[2]) 
      mtext(paste("p1=(",round(p1[ind1],2),",",round(p1[ind2],2),")  p2=(",round(p2[ind1],2),",",round(p2[ind2],2),")",sep=""),cex=0.6)
    }
    
    par(mfrow=c(1,3))
    plot_it(points,res,ind1=1,ind2=2,xlab="lambda 0",ylab="lambda 1",main="Speciation")
    plot_it(points,res,ind1=3,ind2=4,xlab="mu 0",ylab="mu 1",main="Extinction")
    plot_it(points,res,ind1=5,ind2=6,xlab="q0->1",ylab="q1->0",main="Transition")
    
    title(paste(name,"find.mle starting points"),line=-1, outer = TRUE )
  }
  if (seconed) return(p2)
  return(p1)
}


read_points_diff <- function(name,n.points=100){
  load(paste(name,".RData",sep=""))
  print(name)
  res = matrix(-10^100,ncol=7,nrow=n.points)
  missing_points=0
  for (i in 1:n.points){ 
    p=points[i,]  
    name1 = paste("p",i,sep="")
    
    out.file = paste(name,"/",name1,"/res.RData",sep="")  
    if (file.exists(out.file)){      
      load(out.file)
      res[i,1] = l_min
      res[i,2:7] = p_min
    } else{
      #print(paste("No file:",out.file))
      missing_points = missing_points + 1
    }  
  }
  print(paste("Number of missing points:",missing_points))
  
  res[,1] = round(res[,1])
  max_points=sort(unique(res[,1]),decreasing = TRUE)
  p1=res[res[,1]==max_points[1],2:7]
  p1=p1[1,]
  p2=res[res[,1]==max_points[2],2:7]
  
  if (!is.null(dim(p2)))
    p2=p2[1,]
  
  print(paste("p1:",paste(p1,collapse=" "),"  lik:",max_points[1]))
  print(paste("p2:",paste(p2,collapse=" "),"  lik:",max_points[2]))
  
  diff=-1
  if (!is.na(p2[1])){
    diff = max_points[1]-max_points[2]
  } 
  return(diff)
}


plot_res <- function(trees_dir, genera, ml_dir, dir1, dir2, ANALYSIS=1,factor_vec=c(1,2,4),to_pdf=TRUE){ 
    
  if (to_pdf) pdf(file=paste("res_a",ANALYSIS,".pdf",sep=""))
  
  FACTOR=2  
  for (t in 1:length(genera)){
    
    # file name
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    
    ml_file = paste(ml_dir,"/",genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"/",
                    genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,".ML",sep="")
    
    mle = NULL
    if (file.exists(ml_file)){
      ml = get_ml_results(ml_file)
      mle = ml[[1]]
      delta.AIC = ml[[2]]
      model = names(delta.AIC[which(delta.AIC==0)])
      
      print(paste("Model:",model))
      print(mle)
      print(delta.AIC)
    }
    
    out_file = paste(dir1,"/",genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"/",
                     genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,".mcmc",sep="") 
    out_ml_file = paste(dir1,"/",genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"/",
                     genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,".ML",sep="") 
    if (file.exists(out_file)){
      load(out_file)
      load(out_ml_file)
      plot_mcmc_results2(samples,paste(genus,"- best ML point"),priorrate,starting.p,mle,model)
      
    }

    out_file = paste(dir2,"/",genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"/",
                     genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,".mcmc",sep="") 
    out_ml_file = paste(dir2,"/",genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"/",
                        genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,".ML",sep="") 
    if (file.exists(out_file)){
      load(out_file)
      load(out_ml_file)
      plot_mcmc_results2(samples,paste(genus,"- starting.point.bisse"),priorrate,starting.p,mle,model)
           
    }
    
  }  
  if (to_pdf) dev.off()
}

prepare_to_server_trees2 <- function(trees_dir,states,to_run = "to_run.txt",to_server_dir = "to_server",
                                     run_file = "run_mcmc.R",NSTEPS = 2000,factor_vec=c(1,2,4),sampling.f=c(1,1),
                                     prior.type="starting.point.bisse",start.type="starting.point.bisse",root.p=NULL,N_TREES=20){  
  
  tree_files = dir(trees_dir)
  
  for (t in 1:length(tree_files)){
    
    # file name
    tree_file=tree_files[t]   
    tree_file_path=paste(trees_dir,"/",tree_files[t],sep="") 
    if (file.exists(tree_file_path)){
      
      print(tree_file)    
      phy_list <- read.tree(tree_file_path)
      
      for (tr in 1:N_TREES){  
        
        phy = phy_list[[tr]]
        states = states[phy$tip.label]
        
        for (FACTOR in factor_vec){    
          name = paste(tree_file,"_Factor",FACTOR,"_Tree",tr,sep="")
          
          server_commands(name,to_run,to_server_dir,run_file,NSTEPS)    
          file.create("Data.RData")
          save(phy, states, name, FACTOR,sampling.f,root.p=root.p,
               prior.type=prior.type,start.type=start.type,file="Data.RData")
          file.copy("Data.RData",paste(to_server_dir,"/",name,sep=""))
          file.copy("functions.R",paste(to_server_dir,"/",name,sep=""))  
          
          if (start.type=="best.point"){
            points_name = paste("points/",name,"_points",sep="")            
            best_point=read_points(points_name,seconed=seconed)
            save(best_point,file="best_point.RData")
            file.copy("best_point.RData",paste(to_server_dir,"/",name,sep=""))
          } 
        }
      }
    } 
  }
}

plot_res_trees <- function(genera, dir1, ANALYSIS=1,factor_vec=c(1,2,4),to_pdf=TRUE, N_TREES=20){ 
  
  if (to_pdf) pdf(file=paste(dir1,"res_a",ANALYSIS,"_trees",".pdf",sep=""))
  
  res = matrix(NA, nrow=length(genera), ncol=5*length(factor_vec))
  
  for (t in 1:length(genera)){
    
    # file name
    genus_file=genera[t]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    print(genus)
    
    factor_ind=0
    for (FACTOR in factor_vec){   
      
      all_samples = NULL
      finished_trees=0
      for (tr in 1:N_TREES){  
        
        name = paste(genus,"_Analyais",ANALYSIS,"_Factor",FACTOR,"_Tree",tr,sep="")
        
        out_file = paste(dir1,"/",name,"/",name,".mcmc",sep="") 
        out_ml_file = paste(dir1,"/",name,"/",name,".ML",sep="") 
        if (file.exists(out_file)){
          finished_trees = finished_trees+1
          load(out_file)
          load(out_ml_file)
          if (!is.element("mu1",colnames(samples)))
            samples$mu1=samples$mu0
            
          if (is.null(all_samples)){ 
            all_samples = samples
          }else{
            all_samples = rbind(all_samples,samples)
          }              
        }      
      }
      if (!is.null(all_samples)){
        lam = 100*length(which(all_samples$lambda0>all_samples$lambda1))/length(all_samples$lambda0)
        mu = 100*length(which(all_samples$mu0>all_samples$mu1))/length(all_samples$lambda0)
        div = 100*length(which(all_samples$div0>all_samples$div1))/length(all_samples$lambda0)
        q = 100*length(which(all_samples$q01>all_samples$q10))/length(all_samples$lambda0)
        #print(seq(1,4)+4*factor_ind)
        res[t,seq(1,5)+5*factor_ind] = c(lam,mu,div,q,finished_trees)
        factor_ind = factor_ind + 1  
        #plot_mcmc_results2(all_samples,paste(genus," Analyais:",ANALYSIS," Factor:",FACTOR," - ",finished_trees,"trees"),priorrate,starting.p)
        plot_mcmc_results_new(all_samples,paste(genus," Analyais:",ANALYSIS," Factor:",FACTOR," - ",finished_trees,"trees"))
      } 
    }  
  }
  if (to_pdf) dev.off()
  df = data.frame(genus=genera,round(res,1))
  colnames(df) = c("genus",rep(c("lam","mu","div","q","finished_trees"),length(factor_vec)))
  write.csv(df,file=paste(dir1,"_res.csv",sep=""))
}

plot_res_trees2 <- function(trees_dir, dir1, factor_vec=c(1,2,4),to_pdf=TRUE, N_TREES=20){ 
  
  if (to_pdf) pdf(file=paste(dir1,"res","_trees",".pdf",sep=""))
  
  tree_files = dir(trees_dir)
  
  res = matrix(NA, nrow=length(tree_files), ncol=4*length(factor_vec))
  
  for (t in 1:length(tree_files)){
    
    # file name
    tree_file=tree_files[t]   
    tree_file_path=paste(trees_dir,"/",tree_files[t],sep="") 
    
    print(tree_file)
    
    factor_ind=0
    for (FACTOR in factor_vec){   
      
      all_samples = NULL
      finished_trees=0
      for (tr in 1:N_TREES){  
        
        name = paste(tree_file,"_Factor",FACTOR,"_Tree",tr,sep="")
        
        out_file = paste(dir1,"/",name,"/",name,".mcmc",sep="") 
        out_ml_file = paste(dir1,"/",name,"/",name,".ML",sep="") 
        if (file.exists(out_file)){
          finished_trees = finished_trees+1
          load(out_file)
          load(out_ml_file)
          if (is.null(all_samples)){ 
            all_samples = samples
          }else{
            all_samples = rbind(all_samples,samples)
          }              
        }      
      }
      if (!is.null(all_samples)){
        lam = 100*length(which(samples$lambda0>samples$lambda1))/length(samples$lambda0)
        mu = 100*length(which(samples$mu0>samples$mu1))/length(samples$lambda0)
        div = 100*length(which(samples$div0>samples$div1))/length(samples$lambda0)
        q = 100*length(which(samples$q01>samples$q10))/length(samples$lambda0)
        #print(seq(1,4)+4*factor_ind)
        res[t,seq(1,4)+4*factor_ind] = c(lam,mu,div,q)
        factor_ind = factor_ind + 1  
        plot_mcmc_results2(all_samples,paste(tree_file," Factor:",FACTOR," - ",finished_trees,"trees"),priorrate,starting.p)
      } 
    }  
  }
  if (to_pdf) dev.off()
  df = data.frame(genus=tree_files,round(res,1))
  colnames(df) = c("tree_file",rep(c("lam","mu","div","q"),length(factor_vec)))
  write.csv(df,file=paste(dir1,"_res.csv",sep=""))
}

plot_mcmc_results2 <- function(samples,name,priorrate,p6,mle=NULL,model=NULL) {   
  
  old.par <- par( no.readonly = TRUE ) 
  col1= rgb(0,0,1,0.2)
  col2= rgb(1,0,0,0.2)
  
  par(mar=c(5,6,5,2))
  par(mfrow=c(2,2))
  
  plot_dens(samples$lambda0,samples$lambda1,"Speciation rate","")
  legend('topright',c(expression(lambda[0]),expression(lambda[1])), fill = c(col2,col1), bty = 'n', border = NA)
  mtext("a", side = 3, line = 1, adj = 0, cex = 1)
  x=seq(0,max(samples$lambda0,samples$lambda1),max(samples$lambda0,samples$lambda1)/50)
  lines(x, dexp(x,priorrate[1]),col="red",lw=2)
  abline(v = p6[1], col = col1,lw=4)
  abline(v = p6[2], col = col2,lw=4)
  
  plot_dens(samples$mu0,samples$mu1,"Extinction rate","")
  legend('topright',c(expression(mu[0]),expression(mu[1])), fill = c(col2,col1), bty = 'n', border = NA)
  mtext("b", side = 3, line = 1, adj = 0, cex = 1)
  x=seq(0,max(samples$mu0,samples$mu1),max(samples$mu0,samples$mu1)/50)
  lines(x, dexp(x,priorrate[3]),col="red",lw=2)
  abline(v = p6[3], col = col1,lw=4)
  abline(v = p6[4], col = col2,lw=4)
  
  plot_dens(samples$div0,samples$div1,"Diversification rate","")
  legend('topright',c(expression(r[0]),expression(r[1])), fill = c(col2,col1), bty = 'n', border = NA)
  mtext("c", side = 3, line = 1, adj = 0, cex = 1)
  #   x=seq(0,max(samples$div0,samples$div1),0.0001)
  #   lines(x, dexp(x,priorrate[1]-priorrate[3]))
  
  plot_dens(samples$q01,samples$q10,"Transition rate","")
  legend('topright',c(expression(q[01]),expression(q[10])), fill = c(col2,col1), bty = 'n', border = NA)
  mtext("d", side = 3, line = 1, adj = 0, cex = 1)
  x=seq(0,max(samples$q01,samples$q10),max(samples$q01,samples$q10)/50)
  lines(x, dexp(x,priorrate[5]),col="red",lw=2)
  abline(v = p6[5], col = col1,lw=4)
  abline(v = p6[6], col = col2,lw=4)  
  
  title( name,line=-1.2, outer = TRUE )
  if (!is.null(mle)){
    mtext( paste("Model:",model,"    p:",paste(as.character(round(mle,2)),sep="",collapse="   ")),
           line=-2.7, outer = TRUE,cex=0.6)
    
  } 
  
  par( old.par )
  
}

get_ml_results <- function(ml_file){
  load(ml_file)       
  
  n.param = 2*c(length(fit.equal.lam.mu.q$par),
                length(fit.equal.mu.q$par),
                length(fit.equal.lambda.q$par),
                length(fit.equal.lam.mu$par),
                length(fit.equal.q$par),
                length(fit.equal.mu$par),
                length(fit.equal.lambda$par),                  
                length(fit.p6$par))
  
  AIC = n.param - 2*c(fit.equal.lam.mu.q$lnLik,
                      fit.equal.mu.q$lnLik,
                      fit.equal.lambda.q$lnLik,
                      fit.equal.lam.mu$lnLik,
                      fit.equal.q$lnLik,
                      fit.equal.mu$lnLik,
                      fit.equal.lambda$lnLik,                         
                      fit.p6$lnLik)      
  
  delta.AIC = AIC - min(AIC)
  exp.delta = exp(-0.5*delta.AIC)
  w.AIC = round(exp.delta/sum(exp.delta),2)
  print(n.param)
  print(AIC)
  print(delta.AIC)
  print(exp.delta)
  print(w.AIC)
  if (which(delta.AIC==0)==8) mle = fit.p6$par; 
  
  if (which(delta.AIC==0)==7) mle = c(fit.equal.lambda$par[1],fit.equal.lambda$par); 
  if (which(delta.AIC==0)==6) mle = c(fit.equal.mu$par[1:3],fit.equal.mu$par[3],fit.equal.mu$par[4:5]);    
  if (which(delta.AIC==0)==5) mle = c(fit.equal.q$par,fit.equal.q$par[5]); 
  
  if (which(delta.AIC==0)==4) mle = c(rep(fit.equal.lam.mu$par[1],2),   rep(fit.equal.lam.mu$par[2],2),   fit.equal.lam.mu$par[3:4]);     
  if (which(delta.AIC==0)==3) mle = c(rep(fit.equal.lambda.q$par[1],2), fit.equal.lambda.q$par[2:3],      rep(fit.equal.lambda.q$par[4],2));
  if (which(delta.AIC==0)==2) mle = c(fit.equal.mu.q$par[1:2],          rep(fit.equal.mu.q$par[3],2),     rep(fit.equal.mu.q$par[4],2));
  if (which(delta.AIC==0)==1) mle = c(rep(fit.equal.lam.mu.q$par[1],2), rep(fit.equal.lam.mu.q$par[2],2), rep(fit.equal.lam.mu.q$par[3],2)); 
 
  names(delta.AIC) = c("M0","Ms","Me","Mq","Mse","Msq","Meq","Mseq")
  return(list(mle,delta.AIC))  
}

write_ml_results <- function(genera,file_name="ml_results.csv",from_server_dir="from_server_ML",ANALYSIS=1){
  
  res = matrix(NA, nrow=length(genera), ncol=14)
  res2 = matrix(NA, nrow=length(genera), ncol=6)
  res3 = matrix(NA, nrow=length(genera), ncol=6)
  
  for( i in 1:length(genera) ){
    genus_file=genera[i]
    genus = strsplit(genus_file,"[.]")[[1]]
    genus = genus[1]
    
    name = paste(genus,"_Analyais",ANALYSIS,"_Factor",2,"_Tree",1,sep="")
    
    ml_file = paste(from_server_dir,"/",name,"/",name,".ML",sep="")
    
   if (file.exists(ml_file)){
      load(ml_file)       
      
      n.param = 2*c(length(fit.equal.lam.mu.q$par),
                    length(fit.equal.mu.q$par),
                    length(fit.equal.lambda.q$par),
                    length(fit.equal.lam.mu$par),
                    length(fit.equal.q$par),
                    length(fit.equal.mu$par),
                    length(fit.equal.lambda$par),                  
                    length(fit.p6$par))
      
      AIC = n.param - 2*c(fit.equal.lam.mu.q$lnLik,
                          fit.equal.mu.q$lnLik,
                          fit.equal.lambda.q$lnLik,
                          fit.equal.lam.mu$lnLik,
                          fit.equal.q$lnLik,
                          fit.equal.mu$lnLik,
                          fit.equal.lambda$lnLik,                         
                          fit.p6$lnLik)      
      
      delta.AIC = AIC - min(AIC)
      exp.delta = exp(-0.5*delta.AIC)
      w.AIC = round(exp.delta/sum(exp.delta),2)
      print(genus)
      print(n.param)
      print(AIC)
      print(delta.AIC)
      print(exp.delta)
      print(w.AIC)
      if (which(delta.AIC==0)==8) mle = fit.p6$par      
      
      if (which(delta.AIC==0)==7) mle = c(fit.equal.lambda$par[1],fit.equal.lambda$par)
      if (which(delta.AIC==0)==6) mle = c(fit.equal.mu$par[1:3],fit.equal.mu$par[3],fit.equal.mu$par[4:5])    
      if (which(delta.AIC==0)==5) mle = c(fit.equal.q$par,fit.equal.q$par[5])
      
      if (which(delta.AIC==0)==4) mle = c(rep(fit.equal.lam.mu$par[1],2),   rep(fit.equal.lam.mu$par[2],2),   fit.equal.lam.mu$par[3:4])      
      if (which(delta.AIC==0)==3) mle = c(rep(fit.equal.lambda.q$par[1],2), fit.equal.lambda.q$par[2:3],      rep(fit.equal.lambda.q$par[4],2))
      if (which(delta.AIC==0)==2) mle = c(fit.equal.mu.q$par[1:2],          rep(fit.equal.mu.q$par[3],2),     rep(fit.equal.mu.q$par[4],2))
      if (which(delta.AIC==0)==1) mle = c(rep(fit.equal.lam.mu.q$par[1],2), rep(fit.equal.lam.mu.q$par[2],2), rep(fit.equal.lam.mu.q$par[3],2))
      
      res[i,] = c(round(mle,4),round(delta.AIC,2))  
      res2[i,] = round(fit.p6$par,4) 
      res3[i,] = round(c(fit.equal.mu$par[1:3],fit.equal.mu$par[3],fit.equal.mu$par[4:5]),4) 
    }
  }
  to_keep = which(!is.na(res[,1]))
  df = data.frame(genera[to_keep],res[to_keep,])
  colnames(df) = c("Group","lambda0","lambda1","mu0","mu1","q01","q10","M0","Ms","Me","Mq","Mse","Msq","Meq","Mseq")
  
  write.csv(df, file=file_name) 
  print(df)
  
  df2 = data.frame(genera[to_keep],res2[to_keep,])
  colnames(df2) = c("genus","lambda0","lambda1","mu0","mu1","q01","q10")
  write.csv(df2, file="ml2.csv") 
  print(df2)
  
  df3 = data.frame(genera[to_keep],res3[to_keep,])
  colnames(df3) = c("genus","lambda0","lambda1","mu0","mu1","q01","q10")
  write.csv(df3, file="ml3.csv") 
  
}



plot_hist_results <- function(genera, dir1,dat, ANALYSIS=1){ 
  
  res <-read.csv(paste(dir1,"_res.csv",sep=""), header = TRUE, stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("NA","") )
  
  load("properties.g")
  colnames(properties.g)[2] <- "n.sp.data"
  
  colnames(res)[3:17]=c("lam_f1","mu_f1","div_f1","q_f1","finished_trees_f1",
                        "lam_f2","mu_f2","div_f2","q_f2","finished_trees_f2",
                        "lam_f4","mu_f4","div_f4","q_f4","finished_trees_f4")
  
  all.res = cbind(properties.g,dat[,2:8],res[3:17])
  head(all.res)
  
  all.res$diff.pct = abs(all.res$div_f4-all.res$div_f1)
  
  analyzed = which(!is.na(all.res$lam_f2))
  selected = which(all.res$diff.pct<=20 & !is.na(all.res$lam_f2))  
  excluded = setdiff(analyzed,selected)
  
  h=hist(all.res[selected,"div_f2"],seq(0,100,10),
         main="",ylab="",xlab="",ylim=c(0,11))
        #, xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  # tt=wilcox.test(all.res[selected,"div_f2"],mu=50)
  #mtext(paste("wilcox.test: p =",round(tt$p.value,4)))
  
  # print(dip.test(all.res[selected,"div_f2"]))
  
}


plot_all_results <- function(genera, dir1,dat, ANALYSIS=1,to_pdf=TRUE){ 
  
  
  if (to_pdf) pdf(file=paste("sum_a",ANALYSIS,dir1,".pdf",sep=""))
  
  res <-read.csv(paste(dir1,"_res.csv",sep=""), header = TRUE, stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("NA","") )
  
  load("properties.g")
  
  colnames(properties.g)[2] <- "n.sp.data"
  
  colnames(res)[3:17]=c("lam_f1","mu_f1","div_f1","q_f1","finished_trees_f1",
                        "lam_f2","mu_f2","div_f2","q_f2","finished_trees_f2",
                        "lam_f4","mu_f4","div_f4","q_f4","finished_trees_f4")
  
  all.res = cbind(properties.g,dat[,2:8],res[3:17])
  head(all.res)
  
  all.res$diff.pct = abs(all.res$div_f4-all.res$div_f1)
  
  save(all.res,file=paste("all.res_a",ANALYSIS,dir1,".RData",sep=""))
  
  analyzed = which(!is.na(all.res$lam_f2))
  selected = which(all.res$diff.pct<=20 & !is.na(all.res$lam_f2))  
  excluded = setdiff(analyzed,selected)
  
  print(paste("all:",length(all.res$genus),
              "  analyzed:",length(analyzed),
              "  afetr Factors filter:",length(selected)))

  res = all.res[selected,]
  
  old.par <- par( no.readonly = TRUE ) 
  
  par(mfrow=c(2,2), pty="s")
  
  h=hist(all.res[selected,"div_f2"],seq(0,100,2.5),
         main=paste("selected genera (",length(selected),")",sep=""),
         xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  tt=wilcox.test(all.res[selected,"div_f2"],mu=50)
  dd=dip.test(all.res[selected,"div_f2"])
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  print(dip.test(all.res[selected,"div_f2"]))
  #print(dip.test(all.res[selected,"div_f2"], simulate=TRUE, B=5000))
  
  ylim=get_ylim()
  
  rr = par('usr')
  hist(all.res[excluded,"div_f2"],seq(0,100,2.5),
       main=paste("excluded genera (",length(excluded),")",sep=""),
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")")
       ,ylim=ylim)  
  
  hist(res[,"lam_f2"],seq(0,100,2.5),
         main="Speciation",
         xlab = bquote("PP(lam"[N] ~ "> lam"[D] ~ ")"))
  tt=wilcox.test(res[,"lam_f2"],mu=50)
  dd=dip.test(res[,"lam_f2"])
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  hist(res[,"mu_f2"],seq(0,100,2.5),
       main="Extinction",
       xlab = bquote("PP(mu"[N] ~ "> mu"[D] ~ ")"))
  tt=wilcox.test(res[,"mu_f2"],mu=50)
  dd=dip.test(res[,"mu_f2"])
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  
  par(mfrow=c(1,1), pty="s")
  
  hist(res[,"q_f2"],seq(0,100,2.5),
       main="Transition",
       xlab = "PP(qND > qDN)")
  tt=wilcox.test(res[,"q_f2"],mu=50)
  dd=dip.test(res[,"q_f2"])
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  par(mfrow=c(2,2), pty="s")
   
  plot(res[,"div_f2"], res[,"n.sp"], 
       main="Diversification vs. # of species",
       xlim = c(0,100),
 #      ylim = c(0,250),
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"),
       ylab = "# species")
  
  temp = res[,"div_f2"]
  temp[temp>50] = 100-temp[temp>50]
  plot(temp, res[,"n.sp"], 
       main="Modified Diversification vs. # of species",
       xlim = c(0,50),
#       ylim = c(0,250),
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"),
       ylab = "# species")
  c=cor.test(temp, res[,"n.sp"])
  mtext(paste("r  =  ",round(c$estimate,2),"   p  =  ",round(c$p.value,4),sep=""))
  
  plot(res[,"pct.0"], res[,"div_f2"], 
       main="Diversification vs. % State 0",
       xlim = c(0,100),
       ylim = c(0,100),
       ylab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"),
       xlab = "N (%)")
  c=cor.test(res[,"pct.0"], res[,"div_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"   p  =  ",round(c$p.value,4),sep=""))
  
  plot(res[,"pct.0"], res[,"q_f2"], 
       main="Transition vs. % State 0",
       xlim = c(0,100),
       ylim = c(0,100),
       ylab = "PP(qND > qDN)",
       xlab = "N (%)")
  c=cor.test(res[,"pct.0"], res[,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"   p  =  ",round(c$p.value,4),sep=""))
  
  
  par(mfrow=c(2,2), pty="s")
  
  vec = res[,"div_f2"]
  hist(vec,seq(0,100,2.5),
       main="Diversification  min 10 species",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  tt=wilcox.test(vec,mu=50)
  dd=dip.test(vec)
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  vec = res[which(res$n.sp>=25),"div_f2"]
  hist(vec,seq(0,100,2.5),
       main="Diversification  min 25 species",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  tt=wilcox.test(vec,mu=50)
  dd=dip.test(vec)
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  vec = res[which(res$n.sp>=50),"div_f2"]
  hist(vec,seq(0,100,2.5),
       main="Diversification  min 50 species",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  tt=wilcox.test(vec,mu=50)
  dd=dip.test(vec)
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  vec = res[which(res$n.sp>=75),"div_f2"]
  hist(vec,seq(0,100,2.5),
       main="Diversification  min 75 species",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  tt=wilcox.test(vec,mu=75)
  dd=dip.test(vec)
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  
  
  par(mfrow=c(2,2), pty="s")
  
  plot(res[,"pct.0"], res[,"q_f2"], 
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "N %",
       ylab = "q(N->D) > q(D->N) % ")
  c=cor.test(res[,"pct.0"], res[,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
  
  plot(res[,"div_f2"], res[,"q_f2"], 
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "r(N) > r(D) %",
       ylab = "q(N->D) > q(D->N) %")
  c=cor.test(res[,"div_f2"], res[,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
 
  plot(res[,"lam_f2"], res[,"q_f2"], 
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "lam(N) > lam(D) %",
       ylab = "q(N->D) > q(D->N) %")
  c=cor.test(res[,"lam_f2"], res[,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
  
  plot(res[,"mu_f2"], res[,"q_f2"], 
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "mu(N) > mu(D) %",
       ylab = "q(N->D) > q(D->N) %")
  c=cor.test(res[,"mu_f2"], res[,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
  
  
  
  plot_col(res,"n.sp","# species")
  plot_col(res,"pct.non.na","species with state info (%)")
  #plot_n(res)
  
  # make woody categorical
  res$woody.c = NA
  res$woody.c[res$woody>=70] = "woody"
  res$woody.c[res$woody<=30] = "herb"
  
  # make perennial categorical
  res$perennial.c = NA
  res$perennial.c[res$perennial>=70] = "perennial"
  res$perennial.c[res$perennial<=30] = "annual"
  
  res$H.pct.c = NA
  if (ANALYSIS==1){
    res$H.pct.c[res$H.pct.1>=70] = "H"
    res$H.pct.c[res$H.pct.1<=30] = "M"
  }
  if (ANALYSIS==2){
    res$H.pct.c[res$H.pct.2>=70] = "H"
    res$H.pct.c[res$H.pct.2<=30] = "M"
  }
  
  
  plot_box(res)
  plot_box(res,"lam_f2",param="lam")
  plot_box(res,"mu_f2",param="mu")
  #plot_box2(res)
  
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","PP(r(N)>r(D))")
  plot_Fruit.dispersal(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","PP(r(N)>r(D))")
  plot_Pollination(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","PP(r(N)>r(D))")
  plot_monocot_dicot(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","PP(r(N)>r(D))")
  
  
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,"pct.0","div_f2","N %","PP(r(N)>r(D))")
  plot_Growth.Form(res,"pct.0","div_f2","N %","PP(r(N)>r(D))")  
  plot_Fruit.dispersal(res,"pct.0","div_f2","N %","PP(r(N)>r(D))")
  plot_Pollination(res,"pct.0","div_f2","N %","PP(r(N)>r(D))")
  
  #plot_H.pct(res,"pct.0","div_f2","N %","PP(r(N)>r(D))")
  
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,"pct.0","lam_f2","N %","PP(lam(N)>lam(D))")
  plot_Growth.Form(res,"pct.0","lam_f2","N %","PP(lam(N)>lam(D))")  
  plot_Fruit.dispersal(res,"pct.0","lam_f2","N %","PP(lam(N)>lam(D))")
  plot_Pollination(res,"pct.0","lam_f2","N %","PP(lam(N)>lam(D))")
  #plot_H.pct(res,"pct.0","lam_f2","N %","PP(lam(N)>lam(D))")
  
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,"pct.0","mu_f2","N %","PP(mu(N)>mu(D))")
  plot_Growth.Form(res,"pct.0","mu_f2","N %","PP(mu(N)>mu(D))")  
  plot_Fruit.dispersal(res,"pct.0","mu_f2","N %","PP(mu(N)>mu(D))")
  plot_Pollination(res,"pct.0","mu_f2","N %","PP(mu(N)>mu(D))")
  #plot_H.pct(res,"pct.0","mu_f2","N %","PP(mu(N)>mu(D))")
  
  print(paste("number of genera with enough data:",length(which(!is.na(all.res$div_f2)))))
  print(paste("number of genera filtered out:",length(excluded)))
  print(paste("number of genera remaining:",length(selected)))
  
  print(paste("mean number of taxa per genus:", mean(res$n.sp)))
  
  print(paste("% genera with PP(r(N)>r(D))) > 50%:", 100*length(which(res$div_f2>50))/length(res$div_f2)))
  
  c=cor.test(res$lam_f2,res$div_f2)
  print(paste("corr(lam,div):   r  =  ",round(c$estimate,2),"          p  =  ",c$p.value,sep=""))
  c=cor.test(res$mu_f2,res$div_f2)
  print(paste("corr(mu,div):   r  =  ",round(c$estimate,2),"          p  =  ",c$p.value,sep="")) 
  
  tt = wilcox.test(res$div_f2,mu=50)
  print(paste("wilcox.test(res$div_f2,mu=50):  p = ",tt$p.value,sep=""))
  tt = wilcox.test(res$lam_f2,mu=50)
  print(paste("wilcox.test(res$lam_f2,mu=50):  p = ",tt$p.value,sep=""))
  tt = wilcox.test(res$mu_f2,mu=50)
  print(paste("wilcox.test(res$mu_f2,mu=50):  p = ",tt$p.value,sep=""))
  tt = wilcox.test(res$q_f2,mu=50)
  print(paste("wilcox.test(res$q_f2,mu=50):  p = ",tt$p.value,sep=""))
  
  par( old.par )
  
  if (to_pdf) dev.off()
}



plot_all_results_sim <- function(genera, dir1,dat, ANALYSIS=1,to_pdf=TRUE){ 
  
  
  if (to_pdf) pdf(file=paste("sum_a",ANALYSIS,dir1,".pdf",sep=""))
  
  res <-read.csv(paste(dir1,"_res.csv",sep=""), header = TRUE, stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("NA","") )
  
  load("properties.g")
  
  colnames(properties.g)[2] <- "n.sp.data"
  
  colnames(res)[3:7]=c("lam_f2","mu_f2","div_f2","q_f2","finished_trees_f2")
  
  all.res = cbind(properties.g[which(is.element(properties.g$genus,dat$genus)),],dat[,2:8],res[3:7])
  head(all.res) 
 
  save(all.res,file=paste("all.res_a",ANALYSIS,dir1,".RData",sep=""))
  
  analyzed = which(!is.na(all.res$lam_f2))
  selected = which(!is.na(all.res$lam_f2))  
  excluded = setdiff(analyzed,selected)
  
  print(paste("all:",length(all.res$genus),
              "  analyzed:",length(analyzed),
              "  afetr Factors filter:",length(selected)))
  
  res = all.res[selected,]
  
  old.par <- par( no.readonly = TRUE ) 
  
  par(mfrow=c(2,2), pty="s")
  
  h=hist(all.res[selected,"div_f2"],seq(0,100,2.5),
         main=paste("Diversification (",length(selected),")",sep=""),
         xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  tt=wilcox.test(all.res[selected,"div_f2"],mu=50)
  dd=dip.test(all.res[selected,"div_f2"])
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  print(dip.test(all.res[selected,"div_f2"]))
  #print(dip.test(all.res[selected,"div_f2"], simulate=TRUE, B=5000))
  
  ylim=get_ylim()
  
  hist(res[,"lam_f2"],seq(0,100,2.5),
       main="Speciation",
       xlab = bquote("PP(lam"[N] ~ "> lam"[D] ~ ")"))
  tt=wilcox.test(res[,"lam_f2"],mu=50)
  dd=dip.test(res[,"lam_f2"])
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  hist(res[,"mu_f2"],seq(0,100,2.5),
       main="Extinction",
       xlab = bquote("PP(mu"[N] ~ "> mu"[D] ~ ")"))
  tt=wilcox.test(res[,"mu_f2"],mu=50)
  dd=dip.test(res[,"mu_f2"])
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  hist(res[,"q_f2"],seq(0,100,2.5),
       main="Transition",
       xlab = "PP(qND > qDN)")
  tt=wilcox.test(res[,"q_f2"],mu=50)
  dd=dip.test(res[,"q_f2"])
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
   
  par(mfrow=c(2,2), pty="s")
  
  plot(res[,"div_f2"], res[,"n.sp"], 
       main="Diversification vs. # of species",
       xlim = c(0,100),
       #      ylim = c(0,250),
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"),
       ylab = "# species")
  
  temp = res[,"div_f2"]
  temp[temp>50] = 100-temp[temp>50]
  plot(temp, res[,"n.sp"], 
       main="Modified Diversification vs. # of species",
       xlim = c(0,50),
       #       ylim = c(0,250),
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"),
       ylab = "# species")
  c=cor.test(temp, res[,"n.sp"])
  mtext(paste("r  =  ",round(c$estimate,2),"   p  =  ",round(c$p.value,4),sep=""))
  
  plot(res[,"pct.0"], res[,"div_f2"], 
       main="Diversification vs. % State 0",
       xlim = c(0,100),
       ylim = c(0,100),
       ylab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"),
       xlab = "N (%)")
  c=cor.test(res[,"pct.0"], res[,"div_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"   p  =  ",round(c$p.value,4),sep=""))
  
  plot(res[,"pct.0"], res[,"q_f2"], 
       main="Transition vs. % State 0",
       xlim = c(0,100),
       ylim = c(0,100),
       ylab = "PP(qND > qDN)",
       xlab = "N (%)")
  c=cor.test(res[,"pct.0"], res[,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"   p  =  ",round(c$p.value,4),sep=""))
  
  
  par(mfrow=c(2,2), pty="s")
  
  vec = res[,"div_f2"]
  hist(vec,seq(0,100,2.5),
       main="Diversification  min 10 species",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  tt=wilcox.test(vec,mu=50)
  dd=dip.test(vec)
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  vec = res[which(res$n.sp>=25),"div_f2"]
  hist(vec,seq(0,100,2.5),
       main="Diversification  min 25 species",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  tt=wilcox.test(vec,mu=50)
  dd=dip.test(vec)
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  vec = res[which(res$n.sp>=50),"div_f2"]
  hist(vec,seq(0,100,2.5),
       main="Diversification  min 50 species",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  tt=wilcox.test(vec,mu=50)
  dd=dip.test(vec)
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  vec = res[which(res$n.sp>=75),"div_f2"]
  hist(vec,seq(0,100,2.5),
       main="Diversification  min 75 species",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  tt=wilcox.test(vec,mu=75)
  dd=dip.test(vec)
  mtext(paste("runksum:",round(tt$p.value,4),"  unimodality:",round(dd$p.value,4)),cex=0.7)  
  
  
  
  par(mfrow=c(2,2), pty="s")
  
  plot(res[,"pct.0"], res[,"q_f2"], 
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "N %",
       ylab = "q(N->D) > q(D->N) % ")
  c=cor.test(res[,"pct.0"], res[,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
  
  plot(res[,"div_f2"], res[,"q_f2"], 
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "r(N) > r(D) %",
       ylab = "q(N->D) > q(D->N) %")
  c=cor.test(res[,"div_f2"], res[,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
  
  plot(res[,"lam_f2"], res[,"q_f2"], 
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "lam(N) > lam(D) %",
       ylab = "q(N->D) > q(D->N) %")
  c=cor.test(res[,"lam_f2"], res[,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
  
  plot(res[,"mu_f2"], res[,"q_f2"], 
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "mu(N) > mu(D) %",
       ylab = "q(N->D) > q(D->N) %")
  c=cor.test(res[,"mu_f2"], res[,"q_f2"])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""))
  
  
  
  plot_col(res,"n.sp","# species")
  plot_col(res,"pct.non.na","species with state info (%)")
  #plot_n(res)
  
  # make woody categorical
  res$woody.c = NA
  res$woody.c[res$woody>=70] = "woody"
  res$woody.c[res$woody<=30] = "herb"
  
  # make perennial categorical
  res$perennial.c = NA
  res$perennial.c[res$perennial>=70] = "perennial"
  res$perennial.c[res$perennial<=30] = "annual"
  
  res$H.pct.c = NA
  if (ANALYSIS==1){
    res$H.pct.c[res$H.pct.1>=70] = "H"
    res$H.pct.c[res$H.pct.1<=30] = "M"
  }
  if (ANALYSIS==2){
    res$H.pct.c[res$H.pct.2>=70] = "H"
    res$H.pct.c[res$H.pct.2<=30] = "M"
  }
  
  
  plot_box(res)
  plot_box(res,"lam_f2",param="lam")
  plot_box(res,"mu_f2",param="mu")
  #plot_box2(res)
  
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","PP(r(N)>r(D))")
  plot_Fruit.dispersal(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","PP(r(N)>r(D))")
  plot_Pollination(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","PP(r(N)>r(D))")
  plot_monocot_dicot(res,paste("H.pct.",ANALYSIS,sep=""),"div_f2","Hermaphrodite %","PP(r(N)>r(D))")
  
  
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,"pct.0","div_f2","N %","PP(r(N)>r(D))")
  plot_Growth.Form(res,"pct.0","div_f2","N %","PP(r(N)>r(D))")  
  plot_Fruit.dispersal(res,"pct.0","div_f2","N %","PP(r(N)>r(D))")
  plot_Pollination(res,"pct.0","div_f2","N %","PP(r(N)>r(D))")
  
  #plot_H.pct(res,"pct.0","div_f2","N %","PP(r(N)>r(D))")
  
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,"pct.0","lam_f2","N %","PP(lam(N)>lam(D))")
  plot_Growth.Form(res,"pct.0","lam_f2","N %","PP(lam(N)>lam(D))")  
  plot_Fruit.dispersal(res,"pct.0","lam_f2","N %","PP(lam(N)>lam(D))")
  plot_Pollination(res,"pct.0","lam_f2","N %","PP(lam(N)>lam(D))")
  #plot_H.pct(res,"pct.0","lam_f2","N %","PP(lam(N)>lam(D))")
  
  par(mfrow=c(2,2), pty="s")
  
  plot_trop(res,"pct.0","mu_f2","N %","PP(mu(N)>mu(D))")
  plot_Growth.Form(res,"pct.0","mu_f2","N %","PP(mu(N)>mu(D))")  
  plot_Fruit.dispersal(res,"pct.0","mu_f2","N %","PP(mu(N)>mu(D))")
  plot_Pollination(res,"pct.0","mu_f2","N %","PP(mu(N)>mu(D))")
  #plot_H.pct(res,"pct.0","mu_f2","N %","PP(mu(N)>mu(D))")
  
  print(paste("number of genera with enough data:",length(which(!is.na(all.res$div_f2)))))
  print(paste("number of genera filtered out:",length(excluded)))
  print(paste("number of genera remaining:",length(selected)))
  
  print(paste("mean number of taxa per genus:", mean(res$n.sp)))
  
  c=cor.test(res$lam_f2,res$div_f2)
  print(paste("corr(lam,div):   r  =  ",round(c$estimate,2),"          p  =  ",c$p.value,sep=""))
  c=cor.test(res$mu_f2,res$div_f2)
  print(paste("corr(mu,div):   r  =  ",round(c$estimate,2),"          p  =  ",c$p.value,sep="")) 
  
  tt = wilcox.test(res$div_f2,mu=50)
  print(paste("wilcox.test(res$div_f2,mu=50):  p = ",tt$p.value,sep=""))
  tt = wilcox.test(res$lam_f2,mu=50)
  print(paste("wilcox.test(res$lam_f2,mu=50):  p = ",tt$p.value,sep=""))
  tt = wilcox.test(res$mu_f2,mu=50)
  print(paste("wilcox.test(res$mu_f2,mu=50):  p = ",tt$p.value,sep=""))
  tt = wilcox.test(res$q_f2,mu=50)
  print(paste("wilcox.test(res$q_f2,mu=50):  p = ",tt$p.value,sep=""))
  
  par( old.par )
  
  if (to_pdf) dev.off()
}


prep <- function(config_new,config_old,new_seq_model,fixed,new_tree_model,seq_file,
                 sub_dir="run",to_run="to_run.txt",main_dir="mrbayes",
                 NSTEPS='1000000',server="jekyl",time_out=864000){ # 864000 = 10 days
  
  dir.create(paste(main_dir,"/",sub_dir,sep=""))
  dir.create(paste(main_dir,"/",sub_dir,"/temp",sep=""))
  
  command_file = paste(main_dir,"/",sub_dir,"/command.txt",sep="")
  file.create(command_file)
  
  dir_path = paste("/groups/itay_mayrose/nivsabath/",main_dir,"/",sub_dir,"/",sep="")
  
  
  cat("/share/apps/mrbayes_3.2.2/src/mb /groups/itay_mayrose/nivsabath/",
      file= command_file,append = FALSE)
  cat(paste(main_dir,"/",sub_dir,"/mb_config.nex\tmrbayes",sep=""),
      file= command_file,append = TRUE)  
  
  
  cat("perl /groups/itay_mayrose/nivsabath/Scripts/run_cmds_in_q_WithNameInCmd_jekyl.pl ",
      file= to_run,append = TRUE)    
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("command.txt ",
      file= to_run,append = TRUE)
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("temp/ itaym 1 yes ",
      file= to_run,append = TRUE)
  cat(time_out,
      file= to_run,append = TRUE)
  cat(" r",
      file= to_run,append = TRUE)
  cat(sub_dir,
      file= to_run,append = TRUE)  
  cat("\n",
      file= to_run,append = TRUE)  
  
  
  file.copy(seq_file,paste(main_dir,"/",sub_dir,"/",sep=""))
  
  config_file = paste(main_dir,"/",sub_dir,"/mb_config.nex",sep="")
  file.create(config_file)
  
  g=strsplit(as.character(config_new),split=" ",fixed=TRUE)
  x<-function(x) x[1]
  config_new1=sapply(g, x, simplify = "vector")
  g=strsplit(as.character(config_old),split=" ",fixed=TRUE)
  x<-function(x) x[1]
  config_old1=sapply(g, x, simplify = "vector")
  
  cat(config_new[1:3], file= config_file,append = FALSE,sep = "\n")
  
  cat("execute ", file= config_file,append = TRUE,sep = "")
  local_seq_file = paste(dir_path,"mb_final_seq.nex;",sep="")
  cat(local_seq_file, file= config_file,append = TRUE,sep = "\n")
  
  cat(config_new[which(config_new1=="constraint")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="charset")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="partition")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="set")], file= config_file,append = TRUE,sep = "\n")
  
  cat("\n", file= config_file,append = TRUE)
  
  if (new_seq_model){
    cat(config_new[which(config_new1=="Lset")], file= config_file,append = TRUE,sep = "\n")
    cat(config_new[which(config_new1=="Prset")], file= config_file,append = TRUE,sep = "\n")
  }else{
    cat(config_old[which(config_old1=="Lset")], file= config_file,append = TRUE,sep = "\n")
    cat(config_old[which(config_old1=="Prset")], file= config_file,append = TRUE,sep = "\n")
  }
  cat("\n", file= config_file,append = TRUE)
  
  if (fixed){
    cat("prset treeagepr=fixed(1);\n", file= config_file,append = TRUE,sep = "\n")   
  }
  
  
  if (new_tree_model){
    cat("prset brlenspr=clock:uniform;\n", file= config_file,append = TRUE,sep = "\n")
  }else{
    cat("prset brlenspr=clock:birthdeath;\n", file= config_file,append = TRUE,sep = "\n")
  }
  cat("prset clockvarpr=igr;\nprset topologypr = constraints(ingroup);\n\n", file= config_file,append = TRUE)
  
  
  cat("mcmc nruns=2 nchains=4 ngen=", file= config_file,append = TRUE)  
  
  if (is.numeric(NSTEPS)){
    print("Warning: numeric NSTEPS could be wrongly printed to config file") 
    cat(NSTEPS)
  } 
  
  cat(NSTEPS,     file= config_file,append = TRUE)
  cat(" samplefreq=2000 diagnfreq=10000 mcmcdiagn=Yes relburnin=yes burninfrac=0.25 Stoprule=Yes Stopval=0.01 Checkpoint=Yes Checkfreq=10000  "
      ,file= config_file,append = TRUE)
  cat("file="    ,file= config_file,append = TRUE)  
  cat(dir_path   ,file= config_file,append = TRUE)  
  cat("mb.out;\n\n",file= config_file,append = TRUE)  
  
  cat("sumt relburnin=yes burninfrac=0.25 conformat=simple;\nend;\n", file= config_file,append = TRUE)
}



prep_ss <- function(config_new, seq_file, brlenspr,
                    sub_dir="run",to_run="to_run.txt",main_dir="mrbayes_ss",
                    server="jekyl",time_out=864000){ 
  
  dir.create(paste(main_dir,"/",sub_dir,sep=""))
  dir.create(paste(main_dir,"/",sub_dir,"/temp",sep=""))
  
  command_file = paste(main_dir,"/",sub_dir,"/command.txt",sep="")
  file.create(command_file)
  
  dir_path = paste("/groups/itay_mayrose/nivsabath/",main_dir,"/",sub_dir,"/",sep="")
  
  
  cat("/share/apps/mrbayes_3.2.2/src/mb /groups/itay_mayrose/nivsabath/",
      file= command_file,append = FALSE)
  cat(paste(main_dir,"/",sub_dir,"/mb_config.nex\t",sub_dir,sep=""),
      file= command_file,append = TRUE)  
  
  
  cat("perl /groups/itay_mayrose/nivsabath/Scripts/run_cmds_in_q_WithNameInCmd_jekyl.pl ",
      file= to_run,append = TRUE)    
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("command.txt ",
      file= to_run,append = TRUE)
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("temp/ itaym 1 yes ",
      file= to_run,append = TRUE)
  cat(time_out,
      file= to_run,append = TRUE)
  cat(" r",
      file= to_run,append = TRUE)
  cat(sub_dir,
      file= to_run,append = TRUE)  
  cat("\n",
      file= to_run,append = TRUE)  
  
  
  file.copy(seq_file,paste(main_dir,"/",sub_dir,"/",sep=""))
  
  config_file = paste(main_dir,"/",sub_dir,"/mb_config.nex",sep="")
  file.create(config_file)
  
  g=strsplit(as.character(config_new),split=" ",fixed=TRUE)
  x<-function(x) x[1]
  config_new1=sapply(g, x, simplify = "vector")
  
  cat(config_new[1:3], file= config_file,append = FALSE,sep = "\n")
  
  cat("execute ", file= config_file,append = TRUE,sep = "")
  local_seq_file = paste(dir_path,"mb_final_seq.nex;",sep="")
  cat(local_seq_file, file= config_file,append = TRUE,sep = "\n")
  
  cat(config_new[which(config_new1=="constraint")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="charset")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="partition")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="set")], file= config_file,append = TRUE,sep = "\n")
  
  cat("\n", file= config_file,append = TRUE)
  
  cat(config_new[which(config_new1=="Lset")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="Prset")], file= config_file,append = TRUE,sep = "\n")
  
  cat("\n", file= config_file,append = TRUE)
  
  cat("prset treeagepr=fixed(1);\n", file= config_file,append = TRUE,sep = "\n")   
  cat("prset clockvarpr=igr;\nprset topologypr = constraints(ingroup);\n\n", file= config_file,append = TRUE)
  
  
  cat("prset brlenspr=clock:",brlenspr,";\n", file= config_file,append = TRUE,sep = "\n")
  
  cat("ss Filename=", file= config_file,append = TRUE)  
  cat(dir_path      ,file= config_file,append = TRUE)  
  cat(paste(brlenspr,".ss.out;\n\n",sep=""),file= config_file,append = TRUE)  
  
  cat("\nend;\n", file= config_file,append = TRUE)
}



prep_mrbayes <- function(config_new,seq_file,
                 sub_dir="run",to_run="mrbayes_to_run.txt",main_dir="mrbayes",
                 NSTEPS='1000000',server="jekyl",time_out=864000){ # 864000 = 10 days
  
  dir.create(paste(main_dir,"/",sub_dir,sep=""))
  dir.create(paste(main_dir,"/",sub_dir,"/temp",sep=""))
  
  command_file = paste(main_dir,"/",sub_dir,"/command.txt",sep="")
  file.create(command_file)
  
  dir_path = paste("/groups/itay_mayrose/nivsabath/",main_dir,"/",sub_dir,"/",sep="")
  
  
  cat("/share/apps/mrbayes_3.2.2/src/mb /groups/itay_mayrose/nivsabath/",
      file= command_file,append = FALSE)
  cat(paste(main_dir,"/",sub_dir,"/mb_config.nex\tmb",sep=""),
      file= command_file,append = TRUE)  
  
  
  cat("perl /groups/itay_mayrose/nivsabath/Scripts/run_cmds_in_q_WithNameInCmd_jekyl.pl ",
      file= to_run,append = TRUE)    
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("command.txt ",
      file= to_run,append = TRUE)
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("temp/ itaym 1 yes ",
      file= to_run,append = TRUE)
  cat(time_out,
      file= to_run,append = TRUE)
  cat(" r",
      file= to_run,append = TRUE)
  cat(sub_dir,
      file= to_run,append = TRUE)  
  cat("\n",
      file= to_run,append = TRUE)  
  
  
  file.copy(seq_file,paste(main_dir,"/",sub_dir,"/",sep=""))
  
  config_file = paste(main_dir,"/",sub_dir,"/mb_config.nex",sep="")
  file.create(config_file)
  
  g=strsplit(as.character(config_new),split=" ",fixed=TRUE)
  x<-function(x) x[1]
  config_new1=sapply(g, x, simplify = "vector")
   
  cat(config_new[1:3], file= config_file,append = FALSE,sep = "\n")
  
  cat("execute ", file= config_file,append = TRUE,sep = "")
  local_seq_file = paste(dir_path,"mb_final_seq.nex;",sep="")
  cat(local_seq_file, file= config_file,append = TRUE,sep = "\n")
  
  cat(config_new[which(config_new1=="constraint")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="charset")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="partition")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="set")], file= config_file,append = TRUE,sep = "\n")
  
  cat("\n", file= config_file,append = TRUE)
  
  cat(config_new[which(config_new1=="Lset")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="Prset")], file= config_file,append = TRUE,sep = "\n")
  
  cat("\n", file= config_file,append = TRUE)
  
  
  cat("prset treeagepr=fixed(1);\n", file= config_file,append = TRUE,sep = "\n")   
  
  cat("prset brlenspr=clock:birthdeath;\n", file= config_file,append = TRUE,sep = "\n")
 
  cat("prset clockvarpr=igr;\nprset topologypr = constraints(ingroup);\n\n", file= config_file,append = TRUE)
    
  cat("mcmc nruns=2 nchains=4 ngen=", file= config_file,append = TRUE)  
  
  if (is.numeric(NSTEPS)){
    print("Warning: numeric NSTEPS could be wrongly printed to config file") 
    cat(NSTEPS)
  } 
  
  cat(NSTEPS,     file= config_file,append = TRUE)
  cat(" samplefreq=2000 diagnfreq=10000 mcmcdiagn=Yes relburnin=yes burninfrac=0.25 Stoprule=Yes Stopval=0.01 Checkpoint=Yes Checkfreq=10000  "
      ,file= config_file,append = TRUE)
  cat("file="    ,file= config_file,append = TRUE)  
  cat(dir_path   ,file= config_file,append = TRUE)  
  cat("mb.out;\n\n",file= config_file,append = TRUE)  
  
  cat("sumt relburnin=yes burninfrac=0.25 conformat=simple;\nend;\n", file= config_file,append = TRUE)
}



prep_mrbayes_config <- function(config_new,seq_file, sub_dir="run",main_dir="mrbayes",
                         NSTEPS='1000000',server="jekyl",time_out=864000){ # 864000 = 10 days
    
  dir_path = paste("/groups/itay_mayrose/nivsabath/",main_dir,"/",sub_dir,"/",sep="")
  
  dir.create(paste("/groups/itay_mayrose/nivsabath/",main_dir,"/",sub_dir,sep=""))
  dir.create(paste(dir_path,"temp",sep=""))
  
  file.copy(seq_file,dir_path)
  
  config_file = paste(dir_path,"mb_config.nex",sep="")
  file.create(config_file)
  
  g=strsplit(as.character(config_new),split=" ",fixed=TRUE)
  x<-function(x) x[1]
  config_new1=sapply(g, x, simplify = "vector")
  
  cat(config_new[1:3], file= config_file,append = FALSE,sep = "\n")
  
  cat("execute ", file= config_file,append = TRUE,sep = "")
  local_seq_file = paste(dir_path,"mb_final_seq.nex;",sep="")
  cat(local_seq_file, file= config_file,append = TRUE,sep = "\n")
  
  cat(config_new[which(config_new1=="constraint")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="charset")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="partition")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="set")], file= config_file,append = TRUE,sep = "\n")
  
  cat("\n", file= config_file,append = TRUE)
  
  cat(config_new[which(config_new1=="Lset")], file= config_file,append = TRUE,sep = "\n")
  cat(config_new[which(config_new1=="Prset")], file= config_file,append = TRUE,sep = "\n")
  
  cat("\n", file= config_file,append = TRUE)
  
  
  cat("prset treeagepr=fixed(1);\n", file= config_file,append = TRUE,sep = "\n")   
  
  cat("prset brlenspr=clock:birthdeath;\n", file= config_file,append = TRUE,sep = "\n")
  
  cat("prset clockvarpr=igr;\nprset topologypr = constraints(ingroup);\n\n", file= config_file,append = TRUE)
  
  cat("mcmc nruns=2 nchains=4 ngen=", file= config_file,append = TRUE)  
  
  if (is.numeric(NSTEPS)){
    print("Warning: numeric NSTEPS could be wrongly printed to config file") 
    cat(NSTEPS)
  } 
  
  cat(NSTEPS,     file= config_file,append = TRUE)
  cat(" samplefreq=2000 diagnfreq=10000 mcmcdiagn=Yes relburnin=yes burninfrac=0.25 Checkpoint=Yes Checkfreq=10000  "
      ,file= config_file,append = TRUE)
  cat("file="    ,file= config_file,append = TRUE)  
  cat(dir_path   ,file= config_file,append = TRUE)  
  cat("mb.out;\n\n",file= config_file,append = TRUE)  
  
  cat("sumt relburnin=yes burninfrac=0.25 conformat=simple;\nend;\n", file= config_file,append = TRUE)
}

run_rnr <- function(run_dir,g,N.GEN1,N.GEN2, burnin_frac=0.25){  
  dir_path = paste("/groups/itay_mayrose/nivsabath/",run_dir,"/",g,"/",sep="")
  
  if (file.exists(paste(dir_path,"mb.out.con.tre",sep=""))){    
    
    # prep to another mr bayes run
    nc= readLines(paste(dir_path,"mb_config.nex",sep=""))  
    
    # get trees from mr bayes
    
    t.file1 = paste(dir_path,"mb.out.run1.t",sep="")
    t.file2 = paste(dir_path,"mb.out.run1.t",sep="")
    
    t1=read.nexus(t.file1)
    t2=read.nexus(t.file2)
    
    t1 = t1[(round(length(t1)*burnin_frac)+1):length(t1)]
    t2 = t2[(round(length(t2)*burnin_frac)+1):length(t2)]
    
    all_trees = c(t1,t2)
    
    all_trees = all_trees[sample(length(all_trees),length(all_trees))]
    
    write.tree(all_trees, file = paste(dir_path,"proc-trees.txt",sep=""), append = FALSE)
    
    
    # run RogueNaRok
    
    dir.create(paste(dir_path,"rnr",sep=""))
    
    system(paste("/share/apps/RogueNaRok-master/RogueNaRok -i ",
                 dir_path,"proc-trees.txt -w ",
                 dir_path,"rnr -n ",g,sep=""))
    
    if (length(readLines(paste(dir_path,"rnr/RogueNaRok_droppedRogues.",g,sep="")))>2){
      # remove taxa from alignment and write the new file to rnr dir
      
      phy = read.nexus(paste(dir_path,"mb.out.con.tre",sep=""))[[2]]
      n.taxa = length(phy$tip.label)
      seqs = read.table(paste(dir_path,"mb_final_seq.nex",sep=""),skip =5,nrows=n.taxa)
      mb_final_seq = readLines(paste(dir_path,"mb_final_seq.nex",sep=""))
      seq.lines = mb_final_seq[6:(5+n.taxa)]
      
      to.drop = read.table(paste(dir_path,"rnr/RogueNaRok_droppedRogues.",g,sep=""),skip =2)
      to.drop = to.drop$V3
      
      pp=strsplit(as.character(to.drop),split="_",fixed=TRUE)
      x<-function(x) x[1]
      pp=sapply(pp, x, simplify = "vector")
      
      to.drop = to.drop[which(pp==g)] # remove outgroup from drop list
      index_to_keep = 5+which(!is.element(seqs$V1,to.drop)) 
      
      # update num taxa
      new.n.taxa = length(index_to_keep)
      txt=unlist(strsplit(mb_final_seq[3],split=" "))
      new.txt = paste( txt[1],gsub("[0-9]+",new.n.taxa, txt[2]), txt[3])
      
      
      new_nex = paste(dir_path,"rnr/mb_final_seq.nex",sep="")
      cat(mb_final_seq[1:2], file= new_nex,append = FALSE,sep = "\n")
      cat(new.txt, file= new_nex,append = TRUE,sep = "\n")
      cat(mb_final_seq[4:5], file= new_nex,append = TRUE,sep = "\n")
      
      cat(mb_final_seq[index_to_keep], file= new_nex,append = TRUE,sep = "\n")
      cat(mb_final_seq[(6+n.taxa):(7+n.taxa)], file= new_nex,append = TRUE,sep = "\n")
      
      #update ingroup    
      for (t in 1:length(to.drop)){
        #nc = gsub( paste(" ",to.drop[t],sep=""),"",nc)
        nc = gsub( paste(" ",to.drop[t]," ",sep="")," ",nc)
      }
      
    }else{
      file.copy(paste(dir_path,"mb_final_seq.nex",sep=""),paste(dir_path,"rnr/mb_final_seq.nex",sep=""))
    }
    
    # update input and output file paths
    nc = gsub(paste(dir_path,"mb_final_seq.nex",sep=""),paste(dir_path,"rnr/mb_final_seq.nex",sep=""),nc)
    nc = gsub(paste(dir_path,"mb.out",sep=""),paste(dir_path,"rnr/mb.out",sep=""),nc)
        
    # update number of generations
    nc = gsub(paste("ngen=",N.GEN1,sep=""),paste("ngen=",N.GEN2,sep=""),nc)
    
    cat(nc, file= paste(dir_path,"rnr/mb_config.nex",sep=""),append = FALSE,sep = "\n")
    
  }else{
    print(paste(g,"did not finish"))
  }
}  


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


mb_server_commands <- function(main_dir,sub_dir,to_run="mb_to_run.txt",run_file="run_mb_pipe.R",server="jekyl",time_out=864000){ # 864000 = 10 days
  
  dir.create(paste(main_dir,"/",sub_dir,sep=""))
  dir.create(paste(main_dir,"/",sub_dir,"/temp",sep=""))
  file.copy(run_file,paste(main_dir,"/",sub_dir,sep=""))
  
  command_file = paste(main_dir,"/",sub_dir,"/command.txt",sep="")
  file.create(command_file)
  
  dir_path = paste("/groups/itay_mayrose/nivsabath/",main_dir,"/",sub_dir,"/",sep="")
  
  cat("/share/apps/R301/bin/R CMD BATCH '--args working_dir=\"",
      file= command_file,append = FALSE)
  cat(dir_path,
      file= command_file,append = TRUE)
  cat("\"' ",
      file= command_file,append = TRUE)  
  cat(dir_path,
      file= command_file,append = TRUE)
  cat(run_file,
      file= command_file,append = TRUE)
  cat(" ",
      file= command_file,append = TRUE)
  cat(dir_path,
      file= command_file,append = TRUE)
  cat(run_file,
      file= command_file,append = TRUE)  
  cat("out\tr",
      file= command_file,append = TRUE)
  cat(sub_dir,
      file= command_file,append = TRUE)  
  
  if (server=="jekyl"){
    cat("perl /groups/itay_mayrose/nivsabath/Scripts/run_cmds_in_q_WithNameInCmd_jekyl.pl ",
        file= to_run,append = TRUE)    
  }else{ # lecs
    cat("perl /groups/itay_mayrose/nivsabath/Scripts/run_cmds_in_q_WithNameInCmd.pl ",
        file= to_run,append = TRUE)    
  }
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("command.txt ",
      file= to_run,append = TRUE)
  cat(dir_path,
      file= to_run,append = TRUE)
  cat("temp/ itaym 1 yes ",
      file= to_run,append = TRUE)
  cat(time_out,
      file= to_run,append = TRUE)
  cat(" r",
      file= to_run,append = TRUE)
  cat(sub_dir,
      file= to_run,append = TRUE)  
  cat("\n",
      file= to_run,append = TRUE)  
  
}


plot_tree <- function(phy,states,ti,type = "phylogram"){
   
  statecols <- c("red","blue") 
  
  plot(phy, label.offset=0.02,tip.color=statecols[states+1], 
       cex=(0.2+(10/length(phy$tip.label))), main=ti,type=type ) 
  
  tiplabels(col=statecols[states+1], pch=19, cex=(0.1+(1/length(phy$tip.label))))          
  
  legend("topleft", c("H/M","Dioecy"),col=statecols,pch=19)      
}

drop_outgroup <-function(phy){
  gg=strsplit(as.character(phy$tip.label),split="_",fixed=TRUE)
  x<-function(x) x[1]
  gg=sapply(gg, x, simplify = "vector")
  phy = drop.tip(phy,phy$tip.label[which(gg!=g)])
  
  phy$edge.length = phy$edge.length/sum(phy$edge.length)
  return(phy)
}

dist_matrix <-function(phy,order){
  d=cophenetic.phylo(phy)
  d = d[order,order]  
  vec=c()
  for (i in 1:dim(d)[1]-1){
    vec=c(vec,d[i,(i+1):dim(d)[1]])
  }
  return(vec)
}

# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


plot_strict_broad <- function(r1,r2,file.name="res.pdf"){
  
  pdf(file.name)
  
  par(mfrow=c(1,2), pty="s")
  
  vec = r1[,"div_f2"]
  
  hist(vec,seq(0,100,2.5),
       main=paste("Strict Diocy (",length(vec),")",sep=""),
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"),cex.main=1.1)
  tt=wilcox.test(vec,mu=50)
  dd=dip.test(vec)
  mtext(paste("mean<>50:",round(tt$p.value,4),"  bimodality:",round(dd$p.value,4)),cex=0.8)  
  mtext("a",at=-27,side=3,line=2,cex=1.1)
  
  vec = r2[,"div_f2"]
  
  hist(vec,seq(0,100,2.5),
       main=paste("Broad Diocy (",length(vec),")",sep=""),
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"),cex.main=1.1)
  tt=wilcox.test(vec,mu=50)
  dd=dip.test(vec)
  mtext(paste("mean<>50:",round(tt$p.value,4),"  bimodality:",round(dd$p.value,4)),cex=0.8)  
  mtext("b",at=-27,side=3,line=2,cex=1.1)
  
  # make woody categorical
  r1$woody.c = NA
  r1$woody.c[r1$woody>=70] = "woody"
  r1$woody.c[r1$woody<=30] = "herb"
  # make perennial categorical
  r1$perennial.c = NA
  r1$perennial.c[r1$perennial>=70] = "perennial"
  r1$perennial.c[r1$perennial<=30] = "annual"
  
  # make woody categorical
  r2$woody.c = NA
  r2$woody.c[r2$woody>=70] = "woody"
  r2$woody.c[r2$woody<=30] = "herb"
  # make perennial categorical
  r2$perennial.c = NA
  r2$perennial.c[r2$perennial>=70] = "perennial"
  r2$perennial.c[r2$perennial<=30] = "annual"
  
  r1$H.pct.c = NA
  r1$H.pct.c[r1$H.pct.1>=70] = "Mostly H"
  r1$H.pct.c[r1$H.pct.1<=30] = "Mostly M"
  r2$H.pct.c = NA
  r2$H.pct.c[r2$H.pct.2>=70] = "Mostly H"
  r2$H.pct.c[r2$H.pct.2<=30] = "Mostly M"
  
  par(mfrow=c(2,2), pty="s")
  plot_trait_box(r1,r2,"woody.c","herb","woody","Growth Form","a")
  plot_trait_box(r1,r2,"trop.dist","temp","trop","Geographical Distribution","b")
  plot_trait_box(r1,r2,"Fruit.dispersal","animals","abiotic","Fruit Dispersal","c")
  plot_trait_box(r1,r2,"Pollination","animals","abiotic","Pollination","d")
  par(mfrow=c(2,2), pty="s")
  plot_trait_box(r1,r2,"H.pct.c","Mostly H","Mostly M","Hermaphrodite %","")
  
  dev.off()
}


plot_trait_box <-function(r1,r2,col,state0,state1,main,letter){
  
  ylab="PP(r(N) > r(D))"
  
  par(cex.axis=0.6)
  
  res=r1
  labels =c(paste(state0 ,"(",length(which(res[,col]==state0)),")",sep=""),
            paste(state1 ,"(",length(which(res[,col]==state1)),")",sep=""))  
  boxplot(res$div_f2[which(res[,col]==state0)],
          res$div_f2[which(res[,col]==state1)],
          main=main,las=0,
          boxwex = 0.25,ylab=ylab,at=1:2 - 0.2,xlim = c(0.5, 4.5),range =0) 
  text(1:2 - 0.2,-14 , srt = 25, adj = 1, labels = labels, xpd = TRUE,cex=0.7)
  z <- wilcox.test(res$div_f2[which(res[,col]==state0)],res$div_f2[which(res[,col]==state1)]) 
  mtext("Strict" ,cex=0.9,at=c(1.3,110),side=1,line=3)
  mtext(paste("p = ",round(z$p.value,4),sep=""),cex=0.7,at=c(1.3,110),side=1,line=4)
  
  res=r2
  labels =c(paste(state0 ,"(",length(which(res[,col]==state0)),")",sep=""),
            paste(state1 ,"(",length(which(res[,col]==state1)),")",sep=""))  
  boxplot(res$div_f2[which(res[,col]==state0)],
          res$div_f2[which(res[,col]==state1)],
          add=TRUE,
          boxwex = 0.25,ylab=ylab,at=3:4 + 0.2,range =0) 
  text(3:4 + 0.2,-14 , srt = 25, adj = 1, labels = labels, xpd = TRUE,cex=0.7)
  z <- wilcox.test(res$div_f2[which(res[,col]==state0)],res$div_f2[which(res[,col]==state1)]) 
  mtext("Broad" ,cex=0.9,at=c(3.7,110),side=1,line=3)
  mtext(paste("p = ",round(z$p.value,4),sep=""),cex=0.7,at=c(3.7,110),side=1,line=4)
  mtext(letter,at=-1,side=3,line=2,cex=1.1)
  
}


run_musse <- function(phy,states,sampling.f,nsteps=1000){
  lik <- make.musse(phy,states, 3,sampling.f = sampling.f)
  p <- starting.point.musse(phy, 3)
  
  if (p[4]==0 || p[1]/p[4]>10){
    p[4:6] = p[1]/10
  }
  
  priorrate = 1/(2*p)
  prior = make.prior.exponential(priorrate)
  
  samples <- mcmc(lik, p, nsteps=nsteps, w=1, lower=0,
                  prior=prior,print.every=1, control=list(backend="CVODES")) 
  # remove burnin
  samples <- subset(samples, i > length(samples[,1])/10)
  samples$div1 = samples$lambda1-samples$mu1
  samples$div2 = samples$lambda2-samples$mu2
  samples$div3 = samples$lambda3-samples$mu3
  return(samples)
}


##################################################################

plot_musse_results <- function(samples,state_names = c("SI","SC","Dio")){
  col1 <- c("red", "orange", "blue")
  par(mfrow=c(3,1), mar=c(4, 4, 2, 1))
  
  profiles.plot(samples[2:4], col.line=col1, xlab="", ylab="prob. density",main="Speciation")
  legend("topright", state_names, col=col1, lty=1)
  profiles.plot(samples[5:7], col.line=col1, xlab="", ylab="prob. density",main="Extinction")
  profiles.plot(samples[2:4]-samples[5:7], col.line=col1, xlab="rate", ylab="prob. density",main="Diversification")
  
  print(paste("PP(lam1>lam2):",round(100*length(which(samples$lambda1>samples$lambda2))/length(samples$lambda1))))
  print(paste("PP(lam1>lam3):",round(100*length(which(samples$lambda1>samples$lambda3))/length(samples$lambda1))))
  print(paste("PP(lam2>lam3):",round(100*length(which(samples$lambda2>samples$lambda3))/length(samples$lambda1))))
  
  print(paste("PP(mu1>mu2):",round(100*length(which(samples$mu1>samples$mu2))/length(samples$mu1))))
  print(paste("PP(mu1>mu3):",round(100*length(which(samples$mu1>samples$mu3))/length(samples$mu1))))
  print(paste("PP(mu2>mu3):",round(100*length(which(samples$mu2>samples$mu3))/length(samples$mu1))))
  
  print(paste("PP(div1>div2):",round(100*length(which(samples$div1>samples$div2))/length(samples$div1))))
  print(paste("PP(div1>div3):",round(100*length(which(samples$div1>samples$div3))/length(samples$div1))))
  print(paste("PP(div2>div3):",round(100*length(which(samples$div2>samples$div3))/length(samples$div1))))
  
}


comp_dirs <- function(dir1,dir2,genera,FACTOR=2,N_TREES=100){
  
  pdf(paste("compare_",dir1,"_",dir2,"_hists.pdf",sep=""))
  par(mfrow=c(2,3), mar=c(4, 4, 4, 1), pty="s")  
  
  col <- c("red", "blue")
  
  v1 = matrix(NA,ncol=3,nrow=length(genera))
  v2 = NULL
  for (i in 1:length(genera)){  
    genus = genera[i]
    print(i)
    v_temp = matrix(NA,ncol=2,nrow=N_TREES)
    all_samples1 = NULL
    all_samples2 = NULL
    for (tr in 1:N_TREES){  
      
      name = paste(genus,"_Analyais",1,"_Factor",FACTOR,"_Tree",tr,sep="")
      
      out_file1 = paste(dir1,"/",name,"/",name,".mcmc",sep="") 
      data_file1 = paste(dir1,"/",name,"/Data.RData",sep="") 
      out_file2 = paste(dir2,"/",name,"/",name,".mcmc",sep="") 
      data_file2 = paste(dir2,"/",name,"/Data.RData",sep="") 
      
      if (file.exists(out_file1)){
        load(out_file1)
        samples1 = samples
        rm(samples)
        
        div1 = round(100*length(which(samples1$div0>samples1$div1))/length(samples1$lambda0),1)      
        v_temp[tr,1] = div1
        
        if (is.null(all_samples1)){ 
          all_samples1 = samples1
        }else{
          all_samples1 = rbind(all_samples1,samples1)
        } 
      }
      if (file.exists(out_file2)){
        load(out_file2)
        samples2 = samples 
        rm(samples)
        
        div2 = round(100*length(which(samples2$div0>samples2$div1))/length(samples2$lambda0),1)           
        
        if (is.null(all_samples2)){ 
          all_samples2 = samples2
        }else{
          all_samples2 = rbind(all_samples2,samples2)
        } 
        
        v_temp[tr,2] = div2
        
#         if (tr<4){
#           par(mfrow=c(1,2), mar=c(2, 2, 2, 1))
#           
#           rm(phy)
#           rm(states)
#           load(data_file1)
#           plot(phy, cex=0.6, font=1,label.offset=0.02,main=paste(genus,div1))
#           tiplabels(col=col[states+1], pch=19)          
#           legend("topleft", c("N","D"),col=col,pch=19) 
#           
#           rm(phy)
#           rm(states)
#           load(data_file2)
#           plot(phy, cex=0.6, font=1,label.offset=0.02,main=paste("simulated",tr,":",div2))
#           tiplabels(col=col[states+1], pch=19)          
#           legend("topleft", c("N","D"),col=col,pch=19) 
#         }
      }
    }
    if (length(which(!is.na(v_temp[,2])))>3){
      div1 = round(100*length(which(all_samples1$div0>all_samples1$div1))/length(all_samples1$lambda0),1)      
      div2 = round(100*length(which(all_samples2$div0>all_samples2$div1))/length(all_samples2$lambda0),1)           
      div3 = min(round(length(which(div1>=v_temp[,2]))/length(v_temp[,2]),3),
                 round(length(which(div1<=v_temp[,2]))/length(v_temp[,2]),3))
      v1[i,] = c(div1,div2,div3)
      
      if (is.null(v2)){ 
        v2 = v_temp
      }else{
        v2 = rbind(v2,v_temp)
      }     
      
      hist(v_temp[,2],seq(0,100,2.5),
           main=genus,
           xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
      points(div1,0,col="red",pch=17,cex=2)
      mtext(paste("p =",div3), line=-1, cex=0.8)
      
    }
  } 
  
  dev.off()
  
  pdf(paste("compare_",dir1,"_",dir2,".pdf",sep=""))
  par(mfrow=c(3,2), mar=c(4, 4, 4, 1), pty="s")
  
  xx = v1[,1]
  hist(xx,seq(0,100,2.5),
       main="True Data",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  d.t <- dip.test(xx)
  mtext(paste("dip test: Dn=",round(d.t$statistic,2),"   p=",round(d.t$p.value,4),sep=""),cex=0.7)
  
  xx = v1[,2]
  hist(xx,seq(0,100,2.5),
       main="Simulated Data",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  d.t <- dip.test(xx)
  #mtext(paste("dip test (Dn,p):",d.t$statistic,d.t$p.value))
  mtext(paste("dip test: Dn=",round(d.t$statistic,2),"   p=",round(d.t$p.value,4),sep=""),cex=0.7)
  xx = v2[,1]
  hist(xx,seq(0,100,2.5),
       main="true data all trees",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  xx = v2[,2]
  hist(xx,seq(0,100,2.5),
       main="simulated data all trees",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  xx = v1[,3]
  hist(xx,seq(0,1,0.025),
       main="p using simulated data",
       ylab = "Count",
       xlab = "p - value")
  
  plot(v1[,1], v1[,2],      
       xlim = c(0,100),
       ylim = c(0,100),
       main = bquote("compare PP(r"[N] ~ "> r"[D] ~ ")"),
       xlab = "PP(True Data)",
       ylab = "PP(Simulated Data)")
  c=cor.test(v1[,1],v1[,2])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""),cex=0.7)
  points(v1[which(v1[,3]<0.025),1], v1[which(v1[,3]<0.025),2],col="red",lwd=2)
  # plot(v2[,1], v2[,2],      
  #      xlim = c(0,100),
  #      ylim = c(0,100),
  #      main="individual trees",
  #      xlab = "True Data",
  #      ylab = "Simulated Data",
  #      pch = 20)
  res=data.frame(genera,v1)
  colnames(res) = c("genus","")
  write.csv(res,file=paste(dir2,"_p.csv",sep=""))
  dev.off()
  
  save(v1,v2,res, file = paste("compare_",dir1,"_",dir2,".RData",sep=""))
  
}



comp_dirs2 <- function(dir1,dir2,genera,FACTOR=2,N_TREES=100){
  
  pdf(paste("compare2_",dir1,"_",dir2,"_hist.pdf",sep=""))
  
  v1 = matrix(NA,ncol=3,nrow=length(genera))
  v2 = NULL
  for (i in 1:length(genera)){  
    genus = genera[i]
    print(i)
    v_temp = matrix(NA,ncol=2,nrow=N_TREES)
    all_samples1 = NULL
    all_samples2 = NULL
    for (tr in 1:N_TREES){  
      
      name = paste(genus,"_Analyais",1,"_Factor",FACTOR,"_Tree",tr,sep="")
      
      out_file1 = paste(dir1,"/",name,"/",name,".mcmc",sep="") 
      data_file1 = paste(dir1,"/",name,"/Data.RData",sep="") 
      out_file2 = paste(dir2,"/",name,"/",name,".mcmc",sep="") 
      data_file2 = paste(dir2,"/",name,"/Data.RData",sep="") 
      
      if (file.exists(out_file1)){
        load(out_file1)
        samples1 = samples
        rm(samples)
        
        div1 = mean(samples1$lambda0/(samples1$lambda0+samples1$lambda1))
        v_temp[tr,1] = div1
        
        if (is.null(all_samples1)){ 
          all_samples1 = samples1
        }else{
          all_samples1 = rbind(all_samples1,samples1)
        } 
      }
      if (file.exists(out_file2)){
        load(out_file2)
        samples2 = samples 
        rm(samples)
        
        div2 = mean(samples2$lambda0/(samples2$lambda0+samples2$lambda1))          
        
        if (is.null(all_samples2)){ 
          all_samples2 = samples2
        }else{
          all_samples2 = rbind(all_samples2,samples2)
        } 
        
        v_temp[tr,2] = div2        
       
      }
    }
    if (length(which(!is.na(v_temp[,2])))>3){
      div1 = mean(all_samples1$lambda0/(all_samples1$lambda0+all_samples1$lambda1))      
      div2 = mean(all_samples2$lambda0/(all_samples2$lambda0+all_samples2$lambda1))           
      div3 = min(round(length(which(div1>=v_temp[,2]))/length(v_temp[,2]),3),
                 round(length(which(div1<=v_temp[,2]))/length(v_temp[,2]),3))
      v1[i,] = c(div1,div2,div3)
      
      if (is.null(v2)){ 
        v2 = v_temp
      }else{
        v2 = rbind(v2,v_temp)
      }  
      
      hist(v_temp[,2],seq(0,1,0.025),
           main=genus,
           xlab = "lamN/(lamN+lamD)")
      points(div1,0,col="red",pch=17,cex=2)
      
    }
  } 
   
  dev.off()
  
  pdf(paste("compare2_",dir1,"_",dir2,".pdf",sep=""))
  par(mfrow=c(3,2), mar=c(4, 4, 4, 1), pty="s")
  
  xx = v1[,1]
  hist(xx,seq(0,1,0.025),
       main="True Data",
       xlab = "mean(lamN/(lamN+lamD))")
  
  xx = v1[,2]
  hist(xx,seq(0,1,0.025),
       main="Simulated Data",
       xlab = "mean(lamN/(lamN+lamD))")
  
  xx = v2[,1]
  hist(xx,seq(0,1,0.025),
       main="true data all trees",
       xlab = "mean(lamN/(lamN+lamD))")
  
  xx = v2[,2]
  hist(xx,seq(0,1,0.025),
       main="simulated data all trees",
       xlab = "mean(lamN/(lamN+lamD))")
  
  xx = v1[,3]
  hist(xx,seq(0,1,0.025),
       main="p using simulated data",
       ylab = "Count",
       xlab = "p - value")
  
  plot(v1[,1], v1[,2],      
       xlim = c(0,1),
       ylim = c(0,1),
       main = "compare mean(lamN/(lamN+lamD))",
       xlab = "True Data",
       ylab = "Simulated Data")
  c=cor.test(v1[,1],v1[,2])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""),cex=0.7)
  points(v1[which(v1[,3]<0.025),1], v1[which(v1[,3]<0.025),2],col="red",lwd=2)
  
  res=data.frame(genera,v1)
  colnames(res) = c("genus","")
  write.csv(res,file=paste(dir2,"_p.csv",sep=""))
  dev.off()
  
  save(v1,v2,res, file = paste("compare2_",dir1,"_",dir2,".RData",sep=""))
  
}


comp_dirs_multi <- function(dir1,dir2,genera,FACTOR=2,N_TREES=20,N_SIM=20){
  
  pdf(paste("compare_",dir1,"_",dir2,"_trees.pdf",sep=""))
  
  col <- c("red", "blue")
  
  v1 = matrix(NA,ncol=2,nrow=length(genera))
  v2 = matrix(NA,ncol=1,nrow=length(genera)*N_TREES)
  v3 = matrix(NA,ncol=1,nrow=length(genera)*N_TREES*N_SIM)
  ind2=1
  ind3=1
  for (i in 1:length(genera)){  
    genus = genera[i]
    print(i)

    all_samples1 = NULL
    
    for (tr in 1:N_TREES){  
      
      name = paste(genus,"_Analyais",1,"_Factor",FACTOR,"_Tree",tr,sep="")
      
      out_file1 = paste(dir1,"/",name,"/",name,".mcmc",sep="") 
      data_file1 = paste(dir1,"/",name,"/Data.RData",sep="") 
            
      if (file.exists(out_file1)){
        load(out_file1)
        samples1 = samples
        rm(samples)
        
        div1 = round(100*length(which(samples1$div0>samples1$div1))/length(samples1$lambda0),1)      
        v2[ind2] = div1
        ind2=ind2+1
        
        if (is.null(all_samples1)){ 
          all_samples1 = samples1
        }else{
          all_samples1 = rbind(all_samples1,samples1)
        } 
      }
      
      all_samples2 = NULL
      for (tr2 in 1:N_SIM){  
        
        name2 = paste(genus,"_Analyais",1,"_Factor",FACTOR,"_Tree",tr,"_sim",tr2,sep="")
        
        out_file2 = paste(dir2,"/",name2,"/",name2,".mcmc",sep="") 
        data_file2 = paste(dir2,"/",name2,"/Data.RData",sep="") 
        
        if (file.exists(out_file2)){
          load(out_file2)
          samples2 = samples 
          rm(samples)
          
          div2 = round(100*length(which(samples2$div0>samples2$div1))/length(samples2$lambda0),1)           
          v3[ind3] = div2
          ind3=ind3+1
          
          if (is.null(all_samples2)){ 
            all_samples2 = samples2
          }else{
            all_samples2 = rbind(all_samples2,samples2)
          } 
           
          if (tr<4 && tr2==1){
            par(mfrow=c(1,2), mar=c(2, 2, 2, 1))
            
            rm(phy)
            rm(states)
            load(data_file1)
            plot(phy, cex=0.6, font=1,label.offset=0.02,main=paste(genus,div1))
            tiplabels(col=col[states+1], pch=19)          
            legend("topleft", c("N","D"),col=col,pch=19) 
            
            rm(phy)
            rm(states)
            load(data_file2)
            plot(phy, cex=0.6, font=1,label.offset=0.02,main=paste("simulated",tr,":",div2))
            tiplabels(col=col[states+1], pch=19)          
            legend("topleft", c("N","D"),col=col,pch=19) 
          }
        }
      }
    }
    
    div1 = round(100*length(which(all_samples1$div0>all_samples1$div1))/length(all_samples1$lambda0),1)      
    div2 = round(100*length(which(all_samples2$div0>all_samples2$div1))/length(all_samples2$lambda0),1)           
    v1[i,] = c(div1,div2)    
  } 
  
  dev.off()
  
  pdf(paste("compare_",dir1,"_",dir2,".pdf",sep=""))
  par(mfrow=c(3,2), mar=c(4, 4, 4, 1), pty="s")
  
  xx = v1[,1]
  hist(xx,seq(0,100,2.5),
       main="True Data",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  xx = v1[,2]
  hist(xx,seq(0,100,2.5),
       main="Simulated Data",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  xx = v2
  hist(xx,seq(0,100,2.5),
       main="true data all trees",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  xx = v3
  hist(xx,seq(0,100,2.5),
       main="simulated data all trees",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
   
  plot(v1[,1], v1[,2],      
       xlim = c(0,100),
       ylim = c(0,100),
       main = bquote("compare PP(r"[N] ~ "> r"[D] ~ ")"),
       xlab = "PP(True Data)",
       ylab = "PP(Simulated Data)")
  c=cor.test(v1[,1],v1[,2])
  mtext(paste("r  =  ",round(c$estimate,2),"          p  =  ",round(c$p.value,4),sep=""),cex=0.7)
  #points(v1[which(v1[,3]<0.025),1], v1[which(v1[,3]<0.025),2],col="red",lwd=2)

  dev.off()
  
  save(v1,v2,v3, file = paste("compare_",dir1,"_",dir2,".RData",sep=""))
  
}


prep_to_csv <- function(dir1,genera,FACTOR=2,N_TREES=1,N_SIM=100,N_TREES2=20){   
  for (i in 1:length(genera)){  
    genus = genera[i]
    for (tr in 1:N_TREES){  
      for (sim in 1:N_SIM){ 
        print(c(i,tr,sim))
        for (tr2 in 1:N_TREES2){  
          name = paste(genus,"_Analyais",1,"_Factor",FACTOR,"_Tree",tr,"_sim",sim,"_tr",tr2,sep="")          
          server_commands2(name,"prep_to_csv.txt",dir1,"to_csv.R",time_out=3600)    
          file.create("Data.RData")         
          save(name,file="Data.RData")
          file.copy("Data.RData",paste(dir1,"/",name,sep=""))          
        }
      }
    }
  }     
} 


plot_sim_res <- function(dir1,genera,FACTOR=2,N_TREES=1,N_SIM=100,N_TREES2=20){
   
  v1 = matrix(NA,ncol=1,nrow=length(genera))
  #v2 = matrix(NA,ncol=1,nrow=length(genera)*N_TREES)
  v3 = matrix(NA,ncol=N_SIM,nrow=length(genera))
  v4 = matrix(NA,ncol=1,nrow=length(genera)*N_TREES*N_SIM*N_TREES2)
 
  ind4=1
  
  for (i in 1:length(genera)){  
    genus = genera[i]
    all_samples = NULL
    #print(i)
    for (tr in 1:N_TREES){  
      for (sim in 1:N_SIM){        
        samples_sim = NULL
        for (tr2 in 1:N_TREES2){  
          name = paste(genus,"_Analyais",1,"_Factor",FACTOR,"_Tree",tr,"_sim",sim,"_tr",tr2,sep="")          
          out_file = paste(dir1,"/",name,"/",name,".mcmc",sep="") 
          #out_file = paste(dir1,"/",name,"/",name,".mcmc.csv",sep="") 
          if (file.exists(out_file)){
            print(c(i,sim,tr2))
            load(out_file)
            #samples = read.csv(out_file)
            div = round(100*length(which(samples$div0>samples$div1))/length(samples$lambda0),1)      
            v4[ind4] = div
            #print(c(ind4,div))
            ind4=ind4+1
            if (is.null(all_samples)){ 
              all_samples = samples
            }else{
              all_samples = rbind(all_samples,samples)
            }             
            if (is.null(samples_sim)){ 
              samples_sim = samples
            }else{
              samples_sim = rbind(samples_sim,samples)
            }
          }
        }        
        div = round(100*length(which(samples_sim$div0>samples_sim$div1))/length(samples_sim$lambda0),1)      
        v3[i,sim] = div        
      }
    }
    
    div1 = round(100*length(which(all_samples$div0>all_samples$div1))/length(all_samples$lambda0),1)      
    v1[i] = div1    
  } 
  
  dip = matrix(NA,ncol=2,nrow=N_SIM)
  for (s in 1:N_SIM){
    d.t <- dip.test(v3[,s])
    dip[s,1] = d.t$statistic
    h=hist(v3[,s],seq(0,100,2.5),plot=FALSE)
    dip[s,2] = h$density[1]+h$density[length(h$density)]
  }
  
  save(v1,v3,v4, file = paste("plot_sim_res_",dir1,".RData",sep=""))
  
  
  pdf(paste("plot_sim_res_",dir1,".pdf",sep=""))
  par(mfrow=c(2,2), mar=c(4, 4, 4, 1), pty="s")  
  
  xx = v4
  hist(xx,seq(0,100,2.5),
       main="per tree",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  dd=dim(v3)
  xx = matrix(v3,1,dd[1]*dd[2])
  hist(xx,seq(0,100,2.5),
       main="per simulation",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
 
  xx = dip[,1]
  hist(xx,seq(0,0.5,0.001),
       main="",
       xlab="dip statistic")
  
  xx = dip[,2]*100
  hist(xx,seq(0,100,2.5),
       main="",
       xlab = bquote("% PP(r"[N] ~ "> r"[D] ~ ") in 5% quantile"))
  
  dev.off()
  
  
}


plot_sim_res3 <- function(dir1,genera,file.name="",FACTOR=2,N_TREES=1,N_SIM=100,N_TREES2=20){
  
  v1 = matrix(NA,ncol=1,nrow=length(genera))
  #v2 = matrix(NA,ncol=1,nrow=length(genera)*N_TREES)
  v3 = matrix(NA,ncol=N_SIM,nrow=length(genera))
  v4 = matrix(NA,ncol=1,nrow=length(genera)*N_TREES*N_SIM*N_TREES2)
  
  ind4=1
  
  for (i in 1:length(genera)){  
    genus = genera[i]
    all_samples = NULL
    #print(i)
    for (tr in 1:N_TREES){  
      for (sim in 1:N_SIM){        
        samples_sim = NULL
        for (tr2 in 1:N_TREES2){  
          name = paste(genus,"_Analyais",1,"_Factor",FACTOR,"_Tree",tr,"_sim",sim,"_tr",tr2,sep="")          
          #out_file = paste(dir1,"/",name,"/",name,".mcmc",sep="") 
          out_file = paste(dir1,"/",name,"/",name,".mcmc.csv",sep="") 
          if (file.exists(out_file)){
            print(c(i,sim,tr2))
            #load(out_file)
            samples = read.csv(out_file)
            div = round(100*length(which(samples$div0>samples$div1))/length(samples$lambda0),1)      
            v4[ind4] = div
            #print(c(ind4,div))
            ind4=ind4+1
            if (is.null(all_samples)){ 
              all_samples = samples
            }else{
              all_samples = rbind(all_samples,samples)
            }             
            if (is.null(samples_sim)){ 
              samples_sim = samples
            }else{
              samples_sim = rbind(samples_sim,samples)
            }
          }
        }        
        div = round(100*length(which(samples_sim$div0>samples_sim$div1))/length(samples_sim$lambda0),1)      
        v3[i,sim] = div        
      }
    }
    
    div1 = round(100*length(which(all_samples$div0>all_samples$div1))/length(all_samples$lambda0),1)      
    v1[i] = div1    
  } 
  
  dip = matrix(NA,ncol=2,nrow=N_SIM)
  for (s in 1:N_SIM){
    d.t <- dip.test(v3[,s])
    dip[s,1] = d.t$statistic
    h=hist(v3[,s],seq(0,100,2.5),plot=FALSE)
    dip[s,2] = 100*(h$counts[1]+h$counts[length(h$counts)])/sum(h$counts)
  }
  
  save(v1,v3,v4, file = paste("plot_sim_res3_",dir1,file.name,".RData",sep=""))
  
  
  pdf(paste("plot_sim_res_",dir1,file.name,".pdf",sep=""))
  par(mfrow=c(2,2), mar=c(4, 4, 4, 1), pty="s")  
  
  xx = v4
  hist(xx,seq(0,100,2.5),
       main="per tree",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  dd=dim(v3)
  xx = matrix(v3,1,dd[1]*dd[2])
  hist(xx,seq(0,100,2.5),
       main="per simulation",
       xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
  
  xx = dip[,1]
  hist(xx,seq(0,0.2,0.005),
       main="",
       xlab="dip statistic")
  
  xx = dip[,2]*100
  hist(xx,seq(0,30,1),
       main="",
       xlab = bquote("% PP(r"[N] ~ "> r"[D] ~ ") outside 95% quantile"))
  
  dev.off()  
}






plot_tree_hist <- function(dir1,genera,FACTOR=2,N_TREES=100){
  
  pdf(paste("plot_tree_hist",dir1,".pdf",sep=""))
  
  par(mfrow=c(2,3), mar=c(4, 4, 4, 1), pty="s")  
  
  v1 = matrix(NA,ncol=3,nrow=length(genera))
  v2 = NULL
  for (i in 1:length(genera)){  
    genus = genera[i]
    print(i)
    v_temp = matrix(NA,ncol=1,nrow=N_TREES)
    all_samples1 = NULL
    
    for (tr in 1:N_TREES){  
      
      name = paste(genus,"_Analyais",1,"_Factor",FACTOR,"_Tree",tr,sep="")
      
      out_file1 = paste(dir1,"/",name,"/",name,".mcmc",sep="") 
      data_file1 = paste(dir1,"/",name,"/Data.RData",sep="") 
      
      if (file.exists(out_file1)){
        load(out_file1)
        samples1 = samples
        rm(samples)
        
        div1 = round(100*length(which(samples1$div0>samples1$div1))/length(samples1$lambda0),1)      
        v_temp[tr] = div1
        
        if (is.null(all_samples1)){ 
          all_samples1 = samples1
        }else{
          all_samples1 = rbind(all_samples1,samples1)
        } 
      }     
    }    
    hist(v_temp,seq(0,100,2.5),
         main=genus,
         xlab = bquote("PP(r"[N] ~ "> r"[D] ~ ")"))
    
  }   
  dev.off()    
}
