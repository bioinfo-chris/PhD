################################################################
###                   Plot Tomahawk LD Output                ###
################################################################
packages <- c("data.table","ggplot2","ggpubr","viridis","cowplot","magrittr")
invisible(lapply(packages, function(x) suppressPackageStartupMessages(require(x, character=T))))
setDTthreads(threads=10)

# function to get path of this file when executed from Rscripts
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly=FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
home.dir <- dirname(thisFile())

# get arguments 
args <- commandArgs(trailingOnly=T)

# arg 1 path to reference directory releative to where script is
arg.path <- args[1]
ref.path <- paste0(home.dir,"/",arg.path)
setwd(ref.path)
ref.dir <- gsub("^.*/", "", ref.path)

# get project directory for plot titles from arg path
proj.title <- gsub("_", " ", strsplit(arg.path, "/")[[1]][1])

# get frequency directories (exclude all_vars)
var.dirs <- list.files(paste0(ref.path,"/tomahawk"))
freq.dirs <- grep("vars_gt",var.dirs,value=T)
#freq.dirs <- grep("vars_gt0.23",var.dirs,value=T) # specific freq set

# plots, regressions and permutation tests
for (freq.dir in freq.dirs){
  freqn <- gsub("vars_","",freq.dir)
  freq <- strsplit(freqn, "_")[[1]][1] %>% gsub("^.{2}","",.)
  ref <- gsub("ref_","",ref.dir)
  message ("processing variants at frequency > ", freq)
  ld.table <- read.table(paste0(ref.path,"/tomahawk/",freq.dir,"/all_isolates_core_",ref.dir,"_reordered_WGS_positions_",freqn,"_diploid_LD_table.txt"), header=T, sep="\t")
  ld.table$dist <- ld.table$posB-ld.table$posA
  ld.table$dist.pos <- as.numeric(abs(ld.table$dist))
  
  # linear models
  lm.r2 <- summary(lm(R2~dist.pos, data=ld.table))
  capture.output(lm.r2, file=paste0(ref.path,"/tomahawk/",freq.dir,"/all_isolates_core_",ref.dir,"_reordered_WGS_positions_",freqn,"_r2-dist_lm.txt"))
  lm.dprime <- summary(lm(Dprime~dist.pos, data=ld.table))
  capture.output(lm.dprime, file=paste0(ref.path,"/tomahawk/",freq.dir,"/all_isolates_core_",ref.dir,"_reordered_WGS_positions_",freqn,"_dprime-dist_lm.txt"))
  
  # permutation test
  obs.coef <- lm.r2$coefficients[2,1]
  n <- 1000
  pvalue <- 1/n
  perm.test <- data.frame(it=seq(1:n), coef=rep(0, n), r2=rep(0, n))
  for (i in 1:n){
    dist.pos2 <- sample(ld.table$dist.pos, replace=F)
    lm.r2.perm <- summary(lm(ld.table$R2 ~ dist.pos2))
    perm.test[perm.test$it==i,"r2"] <- as.numeric(lm.r2.perm$r.squared)
    perm.test[perm.test$it==i,"coef"] <- as.numeric(lm.r2.perm$coefficients[2,1])
    if(perm.test$coef[i] <= obs.coef){
      pvalue <- pvalue + 1/n
    }
  }
  pvalue <- signif(pvalue, digits=3)
  saveRDS(perm.test, paste0(ref.path,"/tomahawk/",freq.dir,"/all_isolates_core_",ref.dir,"_reordered_WGS_positions_",freqn,"_permutation_test_R_data.RDS"))
  pvalue.res <- paste0(proj.title," permutation P-value: ", pvalue)
  write.table(pvalue.res, file=paste0(ref.path,"/tomahawk/",freq.dir,"/",ref.dir,"_",freq.dir,"_permutation_pvalue.txt"), col.names=F, row.names=F, quote=F)
  
  # LD decay curves
  d.plot.title <- paste0(proj.title, " Linkage Disequilibrium (D) Decay of Variants at Frequency > ", freq)
  plot.subtitle <- paste0("ORF order and direction aligned to reference sequence ", ref, "\n")
  
  d.plot <- ggplot(ld.table,aes(x=dist.pos,y=D)) + theme_classic() +
    geom_point(aes(colour=D), alpha=.2, shape=19, size=3) +
    scale_colour_gradientn(colours=c("#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF", "#FE9F6DFF", "#DE4968FF", "#8C2981FF", "#3B0F70FF", "#000004FF")) +
    scale_x_continuous(breaks=seq(0, 110000, 10000), labels=as.character(seq(0, 110, 10))) +
    xlab("Distance (Kbp)") + ylab("D") +
    ggtitle(d.plot.title, subtitle = plot.subtitle) +
    theme(legend.position = "none",
          plot.title = element_text(size = 20, face="bold",family="Helvetica", colour="azure2", hjust = .5),
          plot.subtitle = element_text(size = 14,family="Helvetica", colour="azure2", hjust = .5),
          axis.title.x = element_text(size=14, vjust=-1, face="bold",family="Helvetica", colour="azure2"),
          axis.title.y = element_text(size=14, vjust=2, face="bold",family="Helvetica", colour="azure2"),
          axis.text = element_text(size=12, colour="azure2"),
          panel.background = element_rect(fill = "black", colour = NA),
          plot.background = element_rect(fill = "black"),
          axis.line = element_line(size = 0.5, linetype = "solid",colour = "azure2"),
          axis.ticks = element_line(size = 0.5, linetype = "solid",colour = "azure2"),
          plot.margin = margin(.8, .8, .8, .8, "cm"))
  ggsave(paste0(ref.path,"/tomahawk/",freq.dir,"/all_isolates_core_",ref.dir,"_reordered_WGS_positions_",freqn,"_D.pdf"), d.plot, height=9, width=16)
  
  # permutation plots
  perm.plot.title <- paste0(proj.title, " LD r2 Decay Regression Permutation Test of Variants at Frequency > ", freq, ": P-value = ", pvalue)
  perm.plot.subtitle <- paste0("Regression result distributions of SNP coordinates randomly permuted 1,000 times, verses observed coordinates (orange bars)\nLD anlaysis conducted using reference sequence ", ref, "\n")
  
  r2.plot <- ggplot(ld.table,aes(x=dist.pos,y=R2)) + theme_bw() +
    geom_point(shape=19, colour="black", alpha=.3, size=3) +
    geom_smooth(method = "lm", formula = y~x, size = 3, color="darkorange") +
    scale_x_continuous(breaks=seq(0, 110000, 10000), labels=as.character(seq(0, 110, 10))) +
    xlab("Distance (Kbp)") + ylab(bquote('r'^2~'')) +
    ggtitle(perm.plot.title) +
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size=14, vjust=-1),
          axis.title.y = element_text(size=14, vjust=2),
          axis.text = element_text(size=12),
          plot.margin = margin(.8, .8, .2, .8, "cm"))
  
  frac <- lm.r2$r.squared/20
  if(lm.r2$r.squared > max(perm.test$r2)){
    max <- lm.r2$r.squared+frac
  } else {
    max <- max(perm.test$r2)+frac
  }
  
  perm.r2 <- ggplot(perm.test, aes(x=r2)) + theme_bw() +
    scale_x_continuous(limits=c(min(perm.test$r2)-frac, max)) +
    geom_histogram(bins=150, colour="white", fill="turquoise4") +
    geom_vline(xintercept=lm.r2$r.squared, colour="darkorange", size=2) +
    geom_hline(yintercept=0, colour="grey", size=.2, linetype="dashed") +
    xlab(bquote('Regression R'^2~'')) + ylab("Frequnecy") +
    theme(axis.title.x = element_text(size=14, vjust=-1),
          axis.title.y = element_text(size=14, vjust=3),
          axis.text = element_text(size=12),
          panel.grid.major   = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_line(colour="grey", size=.2, linetype="dashed"),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(.8, .8, .2, .8, "cm"))
  
  perm.coef <- ggplot(perm.test, aes(x=coef)) + theme_bw() +
    geom_histogram(bins=150, colour="white", fill="turquoise4") +
    geom_vline(xintercept=obs.coef, colour="darkorange", size=2) +
    geom_hline(yintercept=0, colour="grey", size=.2, linetype="dashed") +
    xlab("Regression Coefficient") + ylab("Frequnecy") +
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 14),
          axis.title.x = element_text(size=14, vjust=-1),
          axis.title.y = element_text(size=14, vjust=3),
          axis.text = element_text(size=12),
          panel.grid.major   = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_line(color="grey", size=.2, linetype="dashed"),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(.8, .8, .8, .8, "cm"))
  
  suppressWarnings(comb <- plot_grid(r2.plot, perm.r2, perm.coef, ncol=1, align="v"))
  ggsave(paste0(ref.path,"/tomahawk/",freq.dir,"/all_isolates_core_",ref.dir,"_reordered_WGS_positions_",freqn,"_decay_permutation_test.pdf"), comb, height=9, width=16)
}

################################################################
###                           END                            ###
################################################################
