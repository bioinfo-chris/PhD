################################################################
###             Calculate LD Varient Frequencies             ###
################################################################
packages <- c("data.table")
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

# set greater-than frequencies to filter vcf by
freqs <- read.table(args[2])[[1]]
#freqs <- c(0.03,0.04,0.08,0.16)
#freqs <- c(0.16,0.23)
#freqs <- 0.23 # specific freq

# varient freqs
all.vars <- list.files(paste0(ref.path,"/tomahawk/"))
all.vars <- grep("all_vars", all.vars, value=T)
var <- fread(paste0(ref.path,"/tomahawk/",all.vars,"/all_isolates_core_",ref.dir,"_reordered_WGS_positions_SNP_pres-abs_table.txt"), sep="\t")
snps <- var[,10:ncol(var)]
sample.no <- ncol(snps)
snps$alt.count <- rowSums(snps[,1:sample.no])
snps$alt.freq <- rowSums(snps[,1:sample.no])/sample.no
alt.freq <- data.frame(POS=var$POS, REF=var$REF, ALT=var$ALT, ALT_COUNT=snps$alt.count, ALT_FREQ=snps$alt.freq)
write.table(alt.freq, paste0(ref.path,"/tomahawk/",all.vars,"/all_isolates_core_",ref.dir,"_reordered_WGS_positions_varient_freqs.txt"), quote=F, row.names=F)
for (freq in freqs){
  snps.gt.freq <- alt.freq[alt.freq$ALT_FREQ>freq,]
  vars.to.keep <- snps.gt.freq$POS
  no.vars <- length(vars.to.keep)
  dir.create(paste0(ref.path,"/tomahawk/vars_gt",freq,"_n",no.vars))
  write.table(vars.to.keep, paste0(ref.path,"/tomahawk/vars_gt",freq,"_n",no.vars,"/all_isolates_core_",ref.dir,"_vars_gt",freq,"_n",no.vars,".txt"), quote=F, row.names=F, col.names=F)
}

################################################################
###                           END                            ###
################################################################
