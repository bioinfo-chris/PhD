#####################################################################
### ---===::: Ranavirus Recombination Permuatation Test :::===--- ###
#####################################################################

# permutation test to assess probability of congruent recombination results 
# between observed regions with linkage disequilibrium decay
# and breakpoints detected by GARD in lineages of Ranavirus

# load data
load("PhD/data/rv_recomb_dat.RDS")

# observed proportions of congruent results
bp.has.ld <- sum(rv.recomb.dat$bp$has.ld) / nrow(rv.recomb.dat$bp)
ld.has.bp <- sum(rv.recomb.dat$ld$has.bp) / nrow(rv.recomb.dat$ld)

# set loop variables
rvs <- c("atv","cmtv","fv3","tfv") # ranavirus clades
n <- 1000 # permutation iterations
perm <- data.frame(bp.has.ld.prop=rep(0, n), ld.has.bp.prop=rep(0, n)) # dataframe to record results

# assign p-values and run permutation test
p.bp <- 0 ; p.ld <- 0 ; for(i in 1:n){
  for(rv in rvs){
    
    # get observed results for clade's analyses
    ld.obs <- data.frame(rv.recomb.dat$ld[rv.recomb.dat$ld$clade==rv,c(2:4)]) 
    bp.obs <- rv.recomb.dat$bp[rv.recomb.dat$bp$clade==rv, "bp"]
    
    # set the corresponding random datasets for each iteration
    ld.rand <- data.frame(rv.recomb.dat$core.pos[sample(nrow(rv.recomb.dat$core.pos), nrow(ld.obs)),c(1:3)], has.bp=F)
    bp.rand <- data.frame(bp=sample(c(1:45169), length(bp.obs)), has.ld=F) # 45169 length of alignment
    
    # run the random breakpoints, observed ld
    for(x in 1:nrow(bp.rand)){
      for(y in 1:nrow(ld.obs)){
        bp.rand$has.ld[x] <- ld.obs$start[y] <= bp.rand$bp[x] & ld.obs$end[y] >= bp.rand$bp[x]
        if(bp.rand$has.ld[x]==T){
          break
        }
      }
    }
    assign(paste0(rv,".bp.rand"), bp.rand)
    
    # run the random ld, observed breakpoints
    for(x in 1:nrow(ld.rand)){
      for(y in 1:length(bp.obs)){
        ld.rand$has.bp[x] <- ld.rand$start[x] <= bp.obs[y] & ld.rand$end[x] >= bp.obs[y]
        if(ld.rand$has.bp[x]==T){
          break
        }
      }
    }
    assign(paste0(rv,".ld.rand"), ld.rand)
  }
  
  # record
  rv.bp.rand <- rbind(atv.bp.rand, cmtv.bp.rand, fv3.bp.rand, tfv.bp.rand)
  perm$bp.has.ld.prop[i] <- sum(rv.bp.rand$has.ld)/nrow(rv.bp.rand)
  if(perm$bp.has.ld.prop[i] > bp.has.ld){
    p.bp <- p.bp+(1/n)
  }
  
  rv.ld.rand <- rbind(atv.ld.rand, cmtv.ld.rand, fv3.ld.rand, tfv.ld.rand)
  perm$ld.has.bp.prop[i] <- sum(rv.ld.rand$has.bp)/nrow(rv.ld.rand)
  if(perm$ld.has.bp.prop[i] > ld.has.bp){
    p.ld <- p.ld+(1/n)
  }
} ; p.bp ; p.ld # print p-values


#####################################################################
###                    ---===::: END :::===---                    ###
#####################################################################
