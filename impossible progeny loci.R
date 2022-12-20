# identify impossible loci in mother/seedling data 
# the metadata has to have a "tissue" column where samples are  "mother" or "seedling", and a "families" column where
# all family groups are labelled 
# also has to have columns "mother" and "seedling" -- for each seedling sample these columns contain the mother 
# NSW ID and progeny NSW ID respectively 
# see acacia linifolia for example 

# these are Jason's functions that I've partially annotated 
# they can be found in his Ramet git 
table_genotypes <- function( gtp ) {

  tables <- list()

  igt  <- gtp[1,] # get progeny genotype row of gt
  jgt  <- gtp[2,] # get mother genotype row of gt

  ttt <- table(as.numeric(igt), as.numeric(jgt))  # rotate so progeny and mother are columns
  tcn  <- colnames(ttt) # get the progeny and mother ids
  trn  <- rownames(ttt) # get loci

  g3table <- mat.or.vec(3,3) # make a 3x3 matrix of zeros
  rownames(g3table) <- c("2","1","0") # name rows and columns 1,2,3
  colnames(g3table) <- c("2","1","0")

  for (i in 1:nrow(ttt)) { # for each row (locus) in the mother progeny table...
    for (j in 1:ncol(ttt)) { # for each column
      if ( !(trn[i]=="NA" | tcn[j]=="NA") ) { # if neither mother nor progeny are NA at the locus
        g3table[trn[i],tcn[j]] <- ttt[trn[i],tcn[j]] #
      }

    }
  }

  g2table <- mat.or.vec(2,2)
  g2table[1,1] <- g3table[1,1] + g3table[3,3] + g3table[1,3] + g3table[3,1]
  g2table[1,2] <- g3table[2,1] + g3table[2,3]
  g2table[2,1] <- g3table[1,2] + g3table[3,2]
  g2table[2,2] <- g3table[2,2]

  pair <- list(table9cell=g3table, table4cell=g2table)

  return(pair)
}

summary_of_progeny <- function( gt, fam, tissue ) {

   families <- unique(fam)
   ipro     <- which(tissue=="P") # get index numbers of the progeny
   npro     <- length(ipro) # number of progeny
   v        <- rep(0, npro) # make a vector of zeros the length of progeny
   # make an empty df with the columns "mother" "progeny" and ...
   out <- data.frame(mother=as.character(rep("",npro)), progeny=as.character(rep("",npro)), imposs=v, mhet=v, phet=v, mhetpfix=v, mfixphet=v, mhetphet=v)

   class(out$mother) <- "character"
   class(out$progeny) <- "character"
   c <- 1
   for (i in ipro) { # for each of the progeny

         ifam <- fam[i] # get family of progeny
         j    <- which( fam==ifam & tissue == "M") # get the mother index

         igt  <- which(!is.na(gt[i,])) # get the locations of NAs in the progeny row of gt
         jgt  <- which(!is.na(gt[j,])) # get the locations of NAs in the mother row of gt

         gtp <- gt[c(i, j), ] # make a dataframe with the progeny and mother rows of gt

         plist <- table_genotypes(gtp)
         t9c   <- plist$table9cell

         out[c,1] <- rownames(gt)[j]
         out[c,2] <- rownames(gt)[i]
         out[c,3] <- t9c[1,3] + t9c[3,1]
         out[c,4] <- sum(t9c[,2])
         out[c,5] <- sum(t9c[2,])
         out[c,6] <- t9c[1,2] + t9c[3,2]
         out[c,7] <- t9c[2,1] + t9c[2,3]
         out[c,8] <- t9c[2,2]
         c <- c + 1
   }


   return(out)

}


# Calculations ####

# when i originally ran this with the incorrect families / with intentional mistmatches I found that non mother/offspring pairs had ~1200 impossible loci, whereas true offspring/mother pairs have 300-500, and likely clones have <100 
#------using summary_of_progeny

# source("summary_of_progeny.r") # jason's original file with the function 

mfamp       <- as.vector(dms$meta$analyses[,"families"])

mmatp       <- rep("", length(mfamp))

mmatp[which(dms$meta$analyses[,"tissue"] == "mother")] <- "M"

mmatp[which(dms$meta$analyses[,"tissue"] == "seedling")] <- "P"

sump <- summary_of_progeny(dms$gt,mfamp,mmatp)

sump <- mutate(sump, sumhet = mhet + phet)

fams <- data.frame(dms$sample_names,mfamp,mmatp)

fams <- fams[fams$mmatp %in% "P",]

sump <- data.frame(sump,fams)

# Plots ###

hist(sump[["imposs"]])

plot(sumhet~imposs, data=sump)

text(sumhet~imposs, labels=progeny, data=sump, cex=0.5, pos=1)

text(sumhet~imposs, labels=mfamp, data=sump, cex=0.5, pos=3)
clipr::write_clip(sump)