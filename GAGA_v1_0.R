## GEM and GWAS Automisation Pipeline ##

## This script takes in multiple trait data and runs GEM and GWAS analyses 
## on them, outputting Manhattan plots.
## The outputted data is organised in files with trait data names.
rm(list = ls())

## Navigate to the Data folder containing the full trait dataset in the format
## demonstrated in the documentation.

# setwd("~PATH/TO/YOUR/DIRECTORY")
setwd("~/Documents/GAGA/demo/")
runstats <- TRUE
colno <- 3
rungem <- TRUE
rungwas <- TRUE
gwasmodels <- c('GLM', 'FarmCPU', 'Blink')

## Tick "Source on Save" above and hit save to run!

################################################################################
## Loading packages
################################################################################

print("Loading dependencies...")

#For stats and GEM
library(lmerTest)
library(readxl)
library(pbkrtest)
library(emmeans)
library(broom)
library(plyr)
library(dplyr)
library(EnvStats)
library(qqman)
library(ggplot2)
library(MASS) 
library(gplots)
library(compiler) 
library(tidyverse)
library(data.table) 

# For Plotter
library(shiny) 
library(qvalue) 

# For GAPIT
library(BiocManager)
library(GAPIT3)
library(multtest)
library(scatterplot3d)
library(gplots)
 

################################################################################
## Reading in files
################################################################################

print("Reading in necessary files...")

## Read in delimited files with sequence names and rpkm values
rpkmA.temp <- read.table (file = paste(dir(pattern = "unigenefile_")), header=TRUE,
                          sep= "\t", row.names=1)
meansCutOff = 0.4
delrpkmA = rpkmA.temp[rowMeans(rpkmA.temp) >= meansCutOff, ]
rpkmA <- t(delrpkmA)

## Read B. napus to Arabidopsis matches
codesFile = read.csv(paste(dir(pattern = "marker_")), header=TRUE, na.strings = ".")

## Genotype data for GAPIT
myG <- read.delim(paste(dir(pattern = "snpfile_")), header = FALSE)
## Kinship data for GAPIT
myCV <- read.table(paste(dir(pattern = "qmatrix_")), header = TRUE)

septraits <- read.csv(dir(pattern = "totest_"), header = TRUE)

# Checking all files loaded correctly
expectedfiles <- c("codesFile", "myCV", "myG", "rpkmA", "septraits")
currentfiles <- ls()
filespresent <- as.character(expectedfiles %in% currentfiles)

if(identical(filespresent, c(rep("TRUE", 5)))) {
  print("Files loaded successfully.")
} else {
  stop("Failed to load files. Please check you have titled and formatted all files correctly before proceeding.")
}

print("Sourcing GAPIT scripts...")
source("gapit_functions.txt")

################################################################################
## Seperating traits into individual files and folders.
## These folders will be organised in a master folder named after your trait  
## input file followed by "_Results".
################################################################################

print("Separating traits into individual files and folders...")

## Setting up the master folder.
dataname <- noquote(dir(pattern = "totest_"))
dataname <- gsub("\\..*","",dataname)
dataname <- gsub("totest_","",dataname)

masterfolder <- paste(dataname, "_Results", sep = "")
dir.create(paste(masterfolder))

## Separating traits into individual files

septraits <- septraits[, colSums(is.na(septraits)) != nrow(septraits)]
traitno <- length(colnames(septraits))-colno

for (j in 1:traitno) {
  # Generating trait folder
  dir.create(paste(masterfolder, "/", colnames(septraits)[j+colno], sep = ""))
  setwd(paste(masterfolder, "/", colnames(septraits)[j+colno], sep = ""))
  
  # Generating trait file
  traittosave <- septraits %>% dplyr::select(1:colno, j+colno)
  colnames(traittosave)[colno+1] <- "trait"
  write.table(traittosave, paste("raw_", colnames(septraits)[j+colno], ".txt",
                               sep = ""), sep = "\t", 
              quote = FALSE, row.names = FALSE)
  
  # Running stats test if runstats == TRUE
  if(runstats == TRUE) {
    
    traitname <- colnames(septraits)[j+colno]
    print(paste("Running stats test for ", traitname, "...", sep = ""))

    MyDataframe <- traittosave
    MyDataframe$trait <- as.numeric(MyDataframe$trait)
    
    PValues <- c(as.numeric(shapiro.test(MyDataframe$trait)[2]), 
                 as.numeric(shapiro.test(log10(MyDataframe$trait))[2]), 
                 as.numeric(shapiro.test(sqrt(MyDataframe$trait))[2]))
    
    if (PValues[1] > 0.05 & PValues[2] > 0.05 & PValues[3] > 0.05) {
      StatType <- "None"
    } else {
      Phenotype_Vals <- list(MyDataframe$trait, log10(MyDataframe$trait), sqrt(MyDataframe$trait))[[which.min(replace((0.05-PValues), PValues<=0, NA))]]
      StatType <- list("None", "Logged", "Rooted")[[which.min(replace((0.05-PValues), PValues<=0, NA))]]
      }
    MyDataframe$trait <- Phenotype_Vals
  
    ## Sets up tables to collect LMM p-values and estimate means
    LMMsResults <- data.frame(matrix(ncol = 3, nrow = 1)) # p-values
    colnames(LMMsResults) <- c("Trait", "F.value", "Pr(>F)")
    LMMsResults$Trait <- c(traitname)
    EstimateMeans <- tibble("genotype_id" = 
                              c(unique(MyDataframe$genotype_id))) # estimate means
    EstimateMeans <- arrange(EstimateMeans, genotype_id)
    
    # runs linear mixed model for the trait with location and 
    # genotypes as random factors
    LMMmod<-lmer(trait~(1|location)*genotype_id, data = MyDataframe)
    predgeno <- emmeans(LMMmod, ~ genotype_id)
    trait<-  tidy(predgeno)
    
    # runs ANOVA to get p-values and writes them to a table
    LMManova<- anova(LMMmod, ddf = "K")
    LMMsResults$F.value[1] <- LMManova$`F value`[1]
    LMMsResults$`Pr(>F)`[1] <- LMManova$`Pr(>F)`[1]
    write.csv(as.data.frame(trait), "LMMresult.csv", row.names = F)
    
    # puts estimate means for that trait into a table
    trait <- trait[-c(3:6)]
    colnames(trait) <- c("genotype_id", traitname)
    EstimateMeans <- left_join(EstimateMeans, trait, by = "genotype_id")
    colnames(EstimateMeans)[1] <- "<Trait>"
    EstimateMeans <- filter(EstimateMeans, !is.na(EstimateMeans[2]))
    
    write.table(EstimateMeans,"LMMestimatemeans.txt", sep = "\t", 
                quote = FALSE, row.names = FALSE)
    
    if (length(dir(pattern = "LMMestimatemeans.txt")) == 1) {
      print(paste("Estimate means produced for ", traitname, ".", sep = ""))
    } else {
      stop("Estimate means failed.")
    }

    # 
    # return(DataNormalisation)
    
  } else { # if runstats == FALSE renaming file for GEM and GWAS
    print(paste("No stats testing for ", traitname, "...", sep = ""))
    file.rename(paste("raw_", traitname,".txt", sep=""), "LMMestimatemeans.txt")
  }
  setwd("../../")
}

trait.list <- dir(paste(masterfolder))

# Checking trait folders created match trait names in input file
traitnamesinseptraits <- sort(colnames(septraits[4:length(colnames(septraits))]))

if (identical(traitnamesinseptraits, trait.list)){
  print(paste("Traits successfully separated and organised in the folder, ", 
            masterfolder, ".", sep = ""))
} else {
  stop("There has been a problem setting up the trait folders. Please check the layout of your trait input file.")
}

################################################################################
## GEM ANALYSIS
## GEM analysis. Main output is "GEM-Results-[traitname].txt".
################################################################################

GEManalysis <- function(){
  
  print("Starting GEM Analysis...")
  
  traitname <- getwd()
  traitname <- gsub(".*/","",traitname)
  
  ## Read in delimited file with sequence identifiers and trait values
  traitData <- read.table("LMMestimatemeans.txt", header=TRUE, sep= "\t", 
                          row.names = 1) 
  row.names(traitData) <- gsub("-", "", row.names(traitData))
  
  ## Delete missing varieties from rpkm files
  print("Merging data")
  rpkmmergeA <- merge(traitData, rpkmA, by='row.names')
  rownames(rpkmmergeA) <- rpkmmergeA[,1]
  
  print("Performing linear regression")
  ## Function for apply()ing. Does a linear model between trait data and genes
  performLinearRegression <- function(exp_vector) {
    model <- lm(rpkmmergeA[,2] ~ exp_vector)
    return(as.numeric(anova(model)[1,]))
  }
  
  ## Do the linear models. Each gene is in a column hence apply() over
  ## columns. rpkmmergeA contains lines in column 1 and the trait data
  ## in column 2 hence these are dropped in the apply() call. The beauty
  ## of this way is that it makes the row names the genes in the final
  ## lmResults output.
  lmResults <- t(apply(rpkmmergeA[, -c(1,2)], 2, performLinearRegression))
  colnames(lmResults) <- c("Df", "SumSq", "MeanSq", "Fvalue", "Pvalue")
  
  ## calculate log10P
  log10P <- -log10(lmResults[,"Pvalue"])
  print("Calculated p-values")
  
  ## Get everything together for writing to a file. Made an explicit
  ## rownames column so don't have to use row.names in next regress
  ## plotter step. Not a problem because I use row.names=FALSE in
  ## write.table. Note needed to use as.data.frame otherwise it all came
  ## out as a list of character vectors i.e. in quotes which obviously
  ## can't be order()ed.
  getGeneModels = function(results) {
    unigenes = unlist(strsplit(sub("_", "~#~", results$unigene), "~#~"))
    ArabidopsisHits = codesFile[ codesFile$unigene %in% unigenes, ]
    results$tmp = unigenes
    resultsWithAGI = merge(results, ArabidopsisHits, by.x="tmp", by.y="unigene")  
    ## likely in the future can just merge by="unigene" as identical(as.character(results$unigene), results$tmp) == TRUE
    return(resultsWithAGI)
  }
  
  BnapusUnigene = rownames(lmResults)
  traitAndUnigene = as.data.frame(cbind(trait = traitname, unigene = BnapusUnigene))
  
  print("Producing final results")
  
  finalResults <- cbind(traitAndUnigene, log10P, lmResults)
  finalResults = getGeneModels(finalResults)
  finalResults <- finalResults[order(finalResults$log10P, decreasing = TRUE),]
  
  ## Write these columns for now. In future if the apply() is quick
  ## enough probably combine this script and grapher one together.
  finalResults = finalResults[,c("trait", "unigene",  "A.Chr", "C.Chr", "AGI", 
                                 "log10P", "Df", "SumSq", "MeanSq", "Fvalue", "Pvalue")]
  resultsFile = paste("GEM_results_", traitname, collapse="-", ".txt", sep="")
  write.table(finalResults, resultsFile, quote=FALSE, sep="," ,row.names=FALSE, col.names=TRUE)
  ## Example output in console
  ##     trait      unigene      A.Chr C.Chr    AGI     log10P   Df  SumSq
  ##  minusefvern A_JCVI_14001    A3    C2 AT5G59050.1 14.28305  1 6632.848
  ##  minusefvern A_JCVI_40108    A2    C2 AT1G65480.1 13.99106  1 6532.559
  ##    MeanSq   Fvalue     Pvalue
  ##  6632.848 88.62775 5.211345e-15
  ##  6532.559 85.99291 1.020799e-14

  
  print("Setting up the GEM Manhattan plot data")
  
  GEMinput <- fread(paste("GEM_results_", traitname, ".txt", sep = ""))
  setnames( GEMinput, old = 'unigene', new = 'Marker' )
  setnames( GEMinput, old = 'trait', new = 'Trait' )
  
  GEMinput = GEMinput[ !is.nan( log10P ), ]
  
  GEM.dt = GEMinput[, c('Trait', 'Marker', 'log10P')]
  
  annos.dt.GEM = fread( paste("../../", dir(pattern = "mapping_", "../../"), sep=""))
  
  setnames( annos.dt.GEM, old = 'Unigene', new = 'Marker' )
  annos.dt.GEM = annos.dt.GEM[, c('Marker', 'Chr', 'AGI')]
  annos.dt.GEM = annos.dt.GEM %>% drop_na(Marker)
  annos.dt.GEM = annos.dt.GEM[ !duplicated(annos.dt.GEM), ]
  
  graph.dt.GEM = fread( paste("../../", dir(pattern = "directions_", "../../"), sep=""))
  
  setnames( graph.dt.GEM, old = 'sort.GEM', new = 'sort' )
  setnames( graph.dt.GEM, old = 'Unigene', new = 'Marker' )
  graph.dt.GEM = graph.dt.GEM[, c('Marker', 'Chr', 'sort')]
  graph.dt.GEM = graph.dt.GEM %>% drop_na(Marker, sort)
  graph.dt.GEM = graph.dt.GEM[ !duplicated(graph.dt.GEM), ]
  
  meta.dt.GEM = merge( graph.dt.GEM, annos.dt.GEM, by = c('Marker', 'Chr'))
  
  min_vals_GE = setNames( aggregate( sort ~ Chr, data = meta.dt.GEM, FUN = min ), 
                          nm = c( 'Chr', 'Start' ) )
  max_vals_GE = setNames( aggregate( sort ~ Chr, data = meta.dt.GEM, FUN = max ), 
                          nm = c( 'Chr', 'End' ) )
  
  chromo.dt.GEM = data.table( merge( min_vals_GE, max_vals_GE, by = 'Chr' ) )
  chromo.dt.GEM = separate( data = chromo.dt.GEM,
                            col = 'Chr',
                            into = c( 'Genome', 'ChromoNum' ),
                            sep = c( 1 ),
                            remove = FALSE,
                            convert = TRUE )
  
  
  meta.dt.GEM = merge( meta.dt.GEM, chromo.dt.GEM, by = 'Chr' )
  
  allGEM.dt = merge( GEM.dt, meta.dt.GEM, by = 'Marker', na.rm = TRUE )
  setorder( allGEM.dt, sort )
  
  metaminus <- gsub(":.*", "", meta.dt.GEM$Marker)
  meta.dt.GEM$Marker <- metaminus
  allGEM.dt = merge( GEM.dt, meta.dt.GEM, by = 'Marker', na.rm = TRUE )
  setorder( allGEM.dt, sort )
  
  
  print("Analysing GEM data Q values and Bonferroni")
  
  AllPhenoData <- read.delim(paste("GEM_results_", traitname, ".txt", sep = ""), sep = ",")
  AllPhenoData <- AllPhenoData[order(AllPhenoData$Pvalue),]
  
  pvalues <- AllPhenoData$Pvalue
  qobj <- qvalue(p = pvalues)
  
  qvalues <- qobj$qvalues
  pi0 <- qobj$pi0
  lfdr <- qobj$lfdr
  
  sigSNP <- min(AllPhenoData$unigene[which(lfdr < 0.05)])
  FDRLine <- allGEM.dt$log10P[which(allGEM.dt$Marker == noquote(sigSNP))]
  
  sigSNP2 <- min(AllPhenoData$unigene[which(qvalues < 0.05)])
  QLine <- allGEM.dt$log10P[which(allGEM.dt$Marker == noquote(sigSNP2))]
  
  results <- data.frame(AllPhenoData$unigene, AllPhenoData$A.Chr, 
                        AllPhenoData$C.Chr, AllPhenoData$Pvalue,qvalues, lfdr, pi0)
  colnames(results) <- c("Marker", "A-Chromosome", "C-Chromosome", "P.value",
                         "Q.value", "FDR", "Pi0")
  
  print("Printing results")
  write.csv(results, paste("GEM_fitting_results_", traitname, ".csv", sep = ""))
  
  hist <- hist(qobj)
  png(paste("GEM_fitting_histogram_", traitname, ".png", sep = ""))
  print(hist)
  dev.off()
  
  png(paste("GEM_fitting_qvalues_", traitname, ".png", sep = ""))
  plot(qobj)
  dev.off()
  
  print("Printing the GEM Manhattan plot")
  
  ylim = c(0, max( allGEM.dt$log10P, na.rm = T ) )
  
  # choose A and/or C, if present
  genome_labels = c()
  for (genome_label in c('A', 'C'))
    if ( genome_label %in% allGEM.dt$Genome )
      genome_labels = c( genome_labels, genome_label )
  
  jpeg(
    filename = paste("GEM_Manhattan_", traitname, ".jpeg", sep = ""),
    quality = 1000,
    width = 3000,
    height = 500 * length( genome_labels ),
    units = 'px',
    bg = 'white',
    res = NA
  )
  
  # Sets up a split screen - required when data is for A and C genomes
  par(mfrow = c( length(genome_labels), 1 ) )
  
  for (genome_label in genome_labels) {
    genomeGEM.dt = allGEM.dt[ Genome == genome_label, c('Marker', 'sort', 
                                                        'log10P', 'ChromoNum', 'End')]
    genomeGEM.dt = genomeGEM.dt[ !duplicated(genomeGEM.dt), ]
    
    xlim = c(1, max(genomeGEM.dt$End))
    
    color_pallete_function <- colorRampPalette(
      colors = palette("Tableau"),
      space = "Lab" # Option used when colors do not represent a quantitative scale
    )
    
    num_colors <- length(unique(genomeGEM.dt$ChromoNum))
    Chrom_colors <- color_pallete_function(num_colors)
    
    plot (
      genomeGEM.dt$sort,
      genomeGEM.dt$log10P,
      type = 'p',
      pch = 20,
      xlim = xlim,
      ylim = ylim,
      main = genome_label,
      cex.main = 2,
      ylab = '-log10(p)',
      xlab = ' ',
      cex.lab = 1.5,
      col = Chrom_colors[genomeGEM.dt$ChromoNum],
      xaxt = 'n'
      
    )
    
    abline(h=FDRLine, col="red", lty = 2)
    abline(h=QLine, col="blue", lty = 2)
    
  }
  
  close.screen(all = TRUE)
  dev.off()
  
  print("GEM analysis complete!")
  
}

################################################################################
## GWAS ANALYSIS
## GWAS analysis run in GAPIT. Outputs into folder called "GAPIT".
################################################################################

GWASanalysis <- function(gwasmodels = gwasmodels){
  print("Starting GWAS analysis...")
  
  # Phenotype data
  myY <- read.table("LMMestimatemeans.txt", header = TRUE)
  colnames(myY)[1] <- "Taxa"
  
  dir.create("GAPIToutputs")
  setwd("GAPIToutputs") 

  #Best number of PCA for GWAS
  myGAPIT <- GAPIT(
    Y=myY,
    G=myG,
    CV=myCV,
    PCA.total=0,
    SNP.MAF = 0.05,
    model = 'GLM',
    Geno.View.output = FALSE, 
    Model.selection = TRUE
    )

  # Sets the best number of PCA to be used in the analysis
  PCAresults <- read.csv(paste("GAPIT.GLM.", colnames(myY)[2],
                               ".BIC.Model.Selection.Results.csv", sep = ""))
  BestPCA <- max(PCAresults[,2])
  PCAno <- as.numeric(PCAresults[
    PCAresults$BIC..larger.is.better....Schwarz.1978 == BestPCA,
                      "Number.of.PCs.Covariates"])

  # Run GAPIT across multiple models
  myGAPIT <- GAPIT(
    Y=myY,
    G=myG,
    CV=myCV,
    PCA.total = PCAno,
    SNP.MAF = 0.05,
    model = gwasmodels,
    Geno.View.output = FALSE,
    Multiple_analysis = T
  )
  
  ## Fit Test ##
  
  print("Running fit test and finding significant markers")
  
  BestFit <- data.frame(models = gwasmodels, 
                        FitScore = NA) 
  
  GAPITsig <- data.frame(Model = NA, Marker = NA, P.value = NA, FDR = NA)
  
  for (model in 1:length(gwasmodels)) {
    Testing <- read.csv(paste("GAPIT.Association.GWAS_Results.", BestFit[model,1], ".", colnames(myY)[2], 
                              ".csv", sep = ""))
    # Testing <- filter(Testing, !is.nan(P.value))
    Pvalues <- Testing$P.value
    n <- length(Pvalues)
    Actual_quantiles <- sort(Pvalues)
    Expected_quantiles = seq_along(Pvalues)/n 
    ActualExpected <- data.frame(expected = Expected_quantiles, 
                                 actual = Actual_quantiles)
    
    fitofmodel <- 0
    fitofmodel <- sum((ActualExpected$expected - ActualExpected$actual)^2)
    BestFit$FitScore[model] <- fitofmodel
    
    print(BestFit)
  
    subpvalues <- Testing$P.value
    subqvalues <- qvalue(p = subpvalues)
    fdrvalues <- subqvalues$lfdr
    Testing['FDR'] <- fdrvalues 
    subsetres <- subset(Testing, Testing$FDR <= 0.05)
    
    if (length(subsetres$SNP) >=1) {
      sigonly <- data.frame(Model = BestFit[model,1], Marker = subsetres$SNP, 
                          P.value = subsetres$P.value, FDR = subsetres$FDR)
      GAPITsig <- rbind(GAPITsig, sigonly)
    }
  }
  
  GAPITsig <- GAPITsig[2:length(GAPITsig$Model),]
  GAPITsig <- na.omit(GAPITsig)
  write.csv(GAPITsig, "../GAPIT_Significant_Markers.csv", row.names = F)
  
  print("Saving the best fit model results...")
  
  BestFitModel <- BestFit[which.min(BestFit$FitScore),]
  BestFitModel <- BestFitModel[1,1]
  write.csv(BestFit, "../Bestfit.Model.Analysis.csv", row.names = F)
  
  file.copy(from =  paste("GAPIT.Association.GWAS_Results.", BestFitModel[1], ".", colnames(myY)[2], 
                          ".csv", sep = ""),
            to = paste("../GAPIT.", BestFitModel,".GWAS.Results.csv", 
                       sep = ""))
  
  file.copy(from = dir(pattern = "Circular"), 
            to = "../GAPIT.Circular.Manhattan.pdf")
  
  setwd("../")
  
  print("Setting up the GAPIT Manhattan plot data")
  
  traitname <- colnames(myY)[2]
  gapitinput <- fread(paste("GAPIT.", BestFitModel ,".GWAS.Results.csv", sep = ""))
  setnames( gapitinput, old = 'SNP', new = 'Marker' )
  gapitinput <- cbind(Trait = traitname, gapitinput)
  
  gapitinput = gapitinput[ !is.nan( P.value ), ]
  gapitinput$log10P = -log10(gapitinput$P.value)
  
  gapit.dt = gapitinput[, c('Trait', 'Marker', 'log10P')]
  
  annos.dt.GAPIT = fread( paste("../../", dir(pattern = "mapping_", "../../"), sep=""))
  graph.dt.GAPIT = fread( paste("../../", dir(pattern = "directions_", "../../"), sep=""))
  
  setnames( annos.dt.GAPIT, old = 'Unigene.Marker', new = 'Marker' )
  annos.dt.GAPIT = annos.dt.GAPIT[, c('Marker', 'Chr', 'AGI')]
  annos.dt.GAPIT = annos.dt.GAPIT %>% drop_na(Marker)
  annos.dt.GAPIT = annos.dt.GAPIT[ !duplicated(annos.dt.GAPIT), ]
  
  setnames( graph.dt.GAPIT, old = 'sort.SNP', new = 'sort' )
  setnames( graph.dt.GAPIT, old = 'Unigene.Marker', new = 'Marker' )
  graph.dt.GAPIT = graph.dt.GAPIT[, c('Marker', 'Chr', 'sort')]
  graph.dt.GAPIT = graph.dt.GAPIT %>% drop_na(Marker, sort)
  graph.dt.GAPIT = graph.dt.GAPIT[ !duplicated(graph.dt.GAPIT), ]
  
  meta.dt.GAPIT = merge( graph.dt.GAPIT, annos.dt.GAPIT, by = c('Marker', 'Chr'))
  
  min_vals_GA = setNames( aggregate( sort ~ Chr, data = meta.dt.GAPIT, FUN = min ), 
                          nm = c( 'Chr', 'Start' ) )
  max_vals_GA = setNames( aggregate( sort ~ Chr, data = meta.dt.GAPIT, FUN = max ), 
                          nm = c( 'Chr', 'End' ) )
  
  chromo.dt.GAPIT = data.table( merge( min_vals_GA, max_vals_GA, by = 'Chr' ) )
  chromo.dt.GAPIT = separate( data = chromo.dt.GAPIT,
                              col = 'Chr',
                              into = c( 'Genome', 'ChromoNum' ),
                              sep = c( 1 ),
                              remove = FALSE,
                              convert = TRUE )
  
  
  meta.dt.GAPIT = merge( meta.dt.GAPIT, chromo.dt.GAPIT, by = 'Chr' )
  
  allGAPIT.dt = merge( gapit.dt, meta.dt.GAPIT, by = 'Marker', na.rm = TRUE )
  setorder( allGAPIT.dt, sort )
  
  print("Analysing GAPIT data Q values and Bonferroni")
  
  AllPhenoData <- read.csv(paste("GAPIT.", BestFitModel,".GWAS.Results.csv", sep = ""))
  AllPhenoData <- AllPhenoData[order(AllPhenoData$P.value),]
  
  pvalues <- AllPhenoData$P.value
  qobj <- qvalue(p = pvalues)
  
  qvalues <- qobj$qvalues
  pi0 <- qobj$pi0
  lfdr <- qobj$lfdr
  
  sigSNP <- min(AllPhenoData$SNP[which(lfdr < 0.05)])
  FDRLine <- allGAPIT.dt$log10P[which(allGAPIT.dt$Marker == noquote(sigSNP))]
  
  sigSNP2 <- min(AllPhenoData$SNP[which(qvalues < 0.05)])
  QLine <- allGAPIT.dt$log10P[which(allGAPIT.dt$Marker == noquote(sigSNP2))]
  
  results <- data.frame(AllPhenoData$SNP, qvalues, lfdr, pi0)
  
  print("Printing results")
  
  write.csv(results, paste("GAPIT_fitting_", traitname, ".csv", sep = ""))
  
  hist <- hist(qobj)
  png(paste("GAPIT_fitting_histogram_", traitname, ".png", sep = ""))
  print(hist)
  dev.off()
  
  png(paste("GAPIT_fitting_qplot_", traitname, ".png", sep = ""))
  plot(qobj)
  dev.off()
  
  print("Printing the Manhattan plot")
  
  ylim = c(0, max( allGAPIT.dt$log10P ) )
  
  # choose A and/or C, if present
  genome_labels = c()
  for (genome_label in c('A', 'C'))
    if ( genome_label %in% allGAPIT.dt$Genome )
      genome_labels = c( genome_labels, genome_label )
  
  # Bonf <- p.adjust(gapitinput$P.value, method = "bonferroni")
  # sigSNP <- min(gapitinput$Marker[which(Bonf < 0.05)])
  # bonfLine <- allGAPIT.dt$log10P[which(allGAPIT.dt$Marker == noquote(sigSNP))]
  # 
  # Pvals <- gapitinput$Pvalue
  # sigSNP2 <- min(gapitinput$Marker[which(Pvals < 0.05)])
  # PvalLine <- allGAPIT.dt$log10P[which(allGAPIT.dt$Marker == noquote(sigSNP2))]
  
  jpeg(
    filename =  paste("GAPIT_Manhattan_", traitname, "_", BestFitModel, ".jpeg", sep = ""),
    quality = 1000,
    width = 3000,
    height = 500 * length( genome_labels ),
    units = 'px',
    bg = 'white',
    res = NA
  )
  
  # Sets up a split screen - required when data is for A and C genomes
  par(mfrow = c( length(genome_labels), 1 ) )
  
  for (genome_label in genome_labels) {
    genomeGAPIT.dt = allGAPIT.dt[ Genome == genome_label, c('Marker', 'sort', 
                                                            'log10P', 'ChromoNum', 'End')]
    genomeGAPIT.dt = genomeGAPIT.dt[ !duplicated(genomeGAPIT.dt), ]
    
    xlim = c(1, max(genomeGAPIT.dt$End))
    
    color_pallete_function <- colorRampPalette(
      colors = palette("Tableau"),
      space = "Lab" # Option used when colors do not represent a quantitative scale
    )
    
    num_colors <- length(unique(genomeGAPIT.dt$ChromoNum))
    Chrom_colors <- color_pallete_function(num_colors)
    
    plot (
      genomeGAPIT.dt$sort,
      genomeGAPIT.dt$log10P,
      type = 'p',
      pch = 20,
      xlim = xlim,
      ylim = ylim,
      main = genome_label,
      cex.main = 2,
      ylab = '-log10(p)',
      xlab = ' ',
      cex.lab = 1.5,
      col = Chrom_colors[genomeGAPIT.dt$ChromoNum],
      xaxt = 'n'
      
    )
    
    abline(h=FDRLine, col="red", lty = 2)
    abline(h=QLine, col="blue", lty = 2)
    
  }
  
  close.screen(all = TRUE)
  dev.off()
  
  print("GWAS analysis complete!")


}

################################################################################
## Running the pipeline.
################################################################################

print("Running pipeline...")

setwd(paste(masterfolder, sep = ""))

for (trait in 1:length(trait.list)) {
  
  print(paste("GEM analysis for ", trait.list[trait], " which is number ", trait, 
              " of ", length(trait.list), " traits.", sep = ""))
  
  setwd(paste(trait.list[trait], sep = ""))
  
  if (rungem == TRUE) {
  GEManalysis()
  }
  setwd("../")
}

if (rungem == TRUE) {
  print("Summarising GEM analysis results...")
  
  GEMsig <- data.frame(Trait = NA, Number_of_significant_markers_by_P_value = NA)
  
  for (sigfind2 in 1:length(trait.list)) {
    
    sigfinder2 <- read.csv(paste(trait.list[sigfind2], "/GEM_results_",trait.list[sigfind2], ".txt", sep=""))
    
    if(length(sigfinder2$Pvalue<=0.05) > 0) {
      
      noofsigs <- data.frame(Trait = trait.list[sigfind2],
                      Number_of_significant_markers_by_P_value = length(sigfinder2$Pvalue<=0.05)) 
      
      GEMsig <- rbind(GEMsig, noofsigs)
      
    } else {
      
    }
    
  }
  
  GEMsig <- na.omit(GEMsig)
  
  write.csv(GEMsig, "GEM_traits_with_significant_markers.csv", row.names = F)
  print("GEM summary complete!")
}


for (trait in 1:length(trait.list)) {
  
  print(paste("GWAS analysis for ", trait.list[trait], " which is number ", trait, 
              " of ", length(trait.list), " traits.", sep = ""))
  
  setwd(paste(trait.list[trait], sep = ""))
  
  if (rungwas == TRUE) {
    GWASanalysis(gwasmodels = gwasmodels)
  }
  
  setwd("../")
  
}

if (rungwas == TRUE) {
  print("Summarising GWAS analysis results...")
  
  GAPITsig <- data.frame(Trait = NA, Model = NA, Marker = NA, P.value = NA, FDR = NA)
  
  for (sigfind in 1:length(trait.list)) {
    
    sigfinder <- read.csv(paste(trait.list[sigfind], "/GAPIT_Significant_Markers.csv", sep=""))
    
    if(length(sigfinder$Model) <= 100) {
    
    colone <- data.frame(Trait = rep(trait.list[sigfind], length(sigfinder$Model)))
    sigfinder <- cbind(colone, sigfinder)
    
    GAPITsig <- rbind(GAPITsig, sigfinder)
    
    } else {
      
    }
    
  }
  
  GAPITsig <- na.omit(GAPITsig)
  
  write.csv(GAPITsig, "GWAS_traits_with_significant_markers.csv", row.names = F)
  
  print("GWAS summary complete!")
  
}

print("Pipeline complete!")
