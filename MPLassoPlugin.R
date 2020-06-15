#####################################################################
# Synthetic graph generation: random, AR4, hub, cluster, scale-free
#
# @author Chieh Lo
# @date 01/26/2017
#####################################################################


graph_select <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15, type = "random"){
  if(type=="random"){
    sample <- random_model(OTUnum, Simulations, Strength);
    #print("random");
  }
  else if(type == "AR4"){
    sample <- AR4_model(OTUnum, Simulations, Strength);
    #print("AR4");
  }
  else if(type == "hub"){
    sample <- hub_model(OTUnum, Simulations, Strength);
    #print("hub");
  }
  else if(type == "cluster"){
    sample <- cluster_model(OTUnum, Simulations, Strength);
    #print("cluster");
  }
  else if(type=="scale-free"){
    sample <- scale_free_model(OTUnum, Simulations, Strength);
    #print("scale-free");
  }
  return(sample);
}


AR4_model <- function(OTUnum = 50, Simulations = 1e2, MaxStrength = 0.15 ){
  sim <- huge.generator(n = Simulations, d = OTUnum, graph = "band", v = 0.3, u = 0 , g= 4, verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}

hub_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.2 ){
  sim <- huge.generator(n = Simulations, d = OTUnum, graph = "hub", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));  
}

cluster_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15 ){
  sim <- huge.generator(n = Simulations, d = OTUnum, graph = "cluster", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}

random_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15 ){
  sim <- huge.generator(n = Simulations, d = OTUnum, graph = "random", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}

scale_free_model <- function(OTUnum = 50, Simulations = 1e2, Strength = 0.15 ){
  sim <- huge.generator(n = Simulations, d = OTUnum, graph = "scale-free", v = 0.3, u = 0 , verbose = FALSE)
  x <- round(exp(sim$data));
  data <- t(clr(x+0.1, 1));
  s <- var(data);
  R <- sim$sigma;
  A <- as.matrix(sim$theta);
  diag(A) <- 1;
  return(list("sample" = x, "fraction" = data , "var" = s, "cor" = R, "adj" = A));
}


#####################################################################
# Graphical Lasso optimal lambda selection
#
# @author Chieh Lo
# @date 01/26/2017
#####################################################################

selection <- function(sample, type = "MPLasso", prior = 1, priorType = TRUE, interactionFlag = FALSE, precision=FALSE, precisionRatio = 0.5, priorInfo){
  p <- ncol(sample$sample);
  n <- nrow(sample$sample);
  
  lambdaMin = 0.01;
  lambdaMax = 1;
  interval = 0.01;
  numLambda = (lambdaMax - lambdaMin)/interval + 1;
  BICResult = matrix(0, numLambda, 1);
  
  for (i in 1:numLambda){
    result <- MPLasso(sample, lambdaMin = (lambdaMin+(i-1)*interval), lambdaMax = 10, prior = prior, priorType = priorType, interactionFlag = interactionFlag, precision=precision, precisionRatio = precisionRatio, priorInfo = priorInfo);
    BICResult[i] <- BIC(sample, result);
  } 
  return(BICResult);
}

BIC <- function(sample, result){
  p <- ncol(sample$sample);
  n <- nrow(sample$sample); 
  gamma <- 0.5;
  likelihood <- n/2*(log(det(result$icov)) - sum(diag(sample$var %*% result$icov)) );
  edge <- (sum(result$adj) - p)/2*log(n);
  return(-2*likelihood + edge );
}#####################################################################
# PLasso implementation 
# You can implement your algorithm here
# @author Chieh Lo
# @date 05/13/2017
#####################################################################

algorithm_select <- function(sample, lambdaMin = 0.01, lambdaMax = 10, prior = 1, type = "MPLasso", priorType = TRUE, interactionFlag = FALSE, precision = FALSE, precisionRatio, priorInfo = priorInfo){
  if(type=="MPLasso"){
    result <- MPLasso(sample, lambdaMin, lambdaMax, prior = prior, priorType = priorType, interactionFlag = interactionFlag, precision = precision, precisionRatio = precisionRatio, priorInfo = priorInfo)
  }
  return(result);
}


MPLasso <- function(sample, lambdaMin = 0.01, lambdaMax = 10, prior = 1, priorType, interactionFlag, precision = FALSE, precisionRatio = precisionRatio, priorInfo = priorInfo){
  s <- sample$var;
  adj <- sample$adj;
  if(interactionFlag==TRUE){
    interaction <- sample$inter;
  }
  p <- nrow(s);
  rho <- lambdaMin*matrix(1, p, p); #gLasso regularizer

  ## no pricision ratio (i.e., all priors are correct)
  ## adj encodes true graph structure. For synthetic data, we use the true graph structure as our prior info. 
  ## For real data, adj encodes the co-occurrence information (i.e., taxa is not associated with each other)
  if (priorType == TRUE && precision == FALSE){
    for (i in 1:p){
      for (j in 1:p){
        if (abs(adj[i,j]) == 0 && runif(1, min=0, max=1) > prior && i>j){
          rho[i,j] <- lambdaMax;
          rho[j,i] <- lambdaMax;
        }
      }
    }
  }

  ## consider pricision ratio (given prior information)
  if (priorType == TRUE && precision == TRUE && prior < 1){
    rhoInter <- priorInfo$rhoInter
    rhoOccur <- priorInfo$rhoOccur
    rho <- lambdaMin*matrix(1, p, p); #gLasso regularizer
    if (precisionRatio==0){
      lambdaMax = 1
    }
    else{
      lambdaMax = lambdaMin/precisionRatio 
    }

    for (i in 1:p){
      for (j in 1:p){
        if(rhoOccur[i,j]==1 && rhoInter[i,j]==1){
          rho[i,j] <- lambdaMax
          rho[j,i] <- lambdaMax
        }
      }
    }
  }
  if(prior == 1){
    rho <- lambda_min*matrix(1, p, p); #gLasso regularizer
  }




  if(interactionFlag == TRUE){
    for (i in 1:p){
      for (j in 1:p){
        if (abs(interaction[i,j]) == 1 && runif(1, min=0, max=1) > prior && i!=j){
          rho[i,j] <- 0.01;
          rho[j,i] <- 0.01;
        }
      }
    }
  }
  a<-glasso(s, rho=rho)
  covOpt <- a$w;
  icovOpt <- a$wi;
  d <- 0;
  for(i in 1:p){
    d[i] <- 1/sqrt(covOpt[i,i]);
  }
  D <- diag(d);
  corrOpt <- D%*%covOpt%*%D;
  
  adjOpt <- matrix(0, p, p);
  for(i in 1:p){
    for (j in i:p){
      if (abs(icovOpt[i,j]) > 0){
        adjOpt[i,j] <- 1;
        adjOpt[j,i] <- 1;
      }
    }
  }
  
  
  return(list("cor" = corrOpt, "cov" = covOpt, "icov" = icovOpt, "adj" = adjOpt ));
}



#####################################################################
# Performance evaluation
#
# @author Chieh Lo
# @date 01/26/2017
#####################################################################

evaluation <- function(sample, result){
  p <- ncol(sample$sample);
  n <- nrow(sample$sample);

  RMSE <- RMSE_l0(sample$cor, result$cor);
  ICOV <- ICOV_l0(ginv(sample$cor), result$icov);
  RMSE_f <- RMSE_F(sample$cor, result$cor);
  ICOV_f <- ICOV_F(ginv(sample$cor), result$icov); 

  return(list("RMSE_l0" = RMSE, "ICOV_l0" = ICOV, "RMSE_F" = RMSE_f, "ICOV_F" = ICOV_f));
}

RMSE_l0 <- function(true_cor, estimate_cor){
  p <- nrow(true_cor);
  Error <- abs(true_cor - estimate_cor);
  RMSE <- sum(sum(Error))/(p*(p-1));
  return(RMSE);
}

ICOV_l0 <- function(true_icov, estimate_icov){
  p <- nrow(true_icov);
  Error <- abs(true_icov - estimate_icov);
  RMSE <- sum(sum(Error))/(p*(p-1));
  return(RMSE);
}


RMSE_F <- function(true_cor, estimate_cor){
  p <- nrow(true_cor);
  Error <- abs(true_cor - estimate_cor)*abs(true_cor - estimate_cor);
  RMSE <- sqrt(sum(sum(Error)));
  return(RMSE);
}

ICOV_F <- function(true_icov, estimate_icov){
  p <- nrow(true_icov);
  Error <- abs(true_icov - estimate_icov)*abs(true_icov - estimate_icov);
  RMSE <- sqrt(sum(sum(Error)));
  return(RMSE);
}


ROC_eval <- function(label, pred){
  
  pred <- prediction(pred, label)
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred, 'auc')
  aucstd <- sd(simplify2array(auc@y.values))
  auc <- mean(simplify2array(auc@y.values))
  return(list("perf" = perf, "auc" = auc, "aucsd" = aucstd));
}

ROC_new <- function(label, pred){
  
  pred <- prediction(pred, label)
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred, 'auc')
  aucstd <- sd(simplify2array(auc@y.values))
  auc <- mean(simplify2array(auc@y.values))
  f <- performance(pred,"prec","rec")
  x <- f@x.values
  y <- f@y.values

  AUCPr = matrix(0,length(x), 1);
  for (i in 1:length(x)){
    xTemp <- as.matrix(x[[i]])
    yTemp <- as.matrix(y[[i]])
    yTemp[1] <- yTemp[2];
    AUCPr[i] = trapz(xTemp, yTemp)
    f@y.values[[i]][1] <- f@y.values[[i]][2];
  }
  aucpr <- mean(AUCPr)
  aucprsd <- sd(AUCPr)

  return(list("perf" = perf, "auc" = auc, "aucsd" = aucstd, "f" = f, "aucpr" = aucpr, "aucprsd" = aucprsd));
}


acc_eval <- function(label, pred){
  pred <- prediction(pred, label)
  perf <- performance(pred,"acc")
  acc <- matrix(0, length(label)) 
  for (i in 1:length(label)){
    index <- which.max(perf@y.values[[i]])
    acc[i] <- perf@y.values[[i]][index];
  }
  return(acc);
}

repro_eval <- function(true_adj, test_adj){
  row_sum <- rowSums(true_adj);
  degree <- sort(row_sum, index.return = TRUE, decreasing=TRUE);
  degree_index <- degree$ix;
  num_taxa <- nrow(true_adj);
  percent = c(0.25, 0.5, 0.75, 1);
  reproError = matrix(0, 1, length(percent))
  for(i in 1:length(percent)){
    temp_index <- degree_index[1:round(num_taxa*percent[i])]
    #print(temp_index)
    num_edge_true <- sum(true_adj[temp_index,]) - round(num_taxa*percent[i]) 
    num_edge_test <- sum(test_adj[temp_index,]) - round(num_taxa*percent[i])
    error <- true_adj[temp_index,] -test_adj[temp_index,];
    #print(nrow(error))
    reproError[i] <- sum(abs(error))/round((num_taxa*num_taxa*percent[i]))
    #print(sum(abs(error)))
    #print(round(num_taxa*percent[i]))
  }

  return(reproError);
  
}
#####################################################################
# Calculate the prior information
#
# @author Chieh Lo
# @date 01/26/2017
#####################################################################

occurance <- function(occuPath, totalPaper = 615268){
  a <- read.csv(occuPath, sep = "\t", header=FALSE)
  countTable <- as.matrix(a)
  n = nrow(countTable);
  pValue = matrix(0, n, n);
  CILow = matrix(0, n, n);
  CIHigh = matrix(0, n, n);
  totalPaper = sum(diag(countTable))
  for (i in 1:n){
    for (j in 1:n){
      if (j>i){
        contingency <- matrix(c(countTable[i,j], countTable[j,j] , countTable[i, i], totalPaper), nrow = 2)
        temp <- fisher.test(contingency, alternative = "two.sided");
        pValue[i,j] <- temp$p.value;
        CILow[i,j] <- temp$conf.int[[1]];
        CIHigh[i,j] <- temp$conf.int[[2]];
      }
    }
  }
  
  noAssociation <- matrix(0,n,n);
  adjust <- round(p.adjust(pValue, method = "bonferroni"), 5)
  corrected <- matrix(adjust, ncol=n, nrow = n);
  for (i in 1:n){
    for (j in 1:n){
      if (j>i && (CIHigh[i,j] > 1 && CILow[i,j] < 1) && corrected[i,j] == 1){
        noAssociation[i,j] <- 1;
        noAssociation[j,i] <- 1;
      }
    }
  }
  return(noAssociation)
}

interaction <- function(interactionPath){
  a <- read.csv(interactionPath, sep = "\t", header=FALSE)
  interaction <- as.matrix(a)
  return(interaction)
}
#####################################################################
# Real data preprocessing
#
# @author Chieh Lo
# @date 01/26/2017
#####################################################################

read_OTU <- function(level = "genus", countPath, association = NULL, interaction = NULL){
  if(level == "genus"){
    rawData <- read.csv(countPath, sep = "\t", header=FALSE)
    rawData <- as.matrix(rawData)
    sample = count_data_process(t(rawData), association);
  }
  else if(level == "species"){
    rawData <- t(as.matrix(read.table(countPath)));
    sample = fraction_data_process(rawData, association, interaction);
  }
  return(sample);
}

count_data_process <- function(rawData, association){
  p <- ncol(rawData);
  data <- t(clr(rawData+0.1, 1));
  s <- var(data);
  d <- 0;
  for(i in 1:p){
    d[i] <- 1/sqrt(s[i,i]);
  }
  D <- diag(d);
  R <- D%*%s%*%D;
  A <- association;
  return(list("sample" = rawData, "fraction" = data , "var" = s, "cor" = R, "adj" = association));
}

fraction_data_process <- function(rawData, association, interaction){
  p <- ncol(rawData);
  data <- t(clr(rawData+0.001, 1));
  s <- var(data);
  d <- 0;
  for(i in 1:p){
    d[i] <- 1/sqrt(s[i,i]);
  }
  D <- diag(d);
  R <- D%*%s%*%D;
  A <- association;
  return(list("sample" = rawData, "fraction" = data , "var" = s, "cor" = R, "adj" = association, "inter" = interaction));
}
#####################################################################
# Log-ratio transform
#
# Adopted from Zachary D. Kurtz et. al., Sparse and Compositionally Robust Inference of Microbial Ecological Networks, Plos One, 2015, http://dx.doi.org/10.1371/journal.pcbi.1004226
#####################################################################

#' Centered log-ratio functions
#' @export
clr <- function(x, ...) {
    UseMethod('clr', x)
}

#' @method clr default
#' @export
clr.default <- function(x.f, base=exp(1), tol=.Machine$double.eps) {
    nzero <- (x.f >= tol)
    LOG <- log(ifelse(nzero, x.f, 1), base)
    ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

#' @method clr matrix
#' @export
clr.matrix <- function(x.f, mar=2, ...) {
    apply(x.f, mar, clr, ...)
}
#####################################################################
# PLasso implementation 
# You can implement your algorithm here
# @author Chieh Lo
# @date 05/13/2017
#####################################################################
oldlog <- log
dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")
log <- oldlog

# library required
library(MASS)
library(Matrix)
library(matrixcalc)
library(glasso)
library(ROCR)
library(huge)
library(MCMCpack)
library(glmnet)
library(parallel)
library(matrixStats)


#source files required
#scriptDir <- "./";#dirname(sys.frame(1)$ofile)
#source(paste(scriptDir, '/R/syn_data_generate.R', sep = ''));
#source(paste(scriptDir, '/R/selection.R', sep = ''));
#source(paste(scriptDir, '/R/algorithm.R', sep = ''));
#source(paste(scriptDir, '/R/evaluation.R', sep = ''));
#source(paste(scriptDir, '/R/calculate_co_occurrence.R', sep = ''))
#source(paste(scriptDir, '/R/real_data.R', sep = ''))
#source(paste(scriptDir, '/R/normalization.R', sep = ''))



# five body sites considered in this demo
#bodySites <- c('Anterior_nares', 'Buccal_mucosa', 'Stool', 'Supragingival_plaque', 'Tongue_dorsum')
#bodySites <- c('Anterior_nares')
# flag for evaluating different datasets

input <- function(inputfile) {
# reproducibility evaluation
repro <<- TRUE;

parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
  #print("READING INPUT FILES...");
  #print(interactionFile);
  #for (i in 1:length(bodySites)){
    #print("HMASM:")
    #print(bodySites[i]);
    #bodySite = bodySites[i];
    #prePath = paste(scriptDir, "/data/HMASM/", sep = '');
    prePath = prefix();
    #interactionFile = "/example/species_interaction.csv";
    #print("PREFIX")
    #print(prefix())
  interactionFile <<- parameters["interactionFile",2]
  #print(interactionFile)
    interactionPath = paste(prePath, interactionFile, sep= '/'); 
    #print(interactionPath)
    interactionMatrix <<- interaction(interactionPath); # read interaction file 
    #print(interactionMatrix);
    #print("DONE INTER MATRIX")  
  occurrenceFile <<- parameters["occurrenceFile",2];
  #print(occurrenceFile)
    occuPath = paste(prePath, occurrenceFile, sep = '/')
  #print(occuPath)
    #occuPath = paste(prePath, "/example/association.csv", sep = '');
    NoAssociation <<- occurance(occuPath); # read occurance file
    #print(NoAssociation)
    #print("DONE ASSOC")
    countFile <<- toString(parameters["countFile",2])
    countPath = paste(prePath, countFile, sep = '/');
    #countPath = paste(prePath, "/example/SpeciesCount.txt", sep = '');
    samples <<- read_OTU(level = "species", countPath = countPath, 1-NoAssociation, interactionMatrix); # read OTU
    #print(samples)
    #print("DONE SAMPLES")
 }

run <- function() {
    prior <<- 0;
    interactionFlag <<- TRUE;
    #print(samples);
    selectLambda <<- selection(samples, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
    lambdaOptimal <<- (which.min(selectLambda)-1)*0.01 + 0.01;
    resultPL <<- algorithm_select(samples, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
    #print(resultPL);
}

output <- function(outputfile) {
    write.csv(resultPL$cor, outputfile)#"output.csv")
    nInt <- 2;
    reproduceErrorArrayPL = matrix(0,nInt, 4); 
    if (repro == TRUE){
      for (k in 1:nInt){
        numSample <- nrow(samples$sample)
        subSample <- samples;
        sampleIndex <- sort(sample.int(numSample, round(0.5*numSample)))
        subSample$sample <- subSample$sample[sampleIndex,]
        lassoSub <- fraction_data_process(subSample$sample, 1 - NoAssociation, interactionMatrix)
        selectLambda <<- selection(lassoSub, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
        lambdaOptimal <<- (which.min(selectLambda)-1)*0.01 + 0.01;
        resultPLSample <- algorithm_select(lassoSub, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
        reproduceErrorPL <- repro_eval(resultPL$adj, resultPLSample$adj)
        reproduceErrorArrayPL[k,] <- reproduceErrorPL;
       
      }
      meanPL <- colMeans(1-reproduceErrorArrayPL)
      stdPL <- colSds(1-reproduceErrorArrayPL)
      print("MPLasso HMASM Evaluation Results: 25% (std), 50% (std), 75% (std), 100% (std)")
      print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", meanPL[1], stdPL[1], meanPL[2], stdPL[2], meanPL[3], stdPL[3], meanPL[4], stdPL[4]))
    }
}

  #}  

