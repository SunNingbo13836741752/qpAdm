##Code to summarize qpAdm results from 2 or more way 
#Enable %>%
library(dplyr)
library(readr)
source("read_qpAdm.R")
source("read_qpWave.R")

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
print(args[2])
#List all files and gives the number of ways

logfiles <- list.files(args[1],pattern = "*\\.out") ##Becareful about outlier in the namespace
n <- as.numeric(args[2])


#logfiles <- list.files("/home/xiaowei_mao/NE/git/qpAdm/output/distal2/twoway",pattern = "*.out")
#n<- 2

#Summarize all possible outcomes in the final file
final <- NULL

#run the two functions at the bottom before excuting these block
for (i in 1:length(logfiles)){
  print(i)
  file <- logfiles[i]
  log_lines <- readLines(paste0(args[1],file)) %>% .[!stringr::str_detect(., "nodata")]
  result <- read_qpAdm(log_lines)
  target <- result$proportions$target
  source <- colnames(result$proportions)[c(2:(n+1))] ##the library of vctrs should to be updated, because dplyr::vars does not work, but I could not update it in r940, so I just changed the read_qpAdm.R to keep it as characters 
#  prop <- result$proportions[1,c(2:(n+1))]
  prop <- as.numeric(result$proportions[1,c(2:(n+1))])
#  se <- result$proportions[1,c((n+2):(2*n+1))]
  se <- as.numeric(result$proportions[1,c((n+2):(2*n+1))])
#  snps <- result$proportions["nsnps"]
  snps <- as.numeric(result$proportions["nsnps"])
  #propMinusE<- result$proportions[1,c(2:(2+n-1))]-2*se
  #propPlusE <- result$proportions[1,c(2:(2+n-1))]+2*se
  check1 <- sum(prop>=0)==n #coef estimates larger or equal than 0
  check2 <- sum(prop<=1)==n
  #check3 <- sum(propMinusE>=0)==n #coef estimates +/- 2*se larger or euqal than 0
  #check4 <- sum(propPlusE<=1)==n
  check5 <- result$subsets$comment[1]!="infeasible" #make sure the full model is ok
  check6 <- result$subsets$tail[1] >= 0.05  #Make sure the tail pass 0.05
  #if (check1&&check2&&check3&&check4&&check5&&check6){
  if (check1&&check2&&check5&&check6){
    nest <- result$nest
    nest <- nest[nest$comment!="infeasible",]
    nest <- nest[!is.na(nest$nestp),]
    check7 <- nest$nestp>0.05
    if (sum(check7)<1){
      info <- c(target,source,as.matrix(prop),as.matrix(se),result$subsets$tail[1],as.matrix(snps))
      final <- rbind(final,info)
    }
  }
}


if (is.null(final)){
  print("no results")
}else{
  colnames(final) <- c("Target",paste0("source",seq(1:n)),paste0("Coef",seq(1:n)),paste0("SE",seq(1:n)),"Tail","nSNPs")
  final <- as.data.frame(final,row.names = F)
  write.table(final,file=paste0(n,"wayfinal.txt"),quote=F,row.names=F,col.names=T,sep = " ")
}
