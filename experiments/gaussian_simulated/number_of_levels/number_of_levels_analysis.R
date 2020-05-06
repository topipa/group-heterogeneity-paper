
# fix N to 300
# increase number of levels
# black with ARD-SE kernel
# red with linear + ARD-SE kernel

batches <- 50
par(mfrow=c(5,2))
par(mar = c(1,1,1,1))

Ngs <- c(2,3,4,5,6,10,12,15,20,30)
Nps <- c(150,100,75,60,50,30,25,20,15,10)

for (plo in seq(length(Ngs))) {
  Ng <- Ngs[plo]
  Np <- Nps[plo]
  results_list <- readRDS(paste('tests/gaussian_simulated/triton/newest/data/D_10_N_300_Ng_',Ng,'_Np_',Np,'_rep_1_batch_1_seed_1_.rds', sep = ''))
  for (batch in seq(2,batches)) {
    results_list1 <- readRDS(paste('tests/gaussian_simulated/triton/newest/data/D_10_N_300_Ng_',Ng,'_Np_',Np,'_rep_1_batch_',batch,'_seed_1_.rds', sep = ''))
    results_list$relevance_ranks2 <- cbind(results_list$relevance_ranks2,results_list1$relevance_ranks2[,1])
    results_list$relevance_ranks4 <- cbind(results_list$relevance_ranks4,results_list1$relevance_ranks4[,1])
    
  }
  
  repeats <- dim(results_list$relevance_ranks2)[2]
  
  plot(results_list$group_sd,rowMeans(results_list$relevance_ranks2), xlim = c(-0.05,1.1), ylim = c(0.9,9))
  points(results_list$group_sd+0.01,rowMeans(results_list$relevance_ranks4), ylim = c(0.9,10), col = 'red')
  arrows(results_list$group_sd, rowMeans(results_list$relevance_ranks2) - sqrt(matrixStats::rowVars(results_list$relevance_ranks2)/repeats), results_list$group_sd, rowMeans(results_list$relevance_ranks2) + sqrt(matrixStats::rowVars(results_list$relevance_ranks2)/repeats), length=0.05, angle=90, code=3)
  arrows(results_list$group_sd+0.01, rowMeans(results_list$relevance_ranks4) - sqrt(matrixStats::rowVars(results_list$relevance_ranks4)/repeats), results_list$group_sd+0.01, rowMeans(results_list$relevance_ranks4) + sqrt(matrixStats::rowVars(results_list$relevance_ranks4)/repeats), length=0.05, angle=90, code=3, col = 'red')
  if (plo == 2) {
    legend( "topright", c("SE", "LIN+SE"), 
            col=c("black", "red"),
            lty = c(1,1))
  }
}






# fixed N per level to 10
# increase number of levels


batches <- 50
par(mfrow=c(5,2))
par(mar = c(1,1,1,1))


Ngs <- c(2,4,6,8,10,15,20,30,40,50)
Ns <- c(20,40,60,80,100,150,200,300,400,500)

for (plo in seq(length(Ngs))) {
  Ng <- Ngs[plo]
  N <- Ns[plo]
  results_list <- readRDS(paste('tests/gaussian_simulated/triton/newest/data/D_10_N_',N,'_Ng_',Ng,'_Np_10_rep_1_batch_1_seed_1_.rds', sep = ''))
  for (batch in seq(2,batches)) {
    results_list1 <- readRDS(paste('tests/gaussian_simulated/triton/newest/data/D_10_N_',N,'_Ng_',Ng,'_Np_10_rep_1_batch_',batch,'_seed_1_.rds', sep = ''))
    results_list$relevance_ranks2 <- cbind(results_list$relevance_ranks2,results_list1$relevance_ranks2[,1])
    results_list$relevance_ranks4 <- cbind(results_list$relevance_ranks4,results_list1$relevance_ranks4[,1])
    
  }
  
  repeats <- dim(results_list$relevance_ranks2)[2]
  
  plot(results_list$group_sd,rowMeans(results_list$relevance_ranks2), xlim = c(-0.05,1.1), ylim = c(0.9,9))
  points(results_list$group_sd+0.01,rowMeans(results_list$relevance_ranks4), ylim = c(0.9,10), col = 'red')
  arrows(results_list$group_sd, rowMeans(results_list$relevance_ranks2) - sqrt(matrixStats::rowVars(results_list$relevance_ranks2)/repeats), results_list$group_sd, rowMeans(results_list$relevance_ranks2) + sqrt(matrixStats::rowVars(results_list$relevance_ranks2)/repeats), length=0.05, angle=90, code=3)
  arrows(results_list$group_sd+0.01, rowMeans(results_list$relevance_ranks4) - sqrt(matrixStats::rowVars(results_list$relevance_ranks4)/repeats), results_list$group_sd+0.01, rowMeans(results_list$relevance_ranks4) + sqrt(matrixStats::rowVars(results_list$relevance_ranks4)/repeats), length=0.05, angle=90, code=3, col = 'red')
  
}


























