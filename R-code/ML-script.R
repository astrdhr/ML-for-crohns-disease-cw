# Name: AS
# Date: Thursday 10th December 2020
# Can ML reliably diagnose Crohn's disease from patient GCMS data?

# IMPORT LIBRARIES ----

library(classyfire)
library(R.matlab)
library(ggplot2)
library(scatterplot3d)


# IMPORT DATASETS ----

bl <- readMat("blood/BWG_BL_CDvCTRL.mat") # blood samples
br <- readMat("breath/BWG_BR_CDvCTRL.mat") # breath samples
fa <- readMat("faecal/BWG_FA_CDvCTRL.mat") # faecal samples
ur <- readMat("urine/BWG_UR_CDvCTRL.mat") # urine samples


# GET INPUT DATA MATRIX ----

# blood
bl.X <- bl$XTIC
bl.Y <- bl$CLASS

# breath
br.X <- br$XTIC
br.Y <- br$CLASS

# faecal
fa.X <- fa$XTIC
fa.Y <- fa$CLASS

# urine
ur.X <- ur$XTIC
ur.Y <- ur$CLASS


# COPY SAMPLE NAMES FROM SAM TO ROW NAMES OF X ----

rownames(bl.X) <- as.character(unlist(bl$SAM))
rownames(br.X) <- as.character(unlist(br$SAM))
rownames(fa.X) <- as.character(unlist(fa$SAM))
rownames(ur.X) <- as.character(unlist(ur$SAM))


# EXPLORATORY ANALYSIS ----

# X data matrix contains GCMS data for each sample.
# Dataset Y contains disease state, 1 = control and 2 = disease.

# blood
summary(bl.X)
dim(bl.X)
table(bl.X)

summary(bl.Y)
dim(bl.Y)
table(bl.Y)

# breath
summary(br.X)
dim(br.X)
table(br.X)

summary(br.Y)
dim(br.Y)
table(br.Y)

# faecal
summary(fa.X)
dim(fa.X)
table(fa.X)

summary(fa.Y)
dim(fa.Y)
table(fa.Y)

# urine
summary(ur.X)
dim(ur.X)
table(ur.X)

summary(ur.Y)
dim(ur.Y)
table(ur.Y)


# EXPLORE RETENTION TIMES (RT)

# blood
bl.RT <- bl$RT
i = 23 # set index number to sample you want to plot
plot(bl.RT, bl.X[i,], type="l", xlab="Blood RT (min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(bl.X)[i]))

# breath
br.RT <- br$RT
i = 23 # set index number to sample you want to plot
plot(br.RT, br.X[i,], type="l", xlab="Breath RT (min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(br.X)[i]))

# faecal
fa.RT <- fa$RT
i = 23 # set index number to sample you want to plot
plot(fa.RT, fa.X[i,], type="l", xlab="Faecal RT (min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(fa.X)[i]))

# urine
ur.RT <- ur$RT
i = 22 # set index number to sample you want to plot
plot(ur.RT, ur.X[i,], type="l", xlab="Urine RT (min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(ur.X)[i]))


# DATA PRE-PROCESSING ----

# Pre-process data using PCA
bl.Xpca <- prcomp(bl.X, center = TRUE, scale. = TRUE)
br.Xpca <- prcomp(br.X, center = TRUE, scale. = TRUE)
fa.Xpca <- prcomp(fa.X, center = TRUE, scale. = TRUE)
ur.Xpca <- prcomp(ur.X, center = TRUE, scale. = TRUE)

# Get variance summary
bl.summ <- summary(bl.Xpca)
br.summ <- summary(br.Xpca)
fa.summ <- summary(fa.Xpca)
ur.summ <- summary(ur.Xpca)

# Extract proportion of variance for each PC
bl.exp_var = bl.summ$importance[2,] * 100
br.exp_var = br.summ$importance[2,] * 100
fa.exp_var = fa.summ$importance[2,] * 100
ur.exp_var = ur.summ$importance[2,] * 100

# Extract cumulative proportion of variance for each PC
bl.cum_prop <- bl.summ$importance[3,] * 100
br.cum_prop <- br.summ$importance[3,] * 100
fa.cum_prop <- fa.summ$importance[3,] * 100
ur.cum_prop <- ur.summ$importance[3,] * 100

# Plot proportion of variance for each PC
plot(bl.exp_var, type = "l")
plot(br.exp_var, type = "l")
plot(fa.exp_var, type = "l")
plot(ur.exp_var, type = "l")

# Plot cumulative proportion of variance for each PC
plot(bl.cum_prop, type = "l")
plot(br.cum_prop, type = "l")
plot(fa.cum_prop, type = "l")
plot(ur.cum_prop, type = "l")

# Retain PCs that capture greatest variance
bl.Xscores <- bl.Xpca$x[,1:5]
br.Xscores <- br.Xpca$x[,1:3]
fa.Xscores <- fa.Xpca$x[,1:8]
ur.Xscores <- ur.Xpca$x[,1:4]

# Check PC scores
bl.Xscores
br.Xscores
fa.Xscores
ur.Xscores

# Convert PC scores to dataframe for plotting
bl.Xscores.df <- as.data.frame(bl.Xscores)
bl.Xscores.df <- as.data.frame(bl.Xscores)
bl.Xscores.df <- as.data.frame(bl.Xscores)
bl.Xscores.df <- as.data.frame(bl.Xscores)

# PLOT PCA

# Colour PCA by treatment groups
blue <- "#0096FF"
salmon <- "#F8766D"
my_cols <- c(blue, salmon)
names(my_cols) <- c("1", "2") # 1 for control, 2 for CD

# Set plotting frame
par(mfrow = c(2,2))

# Plot PCAs

# blood
ggplot(as.data.frame(bl.Xscores), aes(x = PC1, y = PC2)) +
  geom_point()

plot(bl.Xscores[,1], bl.Xscores[,5], xlab = "PC1 (60.4%)", ylab = "PC2 (15.0%)", pch = 19, 
     col = my_cols[bl.Y], main = "PCA scores (per blood sample)")

# breath
ggplot(as.data.frame(br.Xscores), aes(x = PC1, y = PC2)) +
  geom_point()

plot(br.Xscores[,1], br.Xscores[,3], xlab = "PC1 (82.3%)", ylab = "PC2 (11.4%)", pch = 19, 
     col = my_cols[br.Y], main = "PCA scores (per breath sample)")

# faecal
ggplot(as.data.frame(fa.Xscores), aes(x = PC1, y = PC2)) +
  geom_point()

plot(fa.Xscores[,1], fa.Xscores[,2], xlab = "PC1 (35.1%)", ylab = "PC2 (21.7%)", pch = 19, 
     col = my_cols[fa.Y], main = "PCA scores (per faecal sample)")

# urine
ggplot(as.data.frame(ur.Xscores), aes(x = PC1, y = PC2)) +
  geom_point()

plot(ur.Xscores[,1], ur.Xscores[,2], xlab = "PC1 (75.4%)", ylab = "PC2 (10.5%)", pch = 19, 
     col = my_cols[ur.Y], main = "PCA scores (per urine sample)")

# Add legend
legend("topright", #position
       legend = c("Control", "Crohn's Disease"),
       fill = c(blue, salmon),
       border = NA, bty = "n")


# PLOT 3D PCAs

par(mfrow = c(1,1))

# blood
scatterplot3d(bl.Xscores[,2], bl.Xscores[,3], bl.Xscores[,1],
              color = my_cols[bl.Y], main = "PCA scores (per blood sample)")

# breath
scatterplot3d(br.Xscores[,1], br.Xscores[,2], br.Xscores[,3],
              color = my_cols[br.Y], main = "PCA scores (per breath sample)")

# faecal
scatterplot3d(fa.Xscores[,1], fa.Xscores[,2], fa.Xscores[,3],
              color = my_cols[fa.Y], main = "PCA scores (per faecal sample)")

# urine
scatterplot3d(ur.Xscores[,2], ur.Xscores[,3], ur.Xscores[,1],
              color = my_cols[ur.Y], main = "PCA scores (per urine sample)")


# BUILD AND TEST CLASSIFICATION ENSEMBLE / ML MODELS ----

# blood
bl.mlm1 <- cfBuild(inputData = bl.Xscores,
                 inputClass = bl.Y,
                 bootNum = 50,
                 ensNum = 50,
                 parallel = TRUE,
                 cpus = 4)


bl.mlm2 <- cfBuild(inputData = bl.Xscores,
                   inputClass = bl.Y,
                   bootNum = 100,
                   ensNum = 100,
                   parallel = TRUE,
                   cpus = 4)

bl.mlm3 <- cfBuild(inputData = bl.Xscores,
                   inputClass = bl.Y,
                   bootNum = 150,
                   ensNum = 150,
                   parallel = TRUE,
                   cpus = 4)

bl.mlm4 <- cfBuild(inputData = bl.Xscores,
                   inputClass = bl.Y,
                   bootNum = 200,
                   ensNum = 200,
                   parallel = TRUE,
                   cpus = 4)

bl.mlm5 <- cfBuild(inputData = bl.Xscores,
                   inputClass = bl.Y,
                   bootNum = 125,
                   ensNum = 125,
                   parallel = TRUE,
                   cpus = 4)

bl.mlm6 <- cfBuild(inputData = bl.Xscores,
                   inputClass = bl.Y,
                   bootNum = 75,
                   ensNum = 75,
                   parallel = TRUE,
                   cpus = 4)



# breath
br.mlm1 <- cfBuild(inputData = br.Xscores,
                 inputClass = br.Y,
                 bootNum = 50,
                 ensNum = 50,
                 parallel = TRUE,
                 cpus = 4)

br.mlm2 <- cfBuild(inputData = br.Xscores,
                   inputClass = br.Y,
                   bootNum = 100,
                   ensNum = 100,
                   parallel = TRUE,
                   cpus = 4)

br.mlm3 <- cfBuild(inputData = br.Xscores,
                   inputClass = br.Y,
                   bootNum = 150,
                   ensNum = 150,
                   parallel = TRUE,
                   cpus = 4)

br.mlm4 <- cfBuild(inputData = br.Xscores,
                   inputClass = br.Y,
                   bootNum = 200,
                   ensNum = 200,
                   parallel = TRUE,
                   cpus = 4)

br.mlm5 <- cfBuild(inputData = br.Xscores,
                   inputClass = br.Y,
                   bootNum = 75,
                   ensNum = 75,
                   parallel = TRUE,
                   cpus = 4)

br.mlm6 <- cfBuild(inputData = br.Xscores,
                   inputClass = br.Y,
                   bootNum = 125,
                   ensNum = 125,
                   parallel = TRUE,
                   cpus = 4)

# faecal
fa.mlm1 <- cfBuild(inputData = fa.Xscores,
                   inputClass = fa.Y,
                   bootNum = 50,
                   ensNum = 50,
                   parallel = TRUE,
                   cpus = 4)

fa.mlm2 <- cfBuild(inputData = fa.Xscores,
                   inputClass = fa.Y,
                   bootNum = 100,
                   ensNum = 100,
                   parallel = TRUE,
                   cpus = 4)

fa.mlm3 <- cfBuild(inputData = fa.Xscores,
                   inputClass = fa.Y,
                   bootNum = 150,
                   ensNum = 150,
                   parallel = TRUE,
                   cpus = 4)

fa.mlm4 <- cfBuild(inputData = fa.Xscores,
                   inputClass = fa.Y,
                   bootNum = 200,
                   ensNum = 200,
                   parallel = TRUE,
                   cpus = 4)

fa.mlm5 <- cfBuild(inputData = fa.Xscores,
                   inputClass = fa.Y,
                   bootNum = 75,
                   ensNum = 75,
                   parallel = TRUE,
                   cpus = 4)

fa.mlm6 <- cfBuild(inputData = fa.Xscores,
                   inputClass = fa.Y,
                   bootNum = 125,
                   ensNum = 125,
                   parallel = TRUE,
                   cpus = 4)

# urine
ur.mlm1 <- cfBuild(inputData = ur.Xscores,
                   inputClass = ur.Y,
                   bootNum = 50,
                   ensNum = 50,
                   parallel = TRUE,
                   cpus = 4)

ur.mlm2 <- cfBuild(inputData = ur.Xscores,
                   inputClass = ur.Y,
                   bootNum = 100,
                   ensNum = 100,
                   parallel = TRUE,
                   cpus = 4)

ur.mlm3 <- cfBuild(inputData = ur.Xscores,
                   inputClass = ur.Y,
                   bootNum = 150,
                   ensNum = 150,
                   parallel = TRUE,
                   cpus = 4)

ur.mlm4 <- cfBuild(inputData = ur.Xscores,
                   inputClass = ur.Y,
                   bootNum = 200,
                   ensNum = 200,
                   parallel = TRUE,
                   cpus = 4)

ur.mlm5 <- cfBuild(inputData = ur.Xscores,
                   inputClass = ur.Y,
                   bootNum = 75,
                   ensNum = 75,
                   parallel = TRUE,
                   cpus = 4)

ur.mlm6 <- cfBuild(inputData = ur.Xscores,
                   inputClass = ur.Y,
                   bootNum = 125,
                   ensNum = 125,
                   parallel = TRUE,
                   cpus = 4)


# PERMUTATION TESTING ----
# To determine statistical significance of performance

# blood
bl.perm <- cfPermute(inputData = bl.Xscores,
                     inputClass = bl.Y,
                     bootNum = 50,
                     ensNum = 50,
                     permNum = 50,
                     parallel = TRUE,
                     cpus = 4)

# breath
br.perm <- cfPermute(inputData = br.Xscores,
                     inputClass = br.Y,
                     bootNum = 50,
                     ensNum = 50,
                     permNum = 50,
                     parallel = TRUE,
                     cpus = 4)

# faecal
fa.perm <- cfPermute(inputData = fa.Xscores,
                     inputClass = fa.Y,
                     bootNum = 50,
                     ensNum = 50,
                     permNum = 50,
                     parallel = TRUE,
                     cpus = 4)

# urine
ur.perm <- cfPermute(inputData = ur.Xscores,
                     inputClass = ur.Y,
                     bootNum = 50,
                     ensNum = 50,
                     permNum = 50,
                     parallel = TRUE,
                     cpus = 4)


# GET STATISTICS FOR CLASSIFICATION ENSEMBLE ----

# GET AVERAGE TEST AND TRAIN ACCURACY OF ENSEMBLE

# blood
getAvgAcc(bl.mlm1)$Test
getAvgAcc(bl.mlm1)$Train

getAvgAcc(bl.mlm2)$Test
getAvgAcc(bl.mlm2)$Train

getAvgAcc(bl.mlm3)$Test
getAvgAcc(bl.mlm3)$Train

getAvgAcc(bl.mlm4)$Test
getAvgAcc(bl.mlm4)$Train

getAvgAcc(bl.mlm5)$Test
getAvgAcc(bl.mlm5)$Train

# breath
getAvgAcc(br.mlm1)$Test
getAvgAcc(br.mlm1)$Train

getAvgAcc(br.mlm2)$Test
getAvgAcc(br.mlm2)$Train

getAvgAcc(br.mlm3)$Test
getAvgAcc(br.mlm3)$Train

getAvgAcc(br.mlm4)$Test
getAvgAcc(br.mlm4)$Train

getAvgAcc(br.mlm5)$Test
getAvgAcc(br.mlm5)$Train

getAvgAcc(br.mlm6)$Test
getAvgAcc(br.mlm6)$Train

# faecal
getAvgAcc(fa.mlm1)$Test
getAvgAcc(fa.mlm1)$Train

getAvgAcc(fa.mlm2)$Test
getAvgAcc(fa.mlm2)$Train

getAvgAcc(fa.mlm3)$Test
getAvgAcc(fa.mlm3)$Train

getAvgAcc(fa.mlm4)$Test
getAvgAcc(fa.mlm4)$Train

getAvgAcc(fa.mlm5)$Test
getAvgAcc(fa.mlm5)$Train

# urine
getAvgAcc(ur.mlm1)$Test
getAvgAcc(ur.mlm1)$Train

getAvgAcc(ur.mlm2)$Test
getAvgAcc(ur.mlm2)$Train

getAvgAcc(ur.mlm3)$Test
getAvgAcc(ur.mlm3)$Train

getAvgAcc(ur.mlm4)$Test
getAvgAcc(ur.mlm4)$Train

getAvgAcc(ur.mlm5)$Test
getAvgAcc(ur.mlm5)$Train

getAvgAcc(ur.mlm6)$Test
getAvgAcc(ur.mlm6)$Train


# GET TEST AND TRAIN ACCURACY OF CLASSIFIERS OF SAMPLE

# blood
getAcc(bl.mlm1)$Test
getAcc(bl.mlm1)$Train

getAcc(bl.mlm2)$Test
getAcc(bl.mlm2)$Train

getAcc(bl.mlm3)$Test
getAcc(bl.mlm3)$Train

getAcc(bl.mlm4)$Test
getAcc(bl.mlm4)$Train

getAcc(bl.mlm5)$Test
getAcc(bl.mlm5)$Train

# breath
getAcc(br.mlm1)$Test
getAcc(br.mlm1)$Train

getAcc(br.mlm2)$Test
getAcc(br.mlm2)$Train

getAcc(br.mlm3)$Test
getAcc(br.mlm3)$Train

getAcc(br.mlm4)$Test
getAcc(br.mlm4)$Train

getAcc(br.mlm5)$Test
getAcc(br.mlm5)$Train

getAcc(br.mlm6)$Test
getAcc(br.mlm6)$Train

# faecal
getAcc(fa.mlm1)$Test
getAcc(fa.mlm1)$Train

getAcc(fa.mlm2)$Test
getAcc(fa.mlm2)$Train

getAcc(fa.mlm3)$Test
getAcc(fa.mlm3)$Train

getAcc(fa.mlm4)$Test
getAcc(fa.mlm4)$Train

getAcc(fa.mlm5)$Test
getAcc(fa.mlm5)$Train

# urine
getAcc(ur.mlm1)$Test
getAcc(ur.mlm1)$Train

getAcc(ur.mlm2)$Test
getAcc(ur.mlm2)$Train

getAcc(ur.mlm3)$Test
getAcc(ur.mlm3)$Train

getAcc(ur.mlm4)$Test
getAcc(ur.mlm4)$Train

getAcc(ur.mlm5)$Test
getAcc(ur.mlm5)$Train

getAcc(ur.mlm6)$Test
getAcc(ur.mlm6)$Train


# GET OVERALL CONFUSION MATRIX OF ENSEMBLE
# Summarises performance of ensemble
getConfMatr(bl.mlm4)
getConfMatr(br.mlm4)
getConfMatr(fa.mlm4)
getConfMatr(ur.mlm4)


# GET OPTIMAL SVM HYPERPARAMETERS OF CLASSIFICATION ENSEMBLE
bl.optParam <- getOptParam(bl.mlm4)
plot(bl.optParam)

br.optParam <- getOptParam(br.mlm4)
plot(br.optParam)

fa.optParam <- getOptParam(fa.mlm4)
plot(fa.optParam)

ur.optParam <- getOptParam(ur.mlm4)
plot(ur.optParam)


# GET PERMUTATION STATISTICS ----

# GET DESCRIPTIVE STATISTICS OF PERMUTATION RUNS

# blood
getPerm5Num(bl.perm)
getPerm5Num(bl.perm)$median      
getPerm5Num(bl.perm)$minimum
getPerm5Num(bl.perm)$maximum
getPerm5Num(bl.perm)$upperQ
getPerm5Num(bl.perm)$lowerQ

# breath
getPerm5Num(br.perm)
getPerm5Num(br.perm)$median      
getPerm5Num(br.perm)$minimum
getPerm5Num(br.perm)$maximum
getPerm5Num(br.perm)$upperQ
getPerm5Num(br.perm)$lowerQ

# faecal
getPerm5Num(fa.perm)
getPerm5Num(fa.perm)$median      
getPerm5Num(fa.perm)$minimum
getPerm5Num(fa.perm)$maximum
getPerm5Num(fa.perm)$upperQ
getPerm5Num(fa.perm)$lowerQ

# urine
getPerm5Num(ur.perm)
getPerm5Num(ur.perm)$median      
getPerm5Num(ur.perm)$minimum
getPerm5Num(ur.perm)$maximum
getPerm5Num(ur.perm)$upperQ
getPerm5Num(ur.perm)$lowerQ


# GET AVERAGED ACCURACIES FOR EACH PERMUTATION
bl.perm$avgAcc
br.perm$avgAcc
fa.perm$avgAcc
ur.perm$avgAcc


# GET TOTAL EXECUTION TIMES FOR PERMUTATION
bl.perm$totalTime
br.perm$totalTime
fa.perm$totalTime
ur.perm$totalTime


# VISUALISE RESULTS ----

# VISUALISE AVERAGE TEST ACCURACY OVER ENSEMBLE ITERATIONS

# blood
ggEnsTrend(bl.mlm1, showText  = TRUE) # checked
ggEnsTrend(bl.mlm2)
ggEnsTrend(bl.mlm3) # optimal, re-run
ggEnsTrend(bl.mlm4)
ggEnsTrend(bl.mlm5)
ggEnsTrend(bl.mlm6)

# breath
ggEnsTrend(br.mlm1, showText  = TRUE)
ggEnsTrend(br.mlm2)
ggEnsTrend(br.mlm3)
ggEnsTrend(br.mlm4)
ggEnsTrend(br.mlm5)
ggEnsTrend(br.mlm6)

# faecal
ggEnsTrend(fa.mlm1, showText  = TRUE)
ggEnsTrend(fa.mlm2) # optimal, re-run
ggEnsTrend(fa.mlm3)
ggEnsTrend(fa.mlm4)
ggEnsTrend(fa.mlm5)
ggEnsTrend(fa.mlm6)

# urine
ggEnsTrend(ur.mlm1, showText  = TRUE)
ggEnsTrend(ur.mlm2)
ggEnsTrend(ur.mlm3)
ggEnsTrend(ur.mlm4)
ggEnsTrend(ur.mlm5)
ggEnsTrend(ur.mlm6)


## VISUALISE CLASSIFICATION ACCURACY

# Set plotting frame
par(mfrow = c(2,2))

# VISUALISE TREND OF TEST ACCURACIES OVER ENSEMBLE ITERATIONS - this is repeat, see above. Retain this one, delete the other.
ggEnsTrend(bl.mlm4)
ggEnsTrend(br.mlm4)
ggEnsTrend(fa.mlm4)
ggEnsTrend(ur.mlm4)

# BARPLOT OF % CLASSIFIED AND MISCLASSIFED SAMPLES
ggClassPred(bl.mlm4, position = "stack", displayAll = TRUE, showText = TRUE)
ggClassPred(br.mlm4, position = "stack", displayAll = TRUE, showText = TRUE)
ggClassPred(fa.mlm4, position = "stack", displayAll = TRUE, showText = TRUE)
ggClassPred(ur.mlm4, position = "stack", displayAll = TRUE, showText = TRUE)

# DENSITY PLOT OF TEST ACCURACIES IN ENSEMBLE
ggEnsHist(bl.mlm4, density = TRUE, percentiles=TRUE, mean=TRUE, median=TRUE)
ggEnsHist(br.mlm4, density = TRUE, percentiles=TRUE, mean=TRUE, median=TRUE)
ggEnsHist(fa.mlm4, density = TRUE, percentiles=TRUE, mean=TRUE, median=TRUE)
ggEnsHist(ur.mlm4, density = TRUE, percentiles=TRUE, mean=TRUE, median=TRUE)


## VISUALISE PERMUTATION RUN

# DENSITY PLOT OF PERMUTATION DISTRIBUTION
ggPermHist(bl.perm, density=TRUE, percentiles = TRUE, median = TRUE)
ggPermHist(br.perm, density=TRUE, percentiles = TRUE, median = TRUE)
ggPermHist(fa.perm, density=TRUE, percentiles = TRUE, median = TRUE)
ggPermHist(ur.perm, density=TRUE, percentiles = TRUE, median = TRUE)

# COMPARE ENSEMBLE AND PERMTUATION DISTRIBUTION
ggFusedHist(bl.mlm4, bl.perm)
ggFusedHist(br.mlm4, br.perm)
ggFusedHist(fa.mlm4, fa.perm)
ggFusedHist(ur.mlm4, ur.perm)


## VISUALISE FUSED PERMUTATION HISTOGRAM ----

# blood

bl.fusedMatr <- rbind( cbind(rep("permutation", length(bl.perm$avgAcc)), bl.perm$avgAcc), 
                       cbind(rep("ensemble", length(bl.mlm4$testAcc)), bl.mlm4$testAcc))

bl.fusedMatr <- as.data.frame(bl.fusedMatr)

colnames(bl.fusedMatr) <- c("type", "acc")

ggplot(data = as.data.frame(bl.fusedMatr), aes(x = as.numeric(as.vector(acc)), fill = type, colour = type)) +
  geom_histogram(binwidth=2, alpha=.6, colour = "#999999") + 
  theme_bw() + 
  scale_fill_manual(values = c("red", "turquoise3")) + 
  xlab("Overall Test Accuracies (%CC)") + 
  ylab("Frequency\n") + 
  geom_vline(xintercept = median(bl.perm$avgAcc), color = "dodgerblue1", linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = median(bl.mlm4$testAcc), color = "pink", linetype = "dashed", size = 0.8)

# breath

br.fusedMatr <- rbind( cbind(rep("permutation", length(br.perm$avgAcc)), br.perm$avgAcc), 
                       cbind(rep("ensemble", length(br.mlm4$testAcc)), br.mlm4$testAcc))

br.fusedMatr <- as.data.frame(br.fusedMatr)

colnames(br.fusedMatr) <- c("type", "acc")

ggplot(data = as.data.frame(br.fusedMatr), aes(x = as.numeric(as.vector(acc)), fill = type, colour = type)) + 
  geom_histogram(binwidth=2, alpha=.6, colour = "#999999") + 
  theme_bw() + 
  scale_fill_manual(values = c("red", "turquoise3")) + 
  xlab("Overall Test Accuracies (%CC)") + 
  ylab("Frequency\n") + 
  geom_vline(xintercept = median(br.perm$avgAcc), color = "dodgerblue1", linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = median(br.mlm4$testAcc), color = "pink", linetype = "dashed", size = 0.8)


# faecal

fa.fusedMatr <- rbind( cbind(rep("permutation", length(fa.perm$avgAcc)), fa.perm$avgAcc), 
                       cbind(rep("ensemble", length(fa.mlm4$testAcc)), fa.mlm4$testAcc))

fa.fusedMatr <- as.data.frame(fa.fusedMatr)

colnames(fa.fusedMatr) <- c("type", "acc")

ggplot(data = as.data.frame(fa.fusedMatr), aes(x = as.numeric(as.vector(acc)), fill = type, colour = type)) + 
  geom_histogram(binwidth=2, alpha=.6, colour = "#999999") + 
  theme_bw() + 
  scale_fill_manual(values = c("red", "turquoise3")) + 
  xlab("Overall Test Accuracies (%CC)") + 
  ylab("Frequency\n") + 
  geom_vline(xintercept = median(fa.perm$avgAcc), color = "dodgerblue1", linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = median(fa.mlm4$testAcc), color = "pink", linetype = "dashed", size = 0.8)


# urine

ur.fusedMatr <- rbind( cbind(rep("permutation", length(ur.perm$avgAcc)), ur.perm$avgAcc), 
                       cbind(rep("ensemble", length(ur.mlm4$testAcc)), ur.mlm4$testAcc))

ur.fusedMatr <- as.data.frame(ur.fusedMatr)

colnames(ur.fusedMatr) <- c("type", "acc")

ggplot(data = as.data.frame(ur.fusedMatr), aes(x = as.numeric(as.vector(acc)), fill = type, colour = type)) + 
  geom_histogram(binwidth=2, alpha=.6, colour = "#999999") +
  theme_bw() + 
  scale_fill_manual(values = c("red", "turquoise3")) + 
  xlab("Overall Test Accuracies (%CC)") + 
  ylab("Frequency\n") + 
  geom_vline(xintercept = median(ur.perm$avgAcc), color = "dodgerblue1", linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = median(ur.mlm4$testAcc), color = "pink", linetype = "dashed", size = 0.8)