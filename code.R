rm(list = ls())

setwd('C:\\Users\\uccoo\\Desktop\\학교\\대학원1\\통분방\\과제\\final project')
library(foreign)
library(poLCA)

#### 1. 데이터 전처리 + 리코딩 ####
data_sat <- read.dta("04264-0002-Data.dta")
var <- paste0('V', c(1646:1648, 1650, 1653,
                     1102, 1216, 1254, 1712,
                     1150, 1151, 1167, 1170))
data_sat_1 <- data_sat[,var]
colnames(data_sat_1) <- c(paste0('q', 1:5), paste0('d', 1:4),
                          'gender', 'race', 'political', 'religion')
summary(data_sat_1)
summary(data_sat_1$political)

# 결측치 제거 : covariates의 결측치는 제거해야 함
na_obs <- which(apply(data_sat_1[,10:13], 1, function(x) any(x == 'MISSING')))
length(na_obs) # 1120개의 결측치
data_sat_2 <- data_sat_1[-na_obs,]

str(data_sat_2) # 최종적으로 1443개의 관측치

# 5개의 설문 문항 전처리 : 1~7의 답변, 1 <- 불만족 / 만족 -> 7
# 1-2 : 불만족, 3-5 : 중간, 6-7 : 만족으로 리코딩
# missing value : -9 -> NA로 변경
for(i in 1:5){
  q <- data_sat_2[,i]
  not.sat <- which(q %in% 1:2); data_sat_2[not.sat, i] <- 1
  mid <- which(q %in% 3:5); data_sat_2[mid, i] <- 2
  sat <- which(q %in% 6:7); data_sat_2[sat, i] <- 3
  missing <- which(q == -9); data_sat_2[missing, i] <- NA
}


for(i in 6:9){
  data_sat_2[,i] <- as.numeric(data_sat_2[,i])
  d <- data_sat_2[,i]
  none <- which(d == 2); data_sat_2[none, i] <- 1
  often <- which(d %in% 3:5); data_sat_2[often, i] <- 2
  lot <- which(d %in% 6:8); data_sat_2[lot, i] <- 3
  missing <- which(d == 1); data_sat_2[missing, i] <- NA
}

summary(data_sat_2)

# 4개의 covariates 리코딩
# gender
data_sat_2$gender <- factor(as.numeric(data_sat_2$gender), 
                            levels = 2:3, labels = c('male', 'female'))
# race
data_sat_2$race <- factor(as.numeric(data_sat_2$race), 
                         levels = 2:3, labels = c('white', 'black'))

# political belief
cons <- which(as.numeric(data_sat_2$political) %in% 2:3)
mod <- which(as.numeric(data_sat_2$political) == 4)
lib <- which(as.numeric(data_sat_2$political) %in% 5:7)
none <- which(as.numeric(data_sat_2$political) == 8)

poli <- rep(0, nrow(data_sat_2))
poli[cons] <- 'CONS'
poli[mod] <- 'MOD'
poli[lib] <- 'LIB'
poli[none] <- 'NONE'
data_sat_2$political <- factor(poli, levels = c('NONE', 'CONS', 'MOD', 'LIB'))

# religion importance
not_imp <- which(as.numeric(data_sat_2$religion) %in% 2:3)
imp <- which(as.numeric(data_sat_2$religion) %in% 4:5)

rel <- rep(0, nrow(data_sat_2))
rel[not_imp] <- 'not_imp'
rel[imp] <- 'imp'
data_sat_2$religion <- factor(rel, levels = c('not_imp', 'imp'))

summary(data_sat_2)
data_sat <- data_sat_2

# response Table
lapply(lapply(apply(data_sat, 2, table, useNA = 'ifany'), '/', nrow(data_sat)), round, 3)


#### 2. 카이제곱 독립성 검정 ####
chi.sq.test <- list()
comb <- combn(4, 2)
for(i in 1:ncol(comb)){
  comb.var <- comb[,i]+9
  data <- data_sat[,comb.var]
  
  chi.sq.test[[i]] <- chisq.test(table(data), correct = FALSE)
  names(chi.sq.test)[i] <- paste(colnames(data)[1], '&', colnames(data)[2])
}

chi.sq.test

#### 3. 로그-선형 모형 ####

loglin.data <- matrix(0, 32, 5); loglin.data <- as.data.frame(loglin.data)
colnames(loglin.data) <- c(colnames(data_sat)[10:13], 'count')

iter <- 1
for(i in 1:2){for(j in 1:2){for(k in 1:4){for(l in 1:2){
  x1 <- as.numeric(data_sat$gender) == i; loglin.data[iter, 1] <- levels(data_sat$gender)[i]
  x2 <- as.numeric(data_sat$race) == j; loglin.data[iter, 2] <- levels(data_sat$race)[j]
  x3 <- as.numeric(data_sat$political) == k; loglin.data[iter, 3] <- levels(data_sat$political)[k]
  x4 <- as.numeric(data_sat$religion) == l; loglin.data[iter, 4] <- levels(data_sat$religion)[l]
  
  loglin.data[iter, 5] <- sum(x1 & x2 & x3 & x4)
  iter <- iter+1
}}}}

sum(loglin.data$count) # observation 수 확인

fit_full_0 <- glm(count ~ (gender + race + political + religion)^4, data = loglin.data,
                  family = poisson())

fit_full <- glm(count ~ gender + race + political + religion +
                  gender*political + gender*religion + race*political + 
                  race*religion + political*religion +
                  race*political*religion + gender*political*religion, data = loglin.data,
                family = poisson())

a0 <- anova(fit_full, fit_full_0)
1-pchisq(a0$Deviance[2], a0$Df[2]) # p-value = 0.87, reduced model 선택

summary(fit_full) # 3차 교호작용 항이 유의하지 않아보이므로 reduced model 적합 후 비교

# gender:political:religion만 제거한 reduced model 적합
fit_1 <- glm(count ~ gender + race + political + religion +
                  gender*political + gender*religion + race*political + 
                  race*religion + political*religion +
                  race*political*religion, data = loglin.data, family = poisson())

a1 <- anova(fit_1, fit_full)
1-pchisq(a1$Deviance[2], a1$Df[2]) # p-value < 0.05


# race:political:religion만 제거한 reduced model 적합
fit_2 <- glm(count ~ gender + race + political + religion +
               gender*political + gender*religion + race*political + 
               race*religion + political*religion +
               gender*political*religion, data = loglin.data, family = poisson())


a2 <- anova(fit_2, fit_full)
1-pchisq(a2$Deviance[2], a2$Df[2]) # p-value = 0.08, reduced model 선택

summary(fit_2)

# 3차 교호작용 항 모두 제거
fit_3 <- glm(count ~ gender + race + political + religion +
               gender*political + gender*religion + race*political + 
               race*religion + political*religion, data = loglin.data,
             family = poisson())

a3 <- anova(fit_3, fit_2)
1-pchisq(a3$Deviance[2], a3$Df[2]) # p-value < 0.05

# race:political 교호작용 항 제거
fit_4 <- glm(count ~ gender + race + political + religion +
               gender*political + gender*religion + race*religion + 
               political*religion + gender*political*religion, data = loglin.data,
             family = poisson())

a4 <- anova(fit_4, fit_2)
1-pchisq(a4$Deviance[2], a4$Df[2]) # p-value = 0.12, reduced model 선택

summary(fit_4)

# 3차 교호작용 다시 제거
fit_5 <- glm(count ~ gender + race + political + religion +
               gender*political + gender*religion + race*religion + 
               political*religion, data = loglin.data,
             family = poisson())

a5 <- anova(fit_5, fit_4)
1-pchisq(a5$Deviance[2], a5$Df[2]) # p-value < 0.05

# fit_4를 최종 모델로 선택
summary(fit_4)

#### 4. LCA with covariates ####
# satisfaction model fitting
f <- cbind(q1, q2, q3, q4, q5)~1
idf_check <- list('nclass = 2' = NULL, 
                  'nclass = 3' = NULL, 'nclass = 4' = NULL, 
                  'nclass = 5' = NULL, 'nclass = 6' = NULL)
model_selection <- matrix(0, 5, 8)
colnames(model_selection) <- c('nclass','npar','G^2', 'df', 'p-value',
                               'AIC','BIC','likelihood')

for(i in 1:5){
  set.seed(2021021352)
  model_selection[i,1] <- i+1
  model <- poLCA(f, data = data_sat[,1:5], nclass = i+1, na.rm = F, verbose = F, nrep = 50)
  model_selection[i,2] <- model$npar
  model_selection[i,3] <- model$Gsq
  model_selection[i,4] <- model$resid.df
  model_selection[i,6] <- model$aic
  model_selection[i,7] <- model$bic
  model_selection[i,8] <- model$llik
  
  boot_gsquare <- rep(0, 100)
  for(j in 1:100){
    set.seed(j)
    iclass <- sample(1:(i+1), nrow(data_sat), replace = T, prob = model$P)
    response <- matrix(0, nrow(data_sat), 5)
    
    for (k in 1:nrow(data_sat)){
      for (l in 1:5){
        response[k, l] <- sample(1:3, 1, prob = model$probs[[l]][iclass[k],])
      }
    }
    boot_data <- as.data.frame(response)
    colnames(boot_data) <- paste0('q', 1:5)
    boot_model <- poLCA(f, data = boot_data, nclass = i+1,
                        na.rm = F, verbose = F)
    boot_gsquare[j] <- boot_model$Gsq
  }
  model_selection[i,5] <- mean(boot_gsquare >= model$Gsq)
  
  # identification check
  idf_check[[i]] <- table(round(model$attempts, 1))
}

model_selection
idf_check
sapply(lapply(idf_check, max), '/', 50)

windows(10, 10)
par(mfrow = c(2, 3))
par(mar = c(3, 5, 2 ,2))
for(i in 1:5){
  barplot(idf_check[[i]], horiz = T, las = 1, 
          main = paste('frequeancy of maximum log-likelihood (# of class = ', i+1, ')', sep = ''))
}

set.seed(2021021352)
f.1 <- cbind(q1, q2, q3, q4, q5) ~ gender+race+political+religion
model.sat <- poLCA(f.1, data = data_sat, nclass = 4,
                   na.rm = F, verbose = F, nrep = 200, maxiter = 10000)
model.sat

round(exp(model.sat$coeff),3)

# drug model fitting
f <- cbind(d1, d2, d3, d4)~1
idf_check <- list('nclass = 2' = NULL, 
                  'nclass = 3' = NULL, 'nclass = 4' = NULL, 
                  'nclass = 5' = NULL, 'nclass = 6' = NULL)
model_selection <- matrix(0, 5, 8)
colnames(model_selection) <- c('nclass','npar','G^2', 'df', 'p-value',
                               'AIC','BIC','likelihood')

for(i in 1:5){
  set.seed(2021021352)
  model_selection[i,1] <- i+1
  model <- poLCA(f, data = data_sat[,6:9], nclass = i+1, na.rm = F, verbose = F, nrep = 50)
  model_selection[i,2] <- model$npar
  model_selection[i,3] <- model$Gsq
  model_selection[i,4] <- model$resid.df
  model_selection[i,6] <- model$aic
  model_selection[i,7] <- model$bic
  model_selection[i,8] <- model$llik
  
  boot_gsquare <- rep(0, 100)
  for(j in 1:100){
    set.seed(j)
    iclass <- sample(1:(i+1), nrow(data_sat), replace = T, prob = model$P)
    response <- matrix(0, nrow(data_sat), 4)
    
    for (k in 1:nrow(data_sat)){
      for (l in 1:4){
        response[k, l] <- sample(1:3, 1, prob = model$probs[[l]][iclass[k],])
      }
    }
    boot_data <- as.data.frame(response)
    colnames(boot_data) <- paste0('d', 1:4)
    boot_model <- poLCA(f, data = boot_data, nclass = i+1,
                        na.rm = F, verbose = F)
    boot_gsquare[j] <- boot_model$Gsq
  }
  model_selection[i,5] <- mean(boot_gsquare >= model$Gsq)
  
  # identification check
  idf_check[[i]] <- table(round(model$attempts, 1))
}

model_selection
idf_check
sapply(lapply(idf_check, max), '/', 50)

windows(10, 10)
par(mfrow = c(2, 3))
par(mar = c(3, 5, 2 ,2))
for(i in 1:5){
  barplot(idf_check[[i]], horiz = T, las = 1, 
          main = paste('frequeancy of maximum log-likelihood (# of class = ', i+1, ')', sep = ''))
}

set.seed(2021021352)
f.1 <- cbind(d1, d2, d3, d4) ~ gender+race+political+religion
model.drug <- poLCA(f.1, data = data_sat, nclass = 3,
                   na.rm = F, verbose = F, nrep = 200, maxiter = 10000)

probs.start.new <- poLCA.reorder(model.drug$probs.start, c(2, 1, 3))

set.seed(2021021352) # baseline 변경을 위해 모형 다시 적합
model.drug.1 <- poLCA(f.1, data = data_sat, nclass = 3,
                      na.rm = F, verbose = F, nrep = 200, maxiter = 10000,
                      probs.start=probs.start.new)
model.drug.1
round(exp(model.drug.1$coeff),3)

#### 5. CMH 검정 ####
data_sat$sat.class <- model.sat$predclass
data_sat$drug.class <- model.drug.1$predclass

CMH.test.sat <- list()
comb <- combn(4, 2)
for(i in 1:ncol(comb)){
  comb.var <- c(comb[,i]+9, 14)
  data <- data_sat[,comb.var]
  
  CMH.test.sat[[i]] <- mantelhaen.test(table(data), correct = FALSE)
  names(CMH.test.sat)[i] <- paste(colnames(data)[1], '&', colnames(data)[2])
}

CMH.test.sat


CMH.test.drug <- list()
comb <- combn(4, 2)
for(i in 1:ncol(comb)){
  comb.var <- c(comb[,i]+9, 15)
  data <- data_sat[,comb.var]
  
  CMH.test.drug[[i]] <- mantelhaen.test(table(data), correct = FALSE)
  names(CMH.test.drug)[i] <- paste(colnames(data)[1], '&', colnames(data)[2])
}

CMH.test.drug


#### 6. 삶의 만족도와 약물 사용의 관계 ####
class.ind <- chisq.test(table(data_sat[,c('sat.class', 'drug.class')]), correct = FALSE)
class.ind
round(class.ind$stdres, 3)


