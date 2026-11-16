###################
#  packages
####################
library(haven)
library(dplyr)
library(tidyr)
library(survey)
library(tidyverse) 
library(MASS)
library(brant)
library(broom)
library(ggplot2)
library(ggrepel)
library(caret)
library(VIM)
library(glmnet)


####################
# regression
###################

#model 1

fit<-coxph(Surv(survival_VTE, VTE)~as.factor(isolation_group)+age+Sex,
           data = data)

#model 2

fit<-coxph(Surv(survival_VTE, VTE)~as.factor(isolation_group)
           +age+Sex+BMI+education1+employment+town_group+smoke0+drink0+PA+cancer+cvd,
           data = data)

#model 3

fit<-coxph(Surv(survival_VTE, VTE)~as.factor(isolation_group)
           +age+Sex+BMI+education1+employment+town_group+smoke0+drink0+PA+cancer+cvd
           +PRS+genoarray+n_22009_0_1+n_22009_0_2+n_22009_0_3+   
             n_22009_0_4+n_22009_0_5+n_22009_0_6+n_22009_0_7+n_22009_0_8+n_22009_0_9+n_22009_0_10,
           data = data)

#model 4

fit<-coxph(Surv(survival_VTE, VTE)~as.factor(isolation_group)
           +age+Sex+BMI+education1+employment+town_group+smoke0+drink0+PA+cancer+cvd
           +PRS+genoarray+n_22009_0_1+n_22009_0_2+n_22009_0_3+   
             n_22009_0_4+n_22009_0_5+n_22009_0_6+n_22009_0_7+n_22009_0_8+n_22009_0_9+n_22009_0_10+as.factor(lonely_group),
           data = data)




################
#   PAF
################

fit<-coxph(Surv(survival_VTE, VTE)~as.factor(isolation_group)
           +age+Sex+BMI+education1+employment+town_group+smoke0+drink0+PA+cancer+cvd
           +PRS+genoarray+n_22009_0_1+n_22009_0_2+n_22009_0_3+   
             n_22009_0_4+n_22009_0_5+n_22009_0_6+n_22009_0_7+n_22009_0_8+n_22009_0_9+n_22009_0_10,
           data = data)


mf <- model.frame(fit)
idx <- match(rownames(mf), rownames(data))
stopifnot(!anyNA(idx))
analysis_data <- data[idx, , drop = FALSE]


calc_PAF <- function(d) {
  model <- coxph(
    formula(fit),
    data = d
  )
  HR <- round(exp(coef(model)[
    c("as.factor(isolation_group)2",
      "as.factor(isolation_group)3")
  ]), 2)
  
  if (any(!is.finite(HR)) || any(HR <= 0)) {
    return(NA_real_)
  }
  
 
case_prop <- prop.table(
    table(factor(
      d$isolation_group[d$VTE == 1],
      levels = c(1, 2, 3))))

  PAF <- sum(case_prop[c("2", "3")] * (HR - 1) / HR)
  
  return(PAF)
}

PAF_original <- calc_PAF(analysis_data)



#############################
#  proteins_association
#############################

result_protein <- data.frame()
for (i in 2:1460) {
  
  
  fit <- lm(all[,i]~isolation0
            +age+Sex+BMI+education1+employment+town_group+smoke0+drink0+PA+cancer+cvd+renal_egfr
            ,data = all)
  fit <- summary(fit)
  beta <-round(fit$coefficients[2,1],4)
  mul_PValue <- format(fit$coefficients[,4]) 
  
  LL<-round(fit$coefficients[,1]-1.96*fit$coefficients[,2],4)
  UL<-round(fit$coefficients[,1]+1.96*fit$coefficients[,2],4)
  mul_CI95<-paste(LL,',',UL)
  mul_cox1 <- data.frame("beta" =beta,
                         "CI95" =mul_CI95,
                         "P"=mul_PValue)
  
  mul_cox1$b_95 <- paste(mul_cox1$beta,'(', mul_cox1$CI95,')', sep = '')
  
  mul_cox1$P1<-mul_cox1$P
  mul_cox1 <- mul_cox1[,c(1,3,4)]
  mul_cox1 <- mul_cox1[2,]
  
  result_protein <- rbind(result_protein,mul_cox1)
  
}


b <- p.adjust(result_protein$P,method = "BH") 
result_protein$padjust <- NA
result_protein$padjust <- b

print(result_protein)




############################
#   proteomic signature
############################

record_ISO = data.frame(repNo = NA, NoMetabs = NA)  
pro_coef_ISO = data.frame(HMDB = colnames(train_imputed1[, c(1:1460)]))  
library(glmnet)


for (inrep in 1:100) {   
  
  x = as.matrix(train_imputed1[, c(2:1460)])  
  y = train_imputed1$isolation0 
  
   
  train_imputed1ing_CV = cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = 10, standardize = FALSE, parallel = TRUE)  
  
  lambda_1se_10F = train_imputed1ing_CV$lambda.1se  
  
  current_coef = as.data.frame(as.matrix(coef(train_imputed1ing_CV, s = lambda_1se_10F)))   
  
  pro_coef_ISO[, dim(pro_coef_ISO)[2] + 1] = data.frame(current_coef)  
  names(pro_coef_ISO)[dim(pro_coef_ISO)[2]] = paste("inrep", inrep, sep='')  
  
  record_ISO[inrep, "repNo"] = paste("inrep", inrep, sep='')  
  record_ISO[inrep, "NoMetabs"] = sum(pro_coef_ISO[, dim(pro_coef_ISO)[2]] != 0)  
  
  print(record_ISO[inrep, "repNo"])  
}  
pro_coef_ISO <- pro_coef_ISO[-1, ]  

pro_coef_ISO$inrep <- apply(pro_coef_ISO[, 2:101], 1, function(row) {  
  if (sum(row == 0) >= 50) {  
    return(0)  
  } else {  
    return(mean(row, na.rm = TRUE) ) 
  }  
}) 

sum(pro_coef_ISO$inrep != 0) 
table(record_ISO$NoMetabs) 


identical(colnames(test_imputed1)[2:1460], pro_coef_ISO$HMDB) 
test_imputed1$proscore_ISO=colSums(t(test_imputed1[, 2:1460])*pro_coef_ISO$inrep)

summary(test_imputed1$proscore_ISO)

save(pro_coef_ISO, file = "pro_coef_ISO.RData")  
save(record_ISO, file = "record_ISO.RData")  



############################
#   mediation analysis
############################


test_imputed1$proscore_ISO_group <- as.ordered(test_imputed1$proscore_ISO_group)

exp1 <- cmest(data = test_imputed1, model = "rb", outcome = "survival_VTE",
              event = "VTE",
              exposure = "isolation0", mediator = "proscore_ISO_group",
              basec = c("age", "Sex","education1","employment","BMI",
                        "town_index","smoke0","drink0","PA","cvd","cancer"), 
              EMint = FALSE, mreg = list("ordinal"), yreg = "coxph", 
              astar = 0, a = 1, mval = list("3"), estimation = "imputation", 
              inference = "bootstrap",nboot=5000)  
summary(exp1)




