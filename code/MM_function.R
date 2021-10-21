CP_Moderated_Mediation <- function(nihtbx,target,n_boot){
  
  setwd('/work2/07939/tg872382/stampede2/connectome/stampede2/cognitive_GPS')
  pheno<-read.csv(file='final_phenotype_for_MM.csv')
  pheno<-fastDummies::dummy_cols(.data = pheno,select_columns =c("married","abcd_site",'race.ethnicity') )
  cov_list<-c("age","sex","high.educ" ,"income", "BMI"  ) #remove sex when sex classified modeling
  cov2<-colnames(pheno)[53:84]
  
  married_cov<-cov2[2:6] ##1st(married_1) is reference
  site_cov<-append(cov2[7:21],cov2[23:27]) ##16th(abcd.site_16) is reference
  race_cov<-cov2[29:32] ##1st(race.ethnicity_1; European ancestry) is reference
  
  cov_list<-append(cov_list,married_cov)
  cov_list<-append(cov_list,site_cov)
  cov_list<-append(cov_list,race_cov)
  
  
  
  ##merge brain data and phenotype data
  CP_brain<-read.csv(file='CP_PCA_transformed.csv')
  CP_brain<-subset(CP_brain,select=c('subjectkey','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'))   ######brain target
  CP_final<-merge(pheno,CP_brain,by='subjectkey')
  ##check remaining NA
  CP_final<-na.omit(CP_final)
  #CP_final<-CP_final[CP_final$sex==1,]
  #CP_final<-CP_final[CP_final$sex==2,]
  
  xmInteraction <-CP_final$CP*CP_final[,target]
  new_CP_final <- cbind(CP_final,xmInteraction)
  
  col_list<-colnames(new_CP_final)
  for (i in 1:length(col_list)){
    if (col_list[i] == target){
      names(new_CP_final)[i] <- 'target'
    }  
  }
  
  col_list<-colnames(new_CP_final)
  for (i in 1:length(col_list)){
    if (col_list[i] == nihtbx){
      names(new_CP_final)[i] <- 'nihtbx'
    }  
  }
  
  CPmodel <-"
  # main
    PC1 ~ a1*CP \n
    nihtbx ~ b1*PC1 \n
    nihtbx ~ cp*CP \n
    PC1 ~ ELS*target \n
    PC1 ~ interaction*xmInteraction \n
  
  
  # covariate
    PC1 ~ d11*age \n
    PC1 ~ d12*sex \n
    PC1 ~ d13*high.educ \n
    PC1 ~ d14*income \n
    PC1 ~ d15*BMI \n
    PC1 ~ d16*married_2 \n
    PC1 ~ d17*married_3 \n
    PC1 ~ d18*married_4 \n
    PC1 ~ d19*married_5 \n
    PC1 ~ d110*married_6 \n
    PC1 ~ d111*abcd_site_1 \n
    PC1 ~ d112*abcd_site_2 \n
    PC1 ~ d113*abcd_site_3 \n
    PC1 ~ d114*abcd_site_4 \n
    PC1 ~ d115*abcd_site_5 \n
    PC1 ~ d116*abcd_site_6 \n
    PC1 ~ d117*abcd_site_7 \n
    PC1 ~ d118*abcd_site_8 \n
    PC1 ~ d119*abcd_site_9 \n
    PC1 ~ d120*abcd_site_10 \n
    PC1 ~ d121*abcd_site_11 \n
    PC1 ~ d122*abcd_site_12 \n
    PC1 ~ d123*abcd_site_13 \n
    PC1 ~ d124*abcd_site_14 \n
    PC1 ~ d125*abcd_site_15 \n
    PC1 ~ d126*abcd_site_17 \n
    PC1 ~ d127*abcd_site_18 \n
    PC1 ~ d128*abcd_site_19 \n
    PC1 ~ d129*abcd_site_20 \n
    PC1 ~ d130*abcd_site_21 \n
    PC1 ~ d111*race.ethnicity_2 \n
    PC1 ~ d112*race.ethnicity_3 \n
    PC1 ~ d113*race.ethnicity_4 \n
    PC1 ~ d114*race.ethnicity_5 \n
  
  
    nihtbx ~ f1*age \n
    nihtbx ~ f2*sex \n
    nihtbx ~ f3*high.educ \n
    nihtbx ~ f4*income \n
    nihtbx ~ f5*BMI \n
    nihtbx ~ f6*married_2 \n
    nihtbx ~ f7*married_3 \n
    nihtbx ~ f8*married_4 \n
    nihtbx ~ f9*married_5 \n
    nihtbx ~ f10*married_6 \n
    nihtbx ~ f11*abcd_site_1 \n
    nihtbx ~ f12*abcd_site_2 \n
    nihtbx ~ f13*abcd_site_3 \n
    nihtbx ~ f14*abcd_site_4 \n
    nihtbx ~ f15*abcd_site_5 \n
    nihtbx ~ f16*abcd_site_6 \n
    nihtbx ~ f17*abcd_site_7 \n
    nihtbx ~ f18*abcd_site_8 \n
    nihtbx ~ f19*abcd_site_9 \n
    nihtbx ~ f20*abcd_site_10 \n
    nihtbx ~ f21*abcd_site_11 \n
    nihtbx ~ f22*abcd_site_12 \n
    nihtbx ~ f23*abcd_site_13 \n
    nihtbx ~ f24*abcd_site_14 \n
    nihtbx ~ f25*abcd_site_15 \n
    nihtbx ~ f26*abcd_site_17 \n
    nihtbx ~ f27*abcd_site_18 \n
    nihtbx ~ f28*abcd_site_19 \n
    nihtbx ~ f29*abcd_site_20 \n
    nihtbx ~ f30*abcd_site_21 \n
    nihtbx ~ f31*race.ethnicity_2 \n
    nihtbx ~ f32*race.ethnicity_3 \n
    nihtbx ~ f33*race.ethnicity_4 \n
    nihtbx ~ f34*race.ethnicity_5 \n
    indirect := a1 * b1 \n 
    IMM := interaction * b1 \n  
    bw1 := ELS * b1 \n 
    gw1 := ELS * interaction \n 
    ratio1 := (indirect)/(indirect+cp) \n 
    ratio_tot := (indirect)/(indirect+cp) \n 
  "
  
  
  
  #result
  CP_fit <- lavaan::sem(CPmodel, data = new_CP_final,fixed.x = TRUE,std.lv = TRUE,se = 'bootstrap',bootstrap=n_boot,parallel = c("multicore"))
  #CP_result <- lavaan::standardizedSolution(CP_fit)
  #save.image('MM_cryst_with_3PC.RData')
  #load('MM_cryst_with_3PC.RData')

  #save.image('MM_cryst_with_3PC.RData')
  
  
  return(CP_fit )
  }
  








EA_Moderated_Mediation <- function(nihtbx,target,n_boot){
  
  setwd('/work2/07939/tg872382/stampede2/connectome/stampede2/cognitive_GPS')
  pheno<-read.csv(file='final_phenotype_for_MM.csv')
  pheno<-fastDummies::dummy_cols(.data = pheno,select_columns =c("married","abcd_site",'race.ethnicity') )
  cov_list<-c("age","sex","high.educ" ,"income", "BMI"  ) #remove sex when sex classified modeling
  cov2<-colnames(pheno)[53:84]
  
  married_cov<-cov2[2:6] ##1st(married_1) is reference
  site_cov<-append(cov2[7:21],cov2[23:27]) ##16th(abcd.site_16) is reference
  race_cov<-cov2[29:32] ##1st(race.ethnicity_1; European ancestry) is reference
  
  cov_list<-append(cov_list,married_cov)
  cov_list<-append(cov_list,site_cov)
  cov_list<-append(cov_list,race_cov)
  
  EA_brain<-read.csv(file='EA_PCA_transformed.csv')
  EA_brain<-subset(EA_brain,select=c('subjectkey','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'))   ######brain target
  EA_final<-merge(pheno,EA_brain,by='subjectkey')
  ###check remaining NA
  EA_final<-na.omit(EA_final)
  #EA_final<-EA_final[EA_final$sex==1,]
  #EA_final<-EA_final[EA_final$sex==2,]
  
  
  
  
  #mediator_CP<-'PC1'
  #mediator_EA<-'PC1'
  
  
  
  
  xmInteraction <-EA_final$EA*EA_final[,target]
  new_EA_final <- cbind(EA_final,xmInteraction)
  
  col_list<-colnames(new_EA_final)
  for (i in 1:length(col_list)){
    if (col_list[i] == target){
      names(new_EA_final)[i] <- 'target'
    }  
  }

  col_list<-colnames(new_EA_final)
  for (i in 1:length(col_list)){
    if (col_list[i] == nihtbx){
      names(new_EA_final)[i] <- 'nihtbx'
    }  
  }
  
  EAmodel <-"
  # main
    PC1 ~ a1*EA \n
    nihtbx ~ b1*PC1 \n
    nihtbx ~ cp*EA \n
    PC1 ~ ELS*target \n
    PC1 ~ interaction*xmInteraction \n
  
  # covariate
    PC1 ~ d11*age \n
    PC1 ~ d12*sex \n
    PC1 ~ d13*high.educ \n
    PC1 ~ d14*income \n
    PC1 ~ d15*BMI \n
    PC1 ~ d16*married_2 \n
    PC1 ~ d17*married_3 \n
    PC1 ~ d18*married_4 \n
    PC1 ~ d19*married_5 \n
    PC1 ~ d110*married_6 \n
    PC1 ~ d111*abcd_site_1 \n
    PC1 ~ d112*abcd_site_2 \n
    PC1 ~ d113*abcd_site_3 \n
    PC1 ~ d114*abcd_site_4 \n
    PC1 ~ d115*abcd_site_5 \n
    PC1 ~ d116*abcd_site_6 \n
    PC1 ~ d117*abcd_site_7 \n
    PC1 ~ d118*abcd_site_8 \n
    PC1 ~ d119*abcd_site_9 \n
    PC1 ~ d120*abcd_site_10 \n
    PC1 ~ d121*abcd_site_11 \n
    PC1 ~ d122*abcd_site_12 \n
    PC1 ~ d123*abcd_site_13 \n
    PC1 ~ d124*abcd_site_14 \n
    PC1 ~ d125*abcd_site_15 \n
    PC1 ~ d126*abcd_site_17 \n
    PC1 ~ d127*abcd_site_18 \n
    PC1 ~ d128*abcd_site_19 \n
    PC1 ~ d129*abcd_site_20 \n
    PC1 ~ d130*abcd_site_21 \n
    PC1 ~ d111*race.ethnicity_2 \n
    PC1 ~ d112*race.ethnicity_3 \n
    PC1 ~ d113*race.ethnicity_4 \n
    PC1 ~ d114*race.ethnicity_5 \n
  
    nihtbx ~ f1*age \n
    nihtbx ~ f2*sex \n
    nihtbx ~ f3*high.educ \n
    nihtbx ~ f4*income \n
    nihtbx ~ f5*BMI \n
    nihtbx ~ f6*married_2 \n
    nihtbx ~ f7*married_3 \n
    nihtbx ~ f8*married_4 \n
    nihtbx ~ f9*married_5 \n
    nihtbx ~ f10*married_6 \n
    nihtbx ~ f11*abcd_site_1 \n
    nihtbx ~ f12*abcd_site_2 \n
    nihtbx ~ f13*abcd_site_3 \n
    nihtbx ~ f14*abcd_site_4 \n
    nihtbx ~ f15*abcd_site_5 \n
    nihtbx ~ f16*abcd_site_6 \n
    nihtbx ~ f17*abcd_site_7 \n
    nihtbx ~ f18*abcd_site_8 \n
    nihtbx ~ f19*abcd_site_9 \n
    nihtbx ~ f20*abcd_site_10 \n
    nihtbx ~ f21*abcd_site_11 \n
    nihtbx ~ f22*abcd_site_12 \n
    nihtbx ~ f23*abcd_site_13 \n
    nihtbx ~ f24*abcd_site_14 \n
    nihtbx ~ f25*abcd_site_15 \n
    nihtbx ~ f26*abcd_site_17 \n
    nihtbx ~ f27*abcd_site_18 \n
    nihtbx ~ f28*abcd_site_19 \n
    nihtbx ~ f29*abcd_site_20 \n
    nihtbx ~ f30*abcd_site_21 \n
    nihtbx ~ f31*race.ethnicity_2 \n
    nihtbx ~ f32*race.ethnicity_3 \n
    nihtbx ~ f33*race.ethnicity_4 \n
    nihtbx ~ f34*race.ethnicity_5 \n
    indirect := a1 * b1 \n 
    IMM := interaction * b1 \n  
    bw1 := ELS * b1 \n 
    gw1 := ELS * interaction \n 
    ratio1 := (indirect)/(indirect+cp) \n 
    ratio_tot := (indirect)/(indirect+cp) \n 
  
  "

  #result
  #CP_fit <- lavaan::sem(CPmodel, data = new_CP_final,fixed.x = TRUE,std.lv = TRUE,se = 'bootstrap',bootstrap=1000)
  #CP_result <- lavaan::standardizedSolution(CP_fit)
  #save.image('MM_cryst_with_3PC.RData')
  #load('MM_cryst_with_3PC.RData')
  EA_fit <- lavaan::sem(EAmodel, data = new_EA_final,fixed.x = TRUE,std.lv = TRUE,se = 'bootstrap',bootstrap=n_boot,parallel = c("multicore"))
  
  #EA_result <- lavaan::standardizedSolution(EA_fit)
  
  #save.image('MM_cryst_with_3PC.RData')
  return (EA_fit)
}



CP_Mediated_Moderation <- function(nihtbx,target,n_boot){
  
  setwd('/work2/07939/tg872382/stampede2/connectome/stampede2/cognitive_GPS')
  pheno<-read.csv(file='final_phenotype_for_MM.csv')
  pheno<-fastDummies::dummy_cols(.data = pheno,select_columns =c("married","abcd_site",'race.ethnicity') )
  cov_list<-c("age","sex","high.educ" ,"income", "BMI"  ) #remove sex when sex classified modeling
  cov2<-colnames(pheno)[53:84]
  
  married_cov<-cov2[2:6] ##1st(married_1) is reference
  site_cov<-append(cov2[7:21],cov2[23:27]) ##16th(abcd.site_16) is reference
  race_cov<-cov2[29:32] ##1st(race.ethnicity_1; European ancestry) is reference
  
  cov_list<-append(cov_list,married_cov)
  cov_list<-append(cov_list,site_cov)
  cov_list<-append(cov_list,race_cov)
  
  
  
  ##merge brain data and phenotype data
  CP_brain<-read.csv(file='CP_PCA_transformed.csv')
  CP_brain<-subset(CP_brain,select=c('subjectkey','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'))   ######brain target
  CP_final<-merge(pheno,CP_brain,by='subjectkey')
  ##check remaining NA
  CP_final<-na.omit(CP_final)
  #CP_final<-CP_final[CP_final$sex==1,]
  #CP_final<-CP_final[CP_final$sex==2,]
  
  xmInteraction <-CP_final$CP*CP_final[,target]
  new_CP_final <- cbind(CP_final,xmInteraction)
  
  col_list<-colnames(new_CP_final)
  for (i in 1:length(col_list)){
    if (col_list[i] == target){
      names(new_CP_final)[i] <- 'target'
    }  
  }
  
  col_list<-colnames(new_CP_final)
  for (i in 1:length(col_list)){
    if (col_list[i] == nihtbx){
      names(new_CP_final)[i] <- 'nihtbx'
    }  
  }
  
  CPmodel <-"
  # main
    PC1 ~ a1*CP \n
    nihtbx ~ b1*PC1 \n
    nihtbx ~ cp*CP \n
    PC1 ~ ELS1*target \n
    nihtbx ~ ELS2*target \n
    PC1 ~ indirect_interaction*xmInteraction \n
    nihtbx ~ direct_interaction*xmInteraction \n
  
  
  # covariate
    PC1 ~ d11*age \n
    PC1 ~ d12*sex \n
    PC1 ~ d13*high.educ \n
    PC1 ~ d14*income \n
    PC1 ~ d15*BMI \n
    PC1 ~ d16*married_2 \n
    PC1 ~ d17*married_3 \n
    PC1 ~ d18*married_4 \n
    PC1 ~ d19*married_5 \n
    PC1 ~ d110*married_6 \n
    PC1 ~ d111*abcd_site_1 \n
    PC1 ~ d112*abcd_site_2 \n
    PC1 ~ d113*abcd_site_3 \n
    PC1 ~ d114*abcd_site_4 \n
    PC1 ~ d115*abcd_site_5 \n
    PC1 ~ d116*abcd_site_6 \n
    PC1 ~ d117*abcd_site_7 \n
    PC1 ~ d118*abcd_site_8 \n
    PC1 ~ d119*abcd_site_9 \n
    PC1 ~ d120*abcd_site_10 \n
    PC1 ~ d121*abcd_site_11 \n
    PC1 ~ d122*abcd_site_12 \n
    PC1 ~ d123*abcd_site_13 \n
    PC1 ~ d124*abcd_site_14 \n
    PC1 ~ d125*abcd_site_15 \n
    PC1 ~ d126*abcd_site_17 \n
    PC1 ~ d127*abcd_site_18 \n
    PC1 ~ d128*abcd_site_19 \n
    PC1 ~ d129*abcd_site_20 \n
    PC1 ~ d130*abcd_site_21 \n
    PC1 ~ d111*race.ethnicity_2 \n
    PC1 ~ d112*race.ethnicity_3 \n
    PC1 ~ d113*race.ethnicity_4 \n
    PC1 ~ d114*race.ethnicity_5 \n
  
  
    nihtbx ~ f1*age \n
    nihtbx ~ f2*sex \n
    nihtbx ~ f3*high.educ \n
    nihtbx ~ f4*income \n
    nihtbx ~ f5*BMI \n
    nihtbx ~ f6*married_2 \n
    nihtbx ~ f7*married_3 \n
    nihtbx ~ f8*married_4 \n
    nihtbx ~ f9*married_5 \n
    nihtbx ~ f10*married_6 \n
    nihtbx ~ f11*abcd_site_1 \n
    nihtbx ~ f12*abcd_site_2 \n
    nihtbx ~ f13*abcd_site_3 \n
    nihtbx ~ f14*abcd_site_4 \n
    nihtbx ~ f15*abcd_site_5 \n
    nihtbx ~ f16*abcd_site_6 \n
    nihtbx ~ f17*abcd_site_7 \n
    nihtbx ~ f18*abcd_site_8 \n
    nihtbx ~ f19*abcd_site_9 \n
    nihtbx ~ f20*abcd_site_10 \n
    nihtbx ~ f21*abcd_site_11 \n
    nihtbx ~ f22*abcd_site_12 \n
    nihtbx ~ f23*abcd_site_13 \n
    nihtbx ~ f24*abcd_site_14 \n
    nihtbx ~ f25*abcd_site_15 \n
    nihtbx ~ f26*abcd_site_17 \n
    nihtbx ~ f27*abcd_site_18 \n
    nihtbx ~ f28*abcd_site_19 \n
    nihtbx ~ f29*abcd_site_20 \n
    nihtbx ~ f30*abcd_site_21 \n
    nihtbx ~ f31*race.ethnicity_2 \n
    nihtbx ~ f32*race.ethnicity_3 \n
    nihtbx ~ f33*race.ethnicity_4 \n
    nihtbx ~ f34*race.ethnicity_5 \n
    indirect := a1 * b1 \n 
    IMM := indirect_interaction * b1 \n  
    bw1 := ELS1 * b1 \n 
    gw1 := ELS1 * indirect_interaction \n 
    ratio1 := (indirect)/(indirect+cp) \n 
    ratio_tot := (indirect)/(indirect+cp) \n 
  "
  
  
  
  #result
  CP_fit <- lavaan::sem(CPmodel, data = new_CP_final,fixed.x = TRUE,std.lv = TRUE,se = 'bootstrap',bootstrap=n_boot,parallel = c("multicore"))
  #CP_result <- lavaan::standardizedSolution(CP_fit)
  #save.image('MM_cryst_with_3PC.RData')
  #load('MM_cryst_with_3PC.RData')
  
  #save.image('MM_cryst_with_3PC.RData')
  
  
  return(CP_fit )
}









EA_Mediated_Moderation <- function(nihtbx,target,n_boot){
  
  setwd('/work2/07939/tg872382/stampede2/connectome/stampede2/cognitive_GPS')
  pheno<-read.csv(file='final_phenotype_for_MM.csv')
  pheno<-fastDummies::dummy_cols(.data = pheno,select_columns =c("married","abcd_site",'race.ethnicity') )
  cov_list<-c("age","sex","high.educ" ,"income", "BMI"  ) #remove sex when sex classified modeling
  cov2<-colnames(pheno)[53:84]
  
  married_cov<-cov2[2:6] ##1st(married_1) is reference
  site_cov<-append(cov2[7:21],cov2[23:27]) ##16th(abcd.site_16) is reference
  race_cov<-cov2[29:32] ##1st(race.ethnicity_1; European ancestry) is reference
  
  cov_list<-append(cov_list,married_cov)
  cov_list<-append(cov_list,site_cov)
  cov_list<-append(cov_list,race_cov)
  
  EA_brain<-read.csv(file='EA_PCA_transformed.csv')
  EA_brain<-subset(EA_brain,select=c('subjectkey','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'))   ######brain target
  EA_final<-merge(pheno,EA_brain,by='subjectkey')
  ###check remaining NA
  EA_final<-na.omit(EA_final)
  #EA_final<-EA_final[EA_final$sex==1,]
  #EA_final<-EA_final[EA_final$sex==2,]
  
  
  
  
  #mediator_CP<-'PC1'
  #mediator_EA<-'PC1'
  
  
  
  
  xmInteraction <-EA_final$EA*EA_final[,target]
  new_EA_final <- cbind(EA_final,xmInteraction)
  
  col_list<-colnames(new_EA_final)
  for (i in 1:length(col_list)){
    if (col_list[i] == target){
      names(new_EA_final)[i] <- 'target'
    }  
  }
  
  col_list<-colnames(new_EA_final)
  for (i in 1:length(col_list)){
    if (col_list[i] == nihtbx){
      names(new_EA_final)[i] <- 'nihtbx'
    }  
  }
  
  EAmodel <-"
  # main
    PC1 ~ a1*EA \n
    nihtbx ~ b1*PC1 \n
    nihtbx ~ cp*EA \n
    PC1 ~ ELS1*target \n
    nihtbx ~ ELS2*target \n
    PC1 ~ indirect_interaction*xmInteraction \n
    nihtbx ~ direct_interaction*xmInteraction \n
  
  # covariate
    PC1 ~ d11*age \n
    PC1 ~ d12*sex \n
    PC1 ~ d13*high.educ \n
    PC1 ~ d14*income \n
    PC1 ~ d15*BMI \n
    PC1 ~ d16*married_2 \n
    PC1 ~ d17*married_3 \n
    PC1 ~ d18*married_4 \n
    PC1 ~ d19*married_5 \n
    PC1 ~ d110*married_6 \n
    PC1 ~ d111*abcd_site_1 \n
    PC1 ~ d112*abcd_site_2 \n
    PC1 ~ d113*abcd_site_3 \n
    PC1 ~ d114*abcd_site_4 \n
    PC1 ~ d115*abcd_site_5 \n
    PC1 ~ d116*abcd_site_6 \n
    PC1 ~ d117*abcd_site_7 \n
    PC1 ~ d118*abcd_site_8 \n
    PC1 ~ d119*abcd_site_9 \n
    PC1 ~ d120*abcd_site_10 \n
    PC1 ~ d121*abcd_site_11 \n
    PC1 ~ d122*abcd_site_12 \n
    PC1 ~ d123*abcd_site_13 \n
    PC1 ~ d124*abcd_site_14 \n
    PC1 ~ d125*abcd_site_15 \n
    PC1 ~ d126*abcd_site_17 \n
    PC1 ~ d127*abcd_site_18 \n
    PC1 ~ d128*abcd_site_19 \n
    PC1 ~ d129*abcd_site_20 \n
    PC1 ~ d130*abcd_site_21 \n
    PC1 ~ d111*race.ethnicity_2 \n
    PC1 ~ d112*race.ethnicity_3 \n
    PC1 ~ d113*race.ethnicity_4 \n
    PC1 ~ d114*race.ethnicity_5 \n
  
    nihtbx ~ f1*age \n
    nihtbx ~ f2*sex \n
    nihtbx ~ f3*high.educ \n
    nihtbx ~ f4*income \n
    nihtbx ~ f5*BMI \n
    nihtbx ~ f6*married_2 \n
    nihtbx ~ f7*married_3 \n
    nihtbx ~ f8*married_4 \n
    nihtbx ~ f9*married_5 \n
    nihtbx ~ f10*married_6 \n
    nihtbx ~ f11*abcd_site_1 \n
    nihtbx ~ f12*abcd_site_2 \n
    nihtbx ~ f13*abcd_site_3 \n
    nihtbx ~ f14*abcd_site_4 \n
    nihtbx ~ f15*abcd_site_5 \n
    nihtbx ~ f16*abcd_site_6 \n
    nihtbx ~ f17*abcd_site_7 \n
    nihtbx ~ f18*abcd_site_8 \n
    nihtbx ~ f19*abcd_site_9 \n
    nihtbx ~ f20*abcd_site_10 \n
    nihtbx ~ f21*abcd_site_11 \n
    nihtbx ~ f22*abcd_site_12 \n
    nihtbx ~ f23*abcd_site_13 \n
    nihtbx ~ f24*abcd_site_14 \n
    nihtbx ~ f25*abcd_site_15 \n
    nihtbx ~ f26*abcd_site_17 \n
    nihtbx ~ f27*abcd_site_18 \n
    nihtbx ~ f28*abcd_site_19 \n
    nihtbx ~ f29*abcd_site_20 \n
    nihtbx ~ f30*abcd_site_21 \n
    nihtbx ~ f31*race.ethnicity_2 \n
    nihtbx ~ f32*race.ethnicity_3 \n
    nihtbx ~ f33*race.ethnicity_4 \n
    nihtbx ~ f34*race.ethnicity_5 \n
    indirect := a1 * b1 \n 
    IMM := indirect_interaction * b1 \n  
    bw1 := ELS1 * b1 \n 
    gw1 := ELS1 * indirect_interaction \n 
    ratio1 := (indirect)/(indirect+cp) \n 
    ratio_tot := (indirect)/(indirect+cp) \n 
  
  "
  
  #result
  #CP_fit <- lavaan::sem(CPmodel, data = new_CP_final,fixed.x = TRUE,std.lv = TRUE,se = 'bootstrap',bootstrap=1000)
  #CP_result <- lavaan::standardizedSolution(CP_fit)
  #save.image('MM_cryst_with_3PC.RData')
  #load('MM_cryst_with_3PC.RData')
  EA_fit <- lavaan::sem(EAmodel, data = new_EA_final,fixed.x = TRUE,std.lv = TRUE,se = 'bootstrap',bootstrap=n_boot,parallel = c("multicore"))
  
  #EA_result <- lavaan::standardizedSolution(EA_fit)
  
  #save.image('MM_cryst_with_3PC.RData')
  return (EA_fit)
}



