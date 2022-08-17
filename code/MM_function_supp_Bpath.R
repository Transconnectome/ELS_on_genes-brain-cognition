CP_Moderated_Mediation <- function(nihtbx,target,n_boot){
  
  setwd('/Users/wangheehwan/Desktop/COGGPS_PAPER/only_EU')
  pheno<-read.csv(file='final_phenotype_for_MM_european.csv')
  pheno<-fastDummies::dummy_cols(.data = pheno,select_columns =c("married","abcd_site") )
  cov_list<-c("age","sex","high.educ" ,"income", "BMI"  ) #remove sex when sex classified modeling
  cov2<-colnames(pheno)[43:84]
  
  married_cov<-cov2[2:6] ##1st(married_1) is reference
  site_cov<-append(cov2[7:21],cov2[23:27]) ##16th(abcd.site_16) is reference
  
  
  cov_list<-append(cov_list,married_cov)
  cov_list<-append(cov_list,site_cov)
  
  
  
  ##merge brain data and phenotype data
  CP_brain<-read.csv(file='CP_PCA_transformed.csv')
  CP_brain<-subset(CP_brain,select=c('subjectkey','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'))   ######brain target
  CP_final<-merge(pheno,CP_brain,by='subjectkey')
  ##check remaining NA
  CP_final<-na.omit(CP_final)
  #CP_final<-CP_final[CP_final$sex==1,]
  #CP_final<-CP_final[CP_final$sex==2,]
  
  zmInteraction <-CP_final$PC1*CP_final[,target]
  new_CP_final <- cbind(CP_final,zmInteraction)
  
  
  
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
  
  # covariate interaction
  for (cov_i in cov_list){
    covInteraction <- new_CP_final$CP*new_CP_final[,cov_i] 
    new_CP_final <- cbind(new_CP_final,covInteraction)
    names(new_CP_final)[length(new_CP_final)] <- paste0(cov_i,'.CP')
  }

  
  CPmodel <-"
  # main
    PC1 ~ a1*CP \n
    nihtbx ~ b1*PC1 \n
    nihtbx ~ cp*CP \n
    nihtbx ~ ELS*target \n
    nihtbx ~ interaction*zmInteraction \n
  
  
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

    
    #interaction covariates 
    #PC1 ~ d21*age.CP \n
    PC1 ~ d22*sex.CP \n
    PC1 ~ d23*high.educ.CP \n
    PC1 ~ d24*income.CP \n
    #PC1 ~ d25*BMI.CP \n
    #PC1 ~ d26*married_2.CP \n
    #PC1 ~ d27*married_3.CP \n
    #PC1 ~ d28*married_4.CP \n
    #PC1 ~ d29*married_5.CP \n
    #PC1 ~ d210*married_6.CP \n
    #PC1 ~ d211*abcd_site_1.CP \n
    #PC1 ~ d212*abcd_site_2.CP \n
    #PC1 ~ d213*abcd_site_3.CP \n
    #PC1 ~ d214*abcd_site_4.CP \n
    #PC1 ~ d215*abcd_site_5.CP \n
    #PC1 ~ d216*abcd_site_6.CP \n
    #PC1 ~ d217*abcd_site_7.CP \n
    #PC1 ~ d218*abcd_site_8.CP \n
    #PC1 ~ d219*abcd_site_9.CP \n
    #PC1 ~ d220*abcd_site_10.CP \n
    #PC1 ~ d221*abcd_site_11.CP \n
    #PC1 ~ d222*abcd_site_12.CP \n
    #PC1 ~ d223*abcd_site_13.CP \n
    #PC1 ~ d224*abcd_site_14.CP \n
    #PC1 ~ d225*abcd_site_15.CP \n
    #PC1 ~ d226*abcd_site_17.CP \n
    #PC1 ~ d227*abcd_site_18.CP \n
    #PC1 ~ d228*abcd_site_19.CP \n
    #PC1 ~ d229*abcd_site_20.CP \n
    #PC1 ~ d230*abcd_site_21.CP \n

    
    #nihtbx ~ f21*age.CP \n
    #nihtbx ~ f22*sex.CP \n
    #nihtbx ~ f23*high.educ.CP \n
    #nihtbx ~ f24*income.CP \n
    #nihtbx ~ f25*BMI.CP \n
    #nihtbx ~ f26*married_2.CP \n
    #nihtbx ~ f27*married_3.CP \n
    #nihtbx ~ f28*married_4.CP \n
    #nihtbx ~ f29*married_5.CP \n
    #nihtbx ~ f210*married_6.CP \n
    #nihtbx ~ f211*abcd_site_1.CP \n
    #nihtbx ~ f212*abcd_site_2.CP \n
    #nihtbx ~ f213*abcd_site_3.CP \n
    #nihtbx ~ f214*abcd_site_4.CP \n
    #nihtbx ~ f215*abcd_site_5.CP \n
    #nihtbx ~ f216*abcd_site_6.CP \n
    #nihtbx ~ f217*abcd_site_7.CP \n
    #nihtbx ~ f218*abcd_site_8.CP \n
    #nihtbx ~ f219*abcd_site_9.CP \n
    #nihtbx ~ f220*abcd_site_10.CP \n
    #nihtbx ~ f221*abcd_site_11.CP \n
    #nihtbx ~ f222*abcd_site_12.CP \n
    #nihtbx ~ f223*abcd_site_13.CP \n
    #nihtbx ~ f224*abcd_site_14.CP \n
    #nihtbx ~ f225*abcd_site_15.CP \n
    #nihtbx ~ f226*abcd_site_17.CP \n
    #nihtbx ~ f227*abcd_site_18.CP \n
    #nihtbx ~ f228*abcd_site_19.CP \n
    #nihtbx ~ f229*abcd_site_20.CP \n
    #nihtbx ~ f230*abcd_site_21.CP \n


    
    
    #summary 
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
  
  
  return(CP_fit)
}









EA_Moderated_Mediation <- function(nihtbx,target,n_boot){
  
  setwd('/Users/wangheehwan/Desktop/COGGPS_PAPER/only_EU')
  pheno<-read.csv(file='final_phenotype_for_MM_european.csv')
  pheno<-fastDummies::dummy_cols(.data = pheno,select_columns =c("married","abcd_site") )
  cov_list<-c("age","sex","high.educ" ,"income", "BMI" ) #remove sex when sex classified modeling
  #cov_list<-c("age","high.educ" ,"income", "BMI")
  cov2<-colnames(pheno)[43:84]
  
  married_cov<-cov2[2:6] ##1st(married_1) is reference
  site_cov<-append(cov2[7:21],cov2[23:27]) ##16th(abcd.site_16) is reference
  
  
  cov_list<-append(cov_list,married_cov)
  cov_list<-append(cov_list,site_cov)
  
  
  EA_brain<-read.csv(file='EA_PCA_transformed.csv')
  EA_brain<-subset(EA_brain,select=c('subjectkey','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'))   ######brain target
  EA_final<-merge(pheno,EA_brain,by='subjectkey')
  ###check remaining NA
  EA_final<-na.omit(EA_final)
  
  #EA_final<-EA_final[EA_final$sex==0,]
  #EA_final<-EA_final[EA_final$sex==1,]
  
  
  
  
  
  #mediator_CP<-'PC1'
  #mediator_EA<-'PC1'
  
  
  
  
  zmInteraction <-EA_final$PC1*EA_final[,target]
  new_EA_final <- cbind(EA_final,zmInteraction)
  
  
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
  
  # covariate interaction
  for (cov_i in cov_list){
    covInteraction <- new_EA_final$EA*new_EA_final[,cov_i] 
    new_EA_final <- cbind(new_EA_final,covInteraction)
    names(new_EA_final)[length(new_EA_final)] <- paste0(cov_i,'.EA')
  }

  
  EAmodel <-"
  # main
    PC1 ~ a1*EA \n
    nihtbx ~ b1*PC1 \n
    nihtbx ~ cp*EA \n
    nihtbx ~ ELS*target \n
    nihtbx ~ interaction*zmInteraction \n
  
  # covariate
    PC1 ~ d11*age \n
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

  
    nihtbx ~ f1*age \n
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


    #interaction covariates 
    PC1 ~ d22*sex.EA \n
    PC1 ~ d23*high.educ.EA \n
    PC1 ~ d24*income.EA \n



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

  return (EA_fit)
}

