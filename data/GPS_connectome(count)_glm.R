### extracting connectome(count) brain features 
setwd("~/Desktop/COGGPS_PAPER/")


###data prep
connectome<-read.csv(file ='final_ver/con_aparc_count.csv')
pheno<-read.csv(file='only_EU/final_phenotype_for_MM_european.csv')

MERGED<-merge(pheno,connectome,by='subjectkey')
#MERGED<-na.omit(MERGED)

COV_LIST <- c('age', 'sex','high.educ', 'income', 'married', 'abcd_site', 'BMI') # should be checked
fact_COV_LIST <- c('age','sex' ,'high.educ','factor(married)','factor(abcd_site)', 'BMI') # should be checked
print(fact_COV_LIST)
COVARS <- paste(fact_COV_LIST, collapse=" + ")
GPS_LIST<-c('CP','EA')


names(MERGED)[43]
BRAIN_START = 43
BRAIN_END <- length(MERGED)
NUM_GPS <- length(GPS_LIST)

BETA_LIST <- c( 'beta_brain')
P_LIST <- c( 'p_brain')

###calculating
for (j in 1:2){
  for (gps in GPS_LIST[[j]]){
    print(paste('--------------',gps,'--------------'))
    # initialize all vector
    brain <- vector()
    
    glm.se <- vector()
    glm.beta <- data.frame()
    glm.p <- data.frame()
    
    # DO GLM
    for (i in BRAIN_START:BRAIN_END){
      print(i)
      brain <- c(brain, names(MERGED)[i])
      
      form <- paste(names(MERGED)[i], '~',gps,'+' ,COVARS)
      
      glm.re <- glm(as.formula(form), data=MERGED, family=gaussian())
      cff <- coef(summary(glm.re))
      glm.se <- c(glm.se, round(cff[2, 2], 5))
      glm.beta <- rbind.data.frame(glm.beta, cff[2,1])
      glm.p <- rbind.data.frame(glm.p, cff[2,4])
      
      
    }
    colnames(glm.beta) <- BETA_LIST
    colnames(glm.p) <- P_LIST
    result <- data.frame(brain, glm.beta, glm.p)
    result <- result[-which(duplicated(result$p_brain)),]
    FDR<-p.adjust(result$p_brain,method='fdr')
    result<-cbind(result,FDR)
    BONF<-p.adjust(result$p_brain,method='bonferroni')
    result<-cbind(result,BONF)
    FILE_NAME <- paste0('only_EU/DTI_resutls/GLM_con(count)_cov_', gps, '_prsice2_corrected.csv')
    write.csv(result, FILE_NAME, quote = F,row.names = F)
    
    
  }
  
} 


###list significant brain features


### making significant list 
brainList<-function(target){
  brain_list<-c()
  for (i in 1:nrow(target)){
    if (target[i,'BONF']<0.05){
      brain_list<-rbind(brain_list,target[i,])
    }
  }
  
  file_name<-paste0('',deparse(substitute(target)),'_brain_list.csv')
  
  write.csv(brain_list,file_name,row.names = F,quote = F)
}

CP_related<-read.csv('only_EU/GLM_con(count)_cov_CP_prsice2_corrected.csv')
EA_related<-read.csv('only_EU/GLM_con(count)_cov_EA_prsice2_corrected.csv')

brainList(CP_related)
brainList(EA_related)


