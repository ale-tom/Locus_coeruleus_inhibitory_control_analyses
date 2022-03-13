set.seed(123)#for reproducibility

pacman::p_load('tidyverse','magrittr','ggpubr','patchwork','gtsummary')

# retrieve data-----------------------------------------------------------------

# load model estimates from MRI session
ssrt <- as.data.frame(read.csv(here::here('Data','CC_exGaussian_estimates.csv')))

# load LC contrast measures
LC <-as.data.frame(read.csv(here::here('Data','CC_LC_CNR.csv')))
CNR      <-as.data.frame(read.csv(here::here('Data','CNR_roi.csv')))# LC contrast segregated by sub-regions
SNR      <- as.data.frame(read.csv(here::here('Data','LC_mean_SD_SNR.csv')))# LC signal-to-noise ratio for each sub-region

# load  fMRI and voxel-based-morphometry measures
FC <- as.data.frame(read.csv(here::here('Data','CC_fMRI_FC_ssrt.csv')))#psycho-physiological interaction measures (i.e. connectivity)
Bval <- as.data.frame(read.csv(here::here('Data','CC_fMRI_beta_ssrt.csv')))# Betas from fMRI
GM   <-as.data.frame(read.csv(here::here('Data','CC_VBM_ROI_8mmRad.csv')))# grey-matter



#load cognitive scores
tests2 <- as.data.frame(read.csv(here::here('Data','CC_cognitive_scores.csv')))

# Data and feature engineering ----------------------------------------------------------------------------------------------------------

#tidy up scores
cattel<- tests2 %>% rename('id'= 'SubCCIDc') %>% rename('scoreTot' = 'Cattell_TotalScore') %>% rename('STW'='STW_total') 


#strip prefixes and uniform column names  
FC$id %<>% gsub("[A-z]","",.)
Bval$id%<>% gsub("[A-z]","",.)
GM$ID%<>% gsub("[A-z]","",.)

names(GM)[1] <- 'id' #for consistency
cattel$id %<>% gsub("[A-z]","",.)
names(SNR)[[1]] = 'ccid'


#select only cases present in all datasets
commonID <- Reduce(dplyr::intersect,list(FC$id,Bval$id,GM$id,LC$Ccid))
LC   %<>% filter(Ccid %in% commonID)
FC   %<>% filter(id %in% commonID)
Bval %<>% filter(id %in% commonID)
GM   %<>% filter(id %in% commonID) 
ssrt   %<>% filter(subject %in% commonID)
cattel %<>% filter(id %in% commonID)




# labels indicate contrasts used in Tsvetanov et al.,2014, we are interested in SS_SSe (successful vs failed inhibition)
Bval <- Bval[Bval$contrast == 'SS_SSe',]
FC <- FC[FC$contrast == 'SS_SSe',]




#merge datasets into a single data frame
df<- FC %>%  cbind(as.data.frame(LC[,c(10:15)])) %>% 
    cbind(as.data.frame(Bval[,2:4])) %>% cbind(as.data.frame(GM[,-1])) %>%cbind(as.data.frame(cattel[,c('scoreTot','STW','acer')]))%>%
    cbind(as.data.frame(ssrt[,c('tf','tf_prob','SSRT','go_match')])) 



# Select sample of elders (aged >=50years) and remove participant with excessive head movements
df=df[df$age>=50,] %>% filter(SD_prob5>0)

# uniform CNR and SNR to the main data frame
CNR %<>% filter(ccid %in% df$id)
names(SNR)[[1]] = 'ccid'
SNR  %<>% filter(ccid %in% df$id) 

# center the variables since we are going to be dealing with interactions
df$agec        <-     c(scale(df$age,scale=F))
df$pSMA2SMA    <-     c(scale(df$pSMA2SMA,scale=F))
df$rIFG2pSMA    <-     c(scale(df$rIFG2pSMA,scale=F))
df$scoreTot    <-     c(scale(df$scoreTot))
df$STW          <-     c(scale(df$STW))
df$Discrep       <-      df$STW-df$scoreTot #calculate discrepancy scores


#Orthogonalise CR to deal with multi-collinearity between age and CR
df$CR<-  residuals(lm(SD_prob5 ~ agec, data = df))


# Descriptive statistics for demographics ------------------------------------------------------------------------------

# descriptive table
df %>% select(age, sex, acer) %>%
  droplevels() %>% tbl_summary() %>% bold_labels()

# define custom plotting theme
My_Theme = theme(
  plot.title = element_text(size = 19, hjust = 0.5, face = 'bold'),#hjust 0.5 is center
  axis.title.x = element_text(size = 19,vjust = -3),
  axis.text.x = element_text(size = 16,vjust = 0),#vjust adjust distance from tick to label
  axis.text.y = element_text(size = 16, margin = margin(t = 0, r = 5, b = 0, l = 0)),
  axis.title.y = element_text(size = 19,vjust = 3),
  axis.ticks = element_line(size=1),
  legend.text=element_text(size=19),
  legend.title = element_text(size=0),
  axis.ticks.length = unit(0.2, "cm"),
  plot.margin = margin(1, 1, 1,1, "cm"))

# plot boxplots for age and acer scores
ageplot <- ggplot(df,aes(x= 1,y = age)) + geom_boxplot(width=0.2,lwd  = 1,color = 'black')+geom_point(aes(x=1, y=age),  position = position_jitter(w = 0.08, h = 0),color="black", size=3,alpha=0.3)+
  
  theme_minimal()+My_Theme+theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  
  labs(title = 'Age',x="", y = "Years") + ylim(50,90) + xlim(0.6,1.4)


acerplot<-ggplot(df,aes(x= 1,y = acer)) + geom_boxplot(width=0.2,lwd  = 1,color = 'black')+geom_point(aes(x=1, y=acer),  position = position_jitter(w = 0.08, h = 0),color="black", size=3,alpha=0.3)+
  
  theme_minimal()+My_Theme+theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  
  labs(title = 'ACE-R', x="", y = "Score") + ylim(88,100) + xlim(0.6,1.4) 


ageplot+acerplot


# Do stats
#  t-test (frequentist & Bayesian) for differences between sexes in age and ACE-R scores
df %>% select(sex,age,acer) %>% gather(key = variable, value = value, -sex) %>%
  group_by(sex,variable) %>% summarise(value = list(value)) %>%
  spread(sex,value) %>%  group_by(variable) %>%
  mutate(p_value = t.test(unlist(F),unlist(M))$p.value,
         t_value = t.test(unlist(F),unlist(M))$statistic[[1]],
         BF      = round(as.data.frame(BayesFactor::ttestBF(unlist(F),unlist(M)))$bf,1))



# Hypothesis #1 -----------------------------------------------------------
# Decreases in locus coeruleus signal intensity are associated with reduced inhibitory control, indexed by increased SSRT.
# The relationship is moderated by age

fit1 <- lm(SSRT~CR*agec+sex, data = df)


#Check regression assumptions
gvlma::gvlma(x = fit1) 
sjPlot::plot_model(fit1, type = 'diag')
# Breusch-Pagan test for heteroscedasticity (p>0.05 insufficient evidence for heteroscedasticity)
lmtest::bptest(fit1)
#Check results
summary(fit1)
#calculate effect-size
effectsize::eta_squared(fit1,partial = FALSE)


#Plot results: The independent variable (CR) is continuous so its effect is measured through the slope of the regression line.
#Thus I plot the predicted values of those regression lines. Moderation says that the slope of the regression line is different
#at every value of the moderator.
sjPlot::plot_model(fit1, type = 'pred', terms = c('CR','agec'))


#Do Bayesian testing 
bf1<-BayesFactor::generalTestBF(SSRT~CR*agec+sex, data = df)
#Inclusion Bayes Factors test whether on average models with a given effect are more likely to have produced the observed
#data than the same models without the given effect.
bayestestR::bf_inclusion(bf1,match_models = TRUE)


# Hypothesis #2 -----------------------------------------------------------
# The relationship between locus coeruleus signal and SSRT varies with age, regardless of 
# age-related decline of grey matter volume in pSMA and rIFG regions.

fit2 <- lm(SSRT~CR*agec+sex+roi_rIFG+roi_rpSMA, data = df)

#Check regression assumptions
gvlma::gvlma(x = fit2) 
sjPlot::plot_model(fit2, type = 'diag')
lmtest::bptest(fit2)

#Check results
summary(fit2)
#calculate effect-size
effectsize::eta_squared(fit2,partial = FALSE)

bf2<-BayesFactor::generalTestBF(SSRT~CR*agec+sex+roi_rIFG+roi_rpSMA, data = df)
bayestestR::bf_inclusion(bf2,match_models = TRUE)

#use Robust regression to check whether results above are determined by outliers
fitr <-MASS::rlm(SSRT~CR*agec+sex+roi_rIFG+roi_rpSMA, data = df)
summary(fitr)

# robust F-test to estimate p-values
sfsmisc::f.robftest(fitr,'CR:agec')
sfsmisc::f.robftest(fitr,'sexM')
sfsmisc::f.robftest(fitr,'roi_rIFG')
sfsmisc::f.robftest(fitr,'roi_rpSMA')

# Hypothesis #3 -----------------------------------------------------------
#	The relationship between locus coeruleus signal and SSRT is strongest 
# in the caudal portion of the locus coeruleus.


# Convert CNR into wide format and bind to df
wCNR <- select(CNR,-med) %>% spread(roi,avg)
wCNR %<>% cbind(df) 
# scale values
# Note: scale() changes the class of the variables. It's safer simplifying them to numeric vectors using the dimension-stripping properties of c() 
wCNR$Caudal  %<>% scale() %>% c()
wCNR$Central %<>% scale() %>% c()
wCNR$Rostral %<>% scale() %>% c()
#orthogonalise against age to address multicollinearity
wCNR$Caudal  <-  residuals(lm(Caudal ~ agec, data = wCNR))
wCNR$Central <-  residuals(lm(Central ~ agec, data = wCNR))
wCNR$Rostral <-  residuals(lm(Rostral ~ agec, data = wCNR))

# prepare SNR data for later analyses
wCNR <- SNR %>% select(c('SNR_rostral','SNR_central','SNR_caudal')) %>% cbind(wCNR,.)


# fit linear model to test hypothesis 3
fit3 <- lm(SSRT~Caudal*Central*Rostral*agec+sex+roi_rIFG+roi_rpSMA,wCNR)

#Check regression assumptions
gvlma::gvlma(x = fit3) 

# Assumptions violated, check diagnostic to identify outliers
#plot(fit3)

# Remove outliers to meet lm assumptions
outl<-car::outlierTest(fit3)
outl <- as.numeric(names(outl$rstudent))
wCNR_clean<-wCNR[(is.na(match(rownames(wCNR),outl))),]

# fit model on cleaned data
fit3.1 <- lm(SSRT~Caudal*Central*Rostral*agec+sex+roi_rIFG+roi_rpSMA,wCNR_clean)

# Now all regression assumptions are satisfied
gvlma::gvlma(x = fit3.1) 

# However there is a huge multicollinearity 
sjPlot::plot_model(fit3.1, type = 'diag')
car::vif(fit3.1)


# Regression assumptions are met however there is strong collinearity, so here I deviate from pre-registration and perform
# Bayesian model comparison to  assess which sub-region drives the LC x age interaction:
# Compare models including an interaction term for each of the LC subregions against a model without interaction to establish
# the interaction term that best accounts for the data
Caudal <- BayesFactor::lmBF(formula = SSRT ~ (Caudal)*agec+sex+roi_rIFG+roi_rpSMA, data = wCNR_clean)
Central <- BayesFactor::lmBF(formula = SSRT ~ (Central)*agec+sex+roi_rIFG+roi_rpSMA, data = wCNR_clean)
Rostral <-  BayesFactor::lmBF(formula = SSRT ~ (Rostral)*agec+sex+roi_rIFG+roi_rpSMA, data = wCNR_clean)
No_Inter <- BayesFactor:::lmBF(formula = SSRT ~  agec+roi_rIFG+roi_rpSMA, data = wCNR_clean)

# Compare models
(bfs <- as.data.frame(c(Caudal,Central,Rostral)/No_Inter))

#Model selected: Rostral x age (Strong evidence 10.53 ±1.53%) against Caudal x age (anecdotal evidence 1.16 ±1.8%)  and middle x age (substantial evidence 7.29 ± 1.53%)
#Note: BF values can deviate marginally from the values above given the stochastic nature of the estimation procedure.



#A possible explanation for the anecdotal evidence in the Caudal portion, is that this part of the locus coeruleus has 
#high inter-subjective variability and partial volume effects. 
#Therefore, here I go beyond the pre-registered analysis and compare the signal-to-noise  ratio (SNR) amongst sub-regions
#of the locus coeruleus

# convert data frame into long format
  anovaROI_snr <- wCNR_clean %>% select(c('ccid', 'SNR_caudal','SNR_central','SNR_rostral')) %>%
    rename(Caudal = SNR_caudal, Central = SNR_central, Rostral = SNR_rostral) %>%
    pivot_longer(cols = Caudal:Rostral,names_to = "LC_ROI",values_to = "SNR") %>% mutate(ccid = as.factor(ccid),LC_ROI = as.factor(LC_ROI))
  
# detect outliers  
 outlier<-anovaROI_snr %>%
   group_by(LC_ROI) %>%
   rstatix::identify_outliers(SNR) %>% as_tibble()
 
 
# remove outliers
 anovaROI_snr <- anti_join(anovaROI_snr, outlier, by=c("LC_ROI", "ccid"))
 
# check  normality
 p<-ggplot(anovaROI_snr, aes(sample=SNR))+stat_qq() 
 p+qplot(sample = SNR,data = anovaROI_snr,shape  = LC_ROI)
 
 
 # Bayesian mixed effects ANOVA
 bayes_rm <- BayesFactor::anovaBF(SNR ~ LC_ROI + ccid,
                                  data = anovaROI_snr, whichRandom = "ccid") 
 bayes_rm
 
 # pairwise comparisons
# create an empty vector to store bfs
 bfpaired <- numeric(3)
 
 #calculate the BFs one by one
 pairs_c <- combn(c('Rostral','Central','Caudal'),2)
 for (i in 1:3) {
   
   bfpaired[i] <- anovaROI_snr %>% filter(LC_ROI %in% pairs_c[,i]) %>%
     droplevels() %>%
     as.data.frame() %>%
     BayesFactor::ttestBF(formula = SNR~LC_ROI, data = .) %>%
     BayesFactor::extractBF() %>%
     .$bf
   
 }
 
 


# Hypothesis #4 -----------------------------------------------------------
#The locus coeruleus contrast is specifically associated with SSRT, not Go reaction times.

fit4 <- lm(go_match~CR*agec+sex+roi_rIFG+roi_rpSMA, data = df)

#Check regression assumptions
gvlma::gvlma(x = fit4) 
sjPlot::plot_model(fit4, type = 'diag')
lmtest::bptest(fit4)

#Check results
summary(fit4)
#calculate effect-size
effectsize::eta_squared(fit4,partial = FALSE)

bf4<-BayesFactor::generalTestBF(go_match~CR*agec+sex+roi_rIFG+roi_rpSMA, data = df)
bayestestR::bf_inclusion(bf4,match_models = TRUE)


# Hypothesis #5 -----------------------------------------------------------
#Unsuccessful ageing is accompanied by reduced locus coeruleus signal and poor inhibitory control, indexed by increased SSRT.

fit5<-lm(SSRT~CR*Discrep*agec+sex+roi_rIFG+roi_rpSMA, data = df)

#Check regression assumptions
gvlma::gvlma(x = fit5) 
sjPlot::plot_model(fit5, type = 'diag')
lmtest::bptest(fit5)

#Check results
summary(fit5)
#calculate effect-size
effectsize::eta_squared(fit5,partial = FALSE)

bf5<-BayesFactor::generalTestBF(SSRT~CR*Discrep*agec+sex+roi_rIFG+roi_rpSMA, data = df)
bayestestR::bf_inclusion(bf5,match_models = TRUE)



# Hypothesis #6 -----------------------------------------------------------
#SSRT is related to functional connectivity between pSMA and rIFG, as a function of age .

fit6 <- lm(SSRT~rIFG2pSMA*agec+sex+roi_rIFG+roi_rpSMA, data = df)
#Check regression assumptions
gvlma::gvlma(x = fit6) 
sjPlot::plot_model(fit6, type = 'diag')
lmtest::bptest(fit6)

#Check results
summary(fit6)
#calculate effect-size
effectsize::eta_squared(fit6,partial = FALSE)

bf6<-BayesFactor::generalTestBF(SSRT~rIFG2pSMA*agec+sex+roi_rIFG+roi_rpSMA, data = df)
bayestestR::bf_inclusion(bf6,match_models = TRUE)

#use Robust regression to check whether results above are determined by outliers
fitr <-MASS::rlm(SSRT~rIFG2pSMA*agec+sex+roi_rIFG+roi_rpSMA, data = df)
summary(fitr)

# robust F-test to estimate p-values
sfsmisc::f.robftest(fitr,'rIFG2pSMA:agec')
sfsmisc::f.robftest(fitr,'sexM')
sfsmisc::f.robftest(fitr,'roi_rIFG')
sfsmisc::f.robftest(fitr,'roi_rpSMA')
sfsmisc::f.robftest(fitr,'rIFG2pSMA')

# Hypothesis #7 -----------------------------------------------------------
#locus coeruleus signal is related to SSRT, but this effect is moderated by functional connectivity between pSMA and rIFG regions.

library('processR')
library(lavaan)

#rename to FC for clarity
df %<>% rename('FC'= 'rIFG2pSMA')

#select structural equation model
pmacroModel(15)
#map variables into the model
labels = list(X='CR',M='FC',Y='SSRT',W='age')
pmacroModel(15,labels=labels)
moderator <- list(name = 'age',site = list(c('b','c')))
model = tripleEquation(X = 'CR',M='FC',Y='SSRT',moderator = moderator)
cat(model)

#fit SEM model using lavaan
semfit = sem(model=model,data=df,se='bootstrap', test = 'bootstrap', bootstrap = 10^3)

#the following syntax  indicates that the standard errors should be bias corrected (but not accelerated). 
#This approach will yield similar results to the PROCESS Macro in SPSS with bias-correct standard errors
parameterestimates(semfit, boot.ci.type = "bca.simple", standardized = TRUE)

# Moderated mediation results
summary(semfit, standardized = TRUE)
estimatesTable2(semfit,vanilla=TRUE)



#extract index moderated mediation
#estimates are obtained through stochastic bootstrapping hence can (marginally) change on each iteration
res <- parameterEstimates(semfit)
idx <- which(res$label == "index.mod.med")
modmed <- data.frame(estimate=res$est[idx],lowCI=res$ci.lower[idx],highCI=res$ci.upper[idx],
                     zvalue = res$z[idx], pvalue=res$pvalue[idx], propmediated = (res$est[36]+res$est[37]))




conditionalEffectPlot(semfit,data=df,mod="age")

