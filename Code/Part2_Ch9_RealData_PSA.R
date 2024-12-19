###############################################################################
## 제2부 제9장 실제 데이터 분석실습 및 매칭기법 비교
###############################################################################
library("Hmisc")              #가중평균/가중분산 
library("tidyverse")          #데이터관리 및 변수사전처리 
summarize = dplyr::summarize  #dplyr의 summarize 함수를 이용
library("Zelig")              #비모수접근 95% CI 계산 
library("MatchIt")            #매칭기법 적용 
library("cobalt")             #처치집단과 통제집단 균형성 점검 
library("treatSens")          #민감도 테스트 실시
library("sensitivitymw")      #민감도 테스트 실시
library("sensitivityfull")    #민감도 테스트 실시 

###############################################################################
## 1. 데이터 소개 
###############################################################################
setwd("D:/data")
mydata=read_csv("observational_study_survey.csv")
# 세대변수는 가변수들 집단으로 리코딩 
mydata=mydata %>% 
  mutate(
    gen20=ifelse(gen=='20s',1,0),
    gen30=ifelse(gen=='30s',1,0),
    gen40=ifelse(gen=='40s',1,0),
    gen50=ifelse(gen=='50s',1,0)
  ) 
# 두번째 처치효과인 보수단체 탄핵반대 집회참여여부는 고려하지 않음 
mydata=mydata %>% 
  filter(rally_con==0) %>%
  select(-rally_con)
# ATC 계산위해 촛불집회 참석여부 더미변수 역코딩 
mydata$Rtreat=ifelse(mydata$rally_pro==1,0,1)
# ATT 계산위해 변수이름 재조정 
mydata = mydata %>% 
  rename(y=vote_will,
         treat=rally_pro)
# 기술통계분석: 변수별 기술통계치 비교 
mean_SD_range = function(myvariable){
  myM=format(round(mean(myvariable),2),2)
  mySD=format(round(sd(myvariable),2),2)
  myMn=format(round(min(myvariable),2),2)
  myMx=format(round(max(myvariable),2),2)
  str_c(myM,"\n(",mySD,")","\n[",myMn,", ",myMx,"]")
}
# 기술통계분석: 처치집단과 통제집단 비교 
mean_SD = function(myvariable){
  myM=format(round(mean(myvariable),2),2)
  mySD=format(round(sd(myvariable),2),2)
  str_c(myM,"\n(",mySD,")")
}
# 집단별 종속변수 및 공변량 비교 
mydata %>% group_by(treat) %>% 
  summarize_at(
    vars(y,female,edu:vote_past,gen20:gen50),
    mean_SD
  ) %>% 
  pivot_longer(cols=2:14, names_to="vars") %>% 
  mutate(treat=ifelse(treat==0,"Control","Treated")) %>% 
  pivot_wider(names_from="treat", values_from=value) %>% 
  write.csv("DesStat_Table01.csv",row.names=FALSE)

# 통상적 OLS 회귀분석
# 공식이 길기 때문에 공변량들만 따로 저장 
cvrt="female+gen20+gen30+gen40+gen50+edu+hhinc+libcon+int_eff+
      park_eva_a+good_eco+vote_past"
y_T_cov=as.formula(str_c("y~treat+",cvrt))
OLS_estimands=lm(y_T_cov,mydata) %>% 
  confint("treat") %>% as_tibble()
names(OLS_estimands)=c("LL95","UL95")
OLS_estimands
###############################################################################
## 2. 성향점수분석 기법 적용 
###############################################################################
att_formula=as.formula(str_c("treat~",cvrt))  
atc_formula=as.formula(str_c("Rtreat~",cvrt)) 
## 1) 성향점수가중 기법 
# 성향점수 계산 
mydata$pscore = glm(att_formula,
                    data=mydata,
                    family=binomial(link='logit')) %>% 
  fitted()
# IPTW 계산(ATE, ATT, ATC 추정) 
mydata = mydata %>% 
  mutate(
    Wate = ifelse(treat==1,1/pscore,1/(1-pscore)),
    Watt = ifelse(treat==1,1,1/(1-pscore)),
    Watc = ifelse(treat==1,1/pscore,1)
  )

## 2) 성향점수기반 그리디 매칭(Greedy matching using propensity score)
# ATT
set.seed(1234)  #정확하게 동일한 결과를 원한다면 
greedy_att=matchit(formula=att_formula,
                   data=mydata,  
                   distance="linear.logit", #선형로짓 형태 성향점수 
                   method="nearest",  # 그리디 매칭 알고리즘
                   caliper=0.15,  #성향점수의 0.15표준편차 허용범위 
                   discard='none',  #공통지지영역에서 벗어난 사례도 포함 
                   ratio=2,  #처치집단 사례당 통제집단 사례는 2개 매칭 
                   replace=FALSE) #동일사례 반복표집 매칭을 허용하지 않음
greedy_att$nn
# ATC
set.seed(4321)  #정확하게 동일한 결과를 원한다면 
greedy_atc=matchit(atc_formula, #통제집단을 '처치집단'으로 가정한 후 실행 
                   data=mydata,  
                   distance="linear.logit", #선형로짓 형태 성향점수 
                   method="nearest",  # 그리디 매칭 알고리즘
                   caliper=0.15,  #성향점수의 0.15표준편차 허용범위 
                   discard='none',  #공통지지영역에서 벗어난 사례도 포함 
                   ratio=2,  #처치집단 사례당 통제집단 사례는 2개 매칭 
                   replace=TRUE) #동일사례 반복표집 매칭을 허용
greedy_atc$nn

# 3) 성향점수기반 최적 매칭(Optimal matching using propensity score)
# ATT 
set.seed(1234)  #정확하게 동일한 결과를 원한다면 
optimal_att=matchit(formula=att_formula,
                    data=mydata,  
                    distance="linear.logit", #선형로짓 형태 성향점수 
                    method="optimal",  # 최적 매칭 알고리즘
                    caliper=0.15,  #성향점수의 0.15표준편차 허용범위 
                    discard='none', #공통지지영역에서 벗어난 사례 보존 
                    ratio=2,  #처치집단 사례당 통제집단 사례는 2개 매칭 
                    replace=FALSE) #동일사례 반복표집 매칭을 허용 
optimal_att$nn
# ATC, 추정불가

# 4) 성향점수기반 전체 매칭(Full matching using propensity score)
# ATT
set.seed(1234)  #정확하게 동일한 결과를 원한다면 
full_att=matchit(formula=att_formula,
                 data=mydata,
                 distance="linear.logit", #선형로짓 형태 성향점수 
                 method="full",  # 전체 매칭 알고리즘
                 discard="none")  #공통지지영역에서 벗어난 사례 보존 
full_att$nn
# ATC
set.seed(4321)  #정확하게 동일한 결과를 원한다면 
full_atc=matchit(formula=atc_formula, #통제집단을 '처치집단'으로 가정한 후 실행 
                 data=mydata,
                 distance="linear.logit", #선형로짓 형태 성향점수 
                 method="full",  # 전체 매칭 알고리즘
                 discard="none")  #공통지지영역에서 벗어난 사례 보존 
full_atc$nn

# 아래의 과정은 Google Cloud Platform으로 진행하였음
# 추정결과는 genetic_att_ch9.RData; genetic_atc_ch9.RData로 저장되어 있음. 
# 5) 성향점수기반 유전 매칭(Genetic matching using propensity score)
# install.packages("rgenoud");library("rgenoud") # 패키지 오류가 발생한다면
# ATT
set.seed(1234)  #정확하게 동일한 결과를 원한다면 
genetic_att=matchit(formula=att_formula,
                    data=mydata,
                    method="genetic",
                    distance="linear.logit",
                    pop.size=1000, #디폴트는 100인데, 안정적 추정을 위해 늘림
                    discard='none',
                    ratio=2,
                    distance.tolerance=1e-05, # 거리차이 기준값 디폴트 
                    ties=TRUE)  #복수의 통제집단 사례들 매칭(디폴트)
genetic_att$nn
# ATC
set.seed(4321)  #정확하게 동일한 결과를 원한다면 
genetic_atc=matchit(formula=atc_formula,
                    data=mydata,
                    method="genetic",
                    distance="linear.logit",
                    pop.size=1000, #디폴트는 100인데, 안정적 추정을 위해 늘림
                    discard='none',
                    ratio=2,
                    distance.tolerance=1e-05, # 거리차이 기준값 디폴트 
                    ties=TRUE)  #복수의 통제집단 사례들 매칭(디폴트)
genetic_atc$nn
# genetic_att=readRDS("genetic_att_ch9.RData");genetic_atc=readRDS("genetic_atc_ch9.RData")


# 6) 마할라노비스 거리점수기반 그리디 매칭(Greedy matching using Mahalanobis distance)
# ATT
set.seed(1234)  #정확하게 동일한 결과를 원한다면 
mahala_att=matchit(formula=att_formula,
                   data=mydata, 
                   distance="mahalanobis", #마할라노비스 거리점수 
                   method="nearest",  # 그리디 매칭 알고리즘
                   caliper=0.15,  #성향점수의 0.15표준편차 허용범위 
                   discard='none',  #공통지지영역에서 벗어난 사례도 포함 
                   ratio=2,  #처치집단 사례당 통제집단 사례는 2개 매칭 
                   replace=FALSE) #동일사례 반복표집 매칭을 허용하지 않음
mahala_att$nn
# ATC 
set.seed(4321)  #정확하게 동일한 결과를 원한다면 
mahala_atc=matchit(formula=atc_formula, #통제집단을 '처치집단'으로 가정한 후 실행 
                   data=mydata,  
                   distance="mahalanobis", #마할라노비스 거리점수 
                   method="nearest",  # 그리디 매칭 알고리즘
                   caliper=0.15,  #성향점수의 0.15표준편차 허용범위 
                   discard='none',  #공통지지영역에서 벗어난 사례도 포함 
                   ratio=2,  #처치집단 사례당 통제집단 사례는 2개 매칭 
                   replace=TRUE) #동일사례 반복표집 매칭을 허용
mahala_atc$nn

# 7) 성향점수층화 기법(집단수는 디폴트인 6을 적용함) 
# ATT 추정
mysubclassN=6
set.seed(1234)
class_att=matchit(formula=att_formula,
                  data=mydata,
                  distance="linear.logit",
                  method="subclass",   # 층화 기법 
                  discard='none', 
                  subclass=mysubclassN)
class_att$nn
# ATC 추정
set.seed(4321)
class_atc=matchit(formula=atc_formula,
                  data=mydata,
                  distance="linear.logit",
                  method="subclass",   # 층화 기법 
                  discard='none', 
                  subclass=mysubclassN)
class_atc$nn


# 8) 준정확매칭(Coarsened exact matching)
# 만약 자동화 뭉뚱그리기를 택하는 경우
matchit(formula=att_formula,data=mydata,
        method="cem",discard='none',L1.break='scott')$nn
matchit(formula=atc_formula,data=mydata,
        method="cem",discard='none',L1.break='scott')$nn
# 이용자선택(user-choice) 뭉뚱그리기 
# 분리지점을 수동으로 선택 
# 세대는 20-40/50-60으로(별도의 더미변수 생성)
mydata=mydata %>% 
  mutate(
    gen5060=ifelse(gen=='50s'|gen=='60s',1,0)
  ) 
gen20_cp = c(-0.5,1.5);gen30_cp = c(-0.5,1.5);gen40_cp = c(-0.5,1.5)
gen50_cp = c(-0.5,1.5) # gen5060 변수만 의미를 갖도록
edu_cp = c(0,3.5,4.5,6)
hhinc_cp = c(0,2.5,5.5,16)
libcon_cp = c(-6,-0.5,0.5,6)
int_eff_cp = c(0,2.5,3.5,5.5)
park_eva_a_cp = c(0,3.5,5)
# CEM의 경우 공식재지정
att_formula_cem = update(att_formula, .~.+gen5060)
atc_formula_cem = update(atc_formula, .~.+gen5060)
# ATT 
set.seed(1234)
cem_att=matchit(formula=att_formula,
                data=mydata,
                method="cem",
                discard='none',
                cutpoints=list(gen20=gen20_cp,
                               gen30=gen30_cp,
                               gen40=gen40_cp,
                               gen50=gen50_cp,
                               edu=edu_cp,
                               hhinc=hhinc_cp,
                               libcon=libcon_cp,
                               int_eff=int_eff_cp,
                               park_eva_a=park_eva_a_cp))
cem_att$nn
# ATC 
set.seed(4321)
cem_atc=matchit(formula=atc_formula,
                data=mydata,
                method="cem",
                discard='none',
                cutpoints=list(gen20=gen20_cp,
                               gen30=gen30_cp,
                               gen40=gen40_cp,
                               gen50=gen50_cp,
                               edu=edu_cp,
                               hhinc=hhinc_cp,
                               libcon=libcon_cp,
                               int_eff=int_eff_cp,
                               park_eva_a=park_eva_a_cp))
cem_atc$nn

###############################################################################
## 3. 공변량 균형성 점검
###############################################################################
# 1) 성향점수가중(PSW) 분석
# 앞서 PSW때 사용한 이용자정의 함수 
balance_check_PSW=readRDS("balance_check_PSW.RData")
# 연속형 변수 형태 공변량들만 선별한 데이터 
mydata_cov = mydata %>% 
  select(pscore,edu:good_eco) %>% 
  as.data.frame()
# 모든 공변량들에 대해 가중평균/가중분산 계산 
result_covbal = list()
for (i in 1:dim(mydata_cov)[2]){
  result_covbal = bind_rows(result_covbal,
                            balance_check_PSW(
                              mydata$treat,mydata_cov[,i],mydata$Wate))
}
result_covbal$cov_name=names(mydata_cov)
result_covbal
# 역치를 넘어서는 공변량은? 
result_covbal %>% 
  filter(abs(A_wgt_Mdiff) >= 0.25|  #절댓값 평균차이가 0.25를 넘는 경우 
           A_wgt_Vratio > 2| A_wgt_Vratio < 0.5)  #분산비 0.5미만 혹은 2.0초과 
# 더미변수들의 경우 표준화를 시키는 것이 타당하지 않을 수도 있음
MeanDifferenceRaw=function(myV){
  wgted_M1 = Hmisc::wtd.mean(x=myV[mydata$treat==1],weights=mydata$Wate[mydata$treat==1])
  wgted_M0 = Hmisc::wtd.mean(x=myV[mydata$treat==0],weights=mydata$Wate[mydata$treat==0])
  wgt_MD=abs(wgted_M1-wgted_M0)
}
mydata %>% 
  select(female,gen20:gen50,vote_past) %>%
  summarize_all(
    MeanDifferenceRaw
  ) %>% 
  pivot_longer(cols=female:vote_past,names_to="Var_Raw")

# 2) 성향점수기반 그리디 매칭 
# 러브플롯 
# ATT, 평균차이 
love_D_greedy_att=love.plot(greedy_att,
                            s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                            stat="mean.diffs", # 두 집단간 공변량 평균차이 
                            drop.distance=FALSE, #성향점수도 포함하여 제시 
                            threshold=0.25, # 평균차이 역치 
                            stars="raw",  #표준화점수가 아닌 원점수를 사용한 경우 
                            sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                            themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00))+ggtitle("ATT, covariate balance")
# ATT, 분산비 
love_VR_greedy_att=love.plot(greedy_att,
                             s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                             stat="variance.ratios", # 두 집단간 분산비 
                             drop.distance=FALSE, #성향점수도 포함하여 제시 
                             threshold=2, # 분산비 역치 
                             sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                             themes=theme_bw())+
  coord_cartesian(xlim=c(0.33,3))+ggtitle("ATT, covariate balance")
# ATC, 평균차이 
love_D_greedy_atc=love.plot(greedy_atc,
                            s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                            stat="mean.diffs", # 두 집단간 공변량 평균차이 
                            drop.distance=FALSE, #성향점수도 포함하여 제시 
                            threshold=0.25, # 평균차이 역치 
                            stars="raw",  #표준화점수가 아닌 원점수를 사용한 경우 
                            sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                            themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00))+ggtitle("ATC, covariate balance")
# ATC, 분산비 
love_VR_greedy_atc=love.plot(greedy_atc,
                             s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                             stat="variance.ratios", # 두 집단간 분산비 
                             drop.distance=FALSE, #성향점수도 포함하여 제시 
                             threshold=2, # 분산비 역치 
                             sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                             themes=theme_bw())+
  coord_cartesian(xlim=c(0.33,3))+ggtitle("ATC, covariate balance")
png("Part2_ch9_love_MD_VR_greedy.png",res=600,height=15,width=25,units="cm")
gridExtra::grid.arrange(love_D_greedy_att,love_VR_greedy_att,
                        love_D_greedy_atc,love_VR_greedy_atc,
                        nrow=2)
dev.off()

# 3) 성향점수기반 최적 매칭 
# 러브플롯 
# ATT, 평균차이 
love_D_optimal_att=love.plot(optimal_att,
                            s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                            stat="mean.diffs", # 두 집단간 공변량 평균차이 
                            drop.distance=FALSE, #성향점수도 포함하여 제시 
                            threshold=0.25, # 평균차이 역치 
                            stars="raw",  #표준화점수가 아닌 원점수를 사용한 경우 
                            sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                            themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00))+ggtitle("ATT, covariate balance")
# ATT, 분산비 
love_VR_optimal_att=love.plot(optimal_att,
                             s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                             stat="variance.ratios", # 두 집단간 분산비 
                             drop.distance=FALSE, #성향점수도 포함하여 제시 
                             threshold=2, # 분산비 역치 
                             sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                             themes=theme_bw())+
  coord_cartesian(xlim=c(0.33,3))+ggtitle("ATT, covariate balance")
png("Part2_ch9_love_MD_VR_optimal.png",res=600,height=8,width=25,units="cm")
gridExtra::grid.arrange(love_D_optimal_att,love_VR_optimal_att,nrow=1)
dev.off()

# 4) 성향점수기반 전체 매칭 
# 러브플롯 
# ATT, 평균차이 
love_D_full_att=love.plot(full_att,
                          s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                          stat="mean.diffs", # 두 집단간 공변량 평균차이 
                          drop.distance=FALSE, #성향점수도 포함하여 제시 
                          threshold=0.25, # 평균차이 역치 
                          stars="raw",  #표준화점수가 아닌 원점수를 사용한 경우 
                          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                          themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00))+ggtitle("ATT, covariate balance")
# ATT, 분산비 
love_VR_full_att=love.plot(full_att,
                           s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                           stat="variance.ratios", # 두 집단간 분산비 
                           drop.distance=FALSE, #성향점수도 포함하여 제시 
                           threshold=2, # 분산비 역치 
                           sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                           themes=theme_bw())+
  coord_cartesian(xlim=c(0.33,3))+ggtitle("ATT, covariate balance")
# ATC, 평균차이 
love_D_full_atc=love.plot(full_atc,
                          s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                          stat="mean.diffs", # 두 집단간 공변량 평균차이 
                          drop.distance=FALSE, #성향점수도 포함하여 제시 
                          threshold=0.25, # 평균차이 역치 
                          stars="raw",  #표준화점수가 아닌 원점수를 사용한 경우 
                          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                          themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00))+ggtitle("ATC, covariate balance")
# ATC, 분산비 
love_VR_full_atc=love.plot(full_atc,
                           s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                           stat="variance.ratios", # 두 집단간 분산비 
                           drop.distance=FALSE, #성향점수도 포함하여 제시 
                           threshold=2, # 분산비 역치 
                           sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                           themes=theme_bw())+
  coord_cartesian(xlim=c(0.33,3))+ggtitle("ATC, covariate balance")
png("Part2_ch9_love_MD_VR_full.png",res=600,height=15,width=25,units="cm")
gridExtra::grid.arrange(love_D_full_att,love_VR_full_att,
                        love_D_full_atc,love_VR_full_atc,
                        nrow=2)
dev.off()


# 5) 성향점수기반 유전 매칭 
# 러브플롯 
# ATT, 평균차이 
love_D_genetic_att=love.plot(genetic_att,
                          s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                          stat="mean.diffs", # 두 집단간 공변량 평균차이 
                          drop.distance=FALSE, #성향점수도 포함하여 제시 
                          threshold=0.25, # 평균차이 역치 
                          stars="raw",  #표준화점수가 아닌 원점수를 사용한 경우 
                          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                          themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00))+ggtitle("ATT, covariate balance")
# ATT, 분산비 
love_VR_genetic_att=love.plot(genetic_att,
                           s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                           stat="variance.ratios", # 두 집단간 분산비 
                           drop.distance=FALSE, #성향점수도 포함하여 제시 
                           threshold=2, # 분산비 역치 
                           sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                           themes=theme_bw())+
  coord_cartesian(xlim=c(0.33,3))+ggtitle("ATT, covariate balance")
# ATC, 평균차이 
love_D_genetic_atc=love.plot(genetic_atc,
                          s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                          stat="mean.diffs", # 두 집단간 공변량 평균차이 
                          drop.distance=FALSE, #성향점수도 포함하여 제시 
                          threshold=0.25, # 평균차이 역치 
                          stars="raw",  #표준화점수가 아닌 원점수를 사용한 경우 
                          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                          themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00))+ggtitle("ATC, covariate balance")
# ATC, 분산비 
love_VR_genetic_atc=love.plot(genetic_atc,
                           s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                           stat="variance.ratios", # 두 집단간 분산비 
                           drop.distance=FALSE, #성향점수도 포함하여 제시 
                           threshold=2, # 분산비 역치 
                           sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                           themes=theme_bw())+
  coord_cartesian(xlim=c(0.33,3))+ggtitle("ATC, covariate balance")
png("Part2_ch9_love_MD_VR_genetic.png",res=600,height=15,width=25,units="cm")
gridExtra::grid.arrange(love_D_genetic_att,love_VR_genetic_att,
                        love_D_genetic_atc,love_VR_genetic_atc,
                        nrow=2)
dev.off()


# 6) 마할라노비스 거리점수기반 그리디 매칭 
# 러브플롯 
# ATT, 평균차이 
love_D_mahala_att=love.plot(mahala_att,
                             s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                             stat="mean.diffs", # 두 집단간 공변량 평균차이 
                             drop.distance=TRUE, #성향점수 없음 
                             threshold=0.25, # 평균차이 역치 
                             stars="raw",  #표준화점수가 아닌 원점수를 사용한 경우 
                             sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                             themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00))+ggtitle("ATT, covariate balance")
# ATT, 분산비 
love_VR_mahala_att=love.plot(mahala_att,
                              s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                              stat="variance.ratios", # 두 집단간 분산비 
                              drop.distance=TRUE, #성향점수 없음 
                              threshold=2, # 분산비 역치 
                              sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                              themes=theme_bw())+
  coord_cartesian(xlim=c(0.33,3))+ggtitle("ATT, covariate balance")
# ATC, 평균차이 
love_D_mahala_atc=love.plot(mahala_atc,
                             s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                             stat="mean.diffs", # 두 집단간 공변량 평균차이 
                             drop.distance=TRUE, #성향점수 없음 
                             threshold=0.25, # 평균차이 역치 
                             stars="raw",  #표준화점수가 아닌 원점수를 사용한 경우 
                             sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                             themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00))+ggtitle("ATC, covariate balance")
# ATC, 분산비 
love_VR_mahala_atc=love.plot(mahala_atc,
                              s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                              stat="variance.ratios", # 두 집단간 분산비 
                              drop.distance=TRUE, #성향점수 없음 
                              threshold=2, # 분산비 역치 
                              sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                              themes=theme_bw())+
  coord_cartesian(xlim=c(0.33,3))+ggtitle("ATC, covariate balance")
png("Part2_ch9_love_MD_VR_mahala.png",res=600,height=15,width=25,units="cm")
gridExtra::grid.arrange(love_D_mahala_att,love_VR_mahala_att,
                        love_D_mahala_atc,love_VR_mahala_atc,
                        nrow=2)
dev.off()

# 7) 성향점수층화 기법
# ATT-평균차이/분산비 
bal.tab(class_att,
        continuous="std",s.d.denom="pooled",
        m.threshold=0.25,v.threshold=2)
# ATC-평균차이 /분산비 
bal.tab(class_atc,
        continuous="std",s.d.denom="pooled",
        m.threshold=0.25,v.threshold=2)

# 8) 준정확 매칭 
# 러브플롯 
# ATT, 평균차이 
love_D_cem_att=love.plot(cem_att,
                            s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                            stat="mean.diffs", # 두 집단간 공변량 평균차이 
                            drop.distance=TRUE, #성향점수 없음 
                            threshold=0.25, # 평균차이 역치 
                            stars="raw",  #표준화점수가 아닌 원점수를 사용한 경우 
                            sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                            themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00))+ggtitle("ATT, covariate balance")
# ATT, 분산비 
love_VR_cem_atc=love.plot(cem_att,
                             s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                             stat="variance.ratios", # 두 집단간 분산비 
                             drop.distance=TRUE, #성향점수 없음 
                             threshold=2, # 분산비 역치 
                             sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                             themes=theme_bw())+
  coord_cartesian(xlim=c(0.33,3))+ggtitle("ATT, covariate balance")
# ATC, 평균차이 
love_D_cem_atc=love.plot(cem_atc,
                            s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                            stat="mean.diffs", # 두 집단간 공변량 평균차이 
                            drop.distance=TRUE, #성향점수 없음 
                            threshold=0.25, # 평균차이 역치 
                            stars="raw",  #표준화점수가 아닌 원점수를 사용한 경우 
                            sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                            themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00))+ggtitle("ATC, covariate balance")
# ATC, 분산비 
love_VR_cem_atc=love.plot(cem_atc,
                             s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                             stat="variance.ratios", # 두 집단간 분산비 
                             drop.distance=TRUE, #성향점수 없음 
                             threshold=2, # 분산비 역치 
                             sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                             themes=theme_bw())+
  coord_cartesian(xlim=c(0.33,3))+ggtitle("ATC, covariate balance")
png("Part2_ch9_love_MD_VR_cem.png",res=600,height=15,width=25,units="cm")
gridExtra::grid.arrange(love_D_cem_att,love_VR_cem_att,
                        love_D_cem_atc,love_VR_cem_atc,
                        nrow=2)
dev.off()


###############################################################################
## 3. 처치효과 추정  
###############################################################################
# 1) 성향점수가중(PSW) 분석
# ATT 추정 
set.seed(1234)  # 동일한 결과를 얻고자 한다면.... 
z_model=zelig(y_T_cov,data=mydata,model='ls', 
              weights="Watt",cite=FALSE)
x_1=setx(z_model,treat=1,data=mydata) 
x_0=setx(z_model,treat=0,data=mydata) 
s_1=sim(z_model,x_1,num=10000)
s_0=sim(z_model,x_0,num=10000)
EST_att=get_qi(s_1,"ev") - get_qi(s_0,"ev")
summary_est_att=tibble(
  LL95=quantile(EST_att,p=c(0.025)),
  PEst=quantile(EST_att,p=c(0.500)),
  UL95=quantile(EST_att,p=c(0.975)),
  estimand="ATT",model="Propensity score weighting"
)
# ATC 추정 
set.seed(1234)  # 동일한 결과를 얻고자 한다면.... 
z_model=zelig(y_T_cov,data=mydata,model='ls',
              weights="Watc",cite=FALSE)
x_1=setx(z_model,treat=1,data=mydata) 
x_0=setx(z_model,treat=0,data=mydata) 
s_1=sim(z_model,x_1,num=10000)
s_0=sim(z_model,x_0,num=10000)
EST_atc=get_qi(s_1,"ev") - get_qi(s_0,"ev")
summary_est_atc=tibble(
  LL95=quantile(EST_atc,p=c(0.025)),
  PEst=quantile(EST_atc,p=c(0.500)),
  UL95=quantile(EST_atc,p=c(0.975)),
  estimand="ATC",model="Propensity score weighting"
)
# ATE 추정 
set.seed(1234)  # 동일한 결과를 얻고자 한다면.... 
z_model=zelig(y_T_cov,data=mydata,model='ls',
              weights="Wate",cite=FALSE)
x_1=setx(z_model,treat=1,data=mydata) 
x_0=setx(z_model,treat=0,data=mydata) 
s_1=sim(z_model,x_1,num=10000)
s_0=sim(z_model,x_0,num=10000)
EST_ate=get_qi(s_1,"ev") - get_qi(s_0,"ev")
summary_est_ate=tibble(
  LL95=quantile(EST_ate,p=c(0.025)),
  PEst=quantile(EST_ate,p=c(0.500)),
  UL95=quantile(EST_ate,p=c(0.975)),
  estimand="ATE",model="Propensity score weighting"
)
# 3가지 효과추정치(estimands) 통합 
PSW_estimands = bind_rows(
  summary_est_att,summary_est_atc,summary_est_ate
)
PSW_estimands  

## 2) 성향점수기반 그리디 매칭 기법 
# 이용자정의 함수 불러오기 (제6장에서 이미 사용된 함수들임)
SUMMARY_EST_ATT = readRDS("SUMMARY_EST_ATT.RData")
SUMMARY_EST_ATC = readRDS("SUMMARY_EST_ATC.RData")
# 처치효과 추정 공식 저장
y_T_cov=as.formula(str_c("y~treat+",cvrt))
y_R_cov=as.formula(str_c("y~Rtreat+",cvrt))
# ATT, ATC 추정 
MD_greedy_att=match.data(greedy_att)
MD_greedy_atc=match.data(greedy_atc)
set.seed(1234)
greedy_ATT = SUMMARY_EST_ATT(myformula=y_T_cov,
                             matched_data=MD_greedy_att,
                             n_sim=10000,
                             model_name="Greedy matching using propensity score")
set.seed(4321)
greedy_ATC = SUMMARY_EST_ATC(myformula=y_R_cov,
                             matched_data=MD_greedy_atc,
                             n_sim=10000,
                             model_name="Greedy matching using propensity score")
# ATE 추정을 위한 이용자정의 함수 
SUMMARY_EST_ATE = function(Out_ATT,Out_ATC,variable_treat,my_naming){
  mypi=prop.table(table(variable_treat))[2]
  EST1000=mypi*as.vector(Out_ATT[[1]]) + (1-mypi)*(as.vector(Out_ATC[[1]]))
  summary_est=tibble(
    LL95=quantile(EST1000,p=c(0.025)),
    PEst=quantile(EST1000,p=c(0.500)),
    UL95=quantile(EST1000,p=c(0.975)),
    estimand="ATE",model=my_naming
  )
  summary_est
}
greedy_ATE=SUMMARY_EST_ATE(greedy_ATT,greedy_ATC,mydata$treat,
                          "Greedy matching using propensity score")
# 3가지 효과추정치(estimands) 통합 
greedy_estimands = bind_rows(
  greedy_ATT[[2]],greedy_ATC[[2]],greedy_ATE
)
greedy_estimands  


## 3) 성향점수기반 최적 매칭 기법 
# ATT 추정 
MD_optimal_att=match.data(optimal_att)
set.seed(1234)
optimal_ATT = SUMMARY_EST_ATT(myformula=y_T_cov,
                             matched_data=MD_optimal_att,
                             n_sim=10000,
                             model_name="Optimal matching using propensity score")
set.seed(4321)
# 3가지 효과추정치(estimands) 통합 
optimal_estimands = optimal_ATT[[2]]
optimal_estimands


## 4) 성향점수기반 전체 매칭 기법 
# ATT, ATC 추정 
MD_full_att=match.data(full_att)
MD_full_atc=match.data(full_atc)
set.seed(1234)
full_ATT = SUMMARY_EST_ATT(myformula=y_T_cov,
                           matched_data=MD_full_att,
                           n_sim=10000,
                           model_name="Full matching using propensity score")
set.seed(4321)
full_ATC = SUMMARY_EST_ATC(myformula=y_R_cov,
                           matched_data=MD_full_atc,
                           n_sim=10000,
                           model_name="Full matching using propensity score")
full_ATE=SUMMARY_EST_ATE(full_ATT,full_ATC,mydata$treat,
                         "Full matching using propensity score")
# 3가지 효과추정치(estimands) 통합 
full_estimands = bind_rows(
  full_ATT[[2]],full_ATC[[2]],full_ATE
)
full_estimands  


## 5) 성향점수기반 유전 매칭 기법 
# ATT, ATC 추정 
MD_genetic_att=match.data(genetic_att)
MD_genetic_atc=match.data(genetic_atc)
set.seed(1234)
genetic_ATT = SUMMARY_EST_ATT(myformula=y_T_cov,
                           matched_data=MD_genetic_att,
                           n_sim=10000,
                           model_name="Genetic matching using propensity score")
set.seed(4321)
genetic_ATC = SUMMARY_EST_ATC(myformula=y_R_cov,
                           matched_data=MD_genetic_atc,
                           n_sim=10000,
                           model_name="Genetic matching using propensity score")
genetic_ATE=SUMMARY_EST_ATE(genetic_ATT,genetic_ATC,mydata$treat,
                         "Genetic matching using propensity score")
# 3가지 효과추정치(estimands) 통합 
genetic_estimands = bind_rows(
  genetic_ATT[[2]],genetic_ATC[[2]],genetic_ATE
)
genetic_estimands  


## 6) 마할라노비스 거리점수기반 그리디 매칭 기법 
# ATT, ATC 추정 
MD_mahala_att=match.data(mahala_att)
MD_mahala_atc=match.data(mahala_atc)
set.seed(1234)
mahala_ATT = SUMMARY_EST_ATT(myformula=y_T_cov,
                             matched_data=MD_mahala_att,
                             n_sim=10000,
                             model_name="Greedy matching using Mahala. distance")
set.seed(4321)
mahala_ATC = SUMMARY_EST_ATC(myformula=y_R_cov,
                             matched_data=MD_mahala_atc,
                             n_sim=10000,
                             model_name="Greedy matching using Mahala. distance")
mahala_ATE=SUMMARY_EST_ATE(mahala_ATT,mahala_ATC,mydata$treat,
                           "Greedy matching using Mahala. distance")
# 3가지 효과추정치(estimands) 통합 
mahala_estimands = bind_rows(
  mahala_ATT[[2]],mahala_ATC[[2]],mahala_ATE
)
mahala_estimands  


## 7) 성향점수층화 기법 
# ATT
MD_class_att=match.data(class_att)
MD_class_atc=match.data(class_atc)
set.seed(1234)
total_wgt=sum(MD_class_att$weights)  # 총 데이터 사례들 
all_sub_sim=list()           # 집단별 시뮬레이션 결과 저장 
for(i in 1:mysubclassN){
  temp_md_sub=MD_class_att %>% filter(subclass==i)  #i번째 집단 
  # 효과추정치 추정모형 
  z_sub_att=zelig(y_T_cov,data=temp_md_sub,
                  model='ls',weights="weights",cite=F)  
  # 독립변수 조건 
  x_sub_att0=setx(z_sub_att,treat=0,data=temp_md_sub)
  x_sub_att1=setx(z_sub_att,treat=1,data=temp_md_sub)
  # 시뮬레이션 
  s_sub_att0=sim(z_sub_att,x_sub_att0,num=10000)
  s_sub_att1=sim(z_sub_att,x_sub_att1,num=10000)
  temp_sim=get_qi(s_sub_att1,"ev")-get_qi(s_sub_att0,"ev") #Estimand 
  # 각 집단이 전체 표본에서 차지하는 비중으로 가중 
  temp_sim_portion=sum(temp_md_sub$weights)/sum(MD_class_att$weights)
  # 결과도출 
  sub_sim=tibble(subclass=rep(i,10000),
                 sim_ev=temp_sim*temp_sim_portion,
                 sim_num=1:10000)
  # 집단별로 얻은 결과 합치기 
  all_sub_sim=bind_rows(all_sub_sim,sub_sim) 
}

# park_eva_a 변수 
MD_class_att %>% 
  filter(subclass==6) %>% 
  group_by(treat) %>% 
  summarize_at(
    vars(edu:good_eco),
    function(x){sd(x,na.rm=TRUE)}
  )

# subclass==6 의 경우만 따로 
i=6 
temp_md_sub=MD_class_att %>% filter(subclass==i)
z_sub_att=zelig(y~treat+female+gen20+gen30+gen40+gen50+edu+hhinc+ 
                  libcon+int_eff+good_eco+vote_past,data=temp_md_sub,
                model='ls',weights="weights",cite=F)  
x_sub_att0=setx(z_sub_att,treat=0,data=temp_md_sub)
x_sub_att1=setx(z_sub_att,treat=1,data=temp_md_sub)
s_sub_att0=sim(z_sub_att,x_sub_att0,num=10000)
s_sub_att1=sim(z_sub_att,x_sub_att1,num=10000)
temp_sim=get_qi(s_sub_att1,"ev")-get_qi(s_sub_att0,"ev") #Estimand 
temp_sim_portion=sum(temp_md_sub$weights)/sum(MD_class_att$weights)
sub_sim=tibble(subclass=rep(i,10000),
               sim_ev=temp_sim*temp_sim_portion,
               sim_num=1:10000)
all_sub_sim=bind_rows(all_sub_sim,sub_sim)

# 6개 집단들의 효과추정치 합산
all_sub_sim_agg = all_sub_sim %>% 
  group_by(sim_num) %>% 
  summarize(att=sum(sim_ev))
# 95% 신뢰구간 계산 
SUB_ATT=all_sub_sim_agg$att
myATT = quantile(SUB_ATT,p=c(0.025,0.5,0.975)) %>% 
  data.frame() %>% 
  t() %>% data.frame() 
names(myATT)=c("LL95","PEst","UL95")
myATT$estimand = "ATT"
myATT$model = "Propensity score subclassification (6 groups)"
myATT %>% as_tibble()
# ATC
set.seed(4321)
all_sub_sim=list()           # 집단별 시뮬레이션 결과 저장 
for(i in 1:mysubclassN){
  temp_md_sub=MD_class_atc %>% filter(subclass==i)  #i번째 집단 
  # 효과추정치 추정모형 
  z_sub_atc=zelig(y_R_cov,data=temp_md_sub,
                  model='ls',weights="weights",cite=F)  
  # 독립변수 조건 
  x_sub_atc0=setx(z_sub_atc,Rtreat=1,data=temp_md_sub)
  x_sub_atc1=setx(z_sub_atc,Rtreat=0,data=temp_md_sub)
  # 시뮬레이션 
  s_sub_atc0=sim(z_sub_atc,x_sub_atc0,num=10000)
  s_sub_atc1=sim(z_sub_atc,x_sub_atc1,num=10000)
  temp_sim=get_qi(s_sub_atc1,"ev")-get_qi(s_sub_atc0,"ev") #Estimand 
  # 각 집단이 전체 표본에서 차지하는 비중으로 가중 
  temp_sim_portion=sum(temp_md_sub$weights)/sum(MD_class_atc$weights)
  # 결과도출 
  sub_sim=tibble(subclass=rep(i,10000),
                 sim_ev=temp_sim*temp_sim_portion,
                 sim_num=1:10000)
  # 집단별로 얻은 결과 합치기 
  all_sub_sim=bind_rows(all_sub_sim,sub_sim) 
}
# 6개 집단들의 효과추정치 합산
all_sub_sim_agg = all_sub_sim %>% 
  group_by(sim_num) %>% 
  summarize(atc=sum(sim_ev))
# 95% 신뢰구간 계산 
SUB_ATC=all_sub_sim_agg$atc
myATC = quantile(SUB_ATC,p=c(0.025,0.5,0.975)) %>% 
  data.frame() %>% 
  t() %>% data.frame() 
names(myATC)=c("LL95","PEst","UL95")
myATC$estimand = "ATC"
myATC$model = "Propensity score subclassification (6 groups)"
myATC %>% as_tibble()
# ATE 추정: ATT, ATC, 처치집단비율(pi)을 이용하여 계산
mypi = prop.table(table(mydata$treat))[2]
SUB_ATE = mypi*SUB_ATT + (1-mypi)*SUB_ATC
myATE = quantile(SUB_ATE,p=c(0.025,0.5,0.975))  %>% 
  data.frame() %>% 
  t() %>% data.frame() 
names(myATE)=c("LL95","PEst","UL95")
myATE$estimand = "ATE"
myATE$model = "Propensity score subclassification (6 groups)"
myATE %>% as_tibble() 
# 3가지 효과추정치 통합 
class_estimands = bind_rows(myATT,myATC,myATE) %>% 
  as_tibble()  
class_estimands


## 8) 준-정확매칭 기법 
# ATT, ATC 추정 
MD_cem_att=match.data(cem_att)
MD_cem_atc=match.data(cem_atc)
set.seed(1234)
cem_ATT = SUMMARY_EST_ATT(myformula=y_T_cov,
                          matched_data=MD_cem_att,
                          n_sim=10000,
                          model_name="CEM, user-choice coarsening")
set.seed(4321)
cem_ATC = SUMMARY_EST_ATC(myformula=y_R_cov,
                          matched_data=MD_cem_atc,
                          n_sim=10000,
                          model_name="CEM, user-choice coarsening")
cem_ATE=SUMMARY_EST_ATE(cem_ATT,cem_ATC,mydata$treat,
                        "CEM, user-choice coarsening")
# 3가지 효과추정치(estimands) 통합 
cem_estimands = bind_rows(
  cem_ATT[[2]],cem_ATC[[2]],cem_ATE
)
cem_estimands  

# 효과추정치 비교 시각화 
bind_rows(PSW_estimands,greedy_estimands,optimal_estimands,full_estimands,
          genetic_estimands,mahala_estimands,class_estimands,cem_estimands) %>% 
  mutate(rid=row_number(),
         model=fct_reorder(model,rid)) %>% 
  ggplot(aes(x=estimand,y=PEst,color=model))+
  geom_point(size=3,position=position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=LL95,ymax=UL95),
                width=0.2,lwd=1,
                position=position_dodge(width=0.3))+
  geom_hline(yintercept=OLS_estimands$LL95,lty=2)+
  geom_hline(yintercept=OLS_estimands$UL95,lty=2)+
  geom_hline(yintercept=0,lty=1,color="red")+
  geom_label(x=0.3,y=0.17,
             label="Naive OLS\n95% CI",color="black")+
  labs(x="Estimands",
       y="Estimates, 95% Confidence Interval",
       color="Models")+
  coord_cartesian(xlim=c(0.7,3),ylim=c(-0.05,.3))+
  theme_bw()+theme(legend.position="top")+
  guides(color = guide_legend(nrow=4))
ggsave("Part2_ch9_Comparison_Estimands.png",height=14,width=21,units='cm')

###############################################################################
## 5. 민감도 분석 
###############################################################################
## 1) 성향점수가중 기법 
SA_PSW_ATT = treatSens(formula=y_T_cov,
                       trt.family=binomial(link='probit'),
                       grid.dim=c(7,5),nsim=20,
                       standardize = FALSE,
                       data=mydata,
                       weights=mydata$Watt)
SA_PSW_ATC = treatSens(formula=y_T_cov,
                       trt.family=binomial(link='probit'),
                       grid.dim=c(7,5),nsim=20,
                       standardize = FALSE, 
                       data=mydata,
                       weights=mydata$Watc)
SA_PSW_ATE = treatSens(formula=y_T_cov,
                       trt.family=binomial(link='probit'), 
                       grid.dim=c(7,5),
                       nsim=20,        
                       standardize = FALSE,
                       data=mydata,weights=mydata$Wate) 

png("Part2_ch9_treatSens.png",res=600,height=25,width=12,units="cm")
par(mfrow=c(3,1))
sensPlot(SA_PSW_ATT)
title(main="Carnegie, Harada, Hill's Sensitivity analysis for ATT")
sensPlot(SA_PSW_ATC)
title(main="Carnegie, Harada, Hill's Sensitivity analysis for ATC")
sensPlot(SA_PSW_ATE)
title(main="Carnegie, Harada, Hill's Sensitivity analysis for ATE")
dev.off()


# 제2부 제6장에서 사용한 이용자정의 함수 불러오기 
MATCH_PAIR_GENERATE = readRDS("MATCH_PAIR_GENERATE.RData")
GAMMA_RANGE_SEARCH = readRDS("GAMMA_RANGE_SEARCH.RData")
## 2) 성향점수기반 그리디 매칭 기법 
# ATT 
SA_greedy_att = MATCH_PAIR_GENERATE(greedy_att,mydata) 
SA_greedy_att_gamma_range = GAMMA_RANGE_SEARCH(SA_greedy_att, 1)
tail(SA_greedy_att_gamma_range)
# ATC 
SA_greedy_atc = -1*MATCH_PAIR_GENERATE(greedy_atc,mydata) 
SA_greedy_atc_gamma_range = GAMMA_RANGE_SEARCH(SA_greedy_atc, 1)
tail(SA_greedy_atc_gamma_range)

## 3) 성향점수기반 최적 매칭 기법 
# ATT 
SA_optimal_att = MATCH_PAIR_GENERATE(optimal_att,mydata) 
SA_optimal_att_gamma_range = GAMMA_RANGE_SEARCH(SA_optimal_att, 1)
tail(SA_optimal_att_gamma_range)

# 제2부 제6장에서 사용한 이용자정의 함수 불러오기 
MATCH_PAIR_GENERATE_FULL = readRDS("MATCH_PAIR_GENERATE_FULL.RData")
GAMMA_RANGE_SEARCH_FULL = readRDS("GAMMA_RANGE_SEARCH_FULL.RData")
## 4) 성향점수기반 전체 매칭 기법 
# ATT 
SA_match_data_full_att = MATCH_PAIR_GENERATE_FULL(
  full_att,mydata)
SA_full_att_gamma_range = GAMMA_RANGE_SEARCH_FULL(
  SA_match_data_full_att[[1]],
  SA_match_data_full_att[[2]],1)
tail(SA_full_att_gamma_range)
# ATC 
SA_match_data_full_atc = MATCH_PAIR_GENERATE_FULL(
  full_atc,mydata)
SA_full_atc_gamma_range = GAMMA_RANGE_SEARCH_FULL(
  SA_match_data_full_atc[[1]],
  SA_match_data_full_atc[[2]],1)
tail(SA_full_atc_gamma_range)


## 5) 성향점수기반 유전 매칭 기법 
# ATT 
SA_genetic_att = MATCH_PAIR_GENERATE(genetic_att,mydata) 
SA_genetic_att_gamma_range = GAMMA_RANGE_SEARCH(SA_genetic_att, 1)
tail(SA_genetic_att_gamma_range)
# ATC 
SA_genetic_atc = -1*MATCH_PAIR_GENERATE(genetic_atc,mydata) 
SA_genetic_atc_gamma_range = GAMMA_RANGE_SEARCH(SA_genetic_atc, 1)
tail(SA_genetic_atc_gamma_range)


## 6) 마할라노비스 거리점수기반 그리디 매칭 기법 
# ATT 
SA_mahala_att = MATCH_PAIR_GENERATE(mahala_att,mydata) 
SA_mahala_att_gamma_range = GAMMA_RANGE_SEARCH(SA_mahala_att, 1)
tail(SA_mahala_att_gamma_range)
# ATC 
SA_mahala_atc = -1*MATCH_PAIR_GENERATE(mahala_atc,mydata) 
SA_mahala_atc_gamma_range = GAMMA_RANGE_SEARCH(SA_mahala_atc, 1)
tail(SA_mahala_atc_gamma_range)

## 7) 성향점수층화 기법: 민감도 분석 없음 
## 8) 준정확매칭 기법: 민감도 분석 없음 
