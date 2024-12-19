# 시율레이션 데이터 생성은 Matching_Simulated_Data.r 참조. 시드넘버=1234
###############################################################################
## 제2부 제8장 시뮬레이션 데이터 대상 준-정확매칭(CEM, Coarsened exact matching)  
###############################################################################
library("tidyverse")   #데이터관리 및 변수사전처리 
summarize = dplyr::summarize   # dplyr의 summarize 함수를 이용
library("Zelig")   #비모수접근 95% CI 계산 
library("MatchIt")   #매칭기법 적용 
library("cobalt")   #처치집단과 통제집단 균형성 점검 

# 데이터 소환
setwd("D:/data")
mydata=read_csv("simdata.csv")
mydata 
###############################################################################
## CEM(subclassification)
###############################################################################
###############################################################################
## 2. 준정확매칭 실시
###############################################################################
## 자동화 뭉뚱그리기(automated coarsening)
# ATT, 자동화 뭉뚱그리기(automated coarsening) 
set.seed(1234)
cem_a_att=matchit(formula=treat~V1+V2+V3,
                data=mydata,
                method="cem",
                discard='none',
                L1.breaks='scott') # 자동화 뭉뚱그리기 옵션(스콧의 방법)
cem_a_att
# 점진적 뭉뚱그리기(progressive coarsening)
set.seed(1234)
K = 10
cem_p_att=matchit(formula=treat~V1+V2+V3,
                  data=mydata,
                  method="cem",
                  discard='none',
                  cutpoints=list(V1=K,V2=K,V3=K)) # 각 공변량을 k개 집단으로 
cem_p_att

# 점진적 뭉뚱그리기를 통한 처치집단 사례수 변화 
set.seed(1234)
Tcases_Included=rep(NA,9)
for (K in 2:10){
  m_cem=matchit(formula=treat~V1+V2+V3,
                data=mydata,
                method="cem",
                discard='none',
                cutpoints=list(V1=K,V2=K,V3=K)) # 각 공변량을 k개 집단으로 
  Tcases_Included[K-1]=m_cem$nn[2,2]  #처치집단 개수 
}
# K변화에 따른 매칭작업 포함 처치집단 개수 변화 
tibble(K=2:10,Cases=Tcases_Included) %>% 
  ggplot()+
  geom_point(aes(x=K, y=Cases))+
  geom_line(aes(x=K, y=Cases))+
  scale_x_continuous(breaks=2:10)+
  scale_y_continuous(breaks=Tcases_Included)+
  coord_cartesian(ylim=c(140,160))+
  labs(x="\nK in Progressive coarsening",
       y="Treated cases matched\n")+
  theme_bw()
ggsave("Part2_ch8_pCoarsening_K.jpeg",height=10,width=14,units='cm')

# ATC, 자동화 뭉뚱그리기(automated coarsening) 
set.seed(4321)
cem_a_atc=matchit(formula=Rtreat~V1+V2+V3,
                  data=mydata,
                  method="cem",
                  discard='none',
                  L1.breaks='scott') # 자동화 뭉뚱그리기 옵션(스콧의 방법)
cem_a_atc

###############################################################################
## 3. 공변량 균형성 점검 
###############################################################################
# ATT, 자동화 뭉뚱그리기(automated coarsening), 러브플롯 
# 평균차이 
love_D_cem_att=love.plot(cem_a_att,
                         s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                         stat="mean.diffs", # 두 집단간 공변량 평균차이 
                         drop.distance=TRUE,
                         threshold=0.1, # 평균차이 역치 
                         sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                         themes=theme_bw())+
  coord_cartesian(xlim=c(-0.75,1.00))+ #X축의 범위 
  ggtitle("ATT, Covariate balance")
# 분산비 
love_VR_cem_att=love.plot(cem_a_att,
                          s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                          stat="variance.ratios", # 두 집단간 분산비 
                          drop.distance=TRUE,
                          threshold=2, # 분산비 역치 
                          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                          themes=theme_bw())+
  coord_cartesian(xlim=c(0.3,3))+ #X축의 범위 
  ggtitle("ATT, Covariate balance")
# ATC, 자동화 뭉뚱그리기(automated coarsening), 러브플롯 
# 평균차이 
love_D_cem_atc=love.plot(cem_a_atc,
                         s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                         stat="mean.diffs", # 두 집단간 공변량 평균차이 
                         drop.distance=TRUE, 
                         threshold=0.1, # 평균차이 역치 
                         sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                         themes=theme_bw())+
  coord_cartesian(xlim=c(-0.75,1.00))+ #X축의 범위 
  ggtitle("ATC, Covariate balance")
# 분산비 
love_VR_cem_atc=love.plot(cem_a_atc,
                          s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                          stat="variance.ratios", # 두 집단간 분산비 
                          drop.distance=TRUE,
                          threshold=2, # 분산비 역치 
                          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                          themes=theme_bw())+
  coord_cartesian(xlim=c(0.3,3))+ #X축의 범위 
  ggtitle("ATC, Covariate balance")
png("Part2_ch8_love_MD_VR_cem_att_atc.png",res=600,height=15,width=25,units="cm")
gridExtra::grid.arrange(love_D_cem_att,love_VR_cem_att,
                        love_D_cem_atc,love_VR_cem_atc,
                        nrow=2)
dev.off()


# ATT, 점진적 뭉뚱그리기를 통한 균형성 달성수준 변화 
set.seed(1234)
abs_M_diff=matrix(NA,nrow=9,ncol=3)  #공변량 평균차이
colnames(abs_M_diff)=str_c("Mean_diff_V",1:3)
rownames(abs_M_diff)=str_c("K",2:10)
Tcases_Included=rep(NA,9)            #매칭된 처치집단 사례수 
for (K in 2:10){
  m_cem=matchit(formula=treat~V1+V2+V3,
                data=mydata,
                method="cem",
                discard='none',
                cutpoints=list(V1=K,V2=K,V3=K)) # 각 공변량을 k개 집단으로 
  # 아래는 표준화 평균차이(절대값), MatchIt 패키지 부속함수를 사용하였음 
  abs_M_diff[K-1,]=abs(as.vector(summary(m_cem,stanadardized=TRUE)$sum.matched[2:4,4]))
  # 아래는 매칭된 처치집단 사례수 
  Tcases_Included[K-1]=m_cem$nn[2,2]  #처치집단 개수
}
CEM_optimalK = as_tibble(abs_M_diff)
CEM_optimalK$K=2:10 
CEM_optimalK$Treated=Tcases_Included
# 시각화
Fig_covbal=CEM_optimalK %>% 
  pivot_longer(cols=Mean_diff_V1:Mean_diff_V3, names_to="Mean_diff") %>% 
  mutate(Mean_diff=str_remove(Mean_diff, "Mean_diff_")) %>% 
  ggplot()+
  geom_point(aes(x=K,y=value,color=Mean_diff))+
  geom_line(aes(x=K,y=value,color=Mean_diff))+
  geom_hline(yintercept=0.1,lty=2,color='red')+
  scale_x_continuous(breaks=2:10)+
  coord_cartesian()+
  labs(x="\nLevels of covaiates by progressive coarsening",
       y="Absolute mean difference of covariate \n between treated and control groups\n",
       color="Covariate")+
  theme_bw()+
  theme(legend.position="top")
Fig_TreatedCases=CEM_optimalK %>% 
  ggplot()+
  geom_point(aes(x=K,y=Treated))+
  geom_line(aes(x=K,y=Treated))+
  scale_x_continuous(breaks=2:10)+
  scale_y_continuous(breaks=Tcases_Included)+
  coord_cartesian()+
  labs(x="\nLevels of covaiates by progressive coarsening",
       y="Number of matched treated cases\n\n")+
  theme_bw()
png("Part2_ch8_TreatedCases_CovBal_cem_att.png",res=600,
    height=15,width=20,units="cm")
gridExtra::grid.arrange(Fig_covbal,Fig_TreatedCases,nrow=2)
dev.off()

# ATC, 점진적 뭉뚱그리기를 통한 균형성 달성수준 변화 
set.seed(4321)
abs_M_diff=matrix(NA,nrow=9,ncol=3)  #공변량 평균차이
colnames(abs_M_diff)=str_c("Mean_diff_V",1:3)
rownames(abs_M_diff)=str_c("K",2:10)
Tcases_Included=rep(NA,9)            #매칭된 처치집단 사례수 
for (K in 2:10){
  m_cem=matchit(formula=Rtreat~V1+V2+V3,
                data=mydata,
                method="cem",
                discard='none',
                cutpoints=list(V1=K,V2=K,V3=K)) # 각 공변량을 k개 집단으로 
  # 아래는 표준화 평균차이(절대값), MatchIt 패키지 부속함수를 사용하였음 
  abs_M_diff[K-1,]=abs(as.vector(summary(m_cem,stanadardized=TRUE)$sum.matched[2:4,4]))
  # 아래는 매칭된 처치집단 사례수 
  Tcases_Included[K-1]=m_cem$nn[2,2]  #처치집단 개수
}
CEM_optimalK = as_tibble(abs_M_diff)
CEM_optimalK$K=2:10 
CEM_optimalK$Treated=Tcases_Included
# 시각화
Fig_covbal=CEM_optimalK %>% 
  pivot_longer(cols=Mean_diff_V1:Mean_diff_V3, names_to="Mean_diff") %>% 
  mutate(Mean_diff=str_remove(Mean_diff, "Mean_diff_")) %>% 
  ggplot()+
  geom_point(aes(x=K,y=value,color=Mean_diff))+
  geom_line(aes(x=K,y=value,color=Mean_diff))+
  geom_hline(yintercept=0.1,lty=2,color='red')+
  scale_x_continuous(breaks=2:10)+
  coord_cartesian()+
  labs(x="\nLevels of covaiates by progressive coarsening",
       y="Absolute mean difference of covariate \n between control and treated groups\n",
       color="Covariate")+
  theme_bw()+
  theme(legend.position="top")
Fig_TreatedCases=CEM_optimalK %>% 
  ggplot()+
  geom_point(aes(x=K,y=Treated))+
  geom_line(aes(x=K,y=Treated))+
  scale_x_continuous(breaks=2:10)+
  scale_y_continuous(breaks=Tcases_Included)+
  coord_cartesian()+
  labs(x="\nLevels of covaiates by progressive coarsening",
       y="Number of matched\ntreated (actually control) cases\n")+
  theme_bw()
png("Part2_ch8_TreatedCases_CovBal_cem_atc.png",res=600,
    height=15,width=20,units="cm")
gridExtra::grid.arrange(Fig_covbal,Fig_TreatedCases,nrow=2)
dev.off()

# ATT, 점진적 뭉뚱그리기(progressive coarsening)
K = 5
cem_p_att=matchit(formula=treat~V1+V2+V3,
                  data=mydata,
                  method="cem",
                  discard='none',
                  cutpoints=list(V1=K,V2=K,V3=K)) # 각 공변량을 k개 집단으로 
# ATC, 점진적 뭉뚱그리기(progressive coarsening)
cem_p_atc=matchit(formula=Rtreat~V1+V2+V3,
                  data=mydata,
                  method="cem",
                  discard='none',
                  cutpoints=list(V1=K,V2=K,V3=K)) # 각 공변량을 k개 집단으로 

###############################################################################
## 4. 처치효과 추정 
###############################################################################
# 이용자정의 함수 불러오기 
SUMMARY_EST_ATT=readRDS("SUMMARY_EST_ATT.RData")
SUMMARY_EST_ATC=readRDS("SUMMARY_EST_ATC.RData")
# ATT, 자동화 뭉뚱그리기(automated coarsening)
# 1단계: 매칭 데이터 생성 
MD_cem_a_att=match.data(cem_a_att)
# 2-5단계: ATT 추정 
set.seed(1234)
CEM_a_ATT = SUMMARY_EST_ATT(myformula="y~treat+V1+V2+V3",
                           matched_data=MD_cem_a_att,
                           n_sim=10000,
                           model_name="CEM, automated coarsening")
# ATC, 자동화 뭉뚱그리기(automated coarsening)
# 1)ATC-성향점수기반 전체 매칭
# 1단계: 매칭 데이터 생성 
MD_cem_a_atc=match.data(cem_a_atc)
# 2-5단계: ATC 추정 
set.seed(4321)
CEM_a_ATC = SUMMARY_EST_ATC(myformula="y~Rtreat+V1+V2+V3",
                           matched_data=MD_cem_a_atc,
                           n_sim=10000,
                           model_name="CEM, automated coarsening")
CEM_a_ATT[[2]]  #95% 신뢰구간 
CEM_a_ATC[[2]]  #95% 신뢰구간 

# 6단계: ATE 점추정치와 95% CI
# 표본내 처치집단 비율 
mypi = prop.table(table(mydata$treat))[2]
CEM_a_ATE=list()
CEM_a_ATE[[1]] = mypi*CEM_a_ATT[[1]] + (1-mypi)*CEM_a_ATC[[1]]
CEM_a_ATE[[2]] = tibble(
  LL95=quantile(CEM_a_ATE[[1]],p=c(0.025)),
  PEst=quantile(CEM_a_ATE[[1]],p=c(0.500)),
  UL95=quantile(CEM_a_ATE[[1]],p=c(0.975)),
  estimand="ATE",model="CEM, automated coarsening")
CEM_a_ATE[[2]]


# ATT, K=5, 점진적 뭉뚱그리기(progressive coarsening)
# 1단계: 매칭 데이터 생성 
MD_cem_p_att=match.data(cem_p_att)
# 2-5단계: ATT 추정 
set.seed(1234)
CEM_p_ATT = SUMMARY_EST_ATT(myformula="y~treat+V1+V2+V3",
                            matched_data=MD_cem_p_att,
                            n_sim=10000,
                            model_name="CEM, 5 levels per cov.")
CEM_p_ATT[[2]]  #95% 신뢰구간 

# ATC, K=5, 점진적 뭉뚱그리기(progressive coarsening)
# 1)ATC-성향점수기반 전체 매칭
# 1단계: 매칭 데이터 생성 
MD_cem_p_atc=match.data(cem_p_atc)
# 2-5단계: ATC 추정 
set.seed(4321)
CEM_p_ATC = SUMMARY_EST_ATC(myformula="y~Rtreat+V1+V2+V3",
                            matched_data=MD_cem_p_atc,
                            n_sim=10000,
                            model_name="CEM, 5 levels per cov.")
CEM_p_ATC[[2]]  #95% 신뢰구간 

# 6단계: ATE 점추정치와 95% CI
# 표본내 처치집단 비율 
mypi = prop.table(table(mydata$treat))[2]
CEM_p_ATE=list()
CEM_p_ATE[[1]] = mypi*CEM_p_ATT[[1]] + (1-mypi)*CEM_p_ATC[[1]]
CEM_p_ATE[[2]] = tibble(
  LL95=quantile(CEM_p_ATE[[1]],p=c(0.025)),
  PEst=quantile(CEM_p_ATE[[1]],p=c(0.500)),
  UL95=quantile(CEM_p_ATE[[1]],p=c(0.975)),
  estimand="ATE",model="CEM, 5 levels per cov.")
CEM_p_ATE[[2]]

# 효과추정치들 저장 
CEM_a_estimands=bind_rows(CEM_a_ATT[[2]],CEM_a_ATC[[2]],CEM_a_ATE[[2]])
CEM_p_estimands=bind_rows(CEM_p_ATT[[2]],CEM_p_ATC[[2]],CEM_p_ATE[[2]])
# 앞서 추정한 효과추정치들 불러오기 
OLS_estimands=readRDS("OLS_estimands.RData" )
PSW_estimands=readRDS("PSW_estimands.RData")
greedy_estimands=readRDS("greedy_estimands.RData")
optimal_estimands=readRDS("optimal_estimands.RData")
full_estimands=readRDS("full_estimands.RData")
genetic_estimands=readRDS("genetic_estimands.RData")
mahala_estimands=readRDS("mahala_estimands.RData")
class_estimands=readRDS("class_estimands.RData")
# 효과추정치 시각화 
bind_rows(PSW_estimands,greedy_estimands,
          optimal_estimands,full_estimands,
          genetic_estimands,mahala_estimands,
          class_estimands,CEM_a_estimands,CEM_p_estimands) %>% 
  mutate(rid=row_number(),
         model=fct_reorder(model,rid)) %>% 
  ggplot(aes(x=estimand,y=PEst,color=model))+
  geom_point(size=3,position=position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=LL95,ymax=UL95),
                width=0.2,lwd=1,
                position=position_dodge(width=0.3))+
  geom_hline(yintercept=OLS_estimands$LL95,lty=2)+
  geom_hline(yintercept=OLS_estimands$UL95,lty=2)+
  geom_label(x=0.7,y=0.5*(OLS_estimands$LL95+OLS_estimands$UL95),
             label="Naive OLS\n95% CI",color="black")+
  labs(x="Estimands",
       y="Estimates, 95% Confidence Interval",
       color="Models")+
  coord_cartesian(ylim=c(0.5,2.5))+
  theme_bw()+theme(legend.position="top")+
  guides(color = guide_legend(nrow=5))
ggsave("Part2_ch8_Comparison_Estimands.png",height=17,width=17,units='cm')
