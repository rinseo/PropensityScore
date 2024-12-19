# 매칭기법별 비교 
# 시뮬레이션 데이터 생성은 Matching_Simulated_Data.r 참조. 시드넘버=1234
#################################################################
## 제2부 제5장 성향점수가중(propensity score weighting, PSW) 기법 
#################################################################
library("tidyverse")   #데이터관리 및 변수사전처리 
summarize = dplyr::summarize   # dplyr의 summarize 함수를 이용
library("Hmisc")       #가중평균 및 가중분산 계산 
library("Zelig")   #비모수접근 95% CI 계산 
# 2020년 현재 treatSens 는 CRAN에서 관리하지 않음.
# 아래와 같은 방식으로 깃허브를 통해 인스톨한 후 구동 
# remotes::install_github("vdorie/treatSens")
library("treatSens")

# 데이터 소환
setwd("D:/data")
mydata=read_csv("simdata.csv")
mydata 

# 성향점수 계산 
mydata$pscore = glm(treat~V1+V2+V3,data=mydata,
                    family=binomial(link='logit')) %>% 
  fitted()
# IPTW 계산(ATE, ATT, ATC 추정) 
mydata = mydata %>% 
  mutate(
    Wate = ifelse(treat==1,1/pscore,1/(1-pscore)),
    Watt = ifelse(treat==1,1,1/(1-pscore)),
    Watc = ifelse(treat==1,1/pscore,1)
  )

# 각주: 성향점수가중치 완전절단(truncation)
quantile(mydata$Wate,prob=0.01*(90:100))
# Gurel & Leite (2012) 극단적 IPTW 조정방법 
mydata$WateTrun = ifelse(mydata$Wate>quantile(mydata$Wate,prob=0.99),
                         quantile(mydata$Wate,prob=0.99),
                         mydata$Wate) 
ggplot(mydata, aes(x=Wate,y=WateTrun))+
  geom_point()+
  labs(x="Raw weight",y="Truncated weight")+
  theme_bw()
ggsave("Part2_ch5_weight_truncation.png",width=8,height=6,unit='cm')

# 공변량 균형성 점검
# 공변량 표준화 
mydata$V1S = (mydata$V1-mean(mydata$V1))/sd(mydata$V1)
# IPTW 적용전과 적용후 평균차이 비교 
mean(mydata$V1S[mydata$treat==1]) - mean(mydata$V1S[mydata$treat==0])
wtd.mean(mydata$V1S[mydata$treat==1],mydata$Wate[mydata$treat==1]) - 
  wtd.mean(mydata$V1S[mydata$treat==0],mydata$Wate[mydata$treat==0])
# IPTW 적용전과 적용후 분산비 비교 
var(mydata$V1S[mydata$treat==1])/var(mydata$V1S[mydata$treat==0])
wtd.var(mydata$V1S[mydata$treat==1],mydata$Wate[mydata$treat==1])/
  wtd.var(mydata$V1S[mydata$treat==0],mydata$Wate[mydata$treat==0])

# 균형성 점검을 위한 이용자함수 설정
balance_check_PSW = function(var_treat,var_cov,var_wgt){
  std_var_cov = (var_cov - mean(var_cov))/sd(var_cov)
  simple_M1 = mean(std_var_cov[var_treat==1])
  simple_M0 = mean(std_var_cov[var_treat==0])
  simple_V1 = var(std_var_cov[var_treat==1])
  simple_V0 = var(std_var_cov[var_treat==0])
  wgted_M1 = Hmisc::wtd.mean(x=std_var_cov[var_treat==1],weights=var_wgt[var_treat==1])
  wgted_M0 = Hmisc::wtd.mean(x=std_var_cov[var_treat==0],weights=var_wgt[var_treat==0])
  wgted_V1 = Hmisc::wtd.var(x=std_var_cov[var_treat==1],weights=var_wgt[var_treat==1])
  wgted_V0 = Hmisc::wtd.var(x=std_var_cov[var_treat==0],weights=var_wgt[var_treat==0])
  B_wgt_Mdiff = simple_M1 - simple_M0
  B_wgt_Vratio = simple_V1/simple_V0
  A_wgt_Mdiff = wgted_M1 - wgted_M0
  A_wgt_Vratio = wgted_V1/wgted_V0
  balance_index = tibble(B_wgt_Mdiff,A_wgt_Mdiff,
                         B_wgt_Vratio,A_wgt_Vratio)
  balance_index
}
saveRDS(balance_check_PSW,"balance_check_PSW.RData")
# 성향점수와 공변량의 균형성 점검 
# balance_check_PSW(mydata$treat,mydata$V1,mydata$Wate) # 앞에서 확인함 
balance_check_PSW(mydata$treat,mydata$V2,mydata$Wate)
balance_check_PSW(mydata$treat,mydata$V3,mydata$Wate)
balance_check_PSW(mydata$treat,mydata$pscore,mydata$Wate)

# 모수통계기법으로 효과추정치 추정 
lm(y~treat+V1+V2+V3,mydata,weights=Wate)$coef['treat'] #ATE  
lm(y~treat+V1+V2+V3,mydata,weights=Wate) %>% 
  confint("treat")   # ATE, 95% CI  

# 효과추정치 정리 
ate=c(lm(y~treat+V1+V2+V3,mydata,weights=Wate)$coef['treat'],
      lm(y~treat+V1+V2+V3,mydata,weights=Wate) %>% confint("treat"))
att=c(lm(y~treat+V1+V2+V3,mydata,weights=Watt)$coef['treat'],
      lm(y~treat+V1+V2+V3,mydata,weights=Watt) %>% confint("treat"))
atc=c(lm(y~treat+V1+V2+V3,mydata,weights=Watc)$coef['treat'],
      lm(y~treat+V1+V2+V3,mydata,weights=Watc) %>% confint("treat"))
estimands_psw = data.frame(rbind(ate,att,atc))
names(estimands_psw)=c("PEst","LL95","UL95") 
estimands_psw %>% select(LL95,PEst,UL95) %>% round(3)

# 비모수통계기법으로 효과추정치 추정 
# ATE 추정 
set.seed(1234)  # 동일한 결과를 얻고자 한다면.... 
z_model=zelig(y~treat+V1+V2+V3,data=mydata,model='ls',
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
summary_est_ate

# ATT 추정 
set.seed(1234)  # 동일한 결과를 얻고자 한다면.... 
z_model=zelig(y~treat+V1+V2+V3,data=mydata,model='ls', #'ls'는 OLS를 의미함
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
z_model=zelig(y~treat+V1+V2+V3,data=mydata,model='ls',
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
# 3가지 효과추정치(estimands) 통합 
PSW_estimands = bind_rows(
  summary_est_att,summary_est_atc,summary_est_ate
)
PSW_estimands
saveRDS(PSW_estimands,"PSW_estimands.RData")

# 시각화 
PSW_estimands %>% 
  ggplot(aes(x=estimand,y=PEst))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=LL95,ymax=UL95),width=0.1,lwd=1)+
  labs(x="Estimands",y="Estimates, 95% Confidence Interval")+
  coord_cartesian(ylim=c(0.5,2.5))+
  theme_bw()+
  ggtitle("propensity score weighting")
ggsave("Part2_ch5_PSW_estimands.png",width=10,height=10,unit='cm')

# 민감도 분석(Carnegie, Harada, & Hill, 2006)
library("treatSens")
SA_PSW_ATE = treatSens(y~treat+V1+V2+V3,
                       trt.family=binomial(link='probit'), #이분변수인 경우 현재 유일하게 제공되는 옵션 
                       grid.dim=c(7,5),   # zeta_T를 7개, zeta_Y를 5개(3번째 결과물)
                       nsim=20,           # 재표집 횟수 
                       standardize = FALSE,   # 누락변수의 효과는 비표준화 회귀계수로 
                       data=mydata,weights=mydata$Wate)  # IPTW 부여 
summary(SA_PSW_ATE)
pnorm(0,0,1)-pnorm(-2,0,1) # zeta_T (zeta_Z) 해석 

SA_PSW_ATT = treatSens(y~treat+V1+V2+V3,
                       trt.family=binomial(link='probit'),
                       grid.dim=c(7,5),nsim=20,
                       standardize = FALSE, 
                       data=mydata,
                       weights=mydata$Watt)
summary(SA_PSW_ATT)

SA_PSW_ATC = treatSens(y~treat+V1+V2+V3,
                       trt.family=binomial(link='probit'),
                       grid.dim=c(7,5),nsim=20,
                       standardize = FALSE, 
                       data=mydata,
                       weights=mydata$Watc)
summary(SA_PSW_ATC)

png("Part2_ch5_sensPlot.png",res=600,height=10,width=15,units="cm")
sensPlot(SA_PSW_ATC)
dev.off()
