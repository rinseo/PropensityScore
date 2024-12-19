###############################################################################
## 제2부 제10장 범주형 원인변수 대상 성향점수분석 기법 
###############################################################################
library("Hmisc")              #가중평균/가중분산 
library("tidyverse")          #데이터관리 및 변수사전처리 
summarize = dplyr::summarize  # dplyr의 summarize 함수를 이용
library("Zelig")              #비모수접근 95% CI 계산 
library("nnet")               #다항로지스틱 회귀모형 

###############################################################################
## 1. 데이터 소개 
###############################################################################
setwd("D:/data")
mydata=read_csv("observational_study_survey.csv")
mydata=mydata %>% 
  mutate(
    gen20=ifelse(gen=='20s',1,0),
    gen30=ifelse(gen=='30s',1,0),
    gen40=ifelse(gen=='40s',1,0),
    gen50=ifelse(gen=='50s',1,0)
  ) 
# 탄핵찬성(촛불집회참여자)=2, 탄핵반대(태극기집회참여자)=1, 통제집단(집회미참여자)=0 
mydata=mydata %>% 
  mutate(rally3=factor(2*rally_pro+rally_con))
count(mydata, rally3) %>% 
  mutate(pct=100*n/(sum(n)))

################################################################
# 일반화 성향점수(GPS, generalized propensity score) 가중
################################################################
# GPS 추정을 위한 공식정의
covs="(female+gen20+gen30+gen40+gen50+edu+hhinc+libcon+int_eff+park_eva_a+good_eco+vote_past)"
pred_ps3=as.formula(str_c("rally3~",covs))
# 다항로지스틱 회귀모형을 이용하여 GPS 추정 
mnlogit=multinom(pred_ps3, mydata)
# GPS 저장 
gps = fitted(mnlogit) %>% data.frame() 
names(gps)=c("control","con","pro")
head(gps)
# 공통지지영역 탐색 
fig_data = bind_cols(mydata %>% select(rally3),
                     as_tibble(gps)) %>%
  pivot_longer(cols=control:pro,
               names_to="type_GPS",
               values_to="value_GPS") %>%
  mutate(
    rally3=factor(rally3,labels=c("No rally","Korean Flag","Candlelight")),
    type_GPS=fct_relevel(type_GPS,"control","con"),
    type_GPS=factor(type_GPS,labels=c("No rally","Korean Flag","Candlelight"))
  )
fig_data
# 태극기 집회 참여자의 경우 극소수
fig_data %>% ggplot(aes(x=value_GPS,fill=rally3))+
  geom_density(alpha=0.4)+
  labs(x="Treatment conditions",y="Density",fill="Treatment status")+
  theme_bw()+
  theme(legend.position="top")+
  facet_wrap(~type_GPS)
ggsave("Part2_Ch10_RCommonSupport_3groups.jpeg",unit="cm",width=16,height=13)

# 태극기 집회 참여자의 경우 극소수
fig_data %>% ggplot(aes(x=value_GPS,fill=rally3))+
  geom_density(alpha=0.4)+
  coord_cartesian(ylim=c(0,15))+ #Y축 조정 
  labs(x="Treatment conditions",y="Density",fill="Treatment status")+
  theme_bw()+
  theme(legend.position="top")+
  facet_wrap(~type_GPS)
ggsave("Part2_Ch10_RCommonSupport_3groups_Yaxis.jpeg",unit="cm",width=16,height=13)

#처치역확률 가중치 
mydata$IPTW = ifelse(mydata$rally3==0, 1/gps$control, #미참여자의 경우 
                     ifelse(mydata$rally3==1,1/gps$con,
                            1/gps$pro)) #참여자의 경우 촛불 대 태극기로 

#################################################################
# 공변량 균형성(covariate balance) 점검 
#################################################################
# 표준화변환을 위항 이용자정의 함수 
standardizedF=function(myvar){(myvar-mean(myvar))/sd(myvar)}
# 더미변수 형태의 공변량들 선정(표준화 없이 공변량 균형성을 살펴봄)
covs_data_R = mydata %>% 
  select(female,gen20:gen50,vote_past) %>% 
  as.data.frame()
# 연속형변수 형태의 공변량들 선정(표준화 변환 적용)
covs_data_SD = mydata %>% 
  select(edu:good_eco) %>% 
  mutate_all(
    standardizedF
  ) %>% 
  as.data.frame()
# 공변량 데이터로 합치기(더미 + 연속형 변수)
covs_data=bind_cols(covs_data_R,covs_data_SD)
# 가중평균후 평균차이(절댓값), 분산비 산출위한 이용자정의 함수 
mybalanceF=function(i){
  M0=wtd.mean(covs_data[,i][mydata$rally3==0],mydata$IPTW[mydata$rally3==0])
  M1=wtd.mean(covs_data[,i][mydata$rally3==1],mydata$IPTW[mydata$rally3==1])
  M2=wtd.mean(covs_data[,i][mydata$rally3==2],mydata$IPTW[mydata$rally3==2])
  V0=wtd.var(covs_data[,i][mydata$rally3==0],mydata$IPTW[mydata$rally3==0])
  V1=wtd.var(covs_data[,i][mydata$rally3==1],mydata$IPTW[mydata$rally3==1])
  V2=wtd.var(covs_data[,i][mydata$rally3==2],mydata$IPTW[mydata$rally3==2])
  MD_01=abs(M0-M1)
  MD_02=abs(M0-M2)
  MD_12=abs(M1-M2)
  VR_01=max(V0,V1)/min(V0,V1)
  VR_02=max(V0,V2)/min(V0,V2)
  VR_12=max(V2,V1)/min(V2,V1)
  cov_name=names(covs_data)[i]
  tibble(
    cov_name,MD_01,MD_02,MD_12,VR_01,VR_02,VR_12
  )
}
# mybalanceF() 함수 적용(반복계산) 
mysummary=list()
for (mylocation in 1:dim(covs_data)[2]){
  mysummary = bind_rows(mysummary,mybalanceF(mylocation))
}
mysummary
# 평균차이(절댓값 기준) 허용범위를 넘는 공변량 비율
mysummary %>% 
  select(contains("MD_")) %>% 
  mutate_all(function(x){ifelse(x > .25,1,0)}) %>% 
  summarize_all(
    sum
  ) 
# 분산비 허용범위를 넘는 공변량 
mysummary %>% 
  filter(row_number()>6) %>%  #더미변수 형태 공변량들은 포함하지 않았음
  select(contains("VR_")) %>% 
  mutate_all(function(x){ifelse(x > 2,1,0)}) %>% 
  summarize_all(
    sum
  ) 

# 처치효과 추정
# 공식정의
set.seed(1234)
pred_y3=as.formula(str_c("vote_will~rally_con+rally_pro+",covs))
z_out_ATE=zelig(pred_y3,data=mydata,
                model="ls",weights="IPTW",cite=FALSE)
# 원인변수 조건 지정 
x_control=setx(z_out_ATE,rally_pro=0,rally_con=0,data=mydata) 
x_pro=setx(z_out_ATE,rally_pro=1,rally_con=0,data=mydata)
x_con=setx(z_out_ATE,rally_pro=0,rally_con=1,data=mydata)
# 시뮬레이션
s_control=sim(z_out_ATE,x_control,num=10000)
s_pro=sim(z_out_ATE,x_pro,num=10000)
s_con=sim(z_out_ATE,x_con,num=10000)
# 점추정치 및 95%신뢰구간 계산 
Compare_control_pro=get_qi(s_pro,"ev") - get_qi(s_control,"ev")
Compare_control_con=get_qi(s_con,"ev") - get_qi(s_control,"ev")
Compare_pro_con=get_qi(s_pro,"ev") - get_qi(s_con,"ev")
PSW_compare = rbind(
  quantile(Compare_control_pro,p=c(0.025,0.5,0.975)),
  quantile(Compare_control_con,p=c(0.025,0.5,0.975)),
  quantile(Compare_pro_con,p=c(0.025,0.5,0.975))
) %>% as_tibble()
names(PSW_compare) = c("LL95","PEst","UL95")
PSW_compare$comparison = c("Candlelight-No","Flag-No","Candlelight-Flag")
PSW_compare$model="Propensity score weighting using GPS"
PSW_compare

# 시각화 
PSW_compare %>% 
  ggplot(aes(x=comparison,y=PEst))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LL95,ymax=UL95),
                width=0.1,lwd=0.5,position=position_dodge(width=0.2))+
  geom_hline(yintercept=0,lty=2,lwd=0.1,color='red')+
  labs(x="Comparison between groups",
       y="Estimates, 95% Confidence interval")+
  coord_cartesian(ylim=c(-0.4,0.4))+
  theme_bw()+
  theme(legend.position="top")
ggsave("Part2_Ch10_PSW_GPS.jpeg",unit="cm",width=11,height=9)

# 통상적 회귀분석 추정결과와의 비교 
# 미참여자 vs 촛불집회; 미참여자 vs 태극기집회
OLS_control_pro_cons = lm(vote_will ~ rally_pro + rally_con + female + 
                            gen20 + gen30 + gen40 + gen50 + edu + hhinc + libcon + 
                            int_eff + park_eva_a + good_eco + vote_past, mydata)
OLS_pro_control=c(OLS_control_pro_cons$coefficients["rally_pro"],
                  confint(OLS_control_pro_cons)["rally_pro",])
OLS_con_control=c(OLS_control_pro_cons$coefficients["rally_con"],
                  confint(OLS_control_pro_cons)["rally_con",])
# 촛불집회 vs 태극기집회; 촛불집회 vs 미참여자 
mydata$rally_no = ifelse(mydata$rally3==0,1,0)
OLS_pro_cons_control = lm(vote_will ~ rally_pro + rally_no + female + 
                            gen20 + gen30 + gen40 + gen50 + edu + hhinc + libcon + 
                            int_eff + park_eva_a + good_eco + vote_past, mydata)
OLS_pro_con=c(OLS_pro_cons_control$coefficients["rally_pro"],
              confint(OLS_pro_cons_control)["rally_pro",])
OLS_compare=rbind(OLS_pro_control,OLS_con_control,OLS_pro_con) %>% as_tibble()
names(OLS_compare)=c("PEst","LL95","UL95")
OLS_compare$comparison = c("Candlelight-No","Flag-No","Candlelight-Flag")
OLS_compare$model="Conventional OLS"
OLS_compare = OLS_compare %>% select(names(PSW_compare)) #변수의 순서 맞추기 
OLS_compare

# 시각화 
bind_rows(PSW_compare,OLS_compare) %>% 
  ggplot(aes(x=comparison,y=PEst,shape=model,color=model))+
  geom_point(size=2,position=position_dodge(width=0.2))+
  geom_errorbar(aes(ymin=LL95,ymax=UL95),
                width=0.1,lwd=0.5,position=position_dodge(width=0.2))+
  geom_hline(yintercept=0,lty=2,lwd=0.1,color='red')+
  labs(x="Comparison between groups",
       y="Estimates, 95% Confidence interval",
       shape="Model",color="Model")+
  coord_cartesian(ylim=c(-0.4,0.4))+
  theme_bw()+
  theme(legend.position="top")
ggsave("Part2_Ch10_estimand_compare.jpeg",unit="cm",width=13,height=11)
