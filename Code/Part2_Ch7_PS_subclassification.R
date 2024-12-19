# 시율레이션 데이터 생성은 Matching_Simulated_Data.r 참조. 시드넘버=1234
###############################################################################
## 제2부 제7장 시뮬레이션 데이터 대상 성향점수층화 기법 
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
## 성향점수 층화(subclassification)
###############################################################################
###############################################################################
## 2. 성향점수층화 기법 실시 
###############################################################################
mysubclassN=8  # 집단수에 따라 추정결과가 조금씩 달라짐(디폴트는 6)
# ATT 추정
set.seed(1234)
class_att=matchit(formula=treat~V1+V2+V3,
                  data=mydata,
                  distance="linear.logit",
                  method="subclass",   # 층화 기법 
                  discard='none', 
                  subclass=mysubclassN)
class_att
# ATC 추정
set.seed(4321)
class_atc=matchit(formula=Rtreat~V1+V2+V3,
                  data=mydata,
                  distance="linear.logit",
                  method="subclass",   # 층화 기법 
                  discard='none', 
                  subclass=mysubclassN)
class_atc

###############################################################################
## 3. 공변량 균형성 점검 
###############################################################################
## 공변량 균형성 점검: cobalt 패키지(4.0.0)의 경우 러브플롯이 지원되지 않음 
# ATT-성향점수층화 기법, 평균차이/분산비 
bal.tab(class_att,
        continuous="std",s.d.denom="pooled",
        m.threshold=0.1,v.threshold=2,disp.subclass=T)
# ATC-성향점수층화 기법, 평균차이 /분산비 
bal.tab(class_atc,
        continuous="std",s.d.denom="pooled",
        m.threshold=0.1,v.threshold=2)

###############################################################################
## 4. 처치효과 추정
###############################################################################
# ATT추정 
MD_class_att=match.data(class_att)
head(MD_class_att) %>% round(2)
# 집단별 ATT추정방법(예시) 
MD_class_att %>% filter(subclass==4) %>% #4번 집단의 경우 
  lm(y~treat+V1+V2+V3,weights=weights,.)

# 모수통계기법 
# 효과추정치 for 루프 
class8_att = rep(NA,8)
for (i in 1:8){
  temp=MD_class_att %>% 
    filter(subclass==i) %>% 
    lm(y~treat+V1+V2+V3,weights=weights,.) #처치효과 추정모형 
  class8_att[i] = temp$coef['treat']  # 처치효과 투입 
}
class8_att %>% round(3)

# 가중치 for 루프 
wgt8_prop = rep(NA,8)
wgt_total=sum(MD_class_att$weights)  
for (i in 1:8){
  wgt_class=sum(MD_class_att$weights[MD_class_att$subclass==i]) #비중
  wgt8_prop[i]=wgt_class/wgt_total   #가중치 
}
wgt8_prop %>% round(3)

# 전체표본에서 얻을 것으로 기대되는 처치효과 
tibble(class8_att,wgt8_prop) %>% 
  mutate(att_wgt=class8_att*wgt8_prop) %>% 
  summarize(sum(att_wgt))

# 비모수통계기법 ATT 추정 
set.seed(1234)
all_sub_sim=list()           # 집단별 시뮬레이션 결과 저장 
for(i in 1:mysubclassN){
  temp_md_sub=MD_class_att %>% filter(subclass==i)  #i번째 집단 
  # 효과추정치 추정모형 
  z_sub_att=zelig(y~treat+V1+V2+V3,data=temp_md_sub,
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
# 8개 집단들의 효과추정치 합산
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
myATT$model = "Propensity score subclassification (8 groups)"
myATT %>% as_tibble()


# 비모수통계기법 ATC 추정 
MD_class_atc=match.data(class_atc)
set.seed(1234)
all_sub_sim=list()           # 집단별 시뮬레이션 결과 저장 
for(i in 1:mysubclassN){
  temp_md_sub=MD_class_atc %>% filter(subclass==i)  #i번째 집단 
  # 효과추정치 추정모형 
  z_sub_atc=zelig(y~Rtreat+V1+V2+V3,data=temp_md_sub,
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
# 8개 집단들의 효과추정치 합산
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
myATC$model = "Propensity score subclassification (8 groups)"
myATC %>% as_tibble()

# ATE 추정: ATT, ATC, 처치집단비율(pi)을 이용하여 계산
mypi = prop.table(table(mydata$treat))[2]
SUB_ATE = mypi*SUB_ATT + (1-mypi)*SUB_ATC
myATE = quantile(SUB_ATE,p=c(0.025,0.5,0.975))  %>% 
  data.frame() %>% 
  t() %>% data.frame() 
names(myATE)=c("LL95","PEst","UL95")
myATE$estimand = "ATE"
myATE$model = "Propensity score subclassification (8 groups)"
myATE %>% as_tibble() 

class_estimands = bind_rows(myATT,myATC,myATE) %>% 
  as_tibble() 

# 성향점수층화 기법 추정 효과추정치 저장 
saveRDS(class_estimands,"class_estimands.RData")

# 앞서 추정한 효과추정치들 불러오기 
OLS_estimands=readRDS("OLS_estimands.RData" )
PSW_estimands=readRDS("PSW_estimands.RData")
greedy_estimands=readRDS("greedy_estimands.RData")
optimal_estimands=readRDS("optimal_estimands.RData")
full_estimands=readRDS("full_estimands.RData")
genetic_estimands=readRDS("genetic_estimands.RData")
mahala_estimands=readRDS("mahala_estimands.RData")

# 효과추정치 시각화 
bind_rows(PSW_estimands,greedy_estimands,
          optimal_estimands,full_estimands,
          genetic_estimands,mahala_estimands,
          class_estimands) %>% 
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
  guides(color = guide_legend(nrow=4))
ggsave("Part2_ch7_Comparison_Estimands.png",height=14,width=17,units='cm')
