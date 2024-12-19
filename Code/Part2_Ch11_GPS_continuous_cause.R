######################################################################
## 제2부 제11장 연속형 변수형태 원인변수 대상 성향점수분석 기법 
## Hirano-Imbens generalized propensity score methods 
######################################################################
# GPS : generalized propensity score

# 원인변수가 연속형인 경우  
library("tidyverse")           # 데이터관리 및 시각화 
summarize = dplyr::summarize   # dplyr의 summarize 함수를 이용
library("Zelig")               # 비모수통계접근 

# 데이터소환 
setwd("D:/data")
simdata=read_csv("simulated_continuous_cause.csv")
summary(simdata)

######################################################################
# 2. 일반화성향점수 추정 
######################################################################
# 공변량이 주어졌을 때의 원인변수 
lmGPS = lm(treat~x1+x2,simdata)
# GPS 추정: 분모부분 
simdata$gps = dnorm(simdata$treat,
                    mean=lmGPS$fitted,        #예측된 y의 평균 
                    sd=summary(lmGPS)$sigma)  #오차항의 분산
# 분자부분 
simdata$numerator = dnorm(simdata$treat,
                    mean=mean(simdata$treat),
                    sd=sd(simdata$treat))

# 분자부분 
simdata$numerator = dnorm(simdata$treat,
                          mean=mean(simdata$treat),
                          sd=sd(simdata$treat))
# 혹은 다음 방법과도 동일
lmNaive = lm(treat~1,simdata)
simdata$numerator = dnorm(simdata$treat,
                          mean=lmGPS$fitted,        #예측된 y의 평균 
                          sd=summary(lmGPS)$sigma)  #오차항의 분산
# IPW 계산 
simdata$IPW=simdata$numerator/simdata$gps

# IPW 계산 
simdata$IPW=simdata$numerator/simdata$gps

######################################################################
# 3. 공변량 균형성 점검 
######################################################################
# 표준화 회귀계수를 얻기 위해 표준화 변환 실시 
stddata = simdata %>%
  mutate_at(
    vars(x1,x2,treat),
    function(x){(x-mean(x))/sd(x)}
  )
lm(x1~treat, stddata)$coef %>% round(4)
lm(x2~treat, stddata)$coef %>% round(4)
lm(x1~treat, stddata, weights=IPW)$coef %>% round(4)
lm(x2~treat, stddata, weights=IPW)$coef %>% round(4)
######################################################################
# 4. 처치효과 추정 
######################################################################
# IPW
set.seed(1234)
z_out_ipw=zelig(y~treat+x1+x2,data=simdata,
                model="ls",weights="IPW",
                cite=FALSE)
# 시뮬레이션 
range_treat = quantile(simdata$treat,prob=0.1*(1:9))
Table_Sim10000=data.frame()
for (i in 1:length(range_treat)){
  # 원인변수의 수준 설정 
  X=setx(z_out_ipw,treat=range_treat[i],data=mydata) 
  # 지정된 원인변수 수준에서 10000번 시뮬레이션 
  S=sim(z_out_ipw,X,num=10000)
  # 각 시뮬레이션 단계에서 얻은 기댓값 
  EV=data.frame(t(get_qi(S,"ev")))
  # 정리 
  Table_Sim10000=rbind(Table_Sim10000, EV)
}
names(Table_Sim10000)=str_c("sim",1:10000)
Table_Sim10000$treat=range_treat
# 점추정치, 95% 신뢰구간 계산 
IPW_estimate = Table_Sim10000 %>% 
  pivot_longer(cols=sim1:sim10000,names_to="sim") %>% 
  group_by(treat) %>% 
  summarize(
    LL95=quantile(value,p=0.025),
    PEst=quantile(value,p=0.5),
    UL95=quantile(value,p=0.975)
  )
IPW_estimate
# Simple OLS 
set.seed(1234)
z_out_ols=zelig(y~treat+x1+x2,data=simdata,
                model="ls",cite=FALSE)
# 시뮬레이션
Table_Sim10000=data.frame()
for (i in 1:length(range_treat)){
  # 원인변수의 수준 설정 
  X=setx(z_out_ols,treat=range_treat[i],data=mydata) 
  # 지정된 원인변수 수준에서 10000번 시뮬레이션 
  S=sim(z_out_ols,X,num=10000)
  # 각 시뮬레이션 단계에서 얻은 기댓값 
  EV=data.frame(t(get_qi(S,"ev")))
  # 정리 
  Table_Sim10000=rbind(Table_Sim10000, EV)
}
names(Table_Sim10000)=str_c("sim",1:10000)
Table_Sim10000$treat=range_treat
OLS_estimate = Table_Sim10000 %>% 
  pivot_longer(cols=sim1:sim10000,names_to="sim") %>% 
  group_by(treat) %>% 
  summarize(
    LL95=quantile(value,p=0.025),
    PEst=quantile(value,p=0.5),
    UL95=quantile(value,p=0.975)
  )

# 추정결과 비교 시각화 
bind_rows(OLS_estimate %>% mutate(model="OLS"),
          IPW_estimate %>% mutate(model="IPW")) %>% 
  ggplot(aes(x=treat,y=PEst,fill=model))+
  geom_point(aes(color=model,shape=model),size=2)+
  geom_line(aes(color=model))+
  geom_ribbon(aes(ymin=LL95,ymax=UL95),alpha=0.3)+
  labs(x="X, continuous variable\n(Dose)",
       y="Point estimates with their 95% CI\n(Response)",
       fill="Model",shape="Model",color="Model")+
  scale_x_continuous(breaks=round(IPW_estimate$treat,1))+
  coord_cartesian(ylim=c(1,3))+
  theme_bw()+
  theme(legend.position="top")
ggsave("Part2_Ch11_IPW_OLS.png",unit="cm",width=11,height=11)
