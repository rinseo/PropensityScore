###############################################################################
## 제1부 제2장 성향점수 추정 및 공통지지영역 점검 
###############################################################################
library("tidyverse")   # 시각화를 위해 ggplot2
summarize = dplyr::summarize   # dplyr의 summarize 함수를 이용

# 데이터 소환
setwd("D:/data")
mydata=read_csv("simdata.csv")
mydata 

# GLM, 로지스틱 회귀모형을 이용한 성향점수 추정 
logis_ps = glm(treat ~ V1+V2+V3, data=mydata,
               family=binomial(link="logit"))
summary(logis_ps)

## 각주 
plogis(-4);  plogis(-3);  plogis(-2)
plogis(-2) - plogis(-3)
plogis(-3) - plogis(-4)

# 성향점수 추정결과 
mydata$ps_prob = predict(logis_ps, mydata,
                         type="response") #확률 
mydata$ps_logit = predict(logis_ps, mydata,
                          type="link")   #로짓

# 선형로짓과 추정확률의 차이
plot_logit = mydata %>% 
  ggplot(aes(x=ps_logit))+
  geom_histogram(bins=40)+
  labs(x="선형로짓(logit, -∞ ~ +∞) 형태의 성향점수",y="빈도")+
  theme_bw()
plot_prob = mydata %>% 
  ggplot(aes(x=ps_prob))+
  geom_histogram(bins=40)+
  labs(x="추정확률(probability, 0~1) 형태의 성향점수",y="빈도")+
  theme_bw()
gridExtra::grid.arrange(plot_logit,plot_prob) #그림합치기 

# jpeg("Part1_ch2_propensity_score.jpeg",res=600,height=20,width=15,units="cm")
gridExtra::grid.arrange(plot_logit,plot_prob)
# dev.off()

# 히스토그램 이용 
plot_prob_cs = ggplot(data=mydata,aes(x=ps_prob)) + 
  geom_histogram(data=mydata %>% filter(treat==1),
                 bins=40,fill="red",alpha=0.4)+
  geom_histogram(data=mydata %>% filter(treat==0),
                 bins=40,fill="blue",alpha=0.4)+
  labs(x="추정확률(probability, 0~1) 형태의 성향점수",y="빈도")+
  theme_bw()+
  ggtitle("추정확률 이용시 공통지지영역")
plot_logit_cs = ggplot(data=mydata,aes(x=ps_logit)) + 
  geom_histogram(data=mydata %>% filter(treat==1),
                 bins=40,fill="red",alpha=0.4)+
  geom_histogram(data=mydata %>% filter(treat==0),
                 bins=40,fill="blue",alpha=0.4)+
  labs(x="선형로짓(logit, -∞ ~ +∞) 형태의 성향점수",y="빈도")+
  theme_bw()+
  ggtitle("선형로짓 이용시 공통지지영역")

# jpeg("Part1_ch2_region_common_support_hist.jpeg",res=600,height=13,width=20,units="cm")
gridExtra::grid.arrange(plot_logit_cs,plot_prob_cs,nrow=1)
# dev.off()

# 박스플롯 이용 
boxplot_prob_cs = mydata %>% 
  mutate(treat = ifelse(treat==0,"통제집단","처치집단")) %>% 
  ggplot(aes(y=ps_prob,color=treat))+
  geom_boxplot()+
  labs(x="추정확률(probability, 0~1) 형태의 성향점수",y="분포",
       color="집단구분")+
  theme_bw()+
  ggtitle("추정확률 이용시 공통지지영역")
boxplot_logit_cs = mydata %>% 
  mutate(treat = ifelse(treat==0,"통제집단","처치집단")) %>% 
  ggplot(aes(y=ps_logit,color=treat))+
  geom_boxplot()+
  labs(x="선형로짓(logit, -∞ ~ +∞) 형태의 성향점수",y="빈도",
       color="집단구분")+
  theme_bw()+
  ggtitle("선형로짓 이용시 공통지지영역")

# jpeg("Part1_ch2_region_common_support_box.jpeg",res=600,height=13,width=20,units="cm")
gridExtra::grid.arrange(boxplot_logit_cs,boxplot_prob_cs,nrow=1)
# dev.off()
