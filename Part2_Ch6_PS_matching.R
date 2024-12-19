# 매칭기법별 비교 
# 시율레이션 데이터 생성은 Matching_Simulated_Data.r 참조. 시드넘버=1234
###############################################################################
## 제2부 제6장 시뮬레이션 데이터 분석실습 및 매칭기법 비교
###############################################################################
library("tidyverse")   #데이터관리 및 변수사전처리 
summarize = dplyr::summarize   # dplyr의 summarize 함수를 이용
library("Zelig")   #비모수접근 95% CI 계산 
library("MatchIt")   #매칭기법 적용 
library("cobalt")   #처치집단과 통제집단 균형성 점검 
library("sensitivitymw")   #민감도 테스트 실시
library("sensitivityfull")   #민감도 테스트 실시(전체매칭용)

# 데이터 소환
setwd("D:/data")
mydata=read_csv("simdata.csv")
mydata 

###############################################################################
## 2. 성향점수매칭 기법 
###############################################################################
# 1)ATT-성향점수기반 그리디 매칭(Greedy matching using propensity score)
set.seed(1234)  #정확하게 동일한 결과를 원한다면 
greedy_att=matchit(formula=treat~V1+V2+V3,
              data=mydata,  
              distance="linear.logit", #선형로짓 형태 성향점수 
              method="nearest",  # 그리디 매칭 알고리즘
              caliper=0.15,  #성향점수의 0.15표준편차 허용범위 
              discard='none',  #공통지지영역에서 벗어난 사례도 포함 
              ratio=2,  #처치집단 사례당 통제집단 사례는 2개 매칭 
              replace=FALSE) #동일사례 반복표집 매칭을 허용하지 않음
greedy_att
# 1)ATC-성향점수기반 그리디 매칭(Greedy matching using propensity score)
set.seed(4321)  #정확하게 동일한 결과를 원한다면 
greedy_atc=matchit(formula=Rtreat~V1+V2+V3, #통제집단을 '처치집단'으로 가정한 후 실행 
              data=mydata,  
              distance="linear.logit", #선형로짓 형태 성향점수 
              method="nearest",  # 그리디 매칭 알고리즘
              caliper=0.15,  #성향점수의 0.15표준편차 허용범위 
              discard='none',  #공통지지영역에서 벗어난 사례도 포함 
              ratio=2,  #처치집단 사례당 통제집단 사례는 2개 매칭 
              replace=TRUE) #동일사례 반복표집 매칭을 허용
greedy_atc

# 2)ATT-성향점수기반 최적 매칭(Optimal matching using propensity score)
set.seed(1234)  #정확하게 동일한 결과를 원한다면 
optimal_att=matchit(formula=treat~V1+V2+V3,
                   data=mydata,  
                   distance="linear.logit", #선형로짓 형태 성향점수 
                   method="optimal",  # 최적 매칭 알고리즘
                   caliper=0.15,  #성향점수의 0.15표준편차 허용범위 
                   discard='none', #공통지지영역에서 벗어난 사례 보존 
                   ratio=2,  #처치집단 사례당 통제집단 사례는 2개 매칭 
                   replace=FALSE) #동일사례 반복표집 매칭을 허용않음 
optimal_att
# 2)ATC-성향점수기반 최적 매칭(Optimal matching using propensity score)
# 추정불가 

# 3)ATT-성향점수기반 전체 매칭(Full matching using propensity score)
set.seed(1234)  #정확하게 동일한 결과를 원한다면 
full_att=matchit(formula=treat~V1+V2+V3,
                 data=mydata,
                 distance="linear.logit", #선형로짓 형태 성향점수 
                 method="full",  # 전체 매칭 알고리즘
                 discard="none")  #공통지지영역에서 벗어난 사례 보존 
full_att
# 3)ATC-성향점수기반 전체 매칭(Full matching using propensity score)
set.seed(4321)  #정확하게 동일한 결과를 원한다면 
full_atc=matchit(formula=Rtreat~V1+V2+V3, #통제집단을 '처치집단'으로 가정한 후 실행 
                 data=mydata,
                 distance="linear.logit", #선형로짓 형태 성향점수 
                 method="full",  # 전체 매칭 알고리즘
                 discard="none")  #공통지지영역에서 벗어난 사례 보존 
full_atc

# 아래의 과정은 Google Cloud Platform으로 진행하였음
# 추정결과는 genetic_att.RData; genetic_atc.RData로 저장되어 있음. 
# 4)ATT-성향점수기반 유전 매칭(Genetic matching using propensity score)
# install.packages("rgenoud");library("rgenoud") # 패키지 오류가 발생한다면
set.seed(1234)  #정확하게 동일한 결과를 원한다면 
genetic_att=matchit(formula=treat~V1+V2+V3,data=mydata,
                    method="genetic",
                    distance="linear.logit",
                    pop.size=1000, #디폴트는 100인데, 안정적 추정을 위해 늘림
                    discard='none',
                    ratio=2,
                    distance.tolerance=1e-05, # 거리차이 기준값 디폴트 
                    ties=TRUE)  #복수의 통제집단 사례들 매칭(디폴트)
genetic_att
# 4)ATC-성향점수기반 유전 매칭(Genetic matching using propensity score)
set.seed(4321)  #정확하게 동일한 결과를 원한다면 
genetic_atc=matchit(formula=Rtreat~V1+V2+V3,data=mydata,
                    method="genetic",
                    distance="linear.logit",
                    pop.size=1000, #디폴트는 100인데, 안정적 추정을 위해 늘림
                    discard='none',
                    ratio=2,
                    distance.tolerance=1e-05, # 거리차이 기준값 디폴트 
                    ties=TRUE)  #복수의 통제집단 사례들 매칭(디폴트)
genetic_atc
# genetic_att=readRDS("genetic_att.RData");genetic_atc=readRDS("genetic_atc.RData")

# 5)ATT-마할라노비스 거리점수기반 그리디 매칭(Greedy matching using Mahalanobis distance)
set.seed(1234)  #정확하게 동일한 결과를 원한다면 
mahala_att=matchit(formula=treat~V1+V2+V3,
                   data=mydata, 
                   distance="mahalanobis", #마할라노비스 거리점수 
                   method="nearest",  # 그리디 매칭 알고리즘
                   caliper=0.15,  #성향점수의 0.15표준편차 허용범위 
                   discard='none',  #공통지지영역에서 벗어난 사례도 포함 
                   ratio=2,  #처치집단 사례당 통제집단 사례는 2개 매칭 
                   replace=FALSE) #동일사례 반복표집 매칭을 허용하지 않음
mahala_att
# 5)ATC-마할라노비스 거리점수기반 그리디 매칭(Greedy matching using Mahalanobis distance)
set.seed(4321)  #정확하게 동일한 결과를 원한다면 
mahala_atc=matchit(formula=Rtreat~V1+V2+V3, #통제집단을 '처치집단'으로 가정한 후 실행 
                   data=mydata,  
                   distance="mahalanobis", #마할라노비스 거리점수 
                   method="nearest",  # 그리디 매칭 알고리즘
                   caliper=0.15,  #성향점수의 0.15표준편차 허용범위 
                   discard='none',  #공통지지영역에서 벗어난 사례도 포함 
                   ratio=2,  #처치집단 사례당 통제집단 사례는 2개 매칭 
                   replace=TRUE) #동일사례 반복표집 매칭을 허용
mahala_atc

###############################################################################
## 3. 공변량 균형성 점검 
###############################################################################
# MatchIt 패키지에서는 summary(), plot() 함수를 이용하여 공변량의 균형성 점검 
summary(greedy_att,standardize=TRUE)  #공변량을 표준화시킨 후 균형성 점검 
# 균형성 점검을 위한 시각화 
plot(greedy_att, type="QQ")  #QQ플롯, 디폴트
plot(greedy_att, type="hist")  #(가중치가 적용된)히스토그램
plot(greedy_att, type="jitter")  #지터링 포함된 산점도

# 1)ATT-성향점수기반 그리디 매칭 균형성 점검
# 코발트 패키지
bal.tab(greedy_att,
        continuous="std",  #연속형 형태 변수들을 표준화
        s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
        m.threshold=0.1, #처치집단과 통제집단 공변량 변수의 평균차이(기준:<0.1)
        v.threshold=2)  #처치집단과 통제집단 공변량 변수의 분산비율(기준:<2)
# 시각화: 히스토그램 
bal.plot(greedy_att,
         var.name="distance", #거리점수에 대해 
         which="both",  #처치집단과 통제집단 모두 
         mirror=TRUE,  #두 집단이 서로 마주닿는 형태로 
         type="histogram")
ggsave("Part2_ch6_greedy_att_hist1.png",height=12,width=24,units='cm')
# ggplot2 패키지의 부속함수들을 추가할 수 있음. 
bal.plot(greedy_att,
         var.name="distance", #거리점수에 대해 
         which="both",  #처치집단과 통제집단 모두 
         mirror=TRUE,  #두 집단이 서로 마주닿는 형태로 
         type="histogram", 
         colors=c("grey80","grey20"))+ #막대색을 바꿈
  labs(x="Propensity score, as linear logit",
       y="Proportion",fill="Treatment")+  #라벨 변경 
  ggtitle("Change in distribution of propensity score before/after matching")+
  theme_bw()  #배경을 흑백으로 
ggsave("Part2_ch6_greedy_att_hist2.png",height=12,width=24,units='cm')

# 거리점수와 세 공변량에 대한 조정전후의 히스토그램 비교 
myf1=bal.plot(greedy_att,"V1",which="both",
              mirror=TRUE,type="histogram",
              colors=c("grey80","grey20"))+theme_bw()
myf2=bal.plot(greedy_att,"V2",which="both",
              mirror=TRUE,type="histogram",
              colors=c("grey80","grey20"))+theme_bw()
myf3=bal.plot(greedy_att,"V3",which="both",
              mirror=TRUE,type="histogram",
              colors=c("grey80","grey20"))+theme_bw()
png("Part2_ch6_greedy_att_V123_hist.png",res=600,height=10,width=30,units="cm")
gridExtra::grid.arrange(myf1,myf2,myf3,nrow=1)
dev.off()

# 러브플롯: Thomas E. Love 박사의 성을 붙인 그래프 
# 1)ATT-성향점수기반 그리디 매칭, 평균차이 
love.plot(greedy_att,
          s.d.denom="pooled", #분산은 처치집단과 통제집단 모두에서 
          stat="mean.diffs", # 두 집단간 공변량 평균차이 
          drop.distance=FALSE, #성향점수도 포함하여 제시 
          threshold=0.1, # 평균차이 역치 
          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
          themes=theme_bw())+
  coord_cartesian(xlim=c(-0.15,1.00)) #X축의 범위 
ggsave("Part2_ch6_love_MD_greedy_att.png",height=7,width=13,units='cm')
# 1)ATT-성향점수기반 그리디 매칭, 분산비 
love.plot(greedy_att,
          s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
          stat="variance.ratios", # 두 집단간 분산비 
          drop.distance=FALSE, #성향점수도 포함하여 제시 
          threshold=2, # 분산비 역치 
          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
          themes=theme_bw())+
  coord_cartesian(xlim=c(0.3,3)) #X축의 범위 
ggsave("Part2_ch6_love_VR_greedy_att.png",height=7,width=13,units='cm')

png("Part2_ch6_greedy_att_V123_hist.png",res=600,height=10,width=30,units="cm")
gridExtra::grid.arrange(myf1,myf2,myf3,nrow=1)
dev.off()

# 1)ATC-성향점수기반 그리디 매칭, 평균차이 
love_D_greedy_atc=love.plot(greedy_atc,
                            s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                            stat="mean.diffs", # 두 집단간 공변량 평균차이 
                            drop.distance=FALSE, #성향점수도 포함하여 제시 
                            threshold=0.1, # 평균차이 역치 
                            sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                            themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00)) #X축의 범위 
# 1)ATC-성향점수기반 그리디 매칭, 분산비 
love_VR_greedy_atc=love.plot(greedy_atc,
                             s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                             stat="variance.ratios", # 두 집단간 분산비 
                             drop.distance=FALSE, #성향점수도 포함하여 제시 
                             threshold=2, # 분산비 역치 
                             sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                             themes=theme_bw())+
  coord_cartesian(xlim=c(0.3,3)) #X축의 범위 
png("Part2_ch6_love_MD_VR_greedy_atc.png",res=600,height=7,width=25,units="cm")
gridExtra::grid.arrange(love_D_greedy_atc,love_VR_greedy_atc,nrow=1)
dev.off()


# 2)ATT-성향점수기반 최적 매칭, 평균차이 
love_D_optimal_att=love.plot(optimal_att,
                            s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                            stat="mean.diffs", # 두 집단간 공변량 평균차이 
                            drop.distance=FALSE, #성향점수도 포함하여 제시 
                            threshold=0.1, # 평균차이 역치 
                            sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                            themes=theme_bw())+
  coord_cartesian(xlim=c(-0.15,1.00)) #X축의 범위 
# 2)ATT-성향점수기반 최적 매칭, 분산비 
love_VR_optimal_att=love.plot(optimal_att,
                             s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                             stat="variance.ratios", # 두 집단간 분산비 
                             drop.distance=FALSE, #성향점수도 포함하여 제시 
                             threshold=2, # 분산비 역치 
                             sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                             themes=theme_bw())+
  coord_cartesian(xlim=c(0.3,3)) #X축의 범위 
png("Part2_ch6_love_MD_VR_optimal_att.png",res=600,height=7,width=25,units="cm")
gridExtra::grid.arrange(love_D_optimal_att,love_VR_optimal_att,nrow=1)
dev.off()


# 3)ATT-성향점수기반 전체 매칭, 평균차이 
love_D_full_att=love.plot(full_att,
                          s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                          stat="mean.diffs", # 두 집단간 공변량 평균차이 
                          drop.distance=FALSE, #성향점수도 포함하여 제시 
                          threshold=0.1, # 평균차이 역치 
                          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                          themes=theme_bw())+
  coord_cartesian(xlim=c(-0.15,1.00)) #X축의 범위 
# 3)ATT-성향점수기반 전체 매칭, 분산비 
love_VR_full_att=love.plot(full_att,
                           s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                           stat="variance.ratios", # 두 집단간 분산비 
                           drop.distance=FALSE, #성향점수도 포함하여 제시 
                           threshold=2, # 분산비 역치 
                           sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                           themes=theme_bw())+
  coord_cartesian(xlim=c(0.3,3)) #X축의 범위 
png("Part2_ch6_love_MD_VR_full_att.png",res=600,height=7,width=25,units="cm")
gridExtra::grid.arrange(love_D_full_att,love_VR_full_att,nrow=1)
dev.off()

# 3)ATC-성향점수기반 전체 매칭, 평균차이 
love_D_full_atc=love.plot(full_atc,
                          s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                          stat="mean.diffs", # 두 집단간 공변량 평균차이 
                          drop.distance=FALSE, #성향점수도 포함하여 제시 
                          threshold=0.1, # 평균차이 역치 
                          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                          themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00)) #X축의 범위 
# 3)ATC-성향점수기반 전체 매칭, 분산비 
love_VR_full_atc=love.plot(full_atc,
                           s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                           stat="variance.ratios", # 두 집단간 분산비 
                           drop.distance=FALSE, #성향점수도 포함하여 제시 
                           threshold=2, # 분산비 역치 
                           sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                           themes=theme_bw())+
  coord_cartesian(xlim=c(0.3,3)) #X축의 범위 
png("Part2_ch6_love_MD_VR_full_atc.png",res=600,height=7,width=25,units="cm")
gridExtra::grid.arrange(love_D_full_atc,love_VR_full_atc,nrow=1)
dev.off()


# 4)ATT-성향점수기반 유전 매칭, 평균차이 
love_D_genetic_att=love.plot(genetic_att,
                          s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                          stat="mean.diffs", # 두 집단간 공변량 평균차이 
                          drop.distance=FALSE, #성향점수도 포함하여 제시 
                          threshold=0.1, # 평균차이 역치 
                          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                          themes=theme_bw())+
  coord_cartesian(xlim=c(-0.15,1.00)) #X축의 범위 
# 4)ATT-성향점수기반 유전 매칭, 분산비 
love_VR_genetic_att=love.plot(genetic_att,
                           s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                           stat="variance.ratios", # 두 집단간 분산비 
                           drop.distance=FALSE, #성향점수도 포함하여 제시 
                           threshold=2, # 분산비 역치 
                           sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                           themes=theme_bw())+
  coord_cartesian(xlim=c(0.3,3)) #X축의 범위 
png("Part2_ch6_love_MD_VR_genetic_att.png",res=600,height=7,width=25,units="cm")
gridExtra::grid.arrange(love_D_genetic_att,love_VR_genetic_att,nrow=1)
dev.off()

# 4)ATC-성향점수기반 유전 매칭, 평균차이 
love_D_genetic_atc=love.plot(genetic_atc,
                          s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                          stat="mean.diffs", # 두 집단간 공변량 평균차이 
                          drop.distance=FALSE, #성향점수도 포함하여 제시 
                          threshold=0.1, # 평균차이 역치 
                          sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                          themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00)) #X축의 범위 
# 4)ATC-성향점수기반 유전 매칭, 분산비 
love_VR_genetic_atc=love.plot(genetic_atc,
                           s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                           stat="variance.ratios", # 두 집단간 분산비 
                           drop.distance=FALSE, #성향점수도 포함하여 제시 
                           threshold=2, # 분산비 역치 
                           sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                           themes=theme_bw())+
  coord_cartesian(xlim=c(0.3,3)) #X축의 범위 
png("Part2_ch6_love_MD_VR_genetic_atc.png",res=600,height=7,width=25,units="cm")
gridExtra::grid.arrange(love_D_genetic_atc,love_VR_genetic_atc,nrow=1)
dev.off()


# 5)ATT-마할라노비스 거리점수 그리디 매칭, 평균차이 
love_D_mahala_att=love.plot(mahala_att,
                             s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                             stat="mean.diffs", # 두 집단간 공변량 평균차이 
                             drop.distance=FALSE, #성향점수도 포함하여 제시 
                             threshold=0.1, # 평균차이 역치 
                             sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                             themes=theme_bw())+
  coord_cartesian(xlim=c(-0.15,1.00)) #X축의 범위 
# 5)ATT-마할라노비스 거리점수 그리디 매칭, 분산비 
love_VR_mahala_att=love.plot(mahala_att,
                              s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                              stat="variance.ratios", # 두 집단간 분산비 
                              drop.distance=FALSE, #성향점수도 포함하여 제시 
                              threshold=2, # 분산비 역치 
                              sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                              themes=theme_bw())+
  coord_cartesian(xlim=c(0.3,3)) #X축의 범위 
png("Part2_ch6_love_MD_VR_mahala_att.png",res=600,height=7,width=25,units="cm")
gridExtra::grid.arrange(love_D_mahala_att,love_VR_mahala_att,nrow=1)
dev.off()

# 5)ATC-마할라노비스 거리점수 그리디 매칭, 평균차이 
love_D_mahala_atc=love.plot(mahala_atc,
                             s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                             stat="mean.diffs", # 두 집단간 공변량 평균차이 
                             drop.distance=FALSE, #성향점수도 포함하여 제시 
                             threshold=0.1, # 평균차이 역치 
                             sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                             themes=theme_bw())+
  coord_cartesian(xlim=c(-0.5,1.00)) #X축의 범위 
# 5)ATC-마할라노비스 거리점수 그리디 매칭, 분산비 
love_VR_mahala_atc=love.plot(mahala_atc,
                              s.d.denom="pooled", # 분산은 처치집단과 통제집단 모두에서 
                              stat="variance.ratios", # 두 집단간 분산비 
                              drop.distance=FALSE, #성향점수도 포함하여 제시 
                              threshold=2, # 분산비 역치 
                              sample.names=c("Unmatched", "Matched"), #처치집단/통제집단 표시 
                              themes=theme_bw())+
  coord_cartesian(xlim=c(0.3,3)) #X축의 범위 
png("Part2_ch6_love_MD_VR_mahala_atc.png",res=600,height=7,width=25,units="cm")
gridExtra::grid.arrange(love_D_mahala_atc,love_VR_mahala_atc,nrow=1)
dev.off()

###############################################################################
## 4. 처치효과 추정 
###############################################################################
# 1)ATT-성향점수기반 그리디 매칭
# 1단계: 매칭 데이터 생성 
MD_greedy_att=match.data(greedy_att)
dim(MD_greedy_att)  # 사례수 x 변수의 수  
head(MD_greedy_att) %>% round(2)  #데이터 형태 
# 2단계: 모형추정 
set.seed(1234) #동일한 결과를 얻기 위해 
z_model=zelig(formula=y~treat+V1+V2+V3,
              data=MD_greedy_att,
              model='ls',
              weights="weights",cite=FALSE)
# 3단계: 추정된 모형을 적용할 X변수의 조건상정 
x_0=setx(z_model,treat=0,data=MD_greedy_att) #통제집단가정 상황
x_1=setx(z_model,treat=1,data=MD_greedy_att) #처치집단가정 상황
# 4단계: 1단계와 2단계를 근거로 기댓값(expected value, ev) 시뮬레이션
s_0=sim(z_model,x_0,num=10000)
s_1=sim(z_model,x_1,num=10000)
# 5단계: ATT의 값을 추정한 후 95%신뢰구간 계산 
EST1=get_qi(s_1,"ev") - get_qi(s_0,"ev")
summary_est1=tibble(
  LL95=quantile(EST1,p=c(0.025)),
  PEst=quantile(EST1,p=c(0.500)),
  UL95=quantile(EST1,p=c(0.975)),
  estimand="ATT",model="Greedy matching using propensity score"
)
summary_est1

# ATT추정을 위한 이용자정의 함수 
SUMMARY_EST_ATT=function(myformula,matched_data,n_sim,model_name){
  # 2단계: 모형추정 
  z_model=zelig(as.formula(myformula),
                data=matched_data,
                model='ls',
                weights="weights",cite=FALSE)
  # 3단계: 추정된 모형을 적용할 X변수의 조건상정 
  x_0=setx(z_model,treat=0,data=matched_data) #통제집단가정 상황
  x_1=setx(z_model,treat=1,data=matched_data) #처치집단가정 상황
  # 4단계: 1단계와 2단계를 근거로 기댓값(expected value, ev) 시뮬레이션
  s_0=sim(z_model,x_0,num=n_sim)
  s_1=sim(z_model,x_1,num=n_sim)
  # 5단계: ATT의 값을 추정한 후 95%신뢰구간 계산 
  EST1=get_qi(s_1,"ev") - get_qi(s_0,"ev")
  summary_est1=tibble(
    LL95=quantile(EST1,p=c(0.025)),
    PEst=quantile(EST1,p=c(0.500)),
    UL95=quantile(EST1,p=c(0.975)),
    estimand="ATT",model=model_name
  )
  rm(z_model,x_0,x_1,s_0,s_1)
  list(EST1,summary_est1)
}
set.seed(1234)
greedy_ATT = SUMMARY_EST_ATT(myformula="y~treat+V1+V2+V3",
                        matched_data=MD_greedy_att,
                        n_sim=10000,
                        model_name="Greedy matching using propensity score")
length(greedy_ATT[[1]])  # 시뮬레이션 결과(ATE 계산시 활용)
greedy_ATT[[2]]  #95% 신뢰구간 

# 1)ATC-성향점수기반 그리디 매칭
# 1단계: 매칭 데이터 생성 
MD_greedy_atc=match.data(greedy_atc)
# ATC추정을 위한 이용자정의 함수 
SUMMARY_EST_ATC=function(myformula,matched_data,n_sim,model_name){
  # 2단계: 모형추정 
  z_model=zelig(as.formula(myformula),
                data=matched_data,
                model='ls',
                weights="weights",cite=FALSE)
  # 3단계: 추정된 모형을 적용할 X변수의 조건상정 
  x_0=setx(z_model,Rtreat=1,data=matched_data) #통제집단가정 상황(이 부분 주의)
  x_1=setx(z_model,Rtreat=0,data=matched_data) #처치집단가정 상황(이 부분 주의)
  # 4단계: 1단계와 2단계를 근거로 기댓값(expected value, ev) 시뮬레이션
  s_0=sim(z_model,x_0,num=n_sim)
  s_1=sim(z_model,x_1,num=n_sim)
  # 5단계: ATT의 값을 추정한 후 95%신뢰구간 계산 
  EST1=get_qi(s_1,"ev") - get_qi(s_0,"ev")
  summary_est1=tibble(
    LL95=quantile(EST1,p=c(0.025)),
    PEst=quantile(EST1,p=c(0.500)),
    UL95=quantile(EST1,p=c(0.975)),
    estimand="ATC",model=model_name
  )
  rm(z_model,x_0,x_1,s_0,s_1)
  list(EST1,summary_est1)
}
set.seed(4321)
greedy_ATC = SUMMARY_EST_ATC(myformula="y~Rtreat+V1+V2+V3",
                        matched_data=MD_greedy_atc,
                        n_sim=10000,
                        model_name="Greedy matching using propensity score")
greedy_ATC[[2]]  #95% 신뢰구간 

# ATE-성향점수기반 그리디 매칭
# 표본내 처치집단 비율 
mypi = prop.table(table(mydata$treat))[2]
# 6단계: ATE 점추정치와 95% CI 
greedy_ATE=list()
greedy_ATE[[1]] = mypi*greedy_ATT[[1]] + (1-mypi)*greedy_ATC[[1]] #10000개의 ATE 
greedy_ATE[[2]] = tibble(
  LL95=quantile(greedy_ATE[[1]],p=c(0.025)),
  PEst=quantile(greedy_ATE[[1]],p=c(0.500)),
  UL95=quantile(greedy_ATE[[1]],p=c(0.975)),
  estimand="ATE",model="Greedy matching using propensity score")
greedy_ATE[[2]]  

# 효과추정치 저장 
greedy_estimands=bind_rows(greedy_ATT[[2]],
                           greedy_ATC[[2]],
                           greedy_ATE[[2]])
greedy_estimands


# 2)ATT-성향점수기반 최적 매칭
# 1단계: 매칭 데이터 생성 
MD_optimal_att=match.data(optimal_att)
# 2-5단계: ATT 추정 
set.seed(1234)
optimal_ATT = SUMMARY_EST_ATT(myformula="y~treat+V1+V2+V3",
                             matched_data=MD_optimal_att,
                             n_sim=10000,
                             model_name="Optimal matching using propensity score")
optimal_ATT[[2]]  #95% 신뢰구간 
# 효과추정치 저장 
optimal_estimands=optimal_ATT[[2]]


# 1)ATT-성향점수기반 전체 매칭
# 1단계: 매칭 데이터 생성 
MD_full_att=match.data(full_att)
# 2-5단계: ATT 추정 
set.seed(1234)
full_ATT = SUMMARY_EST_ATT(myformula="y~treat+V1+V2+V3",
                             matched_data=MD_full_att,
                             n_sim=10000,
                             model_name="Full matching using propensity score")
full_ATT[[2]]  #95% 신뢰구간 

# 1)ATC-성향점수기반 전체 매칭
# 1단계: 매칭 데이터 생성 
MD_full_atc=match.data(full_atc)
# 2-5단계: ATC 추정 
set.seed(4321)
full_ATC = SUMMARY_EST_ATC(myformula="y~Rtreat+V1+V2+V3",
                             matched_data=MD_full_atc,
                             n_sim=10000,
                             model_name="Full matching using propensity score")
full_ATC[[2]]  #95% 신뢰구간 

# 6단계: ATE 점추정치와 95% CI 
full_ATE=list()
full_ATE[[1]] = mypi*full_ATT[[1]] + (1-mypi)*full_ATC[[1]]
full_ATE[[2]] = tibble(
  LL95=quantile(full_ATE[[1]],p=c(0.025)),
  PEst=quantile(full_ATE[[1]],p=c(0.500)),
  UL95=quantile(full_ATE[[1]],p=c(0.975)),
  estimand="ATE",model="Full matching using propensity score")
full_ATE[[2]]

# 효과추정치 저장
full_estimands=bind_rows(full_ATT[[2]],
                           full_ATC[[2]],
                           full_ATE[[2]])
full_estimands



# 1)ATT-성향점수기반 유전 매칭
# 1단계: 매칭 데이터 생성 
MD_genetic_att=match.data(genetic_att)
# 2-5단계: ATT 추정 
set.seed(1234)
genetic_ATT = SUMMARY_EST_ATT(myformula="y~treat+V1+V2+V3",
                           matched_data=MD_genetic_att,
                           n_sim=10000,
                           model_name="Genetic matching using propensity score")
genetic_ATT[[2]]  #95% 신뢰구간 

# 1)ATC-성향점수기반 유전 매칭
# 1단계: 매칭 데이터 생성 
MD_genetic_atc=match.data(genetic_atc)
# 2-5단계: ATC 추정 
set.seed(4321)
genetic_ATC = SUMMARY_EST_ATC(myformula="y~Rtreat+V1+V2+V3",
                           matched_data=MD_genetic_atc,
                           n_sim=10000,
                           model_name="Genetic matching using propensity score")
genetic_ATC[[2]]  #95% 신뢰구간 

# 6단계: ATE 점추정치와 95% CI 
genetic_ATE=list()
genetic_ATE[[1]] = mypi*genetic_ATT[[1]] + (1-mypi)*genetic_ATC[[1]]
genetic_ATE[[2]] = tibble(
  LL95=quantile(genetic_ATE[[1]],p=c(0.025)),
  PEst=quantile(genetic_ATE[[1]],p=c(0.500)),
  UL95=quantile(genetic_ATE[[1]],p=c(0.975)),
  estimand="ATE",model="Genetic matching using propensity score")
genetic_ATE[[2]]

# 효과추정치 저장
genetic_estimands=bind_rows(genetic_ATT[[2]],
                            genetic_ATC[[2]],
                            genetic_ATE[[2]])
genetic_estimands


# 1)ATT-마할라노비스 거리점수기반 그리디 매칭
# 1단계: 매칭 데이터 생성
MD_mahala_att=match.data(mahala_att)
# 2-5단계: ATT 추정 
set.seed(1234)
mahala_ATT = SUMMARY_EST_ATT(myformula="y~treat+V1+V2+V3",
                             matched_data=MD_mahala_att,
                             n_sim=10000,
                             model_name="Greedy matching using Mahalanobis distance")
mahala_ATT[[2]]  #95% 신뢰구간 

# 1)ATC-마할라노비스 거리점수기반 그리디 매칭
# 1단계: 매칭 데이터 생성
MD_mahala_atc=match.data(mahala_atc)
# 2-5단계: ATC 추정 
set.seed(4321)
mahala_ATC = SUMMARY_EST_ATC(myformula="y~Rtreat+V1+V2+V3",
                             matched_data=MD_mahala_atc,
                             n_sim=10000,
                             model_name="Greedy matching using Mahalanobis distance")
mahala_ATC[[2]]  #95% 신뢰구간 

# 6단계: ATE 점추정치와 95% CI 
mahala_ATE=list()
mahala_ATE[[1]] = mypi*mahala_ATT[[1]] + (1-mypi)*mahala_ATC[[1]]
# ATE 95% 신뢰구간 
mahala_ATE[[2]] = tibble(
  LL95=quantile(mahala_ATE[[1]],p=c(0.025)),
  PEst=quantile(mahala_ATE[[1]],p=c(0.500)),
  UL95=quantile(mahala_ATE[[1]],p=c(0.975)),
  estimand="ATE",model="Greedy matching using Mahalanobis distance")
mahala_ATE[[2]]

# 효과추정치 저장
mahala_estimands=bind_rows(mahala_ATT[[2]],
                           mahala_ATC[[2]],
                           mahala_ATE[[2]])
mahala_estimands


# 매칭기법별 효과 추정치 비교 
# PSW 효과추정치 불러오기(제2부 제5장)
PSW_estimands=readRDS("PSW_estimands.RData")  
PSW_estimands
# 일반적 OLS 
OLS_estimands=lm(y~treat+V1+V2+V3, mydata) %>% 
  confint("treat") %>% as_tibble()
names(OLS_estimands)=c("LL95","UL95")
OLS_estimands

# 기법별 효과추정치 비교 
bind_rows(PSW_estimands,greedy_estimands,
          optimal_estimands,full_estimands,
          genetic_estimands,mahala_estimands) %>% 
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
  guides(color = guide_legend(nrow=3))
ggsave("Part2_ch6_Comparison_Estimands.png",height=12,width=16,units='cm')

###############################################################################
## 5. 민감도 분석 
###############################################################################
# 1)ATT-성향점수기반 그리디 매칭
# 매칭된 형태의 행렬 추출 
my_match_mat = greedy_att$match.matrix
head(my_match_mat)
# ratio=2 옵션을 지정했기 때문에 통제집단 사례는 2개씩 매칭된 것 확인
my_match_mat = data.frame(my_match_mat) #데이터프레임형태로 변환 
names(my_match_mat) = str_c("controlcase",1:2)  #변수이름 생성 
head(my_match_mat)
# 처치집단 사례생성(행렬의 가로줄 이름이 바로 처치집단 사례 번호)
# 현재 매칭된 사례번호가 수치형 변수가 아니기에 이를 수치형으로 변환
my_match_mat = my_match_mat %>% 
  mutate(
    treatedcase = as.numeric(as.character(row.names(my_match_mat))),
    controlcase1 = as.numeric(as.character(controlcase1)),
    controlcase2 = as.numeric(as.character(controlcase2))
  ) 
summary(my_match_mat)
# 공통지지영역(common support region)을 벗어난 사례들 제거(-1 혹은 NA)
dropCSR = function(myvar){
  mynewvar=ifelse(myvar<0,NA,myvar)} #-1을 결측값으로 
my_match_mat = my_match_mat %>% 
  mutate_all(dropCSR) %>% 
  drop_na()
summary(my_match_mat)   # NA 제거 확인 
# 매칭된 사례의 종속변수값을 찾아 새롭게 매칭된 데이터를 만들어 봅시다
myT = mydata[my_match_mat$treatedcase,][["y"]]
myC1 =  mydata[my_match_mat$controlcase1,][["y"]]
myC2 =  mydata[my_match_mat$controlcase2,][["y"]]
SA_greedy_att = data.frame(myT,myC1,myC2)
head(SA_greedy_att)

# senmw() 함수로 민감도 테스트 실시 
senmw(SA_greedy_att, gamma=2, method="w") 

# 위의 과정을 이용자 함수로 만들어서 반복적으로 사용
MATCH_PAIR_GENERATE = function(obj_matchit, rawdata){
  my_match_mat = obj_matchit$match.matrix
  my_match_mat = data.frame(my_match_mat)[,1:2] #데이터프레임형태로 변환 
  names(my_match_mat) = str_c("controlcase",1:2)  #변수이름 생성 
  my_match_mat = my_match_mat %>% 
    mutate(
      treatedcase = row.names(my_match_mat),
      controlcase1 = as.numeric(as.character(controlcase1)),
      controlcase2 = as.numeric(as.character(controlcase2))
    ) 
  my_match_mat = my_match_mat %>% 
    mutate_all(function(myvar) {mynewvar = ifelse(myvar<0, NA, myvar)}) %>% 
    drop_na()
  myT = rawdata[my_match_mat$treatedcase,][["y"]]
  myC1 =  rawdata[my_match_mat$controlcase1,][["y"]]
  myC2 =  rawdata[my_match_mat$controlcase2,][["y"]]
  data.frame(myT,myC1,myC2)
}
SA_greedy_att = MATCH_PAIR_GENERATE(greedy_att,mydata)
senmw(SA_greedy_att, gamma=2, method="w") 

# 귀무가설을 수용하게 되는 감마 수치를 탐색  
GAMMA_RANGE_SEARCH = function(Matched_Pair_data,gamma_start){
  mygamma = gamma_start 
  pvalue_2tail = 0   # Gamma의 시작값은 1, p-value의 경우 양측검정 기준 
  myresult = data.frame() 
  while (pvalue_2tail < 0.025) { #양측검정으로 H0를 기각하면 중단 
    result_SA = data.frame(senmw(Matched_Pair_data,gamma=mygamma,method="w"))
    pvalue_2tail = result_SA$pval
    temp = data.frame(cbind(mygamma,result_SA[,c("pval")]))
    names(temp)=c("Gamma","p_value")
    myresult = rbind(myresult,temp)
    mygamma = mygamma+0.1
  }
  myresult
}

# 감마값 범위 탐색 
SA_greedy_att_gamma_range = GAMMA_RANGE_SEARCH(SA_greedy_att, 1)
tail(SA_greedy_att_gamma_range)

# 1)ATC-성향점수기반 그리디 매칭
# ATC의 경우 원인처치변수를 역코딩하였기 때문에 -1을 곱하여 줌  
SA_greedy_atc = -1*MATCH_PAIR_GENERATE(greedy_atc,mydata) 
SA_greedy_atc_gamma_range = GAMMA_RANGE_SEARCH(SA_greedy_atc, 1)
tail(SA_greedy_atc_gamma_range)

# 2)ATT-성향점수기반 최적 매칭
SA_optimal_att = MATCH_PAIR_GENERATE(optimal_att,mydata) 
SA_optimal_att_gamma_range = GAMMA_RANGE_SEARCH(SA_optimal_att, 1)
tail(SA_optimal_att_gamma_range)

# 3)ATT-성향점수기반 전체 매칭
MD_full_att=match.data(full_att)
head(MD_full_att) %>% round(2) # 데이터 형태 
MD_full_att %>% count(subclass)  # 몇개의 집단? 
MD_full_att %>%
  filter(subclass<4) %>%   #1,2,3번 집단만 선별 
  xtabs(~treat+subclass,data=.)    #교차표 구하기 

MD_full_att %>%
  filter(subclass<4) %>%   #1,2,3번 집단만 선별 
  select(treat,subclass,y) %>% 
  arrange(subclass,desc(treat)) %>% round(2) %>% 
  mutate(treat=ifelse(treat==1,"Treat","Control")) 

# 각주 
MD_full_att %>%
  group_by(subclass) %>%
  summarize(lengthT=sum(treat),lengthC=length(treat)-lengthT) %>%
  count(lengthT,lengthC) %>%
  arrange(-lengthT)

g1=MD_full_att %>% filter(subclass==1) %>% arrange(desc(treat))
round(g1$y,2)
g2=MD_full_att %>% filter(subclass==2) %>% arrange(desc(treat))
round(g2$y,2)
g3=MD_full_att %>% filter(subclass==3) %>% arrange(desc(treat))
round(g3$y,2)

# 3)ATT-성향점수기반 전체 매칭
# 민감도 분석 : senfm() 함수를 이용해야 함. 상대적으로 더 복잡 
MATCH_PAIR_GENERATE_FULL = function(full_matched_data, raw_data){
  mynewdata = raw_data 
  mynewdata$full_subclass = full_matched_data$subclass
  temp_row_count = mynewdata %>% count(full_subclass)
  temp_col_count = mynewdata %>% 
    group_by(full_subclass) %>% 
    mutate(rid=row_number()) %>% ungroup()
  ROW_NUMBER_MATRIX_FULL = dim(temp_row_count)[1]
  COL_NUMBER_MATRIX_FULL = max(temp_col_count$rid)
  mymatrix = matrix(NA,nrow=ROW_NUMBER_MATRIX_FULL,
                    ncol=COL_NUMBER_MATRIX_FULL)
  myTCstatus = rep(NA,ROW_NUMBER_MATRIX_FULL)
  for (i in 1:ROW_NUMBER_MATRIX_FULL){
    temp_subclass = mynewdata %>% filter(full_subclass==i)
    y_trt = temp_subclass$y[temp_subclass$treat==1]
    y_ctrl = temp_subclass$y[temp_subclass$treat==0]
    CtoT = c(y_ctrl,y_trt)
    TtoC = c(y_trt,y_ctrl)
    lengthC = rep(length(y_ctrl),length(CtoT))
    lengthT = rep(length(y_trt),length(TtoC))
    y_row = ifelse(lengthC>=lengthT, TtoC, CtoT)
    WhichFirst = ifelse(length(y_ctrl)>=length(y_trt),TRUE,FALSE)
    mymatrix[i,1:length(y_row)] = y_row 
    myTCstatus[i] = WhichFirst 
  }
  final_result = list(mymatrix, myTCstatus)
  final_result
}

SA_match_data_full_att = MATCH_PAIR_GENERATE_FULL(full_att, mydata)
senfm(SA_match_data_full_att[[1]],SA_match_data_full_att[[2]],gamma=1)

GAMMA_RANGE_SEARCH_FULL = function(Matched_Matrix, TC_status, gamma_start) {
  mygamma = gamma_start 
  pvalue_2tail = 0   # Gamma의 시작값은 1, p-value의 경우 양측검정 기준 
  myresult = data.frame() 
  while (pvalue_2tail < 0.025) { #양측검정으로 H0를 기각하면 중단 
    result_SA = data.frame(senfm(Matched_Matrix, 
                                 TC_status,gamma=mygamma))
    pvalue_2tail = result_SA$pval
    temp = data.frame(cbind(mygamma,result_SA[,c("pval")]))
    names(temp)=c("Gamma","p_value")
    myresult = rbind(myresult,temp)
    mygamma = mygamma+0.1
  }
  myresult
}

SA_full_att_gamma_range = GAMMA_RANGE_SEARCH_FULL(
  SA_match_data_full_att[[1]],
  SA_match_data_full_att[[2]],1)
tail(SA_full_att_gamma_range)

# 3)ATC-성향점수기반 전체 매칭
SA_match_data_full_atc = MATCH_PAIR_GENERATE_FULL(
  full_atc,mydata)
SA_full_atc_gamma_range = GAMMA_RANGE_SEARCH_FULL(
  SA_match_data_full_atc[[1]],
  SA_match_data_full_atc[[2]],1)
tail(SA_full_atc_gamma_range)


# 4)ATT-성향점수기반 유전 매칭
SA_genetic_att = MATCH_PAIR_GENERATE(genetic_att,mydata) 
SA_genetic_att_gamma_range = GAMMA_RANGE_SEARCH(SA_genetic_att, 1)
tail(SA_genetic_att_gamma_range)

# 4)ATC-성향점수기반 유전 매칭
SA_genetic_atc = -1*MATCH_PAIR_GENERATE(genetic_atc,mydata) 
SA_genetic_atc_gamma_range = GAMMA_RANGE_SEARCH(SA_genetic_atc, 1)
tail(SA_genetic_atc_gamma_range)


# 5)ATT-마할라노비스 거리점수기반 그리디 매칭
SA_mahala_att = MATCH_PAIR_GENERATE(mahala_att,mydata) 
SA_mahala_att_gamma_range = GAMMA_RANGE_SEARCH(SA_mahala_att, 1)
tail(SA_mahala_att_gamma_range)

# 5)ATC-마할라노비스 거리점수기반 그리디 매칭
SA_mahala_atc = -1*MATCH_PAIR_GENERATE(mahala_atc,mydata) 
SA_mahala_atc_gamma_range = GAMMA_RANGE_SEARCH(SA_mahala_atc, 1)
tail(SA_mahala_atc_gamma_range)

# 효과추정치와 매칭기법 덧붙이기 이용자정의 함수 
temporary_function=function(obj_gamma_range,
                           names_estimand,
                           name_model){
  obj_gamma_range$estimand=names_estimand
  obj_gamma_range$model=name_model
  obj_gamma_range
}
# 5개 매칭 기법들 민감도분석결과 통합 
G_range = bind_rows(
  temporary_function(SA_greedy_att_gamma_range,
                    "ATT","Greedy matching\nusing propensity score"),
  temporary_function(SA_greedy_atc_gamma_range,
                    "ATC","Greedy matching\nusing propensity score"),
  temporary_function(SA_optimal_att_gamma_range,
                    "ATT","Optimal matching\nusing propensity score"),
  temporary_function(SA_full_att_gamma_range,
                    "ATT","Full matching\nusing propensity score"),
  temporary_function(SA_full_atc_gamma_range,
                    "ATC","Full matching\nusing propensity score"),
  temporary_function(SA_genetic_att_gamma_range,
                    "ATT","Genetic matching\nusing propensity score"),
  temporary_function(SA_genetic_atc_gamma_range,
                    "ATC","Genetic matching\nusing propensity score"),
  temporary_function(SA_mahala_att_gamma_range,
                    "ATT","Greedy matching\nusing Mahalanobis distance"),
  temporary_function(SA_mahala_atc_gamma_range,
                    "ATC","Greedy matching\nusing Mahalanobis distance")
)

# 민감도 분석 결과 종합 비교 
G_range %>% 
  mutate(
    Gamma = ifelse(p_value>0.025,NA,Gamma), #0.025를 넘는 유의도인 경우 결측값 
    rid=row_number(),
    estimand=fct_reorder(estimand,rid),
    model=fct_rev(fct_reorder(model,rid))
  ) %>% drop_na() %>% 
  ggplot(aes(x=model,y=Gamma))+
  geom_point(size=2)+
  coord_flip(ylim=c(5,16))+
  scale_y_continuous(breaks=2*(3:9),)+
  labs(x="Matching technique\n",
       y=expression(paste(Gamma," statistic")))+
  theme_bw()+
  facet_wrap(~estimand)
ggsave("Part2_ch6_sensitivity_analysis_gamma.png",height=8,width=13,units='cm')

# 효과추정치 저장 
saveRDS(OLS_estimands,"OLS_estimands.RData")
saveRDS(greedy_estimands,"greedy_estimands.RData")
saveRDS(optimal_estimands,"optimal_estimands.RData")
saveRDS(full_estimands,"full_estimands.RData")
saveRDS(genetic_estimands,"genetic_estimands.RData")
saveRDS(mahala_estimands,"mahala_estimands.RData")

# 사용대상 함수들 별도저장 
saveRDS(SUMMARY_EST_ATC,"SUMMARY_EST_ATC.RData")
saveRDS(SUMMARY_EST_ATT,"SUMMARY_EST_ATT.RData")
saveRDS(GAMMA_RANGE_SEARCH,"GAMMA_RANGE_SEARCH.RData")
saveRDS(GAMMA_RANGE_SEARCH_FULL,"GAMMA_RANGE_SEARCH_FULL.RData")
saveRDS(MATCH_PAIR_GENERATE,"MATCH_PAIR_GENERATE.RData")
saveRDS(MATCH_PAIR_GENERATE_FULL,"MATCH_PAIR_GENERATE_FULL.RData")
