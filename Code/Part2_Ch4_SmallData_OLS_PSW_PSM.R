####################################################################
# 제2부 제4장 수계산을 통한 매칭기법 이해와 적용
####################################################################
library("tidyverse")  #데이터관리 및 사전처리, 시각화를 위해 
summarize = dplyr::summarize   # dplyr의 summarize 함수를 이용
# 예시데이터 호출
setwd("D:/data")
myd = read_csv("SmallData.csv")  #데이터 호출 
myd   # 데이터 확인
lm(y~t+x,myd)    # 전통적 OLS회귀모형 

# 처치역확률 가중치(IPTW) 계산
# 성향점수 계산 
glm(t~x,myd,family=binomial(link='logit')) %>% fitted()

# 성향점수를 변수로 저장한 후
myd$ps = glm(t~x,myd,family=binomial(link='logit')) %>% fitted() 
# 처치 역확률 가중치(IPTW, inverse probability of treatment weights) 계산
myd = myd %>% 
  mutate(
    iptw=(t/ps)+((1-t)/(1-ps))  #ATE 계산용 가중치
  )
# 데이터 확인 
myd %>% select(id,x,t,y,ps,iptw)

# IPTW 부여후 회귀분석 
lm(y~t+x,myd,weights=iptw)   # ATE

####################################################################
# 성향점수 기반 그리디 매칭 
####################################################################
library('MatchIt')  #매칭기법을 위해 
myd = myd %>% select(id,x,t,y)  #필요변수만 선택 
# ATT 계산을 위한 성향점수 기반 그리디 매칭 
set.seed(12345)  #본서의 매칭결과를 정확하게 반복하기 위해
# 1단계: matchit() 함수 
m_out=matchit(t~x,  # 원인처치 선택 확률추정 함수
              myd,  # 데이터 
              distance='linear.logit',  # 성향점수(선형로짓) 
              method='nearest',  # 그리디 매칭, 최인접사례 매칭
              replace=TRUE)      # 반복추출 허용 
summary(m_out)   # 매칭결과 점검  
# 2단계: match.data() 함수 
dm_out=match.data(m_out)
dm_out

# 각주: ATT 계산 
lm(y~t+x,dm_out,weights=weights) %>% summary() 

# 3단계: Zelig 패키지를 이용하여 모형추정치(여기서는 ATT) 추정 
library("Zelig")
# 1) zelig() 함수: 효과추정치 추정모형 설정 
z_psm_att=zelig(y~t+x,data=dm_out,
                model='ls',  #OLS 모형을 사용한다는 의미
                weights="weights", #가중치를 부여한다는 의미
                cite=FALSE)  #젤리그 인용방식에 대한 안내문이 보이지 않도록 
z_psm_att
# 2) setx() 함수: 추정모형을 적용할 독립변수의 조건상정 
x_psm_att0=setx(z_psm_att,t=0,data=dm_out) #통제집단가정 상황
x_psm_att1=setx(z_psm_att,t=1,data=dm_out) #처치집단가정 상황
# 3단계: 1단계와 2단계를 근거로 기댓값(expected value, ev) 시뮬레이션
# 디폴트는 1000회, 만약 5000번의 시뮬레이션을 원하는 경우 num=5000 옵션추가 
s_psm_att0=sim(z_psm_att,x_psm_att0)  #통제집단의 경우 기댓값
s_psm_att1=sim(z_psm_att,x_psm_att1)  #처치집단의 경우 기댓값
# ATT의 값을 추정한 후 95%신뢰구간 계산 
ev_att0=get_qi(s_psm_att0,"ev")  #1000번의 통제집단 기댓값들
ev_att1=get_qi(s_psm_att1,"ev")  #1000번의 처치집단 기댓값들
PSM_ATT=ev_att1-ev_att0          #1000번 추정된 ATT
quantile(PSM_ATT,p=c(0.025,0.50,0.975)) %>% round(3) #ATT의 95%신뢰구간 

# ATC 계산을 위한 성향점수 기반 그리디 매칭 
set.seed(54321)  #본서의 매칭결과를 정확하게 반복하기 위해
# 1단계: matchit() 함수 
# 처치집단과 통제집단을 역코딩한 새로운 원인변수 생성 
myd$tr = ifelse(myd$t==1,0,1) 
m_out_r=matchit(tr~x,  # 원인처치 선택 확률추정 함수
                myd,  # 데이터 
                distance='linear.logit',  # 성향점수 
                method='nearest',  # 그리디  매칭
                replace=TRUE)      # 반복추출 허용 
summary(m_out_r)   # 매칭결과 점검  
# 2단계: match.data() 함수 
dm_out_r=match.data(m_out_r)
dm_out_r
# 3단계: Zelig 패키지를 이용하여 모형추정치(여기서는 ATC) 추정 
# 1) zelig() 함수: 효과추정치 추정모형 설정 
z_psm_atc=zelig(y~tr+x,data=dm_out_r,
                model='ls',  #OLS 모형을 사용한다는 의미
                weights="weights", #가중치를 부여한다는 의미
                cite=FALSE)  #젤리그 인용방식에 대한 안내문이 보이지 않도록 
z_psm_atc
# 2) setx() 함수: 추정모형을 적용할 독립변수의 조건상정 
x_psm_atc0=setx(z_psm_atc,tr=1,data=dm_out_r) #통제집단가정 상황
x_psm_atc1=setx(z_psm_atc,tr=0,data=dm_out_r) #처치집단가정 상황
# 3단계: 1단계와 2단계를 근거로 기댓값(expected value, ev) 시뮬레이션
s_psm_atc0=sim(z_psm_atc,x_psm_atc0)  #통제집단의 경우 기댓값
s_psm_atc1=sim(z_psm_atc,x_psm_atc1)  #처치집단의 경우 기댓값
# ATC의 값을 추정한 후 95%신뢰구간 계산 
ev_atc0=get_qi(s_psm_atc0,"ev")  #1000번의 통제집단 기댓값들
ev_atc1=get_qi(s_psm_atc1,"ev")  #1000번의 처치집단 기댓값들
PSM_ATC=ev_atc1-ev_atc0          #1000번 추정된 ATC
quantile(PSM_ATC,p=c(0.025,0.50,0.975)) %>% round(3) #ATC의 95%신뢰구간 

# ATE 계산
mypi=as.vector(prop.table(table(myd$t)))[2]
mypi
# 공식에 따라 ATE 계산 
PSM_ATE=mypi*PSM_ATT+(1-mypi)*PSM_ATC
quantile(PSM_ATE,p=c(0.025,0.50,0.975)) %>% round(3) #ATE의 95%신뢰구간 
