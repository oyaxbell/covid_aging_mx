## Database management####
library(readr); library(tidyverse); library(survival); library(mediation); library(ggpubr); library(rms); library(psych); library(smoothHR)
library(survminer); library(haven); library(rsq); library(ResourceSelection); library(ggsci);library(timereg); library(coxme)
library(pROC);library(sf); library(rgdal); library(ggpubr); library(ggsci); library(ggmap); library(scales); library(jtools); library(cowplot)
library(ggstance); library(flextable); library(simPH); library(ggthemes); library(lme4); library(lmerTest); library(prismatic)

setwd("C:/Users/HP-PC/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/COVID - Envejecimiento")

covid <- read_csv("200509COVID19MEXICO.csv")

covid$id<-paste0(str_pad(covid$ENTIDAD_RES, 2,pad = "0"),str_pad(covid$MUNICIPIO_RES,3, pad="0"))
covid1<-covid[,c(14:15,18:30,32:35)]
covid<-covid[,-c(14:15,18:30,32:35)]
covid1[covid1==2]<-0
covid1[covid1==97]<-NA;covid1[covid1==98]<-NA;covid1[covid1==99]<-NA
covid<-as.data.frame(cbind(covid, covid1))
covid$EMBARAZO[is.na(covid$EMBARAZO)]<-0
covid$TIPO_PACIENTE[covid$TIPO_PACIENTE==1]<-0;covid$TIPO_PACIENTE[covid$TIPO_PACIENTE==2]<-1
covid$covid<-NULL; covid$covid[covid$RESULTADO==1]<-1;covid$covid[covid$RESULTADO!=1]<-0
covid$edad60<-NULL;covid$edad60[covid$EDAD>=60]<-1;covid$edad60[covid$EDAD<60]<-0
covid$edad40<-NULL;covid$edad40[covid$EDAD>=40]<-0;covid$edad40[covid$EDAD<40]<-1
covid$diabetes_40<-NULL;covid$diabetes_40[covid$DIABETES==1 & covid$edad40==1]<-1;covid$diabetes_40[covid$DIABETES!=1 & covid$edad40!=1]<-0
covid$Mortalidad<-NULL; covid$Mortalidad[is.na(covid$FECHA_DEF)]<-0;covid$Mortalidad[is.na(covid$FECHA_DEF)==FALSE]<-1
covid$FECHA_DEF[is.na(covid$FECHA_DEF)]<-as.Date(Sys.Date())
covid$FU_time<-as.numeric(as.Date(covid$FECHA_DEF)-as.Date(covid$FECHA_SINTOMAS))
covid$FU_time[covid$FU_time<0]<-1
covid$Latencia<-as.numeric(as.Date(covid$FECHA_INGRESO)-as.Date(covid$FECHA_SINTOMAS))
covid$comorb<-covid$DIABETES+covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_d<-covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_ob<-covid$DIABETES+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_dic[covid$comorb>0]<-1;covid$comorb_dic[covid$comorb==0]<-0
covid$diab_ob<-2*covid$DIABETES+covid$OBESIDAD
covid$HABLA_LENGUA_INDIG[is.na(covid$HABLA_LENGUA_INDIG)]<-0
covid$priv[covid$SECTOR==9]<-1; covid$priv[!covid$SECTOR==9]<-0

marg<-read.csv("marg.csv", stringsAsFactors = FALSE)
marg$id<-as.character(str_pad(marg$Clave, 5,pad = "0"))
marg2<-marg%>% dplyr::select(marg,marg_cat,id)

load("df_mx.rda")
dfe3 <- df_mx %>% 
  filter(AÑO==2020) %>% 
  transmute(
    age_levels = sub("-mm", "+", sub("_","-",substr(EDAD_QUIN,6,10))),
    sex = SEXO %>%substr(1,1),
    id = CLAVE,
    value = POB) %<>%
  mutate(age=case_when(
    age_levels=="00-04" | age_levels== "05-09" | age_levels=="10-14" ~ "0-15",
    age_levels== "15-19" ~ "15-19",
    age_levels=="20-24" | age_levels== "25-29" ~ "20-29",
    age_levels=="30-34" | age_levels== "35-39" ~ "30-39",
    age_levels=="40-44" | age_levels== "45-49" ~ "40-49",
    age_levels=="50-54" | age_levels== "55-59" | age_levels=="60-64"~ "50-64",
    age_levels=="65+"~ "65+",
  ))

dfe3$id<-str_pad(dfe3$id, 5,pad = "0")
df_mx$id<-str_pad(df_mx$CLAVE, 5,pad = "0")
pop_mun<-df_mx %>% filter(AÑO==2020)%>%group_by(id)%>%summarise(pop=sum(POB))

camas<-read.csv("dmu_camas.csv")
camas$id<-str_pad(camas$cve_m, 5,pad = "0")

cov<-camas %>% left_join(marg2, by="id") %>% dplyr::select(marg, marg_cat, dmu, id, tot_camas)
cov<-cov %>%left_join(pop_mun, by="id") %>% mutate(camas_adj=tot_camas/pop*10000) %>%dplyr::select(marg, marg_cat,dmu, id, camas_adj)
names(cov)<-c("marg", "marg_cat","dmu", "id", "camas_adj")
covid<- covid %>%left_join(cov, by = "id")
covid$marg<-as.numeric(covid$marg)
covid$marg_cat[covid$marg_cat=="ND"]<-NA
covid$marg_cat[covid$marg_cat=="ND "]<-NA
covid$marg_cat[covid$marg_cat==""]<-NA
covid$marg_cat<-as.factor(covid$marg_cat)
covid<- covid %>% mutate(comorb_cat = comorb %>% 
                           cut(c(0, 1, 2, Inf)))

comorb_levels <- paste0(c("0", "1", ">2"))

covid <- covid %>% 
  mutate(comorb_cat = comorb_cat %>% lvls_revalue(comorb_levels))

covid60<-covid%>%filter(RESULTADO==1, edad60==1)
covid60$edad_cat<-NULL;covid60$edad_cat[covid60$EDAD<70]<-1;covid60$edad_cat[covid60$EDAD>=70 & covid60$EDAD<80]<-2;covid60$edad_cat[covid60$EDAD>=80]<-3


covid1<-covid %>% filter(RESULTADO==1)

### MAPAS ###

mx_mun <- st_read(dsn="shapes", layer="areas_geoestadisticas_municipales", stringsAsFactors=FALSE) %>% 
  st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
bord <-  st_read(dsn="shapes", layer="areas_geoestadisticas_estatales")  %>% 
  st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% 
  rmapshaper::ms_innerlines()
zm <- st_read(dsn="shapes", layer="ZM_2010", stringsAsFactors=FALSE) %>% 
  st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

mx_mun$id<-paste0(str_pad(mx_mun$CVE_ENT, 2,pad = "0"),str_pad(mx_mun$CVE_MUN,3, pad="0"))

#### Incidence rates and mortality by age strata ####
pop<-read.csv("conapo.csv")
names(pop)<-c("r", "year", "ENTIDAD", "CVE", "EDAD", "SEXO", "pop")

age1<-pop %>% filter(year==2020)%>%group_by(SEXO, EDAD)%>%summarise(pop=sum(pop))
covid_age<- covid1 %>% group_by(SEXO, EDAD)%>%summarise(cases=sum(RESULTADO), mort=sum(Mortalidad))
covid_age$SEXO<-factor(covid_age$SEXO, labels = c("Mujeres", "Hombres"))

covid_age1<-left_join(covid_age, age1, by = c("SEXO" = "SEXO", "EDAD" = "EDAD"))
covid_age1$SEXO<-factor(covid_age$SEXO, labels = c("Female", "Male"))


covid_age1<-covid_age1%>%mutate(rate1=cases/pop*100000, rate2=mort/pop*100000)

g1<-ggplot(covid_age1, aes(x=EDAD, y=rate1, col=SEXO))+geom_point()+geom_smooth()+
  theme_classic()+xlim(0, 98)+ylim(0, 75)+ylab("COVID-19 incidence (cases/100,000 habitants)")+
  xlab("Age (years)")+labs(col="Sex")+theme(legend.position="top")+scale_color_jama()
g2<-ggplot(covid_age1, aes(x=EDAD, y=rate2, col=SEXO))+geom_point()+geom_smooth()+
  theme_classic()+xlim(0, 98)+ylim(0, 20)+ylab("COVID-19 mortality (deaths/100,000 habitants)")+
  xlab("Age (years)")+labs(col="Sex")+theme(legend.position="top")+scale_color_jama()

### Histograma de progresión de casos ###

covid4<-covid%>% dplyr::select(FECHA_SINTOMAS,edad60, Mortalidad)%>%drop_na()
covid4$edad60<-as.factor(covid4$edad60)
levels(covid4$edad60)<-c("<60 years", ">60 years")
covid4$Mortalidad<-as.factor(covid4$Mortalidad)
levels(covid4$Mortalidad)<-c("Non-lethal cases", "Lethal cases")
g3<-ggplot(covid4, aes(x=FECHA_SINTOMAS, fill=edad60))+
  geom_histogram(col="black", binwidth = 6)+
  ylab("New confirmed COVID-19 cases")+
  xlab("Symptom onset (date)")+
  theme_classic()+
  geom_vline(xintercept = 5, size = 1, colour = "#FF3721",linetype = "dashed")+
  facet_wrap(~Mortalidad, scales = "free")+
  labs(fill="Age")+theme(legend.position="top")+scale_fill_jama()


cov_comorb<-covid1 %>% group_by(comorb_dic, edad60)%>%
  summarise(hosp=sum(TIPO_PACIENTE), uci=sum(UCI, na.rm=T), intub=sum(INTUBADO, na.rm=T), mort=sum(Mortalidad),n=n()) %>% 
  mutate(freq1 = hosp /n,freq2 = uci / n,freq3 = intub / n, freq4=mort/n)%>%drop_na()
cov_comorb
cov_comorb$edad60<-factor(cov_comorb$edad60, labels = c("<60 years", ">60 years"))
cov_comorb$comorb_dic<-factor(cov_comorb$comorb_dic, labels =c("No-Comorb", "Comorb"))
age1<-c("<60y, 0C", ">60y, 0C","<60y, 1+C", ">60y, 1+C",
        "<60y, 0C", ">60y, 0C","<60y, 1+C", ">60y, 1+C",
        "<60y, 0C", ">60y, 0C","<60y, 1+C", ">60y, 1+C",
        "<60y, 0C", ">60y, 0C","<60y, 1+C", ">60y, 1+C")
freq<-c(cov_comorb$freq1, cov_comorb$freq2, cov_comorb$freq3, cov_comorb$freq4)
des<-c("Hospitalization","Hospitalization","Hospitalization","Hospitalization",
       "ICU admission", "ICU admission", "ICU admission", "ICU admission",
       "Invasive ventilation","Invasive ventilation","Invasive ventilation","Invasive ventilation",
       "Lethality","Lethality","Lethality","Lethality")
cov1<-data.frame(age1, freq, des)

g4<-ggplot(cov1%>% arrange(), aes(y=freq, group=age1,x=age1,fill=age1)) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme_pubr()+ylab("Frequency (%)")+facet_wrap(~des, scales = "free")+
  scale_y_continuous(labels = scales::percent_format())+xlab("")+
  geom_text(position = position_dodge(width= 0.8),aes(label=round(freq*100,2)), vjust=1.6,color="white", size=3.5)+
  scale_fill_jama()+labs(fill="Age")

fig1<-ggarrange(g1, g2,g3, g4,labels = c("A", "B", "C", "D"), nrow=2, ncol=2)

ggsave(fig1,filename = "Figure1.jpg", 
       width = 35, 
       height = 30,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

#### Mixed effects outcome models ####

models<- covid60 %>% dplyr::select(NEUMONIA, HABLA_LENGUA_INDIG,UCI, INTUBADO, TIPO_PACIENTE, SEXO, EDAD, priv,
                            DIABETES, OBESIDAD, TABAQUISMO, EPOC, INMUSUPR, CARDIOVASCULAR,
                            RENAL_CRONICA, marg, id, HIPERTENSION)
names(models)<-c("Pneumonia", "Indigenous","ICU", "Ventilation", "Hospitalization", "Male Sex", "Age", "Private facility",
                 "Diabetes", "Obesity", "Smoking", "COPD", "Immunosupression", "CVD", "CKD", "Social lag index", "id", "Hypertension")

# Pneumonia
m0<-glmer(ICU~Age+`Male Sex`+Indigenous+(1|id), family="binomial", data=models)
summ(m0, exp=T, confint=T)

m1<-glmer(Pneumonia~Age+`Male Sex`+Indigenous+CVD+CKD+COPD+Immunosupression+Smoking+Diabetes+Obesity+Hypertension+`Social lag index`+(1|id), family="binomial", data=models)
summary(m1)
# Hospitalization
m2<-glmer(Hospitalization~Age+`Male Sex`+Indigenous+CVD+CKD+COPD+Immunosupression+Smoking+Diabetes+Obesity+Hypertension+`Social lag index`+(1|id), family="binomial", data=models)


# ICU admission
m3<-glmer(ICU~Age+`Male Sex`+Indigenous+CVD+CKD+COPD+Immunosupression+Smoking+Diabetes+Obesity+`Social lag index`+Hypertension+(1|id), family="binomial", data=models)
summ(m3, exp=T, confint=T)

m4<-glmer(Ventilation~Age+`Male Sex`+Indigenous+CVD+CKD+COPD+Immunosupression+Smoking+Diabetes+Obesity+`Social lag index`+Hypertension+(1|id), family="binomial", data=models)
summ(m4, exp=T, confint=T)

fig2<-plot_summs(m1,m2, m3,m4, ci_level = 0.95, exp = TRUE, scale=TRUE, 
                 model.names = c("Pneumonia", "Hospitalization", "ICU admission", "Ventilation"))+
  xlab("Mixed effects models: OR, 95%CI")

ggsave(fig2,filename = "Figure2.jpg", 
       width = 20, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

## Older adults without comorbidities ###
comorb<- covid60 %>%filter(comorb==0)%>%dplyr::select(NEUMONIA, HABLA_LENGUA_INDIG,UCI, INTUBADO, TIPO_PACIENTE, 
                                                      Mortalidad, SEXO, EDAD, priv, marg, id, FU_time)
nrow(comorb)

m1<-glmer(NEUMONIA~EDAD*marg+SEXO+(1|id), family="binomial",data=comorb)
se <- sqrt(diag(vcov(m1)));tab <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se);exp(tab)

m2<-glmer(UCI~EDAD+marg+SEXO+(1|id), family="binomial",data=comorb)
se <- sqrt(diag(vcov(m2)));tab <- cbind(Est = fixef(m2), LL = fixef(m2) - 1.96 * se, UL = fixef(m2) + 1.96 *se);exp(tab)

m3<-glmer(TIPO_PACIENTE~EDAD+marg+SEXO+(1|id), family="binomial",data=comorb)
se <- sqrt(diag(vcov(m3)));tab <- cbind(Est = fixef(m3), LL = fixef(m3) - 1.96 * se, UL = fixef(m3) + 1.96 *se);exp(tab)

m4<-glmer(INTUBADO~EDAD+marg+SEXO+(1|id), family="binomial",data=comorb)
se <- sqrt(diag(vcov(m4)));tab <- cbind(Est = fixef(m4), LL = fixef(m4) - 1.96 * se, UL = fixef(m4) + 1.96 *se);exp(tab)


m5<-coxph(Surv(FU_time,Mortalidad)~priv+EDAD+marg+NEUMONIA+strata(SEXO)+frailty(id),data=comorb)
summary(m5)
BIC(m5)

#### Effect of age vs. structural vs. comorbidity ####

### Mortality with comorbidities ###

m0<-coxph(Surv(FU_time,Mortalidad)~NEUMONIA+strata(SEXO)+frailty(id),data=covid60)
b0<-BIC(m0)
m1<-coxph(Surv(FU_time,Mortalidad)~EDAD+NEUMONIA+strata(SEXO)+frailty(id),data=covid60)
b1<-BIC(m1)
m2<-coxph(Surv(FU_time,Mortalidad)~priv+marg+NEUMONIA+strata(SEXO)+frailty(id),data=covid60)
b2<-BIC(m2)
m3<-coxph(Surv(FU_time,Mortalidad)~comorb+NEUMONIA+strata(SEXO)+frailty(id),data=covid60)
b3<-BIC(m3)
m4<-coxph(Surv(FU_time,Mortalidad)~EDAD+comorb+NEUMONIA+strata(SEXO)+frailty(id),data=covid60)
b4<-BIC(m4)
m5<-coxph(Surv(FU_time,Mortalidad)~EDAD+priv+marg+NEUMONIA+strata(SEXO)+frailty(id),data=covid60)
b5<-BIC(m5)
m6<-coxph(Surv(FU_time,Mortalidad)~priv+marg+comorb+NEUMONIA+strata(SEXO)+frailty(id),data=covid60)
b6<-BIC(m6)
m7<-coxph(Surv(FU_time,Mortalidad)~EDAD+priv+marg+comorb+NEUMONIA+strata(SEXO)+frailty(id),data=covid60)
b7<-BIC(m7)


delta1<-b1-b0
delta2<-b2-b0
delta5<-b5-b0
delta3<-b3-b0
delta4<-b4-b0
delta6<-b6-b0
delta7<-b7-b0

## Mortality without comorbidities ##

comorb<- covid60 %>%filter(comorb==0)%>%dplyr::select(NEUMONIA, HABLA_LENGUA_INDIG,UCI, INTUBADO, TIPO_PACIENTE, 
                                                      Mortalidad, SEXO, EDAD, priv, marg, id, FU_time)

m0<-coxph(Surv(FU_time,Mortalidad)~NEUMONIA+strata(SEXO)+frailty(id),data=comorb)
b8<-BIC(m0)
m1<-coxph(Surv(FU_time,Mortalidad)~EDAD+NEUMONIA+strata(SEXO)+frailty(id),data=comorb)
b9<-BIC(m1)
m2<-coxph(Surv(FU_time,Mortalidad)~priv+marg+NEUMONIA+strata(SEXO)+frailty(id),data=comorb)
b10<-BIC(m2)
m6<-coxph(Surv(FU_time,Mortalidad)~EDAD+priv+marg+NEUMONIA+strata(SEXO)+frailty(id),data=comorb)
b11<-BIC(m6)


delta8<-b9-b8
delta9<-b10-b8
delta10<-b11-b8

### Figure ###
deltas<-c(delta1, delta2, delta5, delta3, delta4, delta6, delta7, delta8, delta9, delta10)
models<-c("M1", "M2", "M3", "M4", "M5", "M6", "M7","M1", "M2", "M3")
com<-c(rep(0, 7), rep(1,3))
com<-factor(com, labels = c("Comorbidities", "No comorbidities"))

bic<-data.frame(deltas, models, com)

f1<-ggplot(bic, aes(x=models, y=deltas, fill=models))+geom_bar(stat="identity")+
  xlab("")+ylab("BIC change")+facet_wrap(~com, scales="free", shrink=TRUE)+
  scale_fill_grey()+theme_classic()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+labs(fill="Models")

### Mortality in subjects without comorbidity ##

m6<-coxph(Surv(FU_time,Mortalidad)~priv+marg+NEUMONIA+strata(SEXO)+frailty(id),data=comorb)
summary(m6)
BIC(m6)

m7<-coxph(Surv(FU_time,Mortalidad)~priv+marg+EDAD+NEUMONIA+strata(SEXO)+frailty(id),data=comorb)
summary(m7)
BIC(m7)
#### Cox models with frailty penalty ####

## Univariate structural models ###

u1 <- coxph(Surv(FU_time,Mortalidad)~EDAD+marg+strata(SEXO)+frailty(id), data=covid60, x=T)
summary(u1)

u2 <- coxph(Surv(FU_time,Mortalidad)~EDAD+camas_adj+strata(SEXO)+frailty(id), data=covid60, x=T)
summary(u2)

u3 <- coxph(Surv(FU_time,Mortalidad)~priv+EDAD+strata(SEXO)+frailty(id), data=covid60, x=T)
summary(u3)


## Multivariate models ###

mod1 <- coxph(Surv(FU_time,Mortalidad)~priv+OBESIDAD+DIABETES+RENAL_CRONICA+NEUMONIA+INMUSUPR+factor(edad_cat)+marg+strata(SEXO)+frailty(id), data=covid60, x=T)
summary(mod1)
HR<-as.data.frame(cbind(exp(coef(mod1)),exp(confint(mod1))))
HR$Covariate<-c("Private care","Obesity","Diabetes","CKD", "Pneumonia", "Immunosupression",
                "Age 70-80", "Age >80","Social Lag Index")
colnames(HR)<-c("HR", "ciLow", "ciHigh", "Covariate")
HR$breaks<-seq(1:nrow(HR))
h1 <- ggplot(HR, aes(x = HR, y = breaks)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ciHigh, xmin = ciLow), size = .5, height = .2, color = "black") +
  geom_point(size = 3.5, color = "darkblue") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = HR$breaks, labels = HR$Covariate) +
  ylab("") +
  xlab("COVID-19 lethality in older adults, Hazard ratio (HR, 95%CI)")

mod0 <- coxph(Surv(FU_time,Mortalidad)~priv+OBESIDAD+DIABETES+RENAL_CRONICA+EDAD+NEUMONIA+INMUSUPR+marg+strata(SEXO)+frailty(id), data=covid60, x=T)
summary(mod0)

modc <- coxph(Surv(FU_time,Mortalidad)~OBESIDAD+DIABETES+RENAL_CRONICA+NEUMONIA+INMUSUPR+strata(SEXO)+frailty(id), data=covid60, x=T)
summary(modc)

modcs <- coxph(Surv(FU_time,Mortalidad)~EDAD+priv+marg+OBESIDAD+DIABETES+RENAL_CRONICA+NEUMONIA+INMUSUPR+strata(SEXO)+frailty(id), data=covid60, x=T)
summary(modcs)

mod2 <- coxph(Surv(FU_time,Mortalidad)~priv+marg+OBESIDAD+DIABETES+EDAD+RENAL_CRONICA+NEUMONIA+INMUSUPR+strata(SEXO)+frailty(ENTIDAD_RES), data=covid60)
summary(mod2)


### Coxme ###

mod1.1 <- coxme(Surv(FU_time,Mortalidad)~priv+marg+OBESIDAD+DIABETES+EDAD+EPOC+RENAL_CRONICA+NEUMONIA+INMUSUPR+strata(SEXO)+(1|id), data=covid60)
summary(mod1.1)

mod1.2 <- coxme(Surv(FU_time,Mortalidad)~priv+marg+OBESIDAD+EPOC+DIABETES+EDAD+EPOC+RENAL_CRONICA+NEUMONIA+INMUSUPR+strata(SEXO)+(1|ENTIDAD_RES), data=covid60)
summary(mod1.2)

#### Simulations ####
set.seed(123)
covid60$age_med<-covid60$EDAD-median(covid60$EDAD)
covid1$age_med<-covid1$EDAD-median(covid1$EDAD)

mod1 <- coxph(Surv(FU_time,Mortalidad)~priv+age_med+comorb+marg+strata(SEXO)+frailty(id), data=covid60)
summary(mod1)

Sim1 <- coxsimLinear(mod1, b = "age_med", Xj = seq(-20,20),  Xl = seq(0, 0), qi = "Hazard Ratio",spin = TRUE, ci = 0.95)
p1<-simGG(Sim1, alpha = 0.35, type = "ribbons", xlab="Age from median (years)", lcolour="black", rcolour="darkgray")
Sim2 <- coxsimLinear(mod1, b = "comorb", Xj = seq(0,9),  Xl = seq(1,1), qi = "Hazard Ratio",spin = TRUE, ci = 0.95)
p2<-simGG(Sim2, alpha = 0.35, type = "ribbons", xlab="Number of comorbidities", lcolour="black", rcolour="darkgray")+
  scale_x_discrete(limits=0:9, labels=0:9)
Sim3 <- coxsimLinear(mod1, b = "marg", Xj = seq(-2,2),  Xl = seq(-1,-1), qi = "Hazard Ratio",spin = TRUE, ci = 0.95)
p3<-simGG(Sim3, alpha = 0.35, type = "ribbons", xlab="Social lag index", lcolour="black", rcolour="darkgray")

fig3<-ggarrange(ggarrange(h1,f1, labels=c("A", "B"), ncol=2),
                ggarrange(p1, p2, p3, labels=c("C", "D", "E"), ncol=3), nrow=2)

ggsave(fig3,filename = "Figure3.jpg", 
       width = 35, 
       height = 25,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

#### Mortality risk maps ####
pop<-df_mx %>% filter(AÑO==2020, EDAD_QUIN=="pobm_65_mm")%>% group_by(id) %>%
  summarise(pop=sum(POB))%>%dplyr::select(id, pop)

covid_casos<- covid60 %>% left_join(pop, by="id") 
covid_casos<- covid_casos %>%group_by(id) %>%
  summarise(cases=n(), deaths=sum(Mortalidad), pop=median(pop))%>%
  mutate(rate1=cases/pop*100000, rate2=deaths/pop*100000)%>%
  dplyr::select(id, rate1, rate2, pop)

mun.risk<-(ranef(mod1.1)[[1]])
mun.risk<-data.frame(mun.risk, labels(mun.risk))
mun.risk$mun.risk<-exp(mun.risk$mun.risk)
names(mun.risk)<-c("risk", "id")

state.risk<-(ranef(mod1.2)[[1]])
state.risk<-data.frame(state.risk, labels(state.risk))
state.risk$state.risk<-exp(state.risk$state.risk)
names(state.risk)<-c("risk2", "state")

mx_mun$state<-mx_mun$CVE_ENT
mx_mun <- mx_mun %>% left_join(mun.risk, by="id")
mx_mun<-mx_mun %>%left_join(covid_casos, by="id")
mx_mun<-mx_mun %>% left_join(state.risk, "state")

#Juntamos información de mapa con riesgo
map1<-ggplot() +
  geom_sf(data=mx_mun, aes(fill = rate1, geometry=geometry), color = NA) +
  geom_sf(data = bord, color = "darkgray", size = .3)+
  labs(fill="Incidence rate")+
  guides(fill = guide_legend(reverse = FALSE))+theme_map()+
  scale_fill_viridis_c(trans="log", breaks=c(5, 10, 50, 200)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
map2<-ggplot() +
  geom_sf(data=mx_mun, aes(fill = rate2, geometry=geometry), color = NA) +
  geom_sf(data = bord, color = "darkgray", size = .3)+
  labs(fill="Mortality rate")+
  guides(fill = guide_legend(reverse = FALSE))+theme_map()+
  scale_fill_viridis_c(trans="log", breaks=c(5, 10, 50, 200)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
map3<-ggplot() +
  geom_sf(data=mx_mun, aes(fill = risk, geometry=geometry), color = NA) +
  geom_sf(data = bord, color = "darkgray", size = .3)+
  labs(fill="Hazard ratio")+
  guides(fill = guide_legend(reverse = FALSE))+theme_map()+
  scale_fill_viridis_c(trans="log", breaks=c(0.8,1.2, 1.4, 1.6, 1.8)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
cdmx<-ggplot() +
  geom_sf(data=mx_mun, aes(fill = risk, geometry=geometry), color = NA) +
  geom_sf(data = bord, color = "darkgray", size = .7)+
  coord_sf(xlim = c( -99.453209, -98.833209), ylim = c(19.00608, 19.63608), expand = FALSE) +
  theme_map()+
  scale_fill_viridis_c(trans="log",) +
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

out <- ggdraw()+
  draw_plot(map3)+
  draw_plot(cdmx, x = .10, y =.02 , width = .30, height = .32)

#### Lollipop plot ####
pal <- viridis::viridis(6,direction=-1) %>% rev 
pal.25 <- pal %>% clr_darken(shift = .04)
conapo <- read_csv("conapo.csv")
states<- conapo %>% filter(AÑO=="2020") %>% group_by(ENTIDAD) %>% summarise(state=mean(CVE_GEO))
poptot<- df_mx %>% filter(AÑO=="2020") %>%group_by(id)%>%
  summarise(pop_tot=sum(POB))
states$ENTIDAD<-c("AGS", "BCN", "BCS", "CAM", "CHP", "CHH", "CMX", "COA",
                  "COL", "DUR", "GUA", "GRO", "HID", "JAL", "MEX", "MIC", 
                  "MOR", "NAY", "NLE", "OAX", "PUE", "QUE", "ROO", "SLP",
                  "SIN", "SON", "TAB", "TAM", "TLA", "VER", "YUC", "ZAC")
names(states)<-c("name", "state")
states$state<-str_pad(states$state, 2, pad = "0")
mx_mun$id2<-mx_mun$CVE_ENT
mx_mun<- mx_mun %>% mutate(risk1 = risk %>% 
    cut(c(0.5, 0.8, 1, 1.2, 1.4, 1.6,Inf)))

# calculate the group sizes
new_levels <- paste0(c("0.5-0.8", "0.8-1.0", "1.0-1.2", "1.2-1.4", "1.4-1.6", ">1.6"))

mx_mun <- mx_mun %>% 
  mutate(risk1 = risk1 %>% lvls_revalue(new_levels))

mx_mun<- mx_mun %>% left_join(poptot, "id")
mx_mun<-mx_mun %>% left_join(camas, "id") %>%mutate(camas_adj=tot_camas/pop_tot*1000)

map4<-mx_mun %>% 
  mutate(state = id %>% str_sub(1, 2)) %>% 
  left_join(states) %>% 
  group_by(name) %>% 
  ungroup() %>% 
  drop_na() %>% 
  # arrange by decreasing average proprotion
  arrange(risk2 %>% desc) %>% 
  mutate(name = name %>% fct_inorder %>% fct_rev) %>%
  ggplot(aes(risk, name, col = risk1, size = camas_adj))+
  geom_vline(xintercept = 1, size = 1, 
             color = "#B8B8B8", alpha = .5)+
  geom_point(shape = 1)+
  geom_point(aes(x = risk2), shape = 124, color = "cyan", size = 4)+
  scale_color_manual(values = pal.25, guide = NULL)+
  scale_size_area(max_size = 8, breaks = c(0.5, 5, 10, 25))+
  scale_y_discrete(position = "right")+
  theme_minimal()+
  theme(legend.position = c(.95, .18),
        text = element_text(lineheight = .9, size=7),
        axis.title.x = element_text(hjust = .5),
        legend.spacing.y = unit(0.001, "cm"))+
  labs(x = "COVID-19 mortality risk (HR)",
       y = NULL,
       size = "Beds \nper 1,000 \nhabitants")
maps<-ggarrange(out, map4,labels = c("A", "B"))

ggsave(maps,filename = "Figure4.jpg", 
       width = 30, 
       height = 10,
       units=c("cm"),
       dpi = 600,
       limitsize = FALSE)


#### Metropolitan areas ####
zm$id<-zm$CVE_MUN1
zm2<- zm %>% dplyr::select(id, CVE_SUN) %>%as.data.frame()
metro<-covid60 %>% left_join(zm2, "id")

# All municipalities
m1<-coxph(Surv(FU_time, Mortalidad)~priv+EDAD*marg+comorb+marg+NEUMONIA+strata(SEXO)+frailty(id), data=covid60)
summary(m1) 

#Metropolitan areas
m2<-coxph(Surv(FU_time, Mortalidad)~marg+comorb+priv+EDAD+comorb+NEUMONIA+strata(SEXO)+frailty(CVE_SUN), data=metro)
summary(m2)   

#Non-metropolitan areas
non.zm<- metro%>% filter(is.na(CVE_SUN))

m3<-coxph(Surv(FU_time, Mortalidad)~EDAD+marg+comorb+priv+NEUMONIA+strata(SEXO)+frailty(id), data=non.zm)
summary(m3)

#### Obesity paradox ####
ob<- covid60 %>% filter(comorb_ob==0)

m.ob<-coxph(Surv(FU_time, Mortalidad)~EDAD+OBESIDAD+marg+priv+NEUMONIA+strata(SEXO)+frailty(id), data=ob)
summary(m.ob)

m.ob1<-glmer(INTUBADO~EDAD+SEXO+marg+OBESIDAD+(1|id), family="binomial",data=ob)
se <- sqrt(diag(vcov(m.ob1)));tab <- cbind(Est = fixef(m.ob1), LL = fixef(m.ob1) - 1.96 * se, UL = fixef(m.ob1) + 1.96 *se);exp(tab)

m.ob1<-glmer(UCI~EDAD+SEXO+marg+OBESIDAD+(1|id), family="binomial",data=ob)
se <- sqrt(diag(vcov(m.ob1)));tab <- cbind(Est = fixef(m.ob1), LL = fixef(m.ob1) - 1.96 * se, UL = fixef(m.ob1) + 1.96 *se);exp(tab)

m.ob1<-glmer(TIPO_PACIENTE~EDAD+SEXO+marg+OBESIDAD+(1|id), family="binomial",data=ob)
se <- sqrt(diag(vcov(m.ob1)));tab <- cbind(Est = fixef(m.ob1), LL = fixef(m.ob1) - 1.96 * se, UL = fixef(m.ob1) + 1.96 *se);exp(tab)

m.ob1<-glmer(NEUMONIA~EDAD+SEXO+marg+OBESIDAD+(1|id), family="binomial",data=ob)
se <- sqrt(diag(vcov(m.ob1)));tab <- cbind(Est = fixef(m.ob1), LL = fixef(m.ob1) - 1.96 * se, UL = fixef(m.ob1) + 1.96 *se);exp(tab)

