library(readr); library(tidyverse); library(survival); library(mediation); library(ggpubr); library(rms);library(caret)
library(survminer); library(haven); library(rsq); library(ResourceSelection); library(ggsci);library(timereg);
library(ggplot2);library(gridExtra);library(grid);library(officer);library(flextable)

## Manejo de la base de datos####
setwd("C:/Users/facmed/UNIVERSIDAD NACIONAL AUT?NOMA DE M?XICO/OMAR YAXMEHEN BELLO CHAVOLLA - COVID - Envejecimiento")

covid <- read_csv("200509COVID19MEXICO.csv")

table(covid$RESULTADO)

covid$id<-paste0(str_pad(covid$ENTIDAD_RES, 2,pad = "0"),str_pad(covid$MUNICIPIO_RES,3, pad="0"))
covid1<-covid[,c(14:15,18:30,32:35)]
covid<-covid[,-c(14:15,18:30,32:35)]
covid1[covid1==2]<-0
covid1[covid1==97]<-NA;covid1[covid1==98]<-NA;covid1[covid1==99]<-NA
covid<-as.data.frame(cbind(covid, covid1))
covid$EMBARAZO[is.na(covid$EMBARAZO)]<-0
covid$TIPO_PACIENTE[covid$TIPO_PACIENTE==1]<-0;covid$TIPO_PACIENTE[covid$TIPO_PACIENTE==2]<-1
covid$covid<-NULL; covid$covid[covid$RESULTADO==1]<-1;covid$covid[covid$RESULTADO!=1]<-0
covid$edad65<-NULL;covid$edad65[covid$EDAD>=65]<-1;covid$edad65[covid$EDAD<65]<-0
covid$edad40<-NULL;covid$edad40[covid$EDAD>=40]<-0;covid$edad40[covid$EDAD<40]<-1
covid$diabetes_40<-NULL;covid$diabetes_40[covid$DIABETES==1 & covid$edad40==1]<-1;covid$diabetes_40[covid$DIABETES!=1 & covid$edad40!=1]<-0
covid$Mortalidad<-NULL; covid$Mortalidad[is.na(covid$FECHA_DEF)]<-0;covid$Mortalidad[is.na(covid$FECHA_DEF)==FALSE]<-1
covid$FECHA_DEF[is.na(covid$FECHA_DEF)]<-as.Date(Sys.Date())
covid$FU_time<-as.numeric(as.Date(covid$FECHA_DEF)-as.Date(covid$FECHA_SINTOMAS))
covid$FU_time[covid$FU_time>30]<-30
covid$Latencia<-as.numeric(as.Date(covid$FECHA_INGRESO)-as.Date(covid$FECHA_SINTOMAS))
covid$comorb<-covid$DIABETES+covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_d<-covid$OBESIDAD+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$comorb_ob<-covid$DIABETES+covid$EPOC+covid$ASMA+covid$INMUSUPR+covid$HIPERTENSION+
  covid$CARDIOVASCULAR+covid$RENAL_CRONICA+covid$TABAQUISMO+covid$OTRA_COM
covid$diab_ob<-2*covid$DIABETES+covid$OBESIDAD

covid_am<-covid%>%filter(RESULTADO==1,EDAD>=60)
covid_j<-covid%>%filter(RESULTADO==1, EDAD<60)

names(covid)

#Sexo
J_11.1<-table(covid_j$SEXO)[2];J_11.2<-(table(covid_j$SEXO)/nrow(covid_j))[2]
AM_11.1<-table(covid_am$SEXO)[2];AM_11.2<-(table(covid_am$SEXO)/nrow(covid_am))[2]
#Diabetes
J_1.1<-table(covid_j$DIABETES)[2];J_1.2<-(table(covid_j$DIABETES)/nrow(covid_j))[2]
AM_1.1<-table(covid_am$DIABETES)[2];AM_1.2<-(table(covid_am$DIABETES)/nrow(covid_am))[2]
#EPOC
J_2.1<-table(covid_j$EPOC)[2];J_2.2<-(table(covid_j$EPOC)/nrow(covid_j))[2]
AM_2.1<-table(covid_am$EPOC)[2];AM_2.2<-(table(covid_am$EPOC)/nrow(covid_am))[2]
#ASMA
J_3.1<-table(covid_j$ASMA)[2];J_3.2<-(table(covid_j$ASMA)/nrow(covid_j))[2]
AM_3.1<-table(covid_am$ASMA)[2];AM_3.2<-(table(covid_am$ASMA)/nrow(covid_am))[2]
#INMUSUPR
J_4.1<-table(covid_j$INMUSUPR)[2];J_4.2<-(table(covid_j$INMUSUPR)/nrow(covid_j))[2]
AM_4.1<-table(covid_am$INMUSUPR)[2];AM_4.2<-(table(covid_am$INMUSUPR)/nrow(covid_am))[2]
#HIPERTENSION
J_5.1<-table(covid_j$HIPERTENSION)[2];J_5.2<-(table(covid_j$HIPERTENSION)/nrow(covid_j))[2]
AM_5.1<-table(covid_am$HIPERTENSION)[2];AM_5.2<-(table(covid_am$HIPERTENSION)/nrow(covid_am))[2]
#OTRA_COM
J_6.1<-table(covid_j$OTRA_COM)[2];J_6.2<-(table(covid_j$OTRA_COM)/nrow(covid_j))[2]
AM_6.1<-table(covid_am$OTRA_COM)[2];AM_6.2<-(table(covid_am$OTRA_COM)/nrow(covid_am))[2]
#CARDIOVASCULAR
J_7.1<-table(covid_j$CARDIOVASCULAR)[2];J_7.2<-(table(covid_j$CARDIOVASCULAR)/nrow(covid_j))[2]
AM_7.1<-table(covid_am$CARDIOVASCULAR)[2];AM_7.2<-(table(covid_am$CARDIOVASCULAR)/nrow(covid_am))[2]
#OBESIDAD
J_8.1<-table(covid_j$OBESIDAD)[2];J_8.2<-(table(covid_j$OBESIDAD)/nrow(covid_j))[2]
AM_8.1<-table(covid_am$OBESIDAD)[2];AM_8.2<-(table(covid_am$OBESIDAD)/nrow(covid_am))[2]
#RENAL_CRONICA
J_9.1<-table(covid_j$RENAL_CRONICA)[2];J_9.2<-(table(covid_j$RENAL_CRONICA)/nrow(covid_j))[2]
AM_9.1<-table(covid_am$RENAL_CRONICA)[2];AM_9.2<-(table(covid_am$RENAL_CRONICA)/nrow(covid_am))[2]
#TABAQUISMO
J_10.1<-table(covid_j$TABAQUISMO)[2];J_10.2<-(table(covid_j$TABAQUISMO)/nrow(covid_j))[2]
AM_10.1<-table(covid_am$TABAQUISMO)[2];AM_10.2<-(table(covid_am$TABAQUISMO)/nrow(covid_am))[2]
#NEUMONIA
J_12.1<-table(covid_j$NEUMONIA)[2];J_12.2<-(table(covid_j$NEUMONIA)/nrow(covid_j))[2]
AM_12.1<-table(covid_am$NEUMONIA)[2];AM_12.2<-(table(covid_am$NEUMONIA)/nrow(covid_am))[2]
#HOSPITALIZATION
J_13.1<-table(covid_j$TIPO_PACIENTE)[2];J_13.2<-(table(covid_j$TIPO_PACIENTE)/nrow(covid_j))[2]
AM_13.1<-table(covid_am$TIPO_PACIENTE)[2];AM_13.2<-(table(covid_am$TIPO_PACIENTE)/nrow(covid_am))[2]
#ICU
J_14.1<-table(covid_j$INTUBADO)[2];J_14.2<-(table(covid_j$INTUBADO)/nrow(covid_j))[2]
AM_14.1<-table(covid_am$INTUBADO)[2];AM_14.2<-(table(covid_am$INTUBADO)/nrow(covid_am))[2]
#DEATH
J_15.1<-table(covid_j$Mortalidad)[2];J_15.2<-(table(covid_j$Mortalidad)/nrow(covid_j))[2]
AM_15.1<-table(covid_am$Mortalidad)[2];AM_15.2<-(table(covid_am$Mortalidad)/nrow(covid_am))[2]


AM_1<-paste0(AM_1.1," ","(",round(AM_1.2*100,1),")");J_1<-paste0(J_1.1," ","(",round(J_1.2*100,1),")")
AM_2<-paste0(AM_2.1," ","(",round(AM_2.2*100,1),")");J_2<-paste0(J_2.1," ","(",round(J_2.2*100,1),")")
AM_3<-paste0(AM_3.1," ","(",round(AM_3.2*100,1),")");J_3<-paste0(J_3.1," ","(",round(J_3.2*100,1),")")
AM_4<-paste0(AM_4.1," ","(",round(AM_4.2*100,1),")");J_4<-paste0(J_4.1," ","(",round(J_4.2*100,1),")")
AM_5<-paste0(AM_5.1," ","(",round(AM_5.2*100,1),")");J_5<-paste0(J_5.1," ","(",round(J_5.2*100,1),")")
AM_6<-paste0(AM_6.1," ","(",round(AM_6.2*100,1),")");J_6<-paste0(J_6.1," ","(",round(J_6.2*100,1),")")
AM_7<-paste0(AM_7.1," ","(",round(AM_7.2*100,1),")");J_7<-paste0(J_7.1," ","(",round(J_7.2*100,1),")")
AM_8<-paste0(AM_8.1," ","(",round(AM_8.2*100,1),")");J_8<-paste0(J_8.1," ","(",round(J_8.2*100,1),")")
AM_9<-paste0(AM_9.1," ","(",round(AM_9.2*100,1),")");J_9<-paste0(J_9.1," ","(",round(J_9.2*100,1),")")
AM_10<-paste0(AM_10.1," ","(",round(AM_10.2*100,1),")");J_10<-paste0(J_10.1," ","(",round(J_10.2*100,1),")")
AM_11<-paste0(AM_11.1," ","(",round(AM_11.2*100,1),")");J_11<-paste0(J_11.1," ","(",round(J_11.2*100,1),")")
AM_12<-paste0(AM_12.1," ","(",round(AM_12.2*100,1),")");J_12<-paste0(J_12.1," ","(",round(J_12.2*100,1),")")
AM_13<-paste0(AM_13.1," ","(",round(AM_13.2*100,1),")");J_13<-paste0(J_13.1," ","(",round(J_13.2*100,1),")")
AM_14<-paste0(AM_14.1," ","(",round(AM_14.2*100,1),")");J_14<-paste0(J_14.1," ","(",round(J_14.2*100,1),")")
AM_15<-paste0(AM_15.1," ","(",round(AM_15.2*100,1),")");J_15<-paste0(J_15.1," ","(",round(J_15.2*100,1),")")


P_1<-round(prop.test(x=c(J_1.1,AM_1.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value,3);P_1
P_2<-prop.test(x=c(J_2.1,AM_2.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_3<-prop.test(x=c(J_3.1,AM_3.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_4<-prop.test(x=c(J_4.1,AM_4.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_5<-prop.test(x=c(J_5.1,AM_5.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_6<-prop.test(x=c(J_6.1,AM_6.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_7<-prop.test(x=c(J_7.1,AM_7.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_8<-prop.test(x=c(J_8.1,AM_8.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_9<-prop.test(x=c(J_9.1,AM_9.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_10<-prop.test(x=c(J_10.1,AM_10.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_11<-prop.test(x=c(J_11.1,AM_11.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_12<-prop.test(x=c(J_12.1,AM_12.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_13<-prop.test(x=c(J_13.1,AM_13.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_14<-prop.test(x=c(J_14.1,AM_14.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
P_15<-prop.test(x=c(J_15.1,AM_15.1),n=c(nrow(covid_j),nrow(covid_am)))$p.value
GGG<-base::format.pval(c(P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9,P_10,P_11,P_12,P_13,P_14,P_15), eps = .001, digits = 2)

AM_JOV_COVID<-data.frame("bold(Comorbidities)"=c("Diabetes (%)","COPD (%)","Asthma (%)","Immunosuppression (%)","Hypertension (%)",
                                        "Other (%)","CVD (%)","Obesity","CKD (%)","Smoking (%)","Men(%)", "Pneumonia(%)",
                                        "Hospitalization(%)","ICU admission (%)", "Death (%)"),
                    "Age<60\nN=nrows(covid_am)"=c(J_1,J_2,J_3,J_4,J_5,J_6,J_7,J_8,J_9,J_10,J_11,J_12,J_13,J_14,J_15),
                    "Age>=60\nN=nrows(covid_am)"=c(AM_1,AM_2,AM_3,AM_4,AM_5,AM_6,AM_7,AM_8,AM_9,AM_10,AM_11,AM_12,AM_13,AM_14,AM_15),
                    "p-value"=GGG)

colnames(AM_JOV_COVID)<-c("Comorbidities",paste0("Age<60 \nn=",nrow(covid_j)),
                          paste0("Age>=60 \nn=",nrow(covid_am)),"p-value")

AM_JOV_COVID<-flextable(AM_JOV_COVID,cwidth = 0.5*ncol(AM_JOV_COVID))
AM_JOV_COVID
save_as_docx(AM_JOV_COVID,path = "table1.docx")
