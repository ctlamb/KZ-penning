library(tidyverse)
library(here)
library(survival)
library(lubridate)
library(MuMIn)
library(survminer)

##load data
df <- read.csv(here::here("data","Klinse-zaEackerInputs_200607_CL.csv"))

##keep only 2014+ data
df <- df%>%
  filter(Year>=2014)

##make real dates
df <- df%>%
  mutate(enter=paste(Year,enter+1,15,sep="/")%>%ymd%m-%months(1),
         exit=paste(Year,exit,28,sep="/")%>%ymd)

##%m-%months(3) if real dates but these are already changed

##set penned dates
penned <- tribble(
  ~year,~into, ~out,
  2014, "2014-03-15", "2014-07-05",
  2015, "2015-03-15", "2015-07-24",
  2016, "2016-03-15", "2016-07-15",
  2017, "2017-03-15", "2017-07-27",
  2018, "2018-03-15", "2018-07-31",
  2019, "2019-03-15", "2019-07-29")%>%
  mutate(into=ymd(into)%m-%months(2),
         out=ymd(out)%m-%months(2),
         int=into%--%out)



###Assign in or out of pen and 
surv.season <- data.frame()
for(i in 1:nrow(df)){
  a <- df%>%slice(i)
  col.int <- a$enter%--%a$exit
  pen.int <- penned%>%
    filter(year==a$Year[1])%>%
    pull(int)

  
  if(int_start(col.int)%within%pen.int & int_end(col.int)>int_end(pen.int)){
    
    ints <-int_diff(c(int_start(col.int),int_end(pen.int),int_end(col.int)))
    b <- data.frame(id=a$id[1],
                    year=a$Year[1],
                    event=c(rep(0, times=length(ints)-1),max(a$event)),
                    enter.date=int_start(ints),
                    exit.date=int_end(ints),
                    group=a$HerdCode[1],
                    pen=c("in","out"))
  }

  
  if(a$HerdCode[1]%in%"K-Z-Pen" & int_end(col.int)<int_end(pen.int)){
    
    b <- data.frame(id=a$id[1],
                    year=a$Year[1],
                    event=a$event[1],
                    enter.date=int_start(col.int)[1],
                    exit.date=int_end(col.int)[1],
                    group=a$HerdCode[1],
                    pen=c("in"))
  }
  
  if( !int_start(col.int)%within%pen.int ){
    
    b <- data.frame(id=a$id[1],
                    year=a$Year[1],
                    event=a$event[1],
                    enter.date=int_start(col.int)[1],
                    exit.date=int_end(col.int)[1],
                    group=a$HerdCode[1],
                    pen=c("out"))
  }
  
  
  if( int_end(col.int)%within%pen.int ){
    
    ints <-int_diff(c(int_start(pen.int),int_end(col.int)))
    
    b <- data.frame(id=a$id[1],
                    year=a$Year[1],
                    event=c(rep(0, times=length(ints)-1),max(a$event)),
                    enter.date=int_start(ints),
                    exit.date=int_end(ints),
                    group=a$HerdCode[1],
                    pen=c("in"))
  }


  surv.season <- rbind(surv.season, b)
  rm(b)
  rm(a)

}

##get duration
surv.season<-surv.season%>%
  mutate(time=difftime(exit.date,enter.date, units="days")%>%as.numeric)


##model survival
m1 <- coxph(Surv(time, event)~ pen +
              frailty(id), data = surv.season)
m2 <- coxph(Surv(time, event)~ group +
              frailty(id), data = surv.season)
m3 <- coxph(Surv(time, event)~ group + pen +
              frailty(id), data = surv.season)
m4 <- coxph(Surv(time, event)~ group*pen +
              frailty(id), data = surv.season)

summary(m4)

model.sel(m1,m2,m3,m4,  rank="AICc")
cox.zph(m3)
summary(m4)


m1 <- survfit(Surv(time, event)~ pen, data = surv.season)
m2 <- survfit(Surv(time, event)~ group, data = surv.season)
m3 <- survfit(Surv(time, event)~ group + pen, data = surv.season)


ggsurvplot(m1, data = surv.season)
ggsurvplot(m2, data = surv.season)
ggsurvplot(m3, data = surv.season)



##load data
df <- read.csv(here::here("data","Klinse-zaEackerInputs_200607_CL.csv"))

##keep only 2014+ data
df <- df%>%
  filter(Year>=2014)

##make real dates
df <- df%>%
  mutate(enter=paste(Year,enter+1,1,sep="/")%>%ymd%>%yday%>%as.numeric,
         exit=paste(Year,exit,28,sep="/")%>%ymd%>%yday%>%as.numeric)
  


m1 <- survfit(Surv(time=enter, time2=exit, event=event)~ HerdCode, data = df)
ggsurvplot(m1, data = df)

df[df$enter>df$exit,]
