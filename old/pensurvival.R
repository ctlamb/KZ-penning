library(survival)
library(here)
library(lubridate)
library(MuMIn)
library(survminer)
library(hrbrthemes)
library(tidyverse)
library(tidylog)
source("seasonalsurvival_fn.r")

#############################
####ADULT FEMALE SURVIVAL####
#############################

##load data
df <- read_csv(here::here("data","CollarFates4CL-RS_201230.csv"))%>%
  filter(herd%in%c("K","S"),
         Sex%in%"F")
  

pen.yr <- read_csv(here::here("data","Klinse-za capture and demographics_210406_excerpt.csv"))


##get min age via birthdate
df <- df %>% 
  mutate(rstart_dat=dmy(rstart_dat),
         rend_dat=dmy(rend_dat),
    birth=case_when(agecl%in%"C" & month(rstart_dat) %in% 5:12~ paste(year(rstart_dat),"-06-01",sep="-")%>%ymd,
                    agecl%in%"C" & month(rstart_dat) %in% 1:4~ paste(year(rstart_dat)-1,"-06-01",sep="-")%>%ymd,
                    agecl%in%"Y" & month(rstart_dat) %in% 5:12~ paste(year(rstart_dat)-1,"-06-01",sep="-")%>%ymd,
                    agecl%in%"Y" & month(rstart_dat) %in% 1:4~ paste(year(rstart_dat)-2,"-06-01",sep="-")%>%ymd,
                    agecl%in%"A"~ paste(year(rstart_dat)-3,"-06-01",sep="-")%>%ymd))

##classify dead or alive
#one myop but happens in 2002. Watch for myop pen 2020 (don't see here yet)
#x
unique(df$mortcaus)
df <- df %>% 
  filter(!mortcaus %in%"MYOP")%>%
  mutate(dead=case_when(mortcaus%in%c("DEAD")~1,
                        mortcaus%in%c("FADE")~0))

##get columns named appropriately for function
df <- df %>% 
  select(id=animal_id,
         herd,
         start=rstart_dat,
         end=rend_dat,
         dead,
         birth,
         agecl)

##find issues any date issues
df <- df%>%
  mutate(dur=end-start)

df%>%
  filter(dur<1)

##make daily
surv <- stretch_survival_data(df, '1 day')

nrow(surv)
sum(df$dur) ##should be same

##set penned dates
pen <- tribble(
  ~year,~into, ~out,
  2014, "2014-03-27", "2014-07-05",
  2015, "2015-03-18", "2015-07-24",
  2016, "2016-03-16", "2016-07-15",
  2017, "2017-03-23", "2017-07-27",
  2018, "2018-03-13", "2018-07-31",
  2019, "2019-03-12", "2019-07-29",
  2020, "2020-03-11", "2020-08-05",
  2021, "2021-03-15", NA)%>%
  mutate(free_end=(ymd(dplyr::lead(into,1,order_by=year))-1)%>%as.character,
      penned=ymd(into)%--%ymd(out),
      free=ymd(out)%--%ymd(free_end))

##Assign penned time
surv <- surv%>%
  mutate(dur=(end-start)%>%as.numeric(),
         year=year(end),
         season=case_when(start %within% as.list(pen$penned)~"Pen", TRUE~"Free"),
         mgmt=case_when(year<2014~"pre",year>=2014~"mgmt"))

##clean up id's (get rid of trailing letter denoting period)
surv <- surv %>% 
  mutate(id=str_sub(id,1,-2))

ggplot(surv, aes(x=end,y=id,color=season))+
  geom_point(size=0.01)+
  facet_wrap(vars(season))


##decide who was in the pen
pen.yr <- pen.yr %>% 
  select(id=WimsId, paste0("Capt",14:20))%>%
  pivot_longer(-id)%>%
  mutate(id=paste0("CN", str_sub(id,2,-1)),
         year=paste0("20", str_sub(name,-2,-1))%>%as.numeric)

pen.yr.first <- pen.yr %>% 
  drop_na(value)%>%
  select(id,year)%>%
  rbind(surv%>%
          filter(id%in%pen.yr$id)%>%
          group_by(id)%>%
          summarise(year=min(year(start))))%>%
  group_by(id)%>%
  summarise(min.year=min(year))

pen.yr <- pen.yr%>%
  left_join(pen.yr.first)%>%
  mutate(value=case_when(is.na(value)~"W",
                         TRUE~value))%>%
  filter(year>=min.year)%>%
  select(id, loc=value, year)%>%
  distinct(id,year,.keep_all = TRUE)%>%
  left_join(pen%>%select(year,penned,free), by="year")

write_csv(pen.yr, here::here("data","whopenned.csv"))

#CN343S 2019 both in and out
#x

##keep only 2014 onwards
surv <- surv %>% 
  filter(year>2013)

##bind with surv
surv <- surv %>% 
  left_join(pen.yr%>%select(id,loc,year)%>%distinct(), by=c("id","year"))%>%
  left_join(pen.yr%>%mutate(year=year)%>%select(year,id,free)%>%distinct(), by=c("id","year"))%>%
  left_join(pen.yr%>%mutate(year=year+1)%>%select(year,id,loc_prev=loc,free_prev=free)%>%distinct(), by=c("id","year"))


##find errors, join is missing some animals
surv%>%
  filter(is.na(loc))%>%
  group_by(id)%>%
  count

###if NA, drop
surv <- surv %>% 
  drop_na(loc)



ggplot(surv%>%distinct(id,year,loc), aes(x=year,y=id,color=loc,group=id))+
  geom_path()+
  geom_point()


##fix when pen yr assigned to few months before in pen
surv <- surv%>%
  mutate(loc=case_when(start%within%free_prev  ~ loc_prev,
                       TRUE~loc))



##get age and remove calves and yearlings
surv <-surv %>% 
  mutate(age=(end-birth)%>%time_length("year"))%>%
  filter(age>2)

length(unique(surv$id))


##get to pen-only times
surv <- surv%>%filter(start>=ymd("2014-03-27"))

##make sure it all adds up
##missing two penned animals in 2020
##rest of discrepancy is juveniles. Are they collared?
##wild won't add up as I forced wild when not penned
summary <- surv%>%
  filter(season%in%"Pen")%>%
  group_by(year,loc)%>%
  summarize(n_telem_surv=n_distinct(id))
summary%>%filter(year>2013)%>%arrange(loc)%>%
  left_join(read_csv(here::here("data","Klinse-za capture and demographics_210406_excerpt.csv"))%>%
              select(id=WimsId, paste0("Capt",14:20))%>%
              pivot_longer(-id)%>%
              mutate(id=paste0("CN", str_sub(id,2,-1)),
                     year=paste0("20", str_sub(name,-2,-1))%>%as.numeric)%>%
              select(id, loc=value, year)%>%
              drop_na(loc)%>%
              group_by(year,loc)%>%
              summarize(n_pen.yr=n_distinct(id)))


  
##flag
x
##pen-2016 should be 14 not 13



#2019 11 not 12
#2020 13 not 11

ggplot(surv, aes(x=end,y=id,color=season))+
  geom_point(size=0.01)+
  facet_wrap(vars(season, loc))

##fix year for records just prior to penning
surv <- surv%>%
  left_join(pen%>%
              mutate(prepen=(ymd(paste0(year(ymd(into)),"-01-01"))-1) %--% ymd(into))%>%select(year, prepen))%>%
  mutate(year=case_when(end %within% prepen~year-1, TRUE~year))
  

###calculate time since penned
##idea from here https://stackoverflow.com/questions/26553638/calculate-elapsed-time-since-last-event/26555648
set.seed(12345)
id <- c(rep(1, 9), rep(2, 9), rep(3, 9))
time <- c(seq(from = 0, to = 96, by = 12),
          seq(from = 0, to = 80, by = 10),
          seq(from = 0, to = 112, by = 14))
random <- runif(n = 27)
event <- rep(100, 27)

df <- data.frame(cbind(id, time, event, random))
df$event <- ifelse(df$random < 0.55, 0, df$event)
df <- subset(df, select = -c(random))
#df$event <- ifelse(df$time == 0, 100, df$event)
head(df)
df %>%
  mutate(tmpG = cumsum(c(FALSE, as.logical(diff(event==100)))))%>%
  group_by(id) %>%
  mutate(tmp_a = c(0, diff(time)) * !event) %>%
  group_by(tmpG) %>%
  mutate(tae = cumsum(tmp_a)) %>%
  ungroup() %>%
  select(-c(tmp_a, -tmpG))





surv <- surv%>%
  mutate(inpen=case_when(season%in%"Pen" & loc%in%"P"~1,
                          TRUE~0))%>%
  mutate(tmpG = paste(cumsum(c(FALSE, as.logical(diff(inpen==1)))),id))%>%
  group_by(id) %>%
  mutate(tmp_a = case_when(first(loc)%in%"W" & tmpG==min(tmpG)~0,
                           TRUE~dur * !inpen)) %>%
  group_by(tmpG) %>%
  mutate(timesince = cumsum(tmp_a)) %>%
  ungroup()%>%
  select(id,start,end,dead,dur,year,season,loc,loc_prev,inpen,timesince)




##sum up across individuals and covariates of interest
surv.pen <- surv%>%
  group_by(id,year,season,loc)%>%
  summarise(dur=sum(dur),
            dead=max(dead),
            sincepen=max(timesince)/365)%>%
  rename(PenSeason=season,
         PennedThatYear=loc)


surv.pen%>%
  group_by(PennedThatYear,PenSeason)%>%
  summarize(n=n_distinct(id))


hist(surv.pen$dur)

##model survival

##COXPH
m1 <- coxph(Surv(dur, dead)~ 1 +
              frailty(id), data = surv.pen)
m2 <- coxph(Surv(dur, dead)~ PenSeason +
              frailty(id), data = surv.pen)
m3 <- coxph(Surv(dur, dead)~ PenSeason + PennedThatYear +
              frailty(id), data = surv.pen)
m4 <- coxph(Surv(dur, dead)~ PenSeason*PennedThatYear +
              frailty(id), data = surv.pen)
m5 <- coxph(Surv(dur, dead)~ PenSeason+PennedThatYear + sincepen +
              frailty(id), data = surv.pen)

model.sel(m1,m2,m3,m4,m5,  rank="AICc") ##m4 doesnt seem to converge
cox.zph(m3)
summary(m3)
ggforest(coxph(Surv(dur, dead)~ PenSeason+PennedThatYear, data = surv.pen), data=surv.pen)
cox.zph(m5)
summary(m5)
ggforest(coxph(Surv(dur, dead)~ PenSeason+PennedThatYear + sincepen, data = surv.pen), data=surv.pen)

ggforest(coxph(Surv(dur, dead)~ PenSeason+ sincepen, data = surv.pen%>%filter(sincepen>0)), data=surv.pen%>%filter(sincepen>0))
m5 <- coxph(Surv(dur, dead)~ PenSeason + sincepen +
              frailty(id), data = surv.pen%>%filter(sincepen>0))

m1 <- survfit(Surv(dur, dead)~ PenSeason, data = surv.pen)
m2 <- survfit(Surv(dur, dead)~ PennedThatYear, data = surv.pen)
m3 <- survfit(Surv(dur, dead)~ PenSeason + PennedThatYear, data = surv.pen)


ggsurvplot(m3, data = surv.pen, ylim = c(0.8, 1))

plot.dat <- tibble(
  group=rownames(summary(m3)$table),
  surv=summary(m3, times=240, extend=TRUE)$surv,
  se=summary(m3, times=240, extend=TRUE)$std.err,
  n=summary(m3, times=240, extend=TRUE)$n,
  morts=summary(m3, times=240, extend=TRUE)$n.event
)

ggplot(data=plot.dat, aes(x=group,y=surv, ymin=surv-se, ymax=surv+se, label=paste0("n=",n, " (morts=",morts,")")))+
  geom_point()+
  geom_text(nudge_x=0.2, nudge_y=-0.03, size=3)+
  geom_pointrange()+
  theme_ipsum()+
  coord_flip()+
  labs(y="Survival rate (8 month)", x="Group")+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))











##COXPH INDIVIDUAL
surv.boot <- data.frame()
for(i in 1:1000){
mY <- survfit(Surv(dur, dead)~ PenSeason,
               data = surv.pen%>%dplyr::filter(PennedThatYear=="P")%>%dplyr::group_by(id,PenSeason)%>%dplyr::sample_frac(1, replace=TRUE))

mN <- survfit(Surv(dur, dead)~ PenSeason,
               data = surv.pen%>%dplyr::filter(PennedThatYear=="W")%>%dplyr::group_by(id,PenSeason)%>%dplyr::sample_frac(1, replace=TRUE))


surv.boot <- rbind(surv.boot,
                   rbind(
                     tibble(Penned="Y",
                            season=rownames(summary(mY)$table),
                            surv=summary(mY, times=240, extend=TRUE)$surv,
                            se=summary(mY, times=240, extend=TRUE)$std.err,
                            n=summary(mY, times=240, extend=TRUE)$n,
                            morts=summary(mY, times=240, extend=TRUE)$n.event
                            ),
                     tibble(Penned="N",
                            season=rownames(summary(mN)$table),
                            surv=summary(mN, times=240, extend=TRUE)$surv,
                            se=summary(mN, times=240, extend=TRUE)$std.err,
                            n=summary(mN, times=240, extend=TRUE)$n,
                            morts=summary(mN, times=240, extend=TRUE)$n.event)
                     )%>%
                     dplyr::mutate(iter=i)
                   )

}

##Summary stats
sum.stat <-surv.boot%>%
  drop_na(surv)%>%
  group_by(season,Penned)%>%
  summarise(median=median(surv),
            lower=quantile(surv,0.05),
            upper=quantile(surv,0.95),
            se=sd(surv))
sum.stat

####differences between survival for the out of pen period
surv.boot%>%
  filter(season%in%"PenSeason=Free")%>%
  select(Penned, iter, surv)%>%
  pivot_wider(names_from=Penned, values_from=surv)%>%
  drop_na()%>%
  mutate(dif=Y-N)%>%
  summarise(median=median(dif),
            lower=quantile(dif,0.05),
            upper=quantile(dif,0.95))


####differences between survival for the in pen period
surv.boot%>%
  filter(season%in%"PenSeason=Pen")%>%
  select(Penned, iter, surv)%>%
  pivot_wider(names_from=Penned, values_from=surv)%>%
  drop_na()%>%
  mutate(dif=Y-N)%>%
  summarise(median=median(dif),
            lower=quantile(dif,0.05),
            upper=quantile(dif,0.95))




###PLOT
ggplot(data=surv.boot, aes(x=season,y=surv,color=Penned))+
  geom_point(alpha=0.02, position=position_dodge(width = .5))+
  theme_ipsum()+
  labs(y="Survival (monthly)", x="Season")+
  labs(y="Survival rate (8 month)", x="Group")+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))



  
#####################
####CALF SURVIVAL####
#####################
 
##load data
calf <- read_csv(here::here("data","Klinse-za capture and demographics_201215.csv"))%>%
    drop_na(PenYr)
  
#clean up parturition dates
calf_born <- calf %>% 
  mutate(born=ymd(Born_Date),
         born_plot=born)%>%
  filter(!born%in%ymd("2017-04-16"))%>%
  select(PenYr, born, born_plot)

year(calf_born$born_plot)<-2020

##plot
ggplot(calf_born, aes(x=born_plot))+
  geom_histogram()+
  facet_wrap(vars(PenYr))+
  theme_ipsum()

ggplot(calf_born, aes(x=born_plot))+
  geom_histogram()+
  theme_ipsum()

ggplot(calf_born, aes(x=born_plot))+
  geom_density()+
  theme_ipsum()

ggplot(calf_born, aes(x=born_plot, fill=as.factor(PenYr)))+
  geom_density(alpha=0.2)+
  theme_ipsum()

##clean survival data
#calf_surv$`Alive/Dead`%>%unique()
calf_surv <- calf%>%
  filter(CalfAgeat1Yr>0 & !is.na(CalfAgeat1Yr))%>%
  mutate(start=ymd(CaptDate),
         end=start+CalfAgeat1Yr,
         dur=(end-start)%>%as.numeric(),
         ReleaseAge=as.numeric(ReleaseAge))%>%
  mutate(dead=case_when(dur<365~1,
                        TRUE~0),
         dead2=case_when(`Alive/Dead`%in%c("Dead")~1,
                        `Alive/Dead`%in%c("Alive","Unk.","Unk")~0),
         dur=case_when(dur<=365~dur,
                       dur>365~365))%>%
  select(id=WIMS_ID, Sex, start, end, dead, dur, ReleaseAge)



##model survival
m1 <- survfit(Surv(dur, dead)~ 1, data = calf_surv)
m2 <- survfit(Surv(dur, dead)~ Sex, data = calf_surv)

ggsurvplot(m1, data = calf_surv, xlim=c(0,365))
ggsurvplot(m2, data = calf_surv, xlim=c(0,365))




##COXPH
m1 <- coxph(Surv(dur, dead)~ 1 +
              frailty(id), data = calf_surv%>%drop_na(ReleaseAge))
m2 <- coxph(Surv(dur, dead)~ Sex +
              frailty(id), data = calf_surv%>%drop_na(ReleaseAge))
m3 <- coxph(Surv(dur, dead)~ ReleaseAge +
              frailty(id), data = calf_surv%>%drop_na(ReleaseAge))
m4 <- coxph(Surv(dur, dead)~ ReleaseAge + Sex +
              frailty(id), data = calf_surv%>%drop_na(ReleaseAge))

model.sel(m1,m2,m3,m4,  rank="AICc")
cox.zph(m4)
summary(m4)
ggforest(coxph(Surv(dur, dead)~ ReleaseAge, data = calf_surv%>%drop_na(ReleaseAge)))

