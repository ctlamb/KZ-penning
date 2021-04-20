Klinse-Za Penning Demography + Range Recolonization
================
Clayton Lamb
20 April, 2021

## Load Packages & Data

``` r
library(survival)
library(here)
library(lubridate)
library(MuMIn)
library(survminer)
library(hrbrthemes)
library(tidyverse)
library(tidylog)
source(here::here("functions","seasonalsurvival_fn.r"))

#ADULT FEMALE SURVIVAL
df <- read_csv(here::here("data","CollarFates4CL-RS_201230.csv"))%>%
  filter(herd%in%c("K","S"),
         Sex%in%"F")
  
#ADULT FEMALE LOCATIOn (IN/OUT OF PEN)
pen.yr <- read_csv(here::here("data","Klinse-za capture and demographics_210406_excerpt.csv"))

##PENNING DATES
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

##CALF SURVIVAL
calf <- read_csv(here::here("data","Klinse-za capture and demographics_210406_calf.csv"))%>%
    drop_na(PenYr)
```

## Clean and Prep

``` r
#MIN AGE VIA BIRTHDATE (MAYBE REDUNDANT?)
df <- df %>% 
  mutate(rstart_dat=dmy(rstart_dat),
         rend_dat=dmy(rend_dat),
    birth=case_when(agecl%in%"C" & month(rstart_dat) %in% 5:12~ paste(year(rstart_dat),"-06-01",sep="-")%>%ymd,
                    agecl%in%"C" & month(rstart_dat) %in% 1:4~ paste(year(rstart_dat)-1,"-06-01",sep="-")%>%ymd,
                    agecl%in%"Y" & month(rstart_dat) %in% 5:12~ paste(year(rstart_dat)-1,"-06-01",sep="-")%>%ymd,
                    agecl%in%"Y" & month(rstart_dat) %in% 1:4~ paste(year(rstart_dat)-2,"-06-01",sep="-")%>%ymd,
                    agecl%in%"A"~ paste(year(rstart_dat)-3,"-06-01",sep="-")%>%ymd))

#CLASSIFY AS DEAD OR ALIVE
#one myop but happens in 2002. Watch for myop pen 2020 (don't see here yet)
unique(df$mortcaus)
```

    ## [1] "DEAD" "FADE" "MYOP"

``` r
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
  mutate(dur=(end-start)%>%as.numeric)

df%>%
  filter(dur<1) ## all good



## Get into survival format
##make daily survival
surv <- stretch_survival_data(df, '1 day')

nrow(surv)
```

    ## [1] 67065

``` r
sum(df$dur) ##should be same
```

    ## [1] 67065

## Assign penned animals and seasons

``` r
####ASSIGN SEASONS
##Assign penned time and drop data before this period
surv <- surv%>%
  mutate(dur=(end-start)%>%as.numeric(),
         year=year(end),
         season=case_when(start %within% as.list(pen$penned)~"Pen", 
                          start %within% as.list(pen$free)~"Free",
                          TRUE~NA_character_))%>%
  drop_na(season)
           

##clean up id's (get rid of trailing letter denoting period)
surv <- surv %>% 
  mutate(id=str_sub(id,1,-2))

##plot
ggplot(surv, aes(x=end,y=id,color=season))+
  geom_point(size=0.01)+
  facet_wrap(vars(season))+
  theme_ipsum()+
  labs(title="Season")
```

![](README_files/figure-gfm/assign%20penned-1.png)<!-- -->

``` r
####ASSIGN LOCATIONS IN/OUT
##decide who was in the pen
pen.yr <- pen.yr %>% 
  select(id=WimsId, paste0("Capt",14:20))%>%
  pivot_longer(-id)%>%
  mutate(id=paste0("CN", str_sub(id,2,-1)),
         year=paste0("20", str_sub(name,-2,-1))%>%as.numeric)

##figure out first year animals were penned
pen.yr.first <- pen.yr %>% 
  drop_na(value)%>%
  select(id,year)%>%
  rbind(surv%>%
          filter(id%in%pen.yr$id)%>%
          group_by(id)%>%
          summarise(year=min(year(start))))%>%
  group_by(id)%>%
  summarise(min.year=min(year))

##join pen locs and first yeat together
pen.yr <- pen.yr%>%
  left_join(pen.yr.first)%>%
  mutate(value=case_when(is.na(value)~"W",
                         TRUE~value))%>%
  filter(year>=min.year)%>%
  select(id, loc=value, year)%>%
  distinct(id,year,.keep_all = TRUE)%>%
  left_join(pen%>%select(year,penned,free), by="year")

##export for RSF analysis
write_csv(pen.yr, here::here("data","whopenned.csv"))

##bind with surv
surv <- surv %>% 
  left_join(pen.yr%>%select(id,loc,year)%>%distinct(), by=c("id","year"))%>% ##where each animal was per year
  left_join(pen.yr%>%mutate(year=year)%>%select(year,id,free)%>%distinct(), by=c("id","year"))%>% ##add in "free" dates
  left_join(pen.yr%>%mutate(year=year+1)%>%select(year,id,loc_prev=loc,free_prev=free)%>%distinct(), by=c("id","year")) ##add in lag of where they were year before


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
```

![](README_files/figure-gfm/assign%20penned-2.png)<!-- -->

``` r
##fix when pen yr assigned to few months before in pen
surv <- surv%>%
  mutate(loc=case_when(start%within%free_prev  ~ loc_prev,
                       TRUE~loc))



##get age and remove calves and yearlings
surv <-surv %>% 
  mutate(age=(end-birth)%>%time_length("year"))%>%
  filter(age>2)

length(unique(surv$id))
```

    ## [1] 43

``` r
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



ggplot(surv, aes(x=end,y=id,color=season))+
  geom_point(size=0.01)+
  facet_wrap(vars(season, loc))+
  theme_ipsum()
```

![](README_files/figure-gfm/assign%20penned-3.png)<!-- -->

``` r
##fix year for records just prior to penning
surv <- surv%>%
  left_join(pen%>%
              mutate(prepen=(ymd(paste0(year(ymd(into)),"-01-01"))-1) %--% ymd(into))%>%select(year, prepen))%>%
  mutate(year=case_when(end %within% prepen~year-1, TRUE~year))
```

## Calculate time since penned

``` r
##idea from here https://stackoverflow.com/questions/26553638/calculate-elapsed-time-since-last-event/26555648
# set.seed(12345)
# id <- c(rep(1, 9), rep(2, 9), rep(3, 9))
# time <- c(seq(from = 0, to = 96, by = 12),
#           seq(from = 0, to = 80, by = 10),
#           seq(from = 0, to = 112, by = 14))
# random <- runif(n = 27)
# event <- rep(100, 27)
# 
# df <- data.frame(cbind(id, time, event, random))
# df$event <- ifelse(df$random < 0.55, 0, df$event)
# df <- subset(df, select = -c(random))
# #df$event <- ifelse(df$time == 0, 100, df$event)
# head(df)
# df %>%
#   mutate(tmpG = cumsum(c(FALSE, as.logical(diff(event==100)))))%>%
#   group_by(id) %>%
#   mutate(tmp_a = c(0, diff(time)) * !event) %>%
#   group_by(tmpG) %>%
#   mutate(tae = cumsum(tmp_a)) %>%
#   ungroup() %>%
#   select(-c(tmp_a, -tmpG))


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


##how many individuals per group?
surv.pen%>%
  group_by(PennedThatYear,PenSeason)%>%
  summarize(n=n_distinct(id))


##distribution of durations monitored
hist(surv.pen$dur)
```

![](README_files/figure-gfm/time%20since%20penned-1.png)<!-- -->

# Model survival

## CoxPh

``` r
##COXPH survival
m1 <- coxph(Surv(dur, dead)~ 1 +
              frailty(id), data = surv.pen)
m2 <- coxph(Surv(dur, dead)~ PenSeason +
              frailty(id), data = surv.pen)
m3 <- coxph(Surv(dur, dead)~ PenSeason + PennedThatYear +
              frailty(id), data = surv.pen)
# m4 <- coxph(Surv(dur, dead)~ PenSeason*PennedThatYear +
#               frailty(id), data = surv.pen)  ##doesn't converge

model.sel(m1,m2,m3,  rank="AICc")
```

    ## Model selection table 
    ##    (Int) frl(id) PnS PTY family df  logLik  AICc delta weight
    ## m3     +       +   +   +   (NA)  2 -88.378 181.6  0.00  0.850
    ## m1     +       +           (NA)  4 -85.888 185.9  4.37  0.096
    ## m2     +       +   +       (NA)  5 -84.301 187.1  5.50  0.054
    ## Models ranked by AICc(x)

``` r
cox.zph(m3)
```

    ##                chisq df    p
    ## PenSeason      0.662  1 0.42
    ## PennedThatYear 2.631  1 0.10
    ## GLOBAL         3.290  2 0.19

``` r
summary(m3)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ PenSeason + PennedThatYear + 
    ##     frailty(id), data = surv.pen)
    ## 
    ##   n= 275, number of events= 18 
    ## 
    ##                 coef    se(coef) se2    Chisq DF p   
    ## PenSeasonPen    -1.2866 0.7983   0.7983 2.60  1  0.11
    ## PennedThatYearW  0.5462 0.4751   0.4751 1.32  1  0.25
    ## frailty(id)                             0.00  0  0.84
    ## 
    ##                 exp(coef) exp(-coef) lower .95 upper .95
    ## PenSeasonPen       0.2762     3.6206   0.05777     1.320
    ## PennedThatYearW    1.7266     0.5792   0.68046     4.381
    ## 
    ## Iterations: 5 outer, 19 Newton-Raphson
    ##      Variance of random effect= 5e-05   I-likelihood = -88.4 
    ## Degrees of freedom for terms= 1 1 0 
    ## Concordance= 0.765  (se = 0.066 )
    ## Likelihood ratio test= 4.57  on 2 df,   p=0.1

``` r
ggforest(coxph(Surv(dur, dead)~ PenSeason+PennedThatYear, data = surv.pen), data=surv.pen)
```

![](README_files/figure-gfm/model%20coxph-1.png)<!-- -->

``` r
##COXPH time since
# m5 <- coxph(Surv(dur, dead)~ PenSeason + sincepen +
#               frailty(id), data = surv.pen%>%filter(sincepen>0))
# m6 <- coxph(Surv(dur, dead)~ sincepen +
#               frailty(id), data = surv.pen%>%filter(sincepen>0))
# cox.zph(m5)
# summary(m5)
# cox.zph(m6)
# summary(m6)
# ggforest(coxph(Surv(dur, dead)~ PenSeason + sincepen, data = surv.pen%>%filter(sincepen>0)), data=surv.pen%>%filter(sincepen>0))
```

## Kap Meir

``` r
m1 <- survfit(Surv(dur, dead)~ PenSeason, data = surv.pen)
m2 <- survfit(Surv(dur, dead)~ PennedThatYear, data = surv.pen)
m3 <- survfit(Surv(dur, dead)~ PenSeason + PennedThatYear, data = surv.pen)


ggsurvplot(m3, data = surv.pen, ylim = c(0.8, 1))
```

![](README_files/figure-gfm/model%20km-1.png)<!-- -->

``` r
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
```

![](README_files/figure-gfm/model%20km-2.png)<!-- -->

## CoxPh bootstrap

``` r
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
```

![](README_files/figure-gfm/model%20coxph%20boot-1.png)<!-- -->

## Calf survival

``` r
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
```

![](README_files/figure-gfm/model%20calf-1.png)<!-- -->

``` r
ggplot(calf_born, aes(x=born_plot))+
  geom_histogram()+
  theme_ipsum()
```

![](README_files/figure-gfm/model%20calf-2.png)<!-- -->

``` r
ggplot(calf_born, aes(x=born_plot))+
  geom_density()+
  theme_ipsum()
```

![](README_files/figure-gfm/model%20calf-3.png)<!-- -->

``` r
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
```

![](README_files/figure-gfm/model%20calf-4.png)<!-- -->

``` r
ggsurvplot(m2, data = calf_surv, xlim=c(0,365))
```

![](README_files/figure-gfm/model%20calf-5.png)<!-- -->

``` r
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
```

    ## Model selection table 
    ##    (Int) frl(id) Sex     RlA family df  logLik AICc delta weight
    ## m4     +       +   + -0.4270   (NA) 12 -21.338  6.7  0.00  0.907
    ## m3     +       +     -0.4994   (NA) 20 -17.680 11.2  4.56  0.093
    ## m1     +       +               (NA)  0 -31.207 62.4 55.75  0.000
    ## m2     +       +   +           (NA)  1 -30.737 64.1 57.48  0.000
    ## Models ranked by AICc(x)

``` r
cox.zph(m4)
```

    ##            chisq    df      p
    ## ReleaseAge  9.50  0.78 0.0013
    ## Sex         0.28  0.78 0.4945
    ## GLOBAL      9.50 12.60 0.7053

``` r
summary(m4)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ ReleaseAge + Sex + frailty(id), 
    ##     data = calf_surv %>% drop_na(ReleaseAge))
    ## 
    ##   n= 53, number of events= 8 
    ## 
    ##             coef    se(coef) se2    Chisq DF    p    
    ## ReleaseAge  -0.4270 0.2559   0.2263  2.78  1.00 0.095
    ## SexM        -0.6735 0.8534   0.7520  0.62  1.00 0.430
    ## frailty(id)                         11.26 11.04 0.420
    ## 
    ##            exp(coef) exp(-coef) lower .95 upper .95
    ## ReleaseAge    0.6525      1.533   0.39513     1.077
    ## SexM          0.5099      1.961   0.09575     2.716
    ## 
    ## Iterations: 6 outer, 26 Newton-Raphson
    ##      Variance of random effect= 1.742685   I-likelihood = -29.4 
    ## Degrees of freedom for terms=  0.8  0.8 11.0 
    ## Concordance= 0.912  (se = 0.046 )
    ## Likelihood ratio test= 19.74  on 12.6 df,   p=0.09

``` r
ggforest(coxph(Surv(dur, dead)~ ReleaseAge, data = calf_surv%>%drop_na(ReleaseAge)))
```

![](README_files/figure-gfm/model%20calf-6.png)<!-- -->
