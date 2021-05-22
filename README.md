Klinse-Za Penning Demography
================
Clayton Lamb
22 May, 2021

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
library(infer)
library(readxl)
library(lme4)
library(chisq.posthoc.test)
options(scipen=999)
source(here::here("functions","seasonalsurvival_fn.r"))

#ADULT FEMALE SURVIVAL
df <- read_excel(here::here("data","demography", "AllHerdsStartEndDates_210421.xlsx"))%>%
  filter(herd%in%c("K","S"),
         Sex%in%"F")
  
#ADULT FEMALE LOCATIOn (IN/OUT OF PEN)
ad.cap <- read_excel(here::here("data","demography", "Klinse-za capture and demographics_210521.xlsx"),
                     sheet="Adult-Capt")%>%
      drop_na(CalendarYr)

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
  2021, "2021-03-15", "2021-08-05")%>%
  mutate(free_end=(ymd(dplyr::lead(into,1,order_by=year))-1)%>%as.character,
      penned=ymd(into)%--%ymd(out),
      free=ymd(out)%--%ymd(free_end))

##CALF SURVIVAL
calf <- read_excel(here::here("data","demography", "Klinse-za capture and demographics_210521.xlsx"),
                     sheet="Calf-Capt")%>%
    drop_na("Mom_ID")

##MORTS
# mort <- read_excel(here::here("data","demography", "Klinse-za capture and demographics_210521.xlsx"),
#                      sheet="Mortalities")

##PARTURITION
mm.part <- read_csv(here::here("data","demography", "MovementModelResults.csv"))%>%
  drop_na(Caribou_ID)

##CONDITION
cond <- read_csv(here::here("data","health", "KZONCP_Health_Year2_results_210422.csv"))
```

## Clean and Prep

``` r
#MIN AGE VIA BIRTHDATE (MAYBE REDUNDANT?)
df <- df %>% 
  mutate(rstart_dat=rstart_dat%>%date,
         rend_dat=rend_dat%>%date,
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

    ## [1] 75583

``` r
sum(df$dur) ##should be same
```

    ## [1] 75583

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

##drop penned 2021
surv <- surv%>%
  filter(start<ymd("2021-03-15"))
           

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
pen.yr <- ad.cap %>% 
  select(id=WimsId, paste0("Capt",14:21))%>%
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

    ## [1] 45

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
  left_join(read_excel(here::here("data","demography", "Klinse-za capture and demographics_210521.xlsx"),
                     sheet="Adult-Capt")%>%
      drop_na(CalendarYr)%>%
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


##times penned
surv <- surv%>%
  left_join(surv%>%
  filter(loc=="P")%>%
  distinct(id,year(end),loc)%>%
  group_by(id)%>%
  summarize(times.penned=n()))
  



##sum up across individuals and covariates of interest
surv.pen <- surv%>%
  group_by(id,year,season,loc, times.penned)%>%
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

## Assess if differences in body condition or age between penned vs wild

``` r
##prep body condition data
bod.cond <- pen.yr%>%
  select(id,loc,year)%>%
  left_join(cond%>%
              mutate(id=paste0("CN",str_sub(`WII Animal ID`,2,-1)))%>%
              filter(ymd(Date)%>%month %in%2:4)%>%
              select(id,year=Year, condition_clean, age_clean))%>%
    mutate(age_clean=case_when(age_clean%in%c("mature", "mature (6-7)","mature (7)","mature (4-7)","mature (4-5)")~"mature",
                             age_clean%in%c("old","old (8-9)","old (10-11)","old (8-11)")~"old",
                             age_clean%in%c("young","young (2-3)", "yearling")~"young",
                             age_clean%in%c("unk")~"unk"))%>%
  drop_na(condition_clean)


##plot
bod.cond%>%
  filter(condition_clean%in%c("poor", "fair", "good"))%>%
  group_by(loc, condition_clean)%>%
  summarise(n = n()) %>%
  group_by(loc)%>%
  mutate(freq = n / sum(n),
         condition_clean=factor(condition_clean, levels=c("poor", "fair", "good")))%>%
  ungroup()%>%
  ggplot(aes(x=condition_clean, y=freq, fill=loc))+ 
  geom_col(position = "dodge", width=0.4)+
  theme_ipsum()+
  theme(legend.title=element_text(face = "bold",size = rel(1.3)))+
    labs(y="Proportion", x="Body condition")+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))
```

![](README_files/figure-gfm/body-1.png)<!-- -->

``` r
##chi square
bod.matrix <- bod.cond%>%
  filter(condition_clean%in%c("poor", "fair", "good"))%>%
  mutate(condition_clean=factor(condition_clean, levels=c("poor", "fair", "good")))%>%
  group_by(loc, condition_clean)%>%
  summarise(n = n())%>%
  pivot_wider(names_from="condition_clean",values_from="n")%>%
  column_to_rownames(var="loc")%>%
  as.matrix()
  
chisq.test(bod.matrix)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  bod.matrix
    ## X-squared = 0.78713, df = 2, p-value = 0.6746

``` r
# Shows post-hoc pairwise comparisons using fdr method
chisq.posthoc.test(bod.matrix,
                   method = "bonferroni")




##prep and plot ages
bod.cond%>%
  filter(age_clean%in%c("young", "mature", "old"))%>%
  group_by(loc, age_clean)%>%
  summarise(n = n()) %>%
  group_by(loc)%>%
  mutate(freq = n / sum(n),
         age_clean=factor(age_clean, levels=c("young", "mature", "old")))%>%
  ungroup()%>%
  ggplot(aes(x=age_clean, y=freq, fill=loc))+ 
  geom_col(position = "dodge", width=0.4)+
  theme_ipsum()+
  theme(legend.title=element_text(face = "bold",size = rel(1.3)))+
      labs(y="Proportion", x="Age class")+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))
```

![](README_files/figure-gfm/body-2.png)<!-- -->

``` r
##chi square
age.matrix <- bod.cond%>%
  mutate(age_clean=factor(age_clean, levels=c("young", "mature", "old")))%>%
   filter(age_clean%in%c("young", "mature", "old"))%>%
  group_by(loc, age_clean)%>%
  summarise(n = n())%>%
  pivot_wider(names_from="age_clean",values_from="n")%>%
  column_to_rownames(var="loc")%>%
  as.matrix()
  
chisq.test(age.matrix)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  age.matrix
    ## X-squared = 2.2917, df = 2, p-value = 0.3179

``` r
# Shows post-hoc pairwise comparisons using fdr method
chisq.posthoc.test(age.matrix,
                   method = "bonferroni")
```

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
    ##    (Int) frl(id) PnS PTY family df   logLik  AICc delta weight
    ## m3     +       +   +   +   (NA)  6 -102.117 225.6  0.00  0.919
    ## m1     +       +           (NA)  7 -102.791 230.5  4.91  0.079
    ## m2     +       +   +       (NA)  9 -100.263 237.6 12.07  0.002
    ## Models ranked by AICc(x)

``` r
cox.zph(m3)
```

    ##                 chisq   df    p
    ## PenSeason      0.0216 1.00 0.88
    ## PennedThatYear 1.6084 0.98 0.20
    ## GLOBAL         1.6276 6.74 0.97

``` r
ggcoxzph(cox.zph(m3))
```

![](README_files/figure-gfm/model%20coxph-1.png)<!-- -->

``` r
ggcoxdiagnostics(m3, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
```

![](README_files/figure-gfm/model%20coxph-2.png)<!-- -->

``` r
summary(m3)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ PenSeason + PennedThatYear + 
    ##     frailty(id), data = surv.pen)
    ## 
    ##   n= 299, number of events= 21 
    ## 
    ##                 coef    se(coef) se2    Chisq DF   p    
    ## PenSeasonPen    -1.1011 0.6614   0.6613 2.77  1.00 0.096
    ## PennedThatYearW  0.7919 0.4553   0.4505 3.03  1.00 0.082
    ## frailty(id)                             5.14  4.76 0.370
    ## 
    ##                 exp(coef) exp(-coef) lower .95 upper .95
    ## PenSeasonPen       0.3325      3.008   0.09094     1.216
    ## PennedThatYearW    2.2077      0.453   0.90442     5.389
    ## 
    ## Iterations: 7 outer, 28 Newton-Raphson
    ##      Variance of random effect= 0.2550131   I-likelihood = -106.9 
    ## Degrees of freedom for terms= 1.0 1.0 4.8 
    ## Concordance= 0.784  (se = 0.058 )
    ## Likelihood ratio test= 16.18  on 6.74 df,   p=0.02

``` r
ggforest(m3, data=surv.pen)
```

![](README_files/figure-gfm/model%20coxph-3.png)<!-- -->

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

m7 <- coxph(Surv(dur, dead)~ PenSeason + PennedThatYear + times.penned +
              frailty(id), data = surv.pen)
cox.zph(m7)
```

    ##                 chisq df     p
    ## PenSeason      0.0211  1 0.884
    ## PennedThatYear 0.2001  1 0.655
    ## times.penned   4.9687  1 0.026
    ## GLOBAL         5.0587  3 0.168

``` r
ggcoxzph(cox.zph(m7))
```

![](README_files/figure-gfm/model%20coxph-4.png)<!-- -->

``` r
summary(m7)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ PenSeason + PennedThatYear + 
    ##     times.penned + frailty(id), data = surv.pen)
    ## 
    ##   n= 279, number of events= 16 
    ##    (20 observations deleted due to missingness)
    ## 
    ##                 coef    se(coef) se2    Chisq DF p      
    ## PenSeasonPen    -1.0576 0.6697   0.6697  2.49 1  0.11000
    ## PennedThatYearW  0.5016 0.5024   0.5024  1.00 1  0.32000
    ## times.penned    -0.8544 0.2302   0.2302 13.78 1  0.00021
    ## frailty(id)                              0.00 0  0.95000
    ## 
    ##                 exp(coef) exp(-coef) lower .95 upper .95
    ## PenSeasonPen       0.3473     2.8794   0.09347    1.2904
    ## PennedThatYearW    1.6514     0.6056   0.61691    4.4206
    ## times.penned       0.4255     2.3501   0.27101    0.6681
    ## 
    ## Iterations: 6 outer, 23 Newton-Raphson
    ##      Variance of random effect= 0.0000005   I-likelihood = -72.1 
    ## Degrees of freedom for terms= 1 1 1 0 
    ## Concordance= 0.862  (se = 0.037 )
    ## Likelihood ratio test= 26.68  on 3 df,   p=0.000007

``` r
ggforest(m7, data=surv.pen)
```

![](README_files/figure-gfm/model%20coxph-5.png)<!-- -->

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

## KM bootstrap

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

![](README_files/figure-gfm/model%20KM%20boot-1.png)<!-- -->

## Calf survival

``` r
#clean up parturition dates
calf_born <- calf %>% 
  mutate(born=ymd(BornDat),
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
    theme_ipsum()+
  labs(y="Count", x="Date", title="Parturition date in pen")+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))
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
  filter(!WimsId%in%"NA")%>%
  mutate(Born_Date=ymd(BornDat),
         LastObsAsCalf=ymd(ColLastDat))%>%
  select(WimsId,Born_Date,DeadCalf,ReleaseAge=WksAtRel, LastObsAsCalf,Notes.calf=Notes, Sex)%>%
  arrange(DeadCalf)%>%
  mutate(dur=case_when(as.numeric(LastObsAsCalf-Born_Date)<=366 ~as.numeric(LastObsAsCalf-Born_Date),
                       as.numeric(LastObsAsCalf-Born_Date)>366~365),
         year=year(Born_Date))%>%
  select(id=WimsId, Sex, year,start=Born_Date, end=LastObsAsCalf, dead=DeadCalf, dur, ReleaseAge)

cor(calf_surv$year, calf_surv$ReleaseAge, use="complete.obs")
```

    ## [1] 0.6752328

``` r
##model survival
m1 <- survfit(Surv(dur, dead)~ 1, data = calf_surv%>%filter(dur>0))
m2 <- survfit(Surv(dur, dead)~ Sex, data = calf_surv%>%filter(dur>0))

ggsurvplot(m1, data = calf_surv, xlim=c(0,365))
```

![](README_files/figure-gfm/model%20calf-4.png)<!-- -->

``` r
ggsurvplot(m2, data = calf_surv, xlim=c(0,365))
```

![](README_files/figure-gfm/model%20calf-5.png)<!-- -->

``` r
##what is annual calf survival?
summary(survfit(Surv(dur, dead)~ 1, data = calf_surv%>%filter(dur>0)), times = 365)
```

    ## Call: survfit(formula = Surv(dur, dead) ~ 1, data = calf_surv %>% filter(dur > 
    ##     0))
    ## 
    ##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
    ##   365     37      10    0.841  0.0462        0.756        0.937

``` r
##COXPH
m1 <- coxph(Surv(dur, dead)~ 1, data = calf_surv%>%drop_na(ReleaseAge))
m2 <- coxph(Surv(dur, dead)~ Sex, data = calf_surv%>%drop_na(ReleaseAge))
m3 <- coxph(Surv(dur, dead)~ ReleaseAge, data = calf_surv%>%drop_na(ReleaseAge))
m4 <- coxph(Surv(dur, dead)~ ReleaseAge + Sex, data = calf_surv%>%drop_na(ReleaseAge))
m5 <- coxph(Surv(dur, dead)~ Sex + year, data = calf_surv%>%drop_na(ReleaseAge))
m6 <- coxph(Surv(dur, dead)~   year, data = calf_surv%>%drop_na(ReleaseAge))

model.sel(m1,m2,m3,m4,m5,m6,  rank="AICc")
```

    ## Model selection table 
    ##    (Intrc) Sex   RlsAg    year family      class df  logLik AICc delta weight
    ## m3       +     -0.6430           (NA)      coxph  1 -28.322 59.3  0.00  0.639
    ## m4       +   + -0.5665           (NA)      coxph  2 -27.487 61.4  2.06  0.228
    ## m2       +   +                   (NA)      coxph  1 -31.017 64.7  5.39  0.043
    ## m1       +                       (NA) coxph.null  0 -32.383 64.8  5.45  0.042
    ## m6       +             -0.2621   (NA)      coxph  1 -31.362 65.4  6.08  0.031
    ## m5       +   +         -0.2518   (NA)      coxph  2 -30.095 66.6  7.28  0.017
    ## Models ranked by AICc(x)

``` r
cox.zph(m3)
```

    ##            chisq df      p
    ## ReleaseAge  10.7  1 0.0011
    ## GLOBAL      10.7  1 0.0011

``` r
cox.zph(m6)
```

    ##        chisq df    p
    ## year   0.613  1 0.43
    ## GLOBAL 0.613  1 0.43

``` r
ggcoxzph(m3%>%cox.zph)
```

![](README_files/figure-gfm/model%20calf-6.png)<!-- -->

``` r
ggcoxzph(m6%>%cox.zph)
```

![](README_files/figure-gfm/model%20calf-7.png)<!-- -->

``` r
summary(m3)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ ReleaseAge, data = calf_surv %>% 
    ##     drop_na(ReleaseAge))
    ## 
    ##   n= 63, number of events= 8 
    ## 
    ##               coef exp(coef) se(coef)      z Pr(>|z|)   
    ## ReleaseAge -0.6430    0.5257   0.2396 -2.684  0.00728 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##            exp(coef) exp(-coef) lower .95 upper .95
    ## ReleaseAge    0.5257      1.902    0.3287    0.8408
    ## 
    ## Concordance= 0.694  (se = 0.132 )
    ## Likelihood ratio test= 8.12  on 1 df,   p=0.004
    ## Wald test            = 7.2  on 1 df,   p=0.007
    ## Score (logrank) test = 7.04  on 1 df,   p=0.008

``` r
##proportional hazard violations. Mostly because calves died early (<100 days) early on, then later (>200 days) later on


#####massage into monthly survival so seasons/time can be accomadted



calf.surv.daily <- stretch_survival_data(calf_surv%>%
                                           drop_na(ReleaseAge)%>%
                                           mutate(end=start+dur,
                                                  herd="K",
                                                  birth=ReleaseAge), '1 day')%>%
  rename(ReleaseAge=birth)%>%
  mutate(dur=end-start)


calf.surv.monthly <- calf.surv.daily%>%
      group_by(id)%>%
  mutate(age=end-min(start),
         month=month(start))%>%
    group_by(id, ReleaseAge,month)%>%
    summarise(dur=sum(dur),
              dead=max(dead),
              age=mean(age))%>%
  left_join(calf_surv%>%select(id,year))


###try again
##COXPH
m1 <- coxph(Surv(dur, dead)~ 1, data =calf.surv.monthly)
m2 <- coxph(Surv(dur, dead)~ ReleaseAge, data =calf.surv.monthly)
m3 <- coxph(Surv(dur, dead)~ ReleaseAge + age, data =calf.surv.monthly)
m4 <- coxph(Surv(dur, dead)~ ReleaseAge*age, data =calf.surv.monthly)
m5 <- coxph(Surv(dur, dead)~ year+age, data =calf.surv.monthly)
m6 <- coxph(Surv(dur, dead)~ year*age, data =calf.surv.monthly)

model.sel(m1,m2,m3,m4,m5,m6,  rank="AICc")
```

    ## Model selection table 
    ##    (Int)     RlA       age  age:RlA     yer  age:yer family      class df  logLik  AICc delta weight
    ## m4     + -1.1840 -0.056750 0.005921                    (NA)      coxph  3 -41.210  94.4  0.00  0.601
    ## m2     + -0.7142                                       (NA)      coxph  1 -46.870  96.4  1.99  0.223
    ## m3     + -0.6182 -0.007008                             (NA)      coxph  2 -45.344  97.1  2.67  0.158
    ## m5     +         -0.008282          -0.2537            (NA)      coxph  2 -48.298 103.0  8.58  0.008
    ## m1     +                                               (NA) coxph.null  0 -51.547 103.1  8.67  0.008
    ## m6     +         -6.756000          -0.5663 0.003345   (NA)      coxph  3 -47.265 106.5 12.11  0.001
    ## Models ranked by AICc(x)

``` r
cox.zph(m2)
```

    ##              chisq df    p
    ## ReleaseAge 0.00209  1 0.96
    ## GLOBAL     0.00209  1 0.96

``` r
cox.zph(m4)
```

    ##                 chisq df    p
    ## ReleaseAge     0.0239  1 0.88
    ## age            0.5211  1 0.47
    ## ReleaseAge:age 0.2984  1 0.58
    ## GLOBAL         1.7948  3 0.62

``` r
summary(m4)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ ReleaseAge * age, data = calf.surv.monthly)
    ## 
    ##   n= 647, number of events= 8 
    ## 
    ##                     coef exp(coef)  se(coef)      z Pr(>|z|)    
    ## ReleaseAge     -1.184337  0.305949  0.316390 -3.743 0.000182 ***
    ## age            -0.056749  0.944831  0.015838 -3.583 0.000340 ***
    ## ReleaseAge:age  0.005921  1.005939  0.001657  3.573 0.000352 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                exp(coef) exp(-coef) lower .95 upper .95
    ## ReleaseAge        0.3059     3.2685    0.1646    0.5688
    ## age               0.9448     1.0584    0.9160    0.9746
    ## ReleaseAge:age    1.0059     0.9941    1.0027    1.0092
    ## 
    ## Concordance= 0.805  (se = 0.098 )
    ## Likelihood ratio test= 20.67  on 3 df,   p=0.0001
    ## Wald test            = 26.99  on 3 df,   p=0.000006
    ## Score (logrank) test = 31.69  on 3 df,   p=0.0000006

``` r
summary(m2)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ ReleaseAge, data = calf.surv.monthly)
    ## 
    ##   n= 647, number of events= 8 
    ## 
    ##               coef exp(coef) se(coef)      z Pr(>|z|)   
    ## ReleaseAge -0.7142    0.4896   0.2434 -2.934  0.00335 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##            exp(coef) exp(-coef) lower .95 upper .95
    ## ReleaseAge    0.4896      2.042    0.3038     0.789
    ## 
    ## Concordance= 0.668  (se = 0.144 )
    ## Likelihood ratio test= 9.35  on 1 df,   p=0.002
    ## Wald test            = 8.61  on 1 df,   p=0.003
    ## Score (logrank) test = 8.07  on 1 df,   p=0.004

``` r
calf.pred.boot <- data.frame()
##Boot for plot to show uncertainty
for(i in 1:1000){
  m3.i <- coxph(Surv(dur, dead)~ ReleaseAge + age, data = calf.surv.monthly%>%sample_frac(1,replace=TRUE))
 
  calf.pred.i <- expand.grid(ReleaseAge=5:12, age=seq(0,360,by=20), dur=30, dead=0, iter=i) 
  
  calf.pred.i$pred <- exp(-predict(m3.i, newdata=calf.pred.i, type="expected"))
  calf.pred.i$pred
  calf.pred.boot <- rbind(calf.pred.boot,calf.pred.i)
}


calf.pred.boot%>%
         drop_na()%>%
         mutate(pred=pred^12)%>%
         group_by(ReleaseAge,age)%>%
         summarize(median=median(pred))%>%
  ggplot()+
  geom_tile(aes(y=ReleaseAge, x=age,fill=median))
```

![](README_files/figure-gfm/model%20calf-8.png)<!-- -->

``` r
ggplot(data=calf.pred.boot%>%
         drop_na()%>%
         mutate(pred=pred^12)%>%
         group_by(ReleaseAge)%>%
         summarize(median=median(pred), lower=quantile(pred,0.025), upper=quantile(pred,0.975), se=sd(pred))%>%
         mutate(upper2=case_when(median+se>1~1, TRUE~median+se)), aes(x=ReleaseAge*7,y=median))+
         geom_path()+
  geom_ribbon(aes(ymin=median-se, ymax=upper2),alpha=0.5)+
    theme_ipsum()+
  labs(y="Survival rate (annual)", x="Release age (days)")+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))
```

![](README_files/figure-gfm/model%20calf-9.png)<!-- -->

``` r
ggplot(data=calf.pred.boot%>%
         drop_na()%>%
         mutate(pred=pred^12)%>%
         group_by(ReleaseAge,iter)%>%
         summarize(median=median(pred), lower=quantile(pred,0.025), upper=quantile(pred,0.975), se=sd(pred))%>%
         mutate(upper2=case_when(median+se>1~1, TRUE~median+se)), aes(x=ReleaseAge*7,y=median, group=iter))+
         geom_path(alpha=0.05)+
    theme_ipsum()+
  labs(y="Survival rate (annual)", x="Release age (days)")+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))
```

![](README_files/figure-gfm/model%20calf-10.png)<!-- -->

## Pregnancy

``` r
unique(ad.cap$`Preg if prog >1.2 or PSPB > 0.21`)
```

    ## [1] NA  "Y" "N"

``` r
preg <- ad.cap%>%
  filter(Age%in%"Ad")%>%
  mutate(preg=case_when(`Preg if prog >1.2 or PSPB > 0.21`%in% "Y"~1,
                        `Preg if prog >1.2 or PSPB > 0.21`%in% "N"~0,
                        TRUE~NA_real_),
         year=ymd(CapDate)%>%year(),
        class="Pregnancy")%>%
  drop_na(preg)%>%
  mutate(id=paste0("CN", str_sub(WimsId,2,-1)))%>%
  select(id,class, year,value=preg)%>%
  left_join(pen.yr%>%select(id, year, loc), by=c("id","year"))
```

## Parturition

``` r
mm.part.clean <- mm.part%>%
    mutate(`Calf (Predictions)`=case_when(`Validated?`%in% "False Negative"~1, 
         TRUE~`Calf (Predictions)`),##fix a couple wrong predictions
         id=paste0("CN",str_sub(Caribou_ID,2,-1)))%>%
  select(id,year=Year, part=`Calf (Predictions)`)%>%
  left_join(pen.yr%>%select(id, year, loc), by=c("id","year"))


##load calf data
part.cap <- calf%>%
  filter(!Mom_ID %in%"C360K")%>% ##died before calving
  mutate(CalfAgeat1Yr=case_when(is.na(CalfAgeat1Yr)~(-1), TRUE~CalfAgeat1Yr%>%as.numeric()))%>%
  mutate(part=case_when(CalfAgeat1Yr>0~1,TRUE~0),
         loc="P",
         year=paste0(20,PenYr))%>%
  select(id=Mom_ID, year, loc, part)

part <- rbind(part.cap,mm.part.clean)%>%
  mutate(class="Parturition")%>%
  select(id,class,year,value=part, loc)

####Join with preg
part.preg <- rbind(part,preg)


part.preg%>%
  group_by(class,loc)%>%
  summarize(value=mean(value),n=n())


#BOOT to assess differences
part.preg.summary <- part.preg%>%
rep_sample_n(size = nrow(part.preg), replace = TRUE, reps = 1000)%>%
  group_by(replicate,loc,class)%>%
    summarize(value_mean=mean(value))%>%
  group_by(loc,class)%>%
  summarize(  value=median(value_mean),
              lower=quantile(value_mean,0.025),
              upper=quantile(value_mean,0.975))

part.preg%>%
rep_sample_n(size = nrow(part.preg), replace = TRUE, reps = 1000)%>%
  group_by(replicate,loc,class)%>%
    summarize(value_mean=mean(value))%>%
  ggplot(aes(x=loc,y=value_mean, fill=fct_reorder(class, -value_mean)))+
  facet_wrap(vars(fct_reorder(class, -value_mean)))+
  geom_violin()+
  theme_ipsum()+
  labs(y="Rate",x="Penned or Wild", fill="Parameter")
```

![](README_files/figure-gfm/part-1.png)<!-- -->

``` r
##difference in parturition?
part.preg%>%
rep_sample_n(size = nrow(part.preg), replace = TRUE, reps = 10000)%>%
  group_by(replicate,loc,class)%>%
    summarize(value_mean=mean(value))%>%
  filter(class%in%"Parturition")%>%
  select(loc, replicate, value_mean)%>%
  pivot_wider(names_from=loc, values_from=value_mean)%>%
  mutate(dif=P-W)%>%
  ungroup()%>%
  summarise(median=median(dif),
            lower=quantile(dif,0.05),
            upper=quantile(dif,0.95))



null <-glm(value~1,data=part.preg%>%
  filter(class=="Parturition"), family = "binomial")
m1 <- glm(value~loc,data=part.preg%>%
  filter(class=="Parturition"), family = "binomial")

model.sel(null,m1)
```

    ## Model selection table 
    ##      (Intrc) loc          family df  logLik  AICc delta weight
    ## null  0.8391     binomial(logit)  1 -71.029 144.1  0.00  0.534
    ## m1    0.9963   + binomial(logit)  2 -70.130 144.4  0.27  0.466
    ## Models ranked by AICc(x)

## Get data for Avalanches (not run on git, requires some manual moving of tifs, loaded below)

``` r
#https://avalanche.org/wp-content/uploads/2018/08/06_TAR_McCollisterBirkeland.pdf
#https://www.pc.gc.ca/en/pn-np/mtn/securiteenmontagne-mountainsafety/avalanche/echelle-ratings#scale

##get slope
##canopy
#roughness
library(rgee)
library(sf)
library(raster)
ee_Initialize(email="lamb.eco.research@gmail.com")
st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))

##load collars
mcnay <- read_excel(here::here("data", "telem","KlinseZaAllCollars_210418.xlsx"),
                    sheet="KlinseZaAll_210418")%>%
  drop_na(easting)%>%
  st_as_sf(coords=c("easting", "northing"),
        crs=26910)

mcnay <- mcnay%>%
  mutate(Da=str_sub(date_time,0,2),
         Mo=str_sub(date_time,3,5),
         Yr=paste0(20,str_sub(date_time,6,7)),
         Hr=str_sub(date_time,9,10),
         Mn=str_sub(date_time,12,13))%>%
  mutate(datetime=paste0(paste(Yr,Mo,Da,sep="-")," ",paste(Hr,Mn,"00",sep=":"))%>%ymd_hms(),
         Year_=year(datetime))

##towns to clip from telem data (offices collars left on in)
towns <- tribble(~name,~lat,~long,
                 "FSJ", 56.25161403173337, -120.84545438805684,
                 "Mack", 55.3353229828063, -123.0967602935136,
                 "john2", 56.29403813690339, -121.07640843139947)%>%
  st_as_sf(coords=c("long","lat"), crs=4326)%>%
  st_transform(26910)


##remove any locs within the pen and in town
mcnay <- mcnay%>%
  filter(!locn_id%in%c(106995, 32280,1035040))%>%
  st_erase(towns%>%st_as_sf%>%st_buffer(1000))

poly <- sf_as_ee(st_bbox(mcnay)%>%st_as_sfc%>%st_buffer(10000)%>%st_transform(4326))

Map$addLayer(
  eeObject = poly
)

##set pallete
palette = c("#000004FF","#781C6DFF","#ED6925FF","#FCFFA4FF")

##get layerselev <-
elev <-ee$
  Image('USGS/SRTMGL1_003')$
  select('elevation')$
  clip(poly)

slope <- ee$Terrain$slope(elev)$
  clip(poly)

topo.diver <- ee$
  Image('CSP/ERGo/1_0/Global/SRTM_topoDiversity')$
  select('constant')$
  clip(poly)

cc <- ee$
  ImageCollection('NASA/MEASURES/GFCC/TC/v3')$
  select('tree_canopy_cover')$
  filterDate("2010-01-01","2015-12-31")$
  mean()$
  clip(poly)

cc2 <- ee$
  ImageCollection("MODIS/006/MOD44B")$
  select("Percent_Tree_Cover")$
    filterDate("2014-01-01", "2020-12-31")$
    mean()$
    clip(poly)

##PLOT
Map$setCenter(-122.5, 55.5, 7)
Map$addLayer(
  eeObject = slope,
  list(min=0, max=40, palette = palette)
)

Map$addLayer(
  eeObject = topo.diver,
  list(min=0, max=1, palette = palette)
)

Map$addLayer(
  eeObject = cc,
  list(min=0, max=100, palette = palette)
)

Map$addLayer(
  eeObject = cc2,
  list(min=0, max=100, palette = palette)
)

##download rasters
  lyrs <- list(slope, topo.diver, cc)
  lyrs.names <- c("slope", "topo.diver", "cc")
  
  for(j in 1:length(lyrs)){
    ee_image_to_drive(
      image = lyrs[[j]],
      description = paste0(lyrs.names[j]),
      fileNamePrefix=paste0(lyrs.names[j]),
      timePrefix = FALSE,
      folder = "PenningLayers",
      scale = 30,
      maxPixels = 1.0E13,
      fileFormat= 'GeoTIFF',
      crs="EPSG:26910",
      region=poly,
      skipEmptyTiles= TRUE
    )$start()
  }
ee_monitoring(eeTaskList=TRUE)
```

## Prep avalanche layers

``` r
##load rasters
rastlist <- list.files(path =here::here("data","spatial","avalanche"), pattern='.tif$', all.files=TRUE, full.names=TRUE)
allrasters <- stack(rastlist)
plot(allrasters)
```

![](README_files/figure-gfm/Avalanches%20prep-1.png)<!-- -->

``` r
##fix NA's in cc
values(allrasters[["cc"]])[is.na(values(allrasters[["cc"]]))] <- 0

##make avy risk layer
##start by reclassing
allrasters.reclass <- allrasters

allrasters.reclass[["cc"]] <- reclassify(allrasters.reclass[["cc"]],
                     matrix(c(-Inf, 5, 3,
                              5, 30, 2,
                              30, Inf, 1),
                ncol = 3,
                byrow = TRUE))

#plot(allrasters.reclass[["cc"]])


allrasters.reclass[["slope"]] <- reclassify(allrasters.reclass[["slope"]],
                     matrix(c(-Inf, 25, 1,
                              30, 45, 3,
                              25, 30, 2,
                              45,Inf,2),
                ncol = 3,
                byrow = TRUE))
#plot(allrasters.reclass[["slope"]])

allrasters.reclass[["topo.diver"]] <- reclassify(allrasters.reclass[["topo.diver"]],
                     matrix(c(-Inf, 0.33, 1,
                              0.33, 0.66, 2,
                              0.66, Inf, 3),
                ncol = 3,
                byrow = TRUE))

#plot(allrasters.reclass[["topo.diver"]])


avy <- calc(allrasters.reclass, sum)/3

##make anything with a slope of <3 (lakes etc, a 1)
slope0 <- allrasters[["slope"]]
slope0 <- (slope0<3)*100


avy <- slope0+avy
values(avy)[values(avy)>50] <- 1
```

## Plot Avalanche map

``` r
plot(avy)
```

![](README_files/figure-gfm/Avalanches%20map-1.png)<!-- -->

## Add avalanche data to telemetry and then to survival

``` r
mcnay.avy <- mcnay%>%
    filter(animal_id%in%surv.pen$id,
           Year_>=2014)%>%
    dplyr::mutate(avy=raster::extract(avy, .))%>%
  mutate(PenSeason=NA_character_,
         year=NA_real_)


pen.long <- pen%>%dplyr::select(year,penned,free)%>%
  pivot_longer(-year)%>%
  filter(year<2021)

##summarize across penned and unpenned seasons and years

for(i in 1:nrow(pen.long)){

  mcnay.avy<- mcnay.avy%>%
    dplyr::mutate(PenSeason=dplyr::case_when(datetime %within% pen.long[i,]$value~pen.long[i,]$name,
           TRUE~PenSeason),
           year=case_when(datetime %within% pen.long[i,]$value~pen.long[i,]$year,
           TRUE~year))
}



  mcnay.avy.summary <-   mcnay.avy%>%
    tibble%>%
    group_by(animal_id,PenSeason,year)%>%
    summarise(avy=mean(avy,na.rm=TRUE))%>%
    rename(id=animal_id)%>%
    mutate(PenSeason=case_when(PenSeason%in%"penned"~"Pen", 
                          PenSeason%in%"free"~"Free"))


  ##join with survival data
surv.pen <- surv.pen%>%
  left_join(mcnay.avy.summary, by=c("id","PenSeason","year"))

##one animal with missing avy data
##use avergae for animal in other years 1.5 for 314 and 2.1 for 352

surv.pen <- surv.pen%>%
  mutate(avy=case_when(id%in%"CN352S" & year %in%2015~2.1,
                       id%in%"CN314K" & year %in%2015~1.5,
                       TRUE~avy))


ggplot(surv.pen%>%filter(PenSeason%in%"Free"),aes(x=avy, fill=PennedThatYear))+
  geom_density(alpha=0.6)+
      theme_ipsum()+
  labs(y="Count", x="Avalanche rating")+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15),
        axis.text = element_text(size=10),
        legend.text = element_text(size=13),
        legend.title=element_text(size=15))
```

![](README_files/figure-gfm/Avalanches%20to%20telem%20data-1.png)<!-- -->

## Model survival with avalanche risk

``` r
##COXPH survival
m1 <- coxph(Surv(dur, dead)~ 1 +
              frailty(id), data = surv.pen)
m2 <- coxph(Surv(dur, dead)~ PenSeason +
              frailty(id), data = surv.pen)
m3 <- coxph(Surv(dur, dead)~ PenSeason + PennedThatYear +
              frailty(id), data = surv.pen)
m4 <- coxph(Surv(dur, dead)~ PenSeason + PennedThatYear + avy +
              frailty(id), data = surv.pen)

model.sel(m1,m2,m3,m4,  rank="AICc")
```

    ## Model selection table 
    ##    (Int) frl(id) PnS PTY     avy family df   logLik  AICc delta weight
    ## m4     +       +   +   + -0.4557   (NA)  5  -97.269 209.2  0.00      1
    ## m3     +       +   +   +           (NA)  6 -102.117 225.6 16.38      0
    ## m1     +       +                   (NA)  7 -102.791 230.5 21.29      0
    ## m2     +       +   +               (NA)  9 -100.263 237.6 28.45      0
    ## Models ranked by AICc(x)

``` r
summary(m5)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ year + age, data = calf.surv.monthly)
    ## 
    ##   n= 647, number of events= 8 
    ## 
    ##           coef exp(coef)  se(coef)      z Pr(>|z|)  
    ## year -0.253687  0.775934  0.185456 -1.368   0.1713  
    ## age  -0.008282  0.991753  0.004485 -1.847   0.0648 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##      exp(coef) exp(-coef) lower .95 upper .95
    ## year    0.7759      1.289    0.5395     1.116
    ## age     0.9918      1.008    0.9831     1.001
    ## 
    ## Concordance= 0.715  (se = 0.106 )
    ## Likelihood ratio test= 6.5  on 2 df,   p=0.04
    ## Wald test            = 5.8  on 2 df,   p=0.06
    ## Score (logrank) test = 6.47  on 2 df,   p=0.04

``` r
##COXPH survival
m1 <- coxph(Surv(dur, dead)~ 1 , data = surv.pen%>%filter(PenSeason%in%"Free"))
m2 <- coxph(Surv(dur, dead)~ PennedThatYear , data = surv.pen%>%filter(PenSeason%in%"Free"))
m3 <- coxph(Surv(dur, dead)~ PennedThatYear + avy , data = surv.pen%>%filter(PenSeason%in%"Free"))
m4 <- coxph(Surv(dur, dead)~ PennedThatYear*avy , data = surv.pen%>%filter(PenSeason%in%"Free"))


model.sel(m1,m2,m3,m4,  rank="AICc")
```

    ## Model selection table 
    ##    (Int) PTY    avy avy:PTY family      class df  logLik  AICc delta weight
    ## m4     +   + -7.955       +   (NA)      coxph  3 -81.227 170.2  0.00  0.559
    ## m2     +   +                  (NA)      coxph  1 -84.725 171.7  1.53  0.260
    ## m3     +   + -1.038           (NA)      coxph  2 -84.466 173.7  3.56  0.094
    ## m1     +                      (NA) coxph.null  0 -86.951 173.9  3.74  0.086
    ## Models ranked by AICc(x)

``` r
summary(m3)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ PennedThatYear + avy, data = surv.pen %>% 
    ##     filter(PenSeason %in% "Free"))
    ## 
    ##   n= 150, number of events= 18 
    ## 
    ##                    coef exp(coef) se(coef)      z Pr(>|z|)  
    ## PennedThatYearW  1.0963    2.9930   0.5116  2.143   0.0321 *
    ## avy             -1.0385    0.3540   1.4696 -0.707   0.4798  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                 exp(coef) exp(-coef) lower .95 upper .95
    ## PennedThatYearW     2.993     0.3341   1.09797     8.159
    ## avy                 0.354     2.8249   0.01986     6.309
    ## 
    ## Concordance= 0.673  (se = 0.052 )
    ## Likelihood ratio test= 4.97  on 2 df,   p=0.08
    ## Wald test            = 4.67  on 2 df,   p=0.1
    ## Score (logrank) test = 5.05  on 2 df,   p=0.08

``` r
summary(m4)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ PennedThatYear * avy, data = surv.pen %>% 
    ##     filter(PenSeason %in% "Free"))
    ## 
    ##   n= 150, number of events= 18 
    ## 
    ##                               coef      exp(coef)       se(coef)      z Pr(>|z|)  
    ## PennedThatYearW      -13.531063089    0.000001329    5.962422028 -2.269   0.0232 *
    ## avy                   -7.954821128    0.000350966    3.404610944 -2.336   0.0195 *
    ## PennedThatYearW:avy    8.884464354 7218.946839291    3.723799051  2.386   0.0170 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                          exp(coef)     exp(-coef)        lower .95     upper .95
    ## PennedThatYearW        0.000001329 752429.8804457 0.00000000001118        0.1580
    ## avy                    0.000350966   2849.2786269 0.00000044386918        0.2775
    ## PennedThatYearW:avy 7218.946839291      0.0001385 4.88394028400742 10670317.4974
    ## 
    ## Concordance= 0.711  (se = 0.053 )
    ## Likelihood ratio test= 11.45  on 3 df,   p=0.01
    ## Wald test            = 6.76  on 3 df,   p=0.08
    ## Score (logrank) test = 8.63  on 3 df,   p=0.03

``` r
ggforest(m3, data=surv.pen)
```

![](README_files/figure-gfm/Avalanches%20and%20survival-1.png)<!-- -->

``` r
avy.dat <- expand.grid(PennedThatYear=c("W","P"), avy=seq(1.4,2.2, by=0.1), dur=250, dead=0) 
  
avy.dat$pred <- exp(-predict(m4, newdata=avy.dat, type="expected"))

##run by P vs W alone
m2w <- coxph(Surv(dur, dead)~ avy , data = surv.pen%>%filter(PenSeason%in%"Free" & PennedThatYear%in%"W"))
m2p <- coxph(Surv(dur, dead)~ avy , data = surv.pen%>%filter(PenSeason%in%"Free" & PennedThatYear%in%"P"))

summary(m2w)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ avy, data = surv.pen %>% filter(PenSeason %in% 
    ##     "Free" & PennedThatYear %in% "W"))
    ## 
    ##   n= 67, number of events= 12 
    ## 
    ##       coef exp(coef) se(coef)     z Pr(>|z|)
    ## avy 0.8967    2.4516   1.5142 0.592    0.554
    ## 
    ##     exp(coef) exp(-coef) lower .95 upper .95
    ## avy     2.452     0.4079    0.1261     47.68
    ## 
    ## Concordance= 0.498  (se = 0.095 )
    ## Likelihood ratio test= 0.34  on 1 df,   p=0.6
    ## Wald test            = 0.35  on 1 df,   p=0.6
    ## Score (logrank) test = 0.35  on 1 df,   p=0.6

``` r
summary(m2p)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ avy, data = surv.pen %>% filter(PenSeason %in% 
    ##     "Free" & PennedThatYear %in% "P"))
    ## 
    ##   n= 83, number of events= 6 
    ## 
    ##           coef  exp(coef)   se(coef)      z Pr(>|z|)  
    ## avy -7.4450714  0.0005843  3.2901929 -2.263   0.0236 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##     exp(coef) exp(-coef)    lower .95 upper .95
    ## avy 0.0005843       1711 0.0000009248    0.3692
    ## 
    ## Concordance= 0.78  (se = 0.082 )
    ## Likelihood ratio test= 6.06  on 1 df,   p=0.01
    ## Wald test            = 5.12  on 1 df,   p=0.02
    ## Score (logrank) test = 5.66  on 1 df,   p=0.02

``` r
##does avalanche exposure differ between groups?
m1 <- lmer(avy ~  PennedThatYear + (1|id)+ (1|year), data = surv.pen%>%filter(PenSeason%in%"Free"),
                   REML = FALSE)
summary(m1)
```

    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: avy ~ PennedThatYear + (1 | id) + (1 | year)
    ##    Data: surv.pen %>% filter(PenSeason %in% "Free")
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   -143.4   -128.4     76.7   -153.4      145 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1743 -0.5381  0.0548  0.4246  3.6577 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  id       (Intercept) 0.011441 0.10696 
    ##  year     (Intercept) 0.005332 0.07302 
    ##  Residual             0.013233 0.11504 
    ## Number of obs: 150, groups:  id, 44; year, 7
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error t value
    ## (Intercept)      1.70043    0.03511  48.438
    ## PennedThatYearW  0.08135    0.02066   3.938
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## PenndThtYrW -0.264

``` r
lm(avy ~ PennedThatYear, data = surv.pen%>%filter(PenSeason%in%"Free"))%>%summary()
```

    ## 
    ## Call:
    ## lm(formula = avy ~ PennedThatYear, data = surv.pen %>% filter(PenSeason %in% 
    ##     "Free"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.36129 -0.13090 -0.00987  0.09986  0.48061 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value             Pr(>|t|)    
    ## (Intercept)      1.69623    0.01854  91.474 < 0.0000000000000002 ***
    ## PennedThatYearW  0.09245    0.02775   3.332              0.00109 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1689 on 148 degrees of freedom
    ## Multiple R-squared:  0.06979,    Adjusted R-squared:  0.0635 
    ## F-statistic:  11.1 on 1 and 148 DF,  p-value: 0.001089

``` r
surv.pen%>%filter(PenSeason%in%"Free")%>%
  group_by(PennedThatYear)%>%
  summarise(mean=mean(avy))
```
