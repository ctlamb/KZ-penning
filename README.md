Klinse-Za Penning Demography + Range Recolonization
================
Clayton Lamb
10 May, 2021

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
library(chisq.posthoc.test)
source(here::here("functions","seasonalsurvival_fn.r"))

#ADULT FEMALE SURVIVAL
df <- read_excel(here::here("data","AllHerdsStartEndDates_210421.xlsx"))%>%
  filter(herd%in%c("K","S"),
         Sex%in%"F")
  
#ADULT FEMALE LOCATIOn (IN/OUT OF PEN)
ad.cap <- read_excel(here::here("data","Klinse-za capture and demographics_210426.xlsx"),
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
calf <- read_excel(here::here("data","Klinse-za capture and demographics_210506.xlsx"),
                     sheet="Calf-Capt")%>%
    drop_na("Mom_ID")

##MORTS
mort <- read_excel(here::here("data","Klinse-za capture and demographics_210426.xlsx"),
                     sheet="Mortalities")

##PARTURITION
mm.part <- read_csv(here::here("data","MovementModelResults.csv"))%>%
  drop_na(Caribou_ID)

##CONDITION
cond <- read_csv(here::here("data","KZONCP_Health_Year2_results_210422.csv"))
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
  left_join(read_excel(here::here("data","Klinse-za capture and demographics_210426.xlsx"),
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
    ##      Variance of random effect= 5e-07   I-likelihood = -72.1 
    ## Degrees of freedom for terms= 1 1 1 0 
    ## Concordance= 0.862  (se = 0.037 )
    ## Likelihood ratio test= 26.68  on 3 df,   p=7e-06

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

    ## [1] 0.7337618

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
    ##   365     36      11    0.824  0.0483        0.735        0.925

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
    ## m3       +     -0.4621           (NA)      coxph  1 -29.993 62.7  0.00  0.411
    ## m6       +             -0.3125   (NA)      coxph  1 -30.711 64.1  1.44  0.201
    ## m1       +                       (NA) coxph.null  0 -32.131 64.3  1.61  0.184
    ## m4       +   + -0.4337           (NA)      coxph  2 -29.703 65.8  3.15  0.085
    ## m2       +   +                   (NA)      coxph  1 -31.704 66.1  3.42  0.074
    ## m5       +   +         -0.3087   (NA)      coxph  2 -30.340 67.1  4.43  0.045
    ## Models ranked by AICc(x)

``` r
cox.zph(m3)
```

    ##            chisq df       p
    ## ReleaseAge  13.9  1 0.00019
    ## GLOBAL      13.9  1 0.00019

``` r
cox.zph(m6)
```

    ##        chisq df      p
    ## year    7.86  1 0.0051
    ## GLOBAL  7.86  1 0.0051

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
    ##   n= 62, number of events= 8 
    ## 
    ##               coef exp(coef) se(coef)      z Pr(>|z|)  
    ## ReleaseAge -0.4621    0.6300   0.2337 -1.977   0.0481 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##            exp(coef) exp(-coef) lower .95 upper .95
    ## ReleaseAge      0.63      1.587    0.3984    0.9961
    ## 
    ## Concordance= 0.605  (se = 0.143 )
    ## Likelihood ratio test= 4.27  on 1 df,   p=0.04
    ## Wald test            = 3.91  on 1 df,   p=0.05
    ## Score (logrank) test = 3.78  on 1 df,   p=0.05

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
              age=mean(age))


###try again
##COXPH
m1 <- coxph(Surv(dur, dead)~ 1, data =calf.surv.monthly)
m2 <- coxph(Surv(dur, dead)~ ReleaseAge, data =calf.surv.monthly)
m3 <- coxph(Surv(dur, dead)~ ReleaseAge + age, data =calf.surv.monthly)
m4 <- coxph(Surv(dur, dead)~ ReleaseAge*age, data =calf.surv.monthly)

model.sel(m1,m2,m3,m4,  rank="AICc")
```

    ## Model selection table 
    ##    (Int)     RlA       age  age:RlA family      class df  logLik  AICc delta weight
    ## m4     + -1.0790 -0.051420 0.005464   (NA)      coxph  3 -44.323 100.6  0.00  0.376
    ## m2     + -0.4914                      (NA)      coxph  1 -49.021 100.7  0.06  0.365
    ## m1     +                              (NA) coxph.null  0 -51.313 102.6  1.98  0.140
    ## m3     + -0.4505 -0.004699            (NA)      coxph  2 -48.267 102.9  2.29  0.120
    ## Models ranked by AICc(x)

``` r
cox.zph(m2)
```

    ##            chisq df    p
    ## ReleaseAge 0.473  1 0.49
    ## GLOBAL     0.473  1 0.49

``` r
cox.zph(m3)
```

    ##            chisq df    p
    ## ReleaseAge 0.393  1 0.53
    ## age        0.833  1 0.36
    ## GLOBAL     1.075  2 0.58

``` r
cox.zph(m4)
```

    ##                chisq df    p
    ## ReleaseAge     0.138  1 0.71
    ## age            0.526  1 0.47
    ## ReleaseAge:age 0.446  1 0.50
    ## GLOBAL         0.574  3 0.90

``` r
summary(m4)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ ReleaseAge * age, data = calf.surv.monthly)
    ## 
    ##   n= 641, number of events= 8 
    ## 
    ##                     coef exp(coef)  se(coef)      z Pr(>|z|)    
    ## ReleaseAge     -1.079292  0.339836  0.318915 -3.384 0.000714 ***
    ## age            -0.051422  0.949877  0.015120 -3.401 0.000672 ***
    ## ReleaseAge:age  0.005464  1.005479  0.001613  3.387 0.000706 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                exp(coef) exp(-coef) lower .95 upper .95
    ## ReleaseAge        0.3398     2.9426    0.1819    0.6349
    ## age               0.9499     1.0528    0.9221    0.9784
    ## ReleaseAge:age    1.0055     0.9946    1.0023    1.0087
    ## 
    ## Concordance= 0.746  (se = 0.102 )
    ## Likelihood ratio test= 13.98  on 3 df,   p=0.003
    ## Wald test            = 19.11  on 3 df,   p=3e-04
    ## Score (logrank) test = 19.1  on 3 df,   p=3e-04

``` r
summary(m3)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ ReleaseAge + age, data = calf.surv.monthly)
    ## 
    ##   n= 641, number of events= 8 
    ## 
    ##                 coef exp(coef)  se(coef)      z Pr(>|z|)  
    ## ReleaseAge -0.450454  0.637338  0.229431 -1.963   0.0496 *
    ## age        -0.004699  0.995312  0.003987 -1.179   0.2385  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##            exp(coef) exp(-coef) lower .95 upper .95
    ## ReleaseAge    0.6373      1.569    0.4065    0.9992
    ## age           0.9953      1.005    0.9876    1.0031
    ## 
    ## Concordance= 0.612  (se = 0.136 )
    ## Likelihood ratio test= 6.09  on 2 df,   p=0.05
    ## Wald test            = 6.07  on 2 df,   p=0.05
    ## Score (logrank) test = 5.98  on 2 df,   p=0.05

``` r
summary(m2)
```

    ## Call:
    ## coxph(formula = Surv(dur, dead) ~ ReleaseAge, data = calf.surv.monthly)
    ## 
    ##   n= 641, number of events= 8 
    ## 
    ##               coef exp(coef) se(coef)      z Pr(>|z|)  
    ## ReleaseAge -0.4914    0.6118   0.2370 -2.074   0.0381 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##            exp(coef) exp(-coef) lower .95 upper .95
    ## ReleaseAge    0.6118      1.635    0.3845    0.9734
    ## 
    ## Concordance= 0.574  (se = 0.148 )
    ## Likelihood ratio test= 4.59  on 1 df,   p=0.03
    ## Wald test            = 4.3  on 1 df,   p=0.04
    ## Score (logrank) test = 4.08  on 1 df,   p=0.04

``` r
calf.pred.boot <- data.frame()
##Boot for plot to show uncertainty
for(i in 1:1000){
  m3.i <- coxph(Surv(dur, dead)~ ReleaseAge*age, data = calf.surv.monthly%>%sample_frac(1,replace=TRUE))
 
  calf.pred.i <- expand.grid(ReleaseAge=5:12, age=seq(0,360,by=20), dur=365, dead=0, iter=i) 
  
  calf.pred.i$pred <- exp(-predict(m3.i, newdata=calf.pred.i, type="expected"))
  calf.pred.i$pred
  calf.pred.boot <- rbind(calf.pred.boot,calf.pred.i)
}


calf.pred.boot%>%
         drop_na()%>%
         group_by(ReleaseAge,age)%>%
         summarize(median=median(pred))%>%
  ggplot()+
  geom_tile(aes(y=ReleaseAge, x=age,fill=median))
```

![](README_files/figure-gfm/model%20calf-8.png)<!-- -->

``` r
calf.pred.boot <- data.frame()
##Boot for plot to show uncertainty
for(i in 1:1000){
  m3.i <- coxph(Surv(dur, dead)~ ReleaseAge, data = calf.surv.monthly%>%sample_frac(1,replace=TRUE))
 
  calf.pred.i <- expand.grid(ReleaseAge=5:14,  dur=365, dead=0, iter=i) 
  
  calf.pred.i$pred <- exp(-predict(m3.i, newdata=calf.pred.i, type="expected"))
  calf.pred.i$pred
  calf.pred.boot <- rbind(calf.pred.boot,calf.pred.i)
}


ggplot(data=calf.pred.boot%>%
         drop_na()%>%
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
    ## m1    1.1790   + binomial(logit)  2 -64.625 133.4  0.00    0.6
    ## null  0.9605     binomial(logit)  1 -66.068 134.2  0.81    0.4
    ## Models ranked by AICc(x)
