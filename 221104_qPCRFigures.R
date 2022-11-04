######################################
##### Date Created  ## 17 Oct 22 #####
##### Date Modified ## 25 Oct 22 #####
######################################

library(readxl)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate) # to fix inconsistent date/time formats imported from excel
library(patchwork)
library(stringr) # working with strings
library(zoo) # for timeseries data maybe?
#library(tseries) # might work better?
library(worms)
library(ggh4x)
library(kableExtra) # makes raw data tables look 
library(ggExtra)
library(viridis)

options(scipen = 999)

qPCRfiles=list.files(path = "/Users/morrisonme/Documents/Projects", 
                     pattern = "opt-eDNA-[0-9]{2}_qPCR_data", 
                     full.names = TRUE)

qPCRfiles

metadata = lapply(qPCRfiles, function(x) read_excel(x, sheet = 2))
samples = lapply(qPCRfiles, function(x) read_excel(x, sheet = 3))


names(metadata) = qPCRfiles
names(samples) = qPCRfiles


# Ensure dates are in unified format across datasets.
# In character format to make graphing possible
for (i in 1:length(metadata)){
  if ((typeof(metadata[[i]]$eventDate) == "double") == TRUE ) {
    metadata[[i]]$newDate = as.character(metadata[[i]]$eventDate)
  } else {
    ifelse( 
      ((nchar(metadata[[i]]$eventDate) == 5) == TRUE),
      (metadata[[i]]$newDate = as.character(as.Date(as.numeric(metadata[[i]]$eventDate), origin = "1899-12-30"))),
      (metadata[[i]]$newDate2 = 
         as.character(as.Date(parse_date_time(metadata[[i]]$eventDate, orders = "d m y")))))
  }
  metadata[[i]]$eventDate = metadata[[i]]$newDate  # merge Date columns
  metadata[[i]]$eventDate[!is.na(metadata[[i]]$newDate2)] = metadata[[i]]$newDate2[!is.na(metadata[[i]]$newDate2)]
}

# match event date to samples. 
# Always have the vector of smaller size as second arg in match()
for (j in 1:length(samples)) {
  metadata[[j]] = metadata[[j]] %>% 
    drop_na(eventDate)
  samples[[j]]$date = metadata[[j]]$eventDate[match(samples[[j]]$eventID, metadata[[j]]$eventID)]
  samples[[j]]$ecoregion = metadata[[j]]$ecoregion[match(samples[[j]]$eventID, metadata[[j]]$eventID)]
  samples[[j]]$site = metadata[[j]]$samplingSite[match(samples[[j]]$eventID, metadata[[j]]$eventID)]
  samples[[j]]$year = year(samples[[j]]$date)
  samples[[j]]$month = month(samples[[j]]$date)
  samples[[j]]$projectID = substr(samples[[j]]$eventID, 1,2) %>%
    str_remove("-")
  samples[[j]] = mutate(samples[[j]], det_prob = case_when(
      occurrenceStatus == "Detected" ~ 1.00,
      occurrenceStatus == "Suspected" ~ 0.66,
      occurrenceStatus == "Inconclusive" ~ 0.33,
      occurrenceStatus == "Not detected" ~ 0.00)) %>%
    mutate(samples[[j]], detected = case_when(
      occurrenceStatus != "Not detected" ~ 1,
      occurrenceStatus == "Not detected" ~ 0)) %>%
    mutate(spp.labs = paste0(substr(samples[[j]]$scientificName, 1,1), '.'," ", 
                             word(samples[[j]]$scientificName,-1)))
  } 

samples = lapply(samples, function(x) x %>% drop_na(date))


## RELATIVE SUCCESS BY SITE - PERCENT POSITIVE DETECTIONS 
# calculates the percentage of positive samples
rel_success.site = lapply(samples, function(x) {
  x %>% 
    group_by(projectID,year,month,spp.labs,site, detected) %>%
    dplyr::summarise(cnt=n()) %>%
    transmute(detected,spp.labs, freq = cnt/sum(cnt)) %>%
    mutate(percent_pos = case_when(
      detected == 1 ~ freq,
      detected == 0 ~ 0)) }
  )

# combine into one df for plotting
succ.site.df = rel_success.site %>%
  map_df(~bind_rows(.x)) %>%
  subset(select=-c(freq))


# combine all sample df into a single df for plotting
#samples.df = samples %>%
#  map_df(~bind_rows(.x)) 
  
#test = na.omit(samples.df)



#na.months = lapply(rel_success, function(x) {
#  x %>%
 #   right_join(spp, tibble("month" = 1:12), by = "month")%>%
#as_tibble() %>%
#    fill(spp, c("projectID", "year",'spp.labs')) #%>%
  #  mutate_at(c('detected','freq','percent_pos'), ~replace_na(.,0))
#}}
#)

#na.months.df = na.months %>%
#  map_df(~bind_rows(.x)) %>%
#  subset(select=-c(freq))

# test.ciona= subset(rel_success_df, spp.labs %in% c("C. intestinalis")) 

## RELATIVE SUCCESS PER SITE ##
# Need to figure out how to not link the points between samples if non-consecutive
ggplot(subset(rel_success_df, spp.labs %in% c("B. violaceus")))+
  geom_point(mapping=aes(as.factor(month),percent_pos, group=site, col=site), 
             size=1.7, position=position_jitter(width=0.1,height=0.1)) +
  facet_grid(year~spp.labs)+
  xlab('Sampling Date')+
  ylab('Relative success \n(Percent positive samples)')+
  labs(col = "Sample Site")+
  #ggtitle(label='A')+
  scale_x_discrete(limits=as.factor(c(1:12)),
                   breaks = c(1,2,3,4,5,6,
                              7,8,9,10,11,12),
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1),
                     labels = c("0%","25%","50%","75%","100%"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),          
        panel.background = element_blank(), panel.border = element_rect(colour="black", fill=NA),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r=10)),
        plot.title = element_text(),
        axis.text.x=element_text(angle=60, hjust=2, vjust = 1.5))#,
        #legend.position="none")


## RELATIVE SUCCESS - PERCENT POSITIVE DETECTIONS 
# calculates the percentage of positive samples
rel_success = lapply(samples, function(x) {
  x %>% 
    group_by(projectID,year,month,spp.labs,detected) %>%
    dplyr::summarise(cnt=n()) %>%
    transmute(detected,spp.labs, freq = cnt/sum(cnt)) %>%
    mutate(percent_pos = case_when(
      detected == 1 ~ freq)) %>%
    drop_na()}
)

# combine into one df for plotting
success.df = rel_success %>%
  map_df(~bind_rows(.x)) %>%
  subset(select=-c(freq))


# RELATIVE SUCCESS BY PROJECT
ggplot(subset(success.df, spp.labs %in% c("B. violaceus")))+
  geom_line(mapping=aes(as.factor(month), percent_pos, group=spp.labs, col=projectID),
            size=1) +
  facet_grid(year~spp.labs, drop=TRUE)+
  xlab('Sampling Date')+
  ylab('Relative success \n(Percent positive samples)')+
  labs(col = "Project")+
  #ggtitle(label='A')+
  scale_x_discrete(limits=as.factor(c(1:12)),
                   breaks = c(1,2,3,4,5,6,
                              7,8,9,10,11,12),
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  scale_y_continuous(limits=c(0,1),
                     breaks = c(0.00,0.25,0.50,0.75,1),
                     labels = c("0%","25%","50%","75%","100%"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),          
        panel.background = element_blank(), panel.border = element_rect(colour="black", fill=NA),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r=10)),
        plot.title = element_text(),
        axis.text.x=element_text(angle=60, hjust=1))



# ABSOLUTE NUMBER OF POSITIVE SAMPLES
success = lapply(samples, function(x) {
  x %>% 
    group_by(projectID,year,month,spp.labs,detected) %>%
  dplyr::summarise(cnt=n()) %>%
  mutate(num_pos = case_when(
    detected == 1 ~ cnt
   )) }
  )



success.df = success %>%
  map_df(~bind_rows(.x)) %>%
  subset(select=-c(cnt)) 

#success.df[is.na(success.df)] <- 0


ggplot(subset(success.df, spp.labs %in% c("C. intestinalis", "P. vitulina")))+
  geom_point(mapping=aes(as.factor(month), num_pos, group=projectID, col=projectID),
            size=1.5) +
  geom_line(mapping=aes(as.factor(month), num_pos, group=spp.labs, col=projectID),
            alpha=0.5)+
  facet_grid(year~spp.labs, drop=TRUE)+
  xlab('Sampling Date')+
  ylab('Absolute number of positive samples')+
  labs(col = "Project")+
  #ggtitle(label='A')+
  scale_x_discrete(limits=as.factor(c(1:12)),
                   breaks = c(1,2,3,4,5,6,
                              7,8,9,10,11,12),
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
 # scale_y_continuous(limits=c(0,1),
 #                    breaks = c(0.00,0.25,0.50,0.75,1),
 #                    labels = c("0%","25%","50%","75%","100%"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),          
        panel.background = element_blank(), panel.border = element_rect(colour="black", fill=NA),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r=10)),
        plot.title = element_text(),
        axis.text.x=element_text(angle=60, hjust=1))


# AVERAGE qPCR DETECTION RATE PER SAMPLE
det.rate = lapply(samples, function(x) 
    subset(x, select=c(projectID,year,month,spp.labs,det_prob))
  )

det.rate.df = det.rate %>%
  map_df(~bind_rows(.x))

ggplot(subset(subset(det.rate.df, year %in% c("2018","2020","2022")),
              spp.labs %in% c("C. intestinalis","S. clava")))+
  geom_point(mapping=aes(as.factor(month),det_prob, group=projectID, col=projectID),
             position=position_jitter(height=0.1, width=0.1)) +
  facet_grid(year~spp.labs, drop=TRUE)+
  #facet_wrap(~scientificName, scales='free_x')+
  xlab('Sampling Date')+
  ylab('Detection Rate (qPCR replicates)')+
  labs(col = "Project")+
  #ggtitle(label='A')+
  scale_x_discrete(limits=as.factor(c(1:12)),
                   breaks = c(1,2,3,4,5,6,
                              7,8,9,10,11,12),
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  scale_y_continuous(breaks=c(0.00,0.33,0.66,1.00), labels=c("0/3",'1/3','2/3','3/3'),
                     limits=c(0,1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),          
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r=10)),
        plot.title = element_text(),
        axis.text.x=element_text(angle=60, hjust=1))


ggplot(subset(subset(det.rate.df, year %in% c("2018","2020","2022")),
              spp.labs %in% c("C. intestinalis","S. clava")))+
  geom_jitter(mapping=aes(as.factor(month),det_prob, group=projectID, col=projectID)) +
  geom_smooth(mapping=aes(as.factor(month),det_prob, group=projectID, col=projectID))+
  facet_grid(year~spp.labs, drop=TRUE)+
  #facet_wrap(~scientificName, scales='free_x')+
  xlab('Sampling Date')+
  ylab('Detection Rate (qPCR replicates)')+
  labs(col = "Project")+
  #ggtitle(label='A')+
  scale_x_discrete(limits=as.factor(c(1:12)),
                   breaks = c(1,2,3,4,5,6,
                              7,8,9,10,11,12),
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  scale_y_continuous(breaks=c(0.00,0.33,0.66,1.00), labels=c("0/3",'1/3','2/3','3/3'),
                     limits=c(0,1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),          
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour="black", fill=NA),
        axis.title.y = element_text(margin = margin(r=10)),
        plot.title = element_text(),
        axis.text.x=element_text(angle=60, hjust=1))


# AVERAGE CONCENTRATION - DUAL Y-AXES FOR DIFFERENT UNITS

avg.conc = lapply(samples, function(x) {
  x %>%
    subset(select=c(projectID,year,month,site, spp.labs,quantityMean,quantityUnits))%>%
  mutate_at(c('quantityMean'), as.numeric)
})


avg.conc.df = avg.conc %>%
  map_df(~bind_rows(.x))

avg.conc.df[is.na(avg.conc.df)] = 0
avg.conc.df$logMean = log(1+avg.conc.df$quantityMean)

# split dataframe by unique quantityUnits
conc.units.df = split(avg.conc.df, f=avg.conc.df$quantityUnits)
list2env(conc.units.df,envir=.GlobalEnv)




ggplot(subset(subset(`copies/rxn`, year %in% c("2018","2020","2022")),
                         spp.labs %in% c("C. intestinalis","S. clava")))+
  geom_boxplot(mapping=aes(as.factor(month), exp(logMean), group=month, col=projectID),
             size=1) +
  geom_boxplot(data= subset(subset(`pg/rxn`, year %in% c("2018","2020","2022")),
              spp.labs %in% c("C. intestinalis","S. clava")), 
              mapping=aes(as.factor(month), exp(logMean), group=month, col=projectID),
              size=1) +
  #geom_line(mapping=aes(as.factor(month), num_pos, group=spp.labs, col=projectID),
 #           alpha=0.5)+
  facet_grid(year~spp.labs, drop=TRUE)+
  xlab('Sampling Date')+
 # ylab('Absolute number of positive samples')+
  labs(col = "Project")+
  #ggtitle(label='A')+
  scale_x_discrete(limits=as.factor(c(1:12)),
                   breaks = c(1,2,3,4,5,6,
                              7,8,9,10,11,12),
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
   scale_y_continuous(trans="log", name="Copies/reaction", 
                      breaks = c(1,10,100, 1000,3000,6000 ),
                      sec.axis = sec_axis(trans=~.,name="pg/reaction",
                                          breaks = c(1,10,100, 1000,3000,6000 )))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),          
        panel.background = element_blank(), panel.border = element_rect(colour="black", fill=NA),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r=10)),
        plot.title = element_text(),
        axis.text.x=element_text(angle=60, hjust=1))

## HEATMAP OF PRESENCE/ABSENCE FOR EACH MONTH/REGION

hm.lst = lapply(samples, function(x) {
  x %>% 
    group_by(projectID,ecoregion,year,month,spp.labs,detected) %>%
    dplyr::summarise(cnt=n()) %>%
    transmute(detected,spp.labs, freq = cnt/sum(cnt)) %>%
    mutate(percent_pos = case_when(
      detected == 1 ~ freq)) }# %>%
    #drop_na() }
)



# combine into one df for plotting
hm.df = hm.lst %>%
  map_df(~bind_rows(.x))

ggplot()+
  geom_tile(subset(hm.df, projectID %in% 6),mapping=aes(as.factor(month),spp.labs,fill=percent_pos, color= "Not detected"),size=0.1) + 
  scale_fill_viridis(name="       Relative success \n(percent positive samples)",
                      guide="colourbar", trans = "reverse", na.value = "gray94",
                      labels=c("","10%","","30%","","50%")) + 
  scale_colour_manual(values="white") +              
  guides(fill = guide_colourbar(order = 1, ticks = FALSE, label.hjust = -0.0001,
                                label.position = "bottom",
                                title.position="left"), 
         colour=guide_legend(NULL,override.aes=list(fill="grey94"), 
                             order=2, label.position="bottom"))+
  facet_grid(year~ecoregion) +
 # scale_y_continuous(trans = "reverse", breaks = as.factor(1:12))) +
  scale_x_discrete(limits=as.factor(c(1:12)),
                   breaks = c(1,2,3,4,5,6,
                              7,8,9,10,11,12),
                   labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  theme_minimal(base_size = 10) +
  labs(title= paste("GRDI Predator eDNA Detection"), x="Month", y="Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),          
        panel.background = element_blank(), panel.border = element_rect(colour="black", fill=NA),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(margin = margin(r=10)),
        axis.title.x = element_text(margin=margin(t=10)),
        plot.title = element_text(hjust=-2),
        axis.text.x=element_text(angle=60, hjust=1),
        axis.text.y=element_text(face = "italic"), legend.position = "bottom",
        strip.text = element_text (margin = margin (10,30,5,5)))+
  removeGrid()#ggExtra

 ggsave("pred_heatmap.jpg", device="jpg", width = 12, height =10, unit = "cm", dpi = 300)
 