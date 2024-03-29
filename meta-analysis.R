#Title:Zoos as sentinels – A meta-analysis of seroprevalence of terrestrial mammalian viruses in zoos
#Authors: Van Leeuwen, Pauline, Sarah Falconer, Jasmine Veitch, Breanna Pyott, Bryan Hughes, Isabelle Zimmerman, Albrecht Schulte-Hostedde
#import data.txt & virus.txt
library(dplyr)
library(reshape2)
library(ggplot2)
library(forcats)
library(meta)
library(metafor)
library(tidyr)
#data manipulation
datafull = data %>% inner_join(virus, by="v_sp") #incorporate virus data and full dataset
#keeping only variables of interest
datasimple=datafull %>%
  select(host_order, host_fam,host_status_2021, host_continent, v_class.x,v_fam.x, v_sp, host_pop, zoo_loc2,same_continent,zoo_type, number_faciliti, collection_firstyear, collection_lastyear, vector, risk_to_human, paper_ID,pub_year, sample_size_class, sample_size, positive_sample, issue)
head(datasimple)
datawildcaptive = subset(datasimple,datasimple$host_pop!="domestic") #subset for wild/captive animals comparison
datacaptive = subset(datasimple, datasimple$host_pop=="captive") #subset for only captive animals
datacaptive = subset(datacaptive, datacaptive$issue=="ok") #subset for articles satisfying criteria
datacaptive = subset(datacaptive, datacaptive$zoo_type!="reserve") #subset for animals not in reserves
summary(datacaptive)
sum(datacaptive$sample_size)

   #1- Overall test for heterogeneity between studies with the double-arcsine transformation: CAPTIVE ONLY dataset
      #FIRST: Overall effect size, no moderators
      ies.da=escalc(xi=positive_sample, ni=sample_size, data=datacaptive, measure="PFT", add=0) #calculation effect sizes
      pes.da=rma(yi, vi, data=ies.da) #random effects model
      print(pes.da) ##I^2 73.94%, high heterogeneity between studies, estimate is 0.4085 ci.lb=0.3874 ci.up= 0.4297 se=0.0108    p<.0001
      pes=predict(pes.da, transf=transf.ipft.hm, targ=list(ni=datacaptive$sample_size))
      print(pes) #transforming back to real values
      funnel(pes.da, atransf=transf.ipft.hm, targ=list(ni=datacaptive$sample_size), yaxis="sei")#assessing distribution (few outliers, looks good)
      #detecting publication bias with Egger's regression test
      regtest(pes.da, model="lm", predictor="sei") 
      #Egger test shows that the funnel plot is not significantly asymmetrical p=0.3977
     
      
      
  #2- CAPTIVE DATASET ONLY: Subgroup analysis based on PUBLICATION VARIABLES
    #2-A 1st subgroup analysis: sample size class
      host_paper_or=datacaptive %>%
        dplyr::count(sample_size_class, paper_ID, sort=TRUE)
      host_paper_or %>%
        dplyr::count(sample_size_class, sort=TRUE)   #nb paper for each category
      datacaptive %>%
        group_by(sample_size_class) %>% 
        summarise(nb = sum(sample_size)) #sample size by each category   
      #subgroup analysis
      model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive,studlab=paper_ID,byvar=sample_size_class,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                               sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
      summary(model_subgroup) #significant p< 0.0001 higher infection in studies with great sample sizes (more than 10)
      forest(model_subgroup,comb.random = T,  overall=T,subgroup=T,study.results=F, label.test.overall.random=T, test.overall	=T,pooled.events=T, pooled.totals=T, prediction=F )
    
    #2-B 2nd subgroup analysis: year of publication
      host_paper_or=datacaptive %>%
        dplyr::count(pub_year, paper_ID, sort=TRUE)
      host_paper_or %>%
        dplyr::count(pub_year, sort=TRUE)   #nb paper for each category
      datacaptive %>%
        group_by( pub_year) %>% 
        summarise(nb = sum(sample_size)) #sample size by each category   
      #subgroup analysis
      model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive,studlab=paper_ID,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                               sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
      model_metareg= metareg(~ pub_year, x = model_subgroup)
      summary(model_metareg)#Test moderator not significant
      #Visual representation
      datapub= datacaptive %>%
        group_by(paper_ID,pub_year) %>%
        summarise(positive_sample = sum(positive_sample),
                  sample_size = sum(sample_size)) %>%
        ungroup()
      datapub$prop=datapub$positive_sample/datapub$sample_size
      head(datapub)
      pes.subgroup=predict(model_metareg,newmods=cbind(seq(from=2002, to=2021, by=1)), transf=transf.ipft.hm,targ=list(ni=datacaptive$sample_size),addx=TRUE)
      preds.df <- data.frame(pes.subgroup$pred, pes.subgroup$ci.lb, pes.subgroup$ci.ub, pes.subgroup$cr.ub)
      preds.df
      sequence=seq(from=2002, to=2021, by=1)
      ggplot() +
        geom_point(data=datapub,  aes(x=pub_year, y=prop, size=sample_size)) +
        geom_line(data=preds.df, aes(x=sequence, y=pes.subgroup$pred))+
        geom_line(data=preds.df, aes(x=sequence, y=pes.subgroup$ci.lb), linetype=2)+
        geom_line(data=preds.df, aes(x=sequence, y=pes.subgroup$ci.ub), linetype=2)+
        labs(x="Year of publication", y="Seropositivity prevalence")+
        labs(size = "Sample Size")
      
    #2-C 3rd subgroup analysis: length of data collection
      datacaptive$length_study = datacaptive$collection_lastyear-datacaptive$collection_firstyear
      host_paper_or=datacaptive %>%
        dplyr::count(length_study, paper_ID, sort=TRUE)
      host_paper_or %>%
        dplyr::count(length_study, sort=TRUE)   #nb paper for each category
      datacaptive %>%
        group_by(length_study) %>% 
        summarise(nb = sum(sample_size)) #sample size by each category   
      datacaptive1=subset(datacaptive, datacaptive$length_study!="NA") #not taking NA values (one study)
      #subgroup analysis
      model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive1,studlab=paper_ID,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                               sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
      model_metareg= metareg(~ length_study, x = model_subgroup)
      summary(model_metareg)#Test moderator p<0.0001. large heterogeneity due to length oh study, weak negative relationship.
      #Visual representation
      datapub= datacaptive1 %>%
        group_by(paper_ID,length_study) %>%
        summarise(positive_sample = sum(positive_sample),
                  sample_size = sum(sample_size)) %>%
        ungroup()
      datapub$prop=datapub$positive_sample/datapub$sample_size
      head(datapub)
      pes.subgroup=predict(model_metareg,newmods=cbind(seq(from=0, to=29, by=1)), transf=transf.ipft.hm,targ=list(ni=datacaptive1$sample_size),addx=TRUE)
      preds.df <- data.frame(pes.subgroup$pred, pes.subgroup$ci.lb, pes.subgroup$ci.ub, pes.subgroup$cr.ub)
      preds.df
      sequence=seq(from=0, to=29, by=1)
      ggplot() +
        geom_point(data=datapub,  aes(x=length_study, y=prop, size=sample_size)) +
        geom_line(data=preds.df, aes(x=sequence, y=pes.subgroup$pred))+
        geom_line(data=preds.df, aes(x=sequence, y=pes.subgroup$ci.lb), linetype=2)+
        geom_line(data=preds.df, aes(x=sequence, y=pes.subgroup$ci.ub), linetype=2)+
        labs(x="Duration of sampling collection (years)", y="Seropositivity prevalence")+
        labs(size = "Sample Size")
      
      #2-D 4th subgroup analysis: number of facility by paper
      #subgroup analysis
      model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive,studlab=paper_ID,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                               sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
      model_metareg= metareg(~ number_faciliti, x = model_subgroup)
      summary(model_metareg)#Test moderator p<0.0001. large heterogeneity due to length oh study, weak negative relationship.
      
      
  #3- CAPTIVE ONLY: Subgroup analysis according to HOST VARIABLES
      #3-A 1st subgroup analysis: zoo type
      host_paper_or=datacaptive %>%
        dplyr::count(zoo_type, paper_ID, sort=TRUE)
      host_paper_or %>%
        dplyr::count(zoo_type, sort=TRUE)   #nb paper for each category
      datacaptive %>%
        group_by( zoo_type) %>% 
        summarise(nb = sum(sample_size)/11300) #sample size by each category   
      #subgroup analysis
      model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive,studlab=paper_ID,byvar=zoo_type,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                               sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
      summary(model_subgroup) #significantly p= 0.0231, high seroprevalence in reserves
      forest(model_subgroup,comb.random = T,  overall=T,subgroup=T,study.results=F, label.test.overall.random=T, test.overall	=T,pooled.events=T, pooled.totals=T, prediction=F )
      
    #3-B 2nd subgroup analysis: zoo continent=host continent (or partially)
      host_paper_or=datacaptive %>%
        dplyr::count(same_continent, paper_ID, sort=TRUE)
      host_paper_or %>%
        dplyr::count(same_continent, sort=TRUE)   #nb paper for each category
      datacaptive %>%
        group_by( same_continent) %>% 
        summarise(nb = sum(sample_size)/11300) #sample size by each category   
        #subgroup analysis
          model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive,studlab=paper_ID,byvar=same_continent,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                               sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
          summary(model_subgroup) #significantly p= 0.0049, high seroprevalence when in same continent
          forest(model_subgroup,comb.random = T,  overall=T,subgroup=T,study.results=F, label.test.overall.random=T, test.overall	=T,pooled.events=T, pooled.totals=T, prediction=F )
          
      #3-C 3nd subgroup analysis: zoo_continent
          datacaptive %>%
            dplyr::count(zoo_loc2, sort=TRUE)   #nb records for each category
          datacaptive %>%
            group_by( zoo_loc2) %>% 
            summarise(nb = sum(sample_size)/11300) #sample size by each category   
        #subgroup analysis
          model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive,studlab=paper_ID,byvar=zoo_loc2,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                                   sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
          summary(model_subgroup) #significantly p= 0.0001 greatest infection in american zoos and lower in Africa and europe (see forest plot).
          forest(model_subgroup,comb.random = T,  overall=T,subgroup=T,study.results=F, label.test.overall.random=T, test.overall	=T,pooled.events=T, pooled.totals=T, prediction=F )
          
       #3-D 4th subgroup analysis: host_status_2021 (red list UICN from june 2021)
          datacaptive %>%
            dplyr::count(host_status_2021, sort=TRUE)   #nb records for each category
          datacaptive %>%
            group_by( host_status_2021) %>% 
            summarise(nb = sum(sample_size)/11300) #sample size by each category  
          #Let's remove NA and low values
          datastatus = subset(datacaptive, datacaptive$host_status_2021!="NA")
          datastatus = subset(datastatus, datastatus$host_status_2021!="DD")
          datastatus %>%
            dplyr::count(host_status_2021, sort=TRUE)   #nb records for each category
          datastatus %>%
            group_by( host_status_2021) %>% 
            summarise(nb = sum(sample_size)) #sample size by each category  
         #subgroup analysis
          model_subgroup<-metaprop(positive_sample, sample_size,data=datastatus,studlab=paper_ID,byvar=host_status_2021,level.comb =0.99, comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                                   sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
          summary(model_subgroup) #significant p= 0.0408 , more in near threatened 
          forest(model_subgroup,comb.random = T,  overall=T,subgroup=T,study.results=F, label.test.overall.random=T, test.overall	=T,pooled.events=T, pooled.totals=T, prediction=F )
          
      #3-E 5th subgroup analysis: host taxonomy order
            datacaptive %>%
              dplyr::count(host_order, sort=TRUE)   #nb records for each category
            host_sampling = datacaptive %>%
              group_by( host_order) %>% 
              summarise(nb = sum(sample_size)) #sample size by each category
            datacaptive2=host_sampling %>% filter(nb>10) %>% print(n=40)
            datacaptive2= datacaptive[(datacaptive$host_order=="Carnivora") | (datacaptive$host_order=="Diprotodontia")| (datacaptive$host_order=="Cetartiodactyla")| (datacaptive$host_order=="Lagomorpha")| (datacaptive$host_order=="Perissodactyla")| (datacaptive$host_order=="Primates")| (datacaptive$host_order=="Proboscidea")| (datacaptive$host_order=="Rodentia"), ]
            
          #subgroup analysis
          model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive2,studlab=paper_ID,byvar=host_order,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                                   sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
          summary(model_subgroup) #significantly p= 0.0001 variation infection according to host order (see forest plot).
          forest(model_subgroup,comb.random = T,  overall=T,subgroup=T,study.results=F, label.test.overall.random=T, test.overall	=T,pooled.events=T, pooled.totals=T, prediction=F )          
          #but investigation where Heterogeneity tests are significant should be investigated further (with virus variables and within families)
            #for: , , , perrissodactyla, proboscidea
          
      #3-F 6th subgroup analysis: host family
          #select from host_order results
         datacaptive %>%
            group_by( host_fam) %>% 
            summarise(nb = sum(sample_size))%>% print(n=70) #sample size by each category
          datacaptive3= datacaptive[(datacaptive$host_order=="Carnivora") | (datacaptive$host_order=="Cetartiodactyla")|(datacaptive$host_order=="Perissodactyla")| (datacaptive$host_order=="Primates")| (datacaptive$host_order=="Proboscidea"), ]
          host_sampling = datacaptive3 %>%
            group_by( host_fam) %>% 
            summarise(nb = sum(sample_size)) #sample size by each category
          host_sampling %>% filter(nb>30) %>% print(n=40)
          
          datacaptive4= datacaptive[(datacaptive$host_fam=="Atelidae") | (datacaptive$host_fam=="Bovidae")| (datacaptive$host_fam=="Callitrichidae")| (datacaptive$host_fam=="Camelidae")| (datacaptive$host_fam=="Canidae")| (datacaptive$host_fam=="Cebidae")| (datacaptive$host_fam=="Cercopithecidae")| (datacaptive$host_fam=="Cervidae")
                                    | (datacaptive$host_fam=="Elephantidae")| (datacaptive$host_fam=="Equidae")| (datacaptive$host_fam=="Felidae")| (datacaptive$host_fam=="Giraffidae")| (datacaptive$host_fam=="Herpestidae")| (datacaptive$host_fam=="Hominidae")| (datacaptive$host_fam=="Hylobatidae")
                                    | (datacaptive$host_fam=="Lemuridae")| (datacaptive$host_fam=="Pitheciidae")| (datacaptive$host_fam=="Procyonidae")| (datacaptive$host_fam=="Rhinocerotidae")| (datacaptive$host_fam=="Suidae")| (datacaptive$host_fam=="Ursidae"), ]
          #datacaptive4=datacaptive4 %>% drop_na()
          
          datacaptive4 %>%
            dplyr::count(host_fam, sort=TRUE)   #nb records for each category
          datacaptive4 %>%
            group_by( host_fam) %>% 
            summarise(nb = sum(sample_size)) %>% print(n=40) #sample size by each category  
          #subgroup analysis
          model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive4,studlab=paper_ID,byvar=host_fam, comb.random = TRUE, method.tau = "DL",prediction=TRUE,
                                   sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T, control = NULL)
          summary(model_subgroup) #significantly p= 0.0001 variation infection according to host family (see forest plot).
          forest(model_subgroup,comb.random = T,studlab=T,  overall=T,sortvar=host_order,subgroup=T,study.results=F, label.test.overall.random=F, test.overall	=T,pooled.events=T, pooled.totals=T, prediction=F,test.subgroup.random	=F )
          #highest in Cebidae, Cercopithecidae, Elephantidae, Equidae, Rhinocerotidae, Felidae, Atelidae
          #low in Ursidae, Tayassuidae, Procyonidae, Mustelidae, Lemuridae, Herpes, Giraffidae, Camelidae, Callithrichidae
          #but investigation where Heterogeneity tests are significant should be investigated further (with virus variables)
          
    #4- CAPTIVE ONLY: Subgroup analysis according to VIRUS VARIABLES
       #4-A 1st subgroup analysis: risk to human (identified zoonotic agents)
          host_paper_or=datacaptive %>%
            dplyr::count(risk_to_human, paper_ID, sort=TRUE)
          host_paper_or %>%
            dplyr::count(risk_to_human, sort=TRUE)   #nb paper for each category
          datacaptive %>%
            group_by( risk_to_human) %>% 
            summarise(nb = sum(sample_size)/11300) #sample size by each category   
          #subgroup analysis
          model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive,studlab=paper_ID,byvar=risk_to_human,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                                   sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
          summary(model_subgroup) # not significant p=0.2735
          forest(model_subgroup,comb.random = T,studlab=T,  overall=T,sortvar=host_order,subgroup=T,study.results=F, label.test.overall.random=F, test.overall	=T,pooled.events=T, pooled.totals=T, prediction=F,test.subgroup.random	=F )
          
       #4-B 2nd subgroup analysis: transmission type
          host_paper_or=datacaptive %>%
            dplyr::count(vector, paper_ID, sort=TRUE)
          host_paper_or %>%
            dplyr::count(vector, sort=TRUE)   #nb paper for each category
          datacaptive %>%
            group_by( vector) %>% 
            summarise(nb = sum(sample_size)/11300) #sample size by each category   
          #subgroup analysis
          model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive,studlab=paper_ID,byvar=vector,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                                   sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
          summary(model_subgroup) # significant p<0.001, higher infection when vector is through direct contact and unknown (=SV40)
          forest(model_subgroup,comb.random = T,  overall=T,subgroup=T,study.results=F, label.test.overall.random=T, test.overall	=T,pooled.events=T, pooled.totals=T, prediction=F )
          
        #4-C 3rd subgroup analysis: virus class (baltimore classification)
          host_paper_or=datacaptive %>%
            dplyr::count(v_class.x, paper_ID, sort=TRUE)
          host_paper_or %>%
            dplyr::count(v_class.x, sort=TRUE)   #nb paper for each category
          datacaptive %>%
            group_by( v_class.x) %>% 
            summarise(nb = sum(sample_size)) #sample size by each category   
          #subgroup analysis
          model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive,studlab=paper_ID,byvar=v_class.x,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                                   sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
          summary(model_subgroup) # significant p<0.001, higher infection of virus class VI = ssRNA(retroviruses) and I = dsDNA (herpes, SV40)
          forest(model_subgroup,comb.random = T,  overall=T,subgroup=T,study.results=F, label.test.overall.random=T, test.overall	=T,pooled.events=T, pooled.totals=T, prediction=F )
          
        #4-D  4th subgroup analysis: VIRUS species
          #only with significant variation in heterogeneity classes and sample size > 30
          vsp_sampling = datacaptive %>%
            group_by( v_sp) %>% 
            summarise(nb = sum(sample_size)) %>% print(n=40)#sample size by each category
          vsp_sampling %>% filter(nb>30) %>% print(n=40)
          
          datacaptive6= subset(datacaptive, v_sp!="foot and mouth disease virus")
          datacaptive6= subset(datacaptive6, v_sp!="caprine arthritis-encephalitis virus")
          datacaptive6= subset(datacaptive6, v_sp!="feline coronavirus")
          datacaptive6= subset(datacaptive6, v_sp!="maedi-visna virus")
          datacaptive6= subset(datacaptive6, v_sp!="peste des petites ruminants virus")
          datacaptive6= subset(datacaptive6, v_sp!="epizootic haemorrhagic disease virus")
          
          host_paper_or=datacaptive6 %>%
            dplyr::count(v_fam.x, paper_ID, sort=TRUE)
          host_paper_or %>%
            dplyr::count(v_fam.x, sort=TRUE)   #nb paper for each category
          datacaptive6 %>%
            group_by( v_class.x,v_fam.x, v_sp) %>% 
            summarise(nb = sum(sample_size)) %>% print(n=40)#sample size by each category   
          #subgroup analysis
              model_subgroup<-metaprop(positive_sample, sample_size,data=datacaptive6,studlab=paper_ID,byvar=v_sp,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                                       sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
              summary(model_subgroup) # significant p<0.001, higher infection of retro and herpers, and polyo (sv40)
              forest(model_subgroup,comb.random = T,  overall=T,subgroup=T,study.results=F, label.test.overall.random=T, test.overall	=T,pooled.events=T, pooled.totals=T, prediction=F )
              
    #5 CAPTIVE ONLY: host family and virus family interaction    
          
               #subgroup analysis at host order level
               model_subgroup<-metaprop(positive_sample, sample_size, data=datacaptive6, studlab=paper_ID,comb.random = TRUE, method.tau = "DL",prediction=TRUE,
                                        sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
               model_metareg= metareg(~ host_order:v_class.x, x = model_subgroup)
               summary(model_metareg)#Test moderator p<0.0001.
               
          
              #for carnivores and class IV
              datacaptive20= datacaptive[(datacaptive$host_order=="Carnivora" & datacaptive$v_class.x=="IV"), ]
              model_subgroup<-metaprop(positive_sample, sample_size, data=datacaptive20, studlab=paper_ID,comb.random = TRUE, method.tau = "DL",prediction=TRUE,
                                       sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
              model_metareg= metareg(~ host_fam:v_sp, x = model_subgroup)
              summary(model_metareg)#Test moderator p<0.001.
               
              #for class I
              datacaptive20= datacaptive[(datacaptive$v_class.x=="I"), ]
              model_subgroup<-metaprop(positive_sample, sample_size, data=datacaptive20, studlab=paper_ID,comb.random = TRUE, method.tau = "DL",prediction=TRUE,
                                       sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
              model_metareg= metareg(~ host_fam:v_sp, x = model_subgroup)
              summary(model_metareg)#Test moderator p<0.001.
              
              #for primates and class VI
              datacaptive20= datacaptive[(datacaptive$host_order=="Primates" & datacaptive$v_class.x=="VI"), ]
              model_subgroup<-metaprop(positive_sample, sample_size, data=datacaptive20, studlab=paper_ID,comb.random = TRUE, method.tau = "DL",prediction=TRUE,
                                       sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
              model_metareg= metareg(~ host_fam:v_sp, x = model_subgroup)
              summary(model_metareg)#Test moderator p<0.001.
              
              
    #6 CAPTIVE ONLY: meta-regression
      #now that we did subgroup analysis to evaluate potential colinearity, we can build a model with multiple moderators
      cor(x=datacaptive$pub_year, y=datacaptive$sample_size, method="pearson")
      #colinearity with length study and pub_year
      datacaptive$length_study = datacaptive$collection_lastyear-datacaptive$collection_firstyear
      datacaptive1=subset(datacaptive, datacaptive$length_study!="NA") #not taking NA values 
      #build general model without moderators, using full sampling
      model_subgroup<-metaprop(positive_sample, sample_size, data=datacaptive1, studlab=paper_ID,comb.random = TRUE, method.tau = "REML",prediction=TRUE,
                               sm="PFT",comb.fixed = F, overall = T,tau.common=FALSE, method.bias="Egger", backtransf=T,overall.hetstat=T, hakn=T)
      #we already looked at virus x host taxonomy interaction so if we only look at the others
      model_metareg1= metareg(~ zoo_loc2:same_continent+vector+sample_size_class+length_study, x = model_subgroup) #testing moderators
      summary(model_metareg1) #AIC=476, i^2 = 71%, r^2= 10%, p<0.0001
      model_metareg2= metareg(~ zoo_loc2:same_continent+vector+risk_to_human+host_status_2021+zoo_type+number_faciliti+sample_size_class+length_study, x = model_subgroup) #testing moderators
      summary(model_metareg2)
      
      
      #heat maps
      #sample size
          host_v=datacaptive %>%
            group_by(host_order, host_fam,v_class.x,v_fam.x, v_sp) %>% 
            summarise(nb = sum(sample_size))
          host_v=host_v %>% filter(nb>0)
          host_v=host_v %>% filter(host_fam != "NA")
          x = c("I", "II", "III", "IV", "V", "VI", "VII")
          host_v %>%
            mutate(host_fam = fct_reorder(host_fam, desc(nb))) %>% mutate(v_class.x =  factor(v_class.x, levels = x)) %>%
            arrange(v_class.x,v_fam.x) %>%
            ggplot( aes(x = host_fam, y = fct_reorder(v_sp, desc(fct_reorder(v_fam.x, desc(v_class.x)))), fill = nb, na.rm = TRUE)) + theme_bw()+scale_fill_gradient2(low = "#71BF50", high = "#A43020", mid = "#F4CE4B",  midpoint = 30, limit = c(0,70),
                                                                                                                                                                      name="Sample size")+
            geom_tile() +ylab("Virus taxonomy")+xlab("Host taxonomy")+ theme(axis.text.x=element_text(angle = 90, hjust = 0))
      #seroprevalence
          host_v=datacaptive %>%
            dplyr::group_by(host_order,host_fam, v_class.x,v_fam.x, v_sp) %>% 
            dplyr::summarise(nb = mean(positive_sample/sample_size))
          x = c("I", "II", "III", "IV", "V", "VI", "VII")
          host_v %>%
            mutate(host_fam = fct_reorder(host_fam, desc(nb))) %>% mutate(v_class.x =  factor(v_class.x, levels = x)) %>%
            arrange(v_class.x,v_fam.x) %>%
            ggplot( aes(x = fct_reorder(host_fam, desc(fct_reorder(host_order, desc(host_order)))), y = fct_reorder(v_sp, desc(fct_reorder(v_fam.x, desc(v_class.x)))), fill = nb, na.rm = TRUE)) + theme_bw()+scale_fill_gradient2(low = "#71BF50", high = "#A43020", mid = "#F4CE4B",  midpoint = 0.5, limit = c(0,1),
                                                                                                                                                                      name="Mean seroprevalence frequency")+
            geom_tile() +ylab("Virus taxonomy")+xlab("Host taxonomy")+ theme(axis.text.x=element_text(angle = 90, hjust = 0))

      # particular plots
          
          library(Rmisc)
          datacaptiveI= datacaptive[(datacaptive$v_class.x=="I"), ]
          seroplotI=datacaptiveI %>%
            dplyr::group_by(host_order, host_fam, v_sp) %>% 
            dplyr::summarise(mean = (positive_sample/sample_size))  
          seroplotI
          seroplotI <- summarySE(seroplotI, measurevar="mean", groupvars=c("host_order", "host_fam", "v_sp"))
          seroplotI
          plotI = seroplotI %>% drop_na() %>% filter(mean>0)%>%
            ggplot(aes(x=host_fam, y=mean, group=v_sp))+ ggtitle("Host taxonomy") + geom_bar(stat="identity",position=position_dodge(), aes(fill = v_sp)) +ylab("Mean seroprevalence +SE") +
            theme_bw() +  facet_grid(. ~ host_order, scales = "free")+theme(axis.text.x=element_text(angle = 0, hjust = 1))+
            geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +ylim(0,1.1)
          plotI
          
          datacaptiveIII= datacaptive[(datacaptive$v_class.x=="III" | datacaptive$v_sp=="schmallenberg virus" | datacaptive$v_sp=="west nile virus"), ]
          seroplotIII=datacaptiveIII %>%
            dplyr::group_by(host_order, host_fam, v_sp) %>% 
            dplyr::summarise(mean = (positive_sample/sample_size))  
          seroplotIII
          seroplotIII <- summarySE(seroplotIII, measurevar="mean", groupvars=c("host_order", "host_fam", "v_sp"))
          seroplotIII
          plotIII = seroplotIII %>% drop_na() %>% filter(mean>0)%>%
            ggplot(aes(x=fct_reorder(host_fam, desc(fct_reorder(host_order, desc(host_order)))), y=mean, group=v_sp))+ geom_bar(stat="identity",position=position_dodge(), aes(fill = v_sp)) +ylab("Mean seroprevalence +SE") +
            theme_bw() +theme(axis.text.x=element_text(angle = 45, hjust = 1))+facet_grid(. ~ v_sp, scales = "free")+
            geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +ylim(0,1.1)
          plotIII
          
          datacaptiveVI= datacaptive[(datacaptive$v_class.x=="VI" &datacaptive$host_order=="Primates"), ]
          seroplotVI=datacaptiveVI %>%
            dplyr::group_by(host_order, host_fam, v_sp) %>% 
            dplyr::summarise(mean = (positive_sample/sample_size))  
          seroplotVI
          seroplotVI <- summarySE(seroplotVI, measurevar="mean", groupvars=c("host_order", "host_fam", "v_sp"))
          seroplotVI
          plotVI = seroplotVI %>% filter(mean>0)%>%
            ggplot(aes(x=fct_reorder(host_fam, desc(fct_reorder(host_order, desc(host_order)))), y=mean, group=v_sp))+ geom_bar(stat="identity",position=position_dodge(), aes(fill = v_sp)) +ylab("Mean seroprevalence +SE") +
            theme_bw() +theme(axis.text.x=element_text(angle = 45, hjust = 1))+facet_grid(. ~ v_sp, scales = "free")+
            geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +ylim(0,1.1)
          plotVI
          
          library("ggpubr")
          ggarrange(plotI, plotIII, plotVI , common.legend = T,
                    labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)
          
          
          