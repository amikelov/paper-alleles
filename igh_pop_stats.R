#Load libraries, set priorities ----
{
  library(rlang)
  library(ComplexHeatmap)
  library(forcats)
  library(treemapify)
  library(scales)
  library(circlize)
  library(conflicted)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(scales)
  library(ggpubr)
  library(forcats)
  library(purrr)
  library(Biostrings)
  library(cowplot)
  library(ggpmisc)
  conflict_prefer("set_names", winner = "magrittr", quiet = FALSE)
  conflict_prefer("filter", winner = "dplyr", quiet = FALSE)
  conflict_prefer("first", winner = "dplyr", quiet = FALSE)
  conflict_prefer("rename", winner = "dplyr", quiet = FALSE)
  conflict_prefer("count", winner = "dplyr", quiet = FALSE)
  conflict_prefer("unique", winner = "base", quiet = FALSE)
  conflicts_prefer(ggplot2::annotate)
  conflicts_prefer(dplyr::desc)
  options(dplyr.summarise.inform = FALSE)
  `%+%` <- function(a, b) paste0(a, b)
}

#Load MiXCR reports generated with exportReportsTable
mixcr_stats<-read_tsv("/home/amikelov/paper_alleles/mixcrReports.tsv",) %>%
  mutate(donor=str_remove(fileName,"SRR.*_") %>% 
           str_remove_all("clones\\.|.*\\/|_.*|\\..*"),
         .after=1)  %>% 
  relocate(fileName,donor,totalClonotypes,readsUsedInClonotypes,.before="MiXCRVersion")

# Filter out donors with low number of clonotypes ----
goodDepthDonors<-mixcr_stats %>%
  group_by(donor) %>% 
  summarise(readsUsedInClonotypes=sum(readsUsedInClonotypes),
            totalReads=sum(totalReads),
            totalClonotypes=sum(totalClonotypes)) %>% 
  filter(totalClonotypes>3000) %>% pull(donor)

# Load findAlleles tsv outputs to determine significant vGenes ----
findAllelesOutput<-sapply(goodDepthDonors, function(don){
  c(dir("/home/amikelov/paper_alleles/boyd_fastq/alleles", pattern = don %+% "\\..*tsv",full.names = T),
    dir("/home/amikelov/paper_alleles/boyd/alleles", pattern = don %+% "\\..*tsv",full.names = T),
    dir("/home/amikelov/paper_alleles/watson/alleles", pattern = don %+% "\\..*tsv",full.names = T),
    dir("/home/amikelov/paper_alleles/PRJEB26509_naiveIGH/alleles", pattern = don %+% "\\..*tsv",full.names = T)
  )
})

findAllelesOutput<-read_tsv(findAllelesOutput,id="fileName") %>% 
  mutate(donor=str_remove_all(fileName,".*\\/|.tsv"),.after = fileName) %>% 
  filter(str_detect(geneName,'IGHV|IGHJ')) %>% 
  select(donor:totalClonesCountForGene,filteredForAlleleSearchClonesCount)

bigVGenes<-findAllelesOutput %>% 
  filter(str_detect(geneName,'IGHV')) %>% 
  group_by(geneName) %>% 
  summarise(meanTotCount=mean(totalClonesCountForGene),
            maxTotCount=max(totalClonesCountForGene)) %>% 
  filter(maxTotCount>500) %>% 
  pull(geneName)

#Load metadata from series of addNewAlleles commands ----
alleles_db <- read_tsv("/home/amikelov/paper_alleles/mergedLib/meta_popStats.tsv") %>% 
  mutate(geneName=str_remove(alleleName,"\\*.*")) %>% 
  filter(
    str_detect(geneName,'IGHV|IGHJ')) %>% 
  left_join(findAllelesOutput %>% 
              select(donor,alleleName,clonesCount,totalClonesCountForGene),
            join_by(donor==donor, alleleNameInSource==alleleName)) %>% 
  mutate(ethnicity=str_replace(ethnicity,"hispanic","Hispanic/Latino") %>% str_to_title()) 


########ALLELE STATISTICS IN POPULATIONS################

ethnicity_colors<-c(African="#A2F66C",
                    Asian="#FFF780",
                    Caucasian="#5DDDA4",
                    `Hispanic/Latino`="#942AAE",
                    Unknown="#D3D7E0",
                    All="#5EA8DF"
)




cohortDemographics<-alleles_db %>% 
  group_by(ethnicity) %>% 
  summarize(nDonors=n_distinct(donor)) %>%
  mutate(nDonorsPercent=percent(nDonors/sum(nDonors),accuracy= 0.1) )

# Cohort by ethnicity ----
## Fig 3A N donors p/ ethnicity ----
# fig3a<-cohortDemographics %>% 
#   #mutate(ethnicity=factor(ethnicity,levels= ))
#   ggplot(aes(area = nDonors, fill = ethnicity,label = ethnicity %+% "\n" %+% "N=" %+% nDonors %+% "\n" %+%nDonorsPercent)) +
#     geom_treemap(alpha=0.7)+
#     geom_treemap_text()+
#   theme(legend.position = "none")+
#   scale_fill_manual(values = ethnicity_colors)+
#   theme(plot.margin = unit(c(0.4,0.6,1,1), "cm"))

fig3a<-waffle::waffle(cohortDemographics %>%
                        pull(nDonors) %>% 
                        set_names(c("African, N=112",
                                    "Asian, N=18",
                                    "Caucasian, N=92",
                                    "Hispanic/Latino, N=6",
                                    "Unknown, N=222")),
                      rows=13,
                      colors=c("#A2F66C","#FFF780","#5DDDA4",
                               "#942AAE","#D3D7E0"))+
  theme(legend.text=element_text(size=14),
        plot.margin = unit(c(0.4,0.6,1,1), "cm"))


fig3a


#Stats by Vgenes ----
vGeneStats<-alleles_db %>%
  filter( str_detect(geneName,"IGHV") ) %>%
  group_by(donor) %>%
  dplyr::mutate(totalClonesInDonor=sum(clonesCount),
                vGeneFrequency=totalClonesCountForGene/totalClonesInDonor) %>%
  ungroup() %>%
  filter(totalClonesCountForGene>50) %>%
  group_by(geneName) %>%
  summarise(nAlleles=n_distinct(alleleName),
            nDonors=n_distinct(donor),
            ethnicities=paste(unique(ethnicity),collapse=","),
            medianFrequency=median(vGeneFrequency)) %>%
  mutate(rank=dense_rank(dplyr::desc(medianFrequency)),
         donorFraction=nDonors/length(goodDepthDonors))

##Limit analysis to frequent V-genes found more than in 5% of donors----
frequentVs<-vGeneStats %>%
  filter(donorFraction>0.05) %>%
  pull(geneName)

##View stats for filtered out vGenes
vGeneStats %>%
  filter(donorFraction<0.05) 

#Stats by Jgenes ----
jGeneStats<-alleles_db %>%
  filter(str_detect(geneName,"IGHJ") ) %>%
  group_by(donor) %>%
  mutate(totalClonesInDonor=sum(clonesCount),
         jGeneFrequency=totalClonesCountForGene/totalClonesInDonor) %>%
  ungroup() %>%
  filter(totalClonesCountForGene>50
  ) %>%
  group_by(geneName) %>%
  summarise(nAlleles=n_distinct(alleleName),
            nDonors=n_distinct(donor),
            ethnicities=paste(unique(ethnicity),collapse=","),
            medianFrequency=median(jGeneFrequency)) %>%
  mutate(rank=dense_rank(dplyr::desc(medianFrequency)),
         donorFraction=nDonors/length(goodDepthDonors))


#N V alleles in populations----
## Fig S3 Total N V alleles p/ Gene & Ethnicity
# figS3<- alleles_db %>%
#   filter(totalClonesCountForGene>50,
#          ethnicity!="Unknown",
#          geneName %in% frequentVs) %>% 
#   group_by(geneName,ethnicity) %>% 
#   summarise(nAlleles=n_distinct(alleleName),
#             nDonors=n_distinct(donor)) %>%
#   bind_rows(
#     alleles_db %>%
#       filter(totalClonesCountForGene>50,
#              geneName %in% frequentVs) %>% 
#       group_by(geneName) %>% 
#       summarise(nAlleles=n_distinct(alleleName),
#                 nDonors=n_distinct(donor)) %>% 
#       mutate(ethnicity="All")
#   ) %>% 
#   inner_join(vGeneStats %>% select(geneName,medianFrequency,rank),by="geneName") %>% 
#   mutate(vFreqCategory=factor(case_when(
#     rank >= 1 & rank < 15 ~ "1",
#     rank >= 15 & rank < 30 ~ "2",
#     rank >= 30 ~ "3"),
#     levels=c("1","2","3"))) %>% 
#   mutate(Ethnicity=factor(ethnicity,levels=c("All","Caucasian","African","Asian", "Hispanic/Latino"))) %>%  
#   ggplot(aes(x=fct_reorder(geneName,rank),y=nAlleles,fill=Ethnicity))+
#   geom_bar(stat="identity",
#            position = position_dodge2(width = 0.9, preserve = "single"),
#            color="grey38",linewidth=0.15)+
#   theme_bw()+
#   rotate_x_text()+
#   facet_wrap(~vFreqCategory,scales="free_x",nrow=3)+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     legend.position = "none"
#   )+
#   labs(x="V gene name", y="# alleles")+
#   scale_fill_manual(values=ethnicity_colors)

# ggsave(filename = "~/alleles/alleles-paper/plots/NumberOfAlleles.pdf", plot= g,device = "pdf",
#        width = 12,height = 8,units = "in",bg="white")
# ggsave(filename = "~/alleles/alleles-paper/plots/NumberOfAlleles.png", plot= g,device = "png",
#        width = 12,height = 8,units = "in",bg="white")


#N J alleles ----
## Fig S4 N J alleles p/ Gene & Ethnicity----
figS4<- alleles_db %>%
  filter(totalClonesCountForGene>50,
         ethnicity!="Unknown",
         str_detect(geneName,"IGHJ")) %>% 
  group_by(geneName,ethnicity) %>% 
  summarise(nAlleles=n_distinct(alleleName),
            nDonors=n_distinct(donor)) %>%
  bind_rows(
    alleles_db %>%
      filter(totalClonesCountForGene>50,
             str_detect(geneName,"IGHJ")) %>% 
      group_by(geneName) %>% 
      summarise(nAlleles=n_distinct(alleleName),
                nDonors=n_distinct(donor)) %>% 
      mutate(ethnicity="All")
  ) %>% 
  inner_join(jGeneStats %>% select(geneName,medianFrequency,rank),by="geneName") %>% 
  mutate(Ethnicity=factor(ethnicity,levels=c("All","Caucasian","African","Asian", "Hispanic/Latino"))
  ) %>% 
  ggplot(aes(x=geneName,y=nAlleles,fill=Ethnicity))+
  geom_bar(stat="identity",
           position = position_dodge2(width = 0.9, preserve = "single"),
           color="grey38",linewidth=0.15)+
  theme_bw()+
  rotate_x_text()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = 'none'
  )+
  labs(x="J gene name", y="# alleles")+
  scale_fill_manual(values=ethnicity_colors)

# ggsave(filename = "~/alleles/alleles-paper/plots/NumberOfAllelesJ.pdf", plot= fig2c,device = "pdf",
#        width = 12,height = 8,units = "in",bg="white")
# ggsave(filename = "~/alleles/alleles-paper/plots/NumberOfAllelesJ.png", plot= fig2c,device = "png",
#        width = 12,height = 8,units = "in",bg="white")

#  N novel alleles ----
novelAlleles<-scan("~/paper_alleles/mergedLib/novelAlleles.tsv", character(), quote = "") # file from comparing resulting and original json libs, see upstream

allelesStats_igh<-alleles_db %>% 
  filter(totalClonesCountForGene>50) %>% 
  mutate(nDonors=n_distinct(donor),
         isNovel=ifelse(alleleName %in% novelAlleles,"Novel","Known") %>% factor(levels=c("Novel","Known")),
         geneType=ifelse(str_detect(alleleName,"V"),"V","J") %>% factor(levels=c("V","J"))) %>% 
  group_by(geneName,geneType,isNovel) %>% 
  summarise(nAlleles=n_distinct(alleleName),
            nAllelesPerDonors=nAlleles/first(nDonors),
            nDonors=first(nDonors)) 

## Fig 2C  IGHV N alleles ----
fig2c<-allelesStats_igh %>% 
  inner_join(vGeneStats %>% select(geneName,medianFrequency,rank),by="geneName") %>% 
  ggplot(aes(x=fct_reorder(geneName,desc(rank)),y=nAlleles,fill=isNovel,label=nAlleles))+
  geom_bar(stat="identity")+
  facet_wrap(~geneType,scales="free",
             labeller = labeller(geneType=c(V="IGH, V genes",J="IGH, J genes")))+
  theme_bw()+
  coord_flip()+
  scale_fill_manual(values=c(Novel="#2D93FA",Known="#99CCFF"))+
  theme(strip.background =element_rect(fill="white"),
        legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.title.y = element_blank())+
  labs(x= "Genes", y="# Detected alleles")+
  geom_text(aes(color = isNovel),
            size = 3, 
            position = position_stack(vjust = 0.5))+
  scale_color_manual(values=c("white","black"))+
  scale_y_continuous(expand = expansion(mult = c(0.005, .1))) +
  annotate('label',x=-Inf,y=Inf,hjust=1,vjust=0,
           label= length(unique(alleles_db$donor)) %+% " individuals" )



## Fig 2F  IGHJ N alleles ----
fig2f<-allelesStats_igh %>% 
  inner_join(jGeneStats %>% select(geneName,medianFrequency,rank),by="geneName") %>% 
  ggplot(aes(x=fct_reorder(geneName,desc(rank)),y=nAlleles,fill=isNovel))+
  geom_bar(stat="identity")+
  facet_wrap(~geneType,scales="free",
             labeller = labeller(geneType=c(V="IGH, V genes",J="IGH, J genes")))+
  theme_bw()+
  coord_flip()+
  scale_fill_manual(values=c(Novel="#2D93FA",Known="#99CCFF"))+
  theme(strip.background =element_rect(fill="white"),
        legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank())+
  labs(x= "Genes", y="# Detected alleles")+
  geom_text(aes(color = isNovel,label=nAlleles),
            size = 3, 
            position = position_stack(vjust = 0.5))+
  scale_color_manual(values=c("white","black"))+
  scale_y_continuous(expand = expansion(mult = c(0.005, .1))) +
  annotate('label',x=-Inf,y=Inf,hjust=1,vjust=0, 
           label= length(unique(alleles_db$donor)) %+% " ind." )


# N novel alleles per Ethnic group----
novelAllelesStats<-alleles_db %>%
  group_by(ethnicity) %>% 
  mutate(nDonorsEthnicity=n_distinct(donor),
         geneType=ifelse(str_detect(alleleName,"IGHV"),"V","J") %>% factor(levels=c("V","J"))) %>% 
  filter(alleleName %in% novelAlleles ) %>%
  group_by(ethnicity,geneType) %>% 
  summarise(nNovel=n_distinct(alleleName),
            nNovelPerDonors=nNovel/first(nDonorsEthnicity)) %>% 
  mutate(ypos=ifelse(geneType=="V",-0.1,-0.002))


#N novel alleles permutation test ----
donors_ethnicities<-alleles_db |> 
  select(donor,ethnicity) |> 
  unique()

ethnicities <- alleles_db |> 
  pull(ethnicity) |> 
  unique() 

novelPermTest<-lapply(1:10000, function(it) {
  ethnicities_reshuffled <- donors_ethnicities |> mutate(ethnicity=sample(ethnicity))
  
  novelPermutedDb <- alleles_db %>%
    select(-ethnicity) |> 
    left_join(ethnicities_reshuffled, by= "donor") |> 
    group_by(ethnicity) %>% 
    mutate(nDonorsEthnicity=n_distinct(donor),
           geneType=ifelse(str_detect(alleleName,"IGHV"),"V","J") %>% factor(levels=c("V","J"))) %>% 
    filter(alleleName %in% novelAlleles ) %>%
    group_by(ethnicity,geneType) %>% 
    summarise(nNovel=n_distinct(alleleName),
              nNovelPerDonors=nNovel/first(nDonorsEthnicity)) 
  
  crossing(ethn1=ethnicities,ethn2=ethnicities) |> 
    filter(ethn1 < ethn2) |> 
    left_join(novelPermutedDb, by= c( "ethn1" = "ethnicity") ) |> 
    left_join(novelPermutedDb, by= c( "ethn2" = "ethnicity") ) |> 
    replace_na(replace = list(nNovelPerDonors.x=0, nNovelPerDonors.y=0)) |> 
    mutate(nNovelPerDonorsDiff=abs(nNovelPerDonors.x-nNovelPerDonors.y),
           iteration = it)
  
}) |> bind_rows()

novelPermTest |> 
  left_join(
    crossing(ethn1=ethnicities,ethn2=ethnicities) |> 
      filter(ethn1 < ethn2) |> 
      left_join(novelAllelesStats, by= c( "ethn1" = "ethnicity") ) |> 
      left_join(novelAllelesStats, by= c( "ethn2" = "ethnicity") ) |> 
      mutate(nNovelPerDonorsDiff=abs(nNovelPerDonors.x-nNovelPerDonors.y)),
    by=c("ethn1","ethn2")) |> 
  mutate(diffMoreInPermuted=nNovelPerDonorsDiff.x>=nNovelPerDonorsDiff.y) |> 
  group_by(ethn1,ethn2) |> 
  summarise(p=sum(diffMoreInPermuted/max(iteration)))




## Fig 3B N Novel p/ Ethnicity ----
fig3b<-novelAllelesStats %>% 
  mutate(ethnicity=factor(ethnicity,
                          levels=c("Unknown",
                                   "Caucasian",
                                   "African",
                                   "Asian", 
                                   "Hispanic/Latino"))
  ) %>% 
  ggplot(aes(x=ethnicity,fill=ethnicity,y=nNovelPerDonors))+
  geom_bar(stat="identity",position = position_dodge(),
           color="grey38",linewidth=0.15)+
  theme_bw()+
  scale_fill_manual(values=ethnicity_colors)+
  scale_x_discrete(labels=c(Unknown="Unknown",
                            Caucasian="Caucas.",
                            African="African",
                            Asian="Asian", 
                            `Hispanic/Latino`="Hisp./Latino"))+
  facet_wrap(~geneType,
             nrow=1,
             labeller=labeller(geneType=c(V="IGHV Genes",J="IGHJ Genes")))+
  labs(x= "Ethnicity",y="# Novel Alleles Per Donor")+
  theme(legend.position = "none",
        strip.background = element_rect(fill="white"),
        axis.title.x = element_blank(),
        plot.margin = unit(c(3,5.5,6.5,6.5), "points"),
        axis.title.y = element_text(size=14),     # Increase y-axis title font size
        axis.text.x = element_text(size=12),      # Increase x-axis labels font size
        axis.text.y = element_text(size=12),      # Increase y-axis labels font size
        strip.text = element_text(size=14))       # Increase facet labels font size





#N alleles Per donor for all Genes ----
totalPerDonorPerGeneAllelesStats<-alleles_db %>%
  group_by(ethnicity) %>% 
  mutate(nDonorsEthnicity=n_distinct(donor),
         geneType=ifelse(str_detect(alleleName,"IGHV"),"V","J") %>% factor(levels=c("V","J"))) %>% 
  group_by(ethnicity,geneName) %>% 
  summarise(nAlleles=n_distinct(alleleName),
            nAllelesPerDonor=nAlleles/first(nDonorsEthnicity)) 

#N alleles per Donor per gene permutation test ----
genes<-alleles_db %>%
  filter(totalClonesCountForGene>50,
         (geneName %in% frequentVs) | str_detect(alleleName,"IGHJ")) |> 
  pull(geneName) |> 
  unique()
  
totalPerDonorPerGenePermTest<-lapply(1:100, function(it) {
  ethnicities_reshuffled <- donors_ethnicities |> mutate(ethnicity=sample(ethnicity))
  
  totalPermutedDb <- alleles_db %>%
    select(-ethnicity) |> 
    left_join(ethnicities_reshuffled, by= "donor") |> 
    group_by(ethnicity) %>% 
    mutate(nDonorsEthnicity=n_distinct(donor),
           geneType=ifelse(str_detect(alleleName,"IGHV"),"V","J") %>% factor(levels=c("V","J"))) %>% 
    group_by(ethnicity,geneName) %>% 
    summarise(nAlleles=n_distinct(alleleName),
              nAllelesPerDonor=nAlleles/first(nDonorsEthnicity)) 
  
  crossing(ethn1=ethnicities,ethn2=ethnicities,geneName=genes) |> 
    filter(ethn1 < ethn2) |> 
    left_join(totalPermutedDb, by= c( "ethn1" = "ethnicity", "geneName") ) |> 
    left_join(totalPermutedDb, by= c( "ethn2" = "ethnicity", "geneName") ) |> 
    replace_na(replace = list(nAllelesPerDonor.x=0, nAllelesPerDonor.y=0)) |> 
    mutate(nTotalPerDonorsDiff=abs(nAllelesPerDonor.x-nAllelesPerDonor.y),
           iteration = it)
  
}) |> bind_rows()



totalPerDonorPerGenePermTest |> 
  left_join(
    crossing(ethn1=ethnicities,ethn2=ethnicities,geneName=genes) |> 
      filter(ethn1 < ethn2) |> 
      left_join(totalPerDonorPerGeneAllelesStats, by= c( "ethn1" = "ethnicity", "geneName") ) |> 
      left_join(totalPerDonorPerGeneAllelesStats, by= c( "ethn2" = "ethnicity", "geneName") ) |> 
      replace_na(replace = list(nAllelesPerDonor.x=0, nAllelesPerDonor.y=0)) |> 
      mutate(nTotalPerDonorsDiff=abs(nAllelesPerDonor.x-nAllelesPerDonor.y)),
    by=c("ethn1","ethn2", "geneName")) |> 
  mutate(diffMoreInPermuted=nTotalPerDonorsDiff.x>=nTotalPerDonorsDiff.y) |> 
  group_by(ethn1,ethn2,geneName) |> 
  summarise(p=sum(diffMoreInPermuted/max(iteration))) |> 
  ungroup() |> 
  mutate(p.adj=p.adjust(p,"fdr")) |> View()




#Downsample for V ----
alleleStats<-alleles_db %>%
  filter(totalClonesCountForGene>50,
         # ethnicity != "Unknown",
         ethnicity %in% c("Caucasian","African"),
         geneName %in% frequentVs) %>%
  bind_rows(
    alleles_db %>%
      filter(totalClonesCountForGene>50,
             geneName %in% frequentVs) %>%
      mutate(ethnicity="All")
  ) %>%
  nest(.by = c(ethnicity,geneName)) %>%
  #mutate(nDonors=map_int(data,~ .x %>% summarize(nDonor=n_distinct(donor)) %>% pull(nDonor))) %>%
  #group_by(geneName)%>%
  #now
  mutate(minNDonorsForGene=92) %>%
  crossing(iteration=1:10 ) %>%
  mutate(nAlleles=map2(data,minNDonorsForGene,~ .x %>%
                         filter(donor %in% (.x %>%
                                              select(donor) %>%
                                              unique() %>%
                                              slice_sample(n=.y) %>%
                                              pull(donor))) %>%
                         summarize(nAlleles=n_distinct(alleleName)) ) ) %>%
  select(-data)%>%
  unnest() %>%
  group_by(ethnicity,geneName) %>%
  summarise(nAlleles=round(mean(nAlleles)),minNDonorsForGene=dplyr::first(minNDonorsForGene))


# N alleles permutation test ----
permutationTest<-
alleles_db %>%
  filter(totalClonesCountForGene>50,
         # ethnicity != "Unknown",
         ethnicity %in% c("Caucasian","African"),
         geneName %in% frequentVs) %>%
  crossing(iteration=1:1000 ) %>%
  group_by(iteration,geneName) |>
  mutate(ethnicity = sample(ethnicity)) |>
  ungroup() |>
  nest(.by = c(ethnicity,geneName,iteration)) %>%
  mutate(minNDonorsForGene=92) %>%
  mutate(nAlleles=map2(data,minNDonorsForGene,~ .x %>%
                         filter(donor %in% (.x %>%
                                              select(donor) %>%
                                              unique() %>%
                                              slice_sample(n=.y) %>%
                                              pull(donor))) %>%
                         summarize(nAlleles=n_distinct(alleleName)) ) ) %>%
  select(-data)%>%
  unnest() %>%
  group_by(geneName,iteration) |>
  summarise(nAllelesDifAbs.permut=max(nAlleles)-min(nAlleles)) |>
  left_join(alleleStats |>
              filter(ethnicity %in% c("Caucasian","African")) |>
              group_by(geneName) |>
              summarise(nAllelesDifAbs.real=max(nAlleles)-min(nAlleles)),by="geneName")  |>
  mutate(permutDifMore=nAllelesDifAbs.permut>=nAllelesDifAbs.real) |>
  ungroup() |>
  group_by(geneName) |>
  summarise(nWhenDifMore=sum(permutDifMore),
            p.value=nWhenDifMore/max(iteration))

permutationTest <- permutationTest |> mutate(p.adj=p.adjust(p.value,"fdr"))

permutationTest |> write_tsv("~/paper_alleles/paper-alleles/permutationTest1K.tsv")

#N alleles per ethnic group ----
#TODO add permutation test
##Fig 3D IGHV N alleles Ethnic ----
fig3d <- alleleStats %>% 
  inner_join(vGeneStats %>% select(geneName,medianFrequency,rank),by="geneName") %>% 
  mutate(Ethnicity=factor(ethnicity,levels=c("All","Caucasian","African","Asian", "Hispanic/Latino")),
         vFreqCategory=factor(case_when(
           rank >= 1 & rank < 15 ~ "1",
           rank >= 15 & rank < 30 ~ "2",
           rank >= 30 ~ "3"),
           levels=c("1","2","3"))) %>% 
  ggplot(aes(x=fct_reorder(geneName,rank),y=nAlleles ,fill=Ethnicity))+
  geom_bar(stat="identity",position = position_dodge2(width = 0.9, preserve = "single"),
           color="grey38",linewidth=0.15)+
  #stat_summary(fun.data = give.n, geom = "text",size=3)+
  theme_bw()+
  rotate_x_text()+
  facet_wrap(~vFreqCategory,scales="free_x",nrow=3)+
  scale_fill_manual(values=ethnicity_colors)+  
  labs(x="V gene name", y="# Alleles (normalized)")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    plot.margin = unit(c(5.5,5.5,5.5,10),"points")
  ) + scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) 

# ggsave(filename = "~/alleles/alleles-paper/plots/NumberOfAllelesPerFixedNDonors.pdf", plot= g,device = "pdf",
#        width = 12,height = 8,units = "in",bg="white")
# ggsave(filename = "~/alleles/alleles-paper/plots/NumberOfAllelesPerFixedNDonors.png", plot= g,device = "png",
#        width = 12,height = 8,units = "in",bg="white")

#Downsample for J ----
alleleStats_J<-alleles_db %>%
  filter(totalClonesCountForGene>50,
         #ethnicity != "Unknown",
         ethnicity %in% c("Caucasian","African"),
         str_detect(alleleName,"IGHJ")) %>%
  bind_rows(
    alleles_db %>%
      filter(totalClonesCountForGene>50,
             str_detect(alleleName,"IGHJ")) %>%
      mutate(ethnicity="All")
  ) %>%  
  nest(.by = c(ethnicity,geneName)) %>% 
  #mutate(nDonors=map_int(data,~ .x %>% summarize(nDonor=n_distinct(donor)) %>% pull(nDonor))) %>% 
  #group_by(geneName)%>%
  #now 
  mutate(minNDonorsForGene=60) %>%
  crossing(iteration=1:10 ) %>%
  mutate(nAlleles=map2(data,minNDonorsForGene,~ .x %>% 
                         filter(donor %in% (.x %>% 
                                              select(donor) %>%
                                              unique() %>% 
                                              slice_sample(n=.y) %>%
                                              pull(donor))) %>% 
                         summarize(nAlleles=n_distinct(alleleName)) ) ) %>% 
  select(-data)%>% 
  unnest() %>% 
  group_by(ethnicity,geneName) %>% 
  summarise(nAlleles=round(mean(nAlleles)),minNDonorsForGene=dplyr::first(minNDonorsForGene)) 
##Fig 3C IGHJ N alleles Ethnic ----
fig3c <- alleleStats_J %>% 
  inner_join(jGeneStats %>% select(geneName,medianFrequency,rank),by="geneName") %>% 
  mutate(Ethnicity=factor(ethnicity,levels=c("All","Caucasian","African","Asian", "Hispanic/Latino"))) %>% 
  ggplot(aes(x=geneName,y=nAlleles,fill=Ethnicity))+
  geom_bar(stat="identity",position = position_dodge2(width = 0.9, preserve = "single"),
           color="grey38",linewidth=0.15)+
  
  #stat_summary(fun.data = give.n, geom = "text",size=3)+
  theme_bw()+
  rotate_x_text()+
  scale_fill_manual(values=ethnicity_colors)+  
  labs(x="J gene name", y="# Alleles (normalized)")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.y = element_text(size=14),     # Increase y-axis title font size
    axis.text.x = element_text(size=12),      # Increase x-axis labels font size
    axis.text.y = element_text(size=12),      # Increase y-axis labels font size
    strip.text = element_text(size=14),
    legend.position = 'none',
    axis.title.x = element_blank(),
    plot.margin = unit(c(10,5.5,10,6.5),"points")
  )

############# FREQUENCIES OF ALLELES ################

alleleFasta<-"/home/amikelov/paper_alleles/mergedLib/IGH_alleles.popEnriched.fasta"

donorAlleleSeq<-readDNAStringSet(alleleFasta)  
alleleSeq <- tibble(alleleName=names(donorAlleleSeq) %>% str_remove("\\|.*| .*"),
                    alleleSeq=as.character(donorAlleleSeq) %>% str_remove_all("\\."),
                    alleleSeqAA=as.character(translate(donorAlleleSeq)) ) %>% 
  filter(str_detect(alleleName,"IGHV")) %>% 
  select(-alleleSeq) 



#V population frequencies  ----

populationFreqs<-alleles_db %>%
  filter(totalClonesCountForGene>50,
         ethnicity!="Unknown",
         geneName %in% frequentVs) %>% 
  bind_rows(
    alleles_db %>%
      filter(totalClonesCountForGene>50,
             geneName %in% frequentVs) %>% 
      mutate(ethnicity="All")
  ) %>%
  group_by(donor,geneName) %>% 
  mutate(zigosity=n(),
         homozigous=ifelse(zigosity==1,"1,1",as.character(zigosity))) %>%
  separate_rows(homozigous,sep=",") %>% 
  mutate(Ethnicity=factor(str_replace(ethnicity,"Hispanic/Latino","Hisp./Lat.") |> 
                            str_replace("Caucasian","Eur.") ,
                          levels=c("All","Eur.","African","Asian", "Hisp./Lat.")))%>% 
  group_by(geneName,Ethnicity) %>% 
  mutate(totalNHaplotypesGene=n()) %>% 
  ungroup() %>% 
  group_by(geneName,alleleName,Ethnicity,totalNHaplotypesGene) %>% 
  summarize(nHaplotypes=n(),
            alleleFreq=nHaplotypes/dplyr::first(totalNHaplotypesGene)) %>%
  ungroup() %>%
  left_join(alleleSeq,by="alleleName")  %>% 
  mutate(alleleNumber=str_remove(alleleName,".*\\*"),
         isNovel=alleleName %in% novelAlleles)

aaGroups<-populationFreqs %>% 
  filter(Ethnicity=="All") %>% 
  group_by(geneName,alleleSeqAA) %>% 
  summarise(alleleAAFreq=sum(alleleFreq),
            nSynonymousAlleles=n()) %>% 
  ungroup() %>% 
  group_by(geneName) %>% 
  mutate(alleleAASeqGroup=data.table::frank(dplyr::desc(alleleAAFreq),ties.method = "random")) %>% 
  ungroup()

populationFreqs %<>% 
  left_join(aaGroups, by=c("geneName","alleleSeqAA"))

##V heatmaps ----

color_fun = colorRamp2(c(0,0.25,0.75,1),c("#FAFAB4","#8DDA7F","#337E90","#1C0F5C"))
g_leg = Legend(col_fun = color_fun, title = "Population frequency")

allVGenes <- populationFreqs %>% 
  filter(Ethnicity=="All") %>% 
  count(geneName) %>% 
  arrange(desc(n)) %>% 
  pull(geneName) 


allHeatmaps<-
  lapply(allVGenes, function(v){
    heatmapMatrixFreqs<-populationFreqs %>% 
      filter(geneName==v) %>% 
      dplyr::group_by(alleleAASeqGroup) %>%
      dplyr::mutate(freqInAll=alleleFreq[Ethnicity=="All"] %>% 
                      sum()) %>%
      ungroup() %>% 
      mutate(alleleAASeqGroup=fct_reorder(factor(alleleAASeqGroup), alleleFreq,.desc = T)) %>%
      dplyr::arrange(desc(freqInAll)) %>%
      pivot_wider(id_cols = Ethnicity, 
                  names_from = c(alleleNumber,alleleAASeqGroup),
                  names_sep = "_",
                  values_from = alleleFreq,
                  id_expand = T)
    
    heatmapMatrixCounts<-populationFreqs %>% 
      filter(geneName==v) %>% 
      dplyr::group_by(alleleAASeqGroup) %>%
      dplyr::mutate(freqInAll=alleleFreq[Ethnicity=="All"] %>% 
                      sum()) %>%
      ungroup() %>% 
      mutate(alleleAASeqGroup=fct_reorder(factor(alleleAASeqGroup), alleleFreq,.desc = T)) %>%
      dplyr::arrange(desc(freqInAll)) %>%
      pivot_wider(id_cols = Ethnicity,
                  names_from = c(alleleNumber,alleleAASeqGroup),
                  names_sep = "_",
                  values_from = nHaplotypes,
                  id_expand = T)   
    
    alleleFeatures<-populationFreqs %>% 
      filter(geneName==v) %>% 
      group_by(alleleName,alleleNumber,isNovel) %>% 
      summarize(alleleFreq=alleleFreq[Ethnicity=="All"]) %>% 
      mutate(novelCol=ifelse(isNovel,mi_bright[4],"black"),
             novelFont=ifelse(isNovel,"bold","plain")) 
    
    alleleFeatures<-alleleFeatures[match(colnames(heatmapMatrixFreqs[,-1]) %>% str_remove("_.*"), 
                                         alleleFeatures$alleleNumber),]
    
    g<-Heatmap(heatmapMatrixFreqs[,-1] ,
               show_row_dend = F,
               show_column_dend = F,
               col = color_fun,
               na_col = "white",
               rect_gp = gpar(col = "darkgrey", lwd = 1),
               row_split = factor(c("A",rep("R",4)),levels=c("A","R")),
               row_gap = unit(0, "mm"),
               row_title = NULL,
               row_labels = heatmapMatrixFreqs[,1]$Ethnicity,
               row_names_side = "left",
               cluster_rows = F,
               column_split = factor(colnames(heatmapMatrixFreqs[,-1]) %>% str_remove(".*_"),
                                     levels=colnames(heatmapMatrixFreqs[,-1]) %>% str_remove(".*_") %>% unique()),
               column_gap = unit(0, "mm"),
               column_labels =  colnames(heatmapMatrixFreqs[,-1]) %>% str_remove("_.*"),
               column_title=v,
               column_title_gp = gpar(fill = "white", col = "black", border = "white",fontface="bold"),
               column_names_gp = gpar(col=alleleFeatures$novelCol,fontface=alleleFeatures$novelFont,fontsize=8),
               column_names_rot = 0,
               column_names_centered = T,
               cluster_column_slices =FALSE,
               border = TRUE,
               border_gp = gpar(lwd=2),
               
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(heatmapMatrixCounts[,-1][i, j] )) {
                   if ( heatmapMatrixFreqs[,-1][i, j] > 0.5) { font_color = "white" } else { font_color = "black"} 
                   grid.text(heatmapMatrixCounts[,-1][i, j], x, y, gp = gpar(fontsize = 10,col=font_color))
                 }
               },
               
               show_heatmap_legend = FALSE
    )
    return(grid.grabExpr(draw(g)))
  })

g<-tibble(x=letters[1:3],y=letters[1:3],val=c(0,0.5,1)) %>% 
  ggplot(aes(x=x ,y=y,fill=val)) +
  geom_tile(color="grey38",linewidth=0.15)+
  scale_fill_gradientn(colors=c("#FAFAB4","#8DDA7F","#337E90","#1C0F5C"),breaks = c(0,0.5,1) ,
                       na.value = "grey99")+labs(fill="Allele\nfrequency")+
  theme(legend.position = "bottom",plot.background = element_rect(fill = "white", colour = NA),
        panel.background=element_rect(fill = "white", colour = NA))

fig4_leg<-cowplot::get_legend(g )



figS5a<-cowplot::plot_grid(plotlist = allHeatmaps[1:6],nrow=6)
figS5b<-cowplot::plot_grid(plotlist = allHeatmaps[7:16],ncol=2)
figS5c<-cowplot::plot_grid(plotlist = allHeatmaps[17:42],ncol=3)

#IGHV1-69, IGHV3-23,IGHV3-48,IGH3-7

## Fig 4 ----
fig4a<-cowplot::plot_grid(plotlist = allHeatmaps[c(1,4)],nrow=2)
fig4b<-cowplot::plot_grid(plotlist = allHeatmaps[c(18,23)],ncol=2)

fig4<-
  plot_grid(fig4a,fig4b,fig4_leg,nrow=3,
            rel_heights = c(1,0.5,0.1))


ggsave(filename = "plots/Fig4.pdf", plot= fig4,device = "pdf",
       width = 10,height = 8,units = "in",bg="white")
ggsave(filename = "plots/Fig4.png", plot= fig4,device = "png",
       width = 10,height = 8,units = "in",bg="white")

## Fig S5 ----

figS5<-
  plot_grid(
    plot_grid(figS5a,figS5b,nrow=2),
    figS4c,
    rel_widths = c(0.7,1),
    nrow=1)

figS5<-ggdraw()+
  draw_plot(figS5,x=0,y=0,width=1,height=1)+
  draw_plot(fig4_leg,x=0.75,y=0.017,width=0.3,height = 0.1)

ggsave(filename = "plots/FigS5.pdf", plot= figS5,device = "pdf",
       width = 30,height = 24,units = "in",bg="white")
ggsave(filename = "plots/FigS5.png", plot= figS5,device = "png",
       width = 30,height = 24,units = "in",bg="white")





#V gene distribution divergence by  allele freq distributions ----

allVGenes <- populationFreqs %>% 
  filter(Ethnicity=="All") %>% 
  count(geneName) %>% 
  arrange(desc(n)) %>% 
  mutate(geneName=fct_reorder(geneName,n)) %>% 
  pull(geneName)



js_df<-lapply(allVGenes, function(v){
  heatmapMatrixFreqs<-populationFreqs %>% 
    filter(geneName==v) %>% 
    dplyr::group_by(alleleAASeqGroup) %>%
    dplyr::mutate(freqInAll=alleleFreq[Ethnicity=="All"] %>% 
                    sum()) %>%
    ungroup() %>% 
    mutate(nHaplotypes=as.numeric(nHaplotypes)) %>% 
    complete(nesting(alleleName,geneName,alleleNumber),Ethnicity, fill = list(nHaplotypes=0)) %>% 
    group_by(Ethnicity) %>% 
    mutate(alleleFreqNorm=(nHaplotypes+1)/sum(nHaplotypes+1)) |> 
    ungroup()
  
  
  ethnPairs<-combn(heatmapMatrixFreqs$Ethnicity %>% unique() %>% as.character,2,simplify = F)
  
  js_df<-lapply(ethnPairs,function(ethnPair){
    df<-inner_join(heatmapMatrixFreqs %>% filter(Ethnicity==ethnPair[[1]]) %>% select(alleleNumber,alleleFreqNorm),
                   heatmapMatrixFreqs %>% filter(Ethnicity==ethnPair[[2]]) %>% select(alleleNumber,alleleFreqNorm),
                   by="alleleNumber")
    m <- 0.5 * (df$alleleFreqNorm.x + df$alleleFreqNorm.y)
    JS <- 0.5 * (sum(df$alleleFreqNorm.x * log2(df$alleleFreqNorm.x / m)) +
                   sum(df$alleleFreqNorm.y * log2(df$alleleFreqNorm.y / m)))
    
    
    bhattacharyya<- -log(sum(sqrt(df$alleleFreqNorm.x*df$alleleFreqNorm.y)))
    hellinger <- sqrt( sum( ( sqrt(df$alleleFreqNorm.x)-sqrt(df$alleleFreqNorm.y) )^2) ) / sqrt(2)
    
    divergence<-tibble(Ethnicity.x = ethnPair[[1]],
           Ethnicity.y = ethnPair[[2]],
           vGene = v,
           JensenShannonDivergence=JS,
           bhattacharyya=bhattacharyya,
           hellinger=hellinger)
    return(divergence)
  }) %>% bind_rows
  
  # js_df<-plyr::ldply(ethnPairs, rbind) %>% 
  #   mutate(js_divergence=js,
  #          vGene=v) %>% 
  #   set_colnames(c("Ethnicity.x","Ethnicity.y","JensenShannonDivergence","vGene"))
  
  
  return(js_df)
}) %>% bind_rows()


## Hellinger / Jensen-Shannon Permutation test ----
job::job({
vjDistPermTest <- lapply(1:1000, function(it) {
  
  lapply(allVGenes, function(v){
    ethnicities_reshuffled <- donors_ethnicities |> mutate(ethnicity=sample(ethnicity))
    
    populationFreqs<-alleles_db %>%
      select(-ethnicity) |> 
      left_join(ethnicities_reshuffled,by="donor") |> 
      filter(totalClonesCountForGene>50,
             ethnicity!="Unknown",
             ethnicity!="Hispanic/Latino",
             geneName %in% frequentVs) %>% 
      group_by(donor,geneName) %>% 
      mutate(zigosity=n(),
             homozigous=ifelse(zigosity==1,"1,1",as.character(zigosity))) %>%
      separate_rows(homozigous,sep=",") %>% 
      mutate(Ethnicity=factor(ethnicity,
                              levels=c("Caucasian","African","Asian")))%>% 
      group_by(geneName,Ethnicity) %>% 
      mutate(totalNHaplotypesGene=n()) %>% 
      ungroup() %>% 
      group_by(geneName,alleleName,Ethnicity,totalNHaplotypesGene) %>% 
      summarize(nHaplotypes=n(),
                alleleFreq=nHaplotypes/dplyr::first(totalNHaplotypesGene)) %>%
      ungroup() %>% 
      mutate(alleleNumber=str_remove(alleleName,".*\\*"))
    
    heatmapMatrixFreqs<-populationFreqs %>% 
      filter(geneName==v) %>% 
      mutate(nHaplotypes=as.numeric(nHaplotypes)) %>% 
      complete(nesting(alleleName,geneName,alleleNumber),Ethnicity, fill = list(nHaplotypes=0)) %>% 
      group_by(Ethnicity) %>% 
      mutate(alleleFreqNorm=(nHaplotypes+1)/sum(nHaplotypes+1)) |> 
      ungroup()
    
    js_df<-lapply(ethnPairs,function(ethnPair){
      df<-inner_join(heatmapMatrixFreqs %>% filter(Ethnicity==ethnPair[[1]]) %>% select(alleleNumber,alleleFreqNorm),
                     heatmapMatrixFreqs %>% filter(Ethnicity==ethnPair[[2]]) %>% select(alleleNumber,alleleFreqNorm),
                     by="alleleNumber")
      m <- 0.5 * (df$alleleFreqNorm.x + df$alleleFreqNorm.y)
      JS <- 0.5 * (sum(df$alleleFreqNorm.x * log2(df$alleleFreqNorm.x / m)) +
                     sum(df$alleleFreqNorm.y * log2(df$alleleFreqNorm.y / m)))
        #  bhattacharyya<- -log(sum(sqrt(df$alleleFreqNorm.x*df$alleleFreqNorm.y)))
      hellinger <- sqrt( sum( ( sqrt(df$alleleFreqNorm.x)-sqrt(df$alleleFreqNorm.y) )^2) ) / sqrt(2)
      
      divergence<-tibble(Ethnicity.x = ethnPair[[1]],
                         Ethnicity.y = ethnPair[[2]],
                         vGene = v,
                         JensenShannonDivergence=JS,
                     #    bhattacharyya=bhattacharyya,
                         hellinger=hellinger)
      return(divergence)
    }) %>% bind_rows() |> 
      mutate(iteration=it)
    return(js_df)
  }) %>% bind_rows() 
    
}) %>% bind_rows()
},import = "auto")

vjDistPermTest |> write_tsv("VJ_permutation.tsv")


vjDistPermTest.pvalue<-vjDistPermTest |> 
  left_join(js_df,
            by = c("Ethnicity.x","Ethnicity.y","vGene"),
            suffix = c('.reshuf','.orig')) |> 
  mutate(hlngrGreaterThanInOrig=hellinger.reshuf >= hellinger.orig ) |> 
  group_by(Ethnicity.x, Ethnicity.y,vGene) |> 
  summarise(p=sum(hlngrGreaterThanInOrig)/max(iteration)) |> 
  ungroup() |> 
  mutate(p.adj=p.adjust(p,"BH"))

contrast <- function(colour) {
  out   <- rep("black", length(colour))
  light <- farver::get_channel(colour, "l", space = "hcl")
  out[light < 50] <- "white"
  out
}  

autocontrast <- aes(colour = after_scale(contrast(fill)))

##Fig S6 Hellinger V distance ----
figS6<-js_df %>% 
  filter(Ethnicity.x!="All",!str_detect(Ethnicity.x,"Hisp"),!str_detect(Ethnicity.y,"Hisp")) %>% 
  bind_rows(js_df %>% 
              filter(Ethnicity.x!="All",!str_detect(Ethnicity.x,"Hisp"),!str_detect(Ethnicity.y,"Hisp")) %>% 
              mutate(Ethnicity=Ethnicity.y,
                     Ethnicity.y=Ethnicity.x,
                     Ethnicity.x=Ethnicity)) %>% 
  left_join(populationFreqs %>% 
              filter(Ethnicity=="All") %>% 
              count(geneName),
            by=join_by(vGene==geneName)) %>% 
  left_join(vjDistPermTest.pvalue |> 
            bind_rows(vjDistPermTest.pvalue |> 
                        mutate(ethn.tmp=Ethnicity.x,
                                 Ethnicity.x=Ethnicity.y,
                                 Ethnicity.y=ethn.tmp) |> 
                        select(-ethn.tmp)
            ), by=c("Ethnicity.x","Ethnicity.y","vGene")) |> 
  mutate(vGene=fct_reorder(vGene,n,.desc = T),
         # Ethnicity.x=str_replace(Ethnicity.x,"Hispanic/Latino","Hisp/Lat"),
         # Ethnicity.y=str_replace(Ethnicity.y,"Hispanic/Latino","Hisp/Lat"),
         Ethnicity.x=str_replace(Ethnicity.x,"Caucasian","Eur."),
         Ethnicity.y=str_replace(Ethnicity.y,"Caucasian","Eur."),
         p.sign=case_when(p.adj > 0.05 ~ "",
                          p.adj>= 0.01 ~ "*",
                          TRUE ~ "**")) |> 
  ggplot(aes(x=Ethnicity.x ,y=Ethnicity.y,fill=hellinger)) +
  facet_wrap( ~ vGene,scales="free_x")+
  theme_classic()+
  geom_tile(color="grey38",linewidth=0.15)+
  geom_text(aes(label=p.sign, !!!autocontrast), size=5, vjust=0.5, hjust=0.5) +
  scale_fill_gradientn(colors=c("#FAFAB4","#8DDA7F","#337E90","#1C0F5C"),na.value = "grey99")+
  labs(fill="Hellinger\ndistance") +
  theme(axis.title = element_blank())
#theme(legend.position = "none")

ggsave(filename = "plots/figS6.pdf", plot= figS6,device = "pdf",
       width = 16,height = 12,units = "in",bg="white")
ggsave(filename = "plots/figS6.png", plot= figS6,device = "png",
       width = 16,height = 12,units = "in",bg="white")



#J population frequencies  ----

populationFreqsJ<-alleles_db %>%
  filter(totalClonesCountForGene>50,
         ethnicity!="Unknown",
         str_detect(geneName,"IGHJ")) %>% 
  bind_rows(
    alleles_db %>%
      filter(totalClonesCountForGene>50,
             str_detect(geneName,"IGHJ")) %>% 
      mutate(ethnicity="All")
  ) %>%
  group_by(donor,geneName) %>% 
  mutate(zigosity=n(),
         homozigous=ifelse(zigosity==1,"1,1",as.character(zigosity))) %>%
  separate_rows(homozigous,sep=",") %>% 
  mutate(Ethnicity=factor(str_replace(ethnicity,"Hispanic/Latino","Hisp./Lat.") |> 
                            str_replace("Caucasian","European") ,
                          levels=c("All","European","African","Asian", "Hisp./Lat.")))%>% 
  group_by(geneName,Ethnicity) %>% 
  mutate(totalNHaplotypesGene=n()) %>% 
  ungroup() %>% 
  group_by(geneName,alleleName,Ethnicity,totalNHaplotypesGene) %>% 
  summarize(nHaplotypes=n(),
            alleleFreq=nHaplotypes/dplyr::first(totalNHaplotypesGene)) %>%
  ungroup() %>% 
  mutate(alleleNumber=str_remove(alleleName,".*\\*"),
         isNovel=alleleName %in% novelAlleles)



## J Heatmaps ----

color_fun = colorRamp2(c(0,0.25,0.75,1),c("#FAFAB4","#8DDA7F","#337E90","#1C0F5C"))
g_leg = Legend(col_fun = color_fun, title = "Population frequency")

allJGenes <- populationFreqsJ %>% 
  filter(Ethnicity=="All") %>% 
  count(geneName) %>% 
  arrange(desc(n)) %>% 
  pull(geneName) 


allHeatmapsJ<-
  lapply(allJGenes, function(j){
    heatmapMatrixFreqs<-populationFreqsJ %>% 
      filter(geneName==j) %>% 
      mutate(alleleAASeqGroup=alleleNumber) %>% 
      dplyr::group_by(alleleAASeqGroup) %>%
      dplyr::mutate(freqInAll=cur_data() %>% 
                      filter(Ethnicity=="All")%>%
                      pull(alleleFreq) %>% 
                      mean()) %>%
      ungroup() %>% 
      mutate(alleleAASeqGroup=fct_reorder(factor(alleleAASeqGroup), alleleFreq,.desc = T)) %>%
      group_by(alleleAASeqGroup) %>% 
      dplyr::arrange(freqInAll,.by_group = TRUE) %>%
      pivot_wider(id_cols = Ethnicity, 
                  names_from = c(alleleNumber,alleleAASeqGroup),
                  names_sep = "_",
                  values_from = alleleFreq,
                  id_expand = T)
    
    heatmapMatrixCounts<-populationFreqsJ %>% 
      filter(geneName==j) %>% 
      mutate(alleleAASeqGroup=alleleNumber) %>% 
      dplyr::group_by(alleleAASeqGroup) %>% dplyr::mutate(freqInAll=cur_data() %>% filter(Ethnicity=="All")%>% pull(alleleFreq) %>% mean()) %>%
      ungroup() %>% 
      mutate(alleleAASeqGroup=fct_reorder(factor(alleleAASeqGroup), alleleFreq,.desc = T)) %>%
      group_by(alleleAASeqGroup) %>% 
      dplyr::arrange(freqInAll,.by_group = TRUE) %>%
      pivot_wider(id_cols = Ethnicity,
                  names_from = c(alleleNumber,alleleAASeqGroup),
                  names_sep = "_",
                  values_from = nHaplotypes,
                  id_expand = T)   
    
    alleleFeatures<-populationFreqsJ %>% 
      filter(geneName==j) %>% 
      group_by(alleleName,isNovel) %>% 
      summarize(alleleFreq=alleleFreq[Ethnicity=="All"]) %>% 
      mutate(novelCol=ifelse(isNovel,mi_bright[4],"black"),
             novelFont=ifelse(isNovel,"bold","plain")) %>% 
      arrange(desc(alleleFreq))
    
    
    g<-Heatmap(heatmapMatrixFreqs[,-1] ,
               column_split = factor(colnames(heatmapMatrixFreqs[,-1]) %>% str_remove(".*_")) %>% 
                 fct_reorder(heatmapMatrixFreqs[1,-1] %>% as.matrix() %>% t() %>% as.vector(),.desc = T),
               show_row_dend = F,
               show_column_dend = F,
               col = color_fun,
               na_col = "white",
               rect_gp = gpar(col = "darkgrey", lwd = 1),
               column_labels =  colnames(heatmapMatrixFreqs[,-1]) %>% str_remove("_.*"),
               row_labels = heatmapMatrixFreqs[,1]$Ethnicity,
               row_names_side = "left",
               column_title=j,
               column_title_gp = gpar(fill = "white", col = "black", border = "white",fontface="bold"),
               column_names_gp = gpar(col=alleleFeatures$novelCol,fontface=alleleFeatures$novelFont),
               
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(heatmapMatrixCounts[,-1][i, j] )) {
                   if ( heatmapMatrixFreqs[,-1][i, j] > 0.5) { font_color = "white" } else { font_color = "black"} 
                   grid.text(heatmapMatrixCounts[,-1][i, j], x, y, gp = gpar(fontsize = 10,col=font_color))
                 }
               },
               column_gap = unit(0, "mm"),
               border = TRUE,
               border_gp = gpar(lwd=2),
               show_heatmap_legend = FALSE,
               cluster_column_slices =FALSE,
               cluster_rows = F,
               column_names_rot = 0,
               column_names_centered = T
    )
    return(grid.grabExpr(draw(g)))
  })

figS7<-plot_grid(
  plot_grid(plotlist = allHeatmapsJ[c(  4,5,2,4,6)],nrow=1),
  plot_grid(plotlist = allHeatmapsJ[1],nrow=1),
  fig4_leg,
  rel_heights = c(1,1,0.1),
  nrow=3)
ggsave(filename = "plots/figS7.pdf", plot= figS7,device = "pdf",
       width = 12,height = 8,units = "in",bg="white")
ggsave(filename = "plots/figS7.png", plot= figS7,device = "png",
       width = 12,height = 8,units = "in",bg="white")




# TCR Population frequencies ----

alleles_tcr<-read_tsv("/home/amikelov/paper_alleles/mergedLib/tcr_tmp_meta_filt.tsv") %>% 
  filter(str_detect(alleleName,"TRA|TRB")) %>% 
  mutate(geneName=str_remove(alleleName,"\\*.*"),
         ethnicity=str_replace(ethnicity,"hispanic","Hispanic/Latino") %>% str_to_title(),
         chain=ifelse(str_detect(alleleName,"TRA"),"TRA","TRB"),
         geneType=ifelse(str_detect(alleleName,"V"),"V","J") %>% factor(levels=c("V","J"))) 
#filter(str_detect(geneName,'IGHV|IGHJ')) %>% 
# left_join(findAllelesOutput %>% 
#             select(donor,alleleName,clonesCount,totalClonesCountForGene),
#           join_by(donor==donor, alleleNameInSource==alleleName)) %>% 


cohortDemographics_tcr<-alleles_tcr %>% 
  group_by(ethnicity) %>% 
  summarize(nDonors=n_distinct(donor)) %>%
  mutate(nDonorsPercent=percent(nDonors/sum(nDonors),accuracy= 0.1) )

cohortDemographics_tcr %>% 
  ggplot(aes(area = nDonors, fill = ethnicity,label = ethnicity %+% "\n" %+% nDonors %+% "\n" %+%nDonorsPercent)) +
  geom_treemap(alpha=0.7)+
  geom_treemap_text()+
  theme(legend.position = "none")+
  scale_fill_manual(values = ethnicity_colors)


allelesStats_tcr<-alleles_tcr %>%
  mutate(nDonors=n_distinct(donor),
         isNovel=ifelse(str_detect(alleleName,"x"),"Novel","Known") %>% factor(levels=c("Novel","Known"))) %>% 
  group_by(chain,geneName,geneType,isNovel) %>% 
  summarise(nAlleles=n_distinct(alleleName),
            nAllelesPerDonors=nAlleles/first(nDonors)) 



findAllelesOutput_tcr <- dir("/home/amikelov/paper_alleles/mergedLib/liege_tcr",
                             pattern=".tsv",
                             full.names = T)

findAllelesOutput_tcr<-read_tsv(findAllelesOutput_tcr,id="fileName") %>% 
  mutate(donor=str_remove_all(fileName,".*\\/|_alleles.tsv"),.after = fileName) %>% 
  filter(str_detect(geneName,'TRBV|TRBJ|TRAV|TRAJ')) %>% 
  select(donor:totalClonesCountForGene)


geneStats_tcr<-findAllelesOutput_tcr %>%
  mutate(chain=ifelse(str_detect(alleleName,"TRA"),"TRA","TRB"),
         geneType=ifelse(str_detect(geneName,"V"),"V","J") %>% factor(levels=c("V","J"))) %>% 
  group_by(donor,chain,geneType) %>% 
  dplyr::mutate(totalClonesInDonor=sum(clonesCount),
                geneFrequency=totalClonesCountForGene/totalClonesInDonor) %>% 
  ungroup() %>% 
  # filter(totalClonesCountForGene>10) %>% 
  group_by(geneName,chain,geneType) %>% 
  summarise(nAlleles=n_distinct(alleleName),
            nDonors=n_distinct(donor),
            medianFrequency=median(geneFrequency)) %>%
  mutate(rank=dense_rank(dplyr::desc(medianFrequency)),
         donorFraction=nDonors/length(goodDepthDonors)) %>% 
  ungroup()

## Fig 2A  TRAV N alleles ----
fig2a<-allelesStats_tcr %>% 
  filter(chain=="TRA",geneType=="V") %>% 
  left_join(geneStats_tcr %>% select(-chain,-geneType,-nAlleles),by=c("geneName")) %>% 
  ggplot(aes(x=fct_reorder(geneName,medianFrequency),y=nAlleles,fill=isNovel))+
  geom_bar(stat="identity")+
  facet_wrap(~chain,
             scales="free",
             labeller = labeller(chain=c(TRA="TRA, V genes",TRB="TRB, V genes")))+
  theme_bw()+
  coord_flip()+
  scale_fill_manual(values=c(Novel="#2D93FA",Known="#99CCFF"))+
  theme(strip.background =element_rect(fill="white"),
        legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.title.y = element_blank())+
  labs(x= "Genes", y="# Detected alleles")+
  geom_text(aes(color = isNovel,label=nAlleles),
            size = 3, 
            position = position_stack(vjust = 0.5))+
  scale_color_manual(values=c("white","black"))+
  scale_y_continuous(expand = expansion(mult = c(0.005, .1))) +
  annotate('label',x=-Inf,y=Inf,hjust=1,vjust=0,
           label= length(unique(alleles_tcr$donor)) %+% " ind." )

## Fig 2B  TRBV N alleles ----
fig2b<-allelesStats_tcr %>% 
  filter(chain=="TRB",geneType=="V") %>% 
  left_join(geneStats_tcr %>% select(-chain,-geneType,-nAlleles),by=c("geneName")) %>% 
  filter(!is.na(medianFrequency)) %>% 
  ggplot(aes(x=fct_reorder(geneName,medianFrequency),y=nAlleles,fill=isNovel))+
  geom_bar(stat="identity")+
  facet_wrap(~chain,
             scales="free",
             labeller = labeller(chain=c(TRA="TRA, V genes",TRB="TRB, V genes")))+
  theme_bw()+
  coord_flip()+
  scale_fill_manual(values=c(Novel="#2D93FA",Known="#99CCFF"))+
  theme(strip.background =element_rect(fill="white"),
        legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.title.y = element_blank())+
  labs(x= "Genes", y="# Detected alleles")+
  geom_text(aes(color = isNovel,label=nAlleles),
            size = 3, 
            position = position_stack(vjust = 0.5))+
  scale_color_manual(values=c("white","black"))+
  scale_y_continuous(expand = expansion(mult = c(0.005, .1))) +
  annotate('label',x=-Inf,y=Inf,hjust=1,vjust=0,
           label= length(unique(alleles_tcr$donor)) %+% " ind." )


## Fig 2D  TRAJ N alleles ----
fig2d<-allelesStats_tcr %>% 
  filter(chain=="TRA",geneType=="J") %>% 
  left_join(geneStats_tcr %>% select(-chain,-geneType,-nAlleles),by=c("geneName")) %>% 
  ggplot(aes(x=fct_reorder(geneName,medianFrequency),y=nAlleles,fill=isNovel))+
  geom_bar(stat="identity")+
  facet_wrap(~chain,
             scales="free",
             labeller = labeller(chain=c(TRA="TRA, J genes",TRB="TRB, J genes")))+
  theme_bw()+
  coord_flip()+
  scale_fill_manual(values=c(Novel="#2D93FA",Known="#99CCFF"))+
  theme(strip.background =element_rect(fill="white"),
        legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank())+
  labs(x= "Genes", y="# Detected alleles")+
  geom_text(aes(color = isNovel,label=nAlleles),
            size = 3, 
            position = position_stack(vjust = 0.5))+
  scale_color_manual(values=c("white","black"))+
  scale_y_continuous(expand = expansion(mult = c(0.005, .05)),
                     breaks = c(0,1,2)) +
  annotate('label',x=-Inf,y=Inf,hjust=1,vjust=0,
           label= length(unique(alleles_tcr$donor)) %+% " ind." )

## Fig 2E  TRBJ N alleles ----
fig2e<-allelesStats_tcr %>% 
  filter(chain=="TRB",geneType=="J") %>% 
  left_join(geneStats_tcr %>% select(-chain,-geneType,-nAlleles,-rank),by=c("geneName")) %>% 
  filter(!is.na(medianFrequency)) %>% 
  ggplot(aes(x=fct_reorder(geneName,medianFrequency),y=nAlleles,fill=isNovel))+
  geom_bar(stat="identity")+
  facet_wrap(~chain,
             scales="free",
             labeller = labeller(chain=c(TRA="TRA, J genes",TRB="TRB, J genes")))+
  theme_bw()+
  coord_flip()+
  scale_fill_manual(values=c(Novel="#2D93FA",Known="#99CCFF"))+
  theme(strip.background =element_rect(fill="white"),
        legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  labs(x= "Genes", y="# Detected alleles")+
  geom_text(aes(color = isNovel,label=nAlleles),
            size = 3, 
            position = position_stack(vjust = 0.5))+
  scale_color_manual(values=c("white","black"))+
  scale_y_continuous(expand = expansion(mult = c(0.005, .05)),
                     breaks = c(0,1,2)) +
  annotate('label',x=-Inf,y=Inf,hjust=1,vjust=0,
           label= length(unique(alleles_tcr$donor)) %+% " ind." )



populationFreqs_tcr<-alleles_tcr%>%
  group_by(donor,geneName) %>% 
  mutate(zigosity=n(),
         homozigous=ifelse(zigosity==1,"1,1",as.character(zigosity))) %>%
  separate_rows(homozigous,sep=",") %>% 
  # mutate(Ethnicity=factor(str_replace(ethnicity,"Hispanic/Latino","Hisp./Lat."),
  #                         levels=c("All","Caucasian","African","Asian", "Hisp./Lat.")))%>% 
  group_by(geneName) %>% 
  mutate(totalNHaplotypesGene=n()) %>% 
  ungroup() %>% 
  group_by(geneName,alleleName,totalNHaplotypesGene) %>% 
  summarize(nHaplotypes=n(),
            alleleFreq=nHaplotypes/dplyr::first(totalNHaplotypesGene)) %>%
  ungroup() %>%
  #left_join(alleleSeq,by="alleleName")  %>% 
  mutate(alleleNumber=str_remove(alleleName,".*\\*"),
         isNovel=str_detect(alleleName,"x"))

# aaGroups<-populationFreqs %>% 
#   filter(Ethnicity=="All") %>% 
#   group_by(geneName,alleleSeqAA) %>% 
#   summarise(alleleAAFreq=sum(alleleFreq),
#             nSynonymousAlleles=n()) %>% 
#   ungroup() %>% 
#   group_by(geneName) %>% 
#   mutate(alleleAASeqGroup=data.table::frank(dplyr::desc(alleleAAFreq),ties.method = "random")) %>% 
#   ungroup()

# populationFreqs %<>% 
#   left_join(aaGroups, by=c("geneName","alleleSeqAA"))

##TRBV heatmaps ----

color_fun = colorRamp2(c(0,0.25,0.75,1),c("#FAFAB4","#8DDA7F","#337E90","#1C0F5C"))
g_leg = Legend(col_fun = color_fun, title = "Population frequency")

gen<-"TRBV"

allGenes <- populationFreqs_tcr %>% 
  filter(str_detect(alleleName,gen)) %>% 
  mutate(vFamily=str_remove(geneName,gen) %>% str_remove("-.*") %>% as.integer()) %>% 
  count(geneName,vFamily) %>% 
  arrange(desc(n)) %>% 
  pull(geneName) 



allHeatmapsTRBV<-
  lapply(allGenes, function(v){
    heatmapMatrixFreqs<-populationFreqs_tcr %>% 
      filter(geneName==v) %>% 
      # mutate(alleleAASeqGroup=alleleNumber) %>% 
      # ungroup() %>% 
      # mutate(alleleAASeqGroup=fct_reorder(factor(alleleAASeqGroup), alleleFreq,.desc = T)) %>%
      # group_by(alleleAASeqGroup) %>% 
      dplyr::arrange(desc(alleleFreq)) %>%
      pivot_wider(id_cols = geneName, 
                  names_from = alleleNumber,
                  values_from = alleleFreq,
                  id_expand = T)
    
    heatmapMatrixCounts<-populationFreqs_tcr %>% 
      filter(geneName==v) %>% 
      # mutate(alleleAASeqGroup=alleleNumber) %>% 
      # ungroup() %>% 
      # mutate(alleleAASeqGroup=fct_reorder(factor(alleleAASeqGroup), alleleFreq,.desc = T)) %>%
      # group_by(alleleAASeqGroup) %>% 
      dplyr::arrange(desc(alleleFreq)) %>%
      pivot_wider(id_cols = geneName, 
                  names_from = alleleNumber,
                  values_from = nHaplotypes,
                  id_expand = T)
    
    alleleFeatures<-populationFreqs_tcr %>% 
      filter(geneName==v) %>% 
      mutate(novelCol=ifelse(isNovel,mi_bright[4],"black"),
             novelFont=ifelse(isNovel,"bold","plain")) %>% 
      arrange(desc(alleleFreq))
    
    
    g<-Heatmap(heatmapMatrixFreqs[,-1] ,
               # column_split = factor(colnames(heatmapMatrixFreqs[,-1]) %>% str_remove(".*_")) %>% 
               #   fct_reorder(heatmapMatrixFreqs[1,-1] %>% as.matrix() %>% t() %>% as.vector(),.desc = T),
               show_row_dend = F,
               show_column_dend = F,
               col = color_fun,
               na_col = "white",
               rect_gp = gpar(col = "darkgrey", lwd = 1),
               column_labels =  colnames(heatmapMatrixFreqs[,-1]) %>% str_remove("_.*"),
               row_labels = heatmapMatrixFreqs[,1]$geneName %>% str_pad(side="right",width = 8),
               row_names_side = "left",
               row_names_gp = gpar(fontfamily = "mono"),
               # column_title=v,
               # column_title_gp = gpar(fill = "white", col = "black", border = "white",fontface="bold",fontsize=2),
               column_names_gp = gpar(col=alleleFeatures$novelCol,fontface=alleleFeatures$novelFont),
               
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(heatmapMatrixCounts[,-1][i, j] )) {
                   if ( heatmapMatrixFreqs[,-1][i, j] > 0.5) { font_color = "white" } else { font_color = "black"} 
                   grid.text(heatmapMatrixCounts[,-1][i, j], x, y, gp = gpar(fontsize = 10,col=font_color))
                 }
               },
               column_gap = unit(0, "mm"),
               border = TRUE,
               border_gp = gpar(lwd=2),
               show_heatmap_legend = FALSE,
               cluster_column_slices =FALSE,
               cluster_rows = F,
               column_names_rot = 0,
               column_names_centered = T
    )
    return(grid.grabExpr(draw(g)))
  })

figS6<-plot_grid(
  cowplot::plot_grid(plotlist = allHeatmapsTRBV[1:29],nrow=29),
  cowplot::plot_grid(plotlist = allHeatmapsTRBV[30:58],nrow=29),
  nrow=1
)

ggsave(filename = "FigS6.pdf", plot= figS6,device = "pdf",
       width = 16,height = 20,units = "in",bg="white")
ggsave(filename = "FigS6.png", plot= figS6,device = "png",
       width = 16,height = 20,units = "in",bg="white")

##TRAV heatmaps ----

color_fun = colorRamp2(c(0,0.25,0.75,1),c("#FAFAB4","#8DDA7F","#337E90","#1C0F5C"))
g_leg = Legend(col_fun = color_fun, title = "Population frequency")

gen<-"TRAV"

allGenes <- populationFreqs_tcr %>% 
  filter(str_detect(alleleName,gen)) %>% 
  mutate(vFamily=str_remove(geneName,gen) %>% str_remove("-.*") %>% as.integer()) %>% 
  count(geneName,vFamily) %>% 
  arrange(desc(n)) %>% 
  pull(geneName) 



allHeatmapsTRAV<-
  lapply(allGenes, function(v){
    heatmapMatrixFreqs<-populationFreqs_tcr %>% 
      filter(geneName==v) %>% 
      # mutate(alleleAASeqGroup=alleleNumber) %>% 
      # ungroup() %>% 
      # mutate(alleleAASeqGroup=fct_reorder(factor(alleleAASeqGroup), alleleFreq,.desc = T)) %>%
      # group_by(alleleAASeqGroup) %>% 
      dplyr::arrange(desc(alleleFreq)) %>%
      pivot_wider(id_cols = geneName, 
                  names_from = alleleNumber,
                  values_from = alleleFreq,
                  id_expand = T)
    
    heatmapMatrixCounts<-populationFreqs_tcr %>% 
      filter(geneName==v) %>% 
      # mutate(alleleAASeqGroup=alleleNumber) %>% 
      # ungroup() %>% 
      # mutate(alleleAASeqGroup=fct_reorder(factor(alleleAASeqGroup), alleleFreq,.desc = T)) %>%
      # group_by(alleleAASeqGroup) %>% 
      dplyr::arrange(desc(alleleFreq)) %>%
      pivot_wider(id_cols = geneName, 
                  names_from = alleleNumber,
                  values_from = nHaplotypes,
                  id_expand = T)
    
    alleleFeatures<-populationFreqs_tcr %>% 
      filter(geneName==v) %>% 
      mutate(novelCol=ifelse(isNovel,mi_bright[4],"black"),
             novelFont=ifelse(isNovel,"bold","plain")) %>% 
      arrange(desc(alleleFreq))
    
    
    g<-Heatmap(heatmapMatrixFreqs[,-1] ,
               # column_split = factor(colnames(heatmapMatrixFreqs[,-1]) %>% str_remove(".*_")) %>% 
               #   fct_reorder(heatmapMatrixFreqs[1,-1] %>% as.matrix() %>% t() %>% as.vector(),.desc = T),
               show_row_dend = F,
               show_column_dend = F,
               col = color_fun,
               na_col = "white",
               rect_gp = gpar(col = "darkgrey", lwd = 1),
               column_labels =  colnames(heatmapMatrixFreqs[,-1]) %>% str_remove("_.*"),
               row_labels = heatmapMatrixFreqs[,1]$geneName %>% str_pad(side="right",width = 11),
               row_names_side = "left",
               row_names_gp = gpar(fontfamily = "mono"),
               # column_title=v,
               # column_title_gp = gpar(fill = "white", col = "black", border = "white",fontface="bold",fontsize=2),
               column_names_gp = gpar(col=alleleFeatures$novelCol,fontface=alleleFeatures$novelFont),
               
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(!is.na(heatmapMatrixCounts[,-1][i, j] )) {
                   if ( heatmapMatrixFreqs[,-1][i, j] > 0.5) { font_color = "white" } else { font_color = "black"} 
                   grid.text(heatmapMatrixCounts[,-1][i, j], x, y, gp = gpar(fontsize = 10,col=font_color))
                 }
               },
               column_gap = unit(0, "mm"),
               border = TRUE,
               border_gp = gpar(lwd=2),
               show_heatmap_legend = FALSE,
               cluster_column_slices =FALSE,
               cluster_rows = F,
               column_names_rot = 0,
               column_names_centered = T
    )
    return(grid.grabExpr(draw(g)))
  })

figS7<-plot_grid(
  cowplot::plot_grid(plotlist = allHeatmapsTRAV[1:23],nrow=23),
  cowplot::plot_grid(plotlist = allHeatmapsTRAV[24:45],nrow=22),
  nrow=1
)

ggsave(filename = "FigS7.pdf", plot= figS7,device = "pdf",
       width = 16,height = 20,units = "in",bg="white")
ggsave(filename = "FigS7.png", plot= figS7,device = "png",
       width = 16,height = 20,units = "in",bg="white")
#Assemble figures ----
##Figure2 ----
library(patchwork)

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  fig2a + theme(legend.position="bottom",
                legend.box.margin = margin(0, 0, 0, 12))+
    scale_colour_discrete(guide = "none")+
    labs(fill="Alleles")
)
# 
# fig2a<-fig2a+theme(axis.title.y = element_text(size=14),     
#         axis.text.x = element_text(size=12),     
#         axis.text.y = element_text(size=12),     
#         strip.text = element_text(size=14))
# fig2b<-fig2b+theme(axis.title.y = element_text(size=14),     
#                    axis.text.x = element_text(size=12),     
#                    axis.text.y = element_text(size=12),     
#                    strip.text = element_text(size=14))
# fig2c<-fig2c+theme(axis.title.y = element_text(size=14),     
#                    axis.text.x = element_text(size=12),     
#                    axis.text.y = element_text(size=12),     
#                    strip.text = element_text(size=14))
# fig2d<-fig2d+theme(axis.title.y = element_text(size=14),     
#                    axis.text.x = element_text(size=12),     
#                    axis.text.y = element_text(size=12),     
#                    strip.text = element_text(size=14))
# fig2e<-fig2e+theme(axis.title.y = element_text(size=14),     
#                    axis.text.x = element_text(size=12),     
#                    axis.text.y = element_text(size=12),     
#                    strip.text = element_text(size=14))
# fig2f<-fig2f+theme(axis.title.y = element_text(size=14),     
#                    axis.text.x = element_text(size=12),     
#                    axis.text.y = element_text(size=12),     
#                    strip.text = element_text(size=14))


fig2<-plot_grid(fig2a,fig2b,fig2c,fig2d,
                plot_grid(fig2e,fig2f, legend,
                          nrow=3,
                          rel_heights =   c(1.3,0.75,0.2),
                          labels=c("","F")),
                rel_widths = c(1,1,1.4,1,1),
                labels=c("A","B","C","D","E"),
                nrow=1)
ggsave(filename = "plots/Fig2.pdf", plot= fig2,device = "pdf",
       width = 12,height = 8,units = "in",bg="white")
ggsave(filename = "plots/Fig2.png", plot= fig2,device = "png",
       width = 12,height = 8,units = "in",bg="white")




##Figure3 ----

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  novelAllelesStats %>% 
    mutate(ethnicity=factor(ethnicity,
                            levels=c("All",
                                     "Unknown",
                                     "Caucasian",
                                     "African",
                                     "Asian", 
                                     "Hispanic/Latino"))
    ) %>% 
    ggplot(aes(x=ethnicity,fill=ethnicity,y=nNovelPerDonors))+
    geom_bar(stat="identity",position = position_dodge(),
             color="grey38",linewidth=0.15)+
    scale_fill_manual(values=ethnicity_colors,drop = FALSE)+
    labs(x= "Ethnicity",y="# Novel Alleles Per Donor")+
    theme(legend.position = "bottom",
          strip.background =element_rect(fill="white"),
          axis.title.x = element_blank(),
          legend.title=element_text(size=14), 
          legend.text=element_text(size=12), 
          plot.margin = unit(c(3,5.5,6.5,6.5), "points")) + theme(legend.position="bottom",
                                                                  legend.box.margin = margin(0, 0, 0, 0))+
    scale_colour_discrete(guide = "none")+
    guides(fill = guide_legend(nrow = 1))+
    labs(fill="Ethnicity    ")
)

fig3<-
  plot_grid(
    plot_grid(
      plot_grid(fig3a,fig3b,fig3c,
                nrow=3,
                labels=c("A","B","C")),
      fig3d,
      nrow=1,
      labels=c("","D")
    ),
    legend,
    nrow=2,
    rel_heights = c(1,0.05)
  )

ggsave(filename = "plots/Fig3.pdf", plot= fig3,device = "pdf",
       width = 16,height = 8,units = "in",bg = "white")
ggsave(filename = "plots/Fig3.png", plot= fig3,device = "png",
       width = 16,height = 8,units = "in",bg = "white")
