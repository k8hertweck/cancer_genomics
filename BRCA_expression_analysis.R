#### analyzing breast cancer expression data from TCGA ####

# load packages
library(dplyr)
library(ggplot2)

# define exclusion subsetting
`%ni%` <- Negate(`%in%`) 

# read in saved data
fpkmGene <- read.table("data/targetGeneBRCA.csv")
# see all metadata
colnames(fpkmGene)
# view untransformed distribution
hist(fpkmGene$BRCA1) # left skewed
ggplot(fpkmGeneNolog, aes(BRCA1)) + 
  geom_histogram() +
  xlab("BRCA1 expression") +
  ggtitle("Distribution of BRCA1 expression") +
  theme_bw() + 
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16))
#ggsave("figures/BRCA1distributionUntrans.jpg")
hist(fpkmGene$BRCA2) # left skewed

# save untransformed data
fpkmGeneNolog <- fpkmGene

# log transform gene expression data
fpkmGene[1:2] <- fpkmGene[1:2] + 1 # add one pseudo count to all counts to remove zeros
fpkmGene[1:2] <- log2(fpkmGene[1:2]) # apply log2 transformation

## visualizing data distribution for variables of interest
# not variable or not reported: classification_of_tumor, last_known_disease_status, tumor_grade, progression_or_recurrence, disease_type
hist(fpkmGene$BRCA1) # normal
ggplot(fpkmGene, aes(BRCA1)) + 
  geom_histogram() +
  xlab("log2 BRCA1 expression") +
  ggtitle("Distribution of BRCA1 expression") +
  theme_bw() + 
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16))
#ggsave("figures/BRCA1distribution.jpg")
hist(fpkmGene$BRCA2) # normal
ggplot(fpkmGene, aes(BRCA2)) + 
  geom_histogram() +
  xlab("log2 BRCA2 expression") +
  ggtitle("Distribution of BRCA2 expression") +
  theme_bw() + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/BRCA2distribution.jpg")
table(fpkmGene$shortLetterCode) 
# 113 NT= normal tissue, 1102 TP= primary tumor, 7 TM= metastatic
ggplot(fpkmGene, aes(shortLetterCode)) +
  geom_bar() +
  xlab("tissue type") +
  ylab("number of samples") +
  ggtitle("Number of samples by tissue type") + 
  scale_x_discrete(labels=c("NT" = "normal", "TM" = "metastatic", "TP" = "primary tumor")) +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16), title = element_text(size = 16))
#ggsave("figures/BRCAtissueType.jpg")

table(fpkmGene$tumor_stage) # same: table(fpkmGene$subtype_Converted.Stage)
table(fpkmGene$subtype_AJCC.Stage)
table(fpkmGene$vital_status) #table(fpkmGene$subtype_Vital.Status) has missing data
table(fpkmGene$morphology)
table(fpkmGene$subtype_Metastasis) # same: table(fpkmGene$subtype_Metastasis.Coded)

# triple negative indicators
table(fpkmGene$subtype_ER.Status)
table(fpkmGene$subtype_PR.Status)
table(fpkmGene$subtype_HER2.Final.Status)

## initial data subsets
# assess number of samples with BRCA1 expression
meta <- fpkmGene %>%
  filter(shortLetterCode == "TM")
norm <- fpkmGene %>%
  filter(shortLetterCode == "NT")
tum <- fpkmGene %>%
  filter(shortLetterCode == "TP")

#### Compare BRCA1 expression between normal and BC patients ####
# all "normal" samples are paired with a Bca sample
normVcancer <- rbind(norm, tum)
# unpaired, all data
t.test(BRCA1 ~ shortLetterCode, data = normVcancer) # p=2.2e-16
table(normVcancer$shortLetterCode) # 113 normal, 1102 tumor
ggplot(normVcancer, aes(shortLetterCode, BRCA1)) + 
  geom_boxplot() +
  ylab("log2 BRCA1 expression") +
  xlab("tissue type (unpaired)") +
  scale_x_discrete(labels=c("NT" = "normal", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() + 
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16))
#ggsave("figures/BRCA1unpaired.jpg")
# paired
# create paired sample dataset
normVcancerPaired <- normVcancer %>%
  filter(bcr_patient_barcode %in%
           norm$bcr_patient_barcode)
# find patients with more than 1 tumor sample to pair
normVcancerPaired %>% 
  group_by(patient) %>%
  tally() %>%
  filter(n > 2)
# find patients with only 1 normal sample
normVcancerPaired %>%
  group_by(patient) %>%
  tally() %>%
  filter(n == 1)
# only one sample: TCGA-BH-A0BS 
filter(normVcancerPaired, patient == "TCGA-BH-A0BS")
grep("TCGA-BH-A0BS", normVcancerPaired$patient) # remove 67 (no tumor to pair) 
# three or four samples: TCGA-A7-A0DB, TCGA-A7-A0DC, TCGA-A7-A13E
filter(normVcancerPaired, patient == "TCGA-A7-A0DB")
grep("TCGA-A7-A0DB", normVcancerPaired$patient) # remove 182 and 216 (extra tumor)
filter(normVcancerPaired, patient == "TCGA-A7-A0DC")
grep("TCGA-A7-A0DC", normVcancerPaired$patient) # remove 145 (tumor missing metadata)
filter(normVcancerPaired, patient == "TCGA-A7-A13E")
grep("TCGA-A7-A13E", normVcancerPaired$patient) # remove 221 and 226 (extra tumor)
# remove extra samples
normVcancerPaired <- normVcancerPaired[-c(67, 182, 216, 145, 221, 226),]
# sort by barcode
normVcancerPaired <- normVcancerPaired[order(normVcancerPaired$patient),] 
# summarize sample counts
table(normVcancerPaired$shortLetterCode) # 112 NT, 112 TP
# perform t test (paired)
t.test(BRCA1 ~ shortLetterCode, paired=TRUE, data = normVcancerPaired) # 1.003e-15
ggplot(normVcancerPaired, aes(shortLetterCode, BRCA1)) + 
  ylab("log2 BRCA1 expression") +
  xlab("tissue type (paired samples)") +
  scale_x_discrete(labels=c("NT" = "normal", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16))
#ggsave("figures/BRCA1paired.jpg")

# is gene expression higher in deceased individuals?
vital <- fpkmGene %>%
  filter(shortLetterCode == "TP") %>%
  filter(!is.na(shortLetterCode))
table(vital$vital_status) # 947 alive, 154 dead
# BRCA1 and vital status 
t.test(BRCA2 ~ vital_status, data=vital, alternative="less") # 0.03678
ggplot(vital, aes(vital_status, BRCA1)) + 
  ylab("log2 BRCA1 expression") +
  xlab("vital status") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16))
#ggsave("figures/BRCA1vital.jpg")

# compare expression between TP and MT
# combine tumor and metastasis samples
tumVmetaAll <- rbind(tum, meta)
table(tumVmetaAll$shortLetterCode) # 7 TM 1102 TP
# perform t test (unpaired data)
t.test(BRCA1 ~ definition, data = tumVmetaAll) # 0.7794
ggplot(tumVmetaAll, aes(definition, BRCA1)) + 
  ylab("log2 BRCA1 expression") +
  xlab("tissue type") +
  scale_x_discrete(labels=c("TM" = "metastatic", "TP" = "primary tumor")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/BRCA1unpairedMetastasis.jpg")

#### Compare BRCA1 expression with clinical stages of BC patients ####
table(fpkmGene$subtype_AJCC.Stage) # missing data are "[Not Available]"
table(fpkmGene$tumor_stage) # missing data are "not reported"
# compare AJCC stages
AJCC <- tum %>%
  filter(subtype_AJCC.Stage != "[Not Available]") %>%
  filter(subtype_AJCC.Stage != "Stage X")
table(AJCC$subtype_AJCC.Stage)
# perform ANOVA
summary(aov(BRCA1 ~ subtype_AJCC.Stage, data = AJCC)) # 0.036
ggplot(AJCC, aes(subtype_AJCC.Stage, BRCA1)) + 
  ylab("log2 BRCA1 expression") +
  xlab(" AJCC stage") +
  scale_x_discrete(labels=c("Stage I" = "I", "Stage IA" = "IA", "Stage IB" = "IB", "Stage II" = "II", "Stage IIA" = "IIA", "Stage IIB" = "IIB", "Stage III" = "III", "Stage IIIA" = "IIIA", "Stage IIIB" = "IIIB", "Stage IIIC" = "IIIC", "Stage IV" = "IV")) +
  geom_boxplot() +
  theme_bw() 
#ggsave("figures/BRCA1.AJCC.jpg")

#### Compare BRCA1/2 expression together in normal vs. tumor ####
normPaired <- normVcancerPaired %>% 
  filter(shortLetterCode == "NT")
tumPaired <- normVcancerPaired %>%
  filter(shortLetterCode == "TP")
# linear regression (only normal)
norm.mod <- lm(BRCA1 ~ BRCA2, data=normPaired)
summary(norm.mod) # p=<2e-16, R2=0.4686
# linear regression (only TNBC)
tum.mod <- lm(BRCA1 ~ BRCA2, data=tumPaired)
summary(tum.mod) # p=1.22e-14, R2=0.1644
# plot both together
ggplot(normVcancerPaired, aes(BRCA1, BRCA2, col=shortLetterCode)) +
  geom_point() +
  ylab("log2 BRCA2 expression") +
  xlab("log2 BRCA1 expression") +
  theme_bw() +
  geom_smooth(data=subset(tumPaired, shortLetterCode == "TP"), method = "lm", se = FALSE) +
  geom_smooth(data=subset(normPaired, shortLetterCode == "NT"), method = "lm", se = FALSE) +
  scale_color_discrete(name="tissue type",
                      breaks=c("NT", "TP"),
                      labels=c("normal", "tumor"))
#ggsave("figures/BRCA12.jpg")