setwd("/mnt/data/BMOHAMED/mutation_analysis/")

library(tidyverse)

df <- read_csv("IMR90_pos34.csv")

df %>%
  filter(aa != "Tyr") %>%
  mutate(group = str_sub(anticodon,1,1),
         wobble = factor(case_when(group == "A" ~ "A",
                                   group != "A" ~ "G/C/T"))) -> df2

df2 %>%
  ggplot(aes(x = wobble, y = pos34, fill = wobble))+
  geom_violin()+
  xlab("Wobble Position")+
  ylab("Mutation frequency")+
  theme_bw(base_size = 16)+
  ggtitle("A34I mutation rate: IMR90 Empty vs Ras")

ggsave(filename = "A34I_mutationrate_IMR90.png", bg = "white", width = 7, height = 7, dpi = 600)  

df2 %>%
  filter(wobble == "A") %>%
  mutate(sample_group = factor(case_when(Sample == "Empty1" | Sample == "Empty2" |Sample == "Empty3" ~ "Empty",
                                         Sample == "Ras1" | Sample == "Ras2" ~ "Ras"))) %>%
  filter(!(is.na(sample_group)))-> df3

df3 %>%
  group_by(anticodon, sample_group) %>%
  summarise(mean_A34I = mean(pos34)) %>%
  ggplot(aes(x = sample_group, y = mean_A34I, fill = sample_group))+
  geom_violin()+
  geom_boxplot() -> issoacceptor_plot

p <- wilcox.test(data = df3, pos34 ~ sample_group)

df3 %>%
  #group_by(anticodon, sample_group, aa) %>%
  #summarise(mean_A34I = mean(pos34)) %>%
  ggplot(aes(x = factor(aa), y = pos34, colour = sample_group))+
  geom_boxplot()+
  xlab("Isoacceptor")+
  ylab("Mutation frequency")+
  theme_bw(base_size = 16)+
  ggtitle("A34I mutation rates at known isoacceptors")+
  geom_point()


ggsave(filename = "A34I_Isoacceptor_IMR90.png", bg = "white", width = 7, height = 7, dpi = 600)  

wilcox.test(df3$pos34[df3$aa == "Ala" & df3$sample_group == "Empty"],
            df3$pos34[df3$aa == "Ala" & df3$sample_group == "Ras"])



##########################################################################################################################


df <- read_csv("FvS_pos34.csv")

df %>%
  filter(aa != "Tyr") %>%
  mutate(group = str_sub(anticodon,1,1),
         wobble = factor(case_when(group == "A" ~ "A",
                                   group != "A" ~ "G/C/T"))) -> df2

df2 %>%
  ggplot(aes(x = wobble, y = pos34, fill = wobble))+
  geom_violin()+
  xlab("Wobble Position")+
  ylab("Mutation frequency")+
  theme_bw(base_size = 16)+
  ggtitle("A34I mutation rate: Fed vs Starved")

ggsave(filename = "A34I_mutationrate_FvS.png", bg = "white", width = 7, height = 7, dpi = 600)  

df2 %>%
  filter(wobble == "A") %>%
  mutate(sample_group = factor(case_when(Sample == "Fed1" | Sample == "Fed2" | Sample == "Fed5" ~ "Fed",
                                         Sample == "Starved1" | Sample == "Starved2" | Sample == "Starved5" ~ "Starved"))) %>%
  filter(!(is.na(sample_group)))-> df3

df3 %>%
  group_by(anticodon, sample_group) %>%
  summarise(mean_A34I = mean(pos34)) %>%
  ggplot(aes(x = sample_group, y = mean_A34I, fill = sample_group))+
  geom_violin()+
  geom_boxplot() -> issoacceptor_plot

p <- wilcox.test(data = df3, pos34 ~ sample_group)

df3 %>%
  #group_by(anticodon, sample_group, aa) %>%
  #summarise(mean_A34I = mean(pos34)) %>%
  ggplot(aes(x = factor(aa), y = pos34, colour = sample_group))+
  geom_boxplot()+
  xlab("Isoacceptor")+
  ylab("Mutation frequency")+
  theme_bw(base_size = 16)+
  ggtitle("A34I mutation rates at known isoacceptors")+
  geom_point()


ggsave(filename = "A34I_Isoacceptor_FvS.png", bg = "white", width = 7, height = 7, dpi = 600)  

wilcox.test(df3$pos34[df3$aa == "Ala" & df3$sample_group == "Fed"],
            df3$pos34[df3$aa == "Ala" & df3$sample_group == "Starved"])

wilcox.test(df3$pos34[df3$aa == "Val" & df3$sample_group == "Fed"],
            df3$pos34[df3$aa == "Val" & df3$sample_group == "Starved"])

##########################################################################################################################


df <- read_csv("bcat_pos34.csv")

df %>%
  filter(aa != "Gly") %>%
  mutate(group = str_sub(anticodon,1,1),
         wobble = factor(case_when(group == "A" ~ "A",
                                   group != "A" ~ "G/C/T"))) -> df2

df2 %>%
  ggplot(aes(x = wobble, y = pos34, fill = wobble))+
  geom_violin()+
  xlab("Wobble Position")+
  ylab("Mutation frequency")+
  theme_bw(base_size = 16)+
  ggtitle("A34I mutation rate: Null vs bcat/myc")

ggsave(filename = "A34I_mutationrate_bcat.png", bg = "white", width = 7, height = 7, dpi = 600)  

df2 %>%
  filter(wobble == "A") %>%
  mutate(sample_group = factor(case_when(Sample == "WT1_null" | Sample == "WT2_null" | Sample == "WT4_null" | Sample == "WT5_null" ~ "Null",
                                         Sample == "bcat_myc2" | Sample == "bcat_myc3" | Sample == "bcat_myc4" | Sample == "bcat_myc5" ~ "bcat_myc"))) %>%
  filter(!(is.na(sample_group)))-> df3

df3 %>%
  group_by(anticodon, sample_group) %>%
  summarise(mean_A34I = mean(pos34)) %>%
  ggplot(aes(x = sample_group, y = mean_A34I, fill = sample_group))+
  geom_violin()+
  geom_boxplot() -> issoacceptor_plot

p <- wilcox.test(data = df3, pos34 ~ sample_group)

df3 %>%
  #group_by(anticodon, sample_group, aa) %>%
  #summarise(mean_A34I = mean(pos34)) %>%
  ggplot(aes(x = factor(aa), y = pos34, colour = sample_group))+
  geom_boxplot()+
  xlab("Isoacceptor")+
  ylab("Mutation frequency")+
  theme_bw(base_size = 16)+
  ggtitle("A34I mutation rates at known isoacceptors")+
  geom_point()


ggsave(filename = "A34I_Isoacceptor_FvS.png", bg = "white", width = 7, height = 7, dpi = 600)  

wilcox.test(df3$pos34[df3$aa == "Ser" & df3$sample_group == "Null"],
            df3$pos34[df3$aa == "Ser" & df3$sample_group == "bcat_myc"])

wilcox.test(df3$pos34[df3$aa == "Val" & df3$sample_group == "Fed"],
            df3$pos34[df3$aa == "Val" & df3$sample_group == "Starved"])

##########################################################################################################################


df <- read_csv("MDM2_pos34.csv")

df %>%
  filter(aa != "Gly") %>%
  mutate(group = str_sub(anticodon,1,1),
         wobble = factor(case_when(group == "A" ~ "A",
                                   group != "A" ~ "G/C/T"))) -> df2

df2 %>%
  ggplot(aes(x = wobble, y = pos34, fill = wobble))+
  geom_violin()+
  xlab("Wobble Position")+
  ylab("Mutation frequency")+
  theme_bw(base_size = 16)+
  ggtitle("A34I mutation rate: Null vs MDM2")

ggsave(filename = "A34I_mutationrate_MDM2.png", bg = "white", width = 7, height = 7, dpi = 600)  

df2 %>%
  filter(wobble == "A") %>%
  mutate(sample_group = factor(case_when(Sample == "M_MDM2_Null1" | Sample == "M_MDM2_Null2" | Sample == "M_MDM2_Null3" | Sample == "M_MDM2_Null4" | Sample == "M_MDM2_Null5" ~ "Null",
                                         Sample == "M_MDM2_Cre3" | Sample == "M_MDM2_Cre4" | Sample == "M_MDM2_Cre5" ~ "MDM2"))) %>%
  filter(!(is.na(sample_group)))-> df3

df3 %>%
  group_by(anticodon, sample_group) %>%
  summarise(mean_A34I = mean(pos34)) %>%
  ggplot(aes(x = sample_group, y = mean_A34I, fill = sample_group))+
  geom_violin()+
  geom_boxplot() -> issoacceptor_plot

p <- wilcox.test(data = df3, pos34 ~ sample_group)

df3 %>%
  #group_by(anticodon, sample_group, aa) %>%
  #summarise(mean_A34I = mean(pos34)) %>%
  ggplot(aes(x = factor(aa), y = pos34, colour = sample_group))+
  geom_boxplot()+
  xlab("Isoacceptor")+
  ylab("Mutation frequency")+
  theme_bw(base_size = 16)+
  ggtitle("A34I mutation rates at known isoacceptors")+
  geom_point()


ggsave(filename = "A34I_Isoacceptor_MDM2.png", bg = "white", width = 7, height = 7, dpi = 600)  

wilcox.test(df3$pos34[df3$aa == "Ser" & df3$sample_group == "Null"],
            df3$pos34[df3$aa == "Ser" & df3$sample_group == "MDM2"])


