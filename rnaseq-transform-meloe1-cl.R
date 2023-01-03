#rm(list = ls())
library(data.table)
library(dplyr)
library(plyr)
library(forcats)
library(reshape2)
library(ggplot2)
library(stringr)

taa_set1 = c('BIRC5','CTAG1B','PRAME','SSX2','WT1');taa_set1
taa_set2 = c('BIRC5','E6_hpv16','E7_hpv16','MAGEA4','PRAME');taa_set2
#'E1_hpv16', 'E2_hpv16', 'E4_hpv16', 'E5_hpv16',
taa_set2.1 = c('BIRC5','DEPDC1','E2F7','PRAME','SOX30');taa_set2.1

all_taa = unique(c(taa_set1,taa_set2,taa_set2.1,'SOX30','NEK2'))

#hkg = c('C1orf43','CHMP2A', 'EMC7', 'GPI', 'PSMB2', 'PSMB4', 'RAB7A', 'REEP5', 'SNRPD3', 'VC', 'VPS29')

hkg = c('ACTB','RPL19','B2M','RPLP0','GAPDH','LDHA', #High
        'PGK1','TUBB','SDHA','CLTC','G6PD','HPRT1', #Medium #'POLR2A' -- Remove it 
        'ABCF1','GUSB','POLR1B','TBP','ALAS1') #Low

##### HKG genes in Shweta

shwlog2 = fread("~/Box/schavan/Vio/Meloe1/celllines/gene_level_TPM_kallisto_2_with_genenames.tsv",stringsAsFactors=FALSE) %>% 
  # filter(DupGeneName == 0) %>% 
  # select(-c(DupGeneName)) %>% 
  mutate_if(., is.numeric, funs(log2(1 + .))) %>% #### LOG2 transformation
  arrange(Gene)
 head(shwlog2)

#shw_trans_level = fread("/Users/schavan/Box/schavan/Cat/CASKI/abundance.tsv") %>% filter(., Gene == "E5_hpv16")
  
  head(shwlog2); dim(shwlog2) #35658
  
  gep = shwlog2 %>% select(-ENST)
  length(unique(gep$Gene)) #35637
  
  gex_genelevel = melt(gep) %>% 
    group_by(variable, Gene) %>% 
    dplyr::summarise(genelevel = sum(value, na.rm = TRUE))
  
  gex_genelevel_dcast = dcast(gex_genelevel, Gene ~ variable, value.var = "genelevel")
  dim(gex_genelevel_dcast) #35637
  length(unique(gex_genelevel_dcast$Gene)) #35637

  #gex_genelevel_dcast_fil = filter(gex_genelevel_dcast, !Gene %like% "\\.")
  #dim(gex_genelevel_dcast_fil) #28443  
  gex_genelevel_dcast_fil = gex_genelevel_dcast
  #names(gex_genelevel_dcast_fil)  

  
#names(gex_genelevel_dcast_fil)  
shw0 = gex_genelevel_dcast_fil  
colnames(shw0)
dim(shw0)

#sort(colnames(shw0))

write.table(shw0,"~/Box/schavan/Vio/Meloe1/celllines/Raw_gene_level_Log2Plus1TPM_kallisto_GEP_Meloe1.tsv",sep="\t", row.names = F, quote = F, append = F)


################################################################################################
##### HKG QC heatmap
shw.hkg = shw0 %>% filter(Gene %in% hkg) %>% melt() %>% arrange(variable); head(shw.hkg); tail(shw.hkg)

shw.hkg$Gene <- factor(shw.hkg$Gene, levels = unique(hkg),ordered=TRUE)
shw.hkg$variable <- factor(shw.hkg$variable, levels = unique(shw.hkg$variable), ordered =TRUE)

hm = shw.hkg %>% select(Gene, variable, value); dim(hm)
ggp <- ggplot(hm, aes(x=Gene, y=variable, fill = value )) + geom_tile(color="white") + 
  scale_fill_gradient(low = "white", high = "navy") + coord_equal() + theme_minimal() + 
  theme(axis.text.x=element_text(angle=90,hjust=1)) + ggtitle('Pre-HKG-Data'); 
ggp

###############

shw.taa = shw0 %>% filter(Gene %in% all_taa) %>%
  melt() %>% group_by(variable) %>% mutate(Group = "shw")

hm = shw.taa %>% select(Gene, variable, value); dim(hm)
ggp <- ggplot(hm, aes(x=Gene, y=variable, fill = value )) + geom_tile(color="white") + 
  scale_fill_gradient(low = "white", high = "navy") + coord_equal() + theme_minimal() + 
  theme(axis.text.x=element_text(angle=90,hjust=1)) + ggtitle('Pre-HKG-Data-TAA'); 
ggp


################################################################################################
##### ACTUAL NORMALIZATION
################################################################################################

shw.hkg.med.ptn = shw0 %>% filter(Gene %in% hkg) %>%
  melt() %>% 
  group_by(variable) %>% 
  dplyr::summarise(ptn.med.hkg = median(value)) %>% arrange(variable)
head(shw.hkg.med.ptn)
dim(shw.hkg.med.ptn)

shw.hkg.norm = shw0 %>% 
  filter(!is.na(Gene)) %>%
  as.data.frame() %>% 
  melt(.) %>% 
  ungroup %>%
  mutate(hkg.norm.denom = plyr::mapvalues(variable,shw.hkg.med.ptn$variable,shw.hkg.med.ptn$ptn.med.hkg)
         )

head(shw.hkg.norm)
tail(shw.hkg.norm)

shw.hkg.norm$hkg.norm.value = as.numeric(as.character(shw.hkg.norm$value))/as.numeric(as.character(shw.hkg.norm$hkg.norm.denom))
head(shw.hkg.norm)
tail(shw.hkg.norm)

hm.norm = shw.hkg.norm %>% select(Gene, variable, hkg.norm.value) %>% filter(Gene %in% all_taa); dim(hm.norm)
ggp <- ggplot(hm.norm, aes(x=Gene, y=variable, fill = hkg.norm.value )) + geom_tile(color="white") + 
  scale_fill_gradient(low = "white", high = "navy") + coord_equal() + theme_minimal() + 
  theme(axis.text.x=element_text(angle=90,hjust=1)) + ggtitle('HKG-Data-TAA'); 
ggp

shw.hkg.norm = shw.hkg.norm %>% select(Gene,variable,hkg.norm.value)
shw.hkg.med.ptn.norm = shw.hkg.norm %>% dcast(Gene~variable,value.var = 'hkg.norm.value')
dim(shw.hkg.med.ptn.norm)
head(shw.hkg.med.ptn.norm)

write.table(shw.hkg.med.ptn.norm,"~/Box/schavan/Vio/Meloe1/celllines/HKG_gene_level_Log2Plus1TPM_kallisto_GEP_Meloe1.tsv",sep="\t", row.names = F, quote = F, append = F)

################################################################################################
##### Make barcharts for TAA expression PRE vs. POST treatment
################################################################################################

ss = read.table("~/Box/schavan/Vio/Meloe1/celllines/SampleSheet.txt",stringsAsFactors = F,header = T)
head(ss)

as.data.frame(ss)
dim(ss)
ss = ss %>% select(patid,cancertype) %>% distinct()
uniqpatid = unique(ss$patid)

pdf("~/Box/schavan/Vio/Meloe1/celllines/Heatmap_Meloe1_All.pdf")
hm = shw.hkg.norm %>% 
  select(Gene, variable, hkg.norm.value) %>% 
  filter(Gene %like% "gi|14124996" | Gene %like% 'lcl'); 
dim(hm)
hm %>% dcast(.,Gene ~ variable)
unique(hm$Gene)

ggp <- ggplot(hm, aes(x=Gene, y=variable, fill = hkg.norm.value )) + geom_tile(color="white") + 
  scale_fill_gradient(low = "white", high = "navy") + coord_equal() + theme_minimal() + 
  theme(axis.text.x=element_text(angle=90,hjust=1)) + ggtitle('HKG-Data-Meloe1'); 
ggp
dev.off()
################################################################################################
##Frequency table
################################################################################################
hmt = hm %>% group_by(Gene) %>% filter(hkg.norm.value >0)
hmtdf = as.data.frame(table(hmt$Gene, hmt$variable)) 
N=4
hmtdf0 = hmtdf %>% group_by(Var1) %>% dplyr::summarise(freq = sum(Freq)) %>% mutate(total = N)
write.table(hmtdf0,"~/Box/schavan/Vio/Meloe1/celllines/FreqTavle_Meloe1_All_ORFs.txt",sep="\t",append = F,quote = F)
#%>% mutate(freq = rowsum(.))


################################################################################################
##### Make barcharts for TAA expression PRE vs. POST treatment
################################################################################################

cipher_hits = c('lcl|1:c594-235','lcl|1:c2090-1965','lcl|1:c2085-1903','lcl|1:c1866-1573','lcl|1:c1494-1339','lcl|1:c1393-1244',
'lcl|1:c137-3','lcl|1:736-912','lcl|1:729-797','lcl|1:644-778','lcl|1:438-689','lcl|1:2125-2154',
'lcl|1:2076-2153','lcl|1:1209-1370','lcl|1:1105-1194','lcl|1:1070-1108')
pdf("~/Box/schavan/Vio/Meloe1/Heatmap_Meloe1_Cipher_All.pdf")
hm.cipher = shw.hkg.norm %>% 
  select(Gene, variable, hkg.norm.value) %>% 
  filter(Gene %like% "gi|14124996" | Gene %in% cipher_hits); 
dim(hm.cipher)
ggp <- ggplot(hm.cipher, aes(x=Gene, y=variable, fill = hkg.norm.value )) + geom_tile(color="white") + 
  scale_fill_gradient(low = "white", high = "navy") + coord_equal() + theme_minimal() + 
  theme(axis.text.x=element_text(angle=90,hjust=1)) + ggtitle('HKG-Data-Meloe1-CIPHER'); 
ggp
dev.off()

