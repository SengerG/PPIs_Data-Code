# Proteomic Data

# Downloading proteomic data used in this study

# CPTAC COAD cohort proteomic data
# Download proteomic data for normal samples
url = "http://linkedomics.org/cptac-colon/Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__PNNL_Normal_TMT_UnsharedLogRatio.cct"
# Specify destination where the data should be saved
destfile = "/path/to/destination/folder/Name.cct"
download.file(url, destfile = destfile)
COAD_cohort_normal = read.delim("/path/to/destination/folder/Name.cct")

# Download proteomic data for tumor samples
url = "http://linkedomics.org/cptac-colon/Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__PNNL_Tumor_TMT_UnsharedLogRatio.cct"
# Specify destination where the data should be saved
destfile = "/path/to/destination/folder/Name.cct"
download.file(url, destfile = destfile)
COAD_cohort_tumor = read.delim("/path/to/destination/folder/Name.cct")
 
# CPTAC HCC cohort proteomic data
url = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867419310037-mmc1.xlsx"
# Specify destination where the data should be saved
destfile = "/path/to/destination/folder/Name.xlsx"
download.file(url, destfile = destfile)
HCC_cohort = read.xlsx("/path/to/destination/folder/Name.xlsx", sheet = 4)


# CPTAC LUAD cohort proteomic data
url = "https://www.cell.com/cms/10.1016/j.cell.2020.06.013/attachment/49a46b71-468b-45d1-826a-721fa734eff0/mmc3.xlsx"
# Specify destination where the data should be saved
destfile = "/path/to/destination/folder/Name.xlsx"
download.file(url, destfile = destfile)
LUAD_cohort = read.xlsx("/path/to/destination/folder/Name.xlsx", sheet = 2, startRow = 3)


# Proteomic data for TCGA projects (for CPTAC COREAD, OV, BRCA cohorts) were directly downloaded via CPTAC Aspera Connect Client


# A - Statistical analysis on proteomic data

# A.1 - Differential proteomic analysis

path = "~/RData/ProteomeData" # Folder that stores proteome quantification data (processed, filtered, normalized) for tumor and normal samples

files = list.files(path = path, pattern = "^Proteome_Normal_vs_Tumor", all.files = TRUE, full.names = TRUE, recursive = FALSE)
for(i in 1:length(files)){
  load(files[i])
  tumor.cols = which(grepl(pattern = "tumor", names(Proteome.df)))
  normal.cols = which(grepl(pattern = "normal", names(Proteome.df)))
  ans <- apply(Proteome.df, 1, function(row) unlist(wilcox.test(as.numeric(row[tumor.cols[1]:tail(tumor.cols, n=1)]),as.numeric(row[normal.cols[1]:tail(normal.cols, n=1)]))[c("p.value","statistic")])) #Wilcoxon test
  ans.t = t(ans)
  ans.t = as.data.frame(ans.t)
  ans.t$Normal.Median = apply(Proteome.df[,normal.cols[1]:tail(normal.cols, n=1)], 1, median, na.rm = T)
  ans.t$Tumor.Median = apply(Proteome.df[,tumor.cols[1]:tail(tumor.cols, n=1)], 1, median, na.rm = T)
  ans.t$LogFC = ans.t$Tumor.Median - ans.t$Normal.Median
  p = as.numeric(ans.t$p.value)
  p.adjust = p.adjust(p, method = "bonferroni", n = length(p)) #multiple testing correction
  ans.t$p.adjust = p.adjust
  ans.t$DE = "NO"
  ans.t$DE[ans.t$LogFC > 1 & ans.t$p.adjust <= 0.05] <- "UP"
  ans.t$DE[ans.t$LogFC < -1 & ans.t$p.adjust <= 0.05] <- "DOWN"
  s.path = "path/to/destination/folder/"
  save(ans.t,file = s.path)}

# A.2 - Association test between the differentially abundant proteins and protein complex subunits

# Libraries required
library(dplyr)
library(reshape2)

load("~/RData/HumanComplexes_subunits.RData") # Human protein complexes and subunits obtained from the CORUM database
path = "~/RData/Wilcoxon/"
files = list.files(path = path, pattern = "^wilcoxon", all.files = TRUE, full.names = TRUE, recursive = FALSE)
for(i in 1:length(files)){
  load(files[i])
  name = gsub(".RData","",strsplit(basename(files[i]),"_")[[1]][5])
  ans.t$Gene.name = as.character(rownames(ans.t))
  ans.t$Complex.Protein = ifelse(ans.t$Gene.name %in% Subunits_Info_pc_genes$Gene.name, "Subunit","Not.Subunit")
  ans.t$DE = "Not.DE"
  ans.t$DE[abs(ans.t$LogFC) > 1 & ans.t$p.adjust <= 0.05] <- "DE"
  sum = rename(count(ans.t, DE, Complex.Protein), Freq = n)
  table = matrix(sum$Freq, ncol = 2)
  colnames(table) = c("DE","Not.DE")
  rownames(table) = c("Not.complex","Complex")
  res = chisq.test(table)
  assign(name,res)}

# A.3 - Standard deviation in protein abundances across tumor and normal samples

path = "~/RData/ProteomeData" # Folder that stores proteome quantification data (processed, filtered, normalized) for tumor and normal samples
files = list.files(path = path, pattern = "^Proteome_Normal_vs_Tumor", all.files = TRUE, full.names = TRUE, recursive = FALSE)
for(i in 1:length(files)){
  load(files[i])
  case = strsplit(gsub(".RData","",basename(files[i])),"_")[[1]][5]
  tumor.cols = which(grepl(pattern = "tumor", names(Proteome.df))) # Tumor samples
  tumor = Proteome.df[,tumor.cols] # protein abundances across tumor samples
  normal.cols = which(grepl(pattern = "normal", names(Proteome.df))) # Normal samples
  normal = Proteome.df[,normal.cols] # protein abundances across normal samples
  # Starndard deviation in protein abundances across tumor and normal samples
  sd_normal = apply(normal,1,function(x) sd(x,na.rm = TRUE))
  sd_tumor = apply(tumor,1,function(x) sd(x,na.rm = TRUE))
  test = t.test(sd_tumor, sd_normal, alternative = "two.sided")
  assign(case,test)}

# A.4 - Protein level Spearman correlations

# Required library
library(Hmisc)

#Function for correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut])}

path = "~/RData/ProteomeData" # Folder that stores proteome quantification data (processed, filtered, normalized) for tumor and normal samples
files = list.files(path = path, pattern = "^Proteome_", all.files = TRUE, full.names = TRUE, recursive = FALSE)
for(i in 1:length(files)){
  assign("data",get(load(files[i])))
  name = gsub(".RData","",strsplit(basename(files[i]),"_")[[1]][length(strsplit(basename(files[i]),"_")[[1]])])
  data = data[,grepl(pattern = "tumor",names(data))]
  data = data[which(rowMeans(!is.na(data)) > 0.5), ] # Remove proteins that do not have quantification value for half of the tumor samples to prevent to produce NaNs
  mm = t(data)
  res2 <- rcorr(as.matrix(mm), type = "spearman")
  Corr.mat = flattenCorrMatrix(res2$r, res2$P)
  path.tosave = "path/to/destination/folder/"
  save(Corr.mat, file = path.tosave)}


# B - Classification of protein-protein interactions and statistical analysis

# B.1 - Calculation of stoichiometric ratio

# Stoichiometric ratio between possible pairs among the proteins covered by each cohort was caclulated separately.
load("~/RData/Stoichiometry_genenames_PDB.RData")
stoi_proteins = unique(as.character(Stoic.table.human$Gene.name))

# Correlation data
path = "~/RData/Correlations" # Path to the folder where correlations across tumor samples are stored
files = list.files(path = path, pattern = "^Correlation", recursive = FALSE, full.names = TRUE)
for (i in 1:length(files)) {
  load(files[i])
  name = gsub(".RData","",strsplit(basename(files[i]),"_")[[1]][length(strsplit(basename(files[i]),"_")[[1]])])
  # Data covered by PDB data
  covered.data = Corr.mat[as.character(Corr.mat$row) %in% stoi_proteins,]
  covered.data = covered.data[as.character(covered.data$column) %in% stoi_proteins,]
  row_proteins = unique(as.character(covered.data$row))
  case.df = c()
  for(r in row_proteins){
    proI.p = paste0("\\b",r,"\\b")
    proI.p.df = covered.data[as.character(covered.data$row) == r,]
    proI.pdb = Stoic.table.human[grepl(proI.p,as.character(Stoic.table.human$Gene.name)),]
    for (j in 1:length(rownames(proI.p.df))) {
      proII.p = paste0("\\b",proI.p.df$column[j],"\\b")
      pdb.ids = merge(proI.pdb,
                      Stoic.table.human[grepl(proII.p,as.character(Stoic.table.human$Gene.name)),], by = "PDB.ID") # PDB complexes both proteins found
      if(length(rownames(pdb.ids)) == 0){proI.p.df$Stoichiometric.Ratio[j] = "No.info"}
      else{
        pdb.ids$ratio = paste(unlist(lapply(pdb.ids$Chain.x, function(x) length(strsplit(x,";")[[1]]))),unlist(lapply(pdb.ids$Chain.y, function(x) length(strsplit(x,";")[[1]]))),sep = ":")
        pdb.ids = pdb.ids[!(duplicated(pdb.ids$ratio)),]
        if(dim(pdb.ids)[1] == 1){proI.p.df$Stoichiometric.Ratio[j] = pdb.ids$ratio}
        else{proI.p.df$Stoichiometric.Ratio[j] = paste(pdb.ids$ratio, collapse = ";")}}}
    case.df = rbind(case.df,proI.p.df)}
  case.df = case.df[case.df$Stoichiometric.Ratio != "No.info",]
  save(case.df, file = "path/to/destination/folder/")}

# B.2 - Calculation of co-occurrence frequency

# The number of complexes two proteins found together was calculated for each possible pair among subunits of complexes obtained from the CORUM database

load("~/RData/HumanComplexes_subunits.RData") # Human protein complexes and subunits obtained from the CORUM database
subunits = unique(Subunits_Info_pc_genes$Gene.name)
Final.df = c()
for(i in 1:length(subunits)){
  p = subunits[i]
  p.comp = human_proteincomplex[grepl(paste0('\\b',p,'\\b'),human_proteincomplex$subunits.Gene.name.),]
  if(i == length(subunits)) next
  for (j in (i+1):length(subunits)) {
    pII = subunits[j]
    pII.comp = human_proteincomplex[grepl(paste0('\\b',pII,'\\b'),human_proteincomplex$subunits.Gene.name.),]
    together = intersect(p.comp$ComplexID,pII.comp$ComplexID) # The number of complexes (A and B found)
    atleastone = union(p.comp$ComplexID,pII.comp$ComplexID) # The number of complexes (A or B found)
    row = c(p,pII,length(intersect(p.comp$ComplexID,pII.comp$ComplexID)), length(union(p.comp$ComplexID,pII.comp$ComplexID)))
    Final.df = rbind(Final.df,row)}}

# B.3 - Defining context-specific and general interactions

# Required libraries
library(readxlsb)

# BioPlex Interactome: Interactions and Baits for 293T and HCT116 Cells can be downloaded as follow
url = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867421004463-mmc1.zip"

# Specify destination where the data should be saved
destfile = "/path/to/destination/folder/Name.zip"
download.file(url, destfile = destfile)

baits_293T = read_xlsb("/path/to/destination/folder/Huttlin_BioPlex3_Table_S1.xlsb", sheet = 2, col_names = TRUE)
baits_HCT116 = read_xlsb("/path/to/destination/folder/Huttlin_BioPlex3_Table_S1.xlsb", sheet = 4, col_names = TRUE)

interactions_merged = read_xlsb("/path/to/destination/folder/Huttlin_BioPlex3_Table_S1.xlsb", sheet = 6, col_names = TRUE)
shared_baits = intersect(baits_293T$Bait.Symbol, baits_HCT116$Bait.Symbol) # the baits that have been targeted in both cell lines

int_A = interactions_merged[interactions_merged$SymbolA %in% shared_baits,]
int_B = interactions_merged[interactions_merged$SymbolB %in% shared_baits,]
all_interactions = rbind(int_A,int_B)
all_interactions$pair = paste(all_interactions$SymbolA,all_interactions$SymbolB, sep = "_")
all_interactions = all_interactions[!(duplicated(all_interactions$pair)),]

general_interactions = all_interactions[all_interactions$Shared == TRUE,]
contextspecific_interactions = all_interactions[all_interactions$Shared == FALSE,]

# B.4 - Defining competitive and cooperative interactions

# Functions
# Position converter for binding sites
position_converter <- function(x){
  bs = as.character(x)
  num_ch = nchar(bs)
  bs = substring(bs, 2, (num_ch-1))
  pos = strsplit(bs, ",")
  if(length(pos[[1]]) != 0){
    new_pos = c()
    for (i in 1:length(pos[[1]])) {
      p = strsplit(pos[[1]][i], "-")[[1]]
      if(length(p) == 1){item = as.character(p[[1]][1])}
      else if(length(p) > 1){item = seq(as.numeric(p[1]), as.numeric(p[2]))}
      new_pos = c(new_pos,item)}}
  else{new_pos = as.character(x)}
  return(new_pos)}

# Protein interaction interface information from the Interactome Insider
url = "http://interactomeinsider.yulab.org/downloads/interfacesHQ/H_sapiens_interfacesHQ.txt"

# After converting Uniprot IDs to gene names
load("~/RData/H_sapiens_interfacesHQ_genename.RData")

# Remove binary interactions that have binding site information only for one protein
HS_interfaces_genenames = HS_interfaces_genenames[as.character(HS_interfaces_genenames$P1_IRES) != "[]",]
HS_interfaces_genenames = HS_interfaces_genenames[as.character(HS_interfaces_genenames$P2_IRES) != "[]",]


proteins = union(HS_interfaces_genenames$P1,HS_interfaces_genenames$P2)

Summary.df.n = c()
for (pr in proteins) {
  print(which(proteins == pr))
  pr_int = HS_interfaces_genenames[HS_interfaces_genenames$P1 == pr | HS_interfaces_genenames$P2 == pr,]
  pr_partners = unique(union(as.character(pr_int$P1), as.character(pr_int$P2)))
  pr_partners = pr_partners[pr_partners != pr] # List of proteins that first protein interacts with
  for (i in 1:length(rownames(pr_int))) {
    B = as.character(unlist(pr_int[i,1:2]))
    B = B[which(B != pr)]
    B_int = HS_interfaces_genenames[HS_interfaces_genenames$P1 == B | HS_interfaces_genenames$P2 == B,]
    B_partners = unique(union(as.character(B_int$P1), as.character(B_int$P2)))
    B_partners = B_partners[B_partners != B] # List of proteins that second protein interacts with
    shared_partners = intersect(pr_partners, B_partners)
    if(length(shared_partners) == 0)next
    # Finding interaction sites
    for(sh in shared_partners){
      pr_sh = HS_interfaces_genenames[(HS_interfaces_genenames$P1 == pr & HS_interfaces_genenames$P2 == sh) |
                                        (HS_interfaces_genenames$P2 == pr & HS_interfaces_genenames$P1 == sh),]
      B_sh = HS_interfaces_genenames[(HS_interfaces_genenames$P1 == B & HS_interfaces_genenames$P2 == sh) |
                                       (HS_interfaces_genenames$P2 == B & HS_interfaces_genenames$P1 == sh),]
      pr_bs = position_converter(pr_sh[1,which(pr_sh[1,] == sh, arr.ind = FALSE) + 2])
      B_bs = position_converter(B_sh[1,which(B_sh[1,] == sh, arr.ind = FALSE) + 2])
      overlap.pos = sum((pr_bs %in% B_bs) == TRUE) # Number of residues in the overlapping binding site
      line = c(pr,B,sh,overlap.pos,(length(pr_bs) + length(B_bs) - overlap.pos),c.info)
      Summary.df.n = rbind(Summary.df.n, line)}}}

colnames(Summary.df.n) = c("row","column","Shared.Partner","Number.Overlap.Positions","Total.Length","Complex.Info")
Summary.df.n = as.data.frame(Summary.df.n, stringasfactor = FALSE)

# Remove duplicated pairs

Summary.df.n$Pair = paste(Summary.df.n$row,Summary.df.n$column,sep = "_")
pairs = as.character(unique(Summary.df.n$Pair))
for (i in 1:length(pairs)) {
  pair = pairs[i]
  pair_reverse = paste(strsplit(pair,"_")[[1]][2],strsplit(pair,"_")[[1]][1], sep = "_")
  reverse = pairs[pairs == pair_reverse]
  if(length(reverse) > 0){pairs = pairs[pairs != pair_reverse]}}
Summary.df.n = Summary.df.n[Summary.df.n$Pair %in% pairs,]
Summary.df.n$Group = paste(Summary.df.n$row,Summary.df.n$column,Summary.df.n$Shared.Partner,sep = "_")
Summary.df.n = Summary.df.n[!(duplicated(Summary.df.n$Group)),]
Summary.df.n$JaccardIndex = Summary.df.n$Number.Overlap.Positions / Summary.df.n$Total.Length

# B.5 - Transient and Permanent interactions

load("~/RData/Transient_Permanent_interactions.RData") # Transient and permanent interactions obtained from Block et al. (2006) and Mintseris & Weng (2003)

data = rbind(Block_covered[,2:4],Mintseris_covered[,c(4,2,3)])
data$Pair = paste(data$ChainI,data$ChainII,sep = "_")
data = data[!duplicated(data$Pair),]
# To remove reverse pairs
pairs = as.character(unique(data$Pair))
for (i in 1:length(pairs)) {
  pair = pairs[i]
  pair_reverse = paste(strsplit(pair,"_")[[1]][2],strsplit(pair,"_")[[1]][1], sep = "_")
  reverse = pairs[pairs == pair_reverse]
  if(length(reverse) > 0){pairs = pairs[pairs != pair_reverse]}}
data = data[data$Pair %in% pairs,]



# B.6 - Linear regression analysis

load("~/RData/Stoichiometricratios_alldata.RData") # Stoichiometric ratio information for each pair covered by PDB and proteomics data of the corresponding cohort(study)
cases = unique(Final.df$CPTACdata)
path = "~/RData/ProteomeData" # Folder that stores proteome quantification data (processed, filtered, normalized) for tumor and normal samples
lm.data = c() # To store the linear regression analysis results
for(case in cases){
  case.df = Final.df[Final.df$CPTACdata == case,]
  case.df = case.df[case.df$Final.Ratio != "Changes",] #Remove protein pairs participating in different complexes with sometimes even sometimes uneven stoichiometric ratio
  pattern = paste0("*_Tumor_",case,".RData")
  file = list.files(path = path, pattern = pattern, all.files = TRUE, full.names = TRUE, recursive = FALSE) 
  assign("data",get(load(file)))
  data = data[,grepl(pattern = "tumor",names(data))]
  for(i in 1:nrow(case.df)){
    protI = as.character(case.df$row[i])
    protII = as.character(case.df$column[i])
    modeldata = as.data.frame(t(data[c(protI,protII),]))
    colnames(modeldata) = c("row","column")
    ratio = as.numeric(strsplit(strsplit(as.character(case.df$Stoichiometric.Ratio[i]),";")[[1]][1],":")[[1]])
    if(ratio[1] == ratio[2]){lmabundance = lm(row~column, data = modeldata)}
    else{
      max.prot = which.max(ratio)
      if(max.prot == 1){lmabundance = lm(column~row, data = modeldata)}
      else{lmabundance = lm(row~column, data = modeldata)}}
    sum = summary(lmabundance)
    case.df$Slope[i] = as.numeric(sum$coefficients[2,1])
    case.df$pvalue[i] = sum$coefficients[2,4]
    case.df$adj.Rsqaured[i] = sum$adj.r.squared}
  lm.data = rbind(lm.data, case.df)}


