##### load library #####

###################################
# ERROR QUARANTINE #
# this will probably need to be changed when adding to nextflow
#  test4 <- test4 %>% mutate(name = str_remove(name, "./out_bam/")) %>%
#  mutate(name = str_remove(name, "_sorted.bam"))


#################################

#### notes ####

# going to need to rename the files at some point in the pipeline in the header of each fasta
# to include the sample name

#see here for genome coordinates:
# https://www.ncbi.nlm.nih.gov/genome/browse/#!/proteins/86693/757732%7CSevere%20acute%20respiratory%20syndrome%20coronavirus%202/

#setwd("C:/Users/m.bassett/Dropbox (UFL)/coding/scripts from max/snpeff/ncbi_dl/snpeff_out/individual")
setwd("C:/Users/m.bassett/Dropbox (UFL)/coding/scripts from max/snpeff/ncbi_dl/snpeff_out/individual")
library("gdata")
list.of.packages <- c("pheatmap", "gridExtra", "viridis", "splitstackshape", "tidyverse", "data.table", "stringr", "RColorBrewer", "fs", "tools", "ggplot2", "BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(trackViewer, GenomicFeatures)

if(!"trackViewer" %in% installed.packages()[,"Package"]){
  BiocManager::install("trackViewer")
}

lapply(c(list.of.packages, "trackViewer"), require, character.only = TRUE)


#### options ####

# see brewer pals
color_set = "GnBu" 
# false or vector with "non-syn", "other", "del", "syn", "frameshift", "ins", "stop" 
# combinations possible
hide_mutations = F
# intersect mutations present in another vcf: FLASE or relpath/*file_name*.vcf
ref_vcf = F
# intersect type: union, diff
intersect_type = "union"
# show nt changes
show_nt = F
# show amino acid changes, overwrites show_nt
show_as = F
# variant frequency threshold
freq_thres = 0.95
# false or provide nt positions c(start, stop)
show_part = F
# false or gene name  - will overwrite show_part
show_gene = F
# print density plot of mutation counts
print_dens = F
# focus only on the spike region?
spike_only = T
# subset to include ____ region? (T = yes) (each one you say T will be included)
# note: can delete later, going to subset each vcf that gets passed
# heres the idea, the earliest domain that contains a mutation will be first
# (i.e. if there is a mutation in orf1a it will appear first, if not, then spike is first)
# so make unique data frames then merge them together afterwards 
#can remove everything after set column because they will be the same for each region
subset_orf1a = T
subset_spike = T
subset_orf3a = T

# plot will be for protein? (T for yes, if not will be for DNA position)
protein_pos = T




### If I want to expand this, look into permutation tables for T/F values
# https://cran.r-project.org/web/packages/permutations/permutations.pdf
# https://www.r-bloggers.com/2021/05/learning-r-creating-truth-tables/
# https://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r


##### script #####

## define features
# read in GFF

#create a function to subset between named columns
select_cols_between <- function(df, col_start, col_end) {
  # Get the index of the starting and ending column
  start_index <- match(col_start, colnames(df))
  end_index <- match(col_end, colnames(df))
  
  # Select all columns between the starting and ending columns
  df[, start_index:end_index]
}

Gff_files = list.files(path = ".",pattern = "*.gff", recursive = F ,full.names = TRUE)
Gff_file = Gff_files[1]
Gff_file_spike = Gff_files[2]

spike_only = T

if (isTRUE(spike_only)){
  # Define the column names for the GFF3 file
  column_names <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  # Load the GFF3 file into a data frame
  gff3_df <- read.table(Gff_file_spike, sep = "\t", col.names = column_names, comment.char = "#", stringsAsFactors = FALSE)
  
  # Split the attributes column into separate columns
  gff3_df$attributes <- as.character(gff3_df$attributes)
  attributes_df <- data.frame(do.call("rbind", strsplit(gff3_df$attributes, ";")))
  colnames(attributes_df) <- c("key", "value", "note_name")
  attributes_df$value <- gsub("\"", "", attributes_df$value)
  attributes_df$value <- gsub(",", ";", attributes_df$value)
  
  # Merge the attributes data frame with the GFF3 data frame
  gff3_df <- cbind(gff3_df[1:8], attributes_df)
  
  gff_dataframe <- data.frame(gff3_df, stringsAsFactors = F)
  
  pattern <- "Note="
  pattern2 <- "Name="
  
  gff_dataframe$collect_note <- apply(gff_dataframe, 1, function(x) paste(grep(pattern, x, value = TRUE), collapse = ", "))
  gff_dataframe$collect_name <- apply(gff_dataframe, 1, function(x) paste(grep(pattern2, x, value = TRUE), collapse = ", "))
  
  
  gff_dataframe_subset <- gff_dataframe %>%
    relocate(collect_note, .before = NA.)
  gff_dataframe_subset <- gff_dataframe_subset %>%
    relocate(collect_name, .before = NA.)
  
  gff_dataframe_subset <- select_cols_between(gff_dataframe_subset, "seqid", "collect_name")
  
  gff_sub_final<- gff_dataframe_subset[c("start", "end", "note_name")]
  
  
  gff_sub_final$note_name <- gsub("Note=", "", gff_sub_final$note_name)
  gff_sub_final$note_name <- gsub("Name=", "", gff_sub_final$note_name)
  #testnew <- gff3_df[str_detect(gff3_df$note_name, "Note="), ]
  # Remove the original attributes column
  #gff3_df$attributes <- NULL
  
}  


# extract genes and relevant ranges
# subset gff based on gene
if (isFALSE(spike_only)){
  # Define the column names for the GFF3 file
  column_names <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  # Load the GFF3 file into a data frame
  gff3_df <- read.table(Gff_file, sep = "\t", col.names = column_names, comment.char = "#", stringsAsFactors = FALSE)
  
  # Split the attributes column into separate columns
  gff3_df$attributes <- as.character(gff3_df$attributes)
  attributes_df <- data.frame(do.call("rbind", strsplit(gff3_df$attributes, ";")))
  colnames(attributes_df) <- c("key", "value")
  attributes_df$value <- gsub("\"", "", attributes_df$value)
  attributes_df$value <- gsub(",", ";", attributes_df$value)
  
  # Merge the attributes data frame with the GFF3 data frame
  gff3_df <- cbind(gff3_df[1:8], attributes_df)
  
  gff_dataframe <- data.frame(gff3_df, stringsAsFactors = F)
  
  pattern <- "Note="
  pattern2 <- "Name="
  pattern3 <- "^gene="
  
  gff_dataframe$collect_note <- apply(gff_dataframe, 1, function(x) paste(grep(pattern, x, value = TRUE), collapse = ", "))
  gff_dataframe$collect_name <- apply(gff_dataframe, 1, function(x) paste(grep(pattern2, x, value = TRUE), collapse = ", "))
  gff_dataframe$collect_gene <- apply(gff_dataframe, 1, function(x) paste(grep(pattern3, x, value = TRUE), collapse = ", "))
  
  
  gff_dataframe_subset <- gff_dataframe %>%
    relocate(collect_note, .before = NA.)
  gff_dataframe_subset <- gff_dataframe_subset %>%
    relocate(collect_name, .before = NA.)
  gff_dataframe_subset <- gff_dataframe_subset %>%
    relocate(collect_gene, .before = NA.)
  
  
  
  gff_dataframe_subset <- select_cols_between(gff_dataframe_subset, "seqid", "collect_gene")
}

# get VCF files: NOTE, THEY ALL NEED TO BE IN THE SAME DIRECTORY
vcf_list = list.files(path = ".",pattern = "*.vcf", recursive = F ,full.names = TRUE)



n = length(vcf_list)
datalist = list()
datalist = vector("list", length = n)


for (i in 1:n) {
  
  # read in vcf
  vcf_file = fread(vcf_list[i])
  # create variant table
  variant_data = data.table()
  variant_data$position = vcf_file$POS
  vd_test <- as.data.frame(variant_data)
  test1 <- vcf_file
  test1 <- as.data.frame(test1)
  test2 <-  test1 %>% mutate(INFO = strsplit(INFO, "\\|")) %>%
    unnest(INFO) %>%
    group_by(POS) %>%
    mutate(row = row_number()) %>%
    mutate(across(where(is.character), ~ na_if(.,""))) %>%
    spread(row, INFO)
  test3 <- test2[colSums(!is.na(test2)) > 0]
  xxx <- sort(colnames(test3))
  xxx <- as.data.frame(xxx)
  xxx2 <- xxx[2,]
  test4 <- cbind(test3, name = xxx2)
  test4 <- test4 %>% mutate(name = str_remove(name, "./out_bam/")) %>%
    mutate(name = str_remove(name, "_sorted.bam"))
  test5 <- test4
  
  #test5 <- test4[!grepl("-", test4)]
  #test5 <- test5[!grepl("-", test5)]
  #checking subsetting to spike only
  
  
  if (isTRUE(subset_spike)){ 
    test6 <- test5 %>% filter(POS > 21562)
    test6 <- test6 %>% filter(POS < 25385)
    colnames(test6)[24] = "drop"
    test6 <- test6 %>%
      relocate(name, .before = drop)
    
    df_test6 <- test6[1:24]
    
  }
  if (isTRUE(subset_orf1a)){ 
    test7 <- test5 %>% filter(POS < 21563)
    colnames(test7)[24] = "drop"
    test7 <- test7 %>%
      relocate(name, .before = drop)
    
    df_test7 <- test7[1:24]
    
  }
  if (isTRUE(subset_orf3a)){ 
    test8 <- test5 %>% filter(POS > 25392)
    test8 <- test8 %>% filter(POS < 26220)
    colnames(test8)[24] = "drop"
    test8 <- test8 %>%
      relocate(name, .before = drop)
    
    df_test8 <- test8[1:24]
    
  }
  testxx = F
  ## was testing the intergenic regions, could be useful (there is one but its not annotated :( 
  if (isTRUE(testxx)){ 
    test6 <- overalldata %>% filter(POS > 21555)
    test6 <- test6 %>% filter(POS < 21563)
  }
  
  if (isTRUE(subset_spike & subset_orf1a)){
    df_test9 <- rbind(df_test6, df_test7)
  }
  if (isTRUE(subset_spike & subset_orf3a)){
    df_test9 <- rbind(df_test6, df_test8)
  }
  if (isTRUE(subset_spike & subset_orf1a & subset_orf3a)){
    df_test9 <- rbind(df_test6, df_test7, df_test8)
  }
  
  
  #finaldat <- rbind(test_blank, if(exists("test6")))
  
  dataframes <- ls(pattern = "^df_")
  completeddd <- do.call(rbind, mget(dataframes))
  
  
  
  
  #### HEY HERES AN IDEA, CAN DO dplyr GROUP BY THE POSITION AND THEN SUM THE CHANGES
  #test4 %>% mutate_all(vars(-group_cols(), case_when(grep("-")))  
  
  
  
  #df %>%  
  # mutate(group = case_when(grepl("Bl", b) ~ "Group1",
  #                         grepl("re", b, ignore.case = TRUE) ~"Group2"))
  datalist[[i]] <- completeddd
  
}

overalldata <- do.call(rbind, datalist)
colnames(overalldata)[19] = "DNA"
colnames(overalldata)[20] = "Protein"
overalldata <- overalldata %>% drop_na(Protein)

muttype <- unique(overalldata[19])
overalldata <- overalldata[1:24]
overalldata2 <- overalldata %>% select_if(~ !any(is.na(.)))
overall3 <- distinct(overalldata2)
colnames(overall3)[10] = "Effect"
colnames(overall3)[9] = "notsure"

rm(test1, test2, test3, test4, test5, test6, test7, test8, overalldata, overalldata2,
   df_test6, df_test7, df_test8, df_test9)

# want to split df?
#totalDF <- split(finDFmetaRNA, finDFmetaRNA$hexolig)

## define variables and colors
mutation_type_user = c(unique(overall3$Effect))
mut_shapes <- length(mutation_type_user)

overall <- as.data.frame(overall3)
colnames(overall)[11] = "DEFCON"
colnames(overall)[12] = "Gene"



# note: if want to add all of the individual reads, could add name, but that makes 1600+ entries
finalDF <- overall %>% 
  group_by(Gene, POS, REF, ALT, Effect, Protein, DNA) %>%
  summarize(n = n())


if (isTRUE(protein_pos)){ 
  finalDF_prot <- finalDF
  finalDF_prot$Protein <- gsub("p\\.", "", finalDF$Protein)
  
  # Define regular expression to split the mutation column
  split_regex <- "([[:alpha:]]{3})(\\d+)([[:alpha:]]*)"
  
  # Split the mutation column into three parts using the regular expression
  finalDF_prot[, c("aa_ref", "Position", "aa_mut")] <- str_match(finalDF_prot$Protein, split_regex)[, -1]
  
}


gff_sub_final_prot <- gff_sub_final[-1,]

gff_sub_final_prot$note_name <- gsub("gbkey=Prot", "SARS-CoV-2 Spike", gff_sub_final_prot$note_name)
gff_sub_final_prot$note_name <- gsub("SARS-CoV-like_Spike_S1_NTD", "NTD", gff_sub_final_prot$note_name)
gff_sub_final_prot$note_name <- gsub("SARS-CoV-2_Spike_S1_RBD", "RBD", gff_sub_final_prot$note_name)

gff_spike_regions <- gff_sub_final_prot[-1,]

gff_spike_regions <- gff_spike_regions[!grepl("\\[polypeptide binding\\]", gff_spike_regions$note_name),]
gff_spike_regions <- gff_spike_regions[!grepl("\\[posttranslational modification\\]", gff_spike_regions$note_name),]
gff_spike_regions <- gff_spike_regions[!grepl("SARS-CoV-like_Spike_SD1-2_S1-S2_S2", gff_spike_regions$note_name),]



colors <- length(unique(gff_spike_regions$note_name))

colors_test <- viridis(colors, alpha = 0.7)
testtt <- as.data.frame(colors_test)

gff_spike_regions <- cbind(gff_spike_regions, testtt)

gff_spike_regions$sort <- (gff_spike_regions$start - gff_spike_regions$end)


finalDF_prot_S <- finalDF_prot[grepl("S", finalDF_prot$Gene),]



 
## assigning values to each of the spike regions based on notation found in gff file and online ncbi 43740568 
finalDF_prot_S$Position <- as.numeric(finalDF_prot_S$Position)

finalDF_prot_S$region <- ifelse(finalDF_prot_S$Position >= 13 & finalDF_prot_S$Position <= 304, 'NTD',
                              ifelse(finalDF_prot_S$Position >= 319 & finalDF_prot_S$Position <= 437, 'RBD',
                                     ifelse(finalDF_prot_S$Position >= 438 & finalDF_prot_S$Position <= 508, 'receptor binding motif',
                                            ifelse(finalDF_prot_S$Position >= 509 & finalDF_prot_S$Position <= 541, 'RBD',
                                                   ifelse(finalDF_prot_S$Position >= 672 & finalDF_prot_S$Position <= 684, 'S1/S2 cleavage region',
                                                          ifelse(finalDF_prot_S$Position >= 685 & finalDF_prot_S$Position <= 686, 'S1/S2 furin cleavage site',
                                                                 ifelse(finalDF_prot_S$Position >= 687 & finalDF_prot_S$Position <= 709, 'S1/S2 cleavage region',
                                                                        ifelse(finalDF_prot_S$Position >= 788 & finalDF_prot_S$Position <= 806, 'fusion peptide',
                                                                               ifelse(finalDF_prot_S$Position >= 816 & finalDF_prot_S$Position <= 833, 'internal fusion peptide',
                                                                                      ifelse(finalDF_prot_S$Position >= 918 & finalDF_prot_S$Position <= 933, 'heptad repeat 1',
                                                                                             ifelse(finalDF_prot_S$Position >= 1162 & finalDF_prot_S$Position <= 1203, 'heptad repeat 2', 'no region'))))))))))) 
                                                                                                                   




col1 <- c(length(unique(finalDF_prot_S$region)))
# because some of the entries are so much larger, going to subset into 2 plots
subset_finalDF <- subset(finalDF, n < 40)

# define color gradient using viridis package
my_palette <- viridis(col1, alpha = 0.7)
#my_palette_sub <- viridis(max(subset_finalDF$n), alpha = 0.7)

library(cowplot)
# library(RColorBrewer)
# n <- c(length(unique(finalDF_prot_S$region)))
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# x <- sample(col_vector, n)
# pie(rep(1,n), col=sample(col_vector, n))

##test remove later

finalDF <- finalDF_prot_S


counts <- data.frame(table(finalDF$n))

col222 <- length(unique(finalDF$region))

dd.col <- rainbow(col222)
names(dd.col) <- unique(finalDF$region)

plot_prot_notation <- ggplot(gff_spike_regions, aes(y= note_name, color = note_name, alpha = 1, size = 4)) +
  # add overlapping horizontal bars for each note
  geom_segment(stat = "identity", aes(x=start, xend=end, y=reorder(note_name, -sort), yend=note_name), 
               size=5, position = position_dodge(width = 100000), alpha = 1) +
  # add labels for each note
  # add x-axis label
  labs(x=element_blank()) +
  scale_color_manual(values=dd.col) +
  # customize theme
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank(),
        legend.position = "none")

plot(plot_prot_notation)



# create a new column with offset values for points at the same position
finalDF1 <- finalDF %>%
  group_by(Position) %>%
  mutate(offset = seq(-1, 1, length.out = n())) %>%
  ungroup()

finalDF1 <- subset(finalDF1, n < 4)




#### 

# TESTING WITH THE HEATMAP
# need to go back and collect the patient and time information again

testing <- overall
### THIS IS JUST AN EXAMPLE, NEED TO MAKE IT SO THAT THIS IS ALREADY IN THE FASTA HEADER
testing$name <- paste("time10", testing$name, sep = "_")
testing$name <- paste("patientX", testing$name, sep = "_")

#something here is fucky, but im not going to poke into it too much right now
for (i in 1:nrow(testing)) {
  testxx <- testing$name[i]
  lenx <- str_count(i, "_")
  if(lenx < 16){
    testing[c('c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'c10', 'c11', 'c12', 'c13', 'c14', 'c15')] <- str_split_fixed(testing$name, "_", 18)
  }
  else if(lenx > 15){
    testing[c('c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'c10', 'c11', 'c12', 'c13', 'c14', 'c15', 'c16')] <- str_split_fixed(testing$name, "_", 18)
  }
  else {
    print("ERROR: COLUMN NUMBERS ARE OFF")
  }
}

# could test for link mutations, but im going to skip this for now

linked_muations <- testing %>% 
  group_by(name) %>%
  summarize(n = n())

testing <- testing %>% mutate(c2 = str_remove(c2, "time"))
names(testing)[names(testing) == 'c2'] <- 'time'
names(testing)[names(testing) == 'c1'] <- 'patient'
names(testing)[names(testing) == '13'] <- 'DNA_pos'
names(testing)[names(testing) == '14'] <- 'AA_pos'


testing_edit <- select(testing, -c("ID", "QUAL", "FILTER", "FORMAT", "notsure", "5", "6", "7",
                                    "8", "9", "12", "name", 'c3', 'c4', 'c5', 'c6', 'c7', 'c8',
                                   'c9', 'c10', 'c11', 'c12', 'c13', 'c14', 'c15', "DNA"))



testing_edit$Protein <- gsub("p\\.", "", testing_edit$Protein)
  
  # Split the mutation column into three parts using the regular expression
testing_edit[, c("aa_ref", "Position", "aa_mut")] <- str_match(testing_edit$Protein, split_regex)[, -1]


# thought this could work but doesnt look like it works
# for (i in 1:nrow(testing_edit)) {
#   S_check <- any(grepl("S", testing_edit$Gene[i]))
#   if(S_check == TRUE){
testing_edit$region <- ifelse(testing_edit$Gene != "S", "Not Spike",
                              ifelse(testing_edit$Position >= 13 & testing_edit$Position <= 304, 'NTD',
                                     ifelse(testing_edit$Position >= 319 & testing_edit$Position <= 437, 'RBD',
                                            ifelse(testing_edit$Position >= 438 & testing_edit$Position <= 508, 'receptor binding motif',
                                                   ifelse(testing_edit$Position >= 509 & testing_edit$Position <= 541, 'RBD',
                                                          ifelse(testing_edit$Position >= 672 & testing_edit$Position <= 684, 'S1/S2 cleavage region',
                                                                 ifelse(testing_edit$Position >= 685 & testing_edit$Position <= 686, 'S1/S2 furin cleavage site',
                                                                        ifelse(testing_edit$Position >= 687 & testing_edit$Position <= 709, 'S1/S2 cleavage region',
                                                                               ifelse(testing_edit$Position >= 788 & testing_edit$Position <= 806, 'fusion peptide',
                                                                                      ifelse(testing_edit$Position >= 816 & testing_edit$Position <= 833, 'internal fusion peptide',
                                                                                             ifelse(testing_edit$Position >= 918 & testing_edit$Position <= 933, 'heptad repeat 1',
                                                                                                    ifelse(testing_edit$Position >= 1162 & testing_edit$Position <= 1203, 'heptad repeat 2', 'no region'))))))))))))


# } else {
#   testing_edit$region[testing_edit$S_check < 1] <- "Not Spike"
# }

finalDF_heatmap <- testing_edit %>% 
  group_by(Gene, POS, REF, ALT, Effect, Protein, patient, time, aa_ref, Position, aa_mut, region) %>%
  summarize(n = n())

testxx3 <- finalDF_heatmap

testxx3$nmax <- print(n)

testxx3$percentage <- testxx3$n/testxx3$nmax

testxx4 <- testxx3



# TESTING WITH A SIMILAR DF
test_14 <- as.data.frame(read.csv("test_day_14_v3.csv"))
test_14$Nucleotide <- paste(test_14$REF, test_14$POS, test_14$ALT)

testxx4$Nucleotide <- paste(testxx4$REF, testxx4$POS, testxx4$ALT)
testrbind <- testxx4
testrbind <- transform(testrbind, time = as.numeric(time))
testrbind <- transform(testrbind, Position = as.numeric(Position))

#test_14$time <- as.character(test_14$time)
## so this is what the df should look like
test_nuclear <- rbind(testrbind, test_14)

test_nuclear_use <- test_nuclear
test_protein_use <- test_nuclear

list_time <- unique(test_nuclear$time)
list_patient <- unique(test_nuclear$patient)

testy <- test_nuclear


temp324 <- as.data.frame(read.csv("df_test.csv"))
dput(temp324)

temp325 <- as.data.frame(read.csv("final_df_test.csv"))
dput(temp325)



df_merged <- temp324 %>%
  unite("Nucleotide", starts_with("Nucleotide"), sep = "_", na.rm = TRUE) %>%
  mutate(Nucleotide = sub(".*_", "", Nucleotide))

# remove rows with missing data

# reshape the dataframe
df_reshaped <- df_merged %>%
  pivot_wider(names_from = time, values_from = c(n_10, n_14))

# remove row names
rownames_to_column(df_reshaped, "Nucleotide")




df_merged <- unite(temp324, col = "Nucleotide", starts_with("Nucleotide"), na.rm = TRUE, remove = TRUE)
df_complete <- complete(df_merged, Nucleotide, n_10, n_14, fill = list(n_10 = NA, n_14 = NA))



df_reshaped <- temp324 %>% 
  rownames_to_column(var = "Nucleotide") %>% 
  pivot_longer(cols = starts_with("n_"), names_to = "time", values_to = "value") %>% 
  mutate(time = gsub("n_", "", time)) %>% 
  pivot_wider(names_from = "time", values_from = "value")

# Remove the row names and add the "Nucleotide" column back as the row names
rownames(df_reshaped) <- df_reshaped$Nucleotide
df_reshaped$Nucleotide <- NULL

# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Collapse the Nucleotide columns
# temp324 <- temp324 %>%
#   mutate(Nucleotide_combined = coalesce(matches("Nucleotide"), matches("Nucleotide\\.\\d+")))
# 
# # Melt the data
# df_melted <- melt(temp324, id.vars = c("Nucleotide", "patient"), 
#                   measure.vars = c(paste0("n_", c("10", "14"))), 
#                   variable.name = "time", value.name = "value")
# 
# # Create the heatmap using ggplot2
# ggplot(df_melted, aes(x = Nucleotide, y = factor(time), fill = value)) +
#   geom_tile(aes(width = 1, height = 1), color = "white") +
#   scale_fill_gradient(low = "blue", high = "red") +
#   facet_grid(rows = vars(time)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# temp324$Nucleotide_combined <- ifelse(is.na(temp324$Nucleotide), temp324$Nucleotide.1, temp324$Nucleotide)
# 
# # Remove the original Nucleotide columns
# temp324 <- temp324 %>%
#   select(-matches("Nucleotide"))
# 
# # Melt the data
# df_melted <- melt(temp324, id.vars = c("Nucleotide_combined", "patient"), 
#                   measure.vars = c(paste0("n_", c("10", "14"))), 
#                   variable.name = "time", value.name = "value")
# 
# # Create the heatmap using ggplot2
# ggplot(df_melted, aes(x = Nucleotide_combined, y = factor(time), fill = value)) +
#   geom_tile(aes(width = 1, height = 1), color = "white") +
#   scale_fill_gradient(low = "blue", high = "red") +
#   facet_grid(rows = vars(time)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# mutate(Nucleotide_combined = coalesce(Nucleotide.1, Nucleotide)) %>%
#   
#   # select only the relevant columns for melting
#   select(Nucleotide_combined, patient, starts_with("n_")) %>%
#   
#   # reshape the data from wide to long format
#   pivot_longer(cols = starts_with("n_"), 
#                names_to = c("time"), 
#                values_to = "value") %>%
#   
#   # filter out missing values
#   filter(!is.na(value)) %>%
#   
#   # create the heatmap using ggplot2
#   ggplot(aes(x = Nucleotide_combined, y = time, fill = value)) +
#   geom_tile(aes(width = 1, height = 1), color = "white") +
#   scale_fill_gradient(low = "blue", high = "red") +
#   facet_grid(rows = vars(time)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# df_melted <- melt(temp324, id.vars = c("Nucleotide", "patient"), 
#                   measure.vars = c(paste0("n_", c("10", "14"))), 
#                   variable.name = "time", value.name = "value")
# 
# nucleotides <- unique(temp324$Nucleotide)
# 
# # Fill in missing values with 0
# df_melted <- tidyr::complete(df_melted, Nucleotide = nucleotides, time = c(10, 14), fill = list(value = 0))
# 
# # Create the heatmap using ggplot2
# ggplot(df_melted, aes(x = Nucleotide, y = factor(time), fill = value)) +
#   geom_tile(aes(width = 1, height = 1), color = "white") +
#   scale_fill_gradient(low = "blue", high = "red") +
#   facet_grid(rows = vars(time)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# 
# 
# 
# df_melted <- melt(temp324, id.vars = c("Nucleotide", "patient"), 
#                   measure.vars = c(paste0("n_", c("10", "14"))), 
#                   variable.name = "time", value.name = "value")
# 
# # Create the heatmap using ggplot2
# ggplot(df_melted, aes(x = Nucleotide, y = factor(time), fill = value)) +
#   geom_tile(aes(width = 1, height = 1), color = "white") +
#   scale_fill_gradient(low = "blue", high = "red") +
#   facet_grid(rows = vars(time)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ### the following was all test AND IS NOT NECESSARY USING SHINY###
# 
# # test_nuc_1 <- split(test_nuclear, list(test_nuclear$time, test_nuclear$patient))
# # 
# # #this worked when the # rows was the same, but wont work with diff rows
# # # newdf <- do.call("cbind", test_nuc_1)
# # 
# # list_test <- list(test_nuc_1)
# # 
# # pat_list <- unique(test_nuclear$patient)
# # 
# # newdf <- do.call("cbindX", test_nuc_1)
# # #this works, but need to subset before this step
# # 
# # cnum <- ncol(newdf) / 16
# # 
# # x<-split.default(newdf, rep(1:cnum, each = 16))
# # 
# # newdf <- merge(test_nuc_1, by=c("time", "patient"), all = T)
# 
# # merge(dfA, dfB, all.x=TRUE, all.y=TRUE,
# #       by.x=c("ID", "initials"),
# #       by.y=c("ID_mod", "initials"))
# # 
# # str(test_nuc_1)
# 
# # merged_df <- Reduce(function(x, y) merge(x, y, by = c("time", "patient")), test_nuc_1)
# # heres the plan: split it up based on the day and patient, then merge back together
# 
# 
# 
# 
# 
# # just_sub_test <- testxx4 %>%
# #   select("Protein", "percentage", everything())
# # just_sub_test <- just_sub_test[, -c(3:16)]
# # 
# # str(just_sub_test)
# # 
# # 
# 
# 
# # # Pivot data frame to wide format
# # df_wide <- testxx4 %>% pivot_wider(names_from = Protein, values_from = percentage)
# # 
# # # Remove unnecessary columns
# # df_wide <- df_wide[, -c(1:6)]
# # 
# # # Convert row names to integers
# # row.names(df_wide) <- as.integer(row.names(df_wide))
# # 
# # # Sort data frame by row names
# # df_wide <- df_wide[order(row.names(df_wide)), ]
# # 
# # # Convert data frame to matrix
# # mat <- as.matrix(df_wide)
# 
# 
# 
# # Pivot the data frame to the wide format
# # df_wide <- testxx3 %>%
# #   
# # df_wide <- testxx4 %>% 
# #   pivot_wider(testxx4,
# #               id_cols = time,
# #               id_expand = TRUE,
# #               names_from = c("Protein", "Nucleotide", "Gene", "Effect", "patient", "region"),
# #               names_glue = sprintf('{%s}_{%s}_{%s}_{%s}_{%s}_{%s}', "Protein", "Nucleotide", "Gene", "Effect", "patient", "region"),
# #               values_from = percentage)
# 
# # 
# # df_wide <- testxx4 %>%
# #   pivot_wider(names_from = time,
# #               values_from = percentage)
# #   # group_by(Gene, POS, REF, ALT,
# #   #          Effect, Protein, patient, time,
# #   #          aa_ref, Position, aa_mut, region,
# #   #          Nucleotide) %>%
# #   #summarise(percentage = mean(Acc, na.rm = TRUE)) %>% na.omit() %>% 
# # 
# # 
# # # Create the heatmap
# # pheatmap(df_wide, cluster_cols = FALSE, cluster_rows = FALSE)
# 
# # str(testxx4)
# # 
# # testxx4.melt <- testxx4 %>%
# #   melt(id.vars = c("time", "patient", "Effect", "region", "POS", "Position",
# #                    "Gene"), measure.vars = c("Protein", "n"))
# # 
# # 
# # library('maditr')
# # df_wide_prot <- dcast(testxx4, time ~ Protein + n, value.var = "percentage")
# # 
# # # Create the heatmap
# # pheatmap(as.matrix(df_wide[,-1]), cluster_cols = FALSE, cluster_rows = FALSE, border_color = NA)
# 
# ## THIS WILL SEP THE FIRST 3
# 
# ## testing first to generate heatmap
# 
# 
# 
# 
# 
