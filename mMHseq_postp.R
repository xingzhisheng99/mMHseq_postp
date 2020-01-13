# clean up workspace
rm(list=ls())

# set directories
work_directory = "/Users/iris/Desktop/mMHseq_post_processing/" # local
# work_directory = "" # server
# setwd(work_directory)

ken_resource_dir = paste0(work_directory, "Ken_90MH_dbSNP_717_annotated.txt") # MH ID retouched based on match file
#ken_res_annotated_dir = paste0(work_directory, "Ken_90MH_dbSNP_717_annotated.txt")
plotdata_complete_dir = paste0(work_directory, "website_plotdata_7xx.txt")
# by_mh_dir = paste0(work_directory, "plotdata_by_MH/") # separate plotdata by MH

final_output_dir = paste0(work_directory, "additional_grey_bars_annotated.txt")

ken_resource <- read.table(ken_resource_dir, header = TRUE)
plotdata_complete <- read.table(plotdata_complete_dir, header = TRUE)

require(sqldf)
# get all mh IDs and sample IDs
MHapID <- sqldf('SELECT DISTINCT mhap FROM ken_resource')
length(MHapID$mhap)

SampleID <- sqldf('SELECT DISTINCT SampleID FROM plotdata_complete')
length(SampleID$SampleID)

# for (i in 1:length(MHapID$mhap)) {
#   query <- paste0('SELECT * FROM plotdata_complete WHERE MHapID = \'', MHapID$mhap[i], '\'')
#   # print(query)
#   by_mh <- sqldf(query)
#   write.table(as.data.frame(by_mh), file = paste0(by_mh_dir, MHapID$mhap[i], ".txt"), quote=F, sep="\t", row.names=F)
# }
# 
# mh137 <- read.table(paste0(by_mh_dir, "mh08KK-137.txt"), sep = "\t", header = TRUE)
# resource137 <- sqldf('SELECT * FROM ken_resource WHERE mhap = \'mh08KK-137\'')

# class(resource137$pos) # integer

###### get new grey bars (missing due to batch effect - called as variant in a batch but never in the other)
###### with population freq, annotated with the script res_annotation.R ######
###### might be integrated into this script already ######
###### check if 'Ken_90MH_dbSNP_717_annotated.txt' is in the res folder ######

new_grey <- NULL

#for(n in 1:length(MHapID$mhap)){
for(n in 1:1){
  plotdata_by_mh <- data.frame()
  resource_by_mh <- data.frame()
  
  query_pd <- paste0('SELECT * FROM plotdata_complete WHERE MHapID = \'', MHapID$mhap[n], '\'')
  plotdata_by_mh <- sqldf(query_pd)
  #plotdata_by_mh <- read.table(paste0(by_mh_dir, MHapID$mhap[n], ".txt"), sep = "\t", header = TRUE)
  
  query_res <- paste0('SELECT * FROM ken_resource WHERE mhap = \'', MHapID$mhap[n], '\'')
  resource_by_mh <- sqldf(query_res)
  
  for(i in 155:length(SampleID$SampleID)) {
    query <- paste0('SELECT * FROM plotdata_by_mh WHERE SampleID = \'', SampleID$SampleID[i], '\'')
    by_sample <- sqldf(query)
    print(SampleID$SampleID[i])
    
    for (j in 1:length(resource_by_mh$snpid)) {
      print(resource_by_mh$snpid[j])
      #print(resource_by_mh$snpid[j] %in% by_sample$snpid)
      # if is novel, check if pos is occupied --- in case there are other novels in the sample
      if(resource_by_mh$snpid[j] == '.') {
        #print(paste0(resource_by_mh$pos[j]))
        #print(resource_by_mh$pos[j] %in% by_sample$end)
        #print(resource_by_mh$pos[j] %in% by_sample$start)
        # if((!(resource_by_mh$pos[j] %in% by_sample$end)) & (!(resource_by_mh$pos[j] %in% by_sample$start))) {
        if(!(resource_by_mh$pos[j] %in% by_sample$end)) {
          print(paste0(resource_by_mh$pos[j]))
          hap1 <- NULL
          hap2 <- NULL
          
          hap1 <- data.frame("MHapID" = MHapID$mhap[n], "SampleID" = SampleID$SampleID[i], "snpid" = resource_by_mh$snpid[j], 
                             "haplotype" = "hap1", "start" = resource_by_mh$pos[j]-1, "end" = resource_by_mh$pos[j], 
                             "status" = "Purple_placeholder", "EAS" = resource_by_mh$EAS[j], "AMR" = resource_by_mh$AMR[j], 
                             "AFR" = resource_by_mh$AFR[j], "EUR" = resource_by_mh$EUR[j], "SAS" = resource_by_mh$SAS[j],
                             "REF" = resource_by_mh$REF[j], "ALT" = resource_by_mh$ALT[j], "BaseCall" = resource_by_mh$REF[j], 
                             "haplotype.1" = "NA", "ReadDepth" = "NA", "HetRatio" = "NA",
                             "RefForward" = "NA", "RefReverse" = "NA", "AltForward" = "NA", "AltReverse" = "NA", 
                             "gene_or_locus" = "NA", "Chr" = resource_by_mh$chr[j],
                             stringsAsFactors = FALSE)
          hap2 <- data.frame("MHapID" = MHapID$mhap[n], "SampleID" = SampleID$SampleID[i], "snpid" = resource_by_mh$snpid[j], 
                             "haplotype" = "hap2", "start" = resource_by_mh$pos[j]-1, "end" = resource_by_mh$pos[j], 
                             "status" = "Purple_placeholder", "EAS" = resource_by_mh$EAS[j], "AMR" = resource_by_mh$AMR[j], 
                             "AFR" = resource_by_mh$AFR[j], "EUR" = resource_by_mh$EUR[j], "SAS" = resource_by_mh$SAS[j],
                             "REF" = resource_by_mh$REF[j], "ALT" = resource_by_mh$ALT[j], "BaseCall" = resource_by_mh$REF[j], 
                             "haplotype.1" = "NA", "ReadDepth" = "NA", "HetRatio" = "NA",
                             "RefForward" = "NA", "RefReverse" = "NA", "AltForward" = "NA", "AltReverse" = "NA", 
                             "gene_or_locus" = "NA", "Chr" = resource_by_mh$chr[j],
                             stringsAsFactors = FALSE)
          
          new_info <- NULL
          new_info <- rbind(hap1, hap2)
          new_grey <- rbind(new_grey, new_info)
        }
      }
      else if(!(resource_by_mh$snpid[j] %in% by_sample$snpid)){
        #if(!((resource_by_mh$pos[j] %in% by_sample$end) | (resource_by_mh$pos[j] %in% by_sample$start))){ # for verification only
        hap1 <- NULL
        hap2 <- NULL
        
        hap1 <- data.frame("MHapID" = MHapID$mhap[n], "SampleID" = SampleID$SampleID[i], "snpid" = resource_by_mh$snpid[j], 
                           "haplotype" = "hap1", "start" = resource_by_mh$pos[j]-1, "end" = resource_by_mh$pos[j], 
                           "status" = "Ancestral_placeholder", "EAS" = resource_by_mh$EAS[j], "AMR" = resource_by_mh$AMR[j], 
                           "AFR" = resource_by_mh$AFR[j], "EUR" = resource_by_mh$EUR[j], "SAS" = resource_by_mh$SAS[j],
                           "REF" = resource_by_mh$REF[j], "ALT" = resource_by_mh$ALT[j], "BaseCall" = resource_by_mh$REF[j], 
                           "haplotype.1" = "NA", "ReadDepth" = "NA", "HetRatio" = "NA",
                           "RefForward" = "NA", "RefReverse" = "NA", "AltForward" = "NA", "AltReverse" = "NA", 
                           "gene_or_locus" = "NA", "Chr" = resource_by_mh$chr[j],
                           stringsAsFactors = FALSE)
        hap2 <- data.frame("MHapID" = MHapID$mhap[n], "SampleID" = SampleID$SampleID[i], "snpid" = resource_by_mh$snpid[j], 
                           "haplotype" = "hap2", "start" = resource_by_mh$pos[j]-1, "end" = resource_by_mh$pos[j], 
                           "status" = "Ancestral_placeholder", "EAS" = resource_by_mh$EAS[j], "AMR" = resource_by_mh$AMR[j], 
                           "AFR" = resource_by_mh$AFR[j], "EUR" = resource_by_mh$EUR[j], "SAS" = resource_by_mh$SAS[j],
                           "REF" = resource_by_mh$REF[j], "ALT" = resource_by_mh$ALT[j], "BaseCall" = resource_by_mh$REF[j], 
                           "haplotype.1" = "NA", "ReadDepth" = "NA", "HetRatio" = "NA",
                           "RefForward" = "NA", "RefReverse" = "NA", "AltForward" = "NA", "AltReverse" = "NA", 
                           "gene_or_locus" = "NA", "Chr" = resource_by_mh$chr[j],
                           stringsAsFactors = FALSE)
        
        new_info <- NULL
        new_info <- rbind(hap1, hap2)
        new_grey <- rbind(new_grey, new_info)
      }
    }
  }
}

### new entries only
write.table(as.data.frame(new_grey), file = final_output_dir, quote= F, sep= "\t", row.names= F)
### new complete plotdata file

# write.table(as.data.frame(new_grey), file = paste0(work_directory, "website_plotdata.txt"), quote= F, sep= "\t", row.names= F, append = TRUE)
plotdata_complete_new <- NULL
plotdata_complete_new <- rbind(plotdata_complete, new_grey)
write.table(as.data.frame(plotdata_complete_new), file = paste0(work_directory, "website_plotdata_new.txt"), quote= F, sep= "\t", row.names= F)



# ###### use if generating plain text
# new_grey <- data.frame()
# new_grey <- rbind(paste("MHapID", "SampleID", "snpid", "haplotype", "start", "end", "status",
#                         "EAS", "AMR", "AFR", "EUR", "SAS",
#                         "REF", "ALT", "BaseCall", "haplotype.1",
#                         "ReadDepth", "HetRatio",
#                         "RefForward", "RefReverse", "AltForward", "AltReverse", "gene_or_locus",
#                         "Chr", sep = '\t'))


# ###### no population info ######
# for(n in 1:1){
#   plotdata_by_mh <- data.frame()
#   resource_by_mh <- data.frame()
#   
#   query_pd <- paste0('SELECT * FROM plotdata_complete WHERE MHapID = \'', MHapID$mhap[n], '\'')
#   plotdata_by_mh <- sqldf(query_pd)
#   #plotdata_by_mh <- read.table(paste0(by_mh_dir, MHapID$mhap[n], ".txt"), sep = "\t", header = TRUE)
#   
#   query_res <- paste0('SELECT * FROM ken_resource WHERE mhap = \'', MHapID$mhap[n], '\'')
#   resource_by_mh <- sqldf(query_res)
#   
#   for(i in 150:length(SampleID$SampleID)) {
#     query <- paste0('SELECT * FROM plotdata_by_mh WHERE SampleID = \'', SampleID$SampleID[i], '\'')
#     by_sample <- sqldf(query)
#     print(SampleID$SampleID[i])
#     
#     for (j in 1:length(resource_by_mh$snpid)) {
#       print(resource_by_mh$snpid[j])
#       print(resource_by_mh$snpid[j] %in% by_sample$snpid)
#       if(!(resource_by_mh$snpid[j] %in% by_sample$snpid) & !((resource_by_mh$pos[j] %in% by_sample$end) | (resource_by_mh$pos[j] %in% by_sample$start))){
#         #if(!((resource_by_mh$pos[j] %in% by_sample$end) | (resource_by_mh$pos[j] %in% by_sample$start))){ # for verification only
#         hap1 <- paste(MHapID$mhap[n], SampleID$SampleID[i], resource_by_mh$snpid[j], "hap1", 
#                       resource_by_mh$pos[j]-1, resource_by_mh$pos[j], "Ancestral_placeholder",
#                       "NA", "NA", "NA", "NA", "NA",  
#                       resource_by_mh$REF[j], resource_by_mh$ALT[j], resource_by_mh$REF[j],
#                       "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
#                       resource_by_mh$chr[j], sep = "\t")
#         hap2 <- paste(MHapID$mhap[n], SampleID$SampleID[i], resource_by_mh$snpid[j], "hap2", 
#                       resource_by_mh$pos[j]-1, resource_by_mh$pos[j], "Ancestral_placeholder",
#                       "NA", "NA", "NA", "NA", "NA",  
#                       resource_by_mh$REF[j], resource_by_mh$ALT[j], resource_by_mh$REF[j],
#                       "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
#                       resource_by_mh$chr[j], sep = "\t")
#         new_info <- data.frame()
#         new_info <- rbind(hap1, hap2)
#         new_grey <- rbind(new_grey, new_info)
#       }
#     }
#   }
# }

###### no population info ######

# colnames(new_grey) <- c("MHapID", "SampleID", "snpid", "haplotype", "start", "end", "status", 
#                         "EAS", "AMR", "AFR", "EUR", "SAS", 
#                         "REF", "ALT", "BaseCall", "haplotype.1", 
#                         "ReadDepth", "HetRatio", 
#                         "RefForward", "RefReverse", "AltForward", "AltReverse", "gene_or_locus", 
#                         "Chr")
# print(colnames(new_grey))