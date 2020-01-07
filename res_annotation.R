# annotate the resource file
library(data.table) # fread() function
# library(pracma) # strcmp() function
# print(strcmp("abc", "abc"))

rm(list=ls())

# set directories
work_directory = "/Users/iris/Desktop/mMHseq_post_processing/" # local
# work_directory = "" # server
setwd(work_directory)

# MH ID retouched based on match file, can use the original though
ken_resource_dir = paste0(work_directory, "Ken_90MH_dbSNP_717.txt")
db_1000g_dir = paste0(work_directory, "mH_region_1000g/")

output_dir = paste0(work_directory, "Ken_90MH_dbSNP_717_annotated.txt")

ken_resource <- read.table(ken_resource_dir, header = TRUE, stringsAsFactors = FALSE) ### notice stringsAsFactors

# get all mh IDs
require(sqldf)
MHapID <- sqldf('SELECT DISTINCT mhap FROM ken_resource')
length(MHapID$mhap)

db_1000g <- data.frame()
for (i in 1:length(MHapID$mhap)) {
  temp = fread(file = paste0(db_1000g_dir, MHapID$mhap[i], ".vcf"), skip = "#CHROM", select = c(1:8)) # notice the 'skip'
  db_1000g <- rbind(db_1000g, temp)
}
info <- sqldf('SELECT ID, REF, ALT, INFO FROM db_1000g')

##### extract allele freq info for every possible SNP #####
combined_freq <- NULL

for (i in 1:length(info$ID)) {
  if (grepl("rs", info$ID[i])){
    temp <- strsplit(db_1000g$INFO[[i]], split = ";")
    #print(temp[[1]][6])
    eas_af=strsplit(temp[[1]][6], split = "=")[[1]][2] # Eastern Asian
    #print(temp[[1]][7])
    amr_af=strsplit(temp[[1]][7], split = "=")[[1]][2] # American
    #print(temp[[1]][8])
    afr_af=strsplit(temp[[1]][8], split = "=")[[1]][2] # African Amer
    #print(temp[[1]][9])
    eur_af=strsplit(temp[[1]][9], split = "=")[[1]][2] # Europe
    #print(temp[[1]][10])
    sas_af=strsplit(temp[[1]][10], split = "=")[[1]][2] # South Asian
    freq <- data.frame("snpid" = info$ID[i], "REF" = info$REF[i], "ALT" = info$ALT[i], "EAS" = eas_af, "AMR" = amr_af, "AFR"  = afr_af, "EUR" = eur_af, "SAS" = sas_af,
                       stringsAsFactors = FALSE)
    # print(combined_freq)
    combined_freq <- rbind(combined_freq, freq)
  }
  else {
    print(paste0(info$ID[i], " --- different type of id"))
  }
}

write.table(as.data.frame(combined_freq), file = paste0(work_directory, "combined_1000g.txt"), quote= F, sep= "\t", row.names= F)


####### check for REF, ALT 
# for (i in 1:length(ken_resource$snpid)) {
#   if (ken_resource$snpid[i] != '.') {
#     for (j in 1:length((combined_freq$snpid))) {
#       if (ken_resource$snpid[i] == combined_freq$snpid[j]) {
#         if (ken_resource$REF[i] != combined_freq$REF[j]) {
#           print(paste0(ken_resource$snpid[i], " REF: ", ken_resource$REF[i], "; ",
#                        combined_freq$snpid[j], " REF: ", combined_freq$REF[j]))
#         }
#         if (ken_resource$ALT[i] != combined_freq$ALT[j]) {
#           print(paste0(ken_resource$snpid[i], " ALT: ", ken_resource$ALT[i], "; ",
#                        combined_freq$snpid[j], " ALT: ", combined_freq$ALT[j]))
#         }
#       }
#     }
#   }
# }


####### annotate resource file
length(combined_freq$snpid)
length(ken_resource$snpid)
num_novel <- 0
num_1000g <- 0
num_out <- 0
# class(as.character(ken_resource$snpid[i]))
out <- NULL

for (i in 1:length(ken_resource$snpid)) {
  temp <- NULL
  if (ken_resource$snpid[i] == '.'){
    #print(paste0("i = ", i))
    #print(paste0(ken_resource$snpid[i], " is novel"))
    #print(paste0(ken_resource$snpid[i], "NA", "NA", "NA", "NA", "NA"))
    temp <- data.frame("mhap" = ken_resource$mhap[i], "snpid" = ken_resource$snpid[i], "chr" = ken_resource$chr[i], "pos" = ken_resource$pos[i], 
                      "REF" = ken_resource$REF[i], "ALT" = ken_resource$ALT[i], "EAS" = "NA", "AMR" = "NA", "AFR"  = "NA", "EUR" = "NA", "SAS" = "NA",
                      stringsAsFactors = FALSE)
    out <- rbind(out, temp)
    num_novel = num_novel + 1
    next()
  }
  else if (ken_resource$snpid[i] != '.' & !(ken_resource$snpid[i] %in% combined_freq$snpid)) {
    temp <- data.frame("mhap" = ken_resource$mhap[i], "snpid" = ken_resource$snpid[i], "chr" = ken_resource$chr[i], "pos" = ken_resource$pos[i], 
                       "REF" = ken_resource$REF[i], "ALT" = ken_resource$ALT[i], "EAS" = "NA", "AMR" = "NA", "AFR"  = "NA", "EUR" = "NA", "SAS" = "NA",
                       stringsAsFactors = FALSE)
    out <- rbind(out, temp)
    num_out = num_out + 1
    next()
  }
  for (j in 1:length(combined_freq$snpid)) {
    #if(FALSE) {
    if (ken_resource$snpid[i] == combined_freq$snpid[j]){
      #print(paste0("i = ", i, ", j = ", j))
      #print(paste0(ken_resource$snpid[i], "---", combined_freq$snpid[j], ": ", 
      #             combined_freq$EAS[j], " ", combined_freq$AMR[j], " ", combined_freq$AFR[j], " ", combined_freq$EUR[j], " ", combined_freq$SAS[j]))
      temp <- data.frame("mhap" = ken_resource$mhap[i], "snpid" = ken_resource$snpid[i], "chr" = ken_resource$chr[i], "pos" = ken_resource$pos[i], 
                         "REF" = ken_resource$REF[i], "ALT" = ken_resource$ALT[i], "EAS" = combined_freq$EAS[j], "AMR" = combined_freq$AMR[j], "AFR"  = combined_freq$AFR[j], "EUR" = combined_freq$EUR[j], "SAS" = combined_freq$SAS[j],
                         stringsAsFactors = FALSE)
      out <- rbind(out, temp)
      num_1000g = num_1000g + 1
      break()
    }
  }
}

extra <- sqldf('SELECT snpid FROM ken_resource EXCEPT SELECT snpid FROM combined_freq')
print(extra)
# num_out <- length(extra$snpid)-1

# print(num_novel)
# print(num_1000g)
# print(num_out)
sum <- num_1000g + num_novel + num_out
print(sum)
write.table(as.data.frame(out), file = output_dir, quote= F, sep= "\t", row.names= F)




