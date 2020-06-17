library(tidyverse)
library(rentrez)
library(xml2)
library(XML)
library(Biostrings)

library(rbenchmark)
library(compiler)
library(RSQLite)
library(Biobase)
library(RCurl)
library(ggthemes)
library(ggrepel)
library(doParallel)
library(gtools)

unescape_html <- function(str){
  xt <- xml_text(read_html(paste0("<x>", str, "</x>")))
  return(paste("<meta>", xt, "</meta>"))
}

simplifycolumn <- function(x) {
  return(sapply(x, paste, collapse =":"))
}

getSummaryInfo <- function(iteration, webhistory, maxrecords) {
  print(paste("processing", iteration, "iteration"))
  sumInfo <- entrez_summary("assembly",web_history = webhistory, retstart = iteration*maxrecords+1, retmax=maxrecords, always_return_list = T)
  sumtable <- data.frame(matrix(extract_from_esummary(sumInfo, names(sumInfo[[1]])), nrow = length(sumInfo), byrow = T))
  colnames(sumtable) <- names(sumInfo[[1]])
  return(sumtable)
}

processmeta <- function(x) {
  #a <- as_list(read_xml(unescape_html(x)))
  a <- xmlToList(x)
  a <- as.data.frame(t(matrix(unlist(a$Stats), nrow = 3)))
  row.names(a) <- paste(a$V2,a$V3,sep=".")
  a <- select(a,V1)
  return(a)
}

##following conditions to be tested
##https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/
##excluded-from-refseq, suppressed_refseq,  can not be excluded from refseq

selectrepresentative <- function(tid) {
  ##if only one genome for the taxid, use it as representative
  if (nrow(virusinfo %>% filter(taxid == tid)) == 1) {
    return(virusinfo %>% filter(taxid == tid) %>% .$uid)
  }
  ##select fromtype == "ICTV species exemplar" first
  else if (nrow(virusinfo %>% filter(taxid == tid & fromtype == "ICTV species exemplar")) >= 1) {
    return((virusinfo %>% filter(taxid == tid & fromtype == "ICTV species exemplar") %>% arrange(desc(asmreleasedate_genbank)) %>% .$uid)[1])
  } 
  ## if there is refseq entry for the genome and it is identical to genbank, then use one of the latest genbank release as the representative
  else if (nrow(virusinfo %>% filter(taxid == tid & synonym.similarity == "identical")) >= 1) {
    return((virusinfo %>% filter(taxid == tid & synonym.similarity == "identical") %>% arrange(desc(asmreleasedate_genbank)) %>% .$uid)[1])
  }
  ## if there is no refseq entry for the genome, then use one of the latest genbank release as the representative
  else if (nrow(virusinfo %>% filter(taxid == tid & synonym.similarity == "identical")) == 0) {
    return((virusinfo %>% filter(taxid == tid) %>% arrange(desc(asmreleasedate_genbank)) %>% .$uid)[1])
  } 
  ## flag errors if something is amiss
  else {
    return("something is wrong")
  }
}

getlocalpaths <- function(x){
  #ncbidatadir <- "/media/dstore/covid19/ncbi-genomes-2020-05-07/"
  localfile <- paste(ncbidatadir, str_match(x, "\\S+\\/(\\S+)")[,2], "_genomic.fna.gz", sep = "")
  genome <- ifelse(file.exists(localfile), localfile, "unavailable")
  localfile <- paste(ncbidatadir, str_match(x, "\\S+\\/(\\S+)")[,2], "_genomic.gff.gz", sep = "")
  gff <- ifelse(file.exists(localfile), localfile, "unavailable")
  localfile <- paste(ncbidatadir, str_match(x, "\\S+\\/(\\S+)")[,2], "_cds_from_genomic.fna.gz", sep = "")
  cds <- ifelse(file.exists(localfile), localfile, "unavailable")
  localfile <- paste(ncbidatadir, str_match(x, "\\S+\\/(\\S+)")[,2], "_protein.faa.gz", sep = "")
  proteins <- ifelse(file.exists(localfile), localfile, "unavailable")
  return(c(genome,gff,cds,proteins))  
}

#workingdir <- "/media/dstore/covid19/"
#genbankdir <- "/media/dstore/covid19/ncbi-genbank/"
#taxonomydir <- "/media/dstore/covid19/ncbi-taxonomy/"

##directory management
ncbidatadir <- "/media/dstore/covid19/ncbidatadir"
workingdir <- "/media/dstore/covid19/viralgenomics"

dir.create(ncbidatadir, showWarnings = F, recursive = T)
dir.create(workingdir, showWarnings = F, recursive = T)

setwd(workingdir)

##load data objects
if (file.exists("esearchresults.RData")) {
  load("esearchresults.RData")
} else {
  esearchresults <- NULL
}

##get a list of viral genome assembly from NCBI
viralgenomes <- entrez_search("assembly", 'txid10239[Organism:exp] & "latest genbank"[filter] & "complete genome"[filter] NOT anomalous[filter] NOT "derived from surveillance project"[filter]', use_history = T)
currentuids <- entrez_fetch("assembly", web_history = viralgenomes$web_history, rettype="uilist", retmode = "text")
currentuids <- unlist(strsplit(currentuids,"\\n"))

##if esearchresults is null
if(is.null(esearchresults)){
  maxrecords <- 500
  iterations <- as.integer(viralgenomes$count/maxrecords)
  for (i in 0:iterations) {
    tmp <- getSummaryInfo(iteration = i, webhistory = viralgenomes$web_history, maxrecords=maxrecords)
    esearchresults <- bind_rows(esearchresults, tmp)
  }
}

##check and process for new and depricated entries
new <- currentuids[! currentuids %in% esearchresults$uid]
depricated <- esearchresults$uid[! esearchresults$uid %in% currentuids]

if (length(depricated) > 0) {
  esearchresults <- esearchresults[! esearchresults$uid %in% depricated,]
}

if (length(new) > 0){
  maxrecords <- 300 ##this limit is tested previously, may need to change it in the future if something changes
  for (i in 0:as.integer(length(new)/maxrecords)){
    start <- i*maxrecords+1
    end <- ifelse(i*maxrecords+maxrecords > length(new), length(new), i*maxrecords+maxrecords)
    sumInfo <- entrez_summary("assembly",id = new[start:end], always_return_list = T)
    sumtable <- data.frame(matrix(extract_from_esummary(sumInfo, names(sumInfo[[1]])), nrow = length(sumInfo), byrow = T))
    colnames(sumtable) <- names(sumInfo[[1]])
    esearchresults <- bind_rows(esearchresults,sumtable)
  }
}

save(esearchresults, file="esearchresults.RData")

##load data objects
if (file.exists("virusinfo.RData")) {
  load("virusinfo.RData")
} else {
  virusinfo <- NULL
}

##process esearch results to attach to virusinfo
if (is.null(virusinfo)){
  y <- esearchresults
} else {
  y <- esearchresults[esearchresults$uid %in% new,]
}

##property list is not ordered so just paste it together
y$propertylist <- sapply(y$propertylist, paste, collapse =":")
##escape for html characters in metadata column
y$meta <- sapply(y$meta, unescape_html)
##process multi item columns
y <- y %>% unnest_wider(col = synonym, names_sep = ".", names_repair = "universal")
y <- y %>% unnest_wider(col = anomalouslist, names_sep = ".", names_repair = "universal")
y <- y %>% unnest_wider(col = biosource, names_sep = ".", names_repair = "universal")
y <- y %>% unnest_wider(col = gb_bioprojects, names_sep = ".", names_repair = "universal")
y <- y %>% unnest_wider(col = rs_bioprojects, names_sep = ".", names_repair = "universal")
if ("biosource.infraspecieslist" %in% colnames(y)) {
  y <- y %>% unnest_wider(col = biosource.infraspecieslist, names_sep = ".", names_repair = "universal")
}
##unlist column items lists
y <- data.frame(apply(y,2, simplifycolumn))
##remove sortorder column
y <- select(y, -sortorder)
##-biosource.infraspecieslist.)
##do type convert for columns
y <- type.convert(y)
##coverage can be character string
y$coverage <- as.character(y$coverage)
##date columns are not formatted correctly using type convert so take care of them
y <- mutate_at(y, colnames(y)[grep("date", colnames(y))], as.Date)
##get filebase for genomes
y$filebase <- str_match(y$ftppath_genbank, "\\S+\\/(\\S+)")[,2]


##may need to remove this as no information is being gathered from here
##process meta column to extract stats information out of it
#x <- do.call("cbind", lapply(y$meta[1:100], processmeta))
#x <- t(x)
#x <- type.convert(x)
#y <- cbind(y, x)


##when dates are written in sqlite table, they are written as real numbers
##to convert them back to legible date format
##as.Date(18354.0, "1970-01-01") ##origin must be specified
##Sys.Date() as todays date
##taxid should be used instead of speciestaxid. 
##(1) There are more unique entries in taxid vs speciestaxid
##(2) TaxID seems to refer to a strain level identification vs speciestaxid as name suggests
##There can be more than one assembly per taxid, probably because of different strains

##get Lineage information
setwd(ncbidatadir)
download.file("https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz", "new_taxdump.tar.gz")
system("tar -xvzf new_taxdump.tar.gz")
##lineage information
system('sed "s/\t//g" rankedlineage.dmp | sed "s/\'//g" >rankedlineage.dmp.tmp')
rankedlineage <- read.table("rankedlineage.dmp.tmp", sep="|", quote = "", fill = T)
colnames(rankedlineage) <- c("taxid", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom", "empty")
##host information
system('sed "s/\t//g" host.dmp >host.dmp.tmp')
host <- read.table("host.dmp.tmp", sep="|", quote ="")
colnames(host) <- c("taxid", "potentialhosts")

##merged lineage information
y <- merge(y, rankedlineage, by.x = "taxid", by.y = "taxid", all.x = T)
y <- merge(y, host, by.x = "taxid", by.y = "taxid", all.x = T)
y <- droplevels(y)
y <- select(y, -empty)
rm(rankedlineage,host)
##consolidate host information
y <- y %>% mutate(potentialhosts = case_when( .$potentialhosts == "vertebrates,human" ~ "human,vertebrates", is.na(.$potentialhosts) ~ "Unknown", .$potentialhosts == "vertebrates,invertebrates,human" ~ "human,invertebrates,vertebrates", .$potentialhosts == "invertebrates,vertebrates" ~ "vertebrates,invertebrates", TRUE ~ as.character(.$potentialhosts))) %>% mutate(potentialhosts = as.factor(potentialhosts))

##merge new table with the old information and remove uid that have depricated
virusinfo <- bind_row(virusinfo, y) %>% filter(uid %in% currentuids)

##identify representative genome for a given tax id
spectab <- data.frame(table(virusinfo$taxid))
colnames(spectab) <- c("taxid", "genomes")
spectab$taxid <- as.integer(as.character(spectab$taxid))
spectab$representative <- sapply(spectab$taxid, selectrepresentative)
virusinfo <- virusinfo %>% mutate(isrep = case_when( .$uid %in% spectab$representative ~ "taxrep", TRUE ~ "nontaxrep"))
save(virusinfo, file=paste(workingdir,"/virusinfo.RData", sep = ""))

##get assemblies for new virus entries in the assembly database
tarfile <- "genome_assemblies_genome_fasta.tar"
extension <- "_genomic.fna.gz"
filespresent <- gsub("\\S+\\/", "", untar(tarfile, list=T))
expectedfiles <- paste(y$filebase,extension,sep="")
toadd <- gsub(extension,"",expectedfiles[! expectedfiles %in% filespresent])
if (length(toadd) > 0) {
  toadd <- y %>% filter(filebase %in% toadd) %>% select(ftppath_genbank,filebase)
  toadd <- paste(paste(toadd$ftppath_genbank,toadd$filebase,sep="/"), extension, sep="")
  for (i in toadd) { download.file(i,str_match(i, "\\S+\\/(\\S+)")[,2], quiet = T) }
  toadd <- str_match(toadd, "\\S+\\/(\\S+)")[,2]
  for (i in toadd) { system(paste("tar --append --file=", tarfile, " ", i, sep = "")) }
  file.remove(toadd)
}

getgenomeseqproperties <- function(x) {
  filepath <- x
  filebase <- gsub("_genomic.fna.gz","",x)
  untar("genome_assemblies_genome_fasta.tar",files=c(x))
  info<-NULL
  s <- readDNAStringSet(filepath)
  l <- data.frame(letterFrequency(s,c("A","C","G","T")))
  info <- data.frame(filebase = filebase, seqid = as.character(names(s)), seqlen = as.integer(width(s)))
  info <- cbind(info, l)
  kfreq <- oligonucleotideFrequency(s,9) ##9 mers obtained
  write.table(apply(kfreq,1,paste,collapse=","), file = kcountseqfile, col.names = F, sep = "\t", row.names = names(s), append = T, quote = F)
  write.table(paste(colSums(kfreq), collapse = ","), file = kcountgenomefile, col.names = F, sep = "\t", row.names = filebase, append = T, quote = F)
  writeXStringSet(s, filepath = mergedfasta, append = T, compress = T, width = 70)
  file.remove(x)
  return(info)
}

kcountseqfile <- paste(ncbidatadir,"/virus.sequence.kmercounts.txt", sep = "")
kcountgenomefile <- paste(ncbidatadir,"/virus.genome.kmercounts.txt", sep = "")
mergedfasta <- paste(ncbidatadir,"/virusgenomesmerged.fa.gz", sep = "")
sequenceinfofile <- paste(workingdir,"/sequenceinfo.RData",sep="")

##load data objects
if (file.exists(sequenceinfofile)) {
  load(sequenceinfofile)
} else {
  sequenceinfo <- NULL
}

filespresent <- gsub("\\S+\\/", "", untar("genome_assemblies_genome_fasta.tar", list=T))
files2process <- filespresent[! filespresent %in% paste(sequenceinfo$filebase,extension,sep="")]

#for (i in 1:files2process/500) {
#  start <- i*maxrecords+1
#  end <- ifelse(start + maxrecords - 1 <= nrow(y),  start + maxrecords - 1, nrow(y))
#  print(paste("processing", i, "iteration", start, end))
#  kinfo <- do.call("rbind",lapply(y$genome[start:end], getkmers))
#  print("Done counting kmers")
#  write.table(kinfo, "kmercounts.csv", col.names = F, quote = F, sep = ",", row.names = y$genome[start:end], append = T)
#  print("Done writing to file")
#}

##non-essential columns
##coverage
##biosource.sex
##biosource.isolate

sequenceinfo <- bind_rows(lapply(files2process,getgenomeseqproperties))
sequenceinfo <- sequenceinfo %>% mutate(gbid=str_extract(seqid, "\\S+"))
save(sequenceinfo, file=sequenceinfofile)

head(bind_rows(virusinfo,y))

##following code will generate the list of kmers at a given width
##this is the order of kmer counts to be used later
#width <- 9
#s <- DNAString("NNNNNNNNNNNNNN")
#kfreq <- oligonucleotideFrequency(s, width)





#setwd(workingdir)
#save(virusinfo, file="virusinfo.RData")


#toremove <- filespresent[! filespresent %in% expectedfiles]
#if (length(toremove)>0){
#  for (i in toremove) { system(paste("tar --delete --file=", tarfile, " ", i, sep = ""))}
#}

##for protein files
tarfile <- "genome_assemblies_prot_fasta.tar"
extension <- "_protein.faa.gz"
filespresent <- gsub("\\S+\\/", "", untar(tarfile, list=T))
expectedfiles <- paste(y$filebase,extension,sep="")
toadd <- gsub(extension,"",expectedfiles[! expectedfiles %in% filespresent])
toadd <- y %>% filter(filebase %in% toadd & grepl("genbank_has_annotation", propertylist)) %>% .$filebase
if (length(toadd) > 0) {
  toadd <- y %>% filter(filebase %in% toadd) %>% select(ftppath_genbank,filebase)
  toadd <- paste(paste(toadd$ftppath_genbank,toadd$filebase,sep="/"), extension, sep="")
  for (i in toadd) { 
    tryCatch(download.file(i,str_match(i, "\\S+\\/(\\S+)")[,2], quiet = T), 
    error = function(e) print(paste(i, 'did not work out'))) 
  }
  toadd <- str_match(toadd, "\\S+\\/(\\S+)")[,2]
  for (i in toadd) { system(paste("tar --append --file=", tarfile, " ", i, sep = "")) }
  file.remove(toadd)
}
toremove <- filespresent[! filespresent %in% expectedfiles]
if (length(toremove)>0){
  for (i in toremove) { system(paste("tar --delete --file=", tarfile, " ", i, sep = ""))}
}

##for CDS files
tarfile <- "genome_assemblies_cds_fasta.tar"
extension <- "_cds_from_genomic.fna.gz"
filespresent <- gsub("\\S+\\/", "", untar(tarfile, list=T))
expectedfiles <- paste(y$filebase,extension,sep="")
toadd <- gsub(extension,"",expectedfiles[! expectedfiles %in% filespresent])
toadd <- y %>% filter(filebase %in% toadd & grepl("genbank_has_annotation", propertylist)) %>% .$filebase
if (length(toadd) > 0) {
  toadd <- y %>% filter(filebase %in% toadd) %>% select(ftppath_genbank,filebase)
  toadd <- paste(paste(toadd$ftppath_genbank,toadd$filebase,sep="/"), extension, sep="")
  for (i in toadd) { 
    tryCatch(download.file(i,str_match(i, "\\S+\\/(\\S+)")[,2], quiet = T), 
             error = function(e) print(paste(i, 'did not work out'))) 
  }
  toadd <- str_match(toadd, "\\S+\\/(\\S+)")[,2]
  for (i in toadd) { system(paste("tar --append --file=", tarfile, " ", i, sep = "")) }
  file.remove(toadd)
}
toremove <- filespresent[! filespresent %in% expectedfiles]
if (length(toremove)>0){
  for (i in toremove) { system(paste("tar --delete --file=", tarfile, " ", i, sep = ""))}
}

##for gff files
tarfile <- "genome_assemblies_genome_gff.tar"
extension <- "_genomic.gff.gz"
filespresent <- gsub("\\S+\\/", "", untar(tarfile, list=T))
expectedfiles <- paste(y$filebase,extension,sep="")
toadd <- gsub(extension,"",expectedfiles[! expectedfiles %in% filespresent])
toadd <- y %>% filter(filebase %in% toadd & grepl("genbank_has_annotation", propertylist)) %>% .$filebase
if (length(toadd) > 0) {
  toadd <- y %>% filter(filebase %in% toadd) %>% select(ftppath_genbank,filebase)
  toadd <- paste(paste(toadd$ftppath_genbank,toadd$filebase,sep="/"), extension, sep="")
  for (i in toadd) { 
    tryCatch(download.file(i,str_match(i, "\\S+\\/(\\S+)")[,2], quiet = T), 
             error = function(e) print(paste(i, 'did not work out'))) 
  }
  toadd <- str_match(toadd, "\\S+\\/(\\S+)")[,2]
  for (i in toadd) { system(paste("tar --append --file=", tarfile, " ", i, sep = "")) }
  file.remove(toadd)
}
toremove <- filespresent[! filespresent %in% expectedfiles]
if (length(toremove)>0){
  for (i in toremove) { system(paste("tar --delete --file=", tarfile, " ", i, sep = ""))}
}

genomefiles <- untar("genome_assemblies_genome_fasta.tar", list=T)



#x <- bind_cols(lapply(y$ftppath_genbank, getlocalpaths))
#x <- t(x)
#colnames(x) <- c("genome","gff","cds","proteins")
#x <- data.frame(x) %>% mutate_all(as.character)
#y <- cbind(y, x)

##
getgenomeseqproperties <- function(x) {
  filepath <- x
  filebase <- gsub("_genomic.fna.gz","",x)
  untar("genome_assemblies_genome_fasta.tar",files=c(x))
  info<-NULL
  s <- readDNAStringSet(filepath)
  l <- data.frame(letterFrequency(s,c("A","C","G","T")))
  info <- data.frame(filebase = filebase, seqid = as.character(names(s)), seqlen = as.integer(width(s)))
  info <- cbind(info, l)
  file.remove(x)
  return(info)
}

seqinfo <- bind_rows(lapply(genomefiles,getgenomeseqproperties))


getkmers <- function(x) {
  width<-9
  s <- readDNAStringSet(x)
  kfreq <- oligonucleotideFrequency(s,width)
  return(colSums(kfreq))
}



y <- y %>% filter(genome!="unavailable")


kinfo <- do.call("rbind",lapply(y$genome[1:maxrecords], getkmers))
write.table(kinfo, "kmercounts.csv", col.names = T, quote = F, sep = ",", row.names = y$genome[1:maxrecords])
for (i in 1:iterations) {
  start <- i*maxrecords+1
  end <- ifelse(start + maxrecords - 1 <= nrow(y),  start + maxrecords - 1, nrow(y))
  print(paste("processing", i, "iteration", start, end))
  kinfo <- do.call("rbind",lapply(y$genome[start:end], getkmers))
  print("Done counting kmers")
  write.table(kinfo, "kmercounts.csv", col.names = F, quote = F, sep = ",", row.names = y$genome[start:end], append = T)
  print("Done writing to file")
}

kmer <- read.big.matrix("kmercounts.subset.csv", type = "short", sep = ",", binarydescriptor = T, has.row.names = T, header = T, backingfile = "kmer.bigmatrix", descriptorfile = "kmer.bigmatrix.desc")
#k <- kmer
kmer <- kmer[1:50,]
kdist <- filebacked.big.matrix(nrow(kmer),nrow(kmer), type = "double", backingfile = "kmerdist.bk", descriptorfile = "kmerdist.desc")

#getdistance <- function(i,j){
#  shared <- kmer[i,]>0 & kmer[j,]>0
#  kdist[i,j] <- kdist[j,i] <- sum(kmer[i,][shared], kmer[j,][shared]) / (sum(kmer[i,]) + sum(kmer[j,]))
#}

for (i in 1:nrow(kmer)) {
  for (j in i:nrow(kmer)) {
    shared <- kmer[i,]>0 & kmer[j,]>0
    kdist[i,j] <- kdist[j,i] <- sum(kmer[i,][shared], kmer[j,][shared]) / (sum(kmer[i,]) + sum(kmer[j,]))
  }
}

lapply(c(1:10), testf, y=c(1:10))

testf <- function(i,j){
  shared <- kmer[i,]>0 & kmer[j,]>0
  kdist[i,j] <- kdist[j,i] <- sum(kmer[i,][shared], kmer[j,][shared]) / (sum(kmer[i,]) + sum(kmer[j,]))
}
checklocalfiles <- function(x) {
  propertylist <- x["propertylist"]
  ftppath <- x["ftppath_genbank"]
  dcmds <- vector()
  #genome file
  genomedir <- "/media/dstore/covid19/ncbi.genbank.genomes/"
  localfile <- paste(genomedir, str_match(ftppath, "\\S+\\/(\\S+)")[,2], "_genomic.fna.gz", sep = "")
  remotefile <- paste(ftppath, "/", str_match(ftppath, "\\S+\\/(\\S+)")[,2], "_genomic.fna.gz" , sep = "")
  if (! file.exists(localfile) || file.info(localfile)$size == 0){
    dcmds <- append(dcmds, paste("wget -q --continue -O", localfile, remotefile))
  }
  #other annotations
  if (grepl("genbank_has_annotation", propertylist)){
    ##cds file
    cdsdir <- "/media/dstore/covid19/ncbi.genbank.cds/"
    localfile <- paste(cdsdir, str_match(ftppath, "\\S+\\/(\\S+)")[,2], "_cds_from_genomic.fna.gz", sep = "")
    remotefile <- paste(ftppath, "/", str_match(ftppath, "\\S+\\/(\\S+)")[,2], "_cds_from_genomic.fna.gz" , sep = "")
    if (! file.exists(localfile) || file.info(localfile)$size == 0){
      dcmds <- append(dcmds, paste("wget -q --continue -O", localfile, remotefile))
    }
    ## protein file
    proteindir <- "/media/dstore/covid19/ncbi.genbank.proteins/"
    localfile <- paste(proteindir, str_match(ftppath, "\\S+\\/(\\S+)")[,2], "_protein.faa.gz", sep = "")
    remotefile <- paste(ftppath, "/", str_match(ftppath, "\\S+\\/(\\S+)")[,2], "_protein.faa.gz" , sep = "")
    if (! file.exists(localfile) || file.info(localfile)$size == 0){
      dcmds <- append(dcmds, paste("wget -q --continue -O", localfile, remotefile))
    }
    ## gff file
    gffdir <- "/media/dstore/covid19/ncbi.genbank.gff/"
    localfile <- paste(gffdir, str_match(ftppath, "\\S+\\/(\\S+)")[,2], "_genomic.gff.gz", sep = "")
    remotefile <- paste(ftppath, "/", str_match(ftppath, "\\S+\\/(\\S+)")[,2], "_genomic.gff.gz" , sep = "")
    if (! file.exists(localfile) || file.info(localfile)$size == 0){
      dcmds <- append(dcmds, paste("wget -q --continue -O", localfile, remotefile))
    }
  }
  return(dcmds)
}

x <- apply(y, 1, checklocalfiles)
x <- unlist(x)
if ( ! is_empty(x)) {
  fileConn<-file("/media/dstore/covid19/genbankgenomedownloadcmds.txt")
  writeLines(x, fileConn)
  close(fileConn)
}
##following is not running. will need to figure it out.
system("parallel -a /media/dstore/covid19/genbankgenomedownloadcmds.txt -j 8")


##read in genome fasta files, seqids, lengths, gc contents, number of sequences
ginfo <- do.call("rbind", apply(y, 1, processfasta))
y <- merge(y,ginfo %>% group_by(uid) %>% summarise(gc = sum(G,C)/sum(A,C,G,T), nseq = n()) %>% select(uid,gc,nseq), by.x = "uid", by.y = "uid")

getkmers <- function(x) {
  width<-9
  s <- readDNAStringSet(x)
  kfreq <- oligonucleotideFrequency(s,width)
  return(kfreq)
}

kmerprobes <- permutations(4,9,c("A","C","G","T"), repeats.allowed = T)
kmerprobes <- apply(kmerprobes, 1, paste, collapse="")
kmerprobes <- DNAStringSet(kmerprobes)





#boxplots to show GC contents for different groupings (family, order, class, potential hosts)
#y %>% filter(isrep=="taxrep") %>% ggplot(aes(x=family, y=gc)) + geom_boxplot() + geom_jitter(alpha=0.2) + theme_bw() + theme(axis.text = element_text(angle=90,hjust=1))


#cores2use <- detectCores() / 2
#cluster <- makeCluster(cores2use, type="FORK")
#registerDoParallel(cluster)
#x <- parLapply(cluster, y$ftppath_genbank, downloadgb)
#stopCluster(cluster)

##create sqlite database for later use
db = dbConnect(SQLite(), dbname="/media/dstore/covid19/ncbi.virus.20200409.sqlite")
dbWriteTable(db, "assembly", z)
