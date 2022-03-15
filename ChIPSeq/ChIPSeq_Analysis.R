#=============================================================================#
# ChiP-Seq Analysis Workflow:                                                 #
#=============================================================================#

#=============================================================================#
# The following performs a ChIP-Seq Analysis on three different Human         #
# receptors. However this script can be easily modified for any ChIP-Seq data #
#=============================================================================#

# Import Libraries:

# if not already installed, uncomment and run these commands
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocStyle")
#biocLite("systempipeR")
#biocLite("ShortRead")
#biocLite("ChIPseeker")
#biocLite("biomaRt")
#biocLite("Biostrings")
#biocLite("seqLogo")
#biocLite("BCRANK")

# note this R script requires two other programs to run
# namely, Bowtie2 and MACS2.
# it also relies upon output from HOMER.

library(systemPipeR)
library(systemPipeRdata)
library(ShortRead)
library(ChIPseeker)
library(biomaRt)
library(Biostrings)
library(seqLogo)
library(BCRANK)

#=============================================================================#
# Prepare the systemPipeR analysis pipeline:

genWorkenvir("chipseq") #to use the directory structure in this pipeline

# set the working directory as the location of the FASTQ files and workflow
system2("cd", args = "/home/will/R/x86_64-pc-linux-gnu-library/3.2/systemPipeRdata/extdata/workflows/chipseq")
system2("make", args = "-B")
setwd("/home/will/R/x86_64-pc-linux-gnu-library/3.2/systemPipeRdata/extdata/workflows/chipseq")

targetspath = system.file("extdata", "chip_targets.txt", package = "systemPipeR")
targets  = read.delim(targetspath, comment.char = "#")

# note on the systemPipeR .param files:
# these are the file required by the systemPipeR workflow in order for many of 
# the command line arguments to run
# some of the path names in these files may need to be altered in order for 
# the bash software to be properly located and run
 
#=============================================================================#
# Perform quality control of reads:

# uncomment and perform any trimming if necessary
#args = systemArgs(sysma="../../param/trim.param", 
#mytargets = "chip_targets.txt")

#filterReads = function(fq, cutoff=20, Nexceptions=0){
#  qcount = rowSums(as(quality(fq), "matrix") <= cutoff)
#  fq[qcount <=Nexceptions]}
#preprocessReads(args=args, Fct = "filterReads(fq, cutoff=20, Nexceptions=0)",
#                batchsize = 100000,  overwrite=TRUE)
#writeTargetsout(x=args, file="targets_chip_trim.txt", overwrite = TRUE)

# note: if the above code was uncommented and run, change the following
# argument: mytargets = "targets_chip_trim.txt"

dev.off() #flush out plots
# create a FASTQ quality report
args = systemArgs(sysma="../../param/bowtieSE.param", 
                  mytargets = "chip_targets.txt")
fqlist = seeFastq(fastq = infile1(args), batchsize = 100000, klength = 8)
pdf("./results/fastqReport.pdf", height = 18, width=4*length(fqlist))

# plot the FASTQ quality information for each read sequentially
for (i in seq(1, length(fqlist), by=1)){
  seeFastqPlot(fqlist[i])}
dev.off()#close the graphics to the pdf file

#=============================================================================#
# Align the fastq files to a reference with Bowtie2:

args = systemArgs(sysma="../../param/bowtieSE.param",
                   mytargets="chip_targets.txt")
sysargs(args)[1] # check command-line parameters for first FASTQ file
#moduleload(modules(args)) # uncomment and run if a module system is used

# if required, uncomment to unzip the indexed genomes
# change file name to the relevant reference indexed genome
# note if the genome is not already indexed then it needs to
# be indexed using the "bowtie2-build" bash command
#system2("unzip", args = "hg19.1.zip")
#system2("unzip", args = "hg19.2.zip")
#system2("unzip", args = "hg19.3.zip")

# change the paths below to the reference genome machine location
# alignment can take some time to run, but only need be done once
runCommandline(args)
writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE)

# check if all BAM files are created
file.exists(outpaths(args))

# create read and alignment statistics
align_info = alignStats(args=args)
write.table(align_info, "results/alignStats.xls", row.names = FALSE,
            quote = FALSE, sep = "\t")

#=============================================================================#
# Perform peak calling with MACS2:

# before peak calling merge the BAM files by replicates
args = systemArgs(sysma = NULL, mytargets = "targets_bam.txt")
arg_merge = mergeBamByFactor(args, overwrite = TRUE)
writeTargetsout(x = arg_merge, file = "targets_mergeBamByFactor.txt", 
                overwrite = TRUE)

# if there are no controls, call peaks without a reference sample
args = systemArgs(sysma = "../../param/macs2_noinput.param", 
                  mytargets = "targets_mergeBamByFactor.txt")
sysargs(args)[1] # check command-line parameters for first BAM file
runCommandline(args)
file.exists(outpaths(args))
writeTargetsout(x = args, file = "targets_macs.txt", overwrite = TRUE)

# else, call peaks with a reference sample input
writeTargetsRef(infile = "targets_mergeBamByFactor.txt", 
                outfile = "targets_bam_ref.txt", silent = FALSE, overwrite = T)
args = systemArgs(sysma = "../../param/macs2.param", 
                  mytargets = "targets_bam_ref.txt")
sysargs(args)[1] # check command-line parameters for first BAM file
runCommandline(args)
file.exists(outpaths(args))

#=============================================================================#
# Perform peak annotation with ChIPseeker:

args = systemArgs(sysma="../../param/annotate_peaks.param", 
                  mytargets="targets_macs.txt")
txdb = loadDb("hg19.sqlite") #load the queryable reference genome
# annotate peaks

dev.off() #flush out plots
pdf("./results/upsetplots.pdf", height = 6, width=8)

for(i in seq(along=args)){
    peakAnno = annotatePeak(infile1(args)[i], TxDb=txdb, verbose=FALSE)
    upsetplot(peakAnno, vennpie=TRUE)
    df = as.data.frame(peakAnno)
    write.table(df, outpaths(args[i]), quote=FALSE, row.names=FALSE, sep="\t")
}
dev.off() #flush out plots

writeTargetsout(x=args, file="targets_peakanno.txt", overwrite=TRUE)

# save some summary annotation plots to a file
# there should be three plots for each peak file.
pdf("./results/peakAnnoReport.pdf", height = 4, width=20)

peak = readPeakFile(infile1(args)[1])
covplot(peak, weightCol="X.log10.pvalue.")
peakHeatmap(outpaths(args)[1], TxDb=txdb, upstream=3000, 
              downstream=3000, color="blue")
plotAvgProf2(outpaths(args)[1], TxDb=txdb, upstream=3000, downstream=3000,
               xlab="Genomic Region (5'->3')",ylab = "Read Count Frequency")
peak = readPeakFile(infile1(args)[2])
covplot(peak, weightCol="X.log10.pvalue.")
peakHeatmap(outpaths(args)[2], TxDb=txdb, upstream=3000, 
            downstream=3000, color="blue")
plotAvgProf2(outpaths(args)[2], TxDb=txdb, upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')",ylab = "Read Count Frequency")

peak = readPeakFile(infile1(args)[3])
covplot(peak, weightCol="X.log10.pvalue.")
peakHeatmap(outpaths(args)[3], TxDb=txdb, upstream=3000, 
            downstream=3000, color="blue")
plotAvgProf2(outpaths(args)[3], TxDb=txdb, upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')",ylab = "Read Count Frequency")

dev.off() #close the graphics to the pdf file

#=============================================================================#
# Count reads overlapping peak regions:

args = systemArgs(sysma="../../param/count_rangesets.param",
                  mytargets="targets_macs.txt")
args_bam = systemArgs(sysma=NULL, mytargets="targets_bam.txt")
bfl = BamFileList(outpaths(args_bam), yieldSize=50000, index=character())
countDFnames = countRangeset(bfl, args, mode="Union", ignore.strand=TRUE)
writeTargetsout(x=args, file="targets_countDF.txt", overwrite=TRUE)

#=============================================================================#
# Perform differential binding/abundance analysis of peaks:

args_diff = systemArgs(sysma="../../param/rundiff.param", 
                       mytargets="targets_countDF.txt")
cmp = readComp(file=args_bam, format="matrix")
dbrlist = runDiff(args=args_diff, diffFct=run_edgeR, 
                  targets=targetsin(args_bam),cmp=cmp[[1]], 
                  independent=TRUE, dbrfilter=c(Fold=2, FDR=1))
writeTargetsout(x=args_diff, file="targets_rundiff.txt", overwrite=TRUE)

#=============================================================================#
# Perform GO term enrichment analysis:

# Two options are given for GO enrichment here:
 
# option 1:
# unzip gene-to-GO mappings file for organism of interest
# can be found on http://geneontology.org/GO.downloads.annotations.shtml
# use the unzipped file to create a CAT database
#system2("gunzip", args = "goa_human.gaf.gz")
#catdb = makeCATdb("goa_human.gaf")

# option 2:
# create catDB using biomaRt:
listMarts(host="uswest.ensembl.org") # 
m = useMart(biomart="ENSEMBL_MART_ENSEMBL")
x = listAttributes(m)
View(x) # view options for data types to download
go = getBM(attributes=c("go_id", "ensembl_gene_id", "namespace_1003"), mart=m)
go = go[go[,3]!="",]
go[,3] = as.character(go[,3])
write.table(go, "GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE,
            col.names=FALSE, sep="\t")

## create catDB can takes some time, but need only be done once
catdb = makeCATdb(myfile="GOannotationsBiomart_mod.txt", lib=NULL, org="", 
                   colno=c(1,2,3), idconv=NULL)

save(catdb, file="catdb.RData") 
load("catdb.RData")

# perform GO enrichment
args = systemArgs(sysma="../../param/macs2.param", 
                 mytargets="targets_bam_ref.txt")
args_anno = systemArgs(sysma="../../param/annotate_peaks.param", 
                      mytargets="targets_macs.txt")
annofiles = outpaths(args_anno)
gene_ids = sapply(names(annofiles),
                  function(x) unique(as.character(read.delim(annofiles[x])[,"geneId"])))
BatchResult = GOCluster_Report(catdb=catdb, setlist=gene_ids, 
                              method="all", id_type="gene", CLSZ=2, cutoff=0.9,
                               gocats=c("MF", "BP", "CC"))

# Write the gene IDs to a file for HOMER GO Enrichment analysis
write.table(gene_ids$A1a, "HARGeneIDs.txt", quote = FALSE, sep = "\n", 
            row.names = F, col.names = F)
write.table(gene_ids$E1a, "HERaGeneIDs.txt", quote = FALSE, sep = "\n", 
            row.names = F, col.names = F)
write.table(gene_ids$E2a, "HERbGeneIDs.txt", quote = FALSE, sep = "\n",
            row.names = F, col.names = F)
write.table(gene_ids$C1a, "HinputGeneIDs.txt", quote = FALSE, sep = "\n", 
            row.names = F, col.names = F)

# Prepare the data frame for each receptor
HAR = read.delim("results/HARGOEnrich/biological_process.txt")
HAR = HAR[HAR[,6]> 3,] # control the threshold for plotting
CLID1 = replicate(nrow(HAR), "HAR")
Ont1 = replicate(nrow(HAR), "BP")
HAR = cbind(CLID1, HAR, Ont1)
colnames(HAR) = c("CLID","GOID", "Term", "Enrichment", "logP", 
                  "Genes.in.Term", "SampleMatch","", "", "", "", "", "Ont")

HERa = read.delim("results/HERaGOEnrich/biological_process.txt")
HERa = HERa[HERa[,6]> 3,] # control the threshold for plotting
CLID2 = replicate(nrow(HERa), "HERa")
Ont2 = replicate(nrow(HERa), "BP")
HERa = cbind(CLID2, HERa, Ont2)
colnames(HERa) = c("CLID","GOID", "Term", "Enrichment", "logP", 
                  "Genes.in.Term", "SampleMatch","", "", "", "", "", "Ont")

HERb = read.delim("results/HERbGOEnrich/biological_process.txt")
HERb = HERb[HERb[,6]> 3,] # control the threshold for plotting
CLID3 = replicate(nrow(HERb), "HERb")
Ont3 = replicate(nrow(HERb), "BP")
HERb = cbind(CLID3, HERb, Ont3)
colnames(HERb) = c("CLID","GOID", "Term", "Enrichment", "logP", 
                  "Genes.in.Term", "SampleMatch","", "", "", "", "", "Ont")

Hinput = read.delim("results/HinputGOEnrich/biological_process.txt")
Hinput = Hinput[Hinput[,6]> 3,] # control the threshold for plotting
CLID4 = replicate(nrow(Hinput), "Hinput")
Ont4 = replicate(nrow(Hinput), "BP")
Hinput = cbind(CLID4, Hinput, Ont4)
colnames(Hinput) = c("CLID","GOID", "Term", "Enrichment", "logP", 
                  "Genes.in.Term", "SampleMatch","", "", "", "", "", "Ont")

pdf("./results/goBarPlot.pdf", height = 5, width=8)
# create a GO plot for most enriched terms
All_groups = rbind(HAR, HERa, HERb, Hinput)
goBarplot(All_groups, gocat = "BP")
dev.off()

#=============================================================================#
# Perform Motif analysis:

args = systemArgs(sysma="../../param/annotate_peaks.param",
                  mytargets="targets_macs.txt")

# motif analysis on the peaks files
rangefiles = infile1(args)
for(i in seq(along=rangefiles)){
    df = read.delim(rangefiles[i], comment="#")
    peaks = as(df, "GRanges")
    names(peaks) = paste0(as.character(seqnames(peaks)), "_", start(peaks), "-", end(peaks))
    peaks = peaks[order(values(peaks)$X.log10.pvalue, decreasing=TRUE)]
    pseq = getSeq(FaFile("hg19.fa"), peaks)
    names(pseq) = names(peaks)
    writeXStringSet(pseq, paste0(rangefiles[i], ".fa"))
}

# find the highest ranked motifs
set.seed(0)
BCRANKout = bcrank(paste0(rangefiles[1], ".fa"), 
                    restarts=25, use.P1=TRUE, use.P2=TRUE)
topMotifs1 = toptable(BCRANKout)
topMotifs1

# save the highest ranked motifs for the first sample and save them 
pdf("results/seqlogo1.pdf")
for (i in seq(1, dim(topMotifs1)[1], by=1)){
    topMotif = toptable(BCRANKout, i)
    weightMatrix = pwm(topMotif, normalize = FALSE)
    weightMatrixNormalized = pwm(topMotif, normalize = TRUE)
    seqLogo(weightMatrixNormalized)
}
dev.off() #close the graphics to the pdf file

# find the highest ranked motifs
set.seed(0)
BCRANKout = bcrank(paste0(rangefiles[2], ".fa"), 
                   restarts=25, use.P1=TRUE, use.P2=TRUE)
topMotifs2 = toptable(BCRANKout)
topMotifs2

# save the highest ranked motifs for the second sample and save them 
pdf("results/seqlogo2.pdf")
for (i in seq(1, dim(topMotifs2)[1], by=1)){
  topMotif = toptable(BCRANKout, i)
  weightMatrix = pwm(topMotif, normalize = FALSE)
  weightMatrixNormalized = pwm(topMotif, normalize = TRUE)
  seqLogo(weightMatrixNormalized)
}
dev.off() #close the graphics to the pdf file

# find the highest ranked motifs
set.seed(0)
BCRANKout = bcrank(paste0(rangefiles[3], ".fa"), 
                   restarts=25, use.P1=TRUE, use.P2=TRUE)
topMotifs3 = toptable(BCRANKout)
topMotifs3

# save the highest ranked motifs for the third sample and save them 
pdf("results/seqlogo3.pdf")
for (i in seq(1, dim(topMotifs3)[1], by=1)){
  topMotif = toptable(BCRANKout, i)
  weightMatrix = pwm(topMotif, normalize = FALSE)
  weightMatrixNormalized = pwm(topMotif, normalize = TRUE)
  seqLogo(weightMatrixNormalized)
}
dev.off() #close the graphics to the pdf file

# find the highest ranked motifs
set.seed(0)
BCRANKout = bcrank(paste0(rangefiles[4], ".fa"), 
                   restarts=25, use.P1=TRUE, use.P2=TRUE)
topMotifs4 = toptable(BCRANKout)
topMotifs4

# save the highest ranked motifs for the fourth sample and save them 
pdf("results/seqlogo4.pdf")
for (i in seq(1, dim(topMotifs4)[1], by=1)){
  topMotif = toptable(BCRANKout, i)
  weightMatrix = pwm(topMotif, normalize = FALSE)
  weightMatrixNormalized = pwm(topMotif, normalize = TRUE)
  seqLogo(weightMatrixNormalized)
}
dev.off() #close the graphics to the pdf file

#=============================================================================#
# Packages used:

toLatex(sessionInfo())
