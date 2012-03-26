#!/usr/bin/Rscript
#
# LOH/CNV analysis using ExomeCNV
# Configured by Brendan Veeneman, 3/2012
#

################################################################################
#  Libraries, Argument Parsing
################################################################################

suppressPackageStartupMessages(library(DNAcopy,quietly=TRUE));
suppressPackageStartupMessages(library(ExomeCNV,quietly=TRUE));
library(methods); #OSX doesn't include this for Rscript for some reason
die = function(x) {
  #die function
  write(x,stderr());
  q(save="no",status=-1);
}
':=' = function(lhs, rhs) {
  frame = parent.frame()
  lhs = as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs = lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs = list(rhs)
  if (length(lhs) > length(rhs))
    rhs = c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL)) }

args = commandArgs(trailingOnly = TRUE);
if(length(args) != 0)
{
  die("Usage: ExomeCNV.r [no args]");
}

################################################################################
#  Parameters - CHANGE STUFF HERE
################################################################################

demo=FALSE;
tumor_content = 0.9;
read_length = 101; #technically - (read1 + read2) / 2, rounded down
cnv_cutoff_coverage = 10;
normal_coverage = "normal_coverage.gatk";
tumor_coverage =  "tumor_coverage.gatk";
demo.chr.list=paste("chr",c(19:21),sep=""); #demo set
chr.list=paste("chr",c(1:22,"X","Y"),sep=""); #process all chromosomes
#chr.list="chr13";
chrSpec = 0.9999;
chrSens = 0.9999;
allSpec = 0.99;
allSens = 0.99;
allSDundo = c(1, 2);      #default - haven't messed with these yet
allAlpha = c(0.05, 0.01); #default - haven't messed with these yet
outputName = "output";

################################################################################
#  SUBROUTINES
################################################################################

loadDemoCoverage = function() {
  #useful for debugging
  #load chromosomes 19-21 from demo online. overwrites global variable
  chr.list=paste("chr",c("19","20","21"),sep="");
  
  prefix_n = "http://genome.ucla.edu/~fah/ExomeCNV/data/normal.";
  prefix_t = "http://genome.ucla.edu/~fah/ExomeCNV/data/tumor.";
  
  cat("downloading demo data for analysis\n");
  normal = read.all.coverage(prefix_n, ".coverage", chr.list, header=T);
  tumor  = read.all.coverage(prefix_t, ".coverage", chr.list, header=T);
  
  return(list(normal,tumor)	);
}
loadCoverage = function(normal_fp,tumor_fp) {
  normal = read.coverage.gatk(normal_fp);
  tumor  = read.coverage.gatk(tumor_fp);
  return(list(normal,tumor));
}
doCNV = function(normal, tumor, logR) {
  demo.eCNV = c();
  #first pass - exon level, high specificity.  per chr
  for (i in 1:length(chr.list)) {
    idx = (normal$chr == chr.list[i]); #chr number
    cat(paste("calling CNV for ",chr.list[i],sep="")); cat("\n");
    print(Sys.time());
    ecnv = classify.eCNV(normal=normal[idx,], #subset chr from normal
        tumor=tumor[idx,],                    #ibid
        logR=logR[idx],                       #ibid
        min.spec=chrSpec,
        min.sens=chrSens,
        option="spec", #optimize for specificity - don't recommend changing
        c=(1 - tumor_content),
        l=read_length);
    
    demo.eCNV = rbind(demo.eCNV, ecnv); #iteratively add calls
  }
  
  cat("performing final combination step\n");
  print(Sys.time());
  #combining step - use exon calls to infer larger chunks
  demo.cnv = multi.CNV.analyze(normal=normal, #use all data
                               tumor=tumor,   #ibid
                               logR=logR,     #ibid
                               all.cnv.ls=list(demo.eCNV),
                               coverage.cutoff=cnv_cutoff_coverage,
                               min.spec=allSpec,
                               min.sens=allSens,
                               option="auc", #optimize area under curve - keep
                               sdundo = allSDundo,
                               alpha = allAlpha,
                               c=(1 - tumor_content),
                               l=read_length);
  
  do.plot.eCNV(demo.cnv,
               lim.quantile=0.99,
               style="bp",
               bg.cnv=demo.eCNV,
               line.plot=T);
  write.output(demo.eCNV, demo.cnv, outputName);
}
doLOH = function() {
#LOH code
#normal = read.delim("http://genome.ucla.edu/~fah/ExomeCNV/data/normal.baf.txt", header=T)
#tumor = read.delim("http://genome.ucla.edu/~fah/ExomeCNV/data/tumor.baf.txt", header=T)
#eLOH = LOH.analyze(normal, tumor, alpha=0.05, method="two.sample.fisher")
#the.loh = multi.LOH.analyze(normal, tumor, all.loh.ls=list(eLOH), test.alpha=0.001, method="variance.f", sdundo=c(0,0), alpha=c(0.5,0.1))
#do.plot.loh(the.loh, normal, tumor, "two.sample.fisher", plot.style="baf")
#write.loh.output(the.loh, "demo.eloh")
#expanded.loh = expand.loh(the.loh, normal)
#write.loh.output(expanded.loh, "demo.all")
}

################################################################################
#  MAIN
################################################################################

cat("start\n");
print(Sys.time());

#do LOH
#LOH -> guesstimate tumor content

if(demo) {
  chr.list = demo.chr.list;
  c(normal,tumor) := loadDemoCoverage();
} else {
  c(normal,tumor) := loadCoverage(normal_coverage,tumor_coverage);
}

logR = calculate.logR(normal, tumor);

doCNV(normal, tumor, logR);

cat("end\n");
print(Sys.time());


#