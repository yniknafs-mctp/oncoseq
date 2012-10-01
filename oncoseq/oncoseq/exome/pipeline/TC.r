#!/usr/bin/Rscript

suppressPackageStartupMessages(library(DNAcopy,quietly=TRUE));
suppressPackageStartupMessages(library(ExomeCNV,quietly=TRUE));
library(methods); #OSX doesn't include this for Rscript for some reason

die = function(x)
{
  #die function
  write(x,stderr());
  q(save="no",status=-1);
}
':=' = function(lhs,rhs)
{
  frame = parent.frame();
  lhs = as.list(substitute(lhs));
  if (length(lhs) > 1)
    lhs = lhs[-1];
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame);
    return(invisible(NULL)); }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs = list(rhs);
  if (length(lhs) > length(rhs))
    rhs = c(rhs, rep(list(NULL), length(lhs) - length(rhs)));
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame);
  return(invisible(NULL));
}

args = commandArgs(trailingOnly = TRUE);
if(length(args) < 3) { die("Usage: TC.r [normal baf] [tumor baf] [output file]"); }
c(normal.baf.fp,tumor.baf.fp,output.fp) := args;

normal.baf = read.delim(normal.baf.fp, header=T)
tumor.baf  = read.delim(tumor.baf.fp,  header=T)

eLOH = LOH.analyze(normal.baf, tumor.baf, alpha=0.01, method="two.sample.fisher")
colnames(eLOH)[colnames(eLOH) == "normal.baf"] = "normal.b.coverage"
colnames(eLOH)[colnames(eLOH) == "tumor.baf"]  = "tumor.b.coverage"

eLOH$normal.baf = eLOH$normal.b.coverage / eLOH$normal.coverage
eLOH$tumor.baf  = eLOH$tumor.b.coverage  / eLOH$tumor.coverage

eLOH$contam = NA
eLOH[eLOH$normal.baf >  eLOH$tumor.baf,"contam"] = 2 *      eLOH[eLOH$normal.baf >  eLOH$tumor.baf,"tumor.baf"]
eLOH[eLOH$normal.baf <= eLOH$tumor.baf,"contam"] = 2 * (1 - eLOH[eLOH$normal.baf <= eLOH$tumor.baf,"tumor.baf"])

write(mean(eLOH[eLOH$LOH == TRUE,"contam"]),output.fp)
#write(median(eLOH[eLOH$LOH == TRUE,"contam"]),stdout())

#par(mfrow=c(3,2))
#hist(eLOH$normal.baf,breaks=20)
#hist(eLOH$tumor.baf,breaks=20)
#hist(eLOH[eLOH$LOH == TRUE,"normal.baf"],breaks=20)
#hist(eLOH[eLOH$LOH == TRUE,"tumor.baf"],breaks=20)
#hist(eLOH[eLOH$LOH == TRUE,"contam"],breaks=20)
