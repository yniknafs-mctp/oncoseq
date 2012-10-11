#!/usr/bin/Rscript
################################################################################
# LOH/CNV analysis using ExomeCNV                                              #
# Configured by Brendan Veeneman, 3/2012                                       #
#                                                                              #
# Usage: ExomeCNV.r [normal coverage]                                          #
#                   [tumor coverage]                                           #
#                   [output_loh]                                               #
#                   [output_cnv]                                               #
#                   [output_plot]                                              #
#                   [R backup file]                                            #
#                   [plot title]                                               #
#                   [read length]                                              #
#                   [tumor content estimate]                                   #
#                   (silent)                                                   #
#                                                                              #
################################################################################
#  Libraries                                                                   #
################################################################################

suppressPackageStartupMessages(library(DNAcopy,quietly=TRUE));
suppressPackageStartupMessages(library(ExomeCNV,quietly=TRUE));
library(methods); #OSX doesn't include this for Rscript for some reason

################################################################################
#  Parameters - CHANGE STUFF HERE                                              #
################################################################################

silent            = FALSE;
coverage.cutoff   = 10;
exon.spec         = 0.9999;
exon.sens         = 0.9999;
all.spec          = 0.99;
all.sens          = 0.99;
display.CN.cutoff = 1; # = Inf (for everything)
all.sdundo        = c(1, 2);       #default - must be length 2
all.alpha         = c(0.05, 0.01); #default - must be length 2
plot.height       = 800;
plot.width        = 1600;
graphing.offset   = 1000000;

chr.list          = paste("chr",c(1:22,"X","Y"),sep=""); #levels()?
input.header      = FALSE;
input.has.chr     = TRUE;
input.skip        = 1;
coverage.cols     = c("probe","chr","probe_start","probe_end",
                      "targeted.base","sequenced.base","coverage",
                      "average.coverage","base.with..10.coverage");

################################################################################
#  SUBROUTINES                                                                 #
################################################################################

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
'%+=%' = function(lhs,rhs) #XXX doesn't work for arrays
{
  frame = parent.frame();
  dest = as.list(substitute(lhs));
  if (length(dest) > 1)
    dest = dest[-1];
  for (i in 1:length(lhs))
    do.call(`=`, list(dest[[i]], lhs[i] + rhs), envir=frame);
  return(invisible(NULL));
}
printTime = function(x)
{
  cat(x);
  cat(": ");
  cat(format(Sys.time(), "%H:%M:%S")); #%x - m/d/y
  cat("\n");
}
dieIfUnsorted = function(x)
{
  sorted = system(paste("perl -ane '
                         if($chr ne $F[1]){ $s = $F[2]; $chr = $F[1]; }
                         elsif($F[2] < $s){ print \"$F[2] !>= $s\n\"; exit 1; }
                         else{ $s = $F[2]; }' ", x));
  
  if(sorted != 0)
  {
    die(paste("Error: Input file:",x,"isn't sorted."));
  }
}
multi.CNV.analyze.b = function (normal, tumor, logR, all.cnv.ls,
                                coverage.cutoff, c, l, sdundo, alpha,
                                min.spec, min.sens, option) 
{
  #unchanged except to return extra table
  stopifnot(length(sdundo) == length(alpha))
  if (is.null(all.cnv.ls)) {
    all.cnv.ls = list()
  }
  for (i in 1:length(sdundo)) {
    cna = CNV.analyze(normal, tumor, logR = logR,
                      coverage.cutoff = coverage.cutoff, c = c,
                      sdundo = sdundo[i], alpha = alpha[i], plot.cnv = FALSE)
    ecnv = classify.eCNV(cna$cnv, cna$cnv, cna$cnv$seg.mean,
                         min.spec = min.spec, min.sens = min.sens,
                         option = option, c = c, l = l)
    all.cnv.ls[[length(all.cnv.ls) + 1]] = ecnv
  }
  c(joined.probes,the.cnv) := combine.CNV.b(all.cnv.ls);
  return(list(joined.probes,the.cnv));
}
combine.CNV.b = function (cnv.ls) 
{
  #this code is generally unchanged, but now saves
  # a copy of the probes with their corresponding segments
  my.cols = c("logR", "ratio", "copy.number", "lower.cutoff", 
              "upper.cutoff", "spec", "sens");
  my.backup.cols = paste("probe.",my.cols,sep="");
  if(!silent){ cat("start\n"); }
  the.cnv = cnv.ls[[1]]; #probes
  the.cnv[,my.backup.cols] = the.cnv[,my.cols];
  for (i in 2:length(cnv.ls)) #iterate through two sets of segments
  {
    cur.cnv = cnv.ls[[i]];
    if(!silent){ cat("set: ", i, "; number of lines: ", nrow(cur.cnv), "; "); }
    for (j in 1:nrow(cur.cnv)) # iterate through segments
    {
      cur.chunk = cur.cnv[j, ]; #grab this segment
      #list of probes in this segment
      cur.probes = (cur.chunk$chr == the.cnv$chr & 
                    cur.chunk$probe_start <= the.cnv$probe_start &
                    cur.chunk$probe_end >= the.cnv$probe_end)
      #list of probes in this segment with CN: 0,2,NA,or =segment
      probe.to.merge = (the.cnv$copy.number %in% c(0, 2, cur.chunk$copy.number) |
                        is.na(the.cnv$copy.number)) & cur.probes 
      the.cnv[probe.to.merge, my.cols] = cur.chunk[my.cols]
    }
    if(!silent){ cat(" number of merging events: ",
                     nrow(cnv.ls[[i - 1]]) - nrow(cnv.ls[[i]]), "\n"); }
  }
  
  #this relies heavily on column positions for expedience
  joined.probes = the.cnv[,c(1:8,16:22,9:15)]; #joined probes
  my.rep.cols = paste("segment.",my.cols,sep="");
  names(joined.probes)[16:22] = my.rep.cols;
  names(joined.probes)[9:15] = my.cols;  
  
  the.cnv = do.merge.cnv.intervals(the.cnv)
  return(list(joined.probes,the.cnv))
}
loadDemoCoverage = function() #useful for debugging
{
  #load chromosomes 19-21 from demo online. overwrites global variable
  chr.list <<- paste("chr",c(19:21),sep=""); # <<- is deep assignment
  
  prefix_n = "http://genome.ucla.edu/~fah/ExomeCNV/data/normal.";
  prefix_t = "http://genome.ucla.edu/~fah/ExomeCNV/data/tumor.";
  
  if(!silent){ cat("downloading demo data for analysis\n"); }
  normal = read.all.coverage(prefix_n, ".coverage", chr.list, header=T);
  tumor  = read.all.coverage(prefix_t, ".coverage", chr.list, header=T);
  
  return(list(normal,tumor)	);
}
loadCoverage = function(normal.fp,tumor.fp)
{
  normal = read.table(normal.fp,header=input.header,skip=input.skip,sep="\t");
  tumor  = read.table(tumor.fp ,header=input.header,skip=input.skip,sep="\t");
  names(normal) = coverage.cols;
  names(tumor)  = coverage.cols;
  if(!input.has.chr) 
  {
    normal[,"chr"] = paste("chr",normal[,"chr"],sep="");
    tumor[,"chr"]  = paste("chr",tumor[,"chr"] ,sep="");
  }
  return(list(normal,tumor));
}
doLOH = function()
{
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
doCNV = function(normal,tumor,logR,read.length,tumor.content)
{
  mono.probes = c();
  
  for (i in 1:length(chr.list))
  {
    #first pass - exon level, per chr, high specificity
    if(!silent){ printTime(chr.list[i]); }
    
    idx = (normal$chr == chr.list[i]); #subsetting into chr
    chr.CNV = classify.eCNV(normal=normal[idx,], tumor=tumor[idx,], logR=logR[idx],
        min.spec=exon.spec, min.sens=exon.sens, option="spec",
        c=(1 - tumor.content), l=read.length);
    
    mono.probes = rbind(mono.probes, chr.CNV); #add calls
  }
  
  if(!silent){ printTime("final combination");}
  #combining step - use exon calls to infer larger chunks
  c(joined.probes,segment.CNV) := multi.CNV.analyze.b(normal=normal, #use all data
                                                      tumor=tumor,
                                                      logR=logR,
                                                      all.cnv.ls=list(mono.probes),
                                                      coverage.cutoff=coverage.cutoff,
                                                      min.spec=all.spec,
                                                      min.sens=all.sens,
                                                      option="auc",
                                                      sdundo = all.sdundo,
                                                      alpha = all.alpha,
                                                      c=(1 - tumor.content),
                                                      l=read.length);
  
  return(list(joined.probes,segment.CNV));
}
plot.CNV.default = function(output.fp,all.CNV,exon.CNV)
{
  pdf(filename = output.fp, width = plot.width, height = plot.height);
  display.quantile = 1;
  
  do.plot.eCNV(all.CNV,
      lim.quantile=display.quantile,
      style="bp",
      bg.cnv=exon.CNV,
      line.plot=TRUE);
  dev.off();
}
plot.CNV.fancy = function(output.fp,all.CNV,exon.CNV,plot.title)
{
  pdf(filename = output.fp, width = plot.width, height = plot.height);
    
  #defining colors dynamically
  max.cn = max(all.CNV$copy.number, na.rm = TRUE);
  reds = rainbow(2^(max.cn - 2), start = 3/4, end = 0)[2^(max.cn:3 - 2)];
  colors = c("blue", "cyan", "green4", reds); #these are consistent
  
  #adjust data so we can plot it all on one axis
  xaxis.ticks = c(0);
  total.size = graphing.offset;
  for (chr in chr.list)
  {
    max = max(exon.CNV[exon.CNV$chr == chr, "probe_end"]);
    all.CNV[all.CNV$chr == chr,"probe_start"]   = all.CNV[all.CNV$chr == chr,"probe_start"] + total.size;
    all.CNV[all.CNV$chr == chr,"probe_end"]     = all.CNV[all.CNV$chr == chr,"probe_end"] + total.size;
    exon.CNV[exon.CNV$chr == chr,"probe_start"] = exon.CNV[exon.CNV$chr == chr,"probe_start"] + total.size;
    exon.CNV[exon.CNV$chr == chr,"probe_end"]   = exon.CNV[exon.CNV$chr == chr,"probe_end"] + total.size;
    total.size = total.size + max + graphing.offset;
    xaxis.ticks[length(xaxis.ticks) + 1] = total.size;
  }
  
  #set up master plot
  plot(c(0,total.size), #x range
      c(-4,4),          #y range
      pch= NA_integer_,
      xlab = "Chromosome",
      ylab = expression(Log[2]~Copy~Number~Ratio),
      xaxs="i",
      axes=FALSE,
      main=plot.title);
  
  #axis(2,at=c(-3,-2,-1,0,1,2,3),labels=c(0.25,0.5,1,2,4,8,16),tick=TRUE);
  axis(2);
  axis(1,at=xaxis.ticks,labels=FALSE);
  xaxis.labels = (xaxis.ticks[1:(length(xaxis.ticks) - 1)] +
                 xaxis.ticks[2:length(xaxis.ticks)]) / 2;
  axis(1,at=xaxis.labels,tick=FALSE,labels = chr.list);
  
  box(which="plot", lty="solid");
  
  #plot data
  bg.col = "darkgray";
  for (chr in chr.list)
  {
    #subsetting
    seg = all.CNV[all.CNV$chr == chr, ];
    bg  = exon.CNV[exon.CNV$chr  == chr &
                   is.finite(exon.CNV$logR) &
                   is.finite(exon.CNV$segment.logR) &
                   !is.na(exon.CNV$segment.copy.number) &
                   abs(exon.CNV$logR - exon.CNV$segment.logR) < display.CN.cutoff, ];
        
    #plot exons/probes
    points(bg$probe_end, bg$logR, pch = 20, col = bg.col);
    if(bg.col == "darkgray") { bg.col = "gray17"; } else { bg.col = "darkgray"; }
    
    #plot segment lines
    for (cn in 0:max(seg$copy.number, na.rm = TRUE))
    {
      seg.cn = seg[seg$copy.number == cn, ];
      if(nrow(seg.cn) != 0)
      {
        sapply(1:nrow(seg.cn),
            function(i) {
              lines(seg.cn[i, c("probe_start", "probe_end")],
                  rep(seg.cn[i,"logR"], 2),
                  lwd = 3,
                  col = colors[cn + 1])
            });
      }
    }
  }
  
  legend.values = sort(unique(all.CNV[is.finite(all.CNV$copy.number) &
                              all.CNV$copy.number > 0,"copy.number"]));
  legend("bottomright",
      #title="Estimated Tumor Copy Number",
      title="Tumor CN",
      as.expression(legend.values),
      #horiz=TRUE,
      col=colors[legend.values + 1],
      pch=19);
    
  graphics.off();
}

################################################################################
#  MAIN                                                                        #
################################################################################

args = commandArgs(trailingOnly = TRUE);
if(length(args) < 9)
{
  die("Usage: ExomeCNV.r [normal coverage]
          [tumor coverage]
          [output_loh]
          [output_cnv]
          [output_plot]
          [R backup file]
          [plot title]
          [read length]
          [tumor content estimate]
          (silent)");
}

c(normal.fp,tumor.fp,output.LOH,output.CNV,output.plot,rdata,plot.title) := args;
read.length = as.numeric(args[8]); #(read1 + read2) / 2, rounded down
tumor.content = as.numeric(args[9]);
if(length(args) > 9){ silent = as.numeric(args[10]); }
if(!silent){ printTime("checking input sortedness"); }
dieIfUnsorted(normal.fp);
dieIfUnsorted(tumor.fp);

if(!silent){ printTime("start"); }

#do LOH

c(normal,tumor) := loadCoverage(normal.fp,tumor.fp);
#c(normal,tumor) := loadDemoCoverage();

logR = calculate.logR(normal, tumor);
c(probe.CNV, segment.CNV) := doCNV(normal, tumor, logR, read.length, tumor.content);

plot.CNV.fancy(output.plot,segment.CNV,probe.CNV,plot.title);

write.table(probe.CNV, file = output.CNV, sep = "\t",
            row.names = FALSE, quote = FALSE);
save(list = ls(),file = rdata);
        
if(!silent){ printTime("end"); }
