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
#                   [read length]                                              #
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

silent           = FALSE;
tumor.content    = 0.7;
coverage.cutoff  = 10;
exon.spec        = 0.9999;
exon.sens        = 0.9999;
all.spec         = 0.99;
all.sens         = 0.99;
display.quantile = 1;
all.sdundo       = c(1, 2);       #default
all.alpha        = c(0.05, 0.01); #default
plot.height      = 800;
plot.width       = 800;
graphing.offset  = 1000000;

chr.list         = paste("chr",c(1:22,"X","Y"),sep=""); #all
input.header     = FALSE;
input.has.chr    = FALSE;
coverage.cols    = c("probe","chr","probe_start","probe_end",
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
'%+=%' = function(lhs,rhs) #doesn't work for arrays, still
{
  #infuriatingly, %+=% is as close to += as I can get
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
loadDemoCoverage = function()
{
  #useful for debugging
  #load chromosomes 19-21 from demo online. overwrites global variable
  chr.list <<- paste("chr",c(19:21),sep=""); #demo, <<- is deep assignment
  
  prefix_n = "http://genome.ucla.edu/~fah/ExomeCNV/data/normal.";
  prefix_t = "http://genome.ucla.edu/~fah/ExomeCNV/data/tumor.";
  
  if(!silent){ cat("downloading demo data for analysis\n"); }
  normal = read.all.coverage(prefix_n, ".coverage", chr.list, header=T);
  tumor  = read.all.coverage(prefix_t, ".coverage", chr.list, header=T);
  
  return(list(normal,tumor)	);
}
loadCoverage = function(normal.fp,tumor.fp)	
{
  normal = read.table(normal.fp,header=input.header,sep="\t");
  tumor  = read.table(tumor.fp ,header=input.header,sep="\t");
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
doCNV = function(normal,tumor,logR,read.length)
{
  exon.CNV = c();
    
  for (i in 1:length(chr.list))
  {
    #first pass - exon level, per chr, high specificity
    if(!silent){ printTime(chr.list[i]); }
    
    idx = (normal$chr == chr.list[i]); #subsetting into chr
    chr.CNV = classify.eCNV(normal=normal[idx,],
                            tumor=tumor[idx,],
                            logR=logR[idx],
                            min.spec=exon.spec,
                            min.sens=exon.sens,
                            option="spec",       #optimize for specificity
                            c=(1 - tumor.content),
                            l=read.length);
    
    exon.CNV = rbind(exon.CNV, chr.CNV); #add calls
  }
  
  if(!silent){ printTime("final combination");}
  #combining step - use exon calls to infer larger chunks
  all.CNV = multi.CNV.analyze(normal=normal, #use all data
                              tumor=tumor,
                              logR=logR,
                              all.cnv.ls=list(exon.CNV),
                              coverage.cutoff=coverage.cutoff,
                              min.spec=all.spec,
                              min.sens=all.sens,
                              option="auc", #optimize area under curve - keep
                              sdundo = all.sdundo,
                              alpha = all.alpha,
                              c=(1 - tumor.content),
                              l=read.length);
  
  return(list(exon.CNV,all.CNV));
}
plot.CNV.default = function(output.fp,all.CNV,exon.CNV)
{
  png(filename = output.fp, width = plot.width, height = plot.height);
  
  do.plot.eCNV(all.CNV,
               lim.quantile=display.quantile,
               style="bp",
               bg.cnv=exon.CNV,
               line.plot=TRUE);
  dev.off();
}
plot.CNV.fancy = function(output.fp,all.CNV,exon.CNV)
{
  png(filename = output.fp, width = plot.width, height = plot.height);
  
  #defining colors dynamically
  max.cn = max(all.CNV$copy.number, na.rm = TRUE);
  reds = rainbow(2^(max.cn - 2), start = 3/4, end = 0)[2^(max.cn:3 - 2)];
  colors = c("blue", "cyan", "green4", reds);
  chr.list = unique(as.character(all.CNV$chr));
  
  #adjust data so we can plot it all on one axis
  total.size = 0;
  for (chr in chr.list)
  {
    max = max(exon.CNV[exon.CNV$chr == chr, "probe_end"]);
    all.CNV[all.CNV$chr == chr,"probe_start"]   = all.CNV[all.CNV$chr == chr,"probe_start"] + total.size;
    all.CNV[all.CNV$chr == chr,"probe_end"]     = all.CNV[all.CNV$chr == chr,"probe_end"] + total.size;
    exon.CNV[exon.CNV$chr == chr,"probe_start"] = exon.CNV[exon.CNV$chr == chr,"probe_start"] + total.size;
    exon.CNV[exon.CNV$chr == chr,"probe_end"]   = exon.CNV[exon.CNV$chr == chr,"probe_end"] + total.size;
    total.size = total.size + max + graphing.offset;
  }
  
  #set up master plot
  plot(c(0,total.size), #x range
       c(-4,4),         #y range
       pch= NA_integer_,
       xlab = "Chromosome",
       ylab = "Log2 Copy Number Ratio",
       axes=FALSE);
  #axis(2,at=c(-3,-2,-1,0,1,2,3),labels=c(0.25,0.5,1,2,4,8,16),tick=TRUE);
  axis(2);
  box(which="plot", lty="solid");

  #plot data
  bg.col = "darkgray";
  for (chr in chr.list)
  {
    #subsetting
    seg = all.CNV[all.CNV$chr == chr, ];
    bg  = exon.CNV[exon.CNV$chr == chr, ];
    
    #quantile limiting the background exons/probes
    lim.logR = quantile(bg$logR[bg$logR != Inf & bg$logR != -Inf], 
                        c(1 - display.quantile, display.quantile),
                        na.rm = TRUE);      
    bg = bg[(bg$logR < lim.logR[2]) & (bg$logR > lim.logR[1]),];
    
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
  
  legend.values = sort(unique(all.CNV[!is.na(all.CNV$copy.number),"copy.number"]));
  legend("bottomright",
      as.expression(legend.values),
      col=colors,
      pch=19);
  
  dev.off();
}

################################################################################
#  MAIN                                                                        #
################################################################################

args = commandArgs(trailingOnly = TRUE);
if(length(args) < 6)
{
  die("Usage: ExomeCNV.r [normal coverage]
                  [tumor coverage]
                  [output_loh]
                  [output_cnv]
                  [output_plot]
                  [read length]
                  (silent)");
}

c(normal.fp,tumor.fp,output.LOH,output.CNV,output.plot) := args;
read.length = as.numeric(args[6]); #(read1 + read2) / 2, rounded down
if(length(args) == 7){ silent = as.numeric(args[7]); }

if(!silent){ printTime("start"); }

#do LOH
#LOH -> guesstimate tumor content

c(normal,tumor) := loadCoverage(normal.fp,tumor.fp);
#c(normal,tumor) := loadDemoCoverage();

logR = calculate.logR(normal, tumor);
c(exon.CNV,all.CNV) := doCNV(normal, tumor, logR, read.length);

plot.CNV.fancy(output.plot,all.CNV,exon.CNV);

write.table(exon.CNV,
            file = output.CNV,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE);

if(!silent){ printTime("end\n"); }
