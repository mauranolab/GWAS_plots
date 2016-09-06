#!/usr/bin/env Rscript
#Author: Matt Maurano

argv <- commandArgs(TRUE)
f <- argv[1]
pthresh <- as.numeric(argv[2])
numLabels <- as.numeric(argv[3])
filename <- argv[4]


x_lim <- -log(pthresh, base=10) *1.1
y_lim <- NULL #c(0,6) #NB ggplot clips line segments connecting points exceeding limits

cat("Reading", f, "\n")
tmp <- read(f, header=F)
colnames(tmp) <- c("Enrichment", "Fisher_P-value", "DHS", "Study", "GWAS_P_threshold", "Expected_overlaps", "Observed_overlaps", "Total_GWAS_SNPs_below_threshold")

cat("Plotting", filename, " at P>", pthresh, "\n")


#Manual specification of legend by regex on sample name
#NB if legendSamples covers everything, need to cmment out grey line in plot below
legendSamples <- list(
#	"heart"=list(samples="Heart", color="red"),
#	"cardiomyocytes"=list(samples="H7_hESC_T14|HCM", color="green"),
#	"fibroblasts: cardiac/pulmonary/vessel"=list(samples="HCF|HCFaa|HPAF|AoAF|HPF", color="blue"),
#	"fibroblasts: other"=list(samples="BJ-|IMR90|AG|HConF|HFF|HMF|HPdLF|HVMF|NHDF|NHLF|Fibroblasts", color="lightblue"),
	"vascular endothelia"=list(samples="HBMEC|HMVEC|HPAEC|HRGEC|HUVEC|HBVSMC|HBVP", color="orange"),
	"muscle"=list(samples="fMuscle|HSMM|Muscle|LHCN_M2_D4", color="yellow"),
#	"thymus"=list(samples="Thymus", color="pink"),
	"immune"=list(samples="CD\\d|GM\\d\\d\\d\\d\\d|Jurkat|hT[HR]|iTH", color="red")
#	"GI"=list(samples="Intestine|Colon", color="blue")
)

tmp$celltype <- "other"
celltype.colors <- NULL
for(cur in sort(c("other", names(legendSamples)))) {
	if(cur!="other") {
		curParams <- legendSamples[[cur]]
		numSamples <- nrow(subset(tmp[grepl(curParams$samples, tmp$DHS, ignore.case=T, perl=T),], GWAS_P_threshold == pthresh))
		
		#cat(cur, "n=", numSamples, "\n")
		
		tmp[grepl(curParams$samples, tmp$DHS, ignore.case=T, perl=T), "celltype"] <- paste(cur, " (n=", numSamples, ")", sep="")
		celltype.colors <- append(celltype.colors, curParams$color)
	} else {
		#ggplot doesnt respect factor order, only alphabetical
		#TODO use scale_fill_discrete(breaks=c(levels...)
		celltype.colors <- append(celltype.colors, "grey")
	}
}
numSamples <- nrow(subset(tmp[tmp$celltype=="other",], GWAS_P_threshold == pthresh))
tmp[tmp$celltype=="other", "celltype"] <- paste("other (n=", numSamples, ")", sep="")
tmp$celltype <- factor(tmp$celltype)


tmp <- subset(tmp, GWAS_P_threshold >= pthresh)
tmp$logp <- -log(tmp$GWAS_P_threshold, base=10)


library(ggplot2)
library(directlabels)
old <- theme_set(theme_bw(base_size=8))
old <- theme_update(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.x=theme_text(size=6, hjust=1, vjust=1))# , axis.ticks.length=unit(0.1, "lines") )
#old <- theme_update(plot.margin=unit(c(0,7,0,0), "lines")) #NB top, rt, bot, left

labels.cex <- 0.5


#Set up labels
#tmp.labelFrame <- subset(tmp, GWAS_P_threshold == pthresh)

#enable labeling points that don't get to max pthresh
tmp.labelFrame <- tmp[ apply(subset(tmp, select=c(DHS, GWAS_P_threshold)), MARGIN=1, FUN=paste, collapse=".") %in% apply(summaryBy(GWAS_P_threshold~DHS, data=tmp, FUN=min), MARGIN=1, FUN=paste, collapse="."), ]

tmp.labelFrame <- tmp.labelFrame[order(-tmp.labelFrame$Enrichment),]
tmp.labelFrame <- tmp.labelFrame[1:numLabels,]


#Plot
doPlot <- function () {
	print(ggplot() + # & Expected_overlaps>20
	geom_line(data=subset(tmp, grepl("^other", celltype)), aes(x=logp, y=Enrichment, group=DHS, color=celltype)) +
	geom_line(data=subset(tmp, !grepl("^other", celltype)), aes(x=logp, y=Enrichment, group=DHS, color=celltype)) +
	geom_dl(data=tmp.labelFrame, aes(x=logp, y=Enrichment, color=celltype, label=gsub("-DS.+$", "", DHS)), list(last.points, cex=labels.cex, hjust=0)) +
guides(color = guide_legend(override.aes = list(size=5, linetype=1))) +#hack to get legend symbol to be square not dot
	scale_color_manual(name="", values=celltype.colors) +
	opts(legend.position = c(0.13, 0.82)) +
	scale_x_continuous("GWAS P-value threshold", limits=c(0, x_lim), expand=c(0.1,0)) + scale_y_continuous("Fold enrichment of SNPs in DHSs", limits=y_lim, expand=c(0,0.05)))
}


pdf(file=paste(filename, ".pdf", sep=""), width=8, height=6)
doPlot()
dev.off()


#resize for png
old <- theme_set(theme_bw(base_size=12))
old <- theme_update(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.x=theme_text(size=10, hjust=1, vjust=1))# , axis.ticks.length=unit(0.1, "lines") )

labels.cex <- 1

png(file=paste(filename, ".png", sep=""), width=1000, height=600)
doPlot()
dev.off()
