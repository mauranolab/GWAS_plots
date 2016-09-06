#!/usr/bin/perl -w
#Author: Eric Haugen

use strict;
use Getopt::Std;



sub usage
{
    die "Usage:  $0  -p [pvalues.hg19.[GWASname].bed5/starch] -s [starch file with DHS] -r [output results dir]\n";
    exit(1);
}

my($opt_string) = 'r:s:p:';
our($opt_s, $opt_r, $opt_p);
getopts( "$opt_string" ) or usage();
usage() if !$opt_r || !$opt_s || !$opt_p;


my $starchive = $opt_s;
unless (-f $starchive) {
    warn "can't find $starchive";
    usage();
}
warn "Using starch $starchive\n"; 


my $resultsdir = $opt_r;
warn "Output results to $resultsdir\n"; 


my $pbed = $opt_p;
unless (-f $pbed) {
    warn "can't find pvalues $pbed";
    usage();
}
warn "Using pvalues in $pbed\n"; 


my $basedir = "..";
#my $rscript = "$basedir/src/unthresholded_fisher_lowPvalueSNPsEnrichedInDHSs.R";
#my $bedopsdir = "$basedir/x86_64";
my $excludebed = "$basedir/hg19/ccdsGene.Exons.bed";


#OLD my ($pbed,$starchive,$rscript,$bedopsdir,$excludebed) = @ARGV;
chomp( my $gwasName = `basename $pbed | cut -f3 -d .` );

unless (-d $resultsdir) {
    system("mkdir -p $resultsdir");
}
my $outputDir = "$resultsdir/$gwasName";
unless (-d $outputDir) {
    system("mkdir -p $outputDir");
}


print "Starting run ".(localtime);
print "\n";


my $cmd = "bedmap --ec --echo --echo-map-id $pbed $starchive ";
if (defined($excludebed)) {
    $cmd .= "| bedops --ec -n -0% - $excludebed "
}
warn "$cmd\n";

my $totalSNPs = 0;
my %binTotals = ();    # total SNPs per p-value bin
my %cellTotals = ();    # total cell DHSs overlapped by any SNP
my %binCellOverlaps = ();  # binCellOverlaps per cell type per p-value bin 

my $LN10 = log(10); # only calculate this once


open IN, "$cmd | cut -f5- |" or die "$!";
while (<IN>) {
    chomp;
    my ($pvalue,@celltypes) = split /[\|\;]/;
    ++$totalSNPs;
    my $pbin = calcbin($pvalue);
    increment( \%binTotals, $pbin );
    foreach my $celltype (@celltypes) {
        increment( \%cellTotals, $celltype );
        increment( \%binCellOverlaps, "$pbin-$celltype" );
    }
}
close IN;


print "Processing results ".(localtime);
print "\n";


my @reportingBins = ( 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -15, -20, -25, -30, -50, -75, -100, -200 );
foreach my $bin (@reportingBins) {
    # some may not have been observed exactly, still need to add zero counts
    unless (defined($binTotals{$bin})) {
        $binTotals{$bin} = 0;
    }
}

my @sortedbins = sort { $a <=> $b } keys %binTotals;
my $bincumulative = 0;
my %cumulativeBinTotals = ();
foreach my $bin (@sortedbins) {
    $bincumulative += $binTotals{$bin};
    $cumulativeBinTotals{$bin} = $bincumulative; 
    #warn "cumulative\t$bin\t$bincumulative\n";
}


foreach my $cell (sort keys %cellTotals) {
    # One pass to calculate all the cumulative totals (including for disjoint bins with no cell overlaps included)
    my $totalOverlap = $cellTotals{$cell};
    my %cellCumulativeOverlaps = ();
    my $cellCumulativeSum = 0;
    foreach my $bin (@sortedbins) {
        my $totalThreshold = $cumulativeBinTotals{$bin};
        if ($totalThreshold > 0) {
            my $key = "$bin-$cell";
            my $disjointCount = defined( $binCellOverlaps{$key} ) ? $binCellOverlaps{$key} : 0;
            $cellCumulativeSum += $disjointCount; 
        }
        $cellCumulativeOverlaps{$bin} = $cellCumulativeSum; 
        #warn "debug\t$cell\t$bin\t$cellCumulativeSum\n";
    }
    # then a second pass in reverse order, using only the "selected" rounded p-value bins.

    #my $table="$resultsdir/enrichment_$cell.$gwasName.txt";
    my $table="$outputDir/enrichment_$cell.$gwasName.txt";
    open OUT, "> $table" or die "$!";
    print OUT "Enrichment\tFisher_P-value\tDHS\tStudy\tGWAS_P_threshold\tExpected_overlaps\tObserved_overlaps\tTotal_GWAS_SNPs_below_threshold\n";
    foreach my $bin (@reportingBins) {
        my $totalThreshold = $cumulativeBinTotals{$bin};
        if ($totalThreshold > 0) {
            my $fractionTotal = $totalOverlap / $totalSNPs;
            my $overlapThreshold = $cellCumulativeOverlaps{$bin};
            #
            my $expected = $totalThreshold * $fractionTotal;
            my $expectedRounded = int(0.5 + ($expected * 10)) / 10;
            #
            my $enrichment = $overlapThreshold / $expected;
            my $enrichmentRounded = int(0.5 + ($enrichment * 100)) / 100;
            #
            my $low_Pvalue_in_DHS = $overlapThreshold;
            my $low_Pvalue_outside = $totalThreshold - $overlapThreshold;
            my $high_Pvalue_in_DHS = $totalOverlap - $overlapThreshold;
            my $high_Pvalue_outside = $totalSNPs - ($high_Pvalue_in_DHS + $low_Pvalue_in_DHS + $low_Pvalue_outside);
            #chomp (my $pvalueEnrichment=`$rscript $low_Pvalue_in_DHS $low_Pvalue_outside $high_Pvalue_in_DHS $high_Pvalue_outside | cut -f2 -d " " ` );
            my $pvalueEnrichment = "NA";
            #
            my $threshold = displaybin( $bin );
            print OUT "$enrichmentRounded\t$pvalueEnrichment\t$cell\t$gwasName\t$threshold\t$expectedRounded\t$overlapThreshold\t$totalThreshold\n"; 
        }
    }
    close OUT;
}


print "Making plot ".(localtime);
print "\n";
my($plot_pthresh) = 1e-20;
my($plot_numLabels) = 50;
my($plotCmd)="tail -q -n +2 $outputDir/* | grep -v -f excluded_samples.txt > \$TMPDIR/t.txt; ./GWASPlot.R \$TMPDIR/t.txt $plot_pthresh $plot_numLabels $outputDir/$gwasName";
system $plotCmd;


print "Finished ".(localtime);
print "\n";


sub increment {
    my ($hashref, $key ) = @_;
    unless (defined( $hashref->{$key} ) ) {
        $hashref->{$key} = 0;
    }
    ++$hashref->{$key};
}

sub calcbin {
    my $pvalue = shift;
    if (0 == $pvalue) {
        # ugh
        warn "warning, p-value of 0, replacing with 1e-200\n";
        $pvalue = 1e-200;
    }
    my $log10 = log( $pvalue ) / $LN10;
    my $roundLog = int( $log10 );
    #my $result = exp( $roundLog * log(10) );
    my $result = ($log10 < $roundLog) ? $roundLog : ($roundLog + 1);
    return $result;
}

sub displaybin {
    my $bin = shift;
    my %displayLookup = (
        0 => "1.0",
        -1 => "0.1",
        -2 => "0.01",
        -3 => "0.001",
        -4 => "0.0001",
        -5 => "0.00001",
        -6 => "0.000001",
        -7 => "0.0000001",
        -8 => "0.00000001",
        -10 => "1e-10",
        -15 => "1e-15",
        -20 => "1e-20",
        -25 => "1e-25",
        -30 => "1e-30",
        -50 => "1e-50",
        -75 => "1e-75",
        -100 => "1e-100",
        -200 => "1e-200"
    );
    my $result = exp( $bin * $LN10 );
    if (length("$result") > 10) {
        $result = 0 + "1e$bin";
    }
    if (defined($displayLookup{$bin})) {
        $result = $displayLookup{$bin};
        my $confirm = calcbin( $result * 0.99 );
        if ($confirm != $bin) {
            warn "Error $confirm != $bin\n";
        }
    }
    return $result;
}


