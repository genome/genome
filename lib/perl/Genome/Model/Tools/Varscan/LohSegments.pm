
package Genome::Model::Tools::Varscan::LohSegments;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Varscan::FilterVariantCalls    Process somatic pipeline output
#
#    AUTHOR:     Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#    CREATED:    12/09/2009 by D.K.
#    MODIFIED:   12/29/2009 by D.K.
#
#    NOTES:
#
#####################################################################################################################################

use strict;
use warnings;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Varscan::LohSegments {
    is => 'Command',
    has => [                                # specify the command's single-value properties (parameters) <---
        variant_file => {
            is => 'Text',
            doc => "File containing LOH and Germline calls in VarScan-annotation format",
            is_optional => 0,
            is_input => 1,
        },
        min_coverage => {
            is => 'Text',
            doc => "Minimum coverage required to include a site",
            is_optional => 0,
            is_input => 1,
            default => 20,
        },
        min_freq_for_het => {
            is => 'Text',
            doc => "Minimum variant allele frequency in normal to consider a heterozygote",
            is_optional => 0,
            is_input => 1,
            default => 40,
        },
        max_freq_for_het => {
            is => 'Text',
            doc => "Maximum variant allele frequency in normal to consider a heterozygote",
            is_optional => 0,
            is_input => 1,
            default => 60,
        },
        output_basename => {
            is => 'Number',
            doc => "Basename for creating output files",
            is_optional => 1,
            is_input => 1,
            default_value => '0.07',
        },
        varscan_cn_basename => {
            is => 'Number',
            doc => "Basename used for VarScan copy number calling, if plotting these is desired",
            is_optional => 1,
            is_input => 1,
        },
    ],
    has_param => [
        lsf_resource => {
            value => 'select[tmp>1000] rusage[tmp=1000]',
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Generate LOH Frequency Plots from VarScan Germline+LOH Data"                 
}

sub help_synopsis {
    return <<EOS
    This command generates LOH Frequency Plots from VarScan Germline+LOH Data
    EXAMPLE:    gmt varscan loh-segments --variant-file varScan.output.LOHandGermline.snp --output-basename varScan.output.LOH.segments ...
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;

    ## Get required parameters ##
    my $variant_file = $self->variant_file;
    my $output_basename = $self->output_basename;

    my $min_freq_for_het = $self->min_freq_for_het;
    my $max_freq_for_het = $self->max_freq_for_het;
    
    ## Correct to percentage values if a fraction frequency provided ##
    $min_freq_for_het = 100 * $min_freq_for_het if($min_freq_for_het > 0 && $min_freq_for_het < 1);
    $max_freq_for_het = 100 * $max_freq_for_het if($max_freq_for_het > 0 && $max_freq_for_het < 1);

    my %stats = ();
    $stats{'num_snps'} = $stats{'num_het'} = $stats{'num_loh'} = $stats{'num_germline'} = $stats{'loh_segments'} = 0;

    ## Convert variant file to infile ##

    open(OUTFILE, ">$output_basename.infile");
    print OUTFILE "chrom\tposition\tnfreq\ttfreq\tstatus\n";

    open(SEGMENTS, ">$output_basename.rawsegments");
    print SEGMENTS "chrom\tchr_start\tchr_stop\tnum_snps\tone\n";

    my $input = new FileHandle ($self->variant_file);
    my $lineCounter = 0;
    
    my $loh_chrom = my $loh_start = my $loh_stop = my $loh_snps = 0;
    
    while (<$input>)
    {
            chomp;
            my $line = $_;
            $lineCounter++;		
    
            my @lineContents = split(/\t/, $line);
            my $chrom = $lineContents[0];
            my $pos = $lineContents[1];
            my $nreads1 = $lineContents[5];
            my $nreads2 = $lineContents[6];
            my $nfreq = $lineContents[7];
            my $treads1 = $lineContents[9];
            my $treads2 = $lineContents[10];
            my $tfreq = $lineContents[11];
            my $status = $lineContents[13];
            my $diff_p_value = $lineContents[15];
            
            my $normal_cov = $nreads1 + $nreads2;
            my $tumor_cov = $treads1 + $treads2;
            
            $stats{'num_snps'}++;

            $nfreq =~ s/\%//;
            $tfreq =~ s/\%//;
            
            if($normal_cov >= $self->min_coverage && $tumor_cov >= $self->min_coverage && $nfreq >= $min_freq_for_het && $nfreq <= $max_freq_for_het)
            {
                $stats{'num_het'}++;
                ## New logic for segmentation ##
                if($status eq "LOH")
                {
                    $status = 1;
                }
                else
                {
                    ## Determine if there's evidence for partial LOH - say 5% frequency difference away from het  ##
                    
                    my $norm_diff_from_het = abs($nfreq - 50);
                    my $tum_diff_from_het = abs($tfreq - 50);
                    
                    ## Determine if tumor is farther from het and if the diff is significant ##
                    
                    if($tum_diff_from_het > $norm_diff_from_het && $diff_p_value < 0.05)
                    {
                        $status = 1;
                    }
                    else
                    {
                        $status = 0;                        
                    }
                    

                }

                print OUTFILE join("\t", $chrom, $pos, $nfreq, $tfreq, $status) . "\n";
                
                ## If the call was LOH ##
                
                if($status eq "LOH")
                {
                    $stats{'num_loh'}++;
                    ## If we have no LOH region, start one ##
                    
                    if(!$loh_snps)
                    {
                        $loh_chrom = $chrom;
                        $loh_start = $pos;
                        $loh_stop = $pos;
                        $loh_snps = 1;
                    }
                    
                    ## If we're on a different chrom, end any LOH regions ##
                    elsif($chrom ne $loh_chrom)
                    {
                        if($loh_snps > 1)
                        {
                            print SEGMENTS join("\t", $loh_chrom, $loh_start, $loh_stop, $loh_snps, -1) . "\n";
                            $stats{'loh_segments'}++;
                        }
                        ## Also start a new region ##
                        $loh_chrom = $chrom;
                        $loh_start = $pos;
                        $loh_stop = $pos;
                        $loh_snps = 1;
                    }
                    else
                    {
                        ## Extend the LOH region ##
                        
                        $loh_stop = $pos;
                        $loh_snps++;
                    }
                    
                }
                elsif($status eq "Germline")
                {
                    $stats{'num_germline'}++;
                    ## A germline heterozygote will end any LOH regions ##
                    
                    if($loh_snps > 1)
                    {
                        print SEGMENTS join("\t", $loh_chrom, $loh_start, $loh_stop, $loh_snps, -1) . "\n";
                        $stats{'loh_segments'}++;
                    }
                    
                    $loh_chrom = $loh_start = $loh_stop = $loh_snps = 0;
                }
                
                ## Otherwise if we h
            }

    }
    
    close($input);	

    ## Process final LOH region if there was one ##

    if($loh_snps > 1)
    {
        print SEGMENTS join("\t", $loh_chrom, $loh_start, $loh_stop, $loh_snps, -1) . "\n";
        $stats{'loh_segments'}++;
    }


    close(OUTFILE);
    close(SEGMENTS);


    ## Open HTML File ##
    open(INDEX, ">$output_basename.index.html") or die "Can't open outfile: $!\n";
    print INDEX "<HTML><BODY><TABLE CELLSPACING=0 CELLPADDING=5 BORDER=0 WIDTH=\"100%\">\n";
    print INDEX "<TR>\n";
    my $num_printed_in_column = 0;

    open(OUTFILE, ">$output_basename.R") or die "Can't open outfile: $!\n";
    print OUTFILE "snp <- read.table(\"$output_basename.infile\", header=TRUE)\n";
    print OUTFILE "lohsegs <- read.table(\"$output_basename.rawsegments\", header=TRUE)\n";
    print OUTFILE "minus1 <- snp\$pos - snp\$pos - 1\n";

    my @cbs_files = ();
    my $cbs_fileCounter = 0;

    for(my $chrCounter = 1; $chrCounter <= 24; $chrCounter++)
    {
        my $chrom = $chrCounter;
        $chrom = "X" if($chrCounter == 23);
        $chrom = "Y" if($chrCounter == 24);

        my $outfile = $output_basename . ".$chrom.png";
        my $outfile_seg = $output_basename . ".$chrom.LOH.seg";
        my $outfile_dist = $output_basename . ".$chrom.dist.png";

        ## Save name of file ##
        $cbs_files[$cbs_fileCounter] = $outfile_seg;
        $cbs_fileCounter++;
        
        print OUTFILE qq{
library(DNAcopy)
CNA.object <- CNA(snp\$status[snp\$chrom=="$chrom"], snp\$chrom[snp\$chrom=="$chrom"], snp\$position[snp\$chrom=="$chrom"], data.type="binary", sampleid=c("Chromosome $chrom"))
if(nrow(CNA.object) != 0) {
    #scatterplot
    png("$outfile", height=300, width=800)
    maxpos <- max(snp\$pos[snp\$chrom=="$chrom"])
    segment.CNA.object <- segment(CNA.object)
    write.table(segment.CNA.object\$output, file="$outfile_seg")
    pvalue.segment.CNA.object <- segment.CNA.object\$output
    detach(package:DNAcopy)
    par(mar=c(4,4,2,2))
    plot(snp\$pos[snp\$chrom=="$chrom"], snp\$nfreq[snp\$chrom=="$chrom"], pch=19, cex=0.25, ylim=c(0,100), xlim=c(1,maxpos), xlab="Position on chr$chrom", ylab="Variant Allele Freq", col="blue")
    points(snp\$pos[snp\$chrom=="$chrom"], snp\$tfreq[snp\$chrom=="$chrom"], pch=19, cex=0.25, ylim=c(0,100), col="green")
    segments(lohsegs\$chr_start[lohsegs\$chrom=="$chrom"], lohsegs\$one[lohsegs\$chrom=="$chrom"], lohsegs\$chr_stop[lohsegs\$chrom=="$chrom"], lohsegs\$one[lohsegs\$chrom=="$chrom"], col="red")
    segments(pvalue.segment.CNA.object\$loc.start[pvalue.segment.CNA.object\$seg.mean>0.10], pvalue.segment.CNA.object\$seg.mean[pvalue.segment.CNA.object\$seg.mean>0.10] * 100, pvalue.segment.CNA.object\$loc.end[pvalue.segment.CNA.object\$seg.mean>0.10], pvalue.segment.CNA.object\$seg.mean[pvalue.segment.CNA.object\$seg.mean>0.10] * 100, col="red", lwd=2)
    hetDiff <- abs(snp\$nfreq - snp\$tfreq)
    points(snp\$pos[snp\$chrom=="$chrom" & snp\$status=="1"], minus1[snp\$chrom=="$chrom" & snp\$status=="1"], pch=19, cex=0.5, ylim=c(0,100), col="red")
    legend("topright", legend = c("Normal","Tumor","LOH"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5), col=c("blue","green","red"))
    dev.off()
    #plot the distribution of VAFs
    png("$outfile_dist", height=300, width=800)
    par(mar=c(4,4,2,2))
    freqDistNormal <- table(cut(snp\$nfreq[snp\$chrom=="$chrom"], seq(0,105,by=5), left=FALSE, right=FALSE))
    rownames(freqDistNormal) <- seq(0,100,by=5)
    freqDistTumor <- table(cut(snp\$tfreq[snp\$chrom=="$chrom"], seq(0,105,by=5), left=FALSE, right=FALSE))
    rownames(freqDistTumor) <- seq(0,100,by=5)
    plot(prop.table(freqDistNormal), col="blue", type="l", ylim=c(0, 0.50), xlab="Variant Allele Frequency", ylab="Fraction of SNPs")
    lines(prop.table(freqDistTumor), col="green", type="l")
    legend("topright", legend = c("Normal","Tumor"), lty=c(1,1), lwd=c(2.5,2.5), col=c("blue","green"))
    dev.off()
}
};

#points(snp\$pos[snp\$chrom=="$chrom" & snp\$status=="1"], hetDiff[snp\$chrom=="$chrom" & snp\$status=="1"], pch=19, cex=0.5, ylim=c(0,100), col="red")


        ## Get Image filename ##
        my @temp = split(/\//, $outfile);
        my $numdirs = @temp;
        my $image_filename = $temp[$numdirs - 1];

        my @temp2 = split(/\//, $outfile_dist);
        $numdirs = @temp2;
        my $image_filename_dist = $temp2[$numdirs - 1];

        if($self->varscan_cn_basename)
        {
            my $copy_number_filename = $self->varscan_cn_basename . ".$chrom.jpg";
            
            print INDEX "<TD style=\"border: 1px solid black; text-align: center\"><strong>CHROMOSOME $chrom</strong><BR><A HREF=\"$copy_number_filename\"><IMG SRC=\"$copy_number_filename\" HEIGHT=240 WIDTH=320 BORDER=0></A><BR><A HREF=\"$image_filename\"><IMG SRC=\"$image_filename\" HEIGHT=240 WIDTH=320 BORDER=0></A><BR><A HREF=\"$image_filename_dist\"><IMG SRC=\"$image_filename_dist\" HEIGHT=240 WIDTH=320 BORDER=0></A></TD>\n";            
        }
        else
        {
            print INDEX "<TD style=\"border: 1px solid black; text-align: center\"><strong>CHROMOSOME $chrom</strong><BR><A HREF=\"$image_filename\"><IMG SRC=\"$image_filename\" HEIGHT=240 WIDTH=320 BORDER=0></A><BR><A HREF=\"$image_filename_dist\"><IMG SRC=\"$image_filename_dist\" HEIGHT=240 WIDTH=320 BORDER=0></A></TD>\n";            
        }


        $num_printed_in_column++;

        if($num_printed_in_column >= 4)
        {
                print INDEX "</TR><TR>\n";
                $num_printed_in_column = 0;
        }

    }
    
    close(OUTFILE);


    system("R --no-save < $output_basename.R");

    print INDEX "</TR></TABLE></BODY></HTML>\n";
    close(INDEX);
    
    ## Open File for CBS-loh segments ##
    open(LOHSEG, ">$output_basename.segments.cbs") or die "Can't open outfile: $!\n";
    print LOHSEG "chrom\tchr_start\tchr_stop\tnum_markers\tseg_mean\n";
    foreach my $cbs_file (@cbs_files)
    {
        if(-e $cbs_file)
        {
            my $segments = parse_segmented_loh($cbs_file);
            print LOHSEG "$segments\n";
        }
    }
    
    close(LOHSEG);

    return 1;
}



################################################################################################
# 
#
################################################################################################

sub parse_segmented_loh
{
	my $FileName = shift(@_);

        my $segments = "";

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

                $line =~ s/\"//g;
                
                if($lineCounter > 1)
                {
        		my ($id, $name, $chrom, $start, $stop, $num_mark, $seg_mean) = split(/\s+/, $line);
                        $segments .= "\n" if($segments);
                        $segments .= join("\t", $chrom, $start, $stop, $num_mark, $seg_mean);
                }
	}
	
	close($input);	

        return($segments);
}

1;
