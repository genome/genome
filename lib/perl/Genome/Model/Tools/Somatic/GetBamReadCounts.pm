package Genome::Model::Tools::Somatic::GetBamReadCounts;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Somatic::GetBamReadCounts {
    is => 'Command',
    has => [

    ## INPUT/OUTPUT OPTIONS ##
        'bam_file' => {
            type => 'String',
            doc => 'The BAM file in which to examine reads ',
            is_input => 1,
        },
        'variant_file' => {
            type => 'String',
            is_input => 1,
            doc => 'List of variant positions in annotation format',
        },
        'output_file' => {
            type => 'String',
            is_input => 1,
            is_output => 1,
            doc => 'Filename for variants that pass filter',
        },
        'reference' => {
            type => 'String',
            default => '/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa',
            is_optional => 1,
            is_input => 1,
            doc => 'Reference sequence to use',
        },
        ## WGS FILTER OPTIONS ##


        ## SHARED OPTIONS ##
        prepend_chr => {
            is => 'Boolean',
            default => '0',
            is_optional => 1,
            is_input => 1,
            doc => 'prepend the string "chr" to chromosome names. This is primarily used for external/imported bam files.',
        },
        verbose => {
            is => 'Boolean',
            default => '0',
            is_optional => 1,
            is_input => 1,
            doc => 'Print the filtering result for each site.',
        },
        # Make workflow choose 64 bit blades
        lsf_resource => {
            is_param => 1,
            default_value => 'rusage[mem=4000,tmp=1000] select[type==LINUX64 && tmp>1000] span[hosts=1]',
        },
        lsf_queue => {
            is_param => 1,
            default_value => 'long',
        },
        skip => {
            is => 'Boolean',
            default => '0',
            is_input => 1,
            is_optional => 1,
            doc => "If set to true... this will do nothing! Fairly useless, except this is necessary for workflow.",
        },
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
        },
        skip_if_MT => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'Enable this flag to skip evaluation (and auto-pass) mitochondrial sites',
        },

        samtools_version => {
            is => 'Text',
            is_optional => 1,
            is_input => 1,
            doc => 'version of samtools to use',
        },
    ]
};

sub help_brief {
    return "This module uses detailed readcount information from bam-readcounts to filter likely false positives";
}

sub help_synopsis {
    return <<EOS
    gmt somatic get-bam-read-counts \
      --variant-file somatic.snvs --bam-file tumor.bam \
      --output-file somatic.snvs.bam-readcounts
EOS
}

sub help_detail {
    return <<EOS
This module retrieves read counts from a BAM file and calculates the allele frequency
For questions, e-mail Dan Koboldt (dkoboldt\@genome.wustl.edu) or Dave Larson (dlarson\@genome.wustl.edu)
EOS
}

sub execute {
    my $self = shift;

    if ($self->skip) {
        $self->status_message("Skipping execution: Skip flag set");
        return 1;
    }

    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    #test architecture to make sure we can run read count program
    unless (`uname -a` =~ /x86_64/) {
       $self->error_message("Must run on a 64 bit machine");
       die;
    }

    #check on BAM file
    unless(-e $self->bam_file) {
        $self->error_message("BAM file: " . $self->bam_file . " does not exist");
        die;
    }

    unless(-e $self->bam_file . ".bai") {
        $self->error_message("BAM must be indexed");
        die;
    }

    ## Run the FP filter. Note that both WGS and capture use the same filter now ##
    $self->run_filter();

}


##########################################################################################
# Capture filter for high-depth, lower-breadth datasets
# Contact: Dan Koboldt (dkoboldt@genome.wustl.edu)
##########################################################################################

sub run_filter {
    my $self = shift(@_);

    ## Reset counters ##

    my %stats = ();
    $stats{'num_variants'}  = $stats{'num_no_readcounts'} = 0;

    ## Open the output file ##

    my $temp_output_file = Genome::Sys->create_temp_file_path;
    my $ofh = Genome::Sys->open_file_for_writing($temp_output_file);
    unless($ofh)
    {
        $self->error_message("Unable to open " . $self->output_file . " for writing.");
        die;
    }

    ## Run BAM readcounts in batch mode to get read counts for all positions in file ##
    my $readcount_file;
    $self->status_message('Running BAM Readcounts...');

    #First, need to create a variant list file to use for generating the readcounts.
    my $input = Genome::Sys->open_file_for_reading($self->variant_file);

    unless($input)
    {
        $self->error_message("Unable to open " . $self->variant_file . ".");
        die;
    }

    ## Build temp file for positions where readcounts are needed ##
    my ($tfh,$temp_path) = Genome::Sys->create_temp_file;
    unless($tfh)
    {
        $self->error_message("Unable to create temporary file.");
        die;
    }
    $temp_path =~ s/\:/\\\:/g;

    ## Print each line to file, prepending chromosome if necessary ##
    $self->status_message('Printing variants to temp file...');

    while(my $line = $input->getline)
    {
        chomp $line;
        my ($chr, $start, $stop) = split /\t/, $line;
        next unless($start =~ /^\d+$/); #header line

        if ($self->prepend_chr)
        {
            $chr = "chr$chr";
            $chr =~ s/MT$/M/;
        };

        if($stop =~ /^\d+$/)
        {
            #annotation format
            print $tfh "$chr\t$start\t$stop\n";
        }
        else
        {
            #varscan format
            print $tfh "$chr\t$start\t$start\n";
        }
    }
    $tfh->close;
    close($input);

    $readcount_file = Genome::Sys->create_temp_file_path;

    print "Running bam-readcounts on tumor BAM...\n";
    my $cmd = $self->readcount_program() . " -b 15 " . $self->bam_file . " -l $temp_path";
    Genome::Sys->shellcmd(
        cmd => "$cmd > $readcount_file 2> /dev/null",
        input_files => [$self->bam_file],
        output_files => [$readcount_file],
    );


    my $readcounts = Genome::Sys->read_file($readcount_file);
    chomp($readcounts) if($readcounts);

    ## Load the results of the readcounts ##

    my %readcounts_by_position = ();

    ## Open the readcounts file ##
    Genome::Sys->copy_file($readcount_file, $self->output_file . ".readcounts");
    
    my @readcounts = split(/\n/, $readcounts);
    
    foreach my $rc_line (@readcounts)
    {
        (my $chrom, my $pos) = split(/\t/, $rc_line);
        $readcounts_by_position{"$chrom\t$pos"} = $rc_line;
    }

    $self->status_message('Readcounts loaded');

    ## Reopen file for parsing ##
    $input = Genome::Sys->open_file_for_reading($self->variant_file);

    ## Parse the variants file ##
    my $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;

        $stats{'num_variants'}++;

            (my $chrom, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
            next unless($chr_start =~ /^\d+$/); #header line
            unless($chr_stop =~ /^\d+$/)
            {
                my @rest;
                ($chrom, $chr_start, $ref, $var, @rest) = split(/\t/, $line); #varscan snp format doesn't list stop separately
                $chr_stop = $chr_start;

                #convert to annotation format
                $line = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var, @rest);
            }

            $ref = uc($ref);
            $var = uc($var);

            my @lineContents = split(/\t/, $line);
            my $numContents = @lineContents;
            my $rest_of_line = "";
            
            for(my $colCounter = 5; $colCounter < $numContents; $colCounter++)
            {
                $rest_of_line .= "\t" if($rest_of_line);
                $rest_of_line .= $lineContents[$colCounter];
            }

            my $query_string = "";

            if($self->prepend_chr)
            {
                $query_string = "chr" . $chrom . ":" . $chr_start . "-" . $chr_stop;
            }
            else
            {
                $query_string = $chrom . ":" . $chr_start . "-" . $chr_stop;
            }

            ## if the variant allele is an IUPAC code, convert it: ##

            if(!($var =~ /[ACGT]/))
            {
                $var = $self->iupac_to_base($ref, $var);
            }

            if($var =~ /[ACGT]/)
            {

                ## Skip MT chromosome sites,w hich almost always pass ##
                if($self->skip_if_MT && ($chrom eq "MT" || $chrom eq "chrMT")) {
                    ## Auto-pass it to increase performance ##
                    $stats{'num_MT_sites_autopassed'}++;
                    print $ofh "$line\n";
#                } elsif($self->prepend_chr && $chrom =~ "random") {
#                    $stats{'num_random_sites_autopassed'}++;
#                    print $ofh "$line\n";
                }
                else
                {
                    ## Run Readcounts ##
#                    my $cmd = $self->readcount_program() . " -b 15 " . $self->bam_file . " $query_string";
#                    my $readcounts = `$cmd`;
#                    chomp($readcounts) if($readcounts);

                    my $tumor_readcounts = "";
                    if($self->prepend_chr) {
                        $tumor_readcounts = $readcounts_by_position{"chr$chrom\t$chr_start"} if($readcounts_by_position{"chr$chrom\t$chr_start"});
                    }
                    else
                    {
                        $tumor_readcounts = $readcounts_by_position{"$chrom\t$chr_start"} if($readcounts_by_position{"$chrom\t$chr_start"});
                    }


                    ## Parse the results for each allele ##

                    my $tumor_ref_result = $self->read_counts_by_allele($tumor_readcounts, $ref);
                    my $tumor_var_result = $self->read_counts_by_allele($tumor_readcounts, $var);
                    if($tumor_ref_result && $tumor_var_result)
                    {
                        my ($tumor_reads1, $tumor_reads2) = compute_readcount_results($tumor_ref_result, $tumor_var_result);
                        my $tumor_coverage = $tumor_reads1 + $tumor_reads2;
                        my $tumor_var_freq = 0.00;

                        $tumor_var_freq = sprintf("%.3f", $tumor_reads2 / $tumor_coverage) if($tumor_coverage);
                        $tumor_var_freq = ($tumor_var_freq * 100) . '%';
                        print $ofh join("\t", $chrom, $chr_start, $chr_stop, $ref, $var, $tumor_reads1, $tumor_reads2, $tumor_var_freq, $rest_of_line) . "\n";                                

                    }
                    else
                    {
                        $stats{'num_no_readcounts'}++;
                    }
                }
            }
            else
            {
                print "$line\tFailed:UnableToDetermineAllele!\n" if($self->verbose);
                $stats{'num_no_allele'}++;
            }
#       }
    }

    close($input);
    close($ofh);


    Genome::Sys->copy_file($temp_output_file, $self->output_file);

    print $stats{'num_variants'} . " variants\n";
    print $stats{'num_no_allele'} . " failed to determine variant allele\n";
    print $stats{'num_no_readcounts'} . " failed to get readcounts for variant allele\n";
  
    return 1;
}



#############################################################
# Compute_readcount_results - calculate freq, strandedness,
# etc.
#############################################################

sub compute_readcount_results
{
    my ($ref_result, $var_result) = @_;


    ## Parse out the bam-readcounts details for each allele. The fields should be: ##
    #num_reads : avg_mapqual : avg_basequal : avg_semq : reads_plus : reads_minus : avg_clip_read_pos : avg_mmqs : reads_q2 : avg_dist_to_q2 : avgRLclipped : avg_eff_3'_dist
    my ($ref_count, $ref_map_qual, $ref_base_qual, $ref_semq, $ref_plus, $ref_minus, $ref_pos, $ref_subs, $ref_mmqs, $ref_q2_reads, $ref_q2_dist, $ref_avg_rl, $ref_dist_3) = split(/\t/, $ref_result);

    my ($var_count, $var_map_qual, $var_base_qual, $var_semq, $var_plus, $var_minus, $var_pos, $var_subs, $var_mmqs, $var_q2_reads, $var_q2_dist, $var_avg_rl, $var_dist_3) = split(/\t/, $var_result) if($var_result);

    my $ref_strandedness = my $var_strandedness = 0.50;
    $ref_dist_3 = 0.5 if(!$ref_dist_3);

    ## Determine ref strandedness ##

    if(($ref_plus + $ref_minus) > 0) {
        $ref_strandedness = $ref_plus / ($ref_plus + $ref_minus);
        $ref_strandedness = sprintf("%.2f", $ref_strandedness);
    }


    ## If there was no variant result, set some harmless values ##
    if(!$var_result)
    {
        $var_count = $var_plus = $var_minus = 0;
        $var_mmqs = $ref_mmqs;
        $var_map_qual = $ref_map_qual;
        $var_avg_rl = $ref_avg_rl;
        $var_strandedness = $ref_strandedness;
        $var_dist_3 = $ref_dist_3;
    }



    ## Use conservative defaults if we can't get mismatch quality sums ##
    $ref_mmqs = 50 if(!$ref_mmqs);
    $var_mmqs = 0 if(!$var_mmqs);
    my $mismatch_qualsum_diff = $var_mmqs - $ref_mmqs;

    ## Determine map qual diff ##

    my $mapqual_diff = $ref_map_qual - $var_map_qual;


    ## Determine difference in average supporting read length ##

    my $readlen_diff = $ref_avg_rl - $var_avg_rl;

    ## Determine var strandedness ##

    if(($var_plus + $var_minus) > 0)
    {
        $var_strandedness = $var_plus / ($var_plus + $var_minus);
        $var_strandedness = sprintf("%.2f", $var_strandedness);
    }

    ## Determine what we should return ##
    return($ref_count, $var_count); #, $ref_strandedness, $var_strandedness);
}

#############################################################
# Read_Counts_By_Allele - parse out readcount info for an allele
#
#############################################################

sub fails_homopolymer_check {
    (my $self, my $reference, my $min_homopolymer, my $chrom, my $chr_start, my $chr_stop, my $ref, my $var) = @_;

    ## Auto-pass large indels ##

    my $indel_size = length($ref);
    $indel_size = length($var) if(length($var) > $indel_size);

    return(0) if($indel_size > 2);

    ## Build strings of homopolymer bases ##
    my $homoRef = $ref x $min_homopolymer;
    my $homoVar = $var x $min_homopolymer;
#    my $homoA = 'A' x $min_homopolymer;
#    my $homoC = 'C' x $min_homopolymer;
#    my $homoG = 'G' x $min_homopolymer;
#    my $homoT = 'T' x $min_homopolymer;

    ## Build a query string for the homopolymer check ##

    my $query_string = "";

    if($self->prepend_chr) {
        $query_string = "chr" . $chrom . ":" . ($chr_start - $min_homopolymer) . "-" . ($chr_stop + $min_homopolymer);
    } else {
        $query_string = $chrom . ":" . ($chr_start - $min_homopolymer) . "-" . ($chr_stop + $min_homopolymer);
    }

    my $samtools_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
    my $sequence = `$samtools_path faidx $reference $query_string | grep -v \">\"`;
    chomp($sequence);

    if($sequence) {
        if($sequence =~ $homoVar) { #$sequence =~ $homoRef || {
            print join("\t", $chrom, $chr_start, $chr_stop, $ref, $var, "Homopolymer: $sequence") . "\n" if($self->verbose);
            return($sequence);
        }
    }

    return(0);
}


##########################################################################################
# WGS filter for uniform-depth, full-breadth datasets
# Contact: Dave Larson (dlarson@genome.wustl.edu)
##########################################################################################

sub wgs_filter {

}


#############################################################
# Readcount_Program - the path to BAM-readcounts
#
#############################################################

sub readcount_program {
    my $self = shift;
    my $reference = $self->reference;
#    return "bam-readcount0.3 -f $reference";
    return "/gscuser/dlarson/src/genome/bam-readcount/bin/bam-readcount -f $reference";
}


#############################################################
# Read_Counts_By_Allele - parse out readcount info for an allele
#
#############################################################

sub read_counts_by_allele {
    my $self = shift;
    (my $line, my $allele) = @_;

    my @lineContents = split(/\t/, $line);
    my $numContents = @lineContents;

    for(my $colCounter = 5; $colCounter < $numContents; $colCounter++) {
        my $this_allele = $lineContents[$colCounter];
        my @alleleContents = split(/\:/, $this_allele);
        if($alleleContents[0] eq $allele) {
            my $numAlleleContents = @alleleContents;

            return("") if($numAlleleContents < 8);

            my $return_string = "";
            my $return_sum = 0;
            for(my $printCounter = 1; $printCounter < $numAlleleContents; $printCounter++) {
                $return_sum += $alleleContents[$printCounter];
                $return_string .= "\t" if($return_string);
                $return_string .= $alleleContents[$printCounter];
            }

            return($return_string);

#            if($return_sum) {
#                return($return_string);
#            } else {
#                return("");
#            }
        }
    }

    return("");
}

sub iupac_to_base {
    my $self = shift;
    (my $allele1, my $allele2) = @_;

    return($allele2) if($allele2 eq "A" || $allele2 eq "C" || $allele2 eq "G" || $allele2 eq "T");

    if($allele2 eq "M") {
        return("C") if($allele1 eq "A");
        return("A") if($allele1 eq "C");
        return("A");    ## Default for triallelic variant
    } elsif($allele2 eq "R") {
        return("G") if($allele1 eq "A");
        return("A") if($allele1 eq "G");
        return("A");     ## Default for triallelic variant
    } elsif($allele2 eq "W") {
        return("T") if($allele1 eq "A");
        return("A") if($allele1 eq "T");
        return("A");    ## Default for triallelic variant
    } elsif($allele2 eq "S") {
        return("C") if($allele1 eq "G");
        return("G") if($allele1 eq "C");
        return("C");    ## Default for triallelic variant
    } elsif($allele2 eq "Y") {
        return("C") if($allele1 eq "T");
        return("T") if($allele1 eq "C");
        return("C");    ## Default for triallelic variant
    } elsif($allele2 eq "K") {
        return("G") if($allele1 eq "T");
        return("T") if($allele1 eq "G");
        return("G");    ## Default for triallelic variant
    }

    return($allele2);
}

1;
