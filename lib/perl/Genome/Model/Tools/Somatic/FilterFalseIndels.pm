package Genome::Model::Tools::Somatic::FilterFalseIndels;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Readonly;
use Genome::Info::IUB;

class Genome::Model::Tools::Somatic::FilterFalseIndels {
    is => 'Command',
    has => [
	
	## INPUT/OUTPUT OPTIONS ##
       'bam_file' => {
            type => 'String',
            doc => 'BAM file in which to examine reads (usually tumor BAM)',
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
       'filtered_file' => {
           type => 'String',
           is_input => 1,
	   is_optional => 1,
           is_output => 1,
           doc => 'Filename for variants that fail filter (optional)',
       },
        'reference' => {
            type => 'String',            
            is_optional => 1,
            is_input => 1,
            doc => 'Reference sequence to use. Build 37: /gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa  Build 36: /gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa',
            example_values => ['/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa', '/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa'],
        },  
       'min_homopolymer' => {
            type => 'String',
            default => '4',
            is_optional => 1,
            is_input => 1,
            doc => 'Minimum length of a homopolymer to auto-filter small (1-2 bp) indels',
       },       
       ## CAPTURE FILTER OPTIONS ##
       'min_strandedness' => {
            type => 'String',
            default => '0.01',
            is_optional => 1,
            is_input => 1,
            doc => 'Minimum representation of variant allele on each strand',
       },
       'min_good_coverage' => {
            type => 'String',
            default => '20',
            is_optional => 1,
            is_input => 1,
            doc => 'Minimum site coverage to apply var_freq, var_count, and strandedness filters',
       },
       'min_var_freq' => {
            type => 'String',
            default => '0.05',
            is_optional => 1,
            is_input => 1,
            doc => 'Minimum variant allele frequency',
       },
       'min_var_count' => {
            type => 'String',
            default => '2',
            is_optional => 1,
            is_input => 1,
            doc => 'Minimum number of variant-supporting reads [2]',
       },
       'min_read_pos' => {
            type => 'String',
            default => '0.10',
            is_optional => 1,
            is_input => 1,
            doc => 'Minimum average relative distance from start/end of read (start/end of read filter) [0.10]',
       },
       'max_mm_qualsum_diff' => {
            type => 'String',
            default => '50',
            is_optional => 1,
            is_input => 1,
            doc => 'Maximum difference of mismatch quality sum between variant and reference reads (paralog filter)',
       },
       'max_var_mmqs' => {
            type => 'String',
            is_optional => 1,
            is_input => 1,
            doc => 'Maximum mismatch quality sum for variant-supporting reads (paralog filter) [opt:100]',
       },
       'max_mapqual_diff' => {
            type => 'String',
            default => '30',
            is_optional => 1,
            is_input => 1,
            doc => 'Maximum difference of mapping quality between variant and reference reads [30]',
       },
       'max_readlen_diff' => {
            type => 'String',
            default => '15',
            is_optional => 1,
            is_input => 1,
            doc => 'Maximum difference of average supporting read length between variant and reference reads [15]',
       },
       'min_var_readlen' => {
            type => 'String',
            is_optional => 1,
            is_input => 1,
            doc => 'Minimum average aligned read length of variant-supporting reads [recommended: 75]',
       },
       'min_var_dist_3' => {
            type => 'String',
            default => '0.20',
            is_optional => 1,
            is_input => 1,
            doc => 'Minimum average distance to effective 3prime end of read (real end or Q2) for variant-supporting reads',
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
            default_value => 'rusage[mem=4000] select[type==LINUX64] span[hosts=1]',
       },
       lsf_queue => {
            is_param => 1,
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
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
    ]
};

sub help_brief {
    return "This module uses detailed readcount information from bam-readcounts to filter likely false positives";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    gmt somatic filter-false-indels --variant-file somatic.indels --tumor-bam tumor.bam --output-file somatic.indels.FilterFalseIndels 
EOS
}

sub help_detail {                           
    return <<EOS 
This module uses detailed readcount information from bam-readcounts to filter likely false positives
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
        $self->error_message("Tumor bam must be indexed");
        die;
    }

    my $reference = $self->reference;
    my $min_homopolymer = $self->min_homopolymer;

    ## Determine the strandedness and read position thresholds ##
    
    my $min_read_pos = $self->min_read_pos;
    my $max_read_pos = 1 - $min_read_pos;
    my $min_var_freq = $self->min_var_freq;
    my $min_var_count = $self->min_var_count;
    
    my $min_strandedness = $self->min_strandedness;
    my $max_strandedness = 1 - $min_strandedness;

    my $max_mm_qualsum_diff = $self->max_mm_qualsum_diff;
    my $max_mapqual_diff = $self->max_mapqual_diff;
    my $max_readlen_diff = $self->max_readlen_diff;
    my $min_var_dist_3 = $self->min_var_dist_3;

    ## Reset counters ##
    
    my %stats = ();
    $stats{'num_variants'}  = $stats{'num_no_readcounts'} = $stats{'num_pass_filter'} = $stats{'num_no_allele'} = 0;
    $stats{'num_fail_varcount'} = $stats{'num_fail_varfreq'} = $stats{'num_fail_strand'} = $stats{'num_fail_pos'} = $stats{'num_fail_mmqs'} = $stats{'num_fail_mapqual'} = $stats{'num_fail_readlen'} = $stats{'num_fail_dist3'} = 0;
    $stats{'num_MT_sites_autopassed'} = 0;

    ## Open the output file ##
    
    my $ofh = IO::File->new($self->output_file, "w");
    unless($ofh) {
        $self->error_message("Unable to open " . $self->output_file . " for writing. $!");
        die;
    }

    my $filtered_file = $self->output_file . ".removed";
    $filtered_file = $self->filtered_file if($self->filtered_file);

    ## Open the variants file ##

    my $input = new FileHandle ($self->variant_file);

    unless($input) {
        $self->error_message("Unable to open " . $self->variant_file . ". $!");
        die;
    }


    ## Build temp file for positions where readcounts are needed ##

    my ($tfh,$temp_path) = Genome::Sys->create_temp_file;
    unless($tfh) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }
    $temp_path =~ s/\:/\\\:/g;

    ## Print each line to file, prepending chromosome if necessary ##
    print "Printing variants to temp file...\n";
    while(my $line = $input->getline) {
	chomp $line;
	my ($chr, $start, $stop, $ref, $var) = split /\t/, $line;
	if ($self->prepend_chr) {
	    $chr = "chr$chr";
	    $chr =~ s/MT$/M/;
	};

	## Determine indel size ##
	
	my $indel_size = length($ref);
	$indel_size = length($var) if(length($var) > $indel_size);
	
	## Based on start/stop, expand out the region to use for readcounts ##
	
	for(my $position = ($start - $indel_size - 2); $position <= ($stop + $indel_size + 2); $position++)
	{
	    print $tfh "$chr\t$position\t$position\n";
	}

    }
    $tfh->close;

    close($input);

    ## Run BAM readcounts in batch mode to get read counts for all positions in file ##

    print "Running BAM Readcounts...\n";
#    my $cmd = readcount_program() . " -b 15 " . $self->bam_file . " -l $temp_path";
    my $cmd = $self->readcount_program . " -b 1 " . $self->bam_file . " -l $temp_path";
    my $readcount_path = $self->output_file . ".readcounts";

    $cmd .= "> $readcount_path 2> /dev/null";
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->bam_file],
        output_files => [$readcount_path],
    );

    ## Load the results of the readcounts ##
    my %readcounts_by_position = ();

    my $readcount_fh = Genome::Sys->open_file_for_reading($readcount_path);
    while(my $rc_line = $readcount_fh->getline)
    {
     chomp($rc_line);
	(my $chrom, my $pos) = split(/\t/, $rc_line);
	$readcounts_by_position{"$chrom\t$pos"} = $rc_line;
   }

    print "Readcounts loaded\n";

    

    ## Open the filtered output file ##
    
    my $ffh = IO::File->new($filtered_file, "w") if($filtered_file);


    ## Reopen file for parsing ##

    $input = new FileHandle ($self->variant_file);
    

    ## Parse the variants file ##

    my $lineCounter = 0;
    
    while (<$input>)
    {
            chomp;
            my $line = $_;
            $lineCounter++;

            $stats{'num_variants'}++;
            
#            if($lineCounter <= 10)
 #           {
                (my $chrom, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
                
		$ref = uc($ref);
		$var = uc($var);
		
		## Determine indel size ##
		
		my $indel_size = length($ref);
		$indel_size = length($var) if(length($var) > $indel_size);


		## Proceed if indel size is valid ##
    
                if($indel_size > 0)
                {


		    if($indel_size <= 2 && $min_homopolymer && fails_homopolymer_check($self, $reference, $min_homopolymer, $chrom, $chr_start, $chr_stop, $ref, $var))
		    {
			$stats{'num_fail_homopolymer'}++;
			print $ffh join("\t", $line, "Homopolymer") . "\n" if($self->filtered_file);	
		    }
		    ## Skip MT chromosome sites,w hich almost always pass ##
	
		    elsif($chrom eq "MT" || $chrom eq "chrMT")
		    {
			## Auto-pass it to increase performance ##
			$stats{'num_MT_sites_autopassed'}++;
			print $ofh "$line\tMT_autopass\n";
		    }
		    else
		    {
			## Obtain readcounts all around this indel ##
			my $readcounts = "";
			for(my $position = ($chr_start - $indel_size - 2); $position <= ($chr_stop + $indel_size + 2); $position++)
			{
			    $readcounts .= "\n" if($readcounts);
			    
			    ## Append this readcount to the line of readcounts ##
			    
			    if($self->prepend_chr)
			    {
				    $readcounts .= "$position\t" . $readcounts_by_position{"chr$chrom\t$position"} if($readcounts_by_position{"chr$chrom\t$position"});			    
			    }
			    else
			    {
				    $readcounts .= "$position\t" . $readcounts_by_position{"$chrom\t$position"} if($readcounts_by_position{"$chrom\t$position"});			    
			    }			    
			}

		

			my ($ref_result, $var_result) = read_counts_for_indel($readcounts, $chr_start, $chr_stop, $ref, $var);

			if($ref_result && $var_result)
			{
			    $stats{'num_got_readcounts'}++;
				## Parse out the bam-readcounts details for each allele. The fields should be: ##
				# num_reads : avg_mapqual : avg_basequal : avg_semq : reads_plus : reads_minus : avg_clip_read_pos : avg_mmqs : reads_q2 : avg_dist_to_q2 : avgRLclipped : avg_eff_3'_dist
				my ($ref_count, $ref_map_qual, $ref_base_qual, $ref_semq, $ref_plus, $ref_minus, $ref_pos, $ref_subs, $ref_mmqs, $ref_q2_reads, $ref_q2_dist, $ref_avg_rl, $ref_dist_3) = split(/\t/, $ref_result);
				my ($var_allele, $var_count, $var_map_qual, $var_base_qual, $var_semq, $var_plus, $var_minus, $var_pos, $var_subs, $var_mmqs, $var_q2_reads, $var_q2_dist, $var_avg_rl, $var_dist_3) = split(/\t/, $var_result);
    
				my $ref_strandedness = my $var_strandedness = 0.50;

				## Use conservative defaults if we can't get mismatch quality sums ##
				$ref_mmqs = 50 if(!$ref_mmqs);
				$var_mmqs = 0 if(!$var_mmqs);
				my $mismatch_qualsum_diff = $var_mmqs - $ref_mmqs;

				## Determine map qual diff ##
				
				my $mapqual_diff = $ref_map_qual - $var_map_qual;


				## Determine difference in average supporting read length ##
				
				my $readlen_diff = $ref_avg_rl - $var_avg_rl;

    
				## Determine ref strandedness ##
				
				if(($ref_plus + $ref_minus) > 0)
				{
				    $ref_strandedness = $ref_plus / ($ref_plus + $ref_minus);
				    $ref_strandedness = sprintf("%.2f", $ref_strandedness);
				}
    
				## Determine var strandedness ##
    
				if(($var_plus + $var_minus) > 0)
				{
				    $var_strandedness = $var_plus / ($var_plus + $var_minus);
				    $var_strandedness = sprintf("%.2f", $var_strandedness);
				}
    
				my $FilterResult = "";
    
				if($var_count && ($var_plus + $var_minus))
				{
				    ## We must obtain variant read counts to proceed ##

				    my $var_freq = $var_count / ($ref_count + $var_count);

				    my $readcount_info = join("\t", $ref_count, $var_count, $var_freq, $ref_pos, $var_pos, $ref_avg_rl, $var_avg_rl, $ref_strandedness, $var_strandedness, $ref_mmqs, $var_mmqs, $mismatch_qualsum_diff);
				    
				    ## FAILURE 1: READ POSITION ##
				    
				    if(($var_pos < $min_read_pos))# || $var_pos > $max_read_pos))
				    {
					$FilterResult = "ReadPos<$min_read_pos";
					$stats{'num_fail_pos'}++;
				    }
				    
				    ## FAILURE 2: Variant is strand-specific but reference is NOT strand-specific ##
				    
				    elsif(($var_strandedness < $min_strandedness || $var_strandedness > $max_strandedness) && ($ref_strandedness >= $min_strandedness && $ref_strandedness <= $max_strandedness))
				    {
					$FilterResult = "Strandedness: Ref=$ref_strandedness Var=$var_strandedness";
					$stats{'num_fail_strand'}++;
				    }

				    ## FAILURE 2b : Variant allele count does not meet minimum ##    

				    elsif($var_count < $min_var_count)
				    {
					$FilterResult = "VarCount:$var_count<$min_var_count";
					$stats{'num_fail_varcount'}++;					
				    }

				    ## FAILURE 2c : Variant allele frequency does not meet minimum ##    

				    elsif($var_freq < $min_var_freq)
				    {
					$FilterResult = "VarFreq:$var_freq<$min_var_freq";
					$stats{'num_fail_varfreq'}++;					
				    }
				    
				    ## FAILURE 3: Paralog filter for sites where variant allele mismatch-quality-sum is significantly higher than reference allele mmqs
				    
				    elsif($mismatch_qualsum_diff> $max_mm_qualsum_diff)
				    {
					$FilterResult = "MMQSdiff:$var_mmqs-$ref_mmqs=$mismatch_qualsum_diff>$max_mm_qualsum_diff";
					$stats{'num_fail_mmqs'}++;					
				    }
				    
				    ## FAILURE 4: Mapping quality difference exceeds allowable maximum ##
				    elsif($mapqual_diff > $max_mapqual_diff)
				    {
					$FilterResult = "MapQual:$ref_map_qual-$var_map_qual=$mapqual_diff>$max_mapqual_diff";
					$stats{'num_fail_mapqual'}++;
				    }
				    
				    ## FAILURE 5: Read length difference exceeds allowable maximum ##
				    elsif($readlen_diff > $max_readlen_diff)
				    {
					$FilterResult = "ReadLen:$ref_avg_rl-$var_avg_rl=$readlen_diff>$max_readlen_diff";
					$stats{'num_fail_readlen'}++;
				    }
				    ## FAILURE 6: Var read len below minimum ##
				    elsif($self->min_var_readlen && $var_avg_rl < $self->min_var_readlen)
				    {
					$FilterResult = "VarReadLen:$var_avg_rl\n" if ($self->verbose);
					$stats{'num_fail_var_readlen'}++;
				    }
				    elsif($self->max_var_mmqs && $var_mmqs > $self->max_var_mmqs)
				    {
					$FilterResult = "VarMMQS:$var_mmqs\n" if ($self->verbose);
					$stats{'num_fail_var_mmqs'}++;					
				    }
				    elsif($ref_mmqs < 30 && $var_mmqs > 50)
				    {
					$FilterResult = "MMQSthreshold:ref=$ref_mmqs,var=$var_mmqs\n" if ($self->verbose);
					$stats{'num_fail_mmqs_threshold'}++;										
				    }
				    ## SUCCESS: Pass Filter ##				
				    else
				    {					
					$FilterResult = "Pass";					
				    }
				    
				    ## Print filter status and result ##
				    
				    if($FilterResult eq "Pass")
				    {
					$stats{'num_pass_filter'}++;
					print $ofh join("\t", $line, $readcount_info) . "\n";
				    }
				    else
				    {
					$stats{'num_fail_filter'}++;
					print $ffh join("\t", $line, $readcount_info, $FilterResult) . "\n" if($self->filtered_file);					
				    }
					
				    print join("\t", $line, $readcount_info) . "\n" if($self->verbose);

				}
				else
				{
				    $stats{'num_no_readcounts'}++;
				}
			}
			else
			{
			    $stats{'num_no_readcounts_autopassed'}++;
			    
			    if($self->verbose)
			    {
				print join("\t", $chrom, $chr_start, $chr_stop, $ref, $var, "NoReadcounts_autopassed") . "\n";
				my $no_rc_alleles = "";
				my @rc_lines = split(/\n/, $readcounts);
				foreach my $rc_line (@rc_lines)
				{
				    my @rcLineContents = split(/\t/, $rc_line);
				    my $position = $rcLineContents[2];
				    my $numRClineContents = @rcLineContents;
				    for(my $colCounter = 5; $colCounter < $numRClineContents; $colCounter++)
				    {
					(my $allele) = split(/\:/, $rcLineContents[$colCounter]);
					$no_rc_alleles .= "\t$position\t$allele\n" if(length($allele) > 1);
				    }
				}
				print "$no_rc_alleles\n" if($no_rc_alleles);
			    }
			    
			    print $ofh join("\t", $line, "\tNoReadcount_autopass") . "\n";
			}			
		    }

                }
		else
		{
		    print $ffh "$line\tno_allele\n" if($self->filtered_file);
		    print "$line\tFailure: Unable to determine allele!\n" if($self->verbose);        
		    $stats{'num_no_allele'}++;
		}
#            }
    }
    $ffh->close;
    $ofh->close;
    close($input);

    print $stats{'num_variants'} . " variants\n";
    print $stats{'num_fail_homopolymer'} . " failed the homopolymer check\n";
    print $stats{'num_MT_sites_autopassed'} . " MT sites were auto-passed\n";
    print $stats{'num_no_readcounts_autopassed'} . " failed to get readcounts and were auto-passed\n";
    print $stats{'num_got_readcounts'} . " other variants got readcounts\n";
    print $stats{'num_pass_filter'} . " passed the strand filter\n";
    print $stats{'num_fail_filter'} . " failed the strand filter\n";
    print "\t" . $stats{'num_no_allele'} . " failed to determine variant allele\n";
    print "\t" . $stats{'num_no_readcounts'} . " failed to get readcounts for variant allele\n";
    print "\t" . $stats{'num_fail_pos'} . " had read position < $min_read_pos\n";
    print "\t" . $stats{'num_fail_strand'} . " had strandedness < $min_strandedness\n";
    print "\t" . $stats{'num_fail_varcount'} . " had var_count < $min_var_count\n";
    print "\t" . $stats{'num_fail_varfreq'} . " had var_freq < $min_var_freq\n";

    print "\t" . $stats{'num_fail_mmqs'} . " had mismatch qualsum difference > $max_mm_qualsum_diff\n";
    print "\t" . $stats{'num_fail_mapqual'} . " had mapping quality difference > $max_mapqual_diff\n";
    print "\t" . $stats{'num_fail_readlen'} . " had read length difference > $max_readlen_diff\n";	
    print "\t" . $stats{'num_fail_dist3'} . " had var_distance_to_3' < $min_var_dist_3\n";

    ## Print new optional filters if they applied ##
    print "\t" . $stats{'num_fail_mmqs_threshold'} . " had ref MMQS < 30 but var MMQS > 50\n" if($stats{'num_fail_mmqs_threshold'});
    print "\t" . $stats{'num_fail_var_mmqs'} . " had var MMQS > " . $self->max_var_mmqs . "\n" if($stats{'num_fail_var_mmqs'});
    print "\t" . $stats{'num_fail_var_readlen'} . " had var readlen < " . $self->min_var_readlen . "\n" if($stats{'num_fail_var_readlen'});
    
    return 1;
}


#############################################################
# Read_Counts_By_Allele - parse out readcount info for an allele
#
#############################################################

sub fails_homopolymer_check
{
    (my $self, my $reference, my $min_homopolymer, my $chrom, my $chr_start, my $chr_stop, my $ref, my $var) = @_;
    
    ## Auto-pass large indels ##
    
    my $indel_size = length($ref);
    $indel_size = length($var) if(length($var) > $indel_size);
    
    return(0) if($indel_size > 2);
    
    ## Build strings of homopolymer bases ##
    my $homoA = 'A' x $min_homopolymer;
    my $homoC = 'C' x $min_homopolymer;
    my $homoG = 'G' x $min_homopolymer;
    my $homoT = 'T' x $min_homopolymer;
    
    ## Build a query string for the homopolymer check ##
    
    my $query_string = "";
    
    if($self->prepend_chr)
    {
	$query_string = "chr" . $chrom . ":" . ($chr_start - $min_homopolymer) . "-" . ($chr_stop + $min_homopolymer);
    }
    else
    {
	$query_string = $chrom . ":" . ($chr_start - $min_homopolymer) . "-" . ($chr_stop + $min_homopolymer);
    }
    
    my $sequence = `samtools faidx $reference $query_string | grep -v \">\"`;
    chomp($sequence);
    
    if($sequence)
    {
	if($sequence =~ $homoA || $sequence =~ $homoC || $sequence =~ $homoG || $sequence =~ $homoT)
	{
	    return($sequence);
	}
    }
    
    return(0);
}


#############################################################
# Read_Counts_By_Allele - parse out readcount info for an allele
#
#############################################################

sub read_counts_for_indel
{
    my ($readcount_lines, $chr_start, $chr_stop, $ref, $var) = @_;

    ## Determine indel type and size ##
    
    my $indel_type = my $indel_size = "";
    
    if($ref eq "0" || $ref eq "-")
    {
	$indel_type = "INS";
	$indel_size = length($var);
    }
    else
    {
	$indel_type = "DEL";
	$indel_size = length($ref);
    }

    if(!$indel_type || !$indel_size)
    {
	warn "Unable to determine indel type/size for $chr_start $chr_stop $ref $var\n";
	return();
    }
    
    my @matching_results = ();
    my $num_matching_results = 0;
    my @lines = split(/\n/, $readcount_lines);
    
    my %readcounts_by_position = ();
    
    foreach my $line (@lines)
    {
	my @lineContents = split(/\t/, $line);
	my $numContents = @lineContents;
	my $position = $lineContents[0];

	## Build String of the Readcounts ##

	my $readcounts = "";
	for(my $colCounter = 1; $colCounter < $numContents; $colCounter++)
	{
	    $readcounts .= "\t" if($readcounts);
	    $readcounts .= $lineContents[$colCounter];
	}

	$readcounts_by_position{$position} = $readcounts;

	my $position_result = read_counts_by_type_size($readcounts, $indel_type, $indel_size);

	if($position_result)
	{
	    ## Calculate min distance to indel ##
	    my $min_distance = get_min_distance($position, $chr_start, $chr_stop);
	    ## Split multiple possible alleles that were found at this pos ##
	    
	    my @position_result = split(/\n/, $position_result);
	    foreach my $allele_result(@position_result)
	    {
		$matching_results[$num_matching_results] = $position . "\t" . $min_distance . "\t" . $allele_result;
		$num_matching_results++;
	    }

	}
	else
	{

	}
    }


    my @top_results = sort byReadcountThenDistance @matching_results;

    if($top_results[0])
    {
	## Parse out the readcount results for the top matching indel ##
	my @lineContents = split(/\t/, $top_results[0]);
	my $numContents = @lineContents;
	my $position = $lineContents[0];
	my $distance = $lineContents[1];
	my $allele = $lineContents[2];
	my $var_result = "";
	for(my $colCounter = 2; $colCounter < $numContents; $colCounter++)	## Note, allele included ##
	{
	    $var_result .= "\t" if($var_result);
	    $var_result .= $lineContents[$colCounter];	   
	}

	## At the same position, get the readcounts for the reference base ##
	
	my $ref_result = read_counts_for_reference($readcounts_by_position{$position});
	
	return($ref_result, $var_result);
    }

    ## Return ref and var readcounts ##
    return(0);
    
}




sub byReadcountThenDistance
{
    my ($pos_a, $dist_a, $allele_a, $count_a) = split(/\t/, $a);
    my ($pos_b, $dist_b, $allele_b, $count_b) = split(/\t/, $b);
    
    $count_b <=> $count_a
    or
    $dist_a <=> $dist_b;
}


#############################################################
# Read_Counts_By_Allele - parse out readcount info for an allele
#
#############################################################

sub get_min_distance
{
    my ($position, $start, $stop) = @_;
    
    my $dist1 = abs($start - $position) + 1;
    my $dist2 = abs($stop - $position) + 1;
    
    return($dist2) if($dist2 < $dist1);
    return($dist1);
}



#############################################################
# Read_Counts_By_Allele - parse out readcount info for an allele
#
#############################################################

sub read_counts_by_type_size
{
	(my $line, my $indel_type, my $indel_size) = @_;
	
	my $result = "";
	
	my @lineContents = split(/\t/, $line);
	
	my $numContents = @lineContents;
	
	for(my $colCounter = 5; $colCounter < $numContents; $colCounter++)
	{
		my $this_allele_contents = $lineContents[$colCounter];
		my @alleleContents = split(/\:/, $this_allele_contents);
		my $numAlleleContents = @alleleContents;
		my $this_allele = $alleleContents[0];

		if(length($this_allele) > 1)
		{
			my $this_indel_type = my $this_indel_size = 0;
			
			## Determine indel type and size ##			
			if(substr($this_allele, 0, 1) eq '+')
			{
				$this_indel_type = "INS";
				$this_indel_size = length($this_allele) - 1;
			}
			elsif(substr($this_allele, 0, 1) eq '-')
			{
				$this_indel_type = "DEL";
				$this_indel_size = length($this_allele) - 1;				
			}

			## Try to match indel type and size ##

			if(match_type_size($indel_type, $indel_size, $this_indel_type, $this_indel_size))
			{
				my $return_string = "$this_allele";
				my $return_sum = 0;
				for(my $printCounter = 1; $printCounter < $numAlleleContents; $printCounter++)
				{
					$return_string .= "\t" if($return_string);
					$return_string .= $alleleContents[$printCounter];
				}
				
				$result .= "\n" if($result);
				$result .= $return_string;
			}
			else
			{
#				print "$this_allele had $this_indel_type $this_indel_size and did not match $indel_type $indel_size\n";
			}

		}

	}
	
	return($result);
}



#############################################################
# Read_Counts_By_Allele - parse out readcount info for an allele
#
#############################################################

sub match_type_size
{
    my ($indel_type, $indel_size, $test_indel_type, $test_indel_size) = @_;
    
    return(1) if($indel_type eq $test_indel_type && $indel_size == $test_indel_size);

    if($indel_type eq $test_indel_type)
    {
	if($indel_size > 4)
	{
	    my $indel_size_diff = abs($indel_size - $test_indel_size);
	    my $diff_fraction = $indel_size_diff / $indel_size;
	    if($diff_fraction <= 0.20)
	    {
		return(1);
	    }
	}
    }

    return(0);
}


#############################################################
# Read_Counts_By_Allele - parse out readcount info for an allele
#
#############################################################

sub read_counts_for_reference
{
	(my $line) = @_;
	
	my @lineContents = split(/\t/, $line);
	my $numContents = @lineContents;
	
	## parse out the reference allele ##
	my $allele = $lineContents[2];
	
	for(my $colCounter = 5; $colCounter < $numContents; $colCounter++)
	{
		my $this_allele = $lineContents[$colCounter];
		my @alleleContents = split(/\:/, $this_allele);
		if($alleleContents[0] eq $allele)
		{
			my $numAlleleContents = @alleleContents;
			
			return("") if($numAlleleContents < 8);
			
			my $return_string = "";
			my $return_sum = 0;
			for(my $printCounter = 1; $printCounter < $numAlleleContents; $printCounter++)
			{
				$return_sum += $alleleContents[$printCounter];
				$return_string .= "\t" if($return_string);
				$return_string .= $alleleContents[$printCounter];
			}
			
                        return($return_string);
                        
#			if($return_sum)
#			{
#				return($return_string);
#			}
#			else
#			{
#				return("");
#			}
		}
	}
	
	return("");
}






#############################################################
# Readcount_Program - the path to BAM-readcounts
#
#############################################################

sub readcount_program {
    my $self = shift;
    return "/usr/bin/bam-readcount0.3 -f ".$self->reference;
}


1;
