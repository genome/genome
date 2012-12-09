package VarScan::ParseCrossmatch;

use warnings;
use strict;

=head1 NAME

VarScan::ParseCrossmatch - Parses output files from the cross_match aligner

=head1 VERSION

    Version 1.03

=cut

our $VERSION = '1.03';


## Define semi-global variables ##
my %stats = ();
my $verbose;

## Define shared variables ##

my $reference_db = my $fasta_db = my $quality_db;
my $ref_dir, my $fasta_file, my $quality_file, my $min_identity, my $min_align_score;
my $primer_trim, my $default_qual_score, my $min_qual_score;
my $output_dir, my $sample, my $output_alignments, my $output_snps, my $output_indels;
my $min_indel_size, my $max_indel_size;


=head1 SYNOPSIS

    This module

=head1 FUNCTIONS

=cut

################################################################################

=head2	parse_alignments - parses the alignments file

=cut
################################################################################

sub parse_alignments
{
    (my $alignments_file) = @_;	#, my $output_file

    ## Reset statistics ##
    $stats{'num_alignments'} = $stats{'num_snps'} = $stats{'INSERTION'} = $stats{'DELETION'} = $stats{'indels_filtered'} = 0;

    ## Get variables that we'll need ##
    
    $fasta_file = $VarScan::fasta_file;
    $quality_file = $VarScan::quality_file;
    $ref_dir = $VarScan::ref_dir;
    $min_identity = $VarScan::min_identity;
    $min_align_score = $VarScan::min_align_score;
    $primer_trim = $VarScan::primer_trim;
    $default_qual_score = $VarScan::default_qual_score;
    $min_qual_score = $VarScan::min_qual_score;
    $min_indel_size = $VarScan::min_indel_size;
    $max_indel_size = $VarScan::max_indel_size;

    $output_dir = $VarScan::output_dir;
    $sample = $VarScan::sample;
    
    $output_alignments = $VarScan::output_alignments;
    $output_snps = $VarScan::output_snps;
    $output_indels = $VarScan::output_indels;


    ## Get fasta file if set up to do so ##
    
    if($fasta_file && !(-e $fasta_file))
    {
	print "\n***ERROR*** FASTA file not found!\n\n";
    }
    
    ## Get quality file if set up to do so ##
    
    if($quality_file && !(-e $quality_file))
    {
	print "\n***ERROR*** Quality file not found!\n\n";
    }

    ############################
    ### BUILD BIO::DB::FASTA ###
    ############################
    use Cwd;
    
    my $cwd = cwd();
#    my $fasta_db = my $quality_db = my $reference_db;
    
    ## Reference Sequence Dir ##
    
    if($ref_dir && -d $ref_dir)
    {
	## Build the Bio::DB::Fasta index ##
	 $reference_db = VarScan::SequenceFile::build_bio_db_fasta($ref_dir);
    }

    ## Build Bio::DB::Fasta Objects ##


    if($fasta_file && -e $fasta_file)
    {
	## Build the fasta index ##
    
	mkdir("$output_dir/fasta_dir") if(!(-d "$output_dir/fasta_dir"));    
	unlink("$output_dir/fasta_dir/$sample.fasta");
	symlink("$cwd/$fasta_file", "$output_dir/fasta_dir/$sample.fasta");
	$fasta_db = VarScan::SequenceFile::build_bio_db_fasta("$output_dir/fasta_dir/");
    }

    if($quality_file && -e $quality_file)
    {
	print "Building index for quality scores...\n";

	## Build the quality index ##
    
	mkdir("$output_dir/quality_dir") if(!(-d "$output_dir/quality_dir"));    
    
	## Convert the file ##
    
	my $quality_outfile = $output_dir . "/quality_dir/$sample.qual.fa";
	if(!(-e $quality_outfile))
	{
	    VarScan::SequenceFile::quality_to_fasta($quality_file, $quality_outfile);
	}
	
    #    if(!(-e "$output_dir/quality_dir/directory.index"))
	
	## Build the Bio::DB::Fasta index ##
	
	$quality_db = VarScan::SequenceFile::build_bio_db_fasta("$output_dir/quality_dir/");
    }




    print "Minimum alignment score: $min_align_score\n" if($min_align_score);
    print "Minimum identity: $min_identity\n" if($min_identity);

    ## Parse out the best alignments ##

    my %read_alignments = get_best_alignments($alignments_file);

    print "$stats{'num_alignments'} alignments were parsed\n";
    print "$stats{'num_reads_aligned'} reads had at least one alignment\n";
    print "$stats{'reads_with_best_align'} reads had a single best alignment\n";


    ## SNP PARSING ##

    ## Open output file if desired ##

    if($output_alignments)
    {
	    open(ALIGNMENTS, ">$output_alignments") or die "Can't open outfile: $!\n";
#	    print ALIGNMENTS "chrom\tposition\tref_allele\tvar_allele\tread_name\tread_pos\talign_strand\tqual_score\n";
    }

    if($output_snps)
    {
	    open(OUTFILE, ">$output_snps") or die "Can't open outfile: $!\n";
	    print OUTFILE "chrom\tposition\tref_allele\tvar_allele\tread_name\tread_pos\talign_strand\tqual_score\n";
    }

    if($output_indels)
    {
	    open(INDELS, ">$output_indels") or die "Can't open outfile: $!\n";
	    print INDELS "chrom\tchr_start\tchr_stop\tindel_type\tindel_size\tread_name\tread_start\tread_stop\talign_strand\tconf_score\tindel_in_context\n";
    }

    ## Go through best read alignments ##
    foreach my $read_name (keys %read_alignments)
    {
	my $best_read_alignment = $read_alignments{$read_name};
	if($output_alignments)
	{
#	    (my $align_score) = split(/\t/, $best_read_alignment);
	    my $unscored_alignment = $best_read_alignment;
#            $unscored_alignment = "ALIGNMENT\t" . $unscored_alignment;
#	    $unscored_alignment =~ s/$align_score\t//;
	    print ALIGNMENTS "$unscored_alignment\n" 
	}

	if($output_snps)
	{
	    my $read_snps = screen_for_snps($best_read_alignment, $quality_db);
	    print OUTFILE "$read_snps" if($output_snps && $read_snps);
	}

	if($output_indels)
	{
#	    my $read_indels = screen_for_indels($best_read_alignment);
#	    print INDELS "$read_indels" if($output_indels && $read_indels);
	}
    }


    ## Close the output file ##

    close(OUTFILE) if($output_snps);
    close(ALIGNMENTS) if($output_alignments);

    ## Print summary statistics ##

    $stats{'INSERTION'} = 0 if(!$stats{'INSERTION'});
    $stats{'DELETION'} = 0 if(!$stats{'DELETION'});
    $stats{'num_snps'} = 0 if(!$stats{'num_snps'});

    print "$stats{'num_snps'} substitution (SNP) events detected\n";
    print "$stats{'indels_filtered'} indel events removed by heuristic filter\n";
    print "$stats{'INSERTION'} insertion events detected\n";
    print "$stats{'DELETION'} deletion events detected\n";
    
    return(0);
}



##########################
# Sorting function 
##########################

sub by_align_score
{
	my @temp = split(/\s+/, $a);
	my $score_a = $temp[1];
	@temp = split(/\s+/, $b);
	my $score_b = $temp[1];		
	$score_b <=> $score_a;
}


################################################################################

=head2	get_best_alignments - retrieve the single best alignment by score

=cut
################################################################################

sub get_best_alignments
{
    my $alignments_file = shift(@_);
    
    ## Parse and score the BLAT alignments ##   
    my %read_alignments = ();
    my %best_align_scores = ();
    my %best_alignments = ();

    # Alignment variables #
    my $read_name = my $read_start = my $read_stop = my $ref_name = my $ref_start = my $ref_stop = "";
    my $align_score = my $align_strand = "";
    my $current_alignment = "";

    ## Pre-parse the file to get the best alignments for each read ##
    my %alignment_discreps = ();

    my $input = new FileHandle ($alignments_file);
    my $lineCounter = 0;

    while (<$input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;

            if($lineCounter >= 1 && $line && substr($line, 0, 9) eq "ALIGNMENT")
            {
                    $stats{'num_alignments'}++;
            
                    my @lineContents = split(/\s+/, $line);
                    my $numContents = @lineContents;
                    $align_score = $lineContents[1];
                    my $pct_snps = $lineContents[2];
                    $read_name = $lineContents[5];

                    if($lineContents[6] > $lineContents[7])
                    {
                            $read_start = $lineContents[7];
                            $read_stop = $lineContents[6];			
                    }
                    else
                    {
                            $read_start = $lineContents[6];
                            $read_stop = $lineContents[7];			
                    }

                    if($lineContents[9] eq 'C')
                    {
                            $align_strand = '-';
                            
                            $ref_name = $lineContents[10];
                            $ref_start = $lineContents[13];
                            $ref_stop = $lineContents[12];
                    }
                    else
                    {
                            $align_strand = '+';
                            $ref_name = $lineContents[9];
                            $ref_start = $lineContents[10];
                            $ref_stop = $lineContents[11];
                    }

                    ## Calculate identity ##
                    
                    my $alignment_identity = 100 - $pct_snps;

                    if($align_score >= $min_align_score)
                    {
                        if($alignment_identity && $alignment_identity >= $min_identity)
                        {
        #                    $current_alignment = "$align_score\t$ref_name\t$ref_start\t$ref_stop\t$align_strand\t$read_name\t$read_start\t$read_stop";
                            $current_alignment = $line;
                            $read_alignments{$read_name} .= "$current_alignment\n";
                        }
                    }
            }
            elsif($lineCounter >= 1 && $line && substr($line, 0, 11) eq "DISCREPANCY")
            {
#                    my @lineContents = split(/\s+/, $line);
                my $discrep_line = $line;
                $discrep_line =~ s/\s+/\t/g;
                $alignment_discreps{$current_alignment} .= "$discrep_line\n";
            }

    }
    
    close($input);  



    ## Identify reads with best alignments ##

    foreach my $read_name (keys %read_alignments)
    {
	    $stats{'num_reads_aligned'}++;
    
	    my @read_hsps = split(/\n/, $read_alignments{$read_name});
	    @read_hsps = sort by_align_score @read_hsps;
	    my $num_read_hsps = @read_hsps;
    
	    ## Get the best alignment ##
	    my $best_read_alignment = $read_hsps[0];
	    
	    ## Check for reads with competing alignments ##
	    if($num_read_hsps > 1)
	    {
		    my $second_read_alignment = $read_hsps[1];		
		    (my $temp, my $best_score) = split(/\s+/, $best_read_alignment);
		    ($temp, my $second_score) = split(/\s+/, $second_read_alignment);
		    if($best_score == $second_score)
		    {
			    ## Read has competing alignment ##
			    $stats{'reads_with_competing_aligns'}++;
			    $best_read_alignment = "";
		    }
	    }
    
	    ## Proceed if read has single best alignment ##
	    
	    if($best_read_alignment)
	    {
		    $stats{'reads_with_best_align'}++;
		    $best_alignments{$read_name} = $best_read_alignment;
                    if($alignment_discreps{$best_read_alignment})
                    {
                        $best_alignments{$read_name} .= "\n" . $alignment_discreps{$best_read_alignment};
                    }

	    }
    }




    %read_alignments = ();
    %best_align_scores = ();
    
    return(%best_alignments);
    
}


################################################################################

=head2	screen_for_snps - screens a read's best scored alignment for SNPs

=cut
################################################################################

sub screen_for_snps
{
    my $alignment_with_discreps = shift(@_);
    my $quality_db = shift(@_);

    my $detected_snps = "";    

    # Alignment variables #
    my $read_name = my $read_start = my $read_stop = my $ref_name = my $ref_start = my $ref_stop = "";
    my $align_score = my $align_strand = "";
    my $current_alignment = "";


    (my $alignment) = split(/\n/, $alignment_with_discreps);
#    (my $temp, $align_score, $ref_name, $ref_start, $ref_stop, $align_strand, $read_name, $read_start, $read_stop) = split(/\s+/, $alignment);

    ## Parse out alignment info ##

    my @lineContents = split(/\s+/, $alignment);
    my $numContents = @lineContents;
    $align_score = $lineContents[1];
    $read_name = $lineContents[5];

    if($lineContents[6] > $lineContents[7])
    {
            $read_start = $lineContents[7];
            $read_stop = $lineContents[6];			
    }
    else
    {
            $read_start = $lineContents[6];
            $read_stop = $lineContents[7];			
    }

    if($lineContents[9] eq 'C')
    {
            $align_strand = '-';
            
            $ref_name = $lineContents[10];
            $ref_start = $lineContents[13];
            $ref_stop = $lineContents[12];
    }
    else
    {
            $align_strand = '+';
            $ref_name = $lineContents[9];
            $ref_start = $lineContents[10];
            $ref_stop = $lineContents[11];
    }


#    my @lineContents = split(/\t/, $alignment);

    my @discreps = split(/\n/, $alignment_with_discreps);
    my $num_discreps = @discreps;
    
    for(my $discrepCounter = 1; $discrepCounter < $num_discreps; $discrepCounter++)
    {
        (my $temp, my $discrep_type, my $discrep_read_position, my $discrep_genotype, my $discrep_ref_position, my $discrep_context) = split(/\s+/, $discreps[$discrepCounter]);        

        if($discrep_type eq 'S')
        {
            (my $variant_allele, my $qual_score) = split(/[\(\)]/, $discrep_genotype);
            
            my $reference_allele = "N";
            ## If ref dir was given, get the actual reference allele ##
            
            if($reference_db)
            {
                $reference_allele = uc($reference_db->seq($ref_name, $discrep_ref_position => $discrep_ref_position));
            }
            
            if($variant_allele =~ /[ACGT]/)
            {
                $detected_snps .= "$ref_name\t$discrep_ref_position\t$reference_allele\t$variant_allele\t$read_name\t$discrep_read_position\t$align_strand\t$qual_score\t$discrep_context\n";
                $stats{'num_snps'}++;
            }
        }
    }

    return($detected_snps);
}




################################################################################

=head2	screen_for_indels - screens a read's best scored alignment for indels

=cut
################################################################################

sub screen_for_indels
{
    my $alignment = shift(@_);
#    my $quality_db = shift(@_);

    my $detected_indels = "";

    my @lineContents = split(/\t/, $alignment);


    return($detected_indels);
}




################################################################################

=head2	get_readcounts - gets read counts for a set of variants

=cut
################################################################################

sub get_readcounts
{
    (my $variants_file, my $alignments_file, my $output_file) = @_;

    my $verbose, my $quality_file;
    my $output_dir = "./";
    my $sample = "sample";

    $quality_file = $VarScan::quality_file;
    $default_qual_score = $VarScan::default_qual_score;
    $min_qual_score = $VarScan::min_qual_score;
    $min_align_score = $VarScan::min_align_score;
    $min_identity = $VarScan::min_identity;

    ## Get quality file if set up to do so ##
    
    if($quality_file && !(-e $quality_file))
    {
	print "\n***ERROR*** Quality file not found!\n\n";
    }

    ## Build Bio::DB::Fasta Objects ##
    
    my $fasta_db = my $quality_db;

    if($quality_file && -e $quality_file)
    {
	print "Building index for quality scores...\n";

	## Build the quality index ##
    
	mkdir("$output_dir/quality_dir") if(!(-d "$output_dir/quality_dir"));    
    
	## Convert the file ##
    
	my $quality_outfile = $output_dir . "/quality_dir/$sample.qual.fa";
	if(!(-e $quality_outfile))
	{
	    VarScan::SequenceFile::quality_to_fasta($quality_file, $quality_outfile);
	}
	
    #    if(!(-e "$output_dir/quality_dir/directory.index"))
	
	## Build the Bio::DB::Fasta index ##
	
	$quality_db = VarScan::SequenceFile::build_bio_db_fasta("$output_dir/quality_dir/");
    }


    my %stats = ();

    my %variant_regions = my %snps_by_position = my %indels_by_position = my %coverage_by_position = ();
    my $variant_type = "";
    
    ## Parse the combined variants ##
    
    my $input = new FileHandle ($variants_file);
    my $lineCounter = 0;
    while (<$input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;
	    
	    my @lineContents = split(/\t/, $line);

	    if($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name")
	    {
		
	    }
	    elsif($lineContents[6] && $lineContents[3] && ($lineContents[3] eq "INSERTION" || $lineContents[3] eq "DELETION")) ## INDELS ##
	    {
		    (my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $indel_context, my $num_reads, my $num_strands, my $avg_conf_score, my $avg_qual) = split(/\t/, $line);

		    my $indel_key = $chrom . ":" . $chr_start;
		    $indels_by_position{$indel_key} .= "$line\n";

		    ## Build a key using chromosome and first 2 bases of position for storing this indel ##
		    for(my $key_num = substr($chr_start, 0, 2); $key_num <= substr($chr_stop, 0, 2); $key_num++)
		    {
        		    my $chrom_key = $chrom . ":" . $key_num;
			    $variant_regions{$chrom_key}++;
		    }

		    $stats{$indel_type}++;
	    }
	    elsif($lineContents[0] && $lineContents[1] && $lineContents[2] && $lineContents[3]) ## SNPS ##
	    {
		    (my $chrom, my $position, my $allele1, my $allele2, my $num_reads, my $avg_qual, my $num_strands) = split(/\t/, $line);
    
		    my $position_key = $chrom . ":" . $position; #substr($position, 0, 2);
		    $snps_by_position{$position_key} += $num_reads;
		    
		    ## Build a key using chromosome and first 2 bases of position for storing this SNP ##
		    my $chrom_key = $chrom . ":" . substr($position, 0, 2);
		    $variant_regions{$chrom_key}++;			
		    
		    $stats{'SNP'}++;
	    }

    }
    
    close($input);
    
    print $stats{'SNP'} . " SNPs loaded\n" if($stats{'SNP'});
    print $stats{'INSERTION'} . " insertions loaded\n" if($stats{'INSERTION'});
    print $stats{'DELETION'} . " deletions loaded\n" if($stats{'DELETION'});
    
    my %num_ref_reads = my %sum_ref_qual = my %ref_strands_seen = my %num_var_reads = my %num_N_reads = ();    

    ## Parse out the best alignments ##
    print "Loading Read Alignments...\n";
    my %read_alignments = get_best_alignments($alignments_file);


    print "Screening Read Alignments...\n";

    ## Go through best read alignments ##
    foreach my $read_name (keys %read_alignments)
    {
	my $alignment_with_discreps = $read_alignments{$read_name};

        # Alignment variables #
        my $read_name = my $read_start = my $read_stop = my $ref_name = my $ref_start = my $ref_stop = "";
        my $align_score = my $align_strand = "";
        my $current_alignment = "";
    
        (my $alignment) = split(/\n/, $alignment_with_discreps);
    #    (my $temp, $align_score, $ref_name, $ref_start, $ref_stop, $align_strand, $read_name, $read_start, $read_stop) = split(/\s+/, $alignment);
    
        ## Parse out alignment info ##
    
        my @lineContents = split(/\s+/, $alignment);
        my $numContents = @lineContents;
        $align_score = $lineContents[1];
        $read_name = $lineContents[5];
    
        if($lineContents[6] > $lineContents[7])
        {
                $read_start = $lineContents[7];
                $read_stop = $lineContents[6];			
        }
        else
        {
                $read_start = $lineContents[6];
                $read_stop = $lineContents[7];			
        }
    
        if($lineContents[9] eq 'C')
        {
                $align_strand = '-';
                
                $ref_name = $lineContents[10];
                $ref_start = $lineContents[13];
                $ref_stop = $lineContents[12];
        }
        else
        {
                $align_strand = '+';
                $ref_name = $lineContents[9];
                $ref_start = $lineContents[10];
                $ref_stop = $lineContents[11];
        }    
  
        ## Get positions of any discrepancies ##

        my %discrep_positions = ();
        my @discreps = split(/\n/, $alignment_with_discreps);
        my $num_discreps = @discreps;
        
        for(my $discrepCounter = 1; $discrepCounter < $num_discreps; $discrepCounter++)
        {
            (my $temp, my $discrep_type, my $discrep_read_position, my $discrep_genotype, my $discrep_ref_position, my $discrep_context) = split(/\s+/, $discreps[$discrepCounter]);        
    
            if($discrep_type eq 'S')    ## SNPS 
            {
                (my $variant_allele, my $qual_score) = split(/[\(\)]/, $discrep_genotype);
                
                if($variant_allele eq 'N')
                {
                    $discrep_positions{$discrep_ref_position} = 'N';
                }
                elsif($variant_allele =~ /[ACGT]/)
                {
                    $discrep_positions{$discrep_ref_position} = $variant_allele;
                }
            }
            elsif($discrep_type =~ /[DI]/)  ## INDELS ##
            {
                $discrep_positions{$discrep_ref_position} = $discrep_type;
            }
        }

	
        ## Go through each base of the alignment, counting reads when necessary ##
        my $offset = 0;
        for(my $position = $ref_start; $position <= $ref_stop; $position++)
        {
            $offset++;
            my $position_key = $ref_name . ":" . $position;
            if($snps_by_position{$position_key})
            {
                ## Get quality score ##
                my $qual_score = $default_qual_score;                

                ## Infer the read position ##
                my $read_position;
                if($align_strand eq "-")
                {
                    $read_position = $read_stop - $offset;
                }
                else
                {
                    $read_position = $read_start + $offset;
                }
                $read_position = 1 if($read_position < 1);


                if($discrep_positions{$position} && $discrep_positions{$position} eq "N")
                {
                    $num_N_reads{$position_key}++;
                }
                elsif($discrep_positions{$position})
                {
                    $num_var_reads{$position_key}++;
                }
                else
                {
                    if($quality_db)
                    {
                        $qual_score = get_base_quality($quality_db, $read_name, $read_position, $read_position);
                    }	                    
                    
                    $ref_strands_seen{$position_key} .= $align_strand if(!$ref_strands_seen{$position_key} ||(length($ref_strands_seen{$position_key}) < 2 && $align_strand ne $ref_strands_seen{$position_key}));
                    $num_ref_reads{$position_key}++;
                    $sum_ref_qual{$position_key} += $qual_score;                    
                }
            }
        }
        
    }

    print "Saving Read Counts...\n";

    ## Open output file if desired ##

    if($output_file)
    {
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

	if($stats{'INSERTION'} || $stats{'DELETION'})
	{
	    print OUTFILE "chrom\tchr_start\tchr_stop\tindel_type\tindel_size\t";
	    print OUTFILE "reads1\treads2\treadsN\tavg_qual1\tavg_qual2\tstrands1\tstrands2\tindel_conf_score\tindel_in_context\n";
	}
	else
	{
	    print OUTFILE "chrom\tposition\tallele1\tallele2\t";
	    print OUTFILE "read_coverage\treads1\treads2\treadsN\t";
	    print OUTFILE "avg_ref_qual\tavg_var_qual\tnum_ref_strands\tnum_var_strands\tsnp_in_context\n";
	}
    }
    
    $input = new FileHandle ($variants_file);
    $lineCounter = 0;
    while (<$input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;
	    
	    my @lineContents = split(/\t/, $line);

	    if($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name")
	    {
		
	    }
	    elsif($lineContents[6] && $lineContents[3] && ($lineContents[3] eq "INSERTION" || $lineContents[3] eq "DELETION")) ## INDELS ##
	    {
		    (my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, , my $indel_context, my $reads2, my $num_var_strands, my $avg_conf_score, my $avg_var_qual) = split(/\t/, $line);
		    my $indel_key = "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size";
		    my $reads1 = $num_ref_reads{$indel_key};
		    my $readsN = $num_N_reads{$indel_key};
		    $reads1 = 0 if(!$reads1);
		    $readsN = 0 if(!$readsN);
		    my $read_coverage = $reads1 + $reads2 + $readsN;

		    my $num_ref_strands = 0;
		    $num_ref_strands = length($ref_strands_seen{$indel_key}) if($ref_strands_seen{$indel_key});

		    my $avg_ref_qual = $sum_ref_qual{$indel_key} / $reads1 if($reads1);
		    $avg_ref_qual = sprintf("%d", $avg_ref_qual) if($avg_ref_qual);
		    $avg_ref_qual = 0 if(!$avg_ref_qual);

		    print OUTFILE "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size\t";
		    print OUTFILE "$read_coverage\t$reads1\t$reads2\t$readsN\t";
		    print OUTFILE "$avg_ref_qual\t$avg_var_qual\t$num_ref_strands\t$num_var_strands\t$avg_conf_score\t$indel_context\n";

	    }
	    elsif($lineContents[0] && $lineContents[1] && $lineContents[2] && $lineContents[3]) ## SNPS ##
	    {
		    (my $chrom, my $position, my $allele1, my $allele2, my $reads2, my $avg_var_qual, my $num_var_strands, my $snp_in_context) = split(/\t/, $line);
		    my $position_key = $chrom . ":" . $position; #substr($position, 0, 2);
		    $num_ref_reads{$position_key} = 0 if(!$num_ref_reads{$position_key});
		    $num_var_reads{$position_key} = 0 if(!$num_var_reads{$position_key});
		    $num_N_reads{$position_key} = 0 if(!$num_N_reads{$position_key});
		    $snp_in_context = "" if(!$snp_in_context);
    
		    # Calculate read numbers ##
		    my $reads1 = $num_ref_reads{$position_key};
		    my $readsN = $num_N_reads{$position_key};				
		    my $read_coverage = $reads1 + $reads2 + $readsN;
    
		    my $avg_ref_qual = $sum_ref_qual{$position_key} / $reads1 if($reads1);
		    $avg_ref_qual = sprintf("%d", $avg_ref_qual) if($avg_ref_qual);
		    $avg_ref_qual = 0 if(!$avg_ref_qual);
    
		    my $num_ref_strands = 0;
		    $num_ref_strands = length($ref_strands_seen{$position_key}) if($ref_strands_seen{$position_key});
    
		    print OUTFILE "$chrom\t$position\t$allele1\t$allele2\t";
		    print OUTFILE "$read_coverage\t$reads1\t$reads2\t$readsN\t";
		    print OUTFILE "$avg_ref_qual\t$avg_var_qual\t$num_ref_strands\t$num_var_strands\t$snp_in_context\n";
		    
	    }
    }
    
    close($input);
    
    if($output_file)
    {
	close(OUTFILE);	
    }
}



################################################################################

=head2	get_base_quality - Retrieves quality score from Bio::DB::fasta index

=cut
################################################################################

sub get_base_quality
{
    (my $qual_db, my $read_name, my $read_start, my $read_stop) = @_;

    ## Get the quality ##
    my $qpiece = $qual_db->seq($read_name, $read_start => $read_stop);			
#    $qpiece = $qual_db->seq($read_name, ($read_position - 1) => ($read_position - 1)) if(!$qpiece);	
    ## Convert to numeric ##
    
    my $qpiece_num = "";
    
    for(my $baseCounter = 0; $baseCounter < length($qpiece); $baseCounter++)
    {
	    my $Qchar = substr($qpiece, $baseCounter, 1);
	    my $Qnum = ord($Qchar) - 33;
	    $qpiece_num .= " " if($qpiece_num);
	    $qpiece_num .= $Qnum;
    }			
    
    ## If multiple values were returned, take their average ##
    
    my @qual_scores = split(/\s+/, $qpiece_num);
    my $num_scores = @qual_scores;
    
    my $qual_score;

   ## Determine qual score ##


    if($num_scores > 1)
    {
	my $score_sum = 0;
	for(my $scoreCounter = 0; $scoreCounter < $num_scores; $scoreCounter++)
	{
	    $score_sum += $qual_scores[$scoreCounter];
	}
	
	$qual_score = sprintf("%d", $score_sum / $num_scores);
    }
    else
    {
	$qual_score = $qpiece_num;
    }
    
     
#    my $qual_score = $qpiece_num;
    
    
    return($qual_score);
}




=head1 AUTHOR

    Daniel C. Koboldt, << <dkoboldt at genome.wustl.edu> >>
    The Genome Center at Washington University School of Medicine
    St. Louis, Missouri, USA

=head1 COPYRIGHT

    Copyright 2009 Daniel C. Koboldt and Washington University
    All rights reserved.

=head1 LICENSE

    This program is free for non-commercial use.

=cut

1; # End of VarScan::ParseCrossmatch
