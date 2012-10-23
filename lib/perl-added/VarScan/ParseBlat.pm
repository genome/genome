package VarScan::ParseBlat;

use warnings;
use strict;
use Getopt::Long;

=head1 NAME

VarScan::ParseBlat - Parses output files from BLAT in PSLX format

=head1 VERSION

    Version 1.03

=head1 SYNOPSIS

    This module parses BLAT alignments in PSLX format

=cut

our $VERSION = '1.03';

## Define semi-global variables ##
my %stats = ();
my $verbose;
my $no_hp_filter;

## Define alignment scoring parameters ##
my $scoreM =                   	1;	# Match points for scoring alignments
my $scoreN =                   	2;      # Mismatch penalty for scoring alignments
my $scoreQ =                   	3;      # Gap penalty for scoring alignments	

## Define shared variables ##

my $reference_db = my $fasta_db = my $quality_db;
my $ref_dir, my $fasta_file, my $quality_file, my $min_identity, my $min_align_score;
my $primer_trim, my $default_qual_score, my $min_qual_score;
my $output_dir, my $sample, my $output_alignments, my $output_snps, my $output_indels;
my $min_indel_size, my $max_indel_size;
    
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
    $no_hp_filter = $VarScan::no_hp_filter;

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
	    (my $align_score) = split(/\t/, $best_read_alignment);
	    my $unscored_alignment = $best_read_alignment;
	    $unscored_alignment =~ s/$align_score\t//;
	    print ALIGNMENTS "$unscored_alignment\n" 
	}

	if($output_snps)
	{
	    my $read_snps = screen_for_snps($best_read_alignment, $quality_db);
	    print OUTFILE "$read_snps" if($output_snps && $read_snps);
	}

	if($output_indels)
	{
	    my $read_indels = screen_for_indels($best_read_alignment);
	    print INDELS "$read_indels" if($output_indels && $read_indels);
	}
    }


    ## Close the output file ##

    close(OUTFILE) if($output_snps);
    close(ALIGNMENTS) if($output_alignments);

    ## Print summary statistics ##

    $stats{'INSERTION'} = 0 if(!$stats{'INSERTION'});
    $stats{'DELETION'} = 0 if(!$stats{'DELETION'});
    $stats{'num_snps'} = 0 if(!$stats{'num_snps'});

    print "$stats{'num_alignments'} alignments were parsed\n";
    print "$stats{'num_reads_aligned'} reads had at least one alignment\n";
    print "$stats{'reads_with_best_align'} reads had a single best alignment\n";
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
	my @temp = split(/\t/, $a);
	my $score_a = $temp[0];
	@temp = split(/\t/, $b);
	my $score_b = $temp[0];		
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

    my $input = new FileHandle ($alignments_file);
    my $lineCounter = 0;

    while (<$input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;
    
	    ## Parse out only lines matching a BLAT-like result pattern ##

	    my @lineContents = split(/\t/, $line);		
	    my $numContents = @lineContents;
	    
#	    if($line && $line=~/\d+\t\d+\d+\t\d+\t\d+\d+\t/)	
	    if($numContents == 23)
	    {

		    my $numContents = @lineContents;
		    
		    my $match_bases = $lineContents[0];
		    my $mismatch_bases = $lineContents[1];
		    my $rep_match_bases = $lineContents[2];
		    my $query_gaps = $lineContents[4];
		    my $subject_gaps = $lineContents[6];
		    
		    my $read_name = $lineContents[9];
		    my $read_length = $lineContents[10];
		    my $read_start = $lineContents[11];
		    my $read_stop = $lineContents[12];				
		    
		    my $ref_name = $lineContents[13];
		    my $num_blocks = $lineContents[17];
		    my $alignment_gaps = $num_blocks - 1;

		    ## Calculate read-proportional length of the alignment ##
		    
		    my $pct_read_aligned = ($read_stop - $read_start) / $read_length * 100;
		    
		    ## Calculate the sequence identity of the alignment ##
		    my $alignment_identity = 0;
		    $alignment_identity = ($match_bases + $rep_match_bases) / ($match_bases + $rep_match_bases + $mismatch_bases) * 100 if($match_bases);			
	    
		    ## Calculate a BLAST-like alignment score ##
		    
		    my $blast_score = (($match_bases + $rep_match_bases) * $scoreM) - ($mismatch_bases * $scoreN) - ($query_gaps * $scoreQ) - ($subject_gaps * $scoreQ);

		    if($blast_score && $blast_score >= $min_align_score)
		    {
			    if($alignment_identity && $alignment_identity >= $min_identity)
			    {
				    ## Save the score along with the read alignment ##
				    $stats{'num_alignments'}++;	

				    if($verbose)
				    {
					    print $stats{'num_alignments'} . " alignments parsed...\n" if(!($stats{'num_alignments'} % 10000))
				    }

				    if(!$read_alignments{$read_name} || $blast_score >= $best_align_scores{$read_name})
				    {
					$read_alignments{$read_name} .= "\n" if($read_alignments{$read_name});
					$read_alignments{$read_name} .= $blast_score . "\t" . $line;
					$best_align_scores{$read_name} = $blast_score;
				    }
			    }
		    }
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
		    (my $best_score) = split(/\t/, $best_read_alignment);
		    (my $second_score) = split(/\t/, $second_read_alignment);
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
    my $alignment = shift(@_);
    my $quality_db = shift(@_);

    my $detected_snps = "";    
    my @lineContents = split(/\t/, $alignment);
    ## Parse various values from the HSP line ##

    my $align_score = $lineContents[0];
    my $mismatch_bases = $lineContents[2];
    my $query_gap_count = $lineContents[5];
    my $query_gap_bases = $lineContents[6];
    my $chrom_gap_count = $lineContents[7];
    my $chrom_gap_bases = $lineContents[8];
    my $align_strand = $lineContents[9];
    my $read_name = $lineContents[10];
    my $read_length = $lineContents[11];
    my $read_start = $lineContents[12];
    my $read_stop = $lineContents[13];
    my $ref_name = $lineContents[14];
    my $ref_start = $lineContents[16];
    my $ref_stop = $lineContents[17];
    my $num_blocks = $lineContents[18];
    my $block_sizes = $lineContents[19];
    my $read_block_starts = $lineContents[20];
    my $ref_block_starts = $lineContents[21];
    my $read_block_seqs = $lineContents[22];
    my $ref_block_seqs = $lineContents[23];

    ## Put block info into arrays ##
		    
    my @blockSizes = split(/,/, $block_sizes);
    my @readBlockStarts = split(/,/, $read_block_starts);
    my @refBlockStarts = split(/,/, $ref_block_starts);					
    my @readBlockSeqs = split(/,/, $read_block_seqs);
    my @refBlockSeqs = split(/,/, $ref_block_seqs);						
    
    
    ## Reset some variables ##
    
    my $read_block_start = my $ref_block_start = my $read_block_stop = my $ref_block_stop = 0;
    my $read_block_seq = my $ref_block_seq;
    my %read_snps = ();

    ## Iterate through each alignment block, looking at the gaps between them for discrepancies ##
    
    for(my $bCounter = 0; $bCounter < $num_blocks; $bCounter++)
    {
	    ## Parse out this block's info ##
	    
	    my $block_size = $blockSizes[$bCounter];
	    
	    $read_block_start = $readBlockStarts[$bCounter];
	    $ref_block_start = $refBlockStarts[$bCounter];
	    $read_block_stop = $read_block_start + $block_size - 1;
	    $ref_block_stop = $ref_block_start + $block_size - 1;		
	    $read_block_seq = $readBlockSeqs[$bCounter];
	    $ref_block_seq = $refBlockSeqs[$bCounter];							

	    ## Look for SNPs ##

	    if($read_block_seqs && $ref_block_seqs && $mismatch_bases > 0)
	    {
		    my $prev_read_base = my $prev_ref_base = "";
				    
		    for(my $baseCounter = 0; $baseCounter < $block_size; $baseCounter++)
		    {
			    ## Get ref and read positions ##
			    my $ref_position = $ref_block_start + $baseCounter + 1;
			    my $read_position = $read_block_start + $baseCounter + 1;
			    $read_position = $read_length - $read_position if($align_strand eq '-');							
			    my $read_base = uc(substr($read_block_seq, $baseCounter, 1));
			    my $ref_base = uc(substr($ref_block_seq, $baseCounter, 1));
			    
			    $stats{'num_bases'}++;
							    
			    if($read_base ne 'N')
			    {
				    if($read_base ne $ref_base)
				    {
					    my $flank_5 = my $flank_3 = "";	
					    $flank_5 = uc(substr($ref_block_seq, $baseCounter - 5, 5));
					    $flank_3 = uc(substr($ref_block_seq, $baseCounter + 1, 5));
					    my $context = lc($flank_5) . $read_base . lc($flank_3);
					    
					    ## Calculate distance from end of read ##
					    my $dist_from_end = $read_length - $read_position;
					    
					    ## Call SNP if not within the primer-trim section ##
					    
					    if(!$primer_trim || ($read_position > $primer_trim && $dist_from_end > $primer_trim))
					    {
						    ## Get quality score ##
						    
						    my $qual_score = $default_qual_score;

						    if($quality_db)
						    {
							$qual_score = get_base_quality($quality_db, $read_name, $read_position, $read_position);
#							print "Qual score was $qual_score\n";
#							exit(0);
						    }						    

#						    print OUTFILE "$ref_name\t$ref_position\t$ref_base\t$read_base\t$read_name\t$read_position\t$align_strand\t$context\n";
						    $detected_snps .= "$ref_name\t$ref_position\t$ref_base\t$read_base\t$read_name\t$read_position\t$align_strand\t$qual_score\t$context\n";

						    $stats{'num_snps'}++;
						    
						    $read_snps{$read_position} = "$read_position\t$ref_base\t$read_base\t$context";
						    
						    ## Check for DNPs ##
						    
						    my $prev_read_position = $read_position - 1;
						    if($read_snps{$prev_read_position})
						    {
							(my $prev_position, my $prev_ref_base, my $prev_read_base, my $prev_context) = split(/\t/, $read_snps{$prev_read_position});
							my $dnp_ref_position = $ref_position - 1;
							my $dnp_ref_bases = $prev_ref_base . $ref_base;
							my $dnp_read_bases = $prev_read_base . $read_base;
							my $dnp_flank_5 = uc(substr($ref_block_seq, $baseCounter - 5 - 1, 5));                                                                            
							my $dnp_context = lc($dnp_flank_5) . $dnp_read_bases . lc($flank_3);
#							    print DNSUBSTITUTIONS "$ref_name\t$dnp_ref_position\t$dnp_ref_bases\t$dnp_read_bases\t$read_name\t$prev_position\t$align_strand\t$dnp_context\n";
							$stats{'num_dnps'}++;
						    }
					    }
				    }
			    }
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

    my $detected_snps = "";    
    my @lineContents = split(/\t/, $alignment);
    ## Parse various values from the HSP line ##

    my $align_score = $lineContents[0];
    my $mismatch_bases = $lineContents[2];
    my $query_gap_count = $lineContents[5];
    my $query_gap_bases = $lineContents[6];
    my $chrom_gap_count = $lineContents[7];
    my $chrom_gap_bases = $lineContents[8];
    my $align_strand = $lineContents[9];
    my $read_name = $lineContents[10];
    my $read_length = $lineContents[11];
    my $read_start = $lineContents[12];
    my $read_stop = $lineContents[13];
    my $ref_name = $lineContents[14];
    my $ref_start = $lineContents[16];
    my $ref_stop = $lineContents[17];
    my $num_blocks = $lineContents[18];
    my $block_sizes = $lineContents[19];
    my $read_block_starts = $lineContents[20];
    my $ref_block_starts = $lineContents[21];
    my $read_block_seqs = $lineContents[22];
    my $ref_block_seqs = $lineContents[23];

    ## Put block info into arrays ##
		    
    my @blockSizes = split(/,/, $block_sizes);
    my @readBlockStarts = split(/,/, $read_block_starts);
    my @refBlockStarts = split(/,/, $ref_block_starts);					
    my @readBlockSeqs = split(/,/, $read_block_seqs);
    my @refBlockSeqs = split(/,/, $ref_block_seqs);						
    
    
    ## Reset some variables ##
    
    my $read_block_start = my $ref_block_start = my $read_block_stop = my $ref_block_stop = 0;
    my $prev_read_block_stop = my $prev_ref_block_stop = 0;
    my $prev_ref_block_seq = "";	
    my $read_block_seq = my $ref_block_seq = "";				
    my $read_indels = my $read_snps = 0;		

    ## Iterate through each alignment block, looking at the gaps between them for discrepancies ##
    
    for(my $bCounter = 0; $bCounter < $num_blocks; $bCounter++)
    {
	## Parse out this block's info ##
	
	my $block_size = $blockSizes[$bCounter];
	
	$read_block_start = $readBlockStarts[$bCounter];
	$ref_block_start = $refBlockStarts[$bCounter];
	$read_block_stop = $read_block_start + $block_size - 1;
	$ref_block_stop = $ref_block_start + $block_size - 1;		
	$read_block_seq = $readBlockSeqs[$bCounter];
	$ref_block_seq = $refBlockSeqs[$bCounter];							
	
	## Look for indels ##

	if($read_block_seqs && $ref_block_seqs && $bCounter > 0)
	{
		## Get prev block size ##
		my $prev_block_size = $blockSizes[$bCounter - 1];
		
		## Look for indels ##							
            print "GIving it  $prev_ref_block_stop\n";
		my @indel = findGapIndels($read_block_start, $ref_block_start, $read_block_stop, $ref_block_stop, $prev_read_block_stop, $prev_ref_block_stop, $align_strand, $read_length, $read_name);
	
		if(@indel && $indel[0] && $indel[1])
		{
			(my $indel_type, my $indel_size, my $ref_indel_start, my $ref_indel_stop, my $read_indel_start, my $read_indel_stop, my $discrep_penalty) = @indel;							
			
			if($align_strand eq '-')
			{
				$read_indel_start = $read_length - $read_indel_start - 1;
				$read_indel_stop = $read_length - $read_indel_stop - 1;
			}
			
			## Adjust zero-based positions ##
			$ref_indel_start++;
			$ref_indel_stop++;
			$read_indel_start++;
			$read_indel_stop++;
			
			## Get the indel sequence ##
			
			my $indel_context = my $indel_seq = my $indel_qual = "";			
			my $reference_allele = my $variant_allele = "";
			#my $filter_status = "";		
			my $indel_in_context = my $indel_flank_5 = my $indel_flank_3 = "";						
			## Get indel-in-context and determine alleles ##
		
			if($indel_type eq "INSERTION")
			{
			    if($reference_db && $fasta_db)
			    {
				my $ref_context;

				if($reference_db)
				{
				    $ref_context = uc($reference_db->seq($ref_name, ($ref_indel_start - 10) => ($ref_indel_stop + 10)));
				    ## Separate out flank 5 and flank 3 ##
				    $indel_flank_5 = substr($ref_context, 0, 10);
				    $indel_flank_3 = substr($ref_context, (length($ref_context) - 10), 10);							
				}

				## Get reverse complement if read aligned to - strand ##
				
				if($align_strand eq '-' && $read_indel_start >= $read_indel_stop)
				{
					$indel_context = uc($fasta_db->seq($read_name, ($read_indel_start + 10) => ($read_indel_stop - 10)));
					$indel_qual = get_base_quality($quality_db, $read_name, $read_indel_stop, $read_indel_start) if($quality_db);
				}
				else
				{
					$indel_context = uc($fasta_db->seq($read_name, ($read_indel_start - 10) => ($read_indel_stop + 10)));					
					$indel_qual = get_base_quality($quality_db, $read_name, $read_indel_start, $read_indel_stop) if($quality_db);
				}
				
				## If indel is at start of read, adjust the parsing of the insertion sequence ##
				my $context_offset = 10;
				$context_offset = $read_indel_start if($read_indel_start <= 10);
			    

				$variant_allele = substr($indel_context, $context_offset, (abs($read_indel_stop - $read_indel_start) + 1)) if($indel_context);
				
				## Special treatment for insertions with ambiguous coordinates ##
				
				if($ref_indel_start ne $ref_indel_stop && $ref_context && $indel_context)
				{
					## Get the ambiguous ref bases ##
					$reference_allele = substr($ref_context, 10, (length($ref_context) - 20));
					while(length($reference_allele) < length($variant_allele))
					{
						$reference_allele .= "-";
					}
				}
				else
				{
					$reference_allele = '-' x length($variant_allele);
				}
				
				$indel_seq = $variant_allele;
			    }
			    else
			    {
				    ## A quick hack to use the block sequences ##
				    $indel_flank_5 = uc(substr($readBlockSeqs[$bCounter - 1], length($readBlockSeqs[$bCounter - 1]) - 10, 10));
				    $indel_flank_3 = uc(substr($readBlockSeqs[$bCounter], 0, 10));
				
					$variant_allele = $indel_size . " bp";
					$reference_allele = "-";
					$indel_seq = $indel_size . " bp";                                                                
			    }
			}
			elsif($indel_type eq "DELETION")
			{
			    if($reference_db)
			    {
				$indel_context = uc($reference_db->seq($ref_name, ($ref_indel_start - 10) => ($ref_indel_stop + 10)));

				## Separate out flank 5 and flank 3 ##
				$indel_flank_5 = substr($indel_context, 0, 10);
				$indel_flank_3 = substr($indel_context, (length($indel_context) - 10), 10);
			    }
			    else
			    {

				## A quick hack to use the block sequences ##
				$indel_flank_5 = uc(substr($readBlockSeqs[$bCounter - 1], length($readBlockSeqs[$bCounter - 1]) - 10, 10));
				$indel_flank_3 = uc(substr($readBlockSeqs[$bCounter], 0, 10));
			    }
			    
			    if($quality_db)
			    {
				if($align_strand eq '-' && $read_indel_start < $read_indel_stop)
				{
					$indel_qual = get_base_quality($quality_db, $read_name, $read_indel_stop, $read_indel_start);
				}
				else
				{
					$indel_qual = get_base_quality($quality_db, $read_name, $read_indel_start, $read_indel_stop);
				}					
			    }	
				## Set a max indel size to obtain sequence for ##
				
				if($reference_db && $indel_size < 150)
				{
					## Special treatment for deletions with ambiguous coordinates ##
					
					if($read_indel_start ne $read_indel_stop)									
					{
						$reference_allele = uc(substr($indel_context, 10, ($ref_indel_stop - $ref_indel_start + 1)));
						$variant_allele = uc($fasta_db->seq($read_name, $read_indel_start => $read_indel_stop));
						## Flip this sequence if align strand was - ##
						
						$variant_allele = reverse_complement($variant_allele) if($align_strand eq '-');
						
						## Add deletion dashes until length is right ##
						while(length($variant_allele) < length($reference_allele))
						{
							$variant_allele .= "-";
						}
					}
					elsif($indel_context)
					{
						$reference_allele = uc(substr($indel_context, 10, ($ref_indel_stop - $ref_indel_start + 1)));
						$variant_allele = '-' x length($reference_allele);
						$indel_seq = $reference_allele;
					}
				}
				else
				{
					$reference_allele = $indel_size . " bp";
					$variant_allele = "-";
					$indel_seq = $indel_size . " bp";
				}
			}

			$indel_in_context =  $indel_flank_5 . '[' . $reference_allele . '/' . $variant_allele . ']' . "$indel_flank_3";
			
			## Determine discrepancy penalty ##
			
			my $discrep_pct;
			
			if($indel_size < 10 && $discrep_penalty < 10)
			{
				$discrep_pct = $discrep_penalty * 10;
			}
			else
			{
				$discrep_pct = sprintf("%d", ($discrep_penalty / ($discrep_penalty + $indel_size) * 100));
			}
			my $conf_score = 100 - $discrep_pct;
			
			$read_indels++;
			
			my $indel_event = "$ref_name\t$ref_indel_start\t$ref_indel_stop\t$indel_type\t$indel_size\t$read_name\t$read_indel_start\t$read_indel_stop\t$align_strand\t$conf_score\t$indel_in_context\t$indel_qual";
			
			my $filter_status = VarScan::VariantCalling::heuristic_indel_filter($indel_event);
			
			if($filter_status && !$no_hp_filter)
			{
				$stats{'indels_filtered'}++;
			#	$FilterCounts{$filter_status}++;
			}
			else
			{
				$stats{$indel_type}++;
				
				$detected_indels .= "$indel_event\n";
							
			}
			
		} # else !@indel

								
	} # else (!$read_block_seqs || !$ref_block_seqs)

            $prev_ref_block_seq = $ref_block_seq;			
            $prev_read_block_stop = $read_block_stop;
            $prev_ref_block_stop = $ref_block_stop;												

    }    

    return($detected_indels);
}




################################################################################################
# findGapIndels - call indels between PSL alignment blocks
#
################################################################################################

sub findGapIndels
{
	(my $read_block_start, my $ref_block_start, my $read_block_stop, my $ref_block_stop, my $prev_read_block_stop, my $prev_ref_block_stop, my $align_strand, my $read_length, my $read_name) = @_;
	#$read_block_start,        $ref_block_start,    $read_block_stop,    $ref_block_stop,    $prev_read_block_stop,    $prev_ref_block_stop,    $align_strand,    $read_length,    $read_name
	my $read_gap = $read_block_start - $prev_read_block_stop - 1;
	my $ref_gap = $ref_block_start - $prev_ref_block_stop - 1;

	my $indel_type = my $indel_size = my $ref_indel_start = my $ref_indel_stop = 0;
	my $read_indel_start = my $read_indel_stop = 0;
	my $discrep_penalty = 0;
print "read block: $read_block_start-$read_block_stop ref block: $ref_block_start-$ref_block_stop read gap: $read_gap ref_gap: $ref_gap\nmy $ref_gap = $ref_block_start - $prev_ref_block_stop - 1;\n";
	## Use gap sizes to determine type of sequence variant ##
	## If read contains a deletion, read gap should be 0, ref gap should be deletion size ##
	if($read_gap  < $ref_gap)
	{
		$indel_type = 'DELETION';
		$indel_size = $ref_gap - $read_gap;
		$ref_indel_start = $prev_ref_block_stop + 1;
		$ref_indel_stop = $ref_block_start - 1;					
		$read_indel_start = $prev_read_block_stop + 1;
		$read_indel_stop = $read_block_start - 1;					
		$read_indel_stop = $read_indel_start if($read_indel_stop < $read_indel_start);
		
		if($align_strand eq "-")
		{
			$read_indel_start = $read_length - $read_indel_start;
			$read_indel_stop = $read_length - $read_indel_stop;
		}
	
		$discrep_penalty = $read_gap;	#5 * $gap_discrep;
	}
	## If read contains an insertion, read gap is indel size and ref gap should be 0 ##					
	elsif($read_gap > $ref_gap)
	{
		$indel_type = 'INSERTION';
		$indel_size = $read_gap - $ref_gap;	
		$ref_indel_start = $prev_ref_block_stop + 1;
		$ref_indel_stop = $ref_block_start - 1;
		$ref_indel_stop = $ref_indel_start if($ref_indel_stop < $ref_indel_start);
		$read_indel_start = $prev_read_block_stop + 1;
		$read_indel_stop = $read_block_start - 1;					
		
		$discrep_penalty = $ref_gap;	#5 * $gap_discrep

	}
	else
	{
		warn "Read gap and ref gap the same! $read_name read_gap=$read_gap ref_gap=$ref_gap start:$read_block_start,$ref_block_start stop:$read_block_stop,$ref_block_stop prev_stop:$prev_read_block_stop,$prev_ref_block_stop\n";
	}
	
	if($indel_type && $indel_size && $indel_size >= $min_indel_size && $indel_size <= $max_indel_size)
	{
		return($indel_type, $indel_size, $ref_indel_start, $ref_indel_stop, $read_indel_start, $read_indel_stop, $discrep_penalty);
	}
	else
	{
		return(0);
	}
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
		    $snps_by_position{$position_key} .= $allele2;
		    
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
	my $alignment = $read_alignments{$read_name};
	my @lineContents = split(/\t/, $alignment);
	## Parse various values from the HSP line ##
    
	my $align_score = $lineContents[0];
	my $mismatch_bases = $lineContents[2];
	my $query_gap_count = $lineContents[5];
	my $query_gap_bases = $lineContents[6];
	my $chrom_gap_count = $lineContents[7];
	my $chrom_gap_bases = $lineContents[8];
	my $align_strand = $lineContents[9];
	my $read_name = $lineContents[10];
	my $read_length = $lineContents[11];
	my $read_start = $lineContents[12];
	my $read_stop = $lineContents[13];
	my $ref_name = $lineContents[14];
	my $ref_start = $lineContents[16];
	my $ref_stop = $lineContents[17];
	my $num_blocks = $lineContents[18];
	my $block_sizes = $lineContents[19];
	my $read_block_starts = $lineContents[20];
	my $ref_block_starts = $lineContents[21];
	my $read_block_seqs = $lineContents[22];
	my $ref_block_seqs = $lineContents[23];
    
	## Put block info into arrays ##
			
	my @blockSizes = split(/,/, $block_sizes);
	my @readBlockStarts = split(/,/, $read_block_starts);
	my @refBlockStarts = split(/,/, $ref_block_starts);					
	my @readBlockSeqs = split(/,/, $read_block_seqs);
	my @refBlockSeqs = split(/,/, $ref_block_seqs);						
	
	
	## Reset some variables ##
	
	my $read_block_start = my $ref_block_start = my $read_block_stop = my $ref_block_stop = 0;
	my $read_block_seq = my $ref_block_seq;
	my %read_snps = ();
    
	## Iterate through each alignment block, looking at the gaps between them for discrepancies ##
	
	for(my $bCounter = 0; $bCounter < $num_blocks; $bCounter++)
	{
		## Parse out this block's info ##
		
		my $block_size = $blockSizes[$bCounter];
		
		$read_block_start = $readBlockStarts[$bCounter];
		$ref_block_start = $refBlockStarts[$bCounter];
		$read_block_stop = $read_block_start + $block_size - 1;
		$ref_block_stop = $ref_block_start + $block_size - 1;		
		$read_block_seq = $readBlockSeqs[$bCounter];
		$ref_block_seq = $refBlockSeqs[$bCounter];							
    
		## Go through bases and count reads for SNPs  ##
    
		if($read_block_seqs && $ref_block_seqs && $stats{'SNP'})
		{
			my $prev_read_base = my $prev_ref_base = "";
					
			for(my $baseCounter = 0; $baseCounter < $block_size; $baseCounter++)
			{
				## Get ref and read positions ##
				my $ref_position = $ref_block_start + $baseCounter + 1;
				my $read_position = $read_block_start + $baseCounter + 1;
				$read_position = $read_length - $read_position if($align_strand eq '-');							
				my $read_base = uc(substr($read_block_seq, $baseCounter, 1));
				my $ref_base = uc(substr($ref_block_seq, $baseCounter, 1));
				
				$stats{'num_bases'}++;
				
				my $position_key = $ref_name . ":" . $ref_position;
				if($snps_by_position{$position_key})
				{
				    ## Calculate distance from end of read ##
				    my $dist_from_end = $read_length - $read_position;
				    
				    ## Call base if not within the primer-trim section ##
				    
				    if(!$primer_trim || ($read_position > $primer_trim && $dist_from_end > $primer_trim))
				    {
					## Get quality score ##
					
					my $qual_score = $default_qual_score;

					if($quality_db)
					{
					    $qual_score = get_base_quality($quality_db, $read_name, $read_position, $read_position);
					}

					if($read_base eq 'N')
					{
					    $num_N_reads{$position_key}++;
					}
					elsif($read_base ne $ref_base)
					{
					    $num_var_reads{$position_key}++;
					}
					else
					{
					    $ref_strands_seen{$position_key} .= $align_strand if(!$ref_strands_seen{$position_key} ||(length($ref_strands_seen{$position_key}) < 2 && $align_strand ne $ref_strands_seen{$position_key}));
					    $num_ref_reads{$position_key}++;
					    $sum_ref_qual{$position_key} += $qual_score;					    
					}
				    }
				}

			}						
		}
		
		## Go through blocks and see if we can validate indels ##
		if($stats{'INSERTION'} || $stats{'DELETION'})
		{
		    for(my $baseCounter = 0; $baseCounter < $block_size; $baseCounter++)
		    {
			    ## Get ref and read positions ##
			    my $ref_position = $ref_block_start + $baseCounter + 1;
			    my $read_position = $read_block_start + $baseCounter + 1;
			    $read_position = $read_length - $read_position if($align_strand eq '-');

#		    for(my $position = $ref_block_start; $position <= $ref_block_stop; $position++)
#		    {
			my $position_key = $ref_name. ":" . $ref_position;

			if($indels_by_position{$position_key})
			{
			    (my $indel_chrom, my $indel_start, my $indel_stop, my $indel_type, my $indel_size) = split(/\t/, $indels_by_position{$position_key});
			    my $indel_key = "$indel_chrom\t$indel_start\t$indel_stop\t$indel_type\t$indel_size";

			    ## If this block spans the indel with no gaps, call it a reference read ##
			    if($ref_block_start <= ($indel_start - 2) && $ref_block_stop >= ($indel_stop + 2))
			    {
				## Get base quality score ##
				
				my $qual_score = $default_qual_score;
    
				if($quality_db)
				{
				    $qual_score = get_base_quality($quality_db, $read_name, $read_position, $read_position);
				}
				
				$sum_ref_qual{$indel_key} += $qual_score;
				$num_ref_reads{$indel_key}++;
				$ref_strands_seen{$indel_key} .= $align_strand if(!$ref_strands_seen{$indel_key} ||(length($ref_strands_seen{$indel_key}) < 2 && $align_strand ne $ref_strands_seen{$indel_key}));
			    }
			    else
			    {
				## Cannot make a determination if this read block supports/refutes indel ##
				$num_N_reads{$indel_key}++;
			    }
			}
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

    return($default_qual_score) if(!$qual_db);

    ## Get the quality ##
    my $qpiece = $qual_db->seq($read_name, $read_start => $read_stop);			
#    $qpiece = $qual_db->seq($read_name, ($read_position - 1) => ($read_position - 1)) if(!$qpiece);	

    return($default_qual_score) if(!$qpiece);

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




##########################
# Reverse complement
##########################

sub reverse_complement
{
	my $seq = shift(@_);
	$seq = uc($seq);
	my $seq_len = length($seq);
		
	my $rcomp_seq = "";
	
	for(my $baseCounter = 0; $baseCounter < $seq_len; $baseCounter++)
	{
		## Complement the current base ##
		my $this_base = substr($seq, $baseCounter, 1);
		## Complement the current base ##	
		my $comp_base = 'N';		
		$comp_base = 'T' if($this_base eq 'A');
		$comp_base = 'G' if($this_base eq 'C');
		$comp_base = 'C' if($this_base eq 'G');			
		$comp_base = 'A' if($this_base eq 'T');	
		## Put the base at the beginning ##
		
		$rcomp_seq = $comp_base . $rcomp_seq;	
	}

	return($rcomp_seq);
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

1; # End of VarScan::ParseBlat
