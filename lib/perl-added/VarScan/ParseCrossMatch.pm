package VarScan::ParseCrossMatch;

use warnings;
use strict;

=head1 NAME

VarScan::ParseCrossMatch - modules for parsing cross_match output

=head1 VERSION

Version 0.01

=cut

our $VERSION = '1.01';

=head1 SYNOPSIS

This module contains subroutines for parsing cross_match output to obtain best alignments, detect variants, etc.

=head1 FUNCTIONS

################################################################################################
=head2 parse_alignments - a subroutine to parse read alignments in Newbler pairt format
################################################################################################
=cut

sub parse_alignments
{
	(my $FileName, my $output_dir, my $sample_name) = @_;

	my %Stats = ();
	my %ReadAlign = ();

	$Stats{'alignments'} = $Stats{'aligned_reads'} = 0;

	# Alignment variables #
	my $read_name = my $read_start = my $read_stop = my $ref_name = my $ref_start = my $ref_stop = "";
	my $align_score = my $align_strand = "";
        my $current_alignment = "";

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

        ## Pre-parse the file to get the best alignments for each read ##
        my %ReadAlignments = my %AlignmentDiscreps = ();
	$input = new FileHandle ($FileName);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter >= 1 && $line && substr($line, 0, 9) eq "ALIGNMENT")
		{
                        $Stats{'alignments'}++;
                
                        my @lineContents = split(/\s+/, $line);
                        my $numContents = @lineContents;
			$align_score = $lineContents[1];
			$read_name = $lineContents[5];
                
                        if(!$ReadAlign{$read_name})
                        {
                                $Stats{'aligned_reads'}++;
                                $ReadAlign{$read_name} = 1;
                        }		


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

                        $current_alignment = "$align_score\t$ref_name\t$ref_start\t$ref_stop\t$align_strand\t$read_name\t$read_start\t$read_stop";
                        $ReadAlignments{$read_name} .= "$current_alignment\n";

                }
		elsif($lineCounter >= 1 && $line && substr($line, 0, 11) eq "DISCREPANCY")
		{
#                    my @lineContents = split(/\s+/, $line);
                    my $discrep_line = $line;
                    $discrep_line =~ s/\s+/\t/g;
                    $AlignmentDiscreps{$current_alignment} .= "$discrep_line\n";
                }
        }
        
        close($input);


	## Open some output files ##
	
	open(ALIGNMENTS, ">$output_dir/$sample_name.best.aligns.tsv") or die "Can't open outfile: $!\n";
	print ALIGNMENTS "ref_name\tref_start\tref_stop\talign_strand\tread_name\tread_start\tread_stop\talign_score\tgapped_id\tungapped_id\n";
	
	open(SCOREDALIGNS, ">$output_dir/$sample_name.best.aligns.cm") or die "Can't open outfile: $!\n";        
#        print SCOREDALIGNS "AlignScore\tQueryAccno\tQueryStart\tQueryEnd\tQueryLength\tSubjAccno\tSubjStart\tSubjEnd\tSubjLength\tNumIdent\tAlignLength\tQueryAlign\tSubjAlign\n";
        
        open(BLOCKS, ">$output_dir/$sample_name.best.blocks.tsv") or die "Can't open outfile: $!\n";
        print BLOCKS "ref_name\tref_start\tref_stop\talign_strand\tread_name\tread_start\tread_stop\tblock_no\n";
        
        ## Iterate through alignments on per-read basis, selecting only the best for the "best" alignments ##

        foreach my $read_name (keys %ReadAlignments)
        {
            my @aligns = split(/\n/, $ReadAlignments{$read_name});
            @aligns = sort byAlignScoreDesc @aligns;
            my @best = split(/\t/, $aligns[0]);
            my @second = split(/\t/, $aligns[1]) if($aligns[1]);
            if(!$aligns[1] || $best[0] > $second[0])
            {
                $Stats{'best_alignments'}++;
                
                my $align_discreps = $AlignmentDiscreps{$aligns[0]};
                print SCOREDALIGNS "ALIGNMENT\t$aligns[0]\n";
                print SCOREDALIGNS $AlignmentDiscreps{$aligns[0]} if($AlignmentDiscreps{$aligns[0]});
                
                my $align_score = $best[0];
                my $ref_name = $best[1];
                my $ref_start = $best[2];
                my $ref_stop = $best[3];
                my $align_strand = $best[4];
                my $read_name = $best[5];
                my $read_start = $best[6];
                my $read_stop = $best[7];

              
                if($align_discreps)
                {
                    ## Retrieve any discrepancies within the alignment ##
                    
                    my @discreps = split(/\n/, $align_discreps);
                    my %indels_by_read_pos = my %indels_by_ref_pos = my %snps_by_read_pos = my %snps_by_ref_pos = ();
                    foreach my $discrep (@discreps)
                    {
                        (my $temp, my $discrep_type, my $discrep_read_position, my $discrep_genotype, my $discrep_ref_position, my $discrep_context) = split(/\t/, $discrep);
                        
                        my $indel_type = my $indel_size = "";
                        
                        if($discrep_type =~ 'I' || $discrep_type =~ 'D')
                        {
                            $indels_by_read_pos{$discrep_read_position} = $discrep_type;
                            $indels_by_ref_pos{$discrep_ref_position} = $discrep_type;
                            
                            if(length($discrep_type) > 1)
                            {
                                ($discrep_type, $indel_size) = split(/\-/, $discrep_type);
                                $indel_type = "INSERTION" if($discrep_type eq "I");
                                $indel_type = "DELETION" if($discrep_type eq "D");
                            }
                            else
                            {
                                $indel_type = "INSERTION" if($discrep_type eq "I");
                                $indel_type = "DELETION" if($discrep_type eq "D");
                                $indel_size = 1;
                            }

                        }
                        else
                        {
                            $snps_by_read_pos{$discrep_read_position} = 1;
                            $snps_by_ref_pos{$discrep_ref_position} = 1;
                        }
                    }
                    
                    ## Prepare to handle block information ##
                    
                    my $block_number = my $block_read_start = my $block_read_stop = my $block_ref_start = my $block_ref_stop = 0;
                    
                    my $match = my $mismatch = my $gaps = 0;

                    if($align_strand eq '+')
                    {
                        $block_read_start = $read_start;
                        $block_read_stop = $read_start;
                        
                        $block_ref_start = $ref_start;
                        $block_ref_stop = $ref_start;
                        
                        for(my $read_position = $read_start; $read_position <= $read_stop; $read_position++)
                        {
                            if($snps_by_read_pos{$read_position})
                            {
                                $mismatch++;
                            }
                            elsif($indels_by_read_pos{$read_position})
                            {
                                ## Block ends here ##
                                $block_number++;
                                $block_read_stop++;
                                $block_ref_stop++;
                                print BLOCKS "$ref_name\t$block_ref_start\t$block_ref_stop\t$align_strand\t$read_name\t$block_read_start\t$block_read_stop\t$block_number\n";
    
                                (my $gap_type, my $gap_size) = split(/\-/, $indels_by_read_pos{$read_position});
                                $gap_size = 1 if(!$gap_size);
                                $gaps += $gap_size;
    
    #                            $block_read_stop++;
     #                           $block_ref_stop++;
    
                                if($gap_type eq "D")
                                {
                                    $block_ref_stop += $gap_size + 1;
                                }
                                else
                                {
                                    $block_read_stop += $gap_size + 1;
                                }
    
                                $block_ref_start = $block_ref_stop;
                                $block_read_start = $block_read_stop;
                            }
                            else
                            {
                                $match++;
                                
                                $block_read_stop++;
                                $block_ref_stop++;
                            }
                        }
                        
                        ## Block ends here ##
                        $block_number++;
                        print BLOCKS "$ref_name\t$block_ref_start\t$block_ref_stop\t$align_strand\t$read_name\t$block_read_start\t$block_read_stop\t$block_number\n";                        
                    }
                    elsif($align_strand eq '-')
                    {
                        $block_read_start = $read_start;
                        $block_read_stop = $read_start;
                        
                        $block_ref_start = $ref_stop;
                        $block_ref_stop = $ref_stop;
                        
                        for(my $read_position = $read_start; $read_position <= $read_stop; $read_position++)
                        {
                            if($snps_by_read_pos{$read_position})
                            {
                                $mismatch++;
                            }
                            elsif($indels_by_read_pos{$read_position})
                            {
                                ## Block ends here ##
                                $block_number++;
                                $block_read_stop++;
#                                $block_ref_start--;
                                print BLOCKS "$ref_name\t$block_ref_start\t$block_ref_stop\t$align_strand\t$read_name\t$block_read_start\t$block_read_stop\t$block_number\n";
    
                                (my $gap_type, my $gap_size) = split(/\-/, $indels_by_read_pos{$read_position});
                                $gap_size = 1 if(!$gap_size);
                                $gaps += $gap_size;
    
                                if($gap_type eq "D")
                                {
                                    $block_ref_start -= ($gap_size + 1);
                                }
                                else
                                {
                                    $block_read_stop += $gap_size + 1;
                                }
    
                                $block_ref_stop = $block_ref_start;
                                $block_read_start = $block_read_stop + 1;
                            }
                            else
                            {
                                $match++;
                                
                                $block_read_stop++;
                                $block_ref_start--;
                            }
                        }
                        
                        ## Block ends here ##
                        $block_number++;
                        print BLOCKS "$ref_name\t$block_ref_start\t$block_ref_stop\t$align_strand\t$read_name\t$block_read_start\t$block_read_stop\t$block_number\n"; 

                    }

                }
                else
                {
                    ## Ungapped alignment, so print as contiguous block ##
                    print BLOCKS "$ref_name\t$ref_start\t$ref_stop\t$align_strand\t$read_name\t$read_start\t$read_stop\t1\n";
                }


                my $best_alignment = "$ref_name\t$ref_start\t$ref_stop\t$align_strand\t$read_name\t$read_start\t$read_stop\t$align_score";
                print ALIGNMENTS "$best_alignment\n";
                
                
            }
            
        }
        
        sub byAlignScoreDesc
        {
            my @temp = split(/\t/, $a);
            my $score_a = $temp[0];
            @temp = split(/\t/, $b);
            my $score_b = $temp[0];
            $score_b <=> $score_a;
        }
        


        close(ALIGNMENTS);
        close(SCOREDALIGNS);
        close(BLOCKS);

        print $Stats{'alignments'} . " alignments in $FileName\n";
        print $Stats{'aligned_reads'} . " unique reads had alignments\n";
        print $Stats{'best_alignments'} . " reads had a single best alignment\n"; 

        return(0);
}




################################################################################################
=head2 detect_variants - a subroutine to parse read alignments in Newbler pairt format
################################################################################################
=cut

sub detect_variants
{
	(my $FileName, my $output_dir, my $sample_name) = @_;

	my %Stats = ();

	$Stats{'insertions'} = $Stats{'deletions'} = $Stats{'substitutions'} = 0;

	# Alignment variables #
	my $read_name = my $read_start = my $read_stop = my $ref_name = my $ref_start = my $ref_stop = "";
	my $align_score = my $align_strand = "";
        my $current_alignment = "";

        ## Open output files ##

	open(INSERTIONS, ">$output_dir/$sample_name.insertions") or die "Can't open outfile: $!\n";	
	print INSERTIONS "ref_name\tref_start\tref_stop\tindel_type\tindel_size\tread_name\tread_start\tread_stop\talign_strand\tconf_score\tcontext\n";
	
	open(DELETIONS, ">$output_dir/$sample_name.deletions") or die "Can't open outfile: $!\n";	
	print DELETIONS "ref_name\tref_start\tref_stop\tindel_type\tindel_size\tread_name\tread_start\tread_stop\talign_strand\tconf_score\tcontext\n";	
	
	open(SUBSTITUTIONS, ">$output_dir/$sample_name.substitutions") or die "Can't open outfile: $!\n";			
	print SUBSTITUTIONS "ref_name\tref_position\tref_allele\tvar_allele\tread_name\tread_position\talign_strand\tcontext\tbase_quality\n";


	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

        ## Pre-parse the file to get the best alignments for each read ##
	$input = new FileHandle ($FileName);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter >= 1 && $line && substr($line, 0, 9) eq "ALIGNMENT")
		{
                        $Stats{'alignments'}++;
                
                        (my $temp, $align_score, $ref_name, $ref_start, $ref_stop, $align_strand, $read_name, $read_start, $read_stop) = split(/\t/, $_);

                }
		elsif($lineCounter >= 1 && $line && substr($line, 0, 11) eq "DISCREPANCY")
		{
                        (my $temp, my $discrep_type, my $discrep_read_position, my $discrep_genotype, my $discrep_ref_position, my $discrep_context) = split(/\s+/, $line);
                        
                        my $indel_type = my $indel_size = my $indel_ref_start = my $indel_ref_stop = my $indel_read_start = my $indel_read_stop = "";
                        
                        if($discrep_type =~ 'I' || $discrep_type =~ 'D')
                        {
                            if(length($discrep_type) > 1)
                            {
                                ($discrep_type, $indel_size) = split(/\-/, $discrep_type);
                                $indel_type = "INSERTION" if($discrep_type eq "I");
                                $indel_type = "DELETION" if($discrep_type eq "D");
                            }
                            else
                            {
                                $indel_type = "INSERTION" if($discrep_type eq "I");
                                $indel_type = "DELETION" if($discrep_type eq "D");
                                $indel_size = 1;
                            }
                            
                            if($align_strand eq "+")
                            {
                                
                            }
                            elsif($align_strand eq "-")
                            {
                                
                            }

                            if($indel_type eq "INSERTION")
                            {
                                $indel_ref_start = $indel_ref_stop = $discrep_ref_position;
                                $indel_read_start = $discrep_read_position;
                                $indel_read_stop = $indel_read_start + $indel_size - 1;
                                print INSERTIONS "$ref_name\t$indel_ref_start\t$indel_ref_stop\t$indel_type\t$indel_size\t$read_name\t$indel_read_start\t$indel_read_stop\t$align_strand\t$discrep_context\n";
                                $Stats{'insertions'}++;
                            }
                            elsif($indel_type eq "DELETION")
                            {
                                $indel_ref_start = $discrep_ref_position;
                                $indel_ref_stop = $discrep_ref_position + $indel_size - 1;
                                $indel_read_start = $indel_read_stop = $discrep_read_position;
                                print DELETIONS "$ref_name\t$indel_ref_start\t$indel_ref_stop\t$indel_type\t$indel_size\t$read_name\t$indel_read_start\t$indel_read_stop\t$align_strand\t$discrep_context\n";
                                $Stats{'deletions'}++;
                            }
                            
                        }
                        else
                        {
                            (my $variant_allele, my $qual_score) = split(/[\(\)]/, $discrep_genotype);
                            print SUBSTITUTIONS "$ref_name\t$discrep_ref_position\tN\t$variant_allele\t$read_name\t$discrep_read_position\t$align_strand\t$discrep_context\t$qual_score\n";
                            $Stats{'substitutions'}++;
                        }
                }
        }
        
        close($input);

        close(INSERTIONS);
        close(DELETIONS);
        close(SUBSTITUTIONS);

        print "$Stats{'insertions'} insertions detected\n";
        print "$Stats{'deletions'} deletions detected\n";
        print "$Stats{'substitutions'} substitutions detected\n";

        return(0);
}








=head1 AUTHOR

Dan Koboldt, C<< <dankoboldt at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-varscan-parsecrossmatch at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=VarScan>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc VarScan

You can also look for information at:

=over 4

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/VarScan>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/VarScan>

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=VarScan>

=item * Search CPAN

L<http://search.cpan.org/dist/VarScan>

=back

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2008 Dan Koboldt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of VarScan::ParseCrossMatch
