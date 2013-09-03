
package Genome::Model::Tools::Capture::FilterMafWithReview;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FilterMafWithReview - Retrieve, reformat, and annotate mutations from a somatic-variation group
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	02/09/2012 by D.K.
#	MODIFIED:	03/15/2012 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##

my %stats = ();


class Genome::Model::Tools::Capture::FilterMafWithReview {
	is => 'Genome::Model::Tools::Capture',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file		=> { is => 'Text', doc => "File of original MAF predictions" , is_optional => 0},
		review_file	=> { is => 'Text', doc => "File of manual review results: model, build, tumor, then 5 review columns" , is_optional => 0},
		output_file	=> { is => 'Text', doc => "Output file for filtered MAF" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Filters predicted mutations based on manual review"                 
}

sub help_synopsis {
    return <<EOS
This command filters predicted mutations based on manual review
EXAMPLE:	gmt capture filter-maf-with-review --maf-file myMaf.tsv --review-file Compiled-Review.tsv --output-file myMaf.filtered.tsv
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

	my $maf_file = $self->maf_file;
        my $output_file = $self->output_file;
        my $review_file = $self->review_file;

        my %review = load_review($review_file);
        my $maf_header = "";

	my $input = new FileHandle ($maf_file);
	my $lineCounter = 0;
	my %column_index = ();
        my %mutation_counts = ();
        
        open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
        print OUTFILE "#version 2.2\n";
        
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);
	
                
        
		if($lineCounter <= 2 && $line =~ "Hugo_Symbol")
		{
			$maf_header = $line;
			print OUTFILE "$maf_header\n";
			my $numContents = @lineContents;
			
			for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
			{
				if($lineContents[$colCounter])
				{
					$column_index{$lineContents[$colCounter]} = $colCounter;
				}
			}
			
#			$header_line = $line;
		}
		elsif($lineCounter > 1 && !%column_index)
		{
			die "No Header in MAF file!\n";
		}
		elsif(%column_index)
		{
			## Build a record for this line, assigning all values to respective fields ##
			
			my %record = ();

			foreach my $column_name (keys %column_index)
			{
				my $column_number = $column_index{$column_name};
				$record{$column_name} = $lineContents[$column_number];
			}			

			my @tempArray = split(/\-/, $record{'Tumor_Sample_Barcode'});
			my $tumor_sample = "TCGA" . "-" . $tempArray[1] . "-" . $tempArray[2] . "-" . $tempArray[3];
                        $mutation_counts{$tumor_sample}++;

                        my $chrom = $record{'Chromosome'};
                        my $chr_start = $record{'Start_position'};
                        my $chr_stop = $record{'End_position'};
                        my $variant_type = $record{'Variant_Type'};
                        my $ref = $record{'Reference_Allele'};
                        my $var = $record{'Tumor_Seq_Allele2'};
                        
                        $stats{'maf_' . $variant_type}++;
                        
                        if($variant_type ne "SNP")
                        {
                                $chr_start-- if($chr_start eq $chr_stop);
                        }
                        
                        my $variant_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
                        my $mutation_key = join("\t", $variant_key, $tumor_sample);                        
                        
                        my $review_code = "";
                        if($review{$mutation_key})
                        {
                                $stats{'maf_' . $variant_type . '_reviewed'}++;
                                $review_code = $review{$mutation_key};
                        }
                        elsif($review{$variant_key})
                        {
                                $stats{'maf_' . $variant_type . '_pos_reviewed'}++;
                                $review_code = $review{$variant_key};
                        }

                        if($review_code)
                        {
                                $stats{'maf_' . $variant_type . '_reviewed_' . $review_code}++;

                                if($review_code eq 'WT' || $review_code eq 'LQ' || $review_code eq 'LOH' || $review_code eq 'G' || $review_code eq 'D' || $review_code eq 'O')
                                {
                                        ## Remove a false positive ##
                                        $stats{'maf_' . $variant_type . '_removed'}++;
                                }
                                else
                                {
                                        print OUTFILE "$line\n";
                                        $stats{'maf_' . $variant_type . '_retained_pass_review'}++;
                                        $stats{'maf_' . $variant_type . '_retained'}++;
                                }

                        }
                        else
                        {
                                $stats{'maf_' . $variant_type . '_retained'}++;
                                print OUTFILE "$line\n";
                        }
                	
		}

	}

	close($input);
        close(OUTFILE);

        foreach my $tumor_sample (sort keys %mutation_counts)
        {
                $stats{'tumor_samples'}++;
                if(!$review{$tumor_sample})
                {
                        print "No review for $tumor_sample " . $mutation_counts{$tumor_sample} . " mutations\n";
                }
                else
                {
                        $stats{'tumor_samples_with_review'}++;
                }
        }

	
	foreach my $key (sort keys %stats)
	{
		print $stats{$key} . "\t" . $key . "\n";
	}
	
        ## Get SNV validation rate ##
        
#       my $snv_val_rate = ($stats{'snvs_reviewed_S'} + $stats{'snvs_reviewed_V'}) / ($stats{'snvs_reviewed_S'} + $stats{'snvs_reviewed_V'} + $stats{'snvs_reviewed_O'} + $stats{'snvs_reviewed_LQ'} + $stats{'snvs_reviewed_LOH'} + $stats{'snvs_reviewed_G'} + $stats{'snvs_reviewed_D'}) * 100;
#        $snv_val_rate = sprintf("%.2f", $snv_val_rate);
#        print "SNV validation rate: $snv_val_rate\%\n";

#        my $insertion_val_rate = ($stats{'insertions_reviewed_S'} + $stats{'insertions_reviewed_V'}) / ($stats{'insertions_reviewed_S'} + $stats{'insertions_reviewed_V'} + $stats{'insertions_reviewed_O'} + $stats{'insertions_reviewed_LQ'} + $stats{'insertions_reviewed_LOH'} + $stats{'insertions_reviewed_G'} + $stats{'insertions_reviewed_D'}) * 100;
#        $insertion_val_rate = sprintf("%.2f", $insertion_val_rate);
#        print "Insertion validation rate: $insertion_val_rate\%\n";
        
#        my $deletion_val_rate = ($stats{'deletions_reviewed_S'} + $stats{'deletions_reviewed_V'}) / ($stats{'deletions_reviewed_S'} + $stats{'deletions_reviewed_V'} + $stats{'deletions_reviewed_O'} + $stats{'deletions_reviewed_LQ'} + $stats{'deletions_reviewed_LOH'} + $stats{'deletions_reviewed_G'} + $stats{'deletions_reviewed_D'}) * 100;
#        $deletion_val_rate = sprintf("%.2f", $deletion_val_rate);
#        print "Deletion validation rate: $deletion_val_rate\%\n";

	return 1;
}



#############################################################
# parse_file - parses the file
#
#############################################################

sub load_review
{
	my $FileName = shift(@_);
        my %review = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	my $lines = "";
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

                my ($model_id, $build_id, $sample_name, $chrom, $chr_start, $chr_stop, $alleles, $alleles2, $review_code, $review_note) = split(/\t/, $line);
		
                if($alleles =~ '/')
                {
                        my $mutation_key = my $variant_key  = "";
                        my @temp = split(/\-/, $sample_name);
                        my $tumor_sample_name = join("-", "TCGA", $temp[1], $temp[2], $temp[3]);
                        my ($ref, $var) = split(/\//, $alleles);
                        if($ref eq '*' || $ref eq '0' || $ref eq '-' || length($var) > 1)
                        {
                                ## INSERTION ##
                                $ref = "-";
                                $chr_start++;
                                $chr_stop++;
                                $chr_start-- if($chr_stop eq $chr_start);
#                                $chr_stop = $chr_start;
                                $variant_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
                                $mutation_key = join("\t", $variant_key, $tumor_sample_name);
                                
#                                print "Saving $variant_key\n" if($chr_start =~ '12488705');
                                $stats{'insertions_reviewed'}++;
#                               $stats{'insertions_reviewed_' . $review_code}++;
                        }
                        elsif($var eq '*' || $var eq '0' || $var eq '-' || length($ref) > 1)
                        {
                                ## DELETION ##
                                $var = "-";
                                $chr_start++;
                                $variant_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
                                $mutation_key = join("\t", $variant_key, $tumor_sample_name);                                
                                $stats{'deletions_reviewed'}++;
#                                $stats{'deletions_reviewed_' . $review_code}++;
                        }
                        else
                        {
                                ## SNV ##
                                $chr_start = $chr_stop if($chr_stop > $chr_start);
                                $variant_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
                                $mutation_key = join("\t", $variant_key, $tumor_sample_name);
                                $stats{'snvs_reviewed'}++;
#                                $stats{'snvs_reviewed_' . $review_code}++;
                        }
                        
                        if($review_code ne "NA" && $review_code ne "A")
                        {
                                $review{$variant_key} = $review_code;
                                $review{$mutation_key} = $review_code;
                                $review{$tumor_sample_name} = 1;
                        }

                }
                else
                {
                        print "No var: $line\n";
                }
	}
	
	close($input);
	
	return(%review);
}

1;

