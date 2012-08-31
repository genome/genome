package Genome::Model::Tools::Vcf::CompareGenotypes;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# CompareGenotypes.pm -    Merge variants in a VCF file with annotation and dbSNP information
#               
#   AUTHOR:      Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#   CREATED:    03/13/2012 by D.K.
#   MODIFIED:   03/13/2012 by D.K.
#
#   NOTES:   
#         
#####################################################################################################################################


use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Vcf::CompareGenotypes {
    is => 'Command',                       

    has => [                                # specify the command's single-value properties (parameters) <--- 
    vcf_file1   => { is => 'Text', doc => "Sequence-based VCF file", is_input => 1},
    vcf_file2   => { is => 'Text', doc => "Validation or secondary sequence-based VCF file", is_input => 1},
    output_file   => { is => 'Text', doc => "Output file to receive results", is_input => 1, is_output => 1},
    verbose   => { is => 'Text', doc => "Print interim progress updates", is_optional => 1 },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merge variants in a VCF file with VEP and dbSNP annotation"                 
}

sub help_synopsis {
    return <<EOS
This command compares genotypes between two VCF files
EXAMPLE:	gmt vcf merge-with-annotation --vcf-file1 myVCF1.vcf --vcf-file2 myVCF2.vcf --output-file comparison.out ...
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
    my $vcf_file1 = $self->vcf_file1;
    my $vcf_file2 = $self->vcf_file2;
    my $output_file = $self->output_file;
    my $verbose = 1 if(defined($self->verbose));
    
    my %stats = ();
    

    ## Open the outfile ##
    if($output_file)
    {
        open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
        open(OUTPOSITIONS, ">$output_file.positions") or die "Can't open outfile: $!\n";
    }

    ## Reset variables for all positions seen all samples ##
    my %positions = ();
    my %samples = ();
    my %positions1 = my %positions2 = my %positions1filt = ();
    
    my %vcf1 = load_vcf($vcf_file1, 0);
    my %unfiltered_vcf1 = load_vcf($vcf_file1, 1);

    foreach my $gt_key (keys %vcf1)
    {
        my ($chrom, $position, $sample) = split(/\t/, $gt_key);
        my $pos_key = join("\t", $chrom, $position);
        $positions{$pos_key} = 1;
        $positions1{$pos_key} = 1;
        $samples{$sample} = 1;
    }

    foreach my $gt_key (keys %unfiltered_vcf1)
    {
        my ($chrom, $position, $sample) = split(/\t/, $gt_key);
        my $pos_key = join("\t", $chrom, $position);
        $positions1filt{$pos_key} = 1;
    }



    my %vcf2 = load_vcf($vcf_file2, 0);
   

    foreach my $gt_key (keys %vcf2)
    {
        my ($chrom, $position, $sample) = split(/\t/, $gt_key);
        my $pos_key = join("\t", $chrom, $position);
        $positions{$pos_key} = 1;
        $positions2{$pos_key} = 1;
        $samples{$sample} = 1;
    }

    foreach my $position_key (keys %positions)
    {
        $stats{'all_positions'}++;
        
        if($positions1{$position_key} && $positions2{$position_key})
        {
            $stats{'all_positions_both'}++;
            print OUTPOSITIONS join("\t", $position_key, "both") . "\n";
        }
        elsif($positions1{$position_key})
        {
            $stats{'all_positions_1only'}++;
            print OUTPOSITIONS join("\t", $position_key, "1only") . "\n";            
        }
        elsif($positions2{$position_key})
        {
            $stats{'all_positions_2only'}++;
            if($positions1filt{$position_key})
            {
                $stats{'all_positions_2only_1filt'}++;
            }
            print OUTPOSITIONS join("\t", $position_key, "2only") . "\n";
        }
        
        foreach my $sample_name (keys %samples)
        {
            my $gt_key = join("\t", $position_key, $sample_name);
#            print OUTFILE join("\t", $gt_key, "") . "\n";
            
            if($vcf1{$gt_key} && $vcf2{$gt_key})
            {
                $stats{'all_vcf_both'}++;
                
                my ($ref1, $var1, $genotype1) = split(/\t/, $vcf1{$gt_key});
                my ($ref2, $var2, $genotype2) = split(/\t/, $vcf2{$gt_key});
                if($ref1 ne $ref2)
                {
                    warn "Warning: $ref1 does not equal $ref2\n";
                }
                else
                {
                    my $gt1 = convert_genotype($ref1, $var1, $genotype1);
                    my $gt2 = convert_genotype($ref2, $var2, $genotype2);
                    
                    if($gt1 eq $gt2)
                    {
                        $stats{'all_vcf_both_match'}++;
                        if(is_variant($genotype2))
                        {
                            ## Variant match ##
                            $stats{'all_vcf_both_match_variant'}++;
                            print OUTFILE join("\t", $gt_key, "both-match-variant") . "\n"; 
                        }
                        else
                        {
                            ## Reference match ##
                            $stats{'all_vcf_both_match_reference'}++;
                            print OUTFILE join("\t", $gt_key, "both-match-reference") . "\n"; 
                        }
                    }
                    else
                    {
                        $stats{'all_vcf_both_mismatch'}++;
                        print join("\t", "Mismatch", $position_key, $sample_name, $vcf1{$gt_key}, $vcf2{$gt_key}) . "\n";
                        if(is_variant($genotype1) && is_reference($genotype2))
                        {
                            ## Apparent false positive ##
                            $stats{'all_vcf_both_mismatch_variant1'}++;
                            print OUTFILE join("\t", $gt_key, "both-mismatch-variant1", $vcf1{$gt_key}, $vcf2{$gt_key}) . "\n"; 
                        }
                        elsif(is_reference($genotype1) && is_variant($genotype2))
                        {
                            ## Apparent false negative ##
                            $stats{'all_vcf_both_mismatch_variant2'}++;
                            print OUTFILE join("\t", $gt_key, "both-mismatch-variant2", $vcf1{$gt_key}, $vcf2{$gt_key}) . "\n"; 
                        }
                        else
                        {
                            $stats{'all_vcf_both_mismatch_variantdiscrep'}++;
                            if(is_heterozygous($genotype1) && is_homozygous($genotype2))
                            {
                                $stats{'all_vcf_both_mismatch_variantdiscrep_het_hom'}++;
                                print OUTFILE join("\t", $gt_key, "both-mismatch-discrep-het-hom", $vcf1{$gt_key}, $vcf2{$gt_key}) . "\n"; 
                            }
                            elsif(is_homozygous($genotype1) && is_heterozygous($genotype2))
                            {
                                $stats{'all_vcf_both_mismatch_variantdiscrep_hom_het'}++;
                                print OUTFILE join("\t", $gt_key, "both-mismatch-discrep-hom-het", $vcf1{$gt_key}, $vcf2{$gt_key}) . "\n"; 
                            }
                            else
                            {
                                $stats{'all_vcf_both_mismatch_variantdiscrep_het_het'}++;
                                print OUTFILE join("\t", $gt_key, "both-mismatch-discrep-het-het", $vcf1{$gt_key}, $vcf2{$gt_key}) . "\n"; 
                            }
                        }

                    }
                }
            }
            elsif($vcf1{$gt_key})
            {
                $stats{'all_vcf_1only'}++;
                my ($ref1, $var1, $genotype1) = split(/\t/, $vcf1{$gt_key});
                my $gt1 = convert_genotype($ref1, $var1, $genotype1);
                if(is_reference($genotype1))
                {
                    $stats{'all_vcf_1only_reference'}++;
                }
                else
                {
                    $stats{'all_vcf_1only_variant'}++;

                    if(is_heterozygous($genotype1))
                    {
                        $stats{'all_vcf_1only_variant_het'}++;
                        print OUTFILE join("\t", $gt_key, "1only-variant-het", $vcf1{$gt_key}) . "\n"; 
                    }
                    elsif(is_homozygous($genotype1))
                    {
                        $stats{'all_vcf_1only_variant_hom'}++;
                        print OUTFILE join("\t", $gt_key, "1only-variant-hom", $vcf1{$gt_key}) . "\n"; 
                    }
                    else
                    {
                        die "Unable to classify $genotype1 as reference, het or hom\n";
                    }
                }

            }
            elsif($vcf2{$gt_key})
            {
                $stats{'all_vcf_2only'}++;
                my ($ref2, $var2, $genotype2) = split(/\t/, $vcf2{$gt_key});
                my $gt2 = convert_genotype($ref2, $var2, $genotype2);
                
                if(is_variant($genotype2))
                {
                    $stats{'all_vcf_2only_variants'}++;

                    ## See if we filtered it ##
                    
                    if($unfiltered_vcf1{$gt_key})
                    {
                        $stats{'all_vcf_2only_variants_filtered_in_vcf1'}++;
                        print OUTFILE join("\t", $gt_key, "2only-variant-filt", $unfiltered_vcf1{$gt_key}, $vcf2{$gt_key}) . "\n"; 
                    }                    
                    elsif(is_heterozygous($genotype2))
                    {
                        $stats{'all_vcf_2only_variants_het'}++;
                        print OUTFILE join("\t", $gt_key, "2only-variant-het", $vcf2{$gt_key}) . "\n"; 
                    }
                    elsif(is_homozygous($genotype2))
                    {
                        $stats{'all_vcf_2only_variants_hom'}++;
                        print OUTFILE join("\t", $gt_key, "2only-variant-hom", $vcf2{$gt_key}) . "\n"; 
                    }
                    
                    
#                    print join("\t", "FalseNegative", $position_key, $sample_name, $ref2, $var2, $genotype2) . "\n";
                }
                else
                {
                    $stats{'all_vcf_2only_reference'}++;
                }
            }
        }
    }
    
    if($output_file)
    {
        close(OUTFILE);
        close(OUTPOSITIONS);
    }

    
#    print $stats{'num_positions'} . " positions\n";
#    print $stats{'all_vcf_1only'} . " genotypes unique to file 1\n";
#    print $stats{'all_vcf_2only'} . " genotypes unique to file 2\n";
#    print $stats{'all_vcf_both'} . " genotypes available in both files\n";
#    print $stats{'all_vcf_both_match'} . " genotypes match\n";
#    print $stats{'all_vcf_both_match_reference'} . " reference genotypes match\n";
#    print $stats{'all_vcf_both_match_variant'} . " variant genotypes match\n";
#    print $stats{'all_vcf_both_mismatch'} . " genotypes do not match\n";
    
    foreach my $key (sort keys %stats)
    {
        print $stats{$key} . "\t" . $key . "\n";
    }
    
    ## Calculate overall rates ##
    
    my $sample_name = "all";
    
    my $false_positive_rate = $stats{$sample_name . '_vcf_both_mismatch_variant1'} / ($stats{$sample_name . '_vcf_both_mismatch_variant1'} + $stats{$sample_name . '_vcf_both_match_variant'}) * 100;
    my $false_negative_rate = ($stats{'all_vcf_2only_variants'} + $stats{$sample_name . '_vcf_both_mismatch_variant2'}) / ($stats{'all_vcf_2only_variants'} + $stats{$sample_name . '_vcf_both_mismatch_variant2'} + $stats{$sample_name . '_vcf_both_match_variant'}) * 100;
    my $reported_false_negative_rate = $stats{$sample_name . '_vcf_both_mismatch_variant2'} / ($stats{$sample_name . '_vcf_both_mismatch_variant2'} + $stats{$sample_name . '_vcf_both_match_variant'}) * 100;

    my $true_positive_rate = $stats{$sample_name . '_vcf_both_match_variant'} / ($stats{$sample_name . '_vcf_both_match_variant'} + $stats{$sample_name . '_vcf_both_mismatch_variant1'}) * 100;    
    my $true_negative_rate = $stats{$sample_name . '_vcf_both_match_reference'} / ($stats{$sample_name . '_vcf_both_match_reference'} + $stats{$sample_name . '_vcf_both_mismatch_variant2'}) * 100;    

    print "SAMPLE $sample_name\n";
    print join("\t", "False positive rate:", $false_positive_rate) . "\n";
    print join("\t", "False negative rate (overall):", $false_negative_rate) . "\n";
    print join("\t", "False negative rate (reported):", $reported_false_negative_rate) . "\n";
    
    print join("\t", "True positive rate:", $true_positive_rate) . "\n";    
    print join("\t", "True negative rate:", $true_negative_rate) . "\n";
    
    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub is_reference
{
    my $vcf_gt = shift(@_);
    $vcf_gt =~ s/\|/\//;
    return(1) if($vcf_gt eq '0/0' || $vcf_gt eq '0|0');
    my ($a1, $a2) = split(/\//, $vcf_gt);
    return(0) if($a1 || $a2);
    return(1);
}


################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub is_variant
{
    my $vcf_gt = shift(@_);
    $vcf_gt =~ s/\|/\//;
    return(0) if($vcf_gt eq '0/0' || $vcf_gt eq '0|0');
    my ($a1, $a2) = split(/\//, $vcf_gt);
    return(1) if($a1 || $a2);
    return(0);
}



################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub is_heterozygous
{
    my $vcf_gt = shift(@_);
    $vcf_gt =~ s/\|/\//;
    my ($a1, $a2) = split(/\//, $vcf_gt);
    my @num = split(/\//, $vcf_gt);
    my $num_num = @num;
    warn "Can't convert $vcf_gt\n" if($num_num < 2);
    return(1) if($a1 ne $a2);
    return(0);
}

################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub is_homozygous
{
    my $vcf_gt = shift(@_);
    $vcf_gt =~ s/\|/\//;
    my ($a1, $a2) = split(/\//, $vcf_gt);
    return(1) if($a1 eq $a2);
    return(0);
}


################################################################################################
# load_vcf(file, 0) - load VCF but ignore genotypes with a non-pass filter status
# load_vcf(file, 1) - load VCF even if genotypes have non-pass filter status
################################################################################################

sub load_vcf
{
    my $vcf_file = shift(@_);
    my $ignore_filter = shift(@_);
    
    my %genotypes = ();
    
    my $input = new FileHandle ($vcf_file);
    my $lineCounter = 0;

    my @sample_names = ();
    my $num_samples = 0;
    
    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;

        my @lineContents = split(/\t/, $line);
        my $numContents = @lineContents;
        
        if(substr($line, 0, 1) eq '#')
        {
            if(substr($line, 0, 6) eq '#CHROM')
            {
                for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
                {
                    $sample_names[$num_samples] = $lineContents[$colCounter];
                    $num_samples++;
                }
                
#                print OUTFILE "\n";
                
                warn "$num_samples samples\n";
            }
        }
        else
        {
            my ($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format) = split(/\t/, $line);
            
            ## Only take sites passing filters ##
            if($ignore_filter || $filter eq '.' || uc($filter) eq 'PASS')
            {                
                ## Only take sites with a variant allele #
                if($ref && $var)
                {
                    my %gt_columns = ();
                    my @formatContents = split(/\:/, $format);
                    my $numContents2 = @formatContents;
                    for(my $colCounter2 = 0; $colCounter2 < $numContents2; $colCounter2++)
                    {
                        my $field_name = $formatContents[$colCounter2];
                        $gt_columns{$field_name} = $colCounter2;
                    }
                    
                    
                    my $key = join("\t", $chrom, $position, $position, $ref, $var);
                    
                    my $num_gt_ref = my $num_gt_het = my $num_gt_hom = my $num_gt_missing = 0;
                    my $genotypes = "";
                    my $sampleOffset = 0;
                    for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
                    {
                        my $sample_name = $sample_names[$sampleOffset];
                        my ($vcf_genotype) = split(/\:/, $lineContents[$colCounter]);
                        my @vcfContents = split(/\:/, $lineContents[$colCounter]);
                        ## Get filter status ##
                        my $gt_filter = ".";
                        my $filter_column = $gt_columns{'FT'};
                        if($filter_column)
                        {
                            $gt_filter = $vcfContents[$filter_column] if($filter_column && $vcfContents[$filter_column]);
                        }

                        if($vcf_genotype ne '.' && ($ignore_filter || $gt_filter eq "." || uc($gt_filter) eq "PASS"))
                        {
                            my $genotype = convert_genotype($ref, $var, $vcf_genotype);
                            
                            my $genotype_key = join("\t", $chrom, $position, $sample_name);
    #                        $genotypes{$genotype_key} = $genotype;
                            $genotypes{$genotype_key} = join("\t", $ref, $var, $vcf_genotype, $lineContents[$colCounter]);
                            
                            if($genotype eq $ref . $ref)
                            {
                                $num_gt_ref++;
                            }
                            elsif($genotype eq $ref . $var)
                            {
                                $num_gt_het++;
                            }
                            elsif($genotype eq $var . $var)
                            {
                                $num_gt_hom++;   
                            }
                            else
                            {
                                $num_gt_missing++;
                            }                            
                        }
                        else
                        {
                            warn "Ignoring genotype $vcf_genotype with filter status $gt_filter\n";
                        }


                        $sampleOffset++;
                    }


#                    return(%genotypes) if($lineCounter > 20000);

                }
            }
        }
        

    }
    
    close($input);
    
    return(%genotypes);

}

################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub convert_genotype
{
	my ($ref, $var, $genotype) = @_;
	
	return("NN") if($genotype eq '.');
	$genotype =~ s/\|/\//;
	return($ref . $ref) if($genotype eq '0/0');
	return($ref . $var) if($genotype eq '0/1' || $genotype eq '1/0');
	return($var . $var) if($genotype eq '1/1');
	
	if($var =~ '\,')
	{
		my @vars = split(/\,/, $var);
		
		my ($gt1, $gt2) = split(/\//, $genotype);
		
		if($gt1 == 0)
		{
			$genotype = $ref;
		}
		else
		{
			$genotype = $vars[$gt1 - 1];
		}

		if($gt2 == 0)
		{
			$genotype .= $ref;
		}
		else
		{
			$genotype .= $vars[$gt2 - 1];
		}
	}
	else
	{
		warn "Unable to convert $ref $var $genotype\n";		
	}

	return($genotype);
}



1;

