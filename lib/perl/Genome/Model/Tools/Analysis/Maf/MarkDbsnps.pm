
package Genome::Model::Tools::Analysis::Maf::MarkDbsnps;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MarkDbsnps - Perform a proximity analysis on mutations in the MAF file.
#
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	08/24/2010 by D.K.
#	MODIFIED:	08/24/2010 by D.K.
#
#	NOTES:
#
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my %stats = ();
my $max_proximity = 0;

class Genome::Model::Tools::Analysis::Maf::MarkDbsnps {
    is => 'Command',

    has => [                                # specify the command's single-value properties (parameters) <---
                                            maf_file	=> { is => 'Text', doc => "Original MAF file", is_input => 1 },
                                            dbsnp_file	=> { is => 'Text', doc => "Tab-delimited dbSNP file in chrom, start, stop, ref, var, rs, valstatus order", is_input => 1 },
                                            output_file	=> { is => 'Text', doc => "Original MAF file", is_output => 1 },
                                            verbose		=> { is => 'Text', doc => "Print verbose output", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Marks dbSNP RS number and validation status in a MAF file"
}

sub help_synopsis {
    return <<EOS
        This command Fixes DNPs in the MAF file, compressing them into single events
EXAMPLE:	gt analysis maf fix-dnps --maf-file original.maf --output-file corrected.maf
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
    my $maf_file = $self->maf_file;
    my $dbsnp_file = $self->dbsnp_file;

    if(!(-e $maf_file))
    {
        die "Error: MAF file not found!\n";
    }

    my $output_file = $self->output_file;

    warn "Loading SNVs...\n";
#    my %snvs_by_position = load_snvs($self);
    


    ## Load dbSNPs ##
    my %dbsnp = ();
    my $input = new FileHandle ($dbsnp_file);
    my $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;

        my ($chrom, $chr_start, $chr_stop, $ref, $var, $rs_number, $val_status) = split(/\t/, $line);
        my $key = join("\t", $chrom, $chr_start);

	if(length($ref) > 1 || length($var) > 1 || $ref eq "-" || $var eq "-")
	{
	    ## Indel key ##
	    $key = join("\t", $chrom, $chr_start, $chr_start);
	}

#        if($snvs_by_position{$key})
 #       {
            ## Determine approved validation status ##

            my $tcga_val_status = "";
            #by1000genomes;by2Hit2Allele; byCluster; byFrequency; byHapMap; byOtherPop; bySubmitter; alternate_allele

            if($val_status =~ 'by-cluster')
            {
                $tcga_val_status .= ";" if($tcga_val_status);
                $tcga_val_status .= "byCluster";
            }
            if($val_status =~ 'by-frequency')
            {
                $tcga_val_status .= ";" if($tcga_val_status);
                $tcga_val_status .= "byFrequency";
            }
            if($val_status =~ 'by-2hit-2allele')
            {
                $tcga_val_status .= ";" if($tcga_val_status);
                $tcga_val_status .= "by2Hit2Allele";
            }
            if($val_status =~ 'by-1000genomes')
            {
                $tcga_val_status .= ";" if($tcga_val_status);
                $tcga_val_status .= "by1000genomes";
            }
            if($val_status =~ 'by-hapmap')
            {
                $tcga_val_status .= ";" if($tcga_val_status);
                $tcga_val_status .= "byHapmap";
            }
            if($val_status =~ 'by-submitter')
            {
                $tcga_val_status .= ";" if($tcga_val_status);
                $tcga_val_status .= "bySubmitter";
            }
	    
	    $tcga_val_status = "unknown" if(!$tcga_val_status);
            
            
            $dbsnp{$key} = join("\t", $ref, $var, $rs_number, $tcga_val_status);			
#        }



    }

    close($input);


    ## Column index for fields in MAF file ##

    my %column_index = ();
    my @columns = ();

    my %snp_is_dnp = ();

    ## open outfile ##

    open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

    $input = new FileHandle ($maf_file);
    $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;

        my @lineContents = split(/\t/, $line);

        if($lineCounter <= 2 && $line =~ "Chrom")
        {

            my $numContents = @lineContents;

            for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
            {
                if($lineContents[$colCounter])
                {
                    $column_index{$lineContents[$colCounter]} = $colCounter;
                }
            }

            foreach my $column (keys %column_index)
            {
                ## Print out the columns as parsed ##
                #print "$column_index{$column}\t$column\n";
                $columns[$column_index{$column}] = $column;	## Save the column order ##
            }

            print OUTFILE "$line\n";
        }
        elsif($lineCounter < 2)
        {

        }
        elsif($lineCounter > 2 && !@columns)
        {
            die "No Header in MAF file!\n";
        }
        elsif($lineCounter > 2 && @columns)
        {
            $stats{'num_mutations'}++;

            ## Build a record for this line, assigning all values to respective fields ##

            my %record = ();

            foreach my $column_name (keys %column_index)
            {
                my $column_number = $column_index{$column_name};
                $record{$column_name} = $lineContents[$column_number];
            }


            ## Here's how to parse out information for this record ##

            my $chromosome = $record{'Chromosome'};
            my $chr_start = $record{'Start_position'};
	    my $chr_stop = $record{'End_position'};
            my $position = $record{'Start_position'};

            my $tumor_sample = $record{'Tumor_Sample_Barcode'};
            my $variant_type = $record{'Variant_Type'};
            my $ref = $record{'Reference_Allele'};
            my $var = $record{'Tumor_Seq_Allele2'};
            $var = $record{'Tumor_Seq_Allele1'} if($ref eq $var);

	    my $dbsnp_hit = "";

	    $stats{'num_mutations'}++;
            $stats{'num_mutations_' . $variant_type}++;
	    
            if($variant_type eq "SNP" or $variant_type eq "DNP")
            {

		for(my $position = $chr_start; $position <= $chr_stop; $position++)
		{
		    my $key = join("\t", $chromosome, $position);
    
		    if($dbsnp{$key})
		    {
			$dbsnp_hit = $dbsnp{$key};
			if($variant_type eq "DNP")
			{
			    if($dbsnp_hit)
			    {
				## DNP has a second dbSNP hit, so let's combined results ##
				my ($first_rs, $first_status) = split(/\t/, $dbsnp_hit);
				$dbsnp_hit = $dbsnp{$key} if($first_status eq "unknown");
			    }
			}

		    }		    
		}

            }
	    elsif($variant_type eq "INS" || $variant_type eq "DEL")
	    {
		my $indel_size = 0;
		if($variant_type eq "INS")
		{
		    $indel_size = length($var);
		}
		else
		{
		    $indel_size = length($ref);
		}

		for(my $position = $chr_start - $indel_size - 1; $position <= $chr_stop + $indel_size + 1; $position++)
		{
		    my $key = join("\t", $chromosome, $position, $position);
    
		    if($dbsnp{$key})
		    {
			my ($dbsnp_ref, $dbsnp_var, $dbsnp_rs, $dbsnp_status) = split(/\t/, $dbsnp{$key});
			my $dbsnp_indel_type = "DEL";
			my $dbsnp_indel_size = length($dbsnp_ref);
			if(length($dbsnp_var) > length($dbsnp_ref) || $dbsnp_ref eq "-")
			{
			    $dbsnp_indel_type = "INS";
			    $dbsnp_indel_size = length($dbsnp_var);
			}
			
			if($dbsnp_indel_type eq $variant_type)
			{
			    my $size_diff = abs($dbsnp_indel_size - $indel_size);
			    if($size_diff == 0)
			    {
	    			$dbsnp_hit = $dbsnp{$key};				
			    }
			    ## Variable length tri-nucleotide repeat ##
			    elsif($variant_type eq "INS" && !($indel_size % 3) && !($dbsnp_indel_size % 3) && ($dbsnp_var =~ $var || $var =~ $dbsnp_var))
			    {
				$dbsnp_hit = $dbsnp{$key};
			    }
			    elsif($variant_type eq "DEL" && !($indel_size % 3) && !($dbsnp_indel_size % 3) && ($dbsnp_ref =~ $ref || $ref =~ $dbsnp_ref))
			    {
				$dbsnp_hit = $dbsnp{$key};
			    }
			    ## Variable length di-nucleotide repeat ##
			    elsif($variant_type eq "INS" && !($indel_size % 2) && !($dbsnp_indel_size % 2) && ($dbsnp_var =~ $var || $var =~ $dbsnp_var))
			    {
				$dbsnp_hit = $dbsnp{$key};
			    }
			    elsif($variant_type eq "DEL" && !($indel_size % 2) && !($dbsnp_indel_size % 2) && ($dbsnp_ref =~ $ref || $ref =~ $dbsnp_ref))
			    {
				$dbsnp_hit = $dbsnp{$key};
			    }
			    else
			    {
#				print "Same type different size indel: $chromosome $chr_start $chr_stop $ref $var $indel_size != $dbsnp_indel_size $dbsnp{$key}\n";				
			    }
			}

		    }		    
		}
		

	    }
	
	    ## Check to see if we have a matching dbSNP ##

	    if($dbsnp_hit)
	    {
		$stats{'num_mutations_' . $variant_type . '_dbsnp'}++;
		
		my ($dbsnp_allele1, $dbsnp_allele2, $dbsnp_rs, $dbsnp_val_status) = split(/\t/, $dbsnp_hit);

		if($dbsnp_val_status =~ '1000genomes')
		{
		    $stats{'num_mutations_' . $variant_type . '_dbsnp_1000g'}++;
		}
		elsif($dbsnp_val_status =~ 'Hapmap')
		{
			$stats{'num_mutations_' . $variant_type . '_dbsnp_HapMap'}++;
		}



		## Make a new maf line ##
		
		my $new_maf_line = "";
		
		foreach my $column_name (@columns)
		{
		    my $value = $record{$column_name};
		    
		    $value = $dbsnp_rs if($column_name eq "dbSNP_RS");
		    if($column_name eq "dbSNP_Val_Status"){
			$value = $dbsnp_val_status;
		    }
		    

		    $new_maf_line .= "\t" if($new_maf_line);
		    $new_maf_line .= $value;
		}
		
		print OUTFILE "$new_maf_line\n";						
	    }
	    else
	    {
		print OUTFILE "$line\n";
	    }


        }

    }

    close($input);

    close(OUTFILE);

    foreach my $key (sort keys %stats)
    {
	print $stats{$key} . "\t$key\n";
    }

    return(1);
}




################################################################################################
# Load all SNV positions
#
################################################################################################

sub load_snvs
{                               # replace with real execution logic.
    my $self = shift;

    ## Get required parameters ##
    my $maf_file = $self->maf_file;

    if(!(-e $maf_file))
    {
        die "Error: MAF file not found!\n";
    }


    my %snvs = ();

    ## Column index for fields in MAF file ##

    my %column_index = ();
    my @columns = ();


    my $input = new FileHandle ($maf_file);
    my $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;

        my @lineContents = split(/\t/, $line);

        if($lineCounter <= 2 && $line =~ "Chrom")
        {

            my $numContents = @lineContents;

            for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
            {
                if($lineContents[$colCounter])
                {
                    $column_index{$lineContents[$colCounter]} = $colCounter;
                }
            }

            foreach my $column (keys %column_index)
            {
                ## Print out the columns as parsed ##
                #print "$column_index{$column}\t$column\n";
                $columns[$column_index{$column}] = $column;	## Save the column order ##
            }
        }
        elsif($lineCounter < 2)
        {

        }
        elsif($lineCounter > 2 && !@columns)
        {
            die "No Header in MAF file!\n";
        }
        elsif($lineCounter > 2 && @columns)
        {

            ## Build a record for this line, assigning all values to respective fields ##

            my %record = ();

            foreach my $column_name (keys %column_index)
            {
                my $column_number = $column_index{$column_name};
                $record{$column_name} = $lineContents[$column_number];
            }


            ## Here's how to parse out information for this record ##

            my $chromosome = $record{'Chromosome'};
            my $position = $record{'Start_Position'};
            my $tumor_sample = $record{'Tumor_Sample_Barcode'};
            my $variant_type = $record{'Variant_Type'};

            if($variant_type eq "SNP")
            {
                my $key = join("\t", $chromosome, $position);
                $snvs{$key} = 1;
            }

        }

    }

    close($input);

    return(%snvs);

}




1;

