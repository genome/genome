package Genome::Model::Tools::Analysis::Maf::FixDnps;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FixDnps - Perform a proximity analysis on mutations in the MAF file.
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

class Genome::Model::Tools::Analysis::Maf::FixDnps {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Original MAF file", is_input => 1 },
		output_file	=> { is => 'Text', doc => "Original MAF file", is_output => 1 },
		verbose		=> { is => 'Text', doc => "Print verbose output", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Fixes DNPs in the MAF file, compressing them into single events"                 
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

	if(!(-e $maf_file))
	{
		die "Error: MAF file not found!\n";
	}

	my $output_file = $self->output_file;

	warn "Loading SNVs...\n";
	my %snvs_by_position = load_snvs($self);

	## Column index for fields in MAF file ##
	
	my %column_index = ();
	my @columns = ();

	my %snp_is_dnp = ();

	## open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

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
			my $position = $record{'Start_position'};
			my $tumor_sample = $record{'Tumor_Sample_Barcode'};
			my $variant_type = $record{'Variant_Type'};
			my $ref_allele = $record{'Reference_Allele'};
			my $allele1 = $record{'Tumor_Seq_Allele1'};
			my $allele2 = $record{'Tumor_Seq_Allele2'};

			if($variant_type eq "SNP")
			{
				$stats{'num_snvs'}++;
				my $next_position = $position + 1;			
				my $key = join("\t", $chromosome, $position, $tumor_sample);
				my $next_key = join("\t", $chromosome, $next_position, $tumor_sample);
				
				if($snvs_by_position{$key} && $snvs_by_position{$next_key})
				{
					$stats{'num_dnps'}++;
					$stats{'num_snvs_were_dnps'}++;
					
					my $dnp_start = $position;
					my $dnp_stop = $position + 1;
#$record{'Match_Norm_Seq_Allele1'}, $record{'Match_Norm_Seq_Allele2'}, $record{'start_WU'}, $record{'stop_WU'}, $record{'reference_WU'}, $record{'variant_WU'}, $record{'c_position_WU'}, $record{'amino_acid_change_WU'}, $record{'ucsc_cons_WU'}, , $record{'domain_WU'}, $record{'all_domains_WU'}, $record{'deletion_substructures_WU'}					
					my ($next_ref_allele, $next_allele1, $next_allele2, $next_normal_allele1, $next_normal_allele2, $next_start_WU, $next_stop_WU, $next_reference_WU, $next_variant_WU, $next_c_position_WU, $next_amino_acid_change_WU) = split(/\t/, $snvs_by_position{$next_key});
					
					my $dnp_ref = $ref_allele . $next_ref_allele;
					my $dnp_allele1 = $allele1 . $next_allele1;
					my $dnp_allele2 = $allele2 . $next_allele2;
					
					my $dnp_normal_allele1 = $record{'Match_Norm_Seq_Allele1'} . $next_normal_allele1;
					my $dnp_normal_allele2 = $record{'Match_Norm_Seq_Allele2'} . $next_normal_allele2;
					my $dnp_start_WU = $record{'start_WU'};
					my $dnp_stop_WU = $next_stop_WU;
					my $dnp_reference_WU = $record{'reference_WU'} . $next_reference_WU if($record{'reference_WU'} && $next_reference_WU);
					my $dnp_variant_WU = $record{'variant_WU'} . $next_variant_WU if($record{'variant_WU'} && $next_variant_WU);
					my $dnp_c_position_WU = $record{'c_position_WU'} . "-" . $next_c_position_WU if($record{'c_position_WU'});
					my $dnp_amino_acid_change_WU = $record{'amino_acid_change_WU'} . "-" . $next_amino_acid_change_WU if($record{'amino_acid_change_WU'});
					
					$snp_is_dnp{$key} = 1;
					$snp_is_dnp{$next_key} = 1;

					## Make a new maf line ##
					
					my $new_maf_line = "";
					
					foreach my $column_name (@columns)
					{
						my $value = $record{$column_name};
						
						$value = $dnp_start if($column_name eq "Start_position");
						$value = $dnp_stop if($column_name eq "End_position");
						$value = "DNP" if($column_name eq "Variant_Type");
						$value = $dnp_ref if($column_name eq "Reference_Allele");
						$value = $dnp_allele1 if($column_name eq "Tumor_Seq_Allele1");
						$value = $dnp_allele2 if($column_name eq "Tumor_Seq_Allele2");
						$value = $dnp_normal_allele1 if($column_name eq "Match_Norm_Seq_Allele1");
						$value = $dnp_normal_allele2 if($column_name eq "Match_Norm_Seq_Allele2");
						$value = $dnp_start_WU if($column_name eq "start_WU");
						$value = $dnp_stop_WU if($column_name eq "stop_WU");
						$value = $dnp_reference_WU if($column_name eq "reference_WU");
						$value = $dnp_variant_WU if($column_name eq "variant_WU");
						$value = $dnp_c_position_WU if($column_name eq "c_position_WU");
						$value = $dnp_amino_acid_change_WU if($column_name eq "amino_acid_change_WU");

						$new_maf_line .= "\t" if($new_maf_line);
						$new_maf_line .= $value;
					}

					print OUTFILE "$new_maf_line\n";
				}
				elsif($snp_is_dnp{$key})
				{
					## Skip this because we already reported it ##
					$stats{'num_snvs_were_dnps'}++;
				}
				else
				{
					print OUTFILE "$line\n";
				}
			}
			else
			{
				## Print non-SNP lines ##
				print OUTFILE "$line\n";
			}

			
		}

	}

	close($input);
	
	close(OUTFILE);

	print $stats{'num_mutations'} . " mutations in MAF file\n";
	print $stats{'num_snvs'} . " were SNVs\n";
	print $stats{'num_snvs_were_dnps'} . " SNVs were compressed into " . $stats{'num_dnps'} . " DNPs\n";
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
			my $position = $record{'Start_position'};
			my $tumor_sample = $record{'Tumor_Sample_Barcode'};
			my $variant_type = $record{'Variant_Type'};
			my $ref_allele = $record{'Reference_Allele'};
			my $allele1 = $record{'Tumor_Seq_Allele1'};
			my $allele2 = $record{'Tumor_Seq_Allele2'};
			my $var = $allele2;
			$var = $allele1 if($var eq $ref_allele);

			if($variant_type eq "SNP")
			{
				my $key = join("\t", $chromosome, $position, $tumor_sample);			
				$snvs{$key} = join("\t", $ref_allele, $allele1, $allele2, $record{'Match_Norm_Seq_Allele1'}, $record{'Match_Norm_Seq_Allele2'}, $position, $position, $ref_allele, $var);				
				if($record{'c_position_WU'})
				{
					$snvs{$key} .= "\t" . join("\t", $record{'c_position_WU'}, $record{'amino_acid_change_WU'}, $record{'ucsc_cons_WU'}, , $record{'domain_WU'}, $record{'all_domains_WU'}, $record{'deletion_substructures_WU'});
				}
			}
			
		}

	}

	close($input);
	
	return(%snvs);

}




1;

