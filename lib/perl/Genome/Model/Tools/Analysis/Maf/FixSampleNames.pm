
package Genome::Model::Tools::Analysis::Maf::FixSampleNames;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FixSampleNames - Perform a proximity analysis on mutations in the MAF file.
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

class Genome::Model::Tools::Analysis::Maf::FixSampleNames {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Original MAF file", is_input => 1 },
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

	if(!(-e $maf_file))
	{
		die "Error: MAF file not found!\n";
	}

	my $output_file = $self->output_file;




	## Column index for fields in MAF file ##
	
	my %column_index = ();
	my @columns = ();


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
			
			my $tumor_sample = $record{'Tumor_Sample_Barcode'};
			my $normal_sample = $record{'Matched_Norm_Sample_Barcode'};

			my @tumorContents = split(/\-/, $tumor_sample);
			my $numTumorContents = @tumorContents;
			
			if($tumorContents[$numTumorContents - 1] eq "1")
			{
				my $new_tumor_sample = substr($tumor_sample, 0, length($tumor_sample) - 2);
				print "Fixing $tumor_sample\n";
				$line =~ s/$tumor_sample/$new_tumor_sample/;
			}

			my @normalContents = split(/\-/, $normal_sample);
			my $numNormalContents = @normalContents;

			if($normalContents[$numNormalContents - 1] eq "1")
			{
				my $new_normal_sample = substr($normal_sample, 0, length($normal_sample) - 2);
				print "Fixing $normal_sample\n";
				$line =~ s/$normal_sample/$new_normal_sample/;
			}

			## Print non-SNP lines ##
			print OUTFILE "$line\n";

			
		}

	}

	close($input);
	
	close(OUTFILE);

	print $stats{'num_mutations'} . " mutations in MAF file\n";
	print $stats{'num_snvs'} . " were SNVs\n";
	print $stats{'num_snvs_dbsnp'} . " were in dbSNP\n";
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
			my $position = $record{'Start_position'};
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

