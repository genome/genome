
package Genome::Model::Tools::Analysis::Maf::MutationSummary;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MutationSummary - Perform a proximity analysis on mutations in the MAF file.
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

class Genome::Model::Tools::Analysis::Maf::MutationSummary {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Original MAF file" },
		verbose		=> { is => 'Text', doc => "Print verbose output", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Summarizes the mutations in a MAF file"                 
}

sub help_synopsis {
    return <<EOS
This command summarizes the mutations in a MAF file
EXAMPLE:	gt analysis maf mutation-summary --maf-file original.maf 
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


	my %mutation_class_counts = my %mutation_type_counts = ();

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
			
			my $hugo_name = $record{'Hugo_Symbol'};
			my $tumor_sample = $record{'Tumor_Sample_Barcode'};
			my $variant_class = $record{'Variant_Classification'};
			my $variant_type = $record{'Variant_Type'};

			$mutation_type_counts{$variant_type}++;
			$mutation_class_counts{$variant_class}++;
			
		}

	}

	close($input);	


	## Print a summary ##
	
	print $stats{'num_mutations'} . " mutations in the MAF file\n";
	
	print "\nMUTATIONS BY TYPE\n";
	foreach my $variant_type (sort keys %mutation_type_counts)
	{
		print $mutation_type_counts{$variant_type} . " " . $variant_type . "\n";
	}

	print "\nMUTATIONS BY CLASS\n";
	foreach my $variant_class (sort keys %mutation_class_counts)
	{
		print $mutation_class_counts{$variant_class} . " " . $variant_class . "\n";
	}
}




1;

