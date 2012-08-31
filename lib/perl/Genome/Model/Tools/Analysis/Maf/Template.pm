
package Genome::Model::Tools::Analysis::Maf::Template;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Template - Perform a proximity analysis on mutations in the MAF file.
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

class Genome::Model::Tools::Analysis::Maf::Template {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Original MAF file" },
		output_file	=> { is => 'Text', doc => "Output file for proximity report", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "A template for MAF analysis with parsing subroutines"                 
}

sub help_synopsis {
    return <<EOS
This is a stub template with basic functions for MAF parsing. DO NOT call directly.
EXAMPLE:	gt analysis maf template --maf-file original.maf --output-file output.tsv
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
		die "Error: MAF file 1 not found!\n";
	}

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
	}

	my %stats = ();

	## Column index for fields in MAF file ##
	
	my %column_index = ();
	my @columns = ();

	## Parse the MAF file ##
	
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
			warn "Parsing MAF header line...\n";

			if($self->output_file)
			{
				print OUTFILE "$line\n";
			}
			
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
		elsif($lineCounter > 2 && !@columns)
		{
			die "No Header in MAF file!\n";
		}
		elsif(@columns)
		{
			## Build a record for this line, assigning all values to respective fields ##
			
			my %record = ();

			foreach my $column_name (keys %column_index)
			{
				my $column_number = $column_index{$column_name};
				$record{$column_name} = $lineContents[$column_number];
			}			


			## Here's how to parse out information for this record ##
			
			my $hugo_name = $record{'Hugo_Symbol'};


		}

	}

	close($input);	
		



	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}






1;

