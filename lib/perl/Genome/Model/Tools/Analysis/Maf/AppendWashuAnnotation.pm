
package Genome::Model::Tools::Analysis::Maf::AppendWashuAnnotation;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# AppendWashuAnnotation - Perform a proximity analysis on mutations in the MAF file.
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

class Genome::Model::Tools::Analysis::Maf::AppendWashuAnnotation {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Original MAF file", is_input => 1 },
		annotation_file	=> { is => 'Text', doc => "Annotation file in native output format", is_input => 1 },		
		output_file	=> { is => 'Text', doc => "Output name for updated MAF file", is_output => 1 },
		verbose		=> { is => 'Text', doc => "Print verbose output", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Updates a MAF with validation status"                 
}

sub help_synopsis {
    return <<EOS
This command updates a MAF file with validation information
EXAMPLE:	gmt analysis maf update-with-validation --maf-file myMaf.tsv --validation-file myVal.tsv --output-file myMafUpdated.tsv
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

	my %annotation = load_annotation($self->annotation_file);

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
		my $numContents = @lineContents;
		
		if($lineCounter <= 2 && $line =~ "Chrom")
		{
			
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
			
			print OUTFILE "$line\t";
			print OUTFILE "chromosome_name_WU\tstart_WU\tstop_WU\treference_WU\tvariant_WU\ttype_WU\tgene_name_WU\ttranscript_name_WU\ttranscript_species_WU\ttranscript_source_WU\ttranscript_version_WU\tstrand_WU\ttranscript_status_WU\ttrv_type_WU\tc_position_WU\tamino_acid_change_WU\tucsc_cons_WU\tdomain_WU\tall_domains_WU\tdeletion_substructures_WU\ttranscript_error\n";
		}
		elsif($lineCounter < 2)
		{
			print OUTFILE "$line\n";
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
			
			my $chrom = $record{'Chromosome'};
			my $chr_start = $record{'Start_position'};
			$chr_start = $record{'Start_Position'} if(!$chr_start);
			my $chr_stop = $record{'End_position'};
			$chr_stop = $record{'End_Position'} if(!$chr_stop);			
			my $tumor_sample = $record{'Tumor_Sample_Barcode'};
			my $variant_type = $record{'Variant_Type'};
			my $ref = $record{'Reference_Allele'};			
			my $var = $record{'Tumor_Seq_Allele2'};
			$var = $record{'Tumor_Seq_Allele1'} if($var eq $ref);

			my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
			
			if($annotation{$key})
			{
				$stats{'num_mutations_with_annotation'}++;
				print OUTFILE "$line\t$annotation{$key}\n";								
			}
			else
			{
				warn "No annotation for $key\n";
				print OUTFILE "$line\n";				
			}



			
		}

	}

	close($input);
	
	close(OUTFILE);

	foreach my $key (sort keys %stats)
	{
		print $stats{$key} . "\t" . $key . "\n";
	}
	
	return(1);
}




################################################################################################
# Load all SNV positions
#
################################################################################################

sub load_annotation
{                               # replace with real execution logic.
	my $FileName = shift(@_);

	my %annotation = ();
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);

		my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
		
		$annotation{$key} = $line;
	}
	
	close($input);
	return(%annotation);

}




1;

