
package Genome::Model::Tools::Analysis::Solexa::BamToFastq;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# BamToFastq - 	Extract unmapped reads from a BAM file using SAMtools view with filters
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	11/23/2009 by D.K.
#	MODIFIED:	11/23/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Solexa::BamToFastq {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		bam_file	=> { is => 'Text', doc => "The BAM file to process", is_optional => 0 },
		output_file	=> { is => 'Text', doc => "File to receive FASTQ output", is_optional => 0 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Extracts unmapped reads from BAM files"                 
}

sub help_synopsis {
    return <<EOS
This command extracts reads from a BAM file and converts to FASTQ
EXAMPLE:	gmt analysis solexa bam-to-fastq

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
	my $bam_file = $self->bam_file;
	my $output_file = $self->output_file;

	## Get the SAM format ##
	
	my $sam = `samtools view $bam_file | cut --fields=1,2,10,11`;
	
	open(FASTQ, ">$output_file") or die "Can't open FASTQ file: $!\n";
	
	my @lines = split(/\n/, $sam);
	my $lineCounter = 0;
	
	foreach my $line (@lines)
	{
		$lineCounter++;
		
		(my $read_name, my $flag, my $read_seq, my $read_qual) = split(/\t/, $line);

		## Adjust flag values ##
		
		$flag -= 1024 if($flag >= 1024);
		$flag -= 512 if($flag >= 512);
		$flag -= 256 if($flag >= 256);

		## Determine read num ##
		
		my $read_num = 0;
		$read_num = 1 if($flag >= 64);
		$read_num = 2 if($flag >= 128);

		$read_name .= "/" . $read_num if($read_num);

		print FASTQ '@' . $read_name . "\n";
		print FASTQ $read_seq . "\n";
		print FASTQ '+' . $read_name . "\n";
		print FASTQ $read_qual . "\n";		
	}

	close(FASTQ);
	
	print "$lineCounter reads converted\n";
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub unmapped_to_fastq
{                               
	(my $unmapped_file, my $fastq_file, my $read_num) = @_;

	## Open the FASTQ file ##
	
	
	
	## Open the infile ##
	
	my $input = new FileHandle ($unmapped_file);
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		
		
		

	}

	close($input);
	
	close(FASTQ);
}


1;

