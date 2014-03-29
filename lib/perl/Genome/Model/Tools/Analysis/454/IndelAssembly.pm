
package Genome::Model::Tools::Analysis::454::IndelAssembly;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# IndelAssembly - Load 454 reads from a sample-SFF tab-delimited file
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/25/2009 by D.K.
#	MODIFIED:	02/25/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::454::IndelAssembly {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		indel_calls	=> { is => 'Text', doc => "File of called indels e.g. tumor.Insertions.tsv" },
		read_indels	=> { is => 'Text', doc => "File of read indel events e.g. tumor.insertions" },
                sff_file	=> { is => 'Text', doc => "Path to SFF file" },
	],
};

#, is_optional => 1

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Run Newbler assembly of indel-supporting 454 reads"                 
}

sub help_synopsis {
    return <<EOS
This command retrieves the locations of unplaced reads for a given genome model
EXAMPLE:	gmt bowtie --query-file s_1_sequence.fastq --output-file s_1_sequence.Hs36.bowtie
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
	my $indel_calls = $self->indel_calls;
	my $read_indels = $self->read_indels;
	my $sff_file = $self->sff_file;

#	my $output_file = $self->output_file;

	if(!(-e $indel_calls))
	{
		die "Error: Indel calls not found!\n";
	}

	if(!(-e $read_indels))
	{
		die "Error: Read indels not found!\n";
	}


	## Define Bowtie Reference (default to Hs36)

#	my $reference = "/gscmnt/sata194/info/sralign/dkoboldt/human_refseq/Hs36_1c_dkoboldt.bowtie";

	my %indel_calls = ();
 
 	my $input = new FileHandle ($indel_calls);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter > 1)
		{
			(my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size) = split(/\t/, $line);
#			print "$chrom\t$chr_start\t$indel_type\t$indel_size\n";
			my $indel_key = "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size";
			$indel_calls{$indel_key} = $line;
		}
	}
	
	close($input);
 

	## Get the offending reads ## 

	my %indel_reads = ();

 	$input = new FileHandle ($read_indels);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter > 1)
		{
			(my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $read_name) = split(/\t/, $line);
			my $indel_key = "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size";
			$indel_reads{$indel_key} .= "$read_name\n";
		}
	}
	
	close($input);
 
 
	foreach my $indel_key (keys %indel_calls)
	{
		(my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size) = split(/\t/, $indel_key);
		
		if($indel_reads{$indel_key})
		{
			print "$indel_key\t$indel_reads{$indel_key}\n";

			## Build an output directory ##
			
			mkdir("indel_assembly");
			
			## Get the sff file ##
			open(OUTFILE, ">indel.reads");
			print OUTFILE $indel_reads{$indel_key};
			close(OUTFILE);
			
			## Get the sub-sff file ##
			
			system("sfffile -o indel.sff -i indel.reads $sff_file");

			## Run the runAssembly ##
			
			system("runAssembly -o indel_assembly indel.sff");

			## Get the number of assembled contigs ##
			
			my $num_contigs = `grep -c \">\" indel_assembly/454AllContigs.fna`;
			chomp($num_contigs);

			print "$num_contigs CONTIGS\n";

			if($num_contigs)
			{
				exit(0);
			}
#			rmdir("indel_assembly");


		}
	}


 
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

