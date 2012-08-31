
package Genome::Model::Tools::Analysis::Solid::UnmappedFromBam;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# UnmappedFromBam - 	Extract unmapped reads from a BAM file using SAMtools view with filters
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

class Genome::Model::Tools::Analysis::Solid::UnmappedFromBam {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		bam_file	=> { is => 'Text', doc => "The BAM file to process", is_optional => 0 },
		output_dir	=> { is => 'Text', doc => "Directory to contain output files", is_optional => 0 },
		output_name	=> { is => 'Text', doc => "Name to use for output files, i.e. lane", is_optional => 0 },
		verbose	=> { is => 'Text', doc => "Print verbose output to STDOUT", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Extracts unmapped reads from BAM files"                 
}

sub help_synopsis {
    return <<EOS
This command searches for Illumina/Solid data using the database
EXAMPLE:	gmt analysis solexa unmapped-from-bam --bam-file myBam.bam --output-dir unmapped --output-name myBam
			==> unmapped/myBam_1_sequence.pe.unmapped
			==> unmapped/myBam_2_sequence.pe.unmapped
			==> unmapped/myBam_1_sequence.se.unmapped
			==> unmapped/myBam_2_sequence.se.unmapped
			==> unmapped/myBam_sequence.se.unmapped			
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
	my $sam_file = $bam_file . ".sam";
	my $output_dir = $self->output_dir;
	my $name = $self->output_name;


	## Extract the SAM file ##
	print "Extracting SAM format...\n";
	my $input_pipe_cmd = "samtools view -f 4 -F 1 $bam_file";

	## Build output file ##
	
	my $output_file = "$output_dir/$name.unmapped.reads";
#	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	open(OUTFASTQ, ">$output_file.csfastq") or die "Can't open outfile: $!\n";

	## Get the input ##
	
	my $lineCounter = 0;

	#my $input = new FileHandle ('samtools view -f 4 -F 1 ' . $bam_file . '|') or die "Can't open SAM piped input\n";
	open INPUT, "samtools view -f 4 -F 1 $bam_file |" or die "Can't open input: $!\n";
	
	while (<INPUT>)
	{
		chomp;
		my $line = $_;
		
		$lineCounter++;
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		my $read_name = $lineContents[0];
		my $read_seq = $lineContents[9];
		my $read_qual = $lineContents[10];
		my $cs_seq = my $cs_qual = "";
		for(my $colCounter = 11; $colCounter < $numContents; $colCounter++)
		{
			my $colValue = $lineContents[$colCounter];
			
			if(!$cs_seq || !$cs_qual)
			{
				if(substr($colValue, 0, 2) eq "CS")
				{
					$cs_seq = substr($colValue, 5, 99);
				}
				elsif(substr($colValue, 0, 2) eq "CQ")
				{
					$cs_qual = substr($colValue, 5, 99);
				}
			}
		}

#		print OUTFILE "$read_name\t$read_seq\t$read_qual\t$cs_seq\t$cs_qual\n";
		print OUTFASTQ "\@$read_name\n$cs_seq\n+\n$cs_qual\n";
		print $lineCounter .  "\t$read_name\t$read_seq\t$read_qual\t$cs_seq\t$cs_qual\n" if($self->verbose && !($lineCounter % 100000));
	}

#	close(OUTFILE);
	close(OUTFASTQ);
	
	print "$lineCounter unmapped reads extracted\n";	


	return 1;


	my $unmapped_pair_read1 = $output_dir . "/" . $name . "_1_sequence.pe.unmapped";
	my $unmapped_pair_read2 = $output_dir . "/" . $name . "_2_sequence.pe.unmapped";
	my $unmapped_frag_read1 = $output_dir . "/" . $name . "_1_sequence.unmapped";
	my $unmapped_frag_read2 = $output_dir . "/" . $name . "_2_sequence.unmapped";
	my $unmapped_single_read = $output_dir . "/" . $name . "_sequence.unmapped";
	
	## Isolate the unmapped pair 1 ##


	print "$unmapped_pair_read1\n";	
	system("samtools view -f 77 $bam_file | cut --fields=1,10,11 >$unmapped_pair_read1");
	unmapped_to_fastq($unmapped_pair_read1, "$unmapped_pair_read1.fastq", 1);

	print "$unmapped_pair_read2\n";
	system("samtools view -f 141 $bam_file | cut --fields=1,10,11 >$unmapped_pair_read2");
	unmapped_to_fastq($unmapped_pair_read2, "$unmapped_pair_read2.fastq", 2);

	print "$unmapped_frag_read1\n";
	system("samtools view -f 68 -F 8 $bam_file | cut --fields=1,10,11 >$unmapped_frag_read1");
	unmapped_to_fastq($unmapped_frag_read1, "$unmapped_frag_read1.fastq", 1);

	print "$unmapped_frag_read2\n";
	system("samtools view -f 132 -F 8 $bam_file | cut --fields=1,10,11 >$unmapped_frag_read2");
	unmapped_to_fastq($unmapped_frag_read2, "$unmapped_frag_read2.fastq", 2);

	print "$unmapped_single_read\n";
	system("samtools view -f 4 -F 1 $bam_file | cut --fields=1,10,11 >$unmapped_single_read");
	unmapped_to_fastq($unmapped_single_read, "$unmapped_single_read.fastq");
	
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
	
	open(FASTQ, ">$fastq_file") or die "Can't open FASTQ file: $!\n";
	
	## Open the infile ##
	
	my $input = new FileHandle ($unmapped_file);
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		
		(my $read_name, my $read_seq, my $read_qual) = split(/\t/, $line);
		
		$read_name .= "/" . $read_num if($read_num);

		print FASTQ '@' . $read_name . "\n";
		print FASTQ $read_seq . "\n";
		print FASTQ '+' . $read_name . "\n";
		print FASTQ $read_qual . "\n";
	}

	close($input);
	
	close(FASTQ);
}


1;

