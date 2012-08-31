
package Genome::Model::Tools::Blat::FilterAlignments;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FilterAlignments.pm - Filter PSLx-format alignments, restricting them to target regions
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	10/20/2008 by D.K.
#	MODIFIED:	10/21/2008 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Blat::FilterAlignments {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		input_file	=> { is => 'Text', doc => "File containing best BLAT alignments in PSL/PSLX format" },
		target_headers	=> { is => 'Text', doc => "Exclude alignments outside regions in these FASTA headers", is_optional => 1},		
		output_file	=> { is => 'Text', doc => "Output file to contain filtered variants" },	
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Filter scored BLAT alignments"                 
}

sub help_synopsis {
    return <<EOS
This command parses a PSL-format BLAT output file, scores each alignment, and reports the best alignment for uniquely placed reads.
EXAMPLE:	gmt blat parse-alignments myBlatOutput.psl
Scored alignments for uniquely placed reads would be output to myBlatOutput.psl.best-alignments.txt
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
	my $input_file = $self->input_file;
	my $output_file = $self->output_file;

	my %TargetsByChrom = ();
	if($self->target_headers)
	{
		%TargetsByChrom = ParseTargetHeaders($self->target_headers);
	}

	## Open the outpput file ##

	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

	print "Parsing PSL file...\n";
	my $input = new FileHandle ($input_file);
	my $lineCounter = my $pslFormatCounter = my $notFilteredCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
	
		## Parse out only lines matching a BLAT-like result pattern ##
		
		if($line && $line=~/\d+\t\d+\d+\t\d+\t\d+\d+\t/)	
		{
			$pslFormatCounter++;
			my @lineContents = split(/\t/, $line);		
			my $numContents = @lineContents;
		
			my $chrom = $lineContents[14];
			my $chr_start = $lineContents[16];
			my $chr_stop = $lineContents[17];
		
			$chrom =~ s/[^0-9XYM]//g;
			
			my $include_flag = 0;
			if($TargetsByChrom{$chrom})
			{
				my @amplicons = split(/\n/, $TargetsByChrom{$chrom});
				my $num_genes = @amplicons;
				
				for(my $gCounter = 0; $gCounter < $num_genes; $gCounter++)
				{
					(my $amplicon_start, my $amplicon_stop, my $amplicon_name) = split(/\t/, $amplicons[$gCounter]);
					if($chr_stop >= $amplicon_start && $chr_start <= $amplicon_stop)
					{
						$include_flag++;
#						$AmpliconReads{$amplicon_name} .= "$read_name\n";
#						$NumAmpliconReads{$amplicon_name}++;
					}
				}
			}		
		
			if($include_flag)
			{
				print OUTFILE "$line\n";	
				$notFilteredCounter++;		
			}

		}
		elsif($lineCounter == 1)
		{
			print OUTFILE "$line\n";
		}
	}	
	
	close($input);
	
	print "$pslFormatCounter PSL-format lines parsed\n";
	print "$notFilteredCounter remained after filter\n";
	
	close(OUTFILE);
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


#############################################################
# ParseTargetHeaders - parse headers of gene targets
#
#############################################################

sub ParseTargetHeaders
{
	my $FileName = shift(@_);

	my %TargetsByChrom = ();
	my $numAmplicons = 0;

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($line && substr($line, 0, 1) eq ">")
		{
			my @lineContents = split(/\s+/, $line);
			my $numContents = @lineContents;
			
			my $amplicon_name = my $amplicon_chrom = my $amplicon_chr_start = my $amplicon_chr_stop = "";
			
			for(my $eCounter = 0; $eCounter < $numContents; $eCounter++)
			{
				if($eCounter == 0)
				{
					$amplicon_name = substr($lineContents[$eCounter], 1, 999);
				}
				if($lineContents[$eCounter] =~ "Chr")
				{
					my @temp = split(/\:/, $lineContents[$eCounter]);
					$amplicon_chrom = $temp[1];
					$amplicon_chrom =~ s/\,//;
				}
				if($lineContents[$eCounter] =~ "Coords")
				{
					($amplicon_chr_start, $amplicon_chr_stop) = split(/\-/, $lineContents[$eCounter + 1]);
					$amplicon_chr_stop =~ s/\,//;
				}
	
			}
	
			if($amplicon_name && $amplicon_chrom && $amplicon_chr_stop && $amplicon_chr_start)
			{
				$TargetsByChrom{$amplicon_chrom} .= "$amplicon_chr_start\t$amplicon_chr_stop\t$amplicon_name\n";
#				$AmpliconCoords{$amplicon_name} = "$amplicon_chrom\t$amplicon_chr_start\t$amplicon_chr_stop";
				$numAmplicons++;
			}
		}

	}	
	
	print "$numAmplicons target coordinate sets parsed from $FileName\n";
	
	return(%TargetsByChrom);
}


1;

