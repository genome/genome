
package Genome::Model::Tools::Analysis::Indels::DbsnpConcord;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SearchRuns - Search the database for runs
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/01/2009 by D.K.
#	MODIFIED:	04/01/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Indels::DbsnpConcord {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		indels_file	=> { is => 'Text', doc => "File of indels to search" },
		position_diff	=> { is => 'Text', doc => "Max bp difference in indel positions [10]", is_optional => 1 },
		output_file	=> { is => 'Text', doc => "Output file to print indels along with matches", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Obtains reads from a flowcell_id in FastQ format"                 
}

sub help_synopsis {
    return <<EOS
This command determines concordance with dbSNP 129 indels
EXAMPLE:	gmt analysis indels dbsnp-concord --indels-file [myfile.tsv]
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
	my $indels_file = $self->indels_file;

	my %stats = ();

	## Load dbSNP indels ##
	print "Loading dbSNP indels...\n";
	my %solexa_indels = load_indels("/gscmnt/sata194/info/sralign/dkoboldt/SNPseek/snp129.txt.indels.varscan.format");

	my $position_diff = 10;
	$position_diff = $self->position_diff if($self->position_diff);

	my $size_diff = 1;

	my $output_file = $self->output_file if($self->output_file);
	
	if($output_file)
	{
		open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	}
	
	my $input = new FileHandle ($indels_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter >= 1)
		{
			(my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size) = split(/\t/, $line);
			$chrom =~ s/[^0-9XYM]//g;
			$chrom = "chr" . $chrom;
			
			if(!($chrom =~ 'chrom' || $chrom =~ 'chromosome' || $chrom =~ 'ref_name'))
			{
				my $chrom_key = $chrom . "\t" . substr($chr_start, 0, 3);
				
				$stats{'num_indels'}++;

				if($solexa_indels{$chrom_key})
				{
					my $dbsnp_support = "";
					
					my @chrom_indels = split(/\n/, $solexa_indels{$chrom_key});
					
					foreach my $chrom_indel (@chrom_indels)
					{
						(my $indel_chrom, my $indel_chr_start, my $indel_chr_stop, my $chrom_indel_type, my $chrom_indel_size, my $chrom_indel_coverage, my $chrom_indel_reads1, my $chrom_indel_reads2) = split(/\t/, $chrom_indel);
						
						if($indel_chrom eq $chrom)
						{
							if($chr_start <= ($indel_chr_stop + $position_diff) && $chr_stop >= ($indel_chr_start - $position_diff))
							{
								if($indel_type eq $chrom_indel_type)
								{
									if(abs($indel_size - $chrom_indel_size) <= $size_diff)
									{
										$stats{'identical_to_dbsnp'}++ if(($indel_chr_start eq $chr_start || $indel_chr_stop eq $chr_stop || $indel_chr_start eq $chr_stop) && $indel_type eq $chrom_indel_type && $indel_size eq $chrom_indel_size);

										$dbsnp_support .= "$indel_chrom\t$indel_chr_start\t$indel_chr_stop\t$chrom_indel_type\t$chrom_indel_size\t$chrom_indel_coverage\t$chrom_indel_reads1\t$chrom_indel_reads2\n";
									}
								}
	
							}
						}
					}
					
					if($dbsnp_support)
					{
						## Get the best support ##
						
						my @dbsnp_support = split(/\n/, $dbsnp_support);
#						@dbsnp_support = sort byReads2 @dbsnp_support;
						
						if($output_file)
						{
							print OUTFILE "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size\n$dbsnp_support[0]\n"; #$dbsnp_support\n";
						}
					
						$stats{'in_dbsnp'}++;					
					}

				}

			}
		}
	}
	
	close($input);
	close(OUTFILE) if($output_file);

	## Calculate pct ##
	
	$stats{'pct_in_dbsnp'} = $stats{'in_dbsnp'} / $stats{'num_indels'} * 100;
	$stats{'pct_in_dbsnp'} = sprintf("%.2f", $stats{'pct_in_dbsnp'}) . '%';

	print $stats{'num_indels'} . " indels loaded\n";
	print $stats{'in_dbsnp'} . " (" . $stats{'pct_in_dbsnp'} . ") closely matched known dbSNP indels\n";
	print $stats{'identical_to_dbsnp'} . " were identical in size and location\n";
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


sub byReads2
{
	my @temp = split(/\t/, $a);
	my $reads2_a = $temp[7];
	
	@temp = ();
	@temp = split(/\t/, $b);
	my $reads2_b = $temp[7];	

	$reads2_b <=> $reads2_a;
}


#############################################################
# Load indels
#
#############################################################

sub load_indels
{
	my $FileName = shift(@_);
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	my %indels = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		(my $chromosome, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $read_coverage, my $reads1, my $reads2) = split(/\t/, $line);

		if($chromosome && !($chromosome =~ "chrom" || $chromosome =~ 'ref_name'))
		{
			## Match chromosome format ##
			$chromosome =~ s/[^0-9XYM]//g;
			$chromosome = "chr" . $chromosome;

			my $chrom_key = $chromosome . "\t" . substr($chr_start, 0, 3);

			$indels{$chrom_key} .= "$chromosome\t$chr_start\t$chr_stop\t$indel_type\t$indel_size\t$read_coverage\t$reads1\t$reads2\n";


#			print "$lineCounter indels loaded...\n" if(!($lineCounter % 100000));
		}
	}
	
	close($input);

	return(%indels);
}

sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

