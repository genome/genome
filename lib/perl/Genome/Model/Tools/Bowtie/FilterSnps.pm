
package Genome::Model::Tools::Bowtie::FilterSnps;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FilterSnps.pm - 	Limit SNPs by various factors, like ROI positions.
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	03/20/2009 by D.K.
#	MODIFIED:	03/20/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Bowtie::FilterSnps {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File containing SNPs" },
		output_file	=> { is => 'Text', doc => "Output file to receive limited SNPs" },
		not_file	=> { is => 'Text', doc => "Output file to receive excluded SNPs", is_optional => 1 },
		verbose	=> { is => 'Text', doc => "Print interim progress updates", is_optional => 1 },
		positions_file	=> { is => 'Text', doc => "File containing ROI chromosome positions (chr1\t939339)", is_optional => 1 },
		min_coverage	=> { is => 'Text', doc => "Minimum read coverage to include SNP", is_optional => 1 },
		min_reads1	=> { is => 'Text', doc => "Minimum ref reads to include SNP", is_optional => 1 },
		min_reads2	=> { is => 'Text', doc => "Minimum var reads to include SNP", is_optional => 1 },
		min_strands2	=> { is => 'Text', doc => "Minimum var strands to include SNP", is_optional => 1 },
		min_avg_qual	=> { is => 'Text', doc => "Minimum average base quality to include SNP", is_optional => 1 },
		min_var_freq	=> { is => 'Text', doc => "Minimum variant frequency to include SNP", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Filter Bowtie SNP readcounts files by coverage, reads2, frequency, etc."                 
}

sub help_synopsis {
    return <<EOS
This command parses alignments from Bowtie.
EXAMPLE:	gmt bowtie parse-alignments --alignments-file bowtie.txt
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
	my $variants_file = $self->variants_file;
        my $positions_file = $self->positions_file if(defined($self->positions_file));
	my $notfile = $self->not_file if(defined($self->not_file));
	my $outfile = $self->output_file;
	my $verbose = 1 if(defined($self->verbose));
	my $min_coverage = $self->min_coverage if(defined($self->min_coverage));
	my $min_reads1 = $self->min_reads1 if(defined($self->min_reads1));
	my $min_reads2 = $self->min_reads2 if(defined($self->min_reads2));
	my $min_strands2 = $self->min_strands2 if(defined($self->min_strands2));
	my $min_avg_qual = $self->min_avg_qual if(defined($self->min_avg_qual));
	my $min_var_freq = $self->min_var_freq if(defined($self->min_var_freq));

	my %stats = ();
	$stats{'total_snps'} = $stats{'limit_snps'} = 0;
	my $input, my $lineCounter;

	print "Parsing ROI positions...\n";

	my %roi_positions = ();

	if($positions_file)
	{
		## Parse the alignment blocks ##
	
		$input = new FileHandle ($positions_file);
		$lineCounter = 0;
		
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;
			
			if($lineCounter >= 1)# && $lineCounter < 50000)
			{
				(my $chrom, my $position) = split(/\t/, $line);
				$roi_positions{$chrom . "\t" . $position} = 1;
			}
		}
		
		close($input);
		
		print "$lineCounter positions loaded\n";
	}

	print "Parsing SNPs...\n";
	## Parse the combined SNPs ##

	## Open the outfile ##
	open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
	
	if($notfile)
	{
		open(NOTFILE, ">$notfile") or die "Can't open outfile: $!\n";
	}

	$input = new FileHandle ($variants_file);
	$lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter == 1)
		{
			print OUTFILE "$line\n";
			print NOTFILE "$line\n" if($notfile);
		}
		if($lineCounter > 1)# && $lineCounter < 1000)
		{
			$stats{'total_snps'}++;
#			(my $chrom, my $position, my $allele1, my $allele2, my $context, my $num_reads, my $reads) = split(/\t/, $line);			
			(my $chrom, my $position, my $allele1, my $allele2, my $coverage, my $reads1, my $reads2, my $readsN, my $avg_ref_qual, my $avg_var_qual, my $num_ref_strands, my $num_var_strands) = split(/\t/, $line);

			my $variant_freq = $reads2 / $coverage if($coverage);

			my $included_flag = 0;

			if(!$positions_file || $roi_positions{$chrom . "\t" . $position})
			{
				if(!$min_coverage || $coverage >= $min_coverage)
				{
					if(!$min_reads2 || $reads2 >= $min_reads2)
					{
						if(!$min_var_freq || $variant_freq >= $min_var_freq)
						{
							if(!$min_avg_qual || $avg_var_qual >= $min_avg_qual)
							{
                                                                if(!$min_strands2 || $num_var_strands >= $min_strands2)
                                                                {
                                                                        $stats{'limit_snps'}++;
                                                                        print OUTFILE "$line\n";
                                                                        $included_flag = 1;
                                                                }
							}
						}
					}
				}
			}
			
			if(!$included_flag)
			{
				print NOTFILE "$line\n" if($notfile);
			}
		}
	}

	close($input);
	
	close(OUTFILE);

	## Get percentage ##
	
	$stats{'limit_pct'} = sprintf("%.2f", ($stats{'limit_snps'} / $stats{'total_snps'} * 100)) . '%' if($stats{'total_snps'});

	print "$stats{'total_snps'} total SNPs\n";
	print "$stats{'limit_snps'} SNPs ($stats{'limit_pct'}) remained after filter\n";
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

