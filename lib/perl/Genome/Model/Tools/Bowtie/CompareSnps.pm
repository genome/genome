
package Genome::Model::Tools::Bowtie::CompareSnps;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# CompareSnps.pm - 	Compare SNP read counts at a specific set of positions
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

class Genome::Model::Tools::Bowtie::CompareSnps {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		positions_file	=> { is => 'Text', doc => "File containing ROI chromosome positions (chr1\t939339)"},
		variants_file1	=> { is => 'Text', doc => "File of containing called readcounts" },
		variants_file2	=> { is => 'Text', doc => "File of containing called readcounts" },
		variants_file3	=> { is => 'Text', doc => "File of containing called readcounts", is_optional => 1 },
		output_file	=> { is => 'Text', doc => "Output file to receive limited SNPs", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Get ReadCounts for a set of SNPs across multiple Bowtie alignments (e.g. Normal-Tumor)"                 
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


my %roi_positions = ();


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	## Get required parameters ##
        my $positions_file = $self->positions_file;
	my $variants_file1 = $self->variants_file1;
	my $variants_file2 = $self->variants_file2;
	my $variants_file3 = $self->variants_file3 if(defined($self->variants_file3));
	my $outfile = $self->output_file if(defined($self->output_file));

        if(!$positions_file || !(-e $positions_file))
        {
                die "\n***Positions file not found!\n";
        }

        if(!$variants_file1 || !(-e $variants_file1))
        {
                die "\n***Variants file 1 not found!\n";
        }

        if(!$variants_file2 || !(-e $variants_file2))
        {
                die "\n***Variants file 2 not found!\n";
        }

#        if(!defined($variants_file3) || !(-e $variants_file3))
 #       {
  #              die "\n***Variants file 3 not found!\n";
   #     }

	my %stats = ();
	$stats{'total_snps'} = $stats{'limit_snps'} = 0;
	my $input, my $lineCounter;

	print "Parsing ROI positions...\n";


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
				$stats{'total_snps'}++;
			}
		}
		
		close($input);
		
		print "$lineCounter positions loaded\n";
	}

	print "Parsing SNPs...\n";
	## Parse the combined SNPs ##

	## Open the outfile ##
	if($outfile)
	{
		open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
	}

	my %readcounts1 = parse_readcounts($variants_file1);
	my %readcounts2 = parse_readcounts($variants_file2);
	my %readcounts3 = parse_readcounts($variants_file3) if($variants_file3 && -e $variants_file3);

	## Go through every position ##
	
	foreach my $variant_key (keys %roi_positions)
	{
		(my $chromosome, my $position) = split(/\t/, $variant_key);
		
		## Get read counts 1 ##
		
		my $results1 = my $results2 = my $results3 = "";
		
		if($readcounts1{$variant_key})
		{
			$results1 = $readcounts1{$variant_key};		
		}
		else
		{
			$results1 = "N\tN\t0\t0";
		}

		if($readcounts2{$variant_key})
		{
			$results2 = $readcounts2{$variant_key};		
		}
		else
		{
			$results2 = "N\tN\t0\t0";
		}

		if($readcounts3{$variant_key})
		{
			$results3 = $readcounts3{$variant_key};		
		}
		else
		{
			if($variants_file3 && -e $variants_file3)
			{			
				$results3 = "N\tN\t0\t0";
			}
		}


		print OUTFILE "$chromosome\t$position\t$results1\t$results2\t$results3\n" if($outfile);
		
	}

	
	close(OUTFILE) if($outfile);

	print "$stats{'total_snps'} total SNPs\n";

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# Parse Readcounts
#
################################################################################################

sub parse_readcounts
{
	my $variants_file = shift(@_);
	my %readcounts = ();

	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter == 1)
		{

		}
		if($lineCounter > 1)# && $lineCounter < 1000)
		{
			(my $chrom, my $position, my $allele1, my $allele2, my $coverage, my $reads1, my $reads2) = split(/\t/, $line);
			
			my $variant_key = "$chrom\t$position";
		
			if($roi_positions{$variant_key})
			{
				$readcounts{$variant_key} = "$allele1\t$allele2\t$reads1\t$reads2";
			}
		}
	}

	close($input);

	return(%readcounts);
}


1;

