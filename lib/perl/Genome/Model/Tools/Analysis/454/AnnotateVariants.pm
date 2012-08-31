
package Genome::Model::Tools::Analysis::454::AnnotateVariants;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# AnnotateVariants - Load 454 reads from a sample-SFF tab-delimited file
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

class Genome::Model::Tools::Analysis::454::AnnotateVariants {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "Variants file in Varscan format" },
		output_file	=> { is => 'Text', doc => "Output file" },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Run annotation pipeline on Varscan variants"                 
}

sub help_synopsis {
    return <<EOS
This command reformats Varscan variants and runs annotation on them
	
EXAMPLE: gmt analysis 454 annotate-variants --variants-file snps.tsv --output-file snps.annotated.tsv
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
	my $output_file = $self->output_file;

	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	my $snpCounter = 0;
	
	my @variants = ();
	my %variant_info = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my @lineContents = split(/\t/, $line);			

		if($lineContents[3])
		{
			my $chrom = $lineContents[0];
			my $position = $lineContents[1];
			my $allele1 = $lineContents[2];
			my $allele2 = $lineContents[3];
			
			$chrom =~ s/[^0-9XYM]//g;
			
			if($chrom && $position)
			{
				$variant_info{$chrom . "\t" . $position} = $line;
				$variants[$snpCounter] = "$chrom\t$position\t$position\t$allele1\t$allele2";
				$snpCounter++;
			}
		}

	}
	
	close($input);


	## Open the sorted outfile ##

	open(OUTFILE, ">$variants_file.formatted");
	## Sort variants ##
	
	@variants = sort byChrPos @variants;

	foreach my $variant (@variants)
	{
		print OUTFILE "$variant\n";
	}

	close(OUTFILE);

	## Run the annotation ##
	print "Running annotation...\n";
    my $var_obj = Genome::Model::Tools::Annotate::TranscriptVariants->create(
        variant_file => "$variants_file.formatted",
        output_file => "$variants_file.formatted.annotated",
    );
    $var_obj->execute;

	open(OUTFILE, ">$output_file") or die "Can't open outfile $!\n";


	## Parse the annotation file ##

	$input = new FileHandle ("$variants_file.formatted.annotated");
	$lineCounter = 0;
	my $aCounter = 0;	
	my %annotations = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
		
		if($lineCounter == 1)
		{
			print "$line\n";
		}
		else
		{
			(my $chromosome, my $position) = split(/\t/, $line);
			
			## Get the variant info ##
			
			if($variant_info{$chromosome . "\t" . $position})
			{
				my @variant = split(/\t/, $variant_info{$chromosome . "\t" . $position});
				my $numContents = @variant;
				
				## Print annotation ##
				print OUTFILE "$line";


				## Print variant info ##
				
				for(my $colCounter = 4; $colCounter < $numContents; $colCounter++)
				{
					print OUTFILE "\t$variant[$colCounter]";
				}
				
				print OUTFILE "\n";
			}
	
			$annotations{$chromosome . "\t" . $position} = "" if(!$annotations{$chromosome . "\t" . $position});
			$annotations{$chromosome . "\t" . $position} .= "\n" if($annotations{$chromosome . "\t" . $position});
			$annotations{$chromosome . "\t" . $position} .= $line;
		}
	}

	close($input);




	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


sub byChrPos
{
    (my $chrom_a, my $pos_a) = split(/\t/, $a);
    (my $chrom_b, my $pos_b) = split(/\t/, $b);

	$chrom_a =~ s/X/23/;
	$chrom_a =~ s/Y/24/;
	$chrom_a =~ s/M/25/;

	$chrom_b =~ s/X/23/;
	$chrom_b =~ s/Y/24/;
	$chrom_b =~ s/M/25/;

    $chrom_a <=> $chrom_b
    or
    $pos_a <=> $pos_b;
    
#    $chrom_a = 23 if($chrom_a =~ 'X');
#    $chrom_a = 24 if($chrom_a =~ 'Y');
    
}


1;

