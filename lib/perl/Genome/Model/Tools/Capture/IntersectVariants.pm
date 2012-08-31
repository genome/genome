
package Genome::Model::Tools::Capture::IntersectVariants;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# IntersectVariantsForAnnotation - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	10/23/2009 by D.K.
#	MODIFIED:	10/23/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my %variants = ();

class Genome::Model::Tools::Capture::IntersectVariants {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_files	=> { is => 'Text', doc => "List of files to intersect, comma-separated", is_optional => 0, is_input => 1 },
		output_file     => { is => 'Text', doc => "Output file to receive formatted lines", is_optional => 0, is_input => 1, is_output => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Intersects variants between two or more files"                 
}

sub help_synopsis {
    return <<EOS
This command intersects variants between two or more files.
EXAMPLE:	gmt capture intersect-variants --variant-files [file1],[file2] --output-file [file]
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
	my $variant_files = $self->variant_files;
	my $output_file = $self->output_file;
	
	my @variant_files = split(/\,/, $variant_files);
	my $num_files = @variant_files;
	
	if($num_files > 1)
	{
		print "$num_files files will be intersected\n";		

		## Load all variant files ##
		
		foreach my $file (@variant_files)
		{
			if($file && -e $file)
			{
				load_variants($file);
			}
		}

		## Open outfile ##
		
		open(OUTFILE, ">$output_file") or warn "Can't open outfile: $!\n";
	
		my $num_variants = my $num_intersected = 0;
	
		foreach my $key (sort byChrPos keys %variants)
		{
			$num_variants++;
			## Print the variant if it was found in all files ##
			
			if($variants{$key} && $variants{$key} >= $num_files)
			{
				$num_intersected++;
				print OUTFILE "$key\n";
			}
		}
		
		close(OUTFILE);
		
		print "$num_variants unique variants\n";
		print "$num_intersected intersected\n";

	}
	else
	{
		warn "Error: Please provide two or more files\n";
	}

	


}



sub load_variants
{
	my $FileName = shift(@_);

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	my @formatted = ();
	my $formatCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		(my $chrom, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
		
		my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
		
		$variants{$key}++;
	}

	close($input);	
}




sub byChrPos
{
    (my $chrom_a, my $pos_a) = split(/\t/, $a);
    (my $chrom_b, my $pos_b) = split(/\t/, $b);

	$chrom_a =~ s/X/23/;
	$chrom_a =~ s/Y/24/;
	$chrom_a =~ s/MT/25/;
	$chrom_a =~ s/M/25/;
	$chrom_a =~ s/[^0-9]//g;

	$chrom_b =~ s/X/23/;
	$chrom_b =~ s/Y/24/;
	$chrom_b =~ s/MT/25/;
	$chrom_b =~ s/M/25/;
	$chrom_b =~ s/[^0-9]//g;

    $chrom_a <=> $chrom_b
    or
    $pos_a <=> $pos_b;
    
#    $chrom_a = 23 if($chrom_a =~ 'X');
#    $chrom_a = 24 if($chrom_a =~ 'Y');
    
}


1;

