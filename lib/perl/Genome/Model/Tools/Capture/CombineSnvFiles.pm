
package Genome::Model::Tools::Capture::CombineSnvFiles;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# CombineSnvFilesForAnnotation - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
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
use Genome::Model::Tools::Capture::Helpers 'iupac_to_base';

class Genome::Model::Tools::Capture::CombineSnvFiles {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file1	=> { is => 'Text', doc => "Variants in annotation format", is_optional => 0, is_input => 1 },
		variant_file2	=> { is => 'Text', doc => "Variants in annotation format", is_optional => 0, is_input => 1 },
		output_file     => { is => 'Text', doc => "Output file to receive annotation-formatted lines", is_optional => 0, is_input => 1, is_output => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Combines variants in two files into a single file"                 
}

sub help_synopsis {
    return <<EOS
This command combines variants in two files into a single file, making columns 1-5 unique, and preserving additional columns (if shared, from file 1)
EXAMPLE:	gmt capture combine_snv_files --variant-file [file] --output-file [file]
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

	my %stats = ();

	## Get required parameters ##
	my $variant_file1 = $self->variant_file1;
	my $variant_file2 = $self->variant_file2;
	my $output_file = $self->output_file;
	
	warn "Loading variants in file 1...\n";
	my %variants1 = load_variants($variant_file1);

	print "Loading variants in file 2...\n";
	my %variants2 = load_variants($variant_file2);
	
	my %all_variants = ();
	
	foreach my $key (keys %variants1)
	{
		$all_variants{$key}++;
	}

	foreach my $key (keys %variants2)
	{
		$all_variants{$key}++;
	}
	
	
	## Open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

	foreach my $key (sort byChrPos keys %all_variants)
	{
		$stats{'total variants'}++;
		
		if($variants1{$key} && $variants2{$key})
		{
			my @line1 = split(/\t/, $variants1{$key});
			my @line2 = split(/\t/, $variants2{$key});
			
			my $numContents = @line1;
			my $rest_of_line = "";
			for (my $colCounter = 2; $colCounter < $numContents; $colCounter++)
			{
				$rest_of_line .= "\t" if($rest_of_line);
				$rest_of_line .= $line1[$colCounter];
			}
			
			my $ref1 = $line1[0];
			my $ref2 = $line2[0];
			
			my $var1 = $line1[1];
			my $var2 = $line2[1];
			
			my $ref = $ref1;
			$ref .= "/$ref2" if(!($ref =~ $ref2));

			my $var = $var1;
			$var .= "/$var2" if(!($var =~ $var2));
			
			print OUTFILE join("\t", $key, $ref, $var, $rest_of_line) . "\n";
			$stats{'variants shared'}++;
		}
		elsif($variants1{$key})
		{
			print OUTFILE join("\t", $key, $variants1{$key}) . "\n";
			$stats{'variants unique to file 1'}++;
		}
		elsif($variants2{$key})
		{
			print OUTFILE join("\t", $key, $variants2{$key}) . "\n";
			$stats{'variants unique to file 2'}++;
		}
	}
	
	close(OUTFILE);

	foreach my $key (sort keys %stats)
	{
		print "$stats{$key} $key\n";
	}
	
	return(1);
}



################################################################################################
# Load variants - parse a variant file 
#
################################################################################################

sub load_variants
{                               # replace with real execution logic.
	my $FileName = shift(@_);	

	my %variants = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	my @formatted = ();
	my $formatCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);

		if(!(lc($lineContents[0]) =~ "chrom" || lc($lineContents[0]) =~ "ref_name"))
		{
			my $chrom = $lineContents[0];
			$chrom = fix_chrom($chrom);
			my $chr_start = $lineContents[1];
			my $chr_stop = $lineContents[2];
			my $allele1= $lineContents[3];
			my $allele2 = $lineContents[4];
			my $numContents = @lineContents;
			
			my $rest_of_line = "";
			for(my $colCounter = 5; $colCounter < $numContents; $colCounter++)
			{
				$rest_of_line .= "\t" if($rest_of_line);
				$rest_of_line .= $lineContents[$colCounter];
			}

			if($chrom && $chr_start && $chr_stop)
			{
				my $key = join("\t", $chrom, $chr_start, $chr_stop);
				$variants{$key} = "$allele1\t$allele2";
				$variants{$key} .= "\t$rest_of_line" if($rest_of_line);
			}
		}
	}

	close($input);
	
	return(%variants);
	
}

#############################################################
# ParseBlocks - takes input file and parses it
#
#############################################################

sub fix_chrom
{
	my $chrom = shift(@_);
	$chrom =~ s/chr// if(substr($chrom, 0, 3) eq "chr");
	$chrom =~ s/[^0-9XYMNT\_random]//g;	

	return($chrom);
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

