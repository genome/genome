
package Genome::Model::Tools::Analysis::454::MergeMafs;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MergeMafs - Align reads with SSAHA2 or other aligner
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

my %file1_results = my %file2_results = ();

class Genome::Model::Tools::Analysis::454::MergeMafs {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file1	=> { is => 'Text', doc => "Original MAF file" },
		maf_file2	=> { is => 'Text', doc => "Modified or partial maf file" },
		output_file		=> { is => 'Text', doc => "Output file for updated MAF" },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merges changes from maf-file2 into maf-file1 "                 
}

sub help_synopsis {
    return <<EOS
This command updates variants in a MAF file with 454 validation status
EXAMPLE:	gt analysis 454 merge-mafs --maf-file1 original.maf --maf-file2 modified.maf --output-file merged.maf
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
	my $maf_file1 = $self->maf_file1;
	my $maf_file2 = $self->maf_file2;
	my $output_file = $self->output_file;

	if(!(-e $maf_file1))
	{
		die "Error: MAF file 1 not found!\n";
	}

	if(!(-e $maf_file2))
	{
		die "Error: MAF file 2 not found!\n";
	}

	%file1_results = load_maf_file($maf_file1);
	my %file2_novel = get_novel_entries($maf_file2);
	
	
	%file2_results = load_maf_file($maf_file2);

	my %stats = ();


	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

	my @output_lines = ();
	my $num_lines = 0;

	## Parse the MAF file ##
	
	my $input = new FileHandle ($maf_file1);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter == 1 && $line =~ "Chrom")
		{
			## Print the MAF header ##
			print OUTFILE "$line\n";
		}
		elsif($lineCounter == 1)
		{
			warn "Warning: No Header in MAF file!\n";
		}
		else
		{
			my @lineContents = split(/\t/, $line);			
			my $chromosome = $lineContents[4];
			my $chr_start = $lineContents[5];
			my $chr_stop = $lineContents[6];
			my $trv_type = $lineContents[8];
			my $var_type = $lineContents[9];
			my $ref_base = $lineContents[10];
			my $var_base = $lineContents[11];
			my $tumor_sample = $lineContents[15];
			my $normal_sample = $lineContents[16];
			my $current_val_status = $lineContents[24];

			$var_base = $lineContents[12] if($lineContents[12] ne $ref_base);

			## Build the key ##
			
			my $key = "$chromosome\t$chr_start\t$chr_stop\t$trv_type\t$var_type\t$tumor_sample\t$normal_sample";

			my $newline = $line;

			if($file2_results{$key})
			{
				$stats{'num_match'}++;
				
				my @file2contents = split(/\t/, $file2_results{$key});
				my $new_val_status = $file2contents[24];
				
				if($new_val_status ne $current_val_status)
				{
					$stats{'num_updated'}++;
					
					$newline = "";
					
					$lineContents[19] = $file2contents[19];
					$lineContents[20] = $file2contents[20];
					$lineContents[21] = $file2contents[21];
					$lineContents[22] = $file2contents[22];
					$lineContents[24] = $file2contents[24];
					$lineContents[28] = $file2contents[28];
					
					my $numContents = @lineContents;
					
					for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
					{
						$newline .= "\t" if($newline);
						$newline .= $lineContents[$colCounter];
					}
				}
			}
			
			## Append new line to output ##
			
			$output_lines[$num_lines] = $newline;
			$num_lines++;
	
		}

	}

	close($input);	
	
	
	## Add novel events ##
	

	foreach my $key (keys %file2_novel)
	{
		$output_lines[$num_lines] = $file2_novel{$key};
		$num_lines++;
		$stats{'novel_in_file2'}++;
	}

	print $stats{'novel_in_file2'} . " entries in file 2 are new\n";
	print $stats{'num_match'} . " records in file 2 match those in file 1\n";
	print $stats{'num_updated'} . " records will be updated\n";	
	
	my @output_sorted = sort byChrPos @output_lines;
	
	foreach my $line (@output_sorted)
	{
		print OUTFILE "$line\n";
	}

	sub byChrPos {
		my @temp = split(/\t/, $a);
		my $chrom_a = $temp[4];
		my $pos_a = $temp[5];

		@temp = ();
		@temp = split(/\t/, $b);
		my $chrom_b = $temp[4];
		my $pos_b = $temp[5];	

		$chrom_a =~ s/X/23/;
		$chrom_a =~ s/Y/24/;
		$chrom_a =~ s/MT/25/;

		$chrom_b =~ s/X/23/;
		$chrom_b =~ s/Y/24/;
		$chrom_b =~ s/MT/25/;

		$chrom_a <=> $chrom_b
		or
		$pos_a <=> $pos_b;
	}
	
	
	close(OUTFILE);



	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



#############################################################
# parse_file - takes input file and parses it
#
#############################################################

sub load_maf_file
{
	my $FileName = shift(@_);

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	my %results = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my @lineContents = split(/\t/, $line);			
		my $chromosome = $lineContents[4];
		my $chr_start = $lineContents[5];
		my $chr_stop = $lineContents[6];
		my $trv_type = $lineContents[8];
		my $var_type = $lineContents[9];
		my $tumor_sample = $lineContents[15];
		my $normal_sample = $lineContents[16];
		my $validation_status = $lineContents[24];
		my $validation_method = $lineContents[28];
		my $validation_alleles = $lineContents[19] . "\t" . $lineContents[20] . "\t" . $lineContents[21] . "\t" . $lineContents[22];
		
		my $key = "$chromosome\t$chr_start\t$chr_stop\t$trv_type\t$var_type\t$tumor_sample\t$normal_sample";
#		my $result = "$validation_status\t$validation_method\t$validation_alleles";

		$results{$key} = $line;


	}

	close($input);
	
	return(%results);	
}


#############################################################
# parse_file - takes input file and parses it
#
#############################################################

sub get_novel_entries
{
	my $FileName = shift(@_);

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	my %results = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my @lineContents = split(/\t/, $line);			
		my $chromosome = $lineContents[4];
		my $chr_start = $lineContents[5];
		my $chr_stop = $lineContents[6];
		my $trv_type = $lineContents[8];
		my $var_type = $lineContents[9];
		my $tumor_sample = $lineContents[15];
		my $normal_sample = $lineContents[16];
		my $validation_status = $lineContents[24];
		my $validation_method = $lineContents[28];
		my $validation_alleles = $lineContents[19] . "\t" . $lineContents[20] . "\t" . $lineContents[21] . "\t" . $lineContents[22];
		
		my $key = "$chromosome\t$chr_start\t$chr_stop\t$trv_type\t$var_type\t$tumor_sample\t$normal_sample";
#		my $result = "$validation_status\t$validation_method\t$validation_alleles";

		if(!$file1_results{$key})
		{
			$results{$key} = $line;			
		}

	}

	close($input);
	
	return(%results);	
}




1;

