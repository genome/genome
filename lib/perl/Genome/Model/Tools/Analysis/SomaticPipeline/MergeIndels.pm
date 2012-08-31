
package Genome::Model::Tools::Analysis::SomaticPipeline::MergeIndels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MergeIndels - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
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

my %known_dbsnps = ();

class Genome::Model::Tools::Analysis::SomaticPipeline::MergeIndels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		glf_indels	=> { is => 'Text', doc => "TSV file of normalID, tumorID, path to SNVs and annotations", is_optional => 0 },
		sammy_indels	=> { is => 'Text', doc => "TSV file of normalID, tumorID, and path to SNVs", is_optional => 0 },
		regions_file	=> { is => 'Text', doc => "Tab-delimited file of target regions", is_optional => 0 },
		dbsnp_file	=> { is => 'Text', doc => "Tab-delimited file containing known dbSNPs", is_optional => 0 },
		output_dir     => { is => 'Text', doc => "Output dir for intermediate files", is_optional => 0 },
		output_file     => { is => 'Text', doc => "Output file to receive merged data", is_optional => 0 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merges variant predictions from multiple algorithms"                 
}

sub help_synopsis {
    return <<EOS
This command merges somatic variant calls between glfSomatic and Varscan/Sammy
EXAMPLE:	gmt analysis sammy merge-variants --glf-snvs [file] --sammy-snvs [file] --regions-file [file] --dbsnp-file [file] --output-file [file]
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
	my $glf_indels = $self->glf_indels;
	my $sammy_indels = $self->sammy_indels;
	my $regions_file = $self->regions_file;
	my $dbsnp_file = $self->dbsnp_file;
	my $output_file = $self->output_file;
	my $output_dir = $self->output_dir;
	
	print "Loading dbSNPs...\n";
	%known_dbsnps = load_dbsnp_rs($dbsnp_file);
	
	
	## Build directory of samples ##
	
	my %normal_for_tumor = ();
	my %glf_files = my %sammy_files = ();

	print "Parsing glfSomatic variant file list...\n";

	my $input = new FileHandle ($glf_indels);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $normal_sample, my $tumor_sample, my $tier1_file, my $annotation_file) = split(/\t/, $line);
		print "$tumor_sample\n";
		## Merge with annotation ##
		
		$normal_for_tumor{$tumor_sample} = $normal_sample;
		
		## Filter the indel file ##
		
		my $filtered_file = $output_dir . "/" . $tumor_sample . ".glfsomatic.indels.filter";		
		system("gmt analysis somatic-pipeline filter-glf-indels --variants-file $tier1_file --output-file $filtered_file --min-var-freq 0.30 --min-reads2 4 --min-coverage 10"); #--min-var-freq 0.25
	
		## Get number of indels passing filter ##
		
#		my $num_pass_filter = `cat $filtered_file | wc -l`;
#		chomp($num_pass_filter);
#		my $min_freq = 0.25;
#		while($num_pass_filter > 4)
#		{
#			$min_freq += 0.05;
#			system("gmt analysis somatic-pipeline filter-glf-indels --variants-file $tier1_file --output-file $filtered_file --min-var-freq $min_freq --min-reads2 4 --min-coverage 20"); #--min-var-freq 0.25
#			$num_pass_filter = `cat $filtered_file | wc -l`;
#			chomp($num_pass_filter);			
#		}
	
		my $merged_file = $output_dir . "/" . $tumor_sample . ".glfsomatic.indels.filter.annotation.merged";
		
#		if(!(-e $merged_file))
#		{
			system("gmt analysis somatic-pipeline merge-indels-with-annotation --variants-file $filtered_file --annotation-file $annotation_file --output-file $merged_file");
#		}

		$glf_files{$tumor_sample} = $merged_file;
	}

	close($input);



	print "Parsing Sammy variant file list...\n";

	$input = new FileHandle ($sammy_indels);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $normal_sample, my $tumor_sample, my $tier1_file) = split(/\t/, $line);
		print "$tumor_sample\n";
		
		$normal_for_tumor{$tumor_sample} = $normal_sample;
		$sammy_files{$tumor_sample} = $tier1_file;

	}

	close($input);

	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	print OUTFILE "Normal_Sample\tTumor_Sample\tChrom\tChr_Start\tChr_Stop\tRef\tConsensus\tVar_Type\tGene\tTranscript\tTrv_Type\tCodon\tAA_Change\tDomain\tSoftware\tdbSNP\n";	
	
	## Go through all patients we encountered ##
	
	foreach my $tumor_sample (sort keys %normal_for_tumor)
	{
		my $normal_sample = $normal_for_tumor{$tumor_sample};
		
		if($glf_files{$tumor_sample} && $sammy_files{$tumor_sample})
		{
			print "$normal_sample\t$tumor_sample\n";
			
			my $compiled_variants = compile_variants($normal_sample, $tumor_sample, $sammy_files{$tumor_sample}, $glf_files{$tumor_sample});
			
			print OUTFILE $compiled_variants;
#			exit(0);
		}
		else
		{
			print "No glf file for $tumor_sample\n" if(!$glf_files{$tumor_sample});
			print "No Sammy file for $tumor_sample\n" if(!$sammy_files{$tumor_sample});
		}
	}
	
	
	
}







#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub compile_variants
{
	(my $normal_sample, my $tumor_sample, my $sammy_file, my $glf_file) = @_;
	
	my %tier1_variants = my %sammy_variants = my %glf_variants = ();
	my %consensus_calls = ();

	my $compiled_variants = "";	

	my $input;
	my $lineCounter;
	
	if($sammy_file && -e $sammy_file)
	{
		## Parse Sammy ##
	
		$input = new FileHandle ($sammy_file);
		$lineCounter = 0;
		
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;
	
			(my $chromosome, my $position) = split(/\t/, $line);
			my $key = "$chromosome\t$position";
			$tier1_variants{$key} = $line;
			$sammy_variants{$key} = 1;		
	
			my @lineContents = split(/\t/, $line);
			
			my $newline = "";
			for(my $counter = 0; $counter <= 12; $counter++)
			{
				$newline .= "\t" if($newline);
				$lineContents[$counter] = "" if(!$lineContents[$counter]);
				$newline .= "$lineContents[$counter]";
			}
			$tier1_variants{$key} = $newline;
		}
		
		close($input);
		print "Loaded $lineCounter Sammy variants\n";
	}

	## Parse glf ##

	$input = new FileHandle ($glf_file);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chromosome, my $position, my $chr_stop, my $ref_base, my $consensus) = split(/\t/, $line);
		my $key = "$chromosome\t$position";
		my @lineContents = split(/\t/, $line);

		## Save the consensus ##
		
		$consensus_calls{$key} = $consensus;

		my $genotype = consensus_to_gt($consensus);
		
		if(substr($genotype, 0, 1) eq $ref_base)
		{
			$consensus = substr($genotype, 1, 1);	
		}
		elsif(substr($genotype, 1, 1) eq $ref_base)
		{
			$consensus = substr($genotype, 0, 1);
		}

		## Use variant allele in place of consensus ##
		$lineContents[4] = $consensus;

		my $newline = "";
		for(my $counter = 0; $counter <= 12; $counter++)
		{
			$newline .= "\t" if($newline);
			$lineContents[$counter] = "" if(!$lineContents[$counter]);
			$newline .= "$lineContents[$counter]";
		}

		$tier1_variants{$key} = $newline;

		if($lineContents[11] || $lineContents[12] || $lineContents[13])
		{
			$glf_variants{$key} = "$lineContents[11]\t$lineContents[12]\t$lineContents[13]"
		}
		else
		{		
			$glf_variants{$key} = 1;
		}
	}
	
	close($input);
	
	print "Loaded $lineCounter glfSomatic variants\n";
	

	## Go through all variants ##
	
	foreach my $key (keys %tier1_variants)
	{
		## Determine detection Method ##
		my $method = "";
		if($sammy_variants{$key} && $glf_variants{$key})
		{
			$method = "glfSomatic,Varscan";
		}
		elsif($sammy_variants{$key})
		{
			$method = "Varscan";
		}
		elsif($glf_variants{$key})
		{
			$method = "glfSomatic";
		}

		## Get dbsnp ##
		
		my $dbSNP = "";
		$dbSNP = $known_dbsnps{$key} if($known_dbsnps{$key});

		$compiled_variants .= "$normal_sample\t$tumor_sample\t" . $tier1_variants{$key} . "\t" . $method . "\t" . $dbSNP;
		$compiled_variants .= "\t" . $glf_variants{$key} if($glf_variants{$key} && length($glf_variants{$key}) > 1);
		$compiled_variants .= "\t" . $consensus_calls{$key} if($consensus_calls{$key});
		$compiled_variants .= "\n";
	}

	return($compiled_variants);

}


#############################################################
# consensus_to_gt - get genotype equivalent
#
#############################################################

sub consensus_to_gt
{
	my $consensus = shift(@_);
	
	return("AA") if($consensus eq "A");
	return("CC") if($consensus eq "C");
	return("GG") if($consensus eq "G");
	return("TT") if($consensus eq "T");

	return("AC") if($consensus eq "M");
	return("AG") if($consensus eq "R");
	return("AT") if($consensus eq "W");
	return("CG") if($consensus eq "S");
	return("CT") if($consensus eq "Y");
	return("GT") if($consensus eq "K");
	
	return("NN");
}






#############################################################
# load_dbsnp_rs - Load dbSNP RS numbers
#
#############################################################

sub load_dbsnp_rs
{
	my %rs = ();
	
	my $FileName = shift(@_);
	
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chrom, my $position, my $rs_number) = split(/\t/, $line);
		$chrom =~ s/[^0-9XYMT]//g;
		
		my $key = "$chrom\t$position";
		$rs{$key} = $rs_number;
	}

	close($input);

	return(%rs);	
}


#####################################################################################
# call_sammy - Call Sammy3
#
#####################################################################################

sub call_sammy
{
	my $classpath = "/gscuser/dkoboldt/Software/Sammy3";
	my $cmd = "java -Xms3000m -Xmx3000m -classpath $classpath Sammy ";
	return($cmd);
}



sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}

1;

