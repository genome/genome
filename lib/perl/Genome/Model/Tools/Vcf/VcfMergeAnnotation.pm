package Genome::Model::Tools::Vcf::VcfMergeAnnotation;
##########################################################################
# Extract the annotations from a vep and merge them into a Vcf.
# Require a pair of Vcf and vep as inputs.
# Output a merged Vcf (does not overwrite orginal Vcf).
#
#       AUTHOR:         Yang Li (yl5682@truman.edu)
#       CREATED:        05/28/2013
#       MODIFIED:
#
#       NOTES:
#
###########################################################################
use strict;
use warnings;
use Genome;
use File::stat;
use IO::File;
use File::Basename;
use Getopt::Long;
use FileHandle;
use POSIX qw(log10);
use POSIX qw(strftime);
use List::MoreUtils qw(firstidx);
use List::MoreUtils qw(uniq);

class Genome::Model::Tools::Vcf::VcfMergeAnnotation {
    is => 'Command',
    has => [
        vcf_file_name => {
            is => 'Text',
            is_input => 1,
            doc => "filtered VCF file",
            is_optional => 0,
        },

        vep_file_name => {
            is => 'Text',
            is_input => 1,
            doc => "mutations in VEP format",
            is_optional => 0,
        },

        output_file_name => {
            is => 'String',
            is_output => 1,
            doc => "a user defined output file name",
            is_optional => 0,            
        },

        threshold => {
            is => 'Int',
            is_output => 0,
            doc => "rank-by-severity threshold (0 - 13); higher = more severe",
            is_optional => 1,            
        },
        ],
};

sub help_brief {# keep this to just a few words
    "Merge annotations from a VEP to a VCF"
}

sub help_synopsis {
<<'HELP';
    Extract annotations from a VEP and append them to the INFO field of a VCF
HELP
}

sub help_detail {# this is what the user will see with the longer version of help.
<<'HELP';
    Extract 'Gene', 'Consequences', 'Protein_position', 'Amino_acids', and 'Extra' fields from a VEP and append them to the 'INFO' field of the corresponding VCF. Each field is separated from another by a '.'.
HELP
}

################################################################################################
# Execute - the main program logic
################################################################################################
sub execute {
    my $self = shift;
    my $vcfFileName = $self->vcf_file_name;
    my $vepFileName = $self->vep_file_name;
	my $outputFileName = $self->output_file_name;
	my $threshold = $self->threshold;
	my $outputFile;
	
	open (my $vepFile, "<$vepFileName") or die "Cannot open $vepFileName. $!.\n"; # Open Vep file
	open (my $vcfFile, "<$vcfFileName") or die "Cannot open $vcfFileName. $!.\n"; # Open Vcf file
	my ($fileName, $directory, $ext) = fileparse($vcfFileName, qr/\.[^.]*/); # Parse file name, directory, and extension
	open ($outputFile, ">$outputFileName") or die "Cannot create output file. $!\n"; # Open / Create the output file
	
	$threshold = 0 if (!$threshold); # Default $threshold to 0 if not provided
	
	# Call two sub programs
	my ($numOfLines, @vep) = $self->parceVep($vepFile); # Parse information from Vep
	$self->annotateVcf($threshold, $vcfFile, $numOfLines, $outputFile, \@vep); # Merge annotation into Vcf
    
    return 1;
}

# vep Prodcedure
sub parceVep{
	my $self = shift;
	my $vepFile = shift;
    my @vep;
	
	my $numOfLines = 0; # Total number of lines in the Vep file
	# Loop through the vep file
	while(my $line = <$vepFile>){
		chomp $line;
		next if((substr $line, 0, 1) eq "#"); # Skip the comments
		# Parse out the fields from a line
		my ($Uploaded_variation,
			$Location,
			$Allele,
			$Gene,
			$Feature,
			$Feature_type,
			$Consequence,
			$cDNA_position,
			$CDS_position,
			$Protein_position,
			$Amino_acids,
			$Codons,
			$Existing_variation,
			$Extra) = split("\t", $line);
		my ($CHROM, $POS) = split(/:/, $Location); # Break the fields used for identifier

		# Create a record that will be stored in an array called @vep
    	my $rec = {};
    	{
			$rec->{CHROM} = $CHROM;
			$rec->{POS} = $POS;
			$rec->{GENE} = $Gene;
			$rec->{CONSEQ} = $Consequence;
			$rec->{PROTEINPOS} = $Protein_position;
			$rec->{AMINOACIDS} = $Amino_acids;
			$rec->{EXTRA} = $Extra;
		}
		push @vep, $rec; # Append a record into the vep array
		$numOfLines++; # Increment the line counter
		
		## Debugging codes
		#print "Chrom: $rec->{CHROM}, Pos: $rec->{POS}\n";
    	#print "Gene: $rec->{GENE}, Conseq: $rec->{CONSEQ}\n";
    	#print "Protein position: $rec->{PROTEINPOS}, Extra: $rec->{EXTRA}\n\n";
	}
	
	# Close the file
	close ($vepFile);
	return ($numOfLines, @vep);
}

# Vcf Prodcedure
sub annotateVcf{
	my $self = shift;
	my $threshold = shift;
	my $vcfFile = shift;
	my $numOfLines = shift;
	my $outputFile = shift;
	my $vepRef=shift;
	my @vep=@{$vepRef}; # Dereference the array reference
	
	## Pre-define a ranking system for VEP annotation, where higher = more severe ##
	my %vep_class_rank = ();
	$vep_class_rank{'-'} =                 0;
	$vep_class_rank{'NMD_TRANSCRIPT'} =         0;
	$vep_class_rank{'INTERGENIC'} =         0;
	$vep_class_rank{'UPSTREAM'} =             1;
	$vep_class_rank{'DOWNSTREAM'} =         2;
	$vep_class_rank{'INTRONIC'} =             3;
	$vep_class_rank{'COMPLEX_INDEL'} =             3;
	$vep_class_rank{'5PRIME_UTR'} =         4;
	$vep_class_rank{'3PRIME_UTR'} =         5;
	$vep_class_rank{'WITHIN_NON_CODING_GENE'} =     6;
	$vep_class_rank{'WITHIN_MATURE_miRNA'} =     7;
	$vep_class_rank{'PARTIAL_CODON'} =         7;
	$vep_class_rank{'CODING_UNKNOWN'} =             7;
	$vep_class_rank{'SYNONYMOUS_CODING'} =         8;
	$vep_class_rank{'STOP_LOST'} =             9;
	$vep_class_rank{'SPLICE_SITE'} =         10;
	$vep_class_rank{'ESSENTIAL_SPLICE_SITE'} =     11;
	$vep_class_rank{'NON_SYNONYMOUS_CODING'} =     12;
	$vep_class_rank{'STOP_GAINED'} =         13;
	$vep_class_rank{'FRAMESHIFT_CODING'} =         13; 
	
	my(	$CHROM,
		$POS,
		$ID,
		$REF,
		$ALT,
		$QUAL,
		$FILTER,
		$INFO,
		$FORMAT,
		$H_MI_WX_0009_1121413,
		$H_MI_WX_0009_2_1120675P,
		$H_MI_WX_0009_3_1120676P);
	
	# Last line counter (reduce repetitive iterations)
	my $counter = 0; 
	
	# Loop through the Vcf file
	while (my $line = <$vcfFile>){
		#print "VCF: $line\n";
		chomp $line;
		if((substr $line, 0, 1) eq "#"){ # Skip the comments
			print $outputFile $line."\n"; # Write the comments into the output file
			next;
		}
		
		# Parse out the fields from a line
		my ($CHROM,
			$POS,
			$ID,
			$REF,
			$ALT,
			$QUAL,
			$FILTER,
			$INFO,
			$FORMAT,
			$H_MI_WX_0009_1121413,
			$H_MI_WX_0009_2_1120675P,
			$H_MI_WX_0009_3_1120676P) = split("\t", $line);
		
		# Loop through the vep file and add annotations if found
		for (my $i=$counter; $i<$numOfLines; $i++){
			if ($CHROM eq $vep[$i]{CHROM}){ # Chromosome match
				if ($POS == $vep[$i]{POS}){ # Position match					
					$INFO .= ".Gene=$vep[$i]{GENE}";
					if ($vep_class_rank{$vep[$i]{CONSEQ}} >= $threshold){
						$INFO .= ".Consequences=$vep[$i]{CONSEQ}";
					}
					$INFO .=".Protein_position=$vep[$i]{PROTEINPOS}.Amino_acids=$vep[$i]{AMINOACIDS}.$vep[$i]{EXTRA}"; #Append annotation to the INFO field
					$INFO =~s/^\.+//; # Remove leading dots
					$counter = $i++; # Increment counter
				}
			}
		}

		# Prepare a line for the output file
		my $result = "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$H_MI_WX_0009_1121413\t$H_MI_WX_0009_2_1120675P\t$H_MI_WX_0009_3_1120676P\n";

		# Write the line into the output file
		print $outputFile $result;
	}
	
	# Close the files
	close ($vcfFile);
	close ($outputFile);
}
