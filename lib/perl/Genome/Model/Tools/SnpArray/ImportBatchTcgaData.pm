
package Genome::Model::Tools::SnpArray::ImportBatchTcgaData;     # rename this when you give the module file a different name <--

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

my %stats = ();
my %sample_arrays = ();
my %array_sample_names = ();

class Genome::Model::Tools::SnpArray::ImportBatchTcgaData {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		level_2_folders	=> { is => 'Text', doc => "Relative path(s) to extracted Level 2 data archives (should contain birdseed files)", is_optional => 0, is_input => 1 },
		mage_tab_files	=> { is => 'Text', doc => "Relative path(s) to the .sdrf file in the mage-tab archive", is_optional => 0, is_input => 1 },
		use_bsub	=> { is => 'Text', doc => "IF set to 1, will bsub the actual import step", is_optional => 1, is_input => 1 },
		import_format	=> { is => 'Text', doc => "Input format of the data. Currently only 'unknown' is allowed", is_optional => 0, is_input => 1, default => "unknown" },
		description	=> { is => 'Text', doc => "Short description of the import data", is_optional => 1, is_input => 1 },
		source_name	=> { is => 'Text', doc => "Source name for the data, e.g. \"Broad Institute\"", is_optional => 1, is_input => 1 },
		reference_sequence_build	=> { is => 'Text', doc => "Reference sequence build [b36=101947881 b37=106942997]", is_optional => 0, is_input => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Performs batch-imports of SNP array data for TCGA projects"                 
}

sub help_synopsis {
    return <<EOS
This command performs batch-imports of SNP array data for TCGA projects
EXAMPLE:	gmt snp-array import-batch-tcga-data --level-2-folders broad.mit.edu_COAD.Genome_Wide_SNP_6.Level_2.28.1002.0 --mage-tab-files broad.mit.edu_COAD.Genome_Wide_SNP_6.mage-tab.1.1002.0/broad.mit.edu_COAD.Genome_Wide_SNP_6.sdrf.txt
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
	my $mage_tab_files = $self->mage_tab_files;
	my $level_2_folders = $self->level_2_folders;

	my @mage_tab_files = split(/\,/, $mage_tab_files);
	my @level_2_folders = split(/\,/, $level_2_folders);

	my $import_format = $self->import_format;
	my $reference_sequence_build = $self->reference_sequence_build;

	## Get list of samples that already have SNP Array data imported ##
	
	my $sample_list = `genome instrument-data list imported --noheaders --show=sample_name --filter=sequencing_platform=\"affymetrix genotype array\"`;
	chomp($sample_list);
	my @sample_lines = split(/\n/, $sample_list);
	my %already_imported = ();
	
	foreach my $sample (@sample_lines)
	{
		$already_imported{$sample} = 1;
	}

	## Load information in all mage-tab files ##
	
	foreach my $mage_tab_file (@mage_tab_files)
	{
		load_mage_tab($mage_tab_file);		
	}


	## Find any Level 2 Birdseed Files ##
	
	foreach my $level_2_folder (@level_2_folders)
	{
		## Get Birdseed data files ##
		
		opendir (DIR, "$level_2_folder") or warn "Unable to open directory: $level_2_folder";

		my @temp = split(/\//, $level_2_folder);
		my $numElements = @temp;
		my $array_name = $temp[$numElements - 1];
		
		foreach my $filename (readdir DIR)
		{
			if($filename =~ 'birdseed\.data\.txt')
			{
				$stats{'num_birdseed_files'}++;
				my $path_to_file = $level_2_folder . "/" . $filename;
				## Process birdseed file ##
				
				my $key = join("\t", $array_name, $filename);
				if($array_sample_names{$key})
				{
					$stats{'have_tcga_name'}++;
					my $sample_name = $array_sample_names{$key};
					my $washu_sample_name = get_washu_sample_name($sample_name);
					
					if($washu_sample_name)
					{
						$stats{'have_washu_name'}++;
						if($already_imported{$washu_sample_name})
						{
							$stats{'num_already_imported'}++;
							print join("\t", $key, $sample_name, $washu_sample_name) . "\talready_imported\n";												
						}
						else
						{
							$stats{'num_imported'}++;
							print "Attempting to import this sample:\n";
							print join("\t", $key, $sample_name, $washu_sample_name) . "\n";						
							## Run the import ##
							
							my $cmd = "";
							
							if($self->use_bsub)
							{
								$cmd = "bsub -q short genome instrument-data import microarray affymetrix-genotype-array --reference-sequence-build $reference_sequence_build --original-data-file=$path_to_file --sample-name=\"$washu_sample_name\" --import-format=\"$import_format\"";
							}
							else
							{
								$cmd = "genome instrument-data import microarray affymetrix-genotype-array --reference-sequence-build $reference_sequence_build --original-data-file=$path_to_file --sample-name=\"$washu_sample_name\" --import-format=\"$import_format\"";
							}

							## Append optional fields ##
							
							if($self->description)
							{
								$cmd .= " --description=\"" . $self->description . "\"";
							}

							if($self->source_name)
							{
								$cmd .= " --import-source-name=\"" . $self->source_name . "\"";
							}

							## Run the import ##

							system($cmd);
						}						
					}


				}
				else
				{
					print "No mage-tab sample name for $key\n";
				}
			}

		}
		
		closedir(DIR);


		print "\n************CURRENT STATS****************\n";
		print $stats{'num_birdseed_files'} . " birdseed genotype files\n";
		print $stats{'have_tcga_name'} . " were matched to a TCGA sample identifier using mage-tab info\n";
		print $stats{'have_washu_name'} . " have a WashU sample name in the database\n";
		print $stats{'num_imported'} . " genotype files skipped because they're already imported\n";
		print $stats{'num_imported'} . " genotype files imported\n";
		print "*****************************************\n";
		
		
	}


}
	

################################################################################################
# Load Genotypes
#
################################################################################################

sub load_mage_tab
{                               # replace with real execution logic.
	my $mage_tab_file = shift(@_);

	my $input = new FileHandle ($mage_tab_file);
	my $lineCounter = 0;
	my $gtCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if($lineCounter > 1)
		{
			my @lineContents = split(/\t/, $line);
			my $sample_name = $lineContents[0];
			my $archive_name = $lineContents[20];
			my $birdseed_filename = $lineContents[27];
			
			my $key = join("\t", $archive_name, $birdseed_filename);
			
			$array_sample_names{$key} = $sample_name;
			$sample_arrays{$sample_name} = $key;
		}

	}
	close($input);



}




################################################################################################
# Load Genotypes
#
################################################################################################

sub get_washu_sample_name
{
	my $sample = shift(@_);
	
	my @sampleContents = split(/\-/, $sample);
	
	my $query_string = join("-", $sampleContents[1], $sampleContents[2], $sampleContents[3], $sampleContents[4]);
	
	my $washu_name = `genome sample list --noheaders --show=name --filter=name~'\%$query_string\%'`;
	chomp($washu_name);
	
	return($washu_name) if($washu_name);
	return;
}


1;

