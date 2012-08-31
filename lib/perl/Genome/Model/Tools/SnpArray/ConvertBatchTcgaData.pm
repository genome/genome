
package Genome::Model::Tools::SnpArray::ConvertBatchTcgaData;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::SnpArray::ConvertBatchTcgaData {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		level_2_folders	=> { is => 'Text', doc => "Relative path(s) to extracted Level 2 data archives (should contain birdseed files)", is_optional => 0, is_input => 1 },
		mage_tab_files	=> { is => 'Text', doc => "Relative path(s) to the .sdrf file in the mage-tab archive", is_optional => 0, is_input => 1 },
		use_bsub	=> { is => 'Text', doc => "IF set to 1, will bsub the actual import step", is_optional => 1, is_input => 1 },
		snp_info_file	=> { is => 'Text', doc => "SNP array information file downloaded from UCSC", is_optional => 0, is_input => 1, default=> "/gscmnt/xp4101/info/medseq/snp_array_data/snpArrayAffy6.txt" },
		output_dir	=> { is => 'Text', doc => "Output directory to hold converted genotype files", is_optional => 0, is_input => 1 },
		report_only	=> { is => 'Text', doc => "IF set to 1, will skip the actual conversion step and report stats", is_optional => 1, is_input => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Performs batch-conversion of SNP array data for TCGA projects"                 
}

sub help_synopsis {
    return <<EOS
This command performs batch-conversion of SNP array data for TCGA projects
EXAMPLE:	gmt snp-array convert-batch-tcga-data --level-2-folders broad.mit.edu_COAD.Genome_Wide_SNP_6.Level_2.28.1002.0 --mage-tab-files broad.mit.edu_COAD.Genome_Wide_SNP_6.mage-tab.1.1002.0/broad.mit.edu_COAD.Genome_Wide_SNP_6.sdrf.txt
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

	my $output_dir = $self->output_dir;
	my $snp_info_file = $self->snp_info_file;


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
#					my $washu_sample_name = get_washu_sample_name($sample_name);
					
#					if($washu_sample_name)
#					{
						$stats{'have_washu_name'}++;

						print "Attempting to convert this sample:\n";
#						print join("\t", $key, $sample_name, $washu_sample_name) . "\n";						
#						my $output_file = $self->output_dir . "/" . $washu_sample_name . ".genotype";						
						print join("\t", $key, $sample_name) . "\n";						
						my $output_file = $self->output_dir . "/" . $sample_name . ".genotype";						

						my $cmd = "gmt snp-array birdseed-to-genotype --birdseed-file $path_to_file --output-file $output_file --snp-info-file $snp_info_file";

						if(!($self->report_only))
						{
							if($self->use_bsub)
							{
								system("bsub -q tcga -oo $output_file.log -R\"select[model != Opteron250 && mem>4000] rusage[mem=4000]\" \"$cmd\"");								
							}
							else
							{
								system($cmd);															
							}

						}
						else
						{
							print "$cmd\n";
						}

						$stats{'num_converted'}++;
#					}


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
#		print $stats{'have_washu_name'} . " have a WashU sample name in the database\n";
		print $stats{'num_converted'} . " genotype files converted\n";
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

