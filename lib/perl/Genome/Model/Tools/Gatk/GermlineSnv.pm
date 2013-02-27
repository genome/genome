
package Genome::Model::Tools::Gatk::GermlineSnv;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# GermlineIndel - Call the GATK germline indel detection pipeline
#					
#	AUTHOR:		Will Schierding
#
#	CREATED:	03-Mar-2011 by W.S.
#	MODIFIED:	03-Mar-2011 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Gatk::GermlineSnv {
	is => 'Genome::Model::Tools::Gatk',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		bam_file	=> { is => 'Text', doc => "BAM File for Sample", is_optional => 0, is_input => 1 },
		vcf_output_file     => { is => 'Text', doc => "Output file to receive GATK vcf format lines", is_optional => 0, is_input => 1, is_output => 1 },
		verbose_output_file     => { is => 'Text', doc => "STDOUT from GATK", is_optional => 1, is_input => 1, is_output => 1 },
		dbSNP_version     => { is => 'Text', doc => "Version of dbSNP bed file to use", is_optional => 1, is_input => 1, is_output => 1, default => 130 },
		gatk_params => { is => 'Text', doc => "Parameters for GATK", is_optional => 1, is_input => 1, is_output => 1, default => "-T UnifiedGenotyper -et NO_ET" },
		reference_fasta => { is => 'Text', doc => "Parameters for GATK", is_optional => 1, is_input => 1, is_output => 1, default => "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa" },
		run_unsafe_mode => { is => 'Text', doc => "Make GATK print errors instead of dying", is_optional => 1, is_input => 1, default => 1 },
	        mb_of_ram => {
	            is => 'Text',
        	    doc => 'The amount of RAM to use, in megabytes -- if you change this higher, must change lsf_resource to match',
        	    default => 5000,
	        },
		skip_if_output_present => { is => 'Text', doc => "Skip if output is present", is_optional => 1, is_input => 1},
	],
    # Make workflow choose 64 bit blades
    has_param => [
        lsf_queue => {
            default_value => 'long'
        }, 
        lsf_resource => {
            default_value => "-R 'rusage[mem=6000] select[type==LINUX64 && model != Opteron250 && mem>6000 && maxtmp>100000] span[hosts=1]' -M 6000000",
        },
    ],

};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Runs the GATK germline indel detection pipeline"                 
}

sub help_synopsis {
    return <<EOS
This command runs the GATK indel detection pipeline
EXAMPLE:	gmt gatk germline-indel bam-file file.bam --output-file GATK.indel --bed-output-file GATK.indel.bed
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


        # note: this crashes when using v. 131, because the 'get' returns undef
#        my $model = Genome::Model->get( name => "dbSNP-human-".$self->dbSNP_version);
	my $dbsnp_rod;
	if ($self->dbSNP_version == 130) {
#		$dbsnp_bed = '/gscuser/dkoboldt/SNPseek/SNPseek2/ucsc/snp130.txt';
#		$dbsnp_rod = '/gscmnt/sata199/info/gatk_read_recalibration/dbsnp_130_b36.rod';
		$dbsnp_rod = '/gscmnt/sata424/info/medseq/Freimer-Boehnke/dbsnp/dbsnp_130_b36_plusNT.rod';
	}
#	elsif ($self->dbSNP_version == 132) {
#		$dbsnp_bed = '/gscuser/dkoboldt/SNPseek/SNPseek2/ucsc/hg19/snp131.txt';
#	}

# [-L targets.interval_list]

	## Run GATK ##
	my $path_to_gatk = $self->gatk_path;
	my $version = $self->version;
	my $gatk_params;
	if ($version le 5500) {
		$gatk_params = $self->gatk_params;
	}
	elsif ($version ge 5500) {
		$gatk_params = $self->gatk_params .  " -glm SNP";
	}
	else {
		die "cannot determine gatk version to set proper parameter names";
	}
	my $reference_fasta = "-R " . $self->reference_fasta;
	my $output_file = "-o " . $self->vcf_output_file;	
	my $bam_input = "-I ".$self->bam_file;
	my $ram = $self->mb_of_ram;
	my $dbsnp = "-D " . $dbsnp_rod;
	my $stdout_call_stats = "-l INFO";
	my $cmd = 'java -Xms'.$ram.'m -Xmx'.$ram.'m -jar ';
	$cmd .= join(" ", $path_to_gatk, $gatk_params, $stdout_call_stats, $reference_fasta, $dbsnp, $bam_input, $output_file);
	
	## Optionally append stdout output file ##

	if($self->verbose_output_file) {
		$cmd .= " -verbose " . $self->verbose_output_file;
	}

	## Optionally run in unsafe mode ##

	if($self->run_unsafe_mode) {
		$cmd .= " -U ALL --validation_strictness SILENT ";
	}


	## Run GATK Command ##
	my $return;
	if($self->skip_if_output_present && -e $output_file)
	{
		
	}
	else
	{
		Genome::Sys->write_file($self->vcf_output_file, ''); # This will create an empty output file to help prevent GATK from crashing
		$return = Genome::Sys->shellcmd(
                           cmd => "$cmd",
                           output_files => [$self->vcf_output_file],
                           skip_if_output_is_present => 0,
                       );
		unless($return) { 
			$self->error_message("Failed to execute GATK: GATK Returned $return");
			die $self->error_message;
		}
	}

	return $return;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

1;

