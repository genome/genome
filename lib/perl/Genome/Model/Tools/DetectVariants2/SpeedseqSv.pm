#SpeedseqSv.pm


package Genome::Model::Tools::DetectVariants2::SpeedseqSv;
use warnings;
use strict;

use Genome;
use File::Basename;
use IPC::System::Simple;
use Genome::File::Tsv;

class Genome::Model::Tools::DetectVariants2::SpeedseqSv {
    is => 'Genome::Model::Tools::DetectVariants2::Detector',
    has_optional => [
        _legend_output => {
            calculate_from => ['_temp_staging_directory'],
            calculate => q{ File::Spec->join($_temp_staging_directory, 'legend.tsv'); },
        },
    ]
};

	my $test_dir = "/gscmnt/gc2801/analytics/mfulton/test1";
	my $temp_directory = Genome::Sys->create_temp_directory();
	my $data_dir = '/gscmnt/gc2801/analytics/mfulton/genome2/lib/perl/Genome/Test/Data.pm.d/NA12878';
	my $bam = "$data_dir/NA12878.20slice.30X.aligned.bam";
	my $split_bam = "$data_dir/NA12878.20slice.30X.splitters.bam";
	my $discordant_bam = "$data_dir/NA12878.20slice.30X.discordants.bam";
	my $reference_fasta = "$data_dir/human_g1k_v37_20_42220611-42542245.fasta";
	
	my $pkg = 'Genome::Model::Tools::Speedseq::Sv';

sub _detect_variants {
	my $self = shift;

	my $cmd = $pkg->create(
		version => 'test',
		temp_directory => $temp_directory,
	   	reference_fasta => $reference_fasta,
   		full_bam_file => $bam,
   		output_prefix => $self->_sv_staging_output,
   		CNVnator_read_depth => 'true',
		split_read_bam_file => $split_bam,
		genotype_svtyper => 1,
		discordant_read_bam_file => $discordant_bam,
	);
	$cmd->execute();
};
