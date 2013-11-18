package Genome::Model::SmallRna::Command::ClusterCoverage;

use strict;
use warnings;

use Genome;
use Workflow;
# added mim depth on Oct 26
# changed path of ClusterCoverage from GMT biosamtools to GMT RefCov

my $DEFAULT_ZENITH = '5';
my $DEFAULT_MIN_DEPTH = '1';

class Genome::Model::SmallRna::Command::ClusterCoverage {
    is => ['Genome::Model::SmallRna::Command::Base'],
    has_input => [
        bam_file => {
        	is => 'Text',
            doc => 'Input file of BAM alignments',
        },
        zenith_depth => {
            is => 'Text',
            is_output=> 1,
            doc => 'Minimum zepth depth cutoff',
            default_value => $DEFAULT_ZENITH,
        },
        minimum_depth => {
            is => 'Text',
            is_output=> 1,
            doc => 'Minimum depth to filter coverage',
            default_value => $DEFAULT_MIN_DEPTH,
        },
        stats_file => {
            is => 'Text',
            is_output=> 1,
            doc => 'Output file of coverage statistics ',
        },
        bed_file => {
            is => 'Text',
            is_output=> 1,
            doc => 'Output "Regions" BED file containing CLUSTERS from MERGED BED entries',
        },
    ],

    has_optional_param => [
        lsf_queue => {
            default_value => 'workflow',
        },
        lsf_resource => {
            default_value => '-R \'select[mem>16000] rusage[mem=16000]\' -M 16000000 ',
        },
    ],
};


sub execute {
    my $self = shift;
    my $bam_file   = $self->bam_file;
    my $stats_file = $self->stats_file;
    my $bed_file   = $self->bed_file;
    
    my $cmd = 'genome-perl5.10 -S gmt ref-cov cluster-coverage --bam-file='. $bam_file .' --minimum-zenith='. $self->zenith_depth .' --minimum-depth='. $self->minimum_depth .' --stats-file='. $stats_file .' --bed-file='. $bed_file;
    
    my $rv = Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files  => [$bam_file],
        skip_if_output_is_present => 0,
    );

    unless ($rv) {
        $self->error_message('Failed to execute command: '. $cmd);
        die $self->error_message;
    }

    my $err_msg = 'Probably caused by the cutoff failure of either zenith_depth or minimum_depth (or both)';

    unless (-s $stats_file) {
        $self->warning_message("Output stats file: $stats_file is not valid. $err_msg");
    }

    unless (-s $bed_file) {
        $self->warning_message("Output bed file: $bed_file is not valid. $err_msg");
    }
    
    return 1;   
}

1;

__END__

genome-perl5.10 -S gmt bio-samtools cluster-coverage 
--bam-file=/gscmnt/sata141/techd/jhundal/miRNA/64LY0AAXX_AML/NEW_FAR/Lane3_filtered.bam 
--minimum-zenith=5 
--bed-file=/gscmnt/sata141/techd/jhundal/miRNA/64LY0AAXX_AML/NEW_FAR/Lane3_zenith5.bed 
--stats-file=/gscmnt/sata141/techd/jhundal/miRNA/64LY0AAXX_AML/NEW_FAR/Lane3_coverage_stats.tsv
