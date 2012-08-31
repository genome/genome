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
    
    
    my $cmd = 'genome-perl5.10 `which gmt` ref-cov cluster-coverage --bam-file='. $self->bam_file .' --minimum-zenith='. $self->zenith_depth .' --minimum-depth='. $self->minimum_depth .' --stats-file='. $self->stats_file .' --bed-file='. $self->bed_file ;
    
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->bam_file],
        output_files => [$self->bed_file,$self->stats_file],
        skip_if_output_is_present => 0,
    );
    
    return 1;   
}

1;

__END__

genome-perl5.10 `which gmt` bio-samtools cluster-coverage 
--bam-file=/gscmnt/sata141/techd/jhundal/miRNA/64LY0AAXX_AML/NEW_FAR/Lane3_filtered.bam 
--minimum-zenith=5 
--bed-file=/gscmnt/sata141/techd/jhundal/miRNA/64LY0AAXX_AML/NEW_FAR/Lane3_zenith5.bed 
--stats-file=/gscmnt/sata141/techd/jhundal/miRNA/64LY0AAXX_AML/NEW_FAR/Lane3_coverage_stats.tsv
