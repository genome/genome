package Genome::Model::Tools::SmrtAnalysis::ConsensusStats;

use strict;
use warnings;

use Genome;

my $DEFAULT_LSF_RESOURCE = "-g /pacbio/smrtanalysis -M 8000000 -R 'select[type==LINUX64 && mem>=8000 && tmp>=40000] rusage[mem=8000,tmp=20000]'";

class Genome::Model::Tools::SmrtAnalysis::ConsensusStats {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        alignment_summary_gff_file => { },
        variants_gff_file => { },
        cmp_hdf5_file => { },
    ],
    has_optional_input => [
        output_gff_file => { is_output => 1 },
    ],
    has_optional_param => [
        lsf_resource => { default_value => $DEFAULT_LSF_RESOURCE },
    ],
};

sub execute {
    my $self = shift;

    my $output_gff_file = $self->output_gff_file;
    unless ($self->output_gff_file) {
        my (undef, $tmp_file) = File::Temp::tempfile('alignment_summary-XXXXXXX', DIR => '/gscmnt/gc2123/production/lsf_shared_dir', SUFFIX => '.gff');
        $output_gff_file = $tmp_file;
    }
    my $cmd = $self->analysis_bin .'/jConsensusStats '. $self->alignment_summary_gff_file .' '. $self->variants_gff_file .' '. $self->cmp_hdf5_file .' > '. $output_gff_file;
    my @input_files = ($self->alignment_summary_gff_file,$self->variants_gff_file,$self->cmp_hdf5_file);
    $self->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => [$output_gff_file],
        skip_if_output_is_present => 0,
    );
    unless ($self->output_gff_file) {
        #There is some sed stuff going on in the PacBio software, but it doesn't appear necessary. let's just copy for now
        unlink($self->alignment_summary_gff_file) || die('Failed to remove existing alignment summary '. $self->alignment_summary_gff_file);
        File::Copy::copy($output_gff_file,$self->alignment_summary_gff_file) || die ('Failed to copy '. $output_gff_file .' over existing alignment summary '. $self->alignment_summary_gff_file);
        $self->output_gff_file($self->alignment_summary_gff_file);
    }
    return 1;
}
