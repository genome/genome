package Genome::Model::Tools::BioSamtools::ProgressionInstance;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::ProgressionInstance {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        bam_files => {
            is => 'Array',
            doc => 'A list of bam files to merge and run ref-cov',
        },
        target_query_file => {
            is => 'Text',
            doc => 'a file of query names and coordinates relative to the bam targets',
        },
        output_directory => {
            is => 'Text',
            doc => 'The base output directory for ref-cov files',
        },
        samtools_version => {
            is => 'Text',
            doc => 'The version of samtools to use',
            default_value => Genome::Model::Tools::Sam->default_samtools_version,
        }
    ],
    has_output => [
        stats_file => {
            is => 'String',
            is_optional => 1,
        },
        bias_basename => {
            is => 'String',
            is_optional => 1,
        },
        instance => {
            is => 'String',
            is_optional => 1,
        },
    ],
    has_param => [
        lsf_resource => {
            is_optional => 1,
            default_value => "-R 'select[type==LINUX64]'",
        },
    ],
};

sub create {
    my $class = shift;
    my %params = @_;
    my $bam_files = delete($params{bam_files});
    my $self = $class->SUPER::create(%params);
    return unless $self;
    $self->bam_files($bam_files);
    return $self;
}

sub execute {
    my $self = shift;

    my @bam_files = @{$self->bam_files};
    my $instance = scalar(@bam_files);
    my $merged_bam = Genome::Sys->create_temp_file_path('merged_'. $instance .'.bam');
    my %params = (
        merger_name => 'samtools',
        is_sorted => 1,
        files_to_merge => \@bam_files,
        merged_file => $merged_bam,
        use_version => $self->samtools_version,
    );
    my $merge = Genome::Model::Tools::Sam::Merge->execute(%params);
    unless ($merge) {
        $self->error_message('Failed to execute bam file merge tool with params '. Data::Dumper::Dumper(%params));
        die($self->error_message);
    }
    Genome::Sys->create_directory($self->output_directory);
    $self->instance(scalar(@bam_files));
    $self->stats_file($self->output_directory .'/STATS_'. $self->instance .'.tsv');
    $self->bias_basename($self->output_directory .'/bias_'.$self->instance);
    my $cmd = Genome::Model::Tools::BioSamtools::RelativeCoverage->execute(
        bam_file => $merged_bam,
        bed_file => $self->target_query_file,
        stats_file => $self->stats_file,
        bias_file => $self->bias_basename,
    );
    unless ($cmd) {
        die('Failed to execute relative reference coverage.');
    }
    return 1;
}

1;
