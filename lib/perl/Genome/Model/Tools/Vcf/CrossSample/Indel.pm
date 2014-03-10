package Genome::Model::Tools::Vcf::CrossSample::Indel;

use strict;
use warnings;

use Genome;
use File::Basename qw(dirname);

class Genome::Model::Tools::Vcf::CrossSample::Indel {
    is => 'Genome::Model::Tools::Vcf::CrossSample::Base',

    has_optional_input => [
        varscan_version => {
            is => 'Text',
            doc => "Version of varscan to use, will be resolved to the latest default if not specified",
        },
        varscan_min_coverage => {
            is => 'Number',
            default => 3,
            doc => "Minimum base coverage to report readcounts",
        },
        varscan_min_avg_qual => {
            is => 'Number',
            default => 20,
            doc => "Minimum base quality to count a read",
        },
        varscan_min_var_freq => {
            is => 'Number',
            default => 0.20,
            doc => "Minimum variant allele frequency to call a variant",
        },
    ],
    has_transient_optional => [
        _template_path => {
            is => 'File',
        },
    ],
};

sub vcf_accessor {
    return "get_indels_vcf";
}

sub _execute {
    my $self = shift;

    $self->_resolve_varscan_version();

    my $process_uri;
    if ($self->roi_list) {
        $process_uri = $self->_execute_region_limited;
    } else {
        $process_uri = $self->_execute_normal;
    }

    return $process_uri;
}

sub _resolve_varscan_version {
    my $self = shift;
    unless (defined $self->varscan_version) {
        $self->varscan_version(Genome::Model::Tools::Varscan->default_version);
    }
    return;
}

sub get_template_path {
    my $filename = shift;
    (my $dirname = __FILE__) =~ s/.pm//g;
    return File::Spec->join($dirname, $filename);
}

sub _execute_region_limited {
    my $self = shift;

    $self->_template_path(get_template_path('RegionLimit.tsv'));
    my $inputs_filename = $self->_generate_inputs_file();

    my $source_path = 'Vcf::CrossSample::Indel::RegionLimit';
    return $self->_run_rex_process($source_path, $inputs_filename);
}

sub _execute_normal {
    my $self = shift;

    $self->_template_path(get_template_path('Normal.tsv'));
    my $inputs_filename = $self->_generate_inputs_file();

    my $source_path = 'Vcf::CrossSample::Indel';
    return $self->_run_rex_process($source_path, $inputs_filename);
}

sub bam_indexes {
    my $self = shift;

    my @bam_indexes;
    for my $bam ($self->bams) {
        my $index = "$bam.bai";
        unless (-s $index) {
            die $self->error_message("Bam index ($index) does not exist or is empty");
        }
        push @bam_indexes, $index;
    }

    return @bam_indexes;
}

sub bams {
    my $self = shift;

    my @bams;
    for my $build ($self->builds) {
        my $bam = $build->whole_rmdup_bam_file;
        unless (-s $bam) {
            die $self->error_message("Bam ($bam) does not exist or is empty");
        }
        push @bams, $bam;
    }

    return @bams;
}

1;
