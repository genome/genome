package Genome::Model::Build::DifferentialExpression;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::DifferentialExpression {
    is => 'Genome::Model::Build',
};

sub transcript_convergence_directory {
    my $self = shift;
    return $self->data_directory . '/transcript_convergence';
}

sub transcript_gtf_prefix {
    my $self = shift;
    if ($self->model->transcript_convergence_name eq 'cuffcompare') {
        return $self->transcript_convergence_directory .'/cuffcompare';
    } else {
        #TODO: Should this ever be called by anything else?
        die();
    }
}

sub transcript_gtf_file_path {
    my $self = shift;
    if ($self->model->transcript_convergence_name eq 'cuffcompare') {
        return $self->transcript_gtf_prefix .'.combined.gtf';
    } else {
        return $self->transcript_convergence_directory .'/merged.gtf';
    }
}

sub differential_expression_directory {
    my $self = shift;
    return $self->data_directory . '/differential_expression';
}

sub summarize_differential_expression_directory {
    my $self = shift;
    return $self->data_directory .'/summarize_differential_expression';
}

sub summary_report_pdf_file {
    my $self = shift;
    return $self->summarize_differential_expression_directory .'/summary.pdf';
}

sub _disk_usage_result_subclass_names {
    my $self = shift;

    return [];
}

1;

