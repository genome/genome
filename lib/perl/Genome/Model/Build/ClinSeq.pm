package Genome::Model::Build::ClinSeq;
use strict;
use warnings;

use Genome;
class Genome::Model::Build::ClinSeq {
    is => ['Genome::Model::Build',
           'Genome::Model::Build::ClinSeq::FileAccessors',
           'Genome::Model::Build::ClinSeq::InputBuilds'
          ],
};

sub special_compare_functions {
    my $self = shift;
    return ($self->SUPER::special_compare_functions,
        qr(circos\.conf$) => sub {!Genome::Model::Build::ClinSeq::diff_circos_conf(@_)} );
}

sub diff_circos_conf {
    my ($first_file, $second_file) = @_;
    Carp::confess('Missing files to diff!') if @_ != 2;
    my $first_md5  = qx(grep -vP '^file\\s+= ' $first_file | md5sum);
    my $second_md5 = qx(grep -vP '^file\\s+= ' $second_file | md5sum);
    return ($first_md5 eq $second_md5 ? 1 : 0);
}

sub should_run_exome_cnv {
    my $self = shift;
    return ($self->exome_build and $self->processing_profile->exome_cnv);
}

sub input_builds {
    my $self = shift;

    my @input_builds = ( $self->tumor_rnaseq_build, $self->normal_rnaseq_build );
    for my $accessor (qw(wgs_build exome_build)) {
        if ($self->$accessor) {
            push @input_builds, $self->$accessor;
            for my $build (qw(tumor_build normal_build)) {
                if ($self->$accessor->can($build)) {
                    push @input_builds, $self->$accessor->$build;
                }
            }
        }
    }

    return grep {defined($_)} @input_builds;
}

sub reference_sequence_build {
    my $self = shift;
    return $self->model->_resolve_reference;
}

1;
