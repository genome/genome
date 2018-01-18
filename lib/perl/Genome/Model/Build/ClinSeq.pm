package Genome::Model::Build::ClinSeq;
use strict;
use warnings;

use Genome;
use List::MoreUtils qw(uniq);

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
        if (my $input_build = $self->$accessor) {
            push @input_builds, $input_build;
            for my $build (qw(tumor_build normal_build)) {
                if ($input_build->can($build)) {
                    push @input_builds, $input_build->$build;
                }
            }
        }
    }

    return grep {defined($_)} @input_builds;
}

sub reference_sequence_build {
    my $self = shift;
    my @references = $self->_infer_references_from_input_builds;
    if (scalar(@references) == 0) {
        $self->fatal_message("No reference builds on input models?");
    }
    return $references[0];
}

sub _infer_references_from_input_builds {
    my $self = shift;

    my @references;
    for my $input_build ($self->input_builds) {
        if ($input_build->can('reference_sequence_build')) {
            push @references, $input_build->reference_sequence_build;
        }
    }
    return sort {$a->id cmp $b->id} uniq @references;
}

sub best_somatic_build {
    my $self = shift;
    if ($self->wgs_build) {
        return $self->wgs_build;
    }
    elsif ($self->exome_build) {
        return $self->exome_build;
    }
}

sub best_somatic_build_subject_common_name {
    my $self = shift;
    return $self->best_somatic_build->subject->common_name;
}

sub has_microarray_build {
    my $self = shift;

    for my $build_accessor (qw(exome_build wgs_build)) {
        if (my $build = $self->$build_accessor) {
            return 1 if ($build->has_microarray_build);
        }
    }

    return 0;
}

sub tumor_microarray_build {
    my $self = shift;

    for my $build_accessor (qw(exome_build wgs_build)) {
        if (my $build = $self->$build_accessor) {
            return $build->tumor_microarray_build if ($build->has_microarray_build);
        }
    }
}

sub normal_microarray_build {
    my $self = shift;

    for my $build_accessor (qw(exome_build wgs_build)) {
        if (my $build = $self->$build_accessor) {
            return $build->normal_microarray_build if ($build->has_microarray_build);
        }
    }
}

sub _disk_usage_result_subclass_names {
    my $self = shift;

    return [qw(
        Genome::Model::ClinSeq::Command::AnnotateSnvsVcf::Result
    )];
}

1;
