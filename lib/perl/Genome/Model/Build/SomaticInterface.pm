package Genome::Model::Build::SomaticInterface;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::SomaticInterface {
    is => 'UR::Object',
    is_abstract => 1,
};

sub reference_sequence_build {
    my $self = shift;
    $self->fatal_message('Abstract: (reference_sequence_build) needs to be defined on class (%s)', $self->class);
}

sub individual {
    my $self = shift;

    my $subject = $self->subject;
    if ($subject->isa('Genome::Sample')) {
        return $self->subject->individual;
    }
    elsif ($subject->isa('Genome::Individual')) {
        return $self->subject;
    }
    else {
        $self->fatal_message("Can't resolve subject for build (%s) of class (%s)", $self->id, $self->class);
    }
}

sub individual_common_name {
    my $self = shift;
    return $self->individual->common_name;
}

sub has_snvs_annotated_variants_vcf_file {
    my $self = shift;
    $self->fatal_message('Abstract: (has_snvs_annotated_variants_vcf_file) needs to be defined on class (%s)', $self->class);
}

sub indels_effects_file {
    my $self = shift;
    my $tier = shift;
    $self->fatal_message('Abstract: (indels_effects_file) needs to be defined on class (%s)', $self->class);
}

sub snvs_effects_file {
    my $self = shift;
    my $tier = shift;
    return File::Spec->join($self->data_directory, 'effects', "snvs.hq.novel.$tier.v2.bed");
}

sub ran_copycat {
    my $self = shift;

    if (not -s File::Spec->join($self->data_directory, 'variants', 'cnvs.hq')
        and glob(File::Spec->join($self->data_directory, 'variants', 'cnv', 'copy-cat*')))
    {
        return 1;
    }
    else {
        return 0;
    }
}

sub has_microarray_build {
    my $self = shift;
    $self->fatal_message('Abstract: (has_microarray_build) needs to be defined on class (%s)', $self->class);
}

sub normal_bam {
    my $self = shift;
    $self->fatal_message('Abstract: (normal_bam) needs to be defined on class (%s)', $self->class);
}

sub tumor_bam {
    my $self = shift;
    $self->fatal_message('Abstract: (tumor_bam) needs to be defined on class (%s)', $self->class);
}

sub build_for_normal_bam {
    my $self = shift;
    $self->fatal_message('Abstract: (build_for_normal_bam) needs to be defined on class (%s)', $self->class);
}

sub build_for_tumor_bam {
    my $self = shift;
    $self->fatal_message('Abstract: (build_for_tumor_bam) needs to be defined on class (%s)', $self->class);
}

1;
