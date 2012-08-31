package Genome::Model::MetagenomicComposition16s;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicComposition16s {
    is => 'Genome::ModelDeprecated',
    has => [
    map({
            $_ => {
                via => 'processing_profile',
            }
        } Genome::ProcessingProfile::MetagenomicComposition16s->params_for_class
    ),
    ],
};

#< Default Processing Profile >#
# RT66900 was 2278045
# Moved from AQID
sub default_processing_profile_id {
    return 2571784;
}

sub default_processing_profile {
    return Genome::ProcessingProfile::MetagenomicComposition16s->get( __PACKAGE__->default_processing_profile_id );
}
#<>#

sub build_subclass_name {
    return 'metagenomic-16s-composition';
}

sub _additional_parts_for_default_model_name {
    my $self = shift;

    my @parts;
    my $subject = $self->subject;
    if ( $subject->isa('Genome::Sample') and defined $subject->tissue_desc ) {
        push @parts, $subject->tissue_desc;
    }

    return @parts;
}

sub is_for_qc {
    my $self = shift;
    return 1 if $self->name =~ /\-qc$/;
    return 0;
}

1;

