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
sub default_processing_profile_ids {
    # RT66900 from 2278045 to 2571784 
    # RT85266 add 2752939
    return ( 
        '2571784',# RDP 2.2 set 6
        '2752939',# RDP 2.5 set 9
    );
}

sub default_processing_profile_id {
    my @default_processing_profile_ids = default_processing_profile_ids();
    return $default_processing_profile_ids[0];
}
#<>#

sub build_subclass_name {
    return 'metagenomic-16s-composition';
}

sub _additional_parts_for_default_name {
    my $self = shift;

    my @parts;
    my $subject = $self->subject;
    if ( $subject->isa('Genome::Sample') and defined $subject->tissue_desc ) {
        push @parts, $subject->tissue_desc;
    }

    push @parts, $self->processing_profile->classifier;

    return @parts;
}

sub is_for_qc {
    my $self = shift;
    return 1 if $self->name =~ /\-qc$/;
    return 0;
}

1;

