package Genome::Model::Tools::AmpliconAssembly::Classify;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Tools::AmpliconAssembly::Classify {
    is => 'Genome::Model::Tools::AmpliconAssembly',
    # TODO classifier and params
};

#< Helps >#
sub help_detail {
    return sprintf(
        'Classifies successfully assembled amplicons (%s) using RDP.',
        Genome::Model::Tools::AmpliconAssembly::Amplicon->successfully_assembled_requirements_as_string,
    );
}

sub help_synopsis {
}

#< Command >#
sub sub_command_sort_position { 41; }

sub execute {
    my $self = shift;

    my $amplicons = $self->get_amplicons
        or return;

    require Genome::Utility::MetagenomicClassifier::Rdp::Version2x1;
    my $classifier = Genome::Utility::MetagenomicClassifier::Rdp::Version2x1->new()
        or return;

    for my $amplicon ( @$amplicons ) {
        my $bioseq = $amplicon->get_bioseq
            or next;
        my $classification = $classifier->classify($bioseq);
        unless ( $classification ) {
            $self->error_message(
                sprintf(
                    'Can\'t get classification from RDP classifier for amplicon (%s)', 
                    $amplicon->get_name,
                )
            );
            next;
        }
        $amplicon->save_classification($classification); # error?
    }

    return 1;
}

1;

#$HeadURL$
#$Id$
