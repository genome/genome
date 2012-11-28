package Genome::SoftwareResult::Command::FillLookupHashes;

use strict;
use warnings;

use Genome;


class Genome::SoftwareResult::Command::FillLookupHashes {
    is => 'Command::V2',

    has_optional => [
        commit_size => {
            is => 'Number',
            default_value => '10000',
        },
    ],
};

sub execute {
    my $self = shift;

    my $i = 0;
    my $iterator = Genome::SoftwareResult->create_iterator(
        'lookup_hash' => undef);

    while (my $sr = $iterator->next()) {
        $sr->lookup_hash($sr->calculate_lookup_hash);
        $i++;
        if ($i >= $self->commit_size) {
            $self->status_message("Committing chunk of $i rows");
            UR::Context->commit();
            $i = 0;
            $iterator = Genome::SoftwareResult->create_iterator(
                'lookup_hash' => undef);
        }
    }

    return 1;
}

1;
