package Genome::InstrumentData::Command::AlignmentResult::Remove;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::AlignmentResult::Remove {
    is => 'Command::V2',
    has_input => [
        alignment_results => {
            is => 'Genome::InstrumentData::AlignmentResult',
            shell_args_position => 1,
            is_many => 1,
            doc => 'instrument data to remove, specified by id or expression'
        },
    ],
    doc => 'remove an alignment result',
};

sub help_synopsis {
    return <<EOS;
genome instrument-data alignment-result remove 12345                                    # by id
genome instrument-data alignment-result remove instrument_data_id=6789                  # everything for some instrument data 
genome instrument-data alignment-result remove aligner_name=bwa,aligner_version=0.5.9   # everything for some aligner and version
EOS
}

sub help_detail {
    my $self = shift;
    return <<EOS;
Removes one or more alignment results from the system, including deallocation of disk.

Note that if the alignment results have dependent builds the deletion and disk deallocation will fail, so old builds may need to be removed for this to work.
EOS
}

sub execute {
    my $self = shift;
    my @i = $self->alignment_results();
    $self->status_message("Removing " . scalar(@i) . " alignment result entries...");
    sleep 5;
    for my $i (@i) {
        $self->status_message("deleting " . $i->__display_name__ . "...");
        $i->delete;
        print "$i\n";
    }
    $self->status_message("deletion complete.");
    return 1;
}

1;

