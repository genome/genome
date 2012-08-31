package Genome::InstrumentData::Command::AlignmentResult::List;
use strict;
use warnings;
use Genome;

# TODO: the whole ::AlignmentResult::Command tree should appear in the right location so we don't need this empty subclass

class Genome::InstrumentData::Command::AlignmentResult::List {
    is => 'Genome::InstrumentData::AlignmentResult::Command::List',
};

1;

