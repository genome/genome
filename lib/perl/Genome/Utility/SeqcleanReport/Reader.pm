package Genome::Utility::SeqcleanReport::Reader;

#REVIEW fdu 11/20/2009
#This can be part of the refactoring that uses G::U::IO::Reader to
#replace G::U::Parser, and place all current XXXX/Reader.pm under
#new namespace G::U::IO::XXXX::Reader

use strict;
use warnings;

use Genome;

my @header_fields = qw(accession pc_undetermined start end length trash_code comments);

class Genome::Utility::SeqcleanReport::Reader {
    is => 'Genome::Utility::Parser',
    has => [
            separator => {
                          default_value => "\t",
                      },
            header => {
                       default_value => '0',
                   },
            ],
};

sub header_fields {
    return \@header_fields;
}


1;
