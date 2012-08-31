package Genome::Model::Tools::Sx::PhredQualReader;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::PhredQualReader {
    is => 'Genome::Model::Tools::Sx::PhredReaderBase',
};

sub read {
    my $self = shift;

    my ($id, $desc, $data) = $self->_parse_io($self->{_qual_file});
    return if not defined $id;

    my $qual;
    for my $line ( split("\n", $data) ){
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        $qual .= join('', map { chr($_ + 33) } split(/\s+/, $line));
    }

    if ( not $qual ) {
        Carp::confess("Could not convert phred quality to sanger: $data");
    }

    return {
        id => $id,
        desc => $desc,
        qual => $qual,
    };
}

1;

