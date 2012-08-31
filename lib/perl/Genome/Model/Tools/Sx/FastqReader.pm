package Genome::Model::Tools::Sx::FastqReader;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::FastqReader {
    is => 'Genome::Model::Tools::Sx::SeqReader',
};

sub read {
    my $self = shift;

    my $line = $self->{_file}->getline;
    return if not defined $line;
    chomp $line;

    my ($id, $desc) = split(/\s+/, $line, 2);
    if ( not $id =~ s/^@// ) {
        Carp::confess('Id line does not start with an "@". On line: '.$line);
    }

    my $seq = $self->{_file}->getline;
    chomp $seq; 

    my $qual_header = $self->{_file}->getline; 
    if ( $qual_header !~ /^\+/ ) {
        Carp::confess('Quality header does not start with a "+" for '.$line);
    }
    
    my $qual = $self->{_file}->getline;
    chomp $qual;

    if ( length $seq != length $qual ) {
        Carp::confess( join("\n", 'Sequence and quality lengths differ:', $line, $seq, $qual) );
    }

    return {
        id => $id,
        desc => $desc,
        seq => $seq,
        qual => $qual,
    };
}

1;

