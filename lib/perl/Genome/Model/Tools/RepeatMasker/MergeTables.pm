package Genome::Model::Tools::RepeatMasker::MergeTables;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RepeatMasker::MergeTables {
    is => 'Genome::Model::Tools::RepeatMasker::TableI',
    has => [
        input_tables => { },
    ],
};

sub execute {
    my $self = shift;
    my %repeats;
    for my $table ( @{$self->input_tables} ) {
        my $table_fh = IO::File->new($table,'r');
        my $family;
        while ( my $line = $table_fh->getline ) {
            chomp($line);
            my @entry = split("\t",$line);
            if (scalar(@entry) == 2) {
                if ($entry[0] eq 'sequences:') {
                    my $total_count = $self->_total_count || 0;
                    $self->_total_count($total_count + $entry[1]);
                } elsif ($entry[0] eq 'total length:') {
                    my $total_bp = $self->_total_bp || 0;
                    $self->_total_bp($total_bp + $entry[1]);
                }
                # new repeat family
            } elsif (scalar(@entry) == 4) {
                $family = $entry[0];
                $family =~ s/://;
                $repeats{$family}{elements} += $entry[1];
                $repeats{$family}{base_pair} += $entry[2];
            } elsif (scalar(@entry) == 5) {
                unless ($family) { die; }
                my $class = $entry[1];
                $class =~ s/://;
                $repeats{$family}{$class}{elements} += $entry[2];
                $repeats{$family}{$class}{base_pair} += $entry[3];
            }
        }
        $table_fh->close;
    }
    $self->print_table_from_hash_ref(\%repeats);
    return 1;
}

1;
