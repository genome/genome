package Genome::Model::Tools::Sx::Split::ByNs;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Split::ByNs {
    is  => 'Genome::Model::Tools::Sx::Base',
    has => {
        number_of_ns => {
            is => 'Integer',
            doc => 'Number of Ns to split sequence.',
        },
    },
};

sub help_brief { 'Split sequences by Ns' }

sub help_detail { 'Split sequences by a set number of Ns.' }
    
sub execute {
    my $self = shift;

    return if not $self->_init;

    my $reader = $self->_reader;
    while ( my $seqs = $reader->read ) {
        for my $seq ( @$seqs ) {
            my $it = $self->iterator_to_split_sequence($seq);
            while ( my $new_seq = $it->() ) { $self->_writer->write([ $new_seq ]); }
        }
    }

    return 1;
}

sub iterator_to_split_sequence {
    my ($self, $seq) = @_;

    my ($pos, $split_bases, $ns);
    my $string_o_ns = 'N' x $self->number_of_ns;
    my $remaining_bases = uc($seq->{seq});
    my $remaining_quals = $seq->{qual};
    my $cnt = 0;
    return sub{
        while ( $remaining_bases && length $remaining_bases ) {
            my $pos = index($remaining_bases, $string_o_ns);
            if ( $pos == -1 ) { # not enough Ns, use full sequence
                $split_bases = $remaining_bases;
                $ns = '';
                $remaining_bases = undef;
            }
            else { # match
                # Remove the split bases using the position of the index match
                $split_bases = substr($remaining_bases, 0, $pos, '');
                # Get the Ns using regex
                ($ns) = $remaining_bases =~ /^(N+)/;
                # Remove Ns
                substr($remaining_bases, 0, length($ns), '');
            }

            my %split_seq = (
                id => join('.', $seq->{id}, ++$cnt),
                seq => $split_bases,
            );
            if ( $remaining_quals ) {
                # Get the split quals
                $split_seq{qual} = substr($remaining_quals, 0, length($split_seq{seq}), '');
                # Remove the quals for the Ns
                substr($remaining_quals, 0, length($ns), '');
            }

            return \%split_seq;
        }
    };
}

1;

