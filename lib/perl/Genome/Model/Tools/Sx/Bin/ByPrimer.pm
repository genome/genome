package Genome::Model::Tools::Sx::Bin::ByPrimer;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use Regexp::Common;

class Genome::Model::Tools::Sx::Bin::ByPrimer {
    is  => 'Genome::Model::Tools::Sx::Base',
    has => [
        primers => {
            is  => 'Text',
            is_many => 1,
            doc => 'Names and primer sequences to bin sequences.',
        },
        remove => {
            is => 'Boolean',
            doc => 'If found, remove the primer from the sequence.',
        },
    ],
};

sub help_brief {
    return 'Bin sequences by primers';
}

sub help_synopsis {
    return <<HELP
HELP
}

sub help_detail {
    return <<HELP 
HELP
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;

    my %primers;
    for my $primer ( $self->primers ) {
        my ($name, $sequence) = split('=', $primer);
        if ( not $name or $name eq '' ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ primers /],
                desc => "Primer ($primer) does not have a name. Primers should have format of name=sequence",
            );
        }
        if ( not $sequence or $sequence eq '' ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ primers /],
                desc => "Primer ($primer) does not have a sequence. Primers should have format of name=sequence",
            );
        }
        push @{$primers{$name}}, $sequence;
    }

    $self->{_primers} = \%primers;

    return @errors;
}

sub execute {
    my $self = shift;

    my $init = $self->_init;
    return if not $init;

    my $reader = $self->_reader;
    my $writer = $self->_writer;

    my $binner = ( $self->remove )
    ? $self->_create_bin_by_primer_and_remove
    : $self->_create_bin_by_primer;

    while ( my $seqs = $reader->read ) {
        for my $seq ( @$seqs ) {
            $binner->($seq);
        }
        $seqs = $self->_update_writer_for_paired_reads( $seqs );
        $writer->write($seqs);
    }

    return 1;
}

sub _update_writer_for_paired_reads {
    my ( $self, $seqs ) = @_;

    return $seqs unless scalar @$seqs == 2;

    #discard paired reads w/ only 1 hit or 2 hits from diff primer sets
    unless( exists @$seqs[0]->{writer_name} and exists @$seqs[1]->{writer_name} and @$seqs[0]->{writer_name} eq @$seqs[1]->{writer_name} ) {
        delete @$seqs[0]->{writer_name};
        delete @$seqs[1]->{writer_name};
        return $seqs;
    }
    #discard paired reads where fwd/rev hits are from same fwd or rev primer set
    if( @$seqs[0]->{primer_name} eq @$seqs[1]->{primer_name} ) {
        delete @$seqs[0]->{writer_name};
        delete @$seqs[1]->{writer_name};
        return $seqs;
    }

    return $seqs;
}

sub _create_bin_by_primer {
    my $self = shift;

    my $primers = $self->{_primers};

    return sub{
        my $seq = shift;
        for my $primer_name ( keys %$primers ) {
            my $writer_name = $primer_name;
            $writer_name =~ s/\.[FR]$//;
            for my $primer ( @{$primers->{$primer_name}} ) {
                if ( $seq->{seq} =~ /^$primer/ ) {
                    $seq->{primer_name} = $primer_name;
                    $seq->{writer_name} = $writer_name;
                    return 1;
                }
            }
        }
    };
}

sub _create_bin_by_primer_and_remove {
    my $self = shift;

    my $primers = $self->{_primers};

    return sub{
        my $seq = shift;
        for my $primer_name ( keys %$primers ) {
            my $writer_name = $primer_name;
            $writer_name =~ s/\.[FR]$//;
            for my $primer ( @{$primers->{$primer_name}} ) {
                if ( $seq->{seq} =~ s/^$primer// ) {
                    substr($seq->{qual}, 0, length($primer), '') if $seq->{qual}; # rm qual
                    $seq->{primer_name} = $primer_name;
                    $seq->{writer_name} = $writer_name;
                    return 1;
                }
            }
        }
    };
}

1;

