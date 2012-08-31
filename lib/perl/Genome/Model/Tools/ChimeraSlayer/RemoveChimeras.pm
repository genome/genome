package Genome::Model::Tools::ChimeraSlayer::RemoveChimeras;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ChimeraSlayer::RemoveChimeras {
    is => 'Command::V2',
    has => [
        sequences => {
            is => 'Text',
            doc => 'The sequences to remove chimeras. Type is determined by extension. For more options, see "gmt sx -h" for format help.',
        },
        chimeras => {
            is => 'Text',
            doc => 'The file of chimeras.',
        },
        output => {
            is => 'Text',
            doc => 'The output to write chimera free sequences. Type is determined by extension. For more options, see "gmt sx -h" for format help.',
        },
        _metrics => { is_optional => 1, is_transient => 1, },
    ],
};

sub help_brief {
    return 'remove chimeric sequences from a file';
}

sub help_detail {
    return <<HELP;
This command takes a sequences file (fasta, fastq, etc) and a chimeras file (CPC, final verdict) and outputs sequences that are not chimeric.
HELP
}

sub execute {
    my $self = shift;

    $self->status_message('Remove chimeras...');

    $self->status_message('Sequences: '.$self->sequences);
    my ($reader_class, $reader_params) = Genome::Model::Tools::Sx::Reader->parse_reader_config($self->sequences);
    return if not $reader_class;
    my $reader = $reader_class->create(%$reader_params);
    return if not $reader;

    $self->status_message('Chimeras: '.$self->chimeras);
    my $chimera_reader = Genome::Model::Tools::ChimeraSlayer::Reader->create(input => $self->chimeras);
    return if not $chimera_reader;

    $self->status_message('Output: '.$self->output);
    my ($writer_class, $writer_params) = Genome::Model::Tools::Sx::Writer->parse_writer_config($self->output);
    my $writer = $writer_class->create(%$writer_params);
    return if not $writer;

    $self->status_message('Read sequences...');
    my %metrics = ( sequences => 0, chimeras => 0 );
    my $chimera = $chimera_reader->read;
    while ( my $seq = $reader->read ) {
        $metrics{sequences}++;
        if ( $seq->{id} eq $chimera->{id} ) {
            my $verdict = $chimera->{verdict}; # store verdict
            $chimera = $chimera_reader->read; # get next chimera
            if ( $verdict eq 'YES' ) { # do not write seq if it is a chimera
                $metrics{chimeras}++;
                next;
            }
        }
        $metrics{output}++;
        $writer->write($seq);
    }
    $self->status_message('Read sequences...OK');

    $self->_metrics(\%metrics);
    $self->status_message('Sequence count: '.$metrics{sequences});
    $self->status_message('Chimera count: '.$metrics{chimeras});
    $self->status_message('Output count: '.$metrics{output});

    $self->status_message('Remove chimeras...OK');

    return 1;
}

1;

