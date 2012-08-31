package Genome::Utility::BioPerl;

use strict;
use warnings;

use Genome;

require Bio::SeqIO;
require Bio::Seq;
require Bio::Seq::Quality;
use Carp 'confess';
use Data::Dumper 'Dumper';

class Genome::Utility::BioPerl {
};

sub _create_bioseq_io {
    my ($self, $fh, $format) = @_;

    my $bioseq_io;
    eval{
        $bioseq_io = Bio::SeqIO->new(
            '-fh' => $fh,
            '-format' => $format,
        ); 
    };
    unless ( $bioseq_io ) {
        confess sprintf(
            "Failed to create %s Bio::SeqIO => %s",
            $format,
            $@,
        );
    }

    return $bioseq_io;
}

sub create_bioseq_writer {
    my ($self, $file, $format) = @_;

    my $fh = Genome::Sys->open_file_for_writing($file);
    unless ( $fh ) {
        confess "Can't open file for writing to create bioseq writer.  See above error.";
    }

    $format = 'fasta' unless defined $format;

    return $self->_create_bioseq_io($fh, $format);
}

sub create_bioseq_reader {
    my ($self, $file, $format) = @_;

    my $fh = Genome::Sys->open_file_for_reading($file);
    unless ( $fh ) {
        confess "Can't open file for reading to create bioseq reader.  See above error.";
    }

    $format = 'fasta' unless defined $format;

    return $self->_create_bioseq_io($fh, $format);
}

sub create_bioseq_from_fasta_and_qual {
    my ($self, %params) = @_;

    $self->validate_fasta_and_qual_bioseq($params{fasta}, $params{qual}); # confesses on error
    
    my $bioseq;
    eval {
        $bioseq = Bio::Seq::Quality->new(
            '-id' => $params{fasta}->id,
            '-desc' => $params{fasta}->desc,
            '-alphabet' => 'dna',
            '-force_flush' => 1,
            '-seq' => $params{fasta}->seq,
            '-qual' => $params{qual}->qual,
        ),
    };

    if ( $@ ) {
        $self->error_message(
            "Can't create combined fasta/qual (".$params{fasta}->id.") bioseq: $@"
        );
        return;
    }

    return $bioseq;
}

sub validate_bioseq {
    my ($self, $bioseq) = @_;

    unless ( $bioseq ) {
        confess $self->class." => No bioseq given to validate.";
    }

    unless ( ref($bioseq) ) {
        confess $self->class." => Bioseq given to validate is not an object: $bioseq";
    }

    unless ( $bioseq->seq =~ /^[ATGCNX]+$/i ) {
        confess sprintf(
            "%s => Bioseq (%s) has illegal characters in sequence:\n%s",
            $self->class,
            $bioseq->id,
            $bioseq->seq,
        );
    }

    return 1 unless $bioseq->can('qual');
    
    unless ( length($bioseq->seq) == scalar(@{$bioseq->qual}) ) {
        confess sprintf(
            '%s => Bioseq (%s) has unequal length for fasta and quality',
            $self->class,
            $bioseq->id,
        );
    }

    return 1;
}

sub validate_fasta_and_qual_bioseq {
    my ($self, $fasta, $qual) = @_;

    unless ( $qual ) {
        confess $self->class." => No qual bioseq given to validate.";
    }
    
    Genome::Utility::BioPerl->validate_bioseq($fasta);

    # validate length of seq v. qual
    unless ( length($fasta->seq) == scalar(@{$qual->qual}) ) {
        confess sprintf(
            '%s => Unequal length for fasta (%s) and quality (%s)',
            $self->class,
            $fasta->id,
            $qual->id,
        );
    }

    return 1;
}

1;

=pod

=head1 Disclaimer

Copyright (C) 2010 Genome Center at Washington University in St. Louis

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$
