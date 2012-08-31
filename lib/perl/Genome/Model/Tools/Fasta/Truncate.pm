package Genome::Model::Tools::Fasta::Truncate;

use strict;
use warnings;

use Genome;

use FASTAParse;

class Genome::Model::Tools::Fasta::Truncate {
    is => 'Command',
    has => [
        input_fasta_file => { is => 'Text', },
        output_fasta_file => { is => 'Text', },
        start => { is => 'Integer', default_value => 1 },
        length => {is => 'Integer', default_value => 50 },
        #TODO: Add report of sequneces removed ie. zero length and less than $self->length
    ],
};

sub execute {
    my $self = shift;

    my $output_fh = Genome::Sys->open_file_for_writing($self->output_fasta_file);
    unless ($output_fh) {
        $self->error_message('Failed to open output file '. $self->output_fasta_file);
        die($self->error_message);
    }
    my $fasta_reader = Genome::Sys->open_file_for_reading($self->input_fasta_file);
    unless ($fasta_reader) {
        $self->error_message('Failed to open fasta file '. $self->input_fasta_file);
        die($self->error_message);
    }
    local $/ = "\n>";
    while (<$fasta_reader>) {
        if ($_) {
            chomp;
            if ($_ =~ /^>/) { $_ =~ s/\>//g }
            my $myFASTA = FASTAParse->new();
            $myFASTA->load_FASTA( fasta => '>' . $_ );
            my $truncated_sequence = substr( $myFASTA->sequence(), $self->start - 1, $self->length );
            if (length( $truncated_sequence ) == 50) {
                print $output_fh '>'. $myFASTA->id() ."\n". $truncated_sequence ."\n";
            }
        }
    }
    $fasta_reader->close;
    $output_fh->close;
    return 1;
}

1;
