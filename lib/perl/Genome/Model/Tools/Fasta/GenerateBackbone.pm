package Genome::Model::Tools::Fasta::GenerateBackbone;

use strict;
use warnings;

use Genome;
use FASTAParse;

class Genome::Model::Tools::Fasta::GenerateBackbone {
    is => 'Genome::Model::Tools::Fasta',
    has => [
            output_file => {
                            is => 'Text',
                            doc => 'The file path name to write backbone',
                        },
        ],
};

sub help_brief {
    "Convert a transcriptome reference into a backbone file.",
}

sub help_detail {
    return <<EOS
        RefCov requires a backbone file, a tab seperated value file, that
contains the name and length of all the genes within the transcriptome reference.
EOS
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }

    unless (Genome::Sys->validate_file_for_reading($self->fasta_file)) {
        $self->error_message('Failed to validate file '. $self->fasta_file .' for reading.');
        return;
    }
    unless (Genome::Sys->validate_file_for_writing($self->output_file)) {
        $self->error_message('Failed to validate file '. $self->output_file .' for writing.');
        return;
    }
    return $self;
}

sub execute {
    my $self = shift;

    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    unless ($output_fh) {
        $self->error_message('Failed to open output file '. $self->output_file);
        return;
    }

    my $fasta_reader = Genome::Sys->open_file_for_reading($self->fasta_file);
    unless ($fasta_reader) {
        $self->error_message('Failed to open fasta file '. $self->fasta_file);
        return;
    }


    local $/ = "\n>";
    while (<$fasta_reader>) {
        if ($_) {
            chomp;
            if ($_ =~ /^>/) { $_ =~ s/\>//g }
            my $myFASTA = FASTAParse->new();
            $myFASTA->load_FASTA( fasta => '>' . $_ );
            my $seqlen = length( $myFASTA->sequence() );
            print $output_fh join("\t", $myFASTA->id(), '1', $seqlen, $myFASTA->id()) . "\n";
        }
    }
    $fasta_reader->close;
    $output_fh->close;

    return 1;
}
