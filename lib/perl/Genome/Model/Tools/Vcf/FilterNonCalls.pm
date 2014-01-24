package Genome::Model::Tools::Vcf::FilterNonCalls;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Reader;
#use Genome::File::Vcf::Entry;

class Genome::Model::Tools::Vcf::FilterNonCalls {
    doc => 'Remove all lines from a VCF file where the alt is "." or "N" with no other alts present.',
    is => 'Command',
    has_input => [
        input_file=> {
            is => 'Text',
            doc => "Input VCF file",
        },
        output_file=> {
            is => 'Text',
            is_output => 1,
            doc => "Output VCF file",
        },
    ],
};

sub help_synopsis {
    <<'HELP';
Remove all lines from a VCF file where the alt is "." or "N" with no other alts present.
HELP
}

sub help_detail {
    <<'HELP';
Remove all lines from a VCF file where the alt is "." or "N" with no other alts present.
HELP
}

sub execute {
    my $self = shift;

    my $reader = $self->_initialize_reader($self->input_file);
    my $output_fh = $self->_initialize_output($reader, $self->output_file);

    while (my $entry = $reader->next) {
        if ($self->_entry_is_valid($entry)) {
            $output_fh->print($entry->to_string . "\n");
        }
    }
    $output_fh->close;

    return 1;
}

# The vcf entry is valid if there are any alt values besides "." or "N"
sub _entry_is_valid {
    my ($self, $entry) = @_;

    my $alts = $entry->{alternate_alleles};
    for my $alt (@$alts) {
        unless ( ($alt eq ".") || (uc($alt) eq "N") ) {
            return 1;
        }
    }
    return;
}

sub _initialize_reader {
    my ($self, $input_file) = @_;
    Genome::Sys->validate_file_for_reading($input_file);
    my $reader = Genome::File::Vcf::Reader->new($input_file);
    unless ($reader) {
        die $self->error_message("Failed to create a Genome::File::Vcf::Reader object with file $input_file");
    }
    return $reader;
}

# Open an output filehandle and copy the header from the input vcf over
sub _initialize_output {
    my ($self, $reader, $output_file) = @_;
    my $output_fh = Genome::Sys->open_gzip_file_for_writing($output_file);
    my $header = $reader->header->to_string;
    unless ($header) {
        die $self->error_message("Header not returned from the vcf reader");
    }
    $output_fh->print($header."\n");
    return $output_fh;
}

1;

