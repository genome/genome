package Genome::Model::Tools::Fasta::Concat;

use Genome;
use IO::File;
use Carp qw/confess/;

class Genome::Model::Tools::Fasta::Concat {
    is => 'Genome::Command::Base',
    has_input => [
        input_files => {
            is => 'Text',
            is_many => 1,
            doc => 'The input fasta files to concatenate',
        },
        output_file => {
            is => 'Text',
            doc => 'The output file to write', 
        }
    ]
};

sub execute {
    my $self = shift;

    my @missing_input_files = grep { ! -s $_ } $self->input_files;
    if (@missing_input_files) {
        confess "Missing input files: \n\t" . join("\n\t", @missing_input_files);
    }

    my %sequences;
    my $output_file = $self->output_file;
    my $out_fh = new IO::File(">$output_file") or
        confess "Failed to open output fasta file: $output_file";

    for my $input ($self->input_files) {
        my $fh = new IO::File("<$input") or
            confess "Failed to open input fasta file: $input";

        while (<$fh>) {
            chomp;
            next if $_ eq "";
            if (/^[>;](.*)/) {
                my $sequence_name = $1;
                if (exists $sequences{$sequence_name}) {
                    confess "Duplicate sequence '$sequence_name' found in $input and $sequences{$sequence_name}";
                }
                $sequences{$sequence_name} = $input;
                $out_fh->print(">$sequence_name\n");
            } else {
                $out_fh->print("$_\n");
            }
        }

        $fh->close();
    }

    return 1;
}

1;
