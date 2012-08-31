package Genome::Model::Tools::Joinx::Union;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Joinx::Union {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_files => {
            is => 'Text',
            doc => 'Sorted bed "A"',
            shell_args_position => 1,
            is_many => 1,
        },
    ],
    has_optional_input => [
        output_file => {
            is => 'Text',
            doc => 'The output file (defaults to stdout)',
        },
    ],
};

sub help_brief {
    "Compute union of 2 bed files."
}

sub help_synopsis {
    my $self = shift;
    "gmt joinx union a.bed b.bed [--output-file=n.bed]"
}

sub execute {
    my $self = shift;
    my $output = "-";
    my $output_file = $self->output_file || '-';
    my @input_files = $self->input_files;
    print "INPUT FILES:\n\t" . join("\n\t", @input_files) . "\n";
    my %params = (
        unique => 1,
        input_files => \@input_files,
        output_file => $output_file,
    );

    my $cmd = Genome::Model::Tools::Joinx::Sort->create(
        \%params
    );
    return $cmd->execute();
}

1;
