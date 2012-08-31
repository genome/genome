package Genome::Model::Tools::Joinx::UnionBedPositions;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::Tools::Joinx::UnionBedPositions {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_files => {
            is => 'Text',
            is_many => 1,
            doc => 'List of bed files to sort',
            shell_args_position => 1,
        },
    ],
    has_optional_input => [
        unique => {
            is => 'Boolean',
            value => '1',
            doc => 'Print only unique entries (which are compared up to and including alleles)',
        },
        merge_only => {
            is => 'Boolean',
            value => '0',
            doc => 'If set, then the pre-sorted input files just merged',
        },
        output_file => {
            is => 'Text',
            is_output => 1,
            is_input => 1,
            doc => 'The output file (defaults to stdout)',
        },
    ],
};

sub help_brief {
    "Sorts one or more bed files."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt joinx sort a.bed [b.bed ...] --output-file sorted.bed
EOS
}

sub flags {
    my $self = shift;
    my @flags = ('--stable');
    push(@flags, "--merge-only") if $self->merge_only;
    push(@flags, "--unique") if $self->unique;
    return @flags;
}

sub execute {
    my $self = shift;
    my $output = "-";
    $output = $self->output_file if (defined $self->output_file);
    # Grep out empty files
    my @inputs = grep { -s $_ } $self->input_files;

    # If all input files are empty, make sure the output file at least exists
    unless (@inputs) {
        if (defined $self->output_file) {
            unless (system("touch $output") == 0) {
                die $self->error_message("Failed to touch $output");
            }
        }
        return 1;
    }

    my $flags = join(" ", $self->flags);
    my $cmd = $self->joinx_path . " sort $flags " .
        join(" ", map { "<(zcat $_)"} @inputs) .
        " | bgzip -c > $output";
    $cmd = "bash -c \"".$cmd."\"";

    my %params = (
        cmd => $cmd,
        # Sometimes these files come in empty in pipelines. We don't want this to choke when this happens, so don't check input file sizes.
        #input_files => \@inputs,
        allow_zero_size_output_files=>1,
    );
    $params{output_files} = [$output] if $output ne "-";
    Genome::Sys->shellcmd(%params);

    return 1;
}

1;
