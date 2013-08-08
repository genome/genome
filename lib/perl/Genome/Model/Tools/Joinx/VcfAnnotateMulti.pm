package Genome::Model::Tools::Joinx::VcfAnnotateMulti;

use Data::Dumper;
use Carp qw/confess/;
use Genome;

use strict;
use warnings;

class Genome::Model::Tools::Joinx::VcfAnnotateMulti {
    is => "Genome::Model::Tools::Joinx",
    has_input => [
        annotation_specs => {
            is => "Genome::Model::Tools::Joinx::VcfAnnotationSpec",
            is_many => 1
        },
        input_file => {
            is => "Text",
            doc => "The input file to annotate",
        },
        output_file => {
            is => "Text",
            is_optional => 1,
        },
        bgzip_output => {
            is => "Flag",
            default => 0,
            doc => "Compress output with bgzip (requires --output-file)",
        },
    ],
};

sub _is_hidden_in_docs { 1 }

sub build_command_line {
    my $self = shift;
    my @cmds = $self->_resolve_input;

    my $jx = $self->joinx_path;

    my @specs = $self->annotation_specs;
    confess $self->error_message("No annotation specs given") unless @specs;
    for my $spec (@specs) {
        my $afile = $spec->annotation_file;
        my @cmdline = qq/$jx vcf-annotate --input-file - --annotation-file $afile/;
        if (!$spec->identifiers) {
            push(@cmdline, "--no-identifiers");
        }

        if (!$spec->info) {
            push(@cmdline, "--no-info");
        }

        elsif (defined $spec->info_fields) {
            my $info_fields = "-I " . join(" -I ", @{$spec->info_fields});
            push(@cmdline, $info_fields);
        }

        push(@cmds, join(" ", @cmdline));
    }

    return @cmds;
}

sub _resolve_input {
    my $self = shift;

    my $input_file = $self->input_file;
    confess $self->error_message("Input file $input_file does not exist or is empty") unless -s $input_file;

    my $program = "cat";
    if (Genome::Sys->file_is_gzipped($input_file)) {
        $program = "bgzip -dc";
    }

    return "$program $input_file";
}

sub _resolve_output {
    my $self = shift;
    return "" unless $self->output_file;

    my $output_file = $self->output_file;

    if ($self->bgzip_output) {
        return " | bgzip -c > $output_file";
    }
    else {
        return " -o $output_file";
    }
}

sub _validate_args {
    my $self = shift;
    if ($self->bgzip_output && !$self->output_file) {
        confess $self->error_message("If bgzip_output is set, output_file must also be set, otherwise binary nonsense will spew forth.");
    }
}

sub execute {
    my $self = shift;

    $self->_validate_args;

    my @cmds = $self->build_command_line;
    my $pipe = join(" | ", @cmds) . $self->_resolve_output;

    my %params = (
        cmd => $pipe,
    );
    Genome::Sys->shellcmd(%params);


    return 1;
}

1;
