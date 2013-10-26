package Genome::Model::Tools::Dindel::MakeDindelWindows;

use strict;
use warnings;

use Genome;
use File::Spec;

class Genome::Model::Tools::Dindel::MakeDindelWindows {
    is => 'Genome::Model::Tools::Dindel::Base',
    has_input => [
        input_dindel_file => {
            is => 'Path',
            doc => 'file of dindel formatted indels to examine. get this from GetCigarIndels or VcfToDindel followed by RealignCandidates',
        },
    ],
    has_calculated_output => [
        window_file => {
            is => 'Path',
            is_calculated => 1,
            calculate =>  q{ $output_prefix . ".txt" },
            calculate_from => ['output_prefix'],
        },
    ],
    has_optional_transient => {
        output_prefix => {
            is_calculated => 1,
            calculate =>  q{ File::Spec->join($output_directory, "dindel_window_file") },
            calculate_from => ['output_directory'],
        },
    },
};

sub help_brief {
    'make window files for dindel'
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
EOS
}


sub execute {
    my $self = shift;

    $self->create_output_directory();

    my @cmd = (
        'python', $self->python_script('makeWindows'),
        '--inputVarFile', $self->input_dindel_file,
        '--windowFilePrefix', $self->output_prefix,
    );

    return $self->run_with_single_output(@cmd);
}

sub run_with_single_output {
    my ($self, @cmd) = @_;

    push @cmd, '--numWindowsPerFile', 9_999_999; # author of .py hard-coded 10_000_000 as starting value.
    my $result = $self->run(@cmd);

    # ensure there is only one output.
    my @output_files = glob($self->output_prefix . "*");
    Genome::Sys->concatenate_files(\@output_files, $self->window_file);

    return $result;
}

sub run {
    my ($self, @cmd) = @_;
    return $self->shellcmd_arrayref(
        cmd => \@cmd,
        input_files => [
            $self->input_dindel_file,
        ],
    );
}

1;
