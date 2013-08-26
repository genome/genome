package Genome::Model::Tools::Dindel::MakeDindelWindows;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dindel::MakeDindelWindows {
    is => 'Genome::Model::Tools::Dindel::Base',
    has_input => [
        input_dindel_file => {
            is => 'Path',
            doc => 'file of dindel formatted indels to examine. get this from getcigarindels or vcftodindel followed by realigncandidates',
        },
        num_windows_per_file =>  {
            is => 'Number',
            is_optional => 1,
        },
        output_directory => {
            is => 'Path',
        },
    ],
    has_output => [
        output_files => {
            is => 'Path',
            is_many => 1,
            is_calculated => 1,
            calculate =>  q{ glob("$output_prefix*") },
            calculate_from => ['output_prefix'],
        },
    ],
    has_transient => {
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

    my @cmd = (
        'python', $self->python_script('makeWindows'),
        '--inputVarFile', $self->input_dindel_file,
        '--windowFilePrefix', $self->output_prefix,
    );

    if (defined($self->num_windows_per_file)) {
        push @cmd, '--numWindowsPerFile', $self->num_windows_per_file;
    }

    return Genome::Sys->shellcmd_arrayref(
        cmd => \@cmd,
        input_files => [
            $self->input_dindel_file,
        ],
    );
}

1;
