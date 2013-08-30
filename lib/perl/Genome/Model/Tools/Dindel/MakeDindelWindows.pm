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
            doc => 'The number of windows dindel will put in each window_file.  If < 1, there will be only one window_file.',
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

    if (defined($self->num_windows_per_file)) {
        if ($self->num_windows_per_file < 1) {
            return $self->run_with_single_output(@cmd);
        } else {
            push @cmd, '--numWindowsPerFile', $self->num_windows_per_file;
            return $self->run(@cmd);
        }
    } else {
        return $self->run(@cmd);
    }
}

sub create_output_directory {
    my $self = shift;
    if (! -d $self->output_directory) {
        Genome::Sys->create_directory($self->output_directory);
    }
    if (! -d $self->output_directory) {
        die "Couldn't create output_directory " . $self->output_directory;
    }
}

sub run_with_single_output {
    my ($self, @cmd) = @_;

    push @cmd, '--numWindowsPerFile', 9_000_000;
    my $result = $self->run(@cmd);

    my @output_files = $self->output_files;
    if (scalar(@output_files) > 1) {
        die sprintf("Found more than one output_file (%s) when only one should have been produced.",
            scalar(@output_files));
    } else {
        return $result;
    }
}

sub run {
    my ($self, @cmd) = @_;
    return Genome::Sys->shellcmd_arrayref(
        cmd => \@cmd,
        input_files => [
            $self->input_dindel_file,
        ],
    );
}

1;
