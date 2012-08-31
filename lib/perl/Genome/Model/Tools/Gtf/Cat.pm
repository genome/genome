package Genome::Model::Tools::Gtf::Cat;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gtf::Cat {
    is => ['Command',],
    has_input => [
        input_files => { },
        output_file => {
            is_output => 1,
        },
        remove_originals => {
            default_value => 1,
            is_optional => 1,
        }
    ],
};

sub execute {
    my $self = shift;
    my @input_files;
    if (ref($self->input_files) eq 'ARRAY') {
        @input_files = @{$self->input_files};
    } else {
        @input_files = split(',',$self->input_files);
    }
    my @files_with_size;
    for my $file (@input_files) {
        if (-e $file) {
            if (-s $file) {
                push @files_with_size, $file;
            }
        } else {
            die('Missing file: '. $file);
        }
    }
    my $cmd = 'cat '. join(' ',@files_with_size) .' > '. $self->output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => \@files_with_size,
        output_files => [$self->output_file],
    );
    if ($self->remove_originals) {
        for my $file (@input_files) {
            unlink($file) || warn 'Failed to remove GTF file: '. $file;
        }
    }
    return 1;
}

1;
