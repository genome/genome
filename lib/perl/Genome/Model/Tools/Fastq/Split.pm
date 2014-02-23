package Genome::Model::Tools::Fastq::Split;

use strict;
use warnings;

require Cwd;
require File::Basename;
require File::Copy;
require Genome;

class Genome::Model::Tools::Fastq::Split {
    is  => 'Command',
    has => [
        fastq_file => {
            type     => 'Text',
            doc      => 'FASTQ file that contains both sequences and quality values',
            is_input => 1,
        },
        split_size => {
            type     => 'Integer',
            doc      => 'Number of fastq sequences for each output file',
            is_input => 1,
        },
    ],
    has_optional => [
        output_directory => {
            is_input => 1,
        },
        split_files => {
            is        => 'Text',
            is_many   => 1,
            is_output => 1,
        },
        show_list => {
            is  => 'Boolean',
            doc => 'show list of split files',
        },
    ],
};

sub help_brief {
    return 'Split fastq into multiple fastqs each with --split-size sequences.'
}


sub help_detail {
    return 'Split fastq into multiple fastqs each with --split-size sequences.'
}

sub execute {
    my $self = shift;

    $self->dump_status_messages($self->show_list);

    my @suffix = qw/\.txt \.fastq \.fq/;
    my ($fastq_basename, $fastq_dirname, $fastq_suffix) = File::Basename::fileparse($self->fastq_file, @suffix);
    unless ($fastq_basename && $fastq_dirname && $fastq_suffix) {
        die('Failed to parse fastq file name ' . $self->fastq_file);
    }
    unless (Genome::Sys->validate_directory_for_read_write_access($fastq_dirname)) {
        $self->error_message('Failed to validate directory ' . $fastq_dirname . " for read/write access:  $!");
        die($self->error_message);
    }

    my $cwd = Cwd::getcwd();
    my $tmp_dir = $self->output_directory || Genome::Sys->base_temp_directory;
    chdir($tmp_dir);
    my $cmd = 'split -l ' . $self->split_size * 4 . ' -a 5 -d ' . $self->fastq_file . ' ' . $fastq_basename . '-';
    Genome::Sys->shellcmd(
        cmd         => $cmd,
        input_files => [$self->fastq_file],
    );
    chdir($cwd);
    #If more than one lane processed in the same output_directory, this becomes a problem
    my @tmp_fastqs = grep { $_ !~ /\.$fastq_suffix$/ } grep { /$fastq_basename-\d+$/ } glob($tmp_dir . '/' . $fastq_basename . '*');

    # User should provide a directory as input, then we can keep output fastqs on tmp
    # and distribute bfqs in a downstream process
    # However, by default write fastqs to the source fastq file dir
    my $output_dir = $self->output_directory || $fastq_dirname;

    for my $tmp_fastq (@tmp_fastqs){
        my ($tmp_fastq_basename,$tmp_fastq_dirname) = File::Basename::fileparse($tmp_fastq);
        my $fastq_file = $output_dir . '/' . $tmp_fastq_basename . $fastq_suffix;
        unless (File::Copy::move($tmp_fastq,$fastq_file,) ) {
            die('Failed to move file ' . $tmp_fastq . ' to ' . $fastq_file . ":  $!");
        }
        $self->add_split_file($fastq_file);
    }

    if ($self->show_list) {
        for my $file ($self->split_files) {
            $self->debug_message($file);
        }
    }

    return 1;
}

1;

