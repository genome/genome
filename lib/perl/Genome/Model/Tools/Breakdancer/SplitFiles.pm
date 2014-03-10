package Genome::Model::Tools::Breakdancer::SplitFiles;

use strict;
use warnings;
use Genome;
use File::Basename 'dirname';
use Carp 'confess';

class Genome::Model::Tools::Breakdancer::SplitFiles {
    is => 'Command',
    has => [
        input_file => {
            is => 'FilePath',
            doc => 'Breakdancer file to be split up',
        },
    ],
    has_optional => [
        output_directory => {
            is => 'DirectoryPath',
            doc => 'Directory into which split files should go, defaults to same directory as input file',
        },
        output_file_template => {
            is => 'String',
            doc => 'Naming scheme that output file should follow, CHR is replaced with chromosome name',
            default => 'breakdancer_CHR',
        },
        output_files => {
            is => 'ARRAY',
            is_transient => 1,
            doc => 'Array containing paths of all created output files',
        },
        split_column => {
            is => 'Number',
            default => 'chr2',
            valid_values => ['chr1', 'chr2'],
            doc => 'Chromosome column to split on',
        },
        create_other => {
            is => "Boolean",
            default => 1,
            doc => 'Whether or not to restrict to canonical chromosomes plus an "other" chromosome',
        },
    ],
};

my @FULL_CHR_LIST = (1..22, 'X', 'Y', 'MT');

sub help_synopsis {
    return 'Splits up a breakdancer output file by chromosome';
}

sub help_brief {
    return 'Splits up a breakdancer output file by chromosome';
}

sub help_detail {
    return 'Splits up a breakdancer output file by chromosome';
}

sub execute {
    my $self = shift;
    confess 'No file at ' . $self->input_file unless -e $self->input_file;

    if (defined $self->output_directory) {
        my $dir = Genome::Sys->create_directory($self->output_directory);
        confess 'Could not find or create output directory ' . $self->output_directory unless defined $dir;
    }
    else {
        my $dir = dirname($self->input_file);
        confess "Could not create output directory $dir!" unless Genome::Sys->create_directory($dir);
        $self->output_directory($dir);
    }

    unless ($self->output_file_template =~ /CHR/) {
        $self->warning_message("Given template without CHR, appending it to the end");
        $self->output_file_template($self->output_file_template . 'CHR');
    }

    $self->debug_message("Split files being written to " . $self->output_directory);

    my $split_index = $self->split_column eq 'chr1' ? 0 : 3;
    my $input_fh = IO::File->new($self->input_file, 'r');
    my %output_handles;
    my $output_fh;
    my $chrom;
    my $header;
    my @files;
    while (my $line = $input_fh->getline) {
        if ($line =~ /^#/) {
            next if defined $header;
            if ($line =~ /^#Chr1/i) {
                $header = $line;
            }
            next;
        }

        my @fields = split("\t", $line);
        my $chr = $fields[$split_index];


        if($self->create_other) {
            unless (grep { $_ eq $chr } @FULL_CHR_LIST) {
                $chr = 'other';
            }
        }

        unless (defined $chrom and $chrom eq $chr) {
            $chrom = $chr;
            if (exists $output_handles{$chrom}) {
                $output_fh = $output_handles{$chrom};
            }
            else {
                my $file_name = $self->output_directory . '/' . $self->output_file_template;
                $file_name =~ s/CHR/$chrom/;
                unlink $file_name if -e $file_name;
                $output_fh = IO::File->new($file_name, 'w');
                confess "Could not get file handled for output file $output_fh!" unless $output_fh;
                $output_handles{$chrom} = $output_fh;
                $output_fh->print($header);
                push @files, $file_name;
                $self->debug_message("Created output file $file_name");
            }
        }

        $output_fh->print($line);
    }

    $input_fh->close;
    map { $output_handles{$_}->close } keys %output_handles;
    $self->output_files(\@files);
    return 1;
}

1;

