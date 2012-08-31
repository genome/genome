package Genome::Model::Tools::Breakdancer::MergeFiles;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::Tools::Breakdancer::MergeFiles {
    is => 'Command',
    has => [
        input_files => {
            is => 'String',
            doc => 'Comma delimited list of breakdancer output files to be merged',
        },
        output_file => {
            is => 'FilePath',
            doc => 'Path to merged file',
        },
    ],
};

sub help_synopsis {
    return 'Merges several breakdancer output files into one file';
}

sub help_brief {
    return 'Merges several breakdancer output files into one file';
}

sub help_detail {
    return 'Merges several breakdancer output files into one file';
}

sub execute {
    my $self = shift;
    #my @files = grep { -e $_ } split(',', $self->input_files);
    #confess 'No files found that for merging!' unless @files;
    my @files = split ',', $self->input_files;

    for my $input_file (@files) {
        unless (-e $input_file) {
            $self->error_message("Input file for merge: $input_file not existing");
            die;
        }
    }

    if (-e $self->output_file) {
        $self->warning_message('Removing existing output file at ' . $self->output_file);
        unlink $self->output_file;
    }
    my $output_fh = IO::File->new($self->output_file, 'w');
    confess 'Could not get file handle for output file ' . $self->output_file unless $output_fh;

    # The header from each input file can be ignored, but we still need the column header
    my $header;
    my $input_fh;
    for my $file (@files) {
        $input_fh->close if $input_fh;
        $input_fh = IO::File->new($file, 'r');
        confess 'Could not get file handle for input file ' . $file unless $input_fh;

        while (my $line = $input_fh->getline) {
            if ($line =~ /^#/) {
                next if defined $header;
                next unless $line =~ /^#Chr1/i; #tigra validation output using CHR1
                $header = $line;
            }
            $output_fh->print($line);
        }
    }

    $output_fh->close;
    return 1;
}

        
1;

