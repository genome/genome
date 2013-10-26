package Genome::Model::Tools::Mutect::MergeOutputFiles;

use Data::Dumper;
use Genome;
use Carp qw/confess/;
use strict;
use warnings;

class Genome::Model::Tools::Mutect::MergeOutputFiles {
    is => ['Command'],
    doc => "Merge together Mutect's native output files",
    has_input => [
        mutect_output_files => {
            is => 'Filename',
            is_many => 1,
            is_optional => 0,
            doc => 'mutect file names to merge together. Should be in order or you should sort afterwards.',
        },
    ],
    has_output => [
        merged_file => {
            is => 'String',
            doc => 'merged file',
        },
    ],
};

sub help_brief {
    return <<"EOS"
merge together mutect output files while respecting their headers
EOS
}

sub help_detail {
    return <<"EOS"
This tool is used to merge together output files from parallel runs of the main
mutect tool. You should pass these output files in order or sort after merging.
Headers from individual files are stripped out and a single header from the
first file is passed forwards as the header.
EOS
}

sub execute {
    my $self = shift;
    # the first two lines of every file should consist of a comment containing the mutect version
    # and also the header line

    my $output_fh = Genome::Sys->open_file_for_writing($self->merged_file);

    my $has_header = 0;
    for my $file ($self->mutect_output_files) {
        my $fh = Genome::Sys->open_file_for_reading($file);
        unless($has_header) {
            print $output_fh $fh->getline, $fh->getline;
            $has_header = 1;
        }
        else {
            #skip the header of this file as we already have one
            $fh->getline;
            $fh->getline;
        }
        while(my $line = $fh->getline) {
            print $output_fh $line;
        }
    }
    return 1;
}

1;
