package Genome::Model::Tools::BedTools::Sort;

use strict;
use warnings;

use Genome;
use File::Copy;

class Genome::Model::Tools::BedTools::Sort {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The input file to sort.',
        },
        output_file => {
            is => 'Text',
            doc => 'The output file to write sorted output.  If not defined the original file is overwritten with the sorted output.',
            is_optional => 1,
            is_output => 1,
        },
    ],
};

sub help_brief {
    "Sorts a feature file in various and useful ways.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools sort --input-file example.bed
EOS
}

sub help_detail {
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/.
EOS
}

sub execute {
    my $self = shift;

    my $input_file = $self->input_file;
    my $output_file = $self->output_file;
    my $copy_over_input = 0;
    unless (defined($output_file)) {
        $output_file = Genome::Sys->create_temp_file_path;
        $copy_over_input = 1;
        $self->output_file($input_file);
    }
    my $cmd = $self->bedtools_path .'/bin/sortBed -i  '. $input_file .' > '. $output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$input_file],
        output_files => [$output_file],
    );
    unless (copy($output_file,$input_file)) {
        die('Failed to copy sorted temp file '. $output_file .' over input file '. $input_file);
    }
    return 1;
}

1;
