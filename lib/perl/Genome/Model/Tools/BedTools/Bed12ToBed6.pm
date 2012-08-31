package Genome::Model::Tools::BedTools::Bed12ToBed6;

use strict;
use warnings;

use Genome;
use File::Copy;

class Genome::Model::Tools::BedTools::Bed12ToBed6 {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_bed12_file => {
            is => 'Text',
            doc => 'The input BED12 file.',
        },
        output_bed6_file => {
            is => 'Text',
            doc => 'The output file to write BED6 format.',
            is_optional => 1,
            is_output => 1,
        },
        force_score => {
            is => 'Boolean',
            doc => 'Force the score to be the (1-based) block number from the BED12.',
            is_optional => 1,
            default_value => 0,
        },
    ],
};

sub help_brief {
    "Splits BED12 features into discrete BED6 features..",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools bed12-to-bed6 --input-bed12-file example.bed --output-bed6-file= output.bed
EOS
}

sub help_detail {
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/.
EOS
}

sub execute {
    my $self = shift;

    my $input_file = $self->input_bed12_file;
    my $output_file = $self->output_bed6_file;
    my $cmd = $self->bedtools_path .'/bin/bed12ToBed6 -i  '. $input_file;
    if ($self->force_score) {
        $cmd .= ' -n ';
    }
    $cmd .= ' > '. $output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$input_file],
        output_files => [$output_file],
    );
    return 1;
}

1;
