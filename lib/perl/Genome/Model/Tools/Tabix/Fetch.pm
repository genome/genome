package Genome::Model::Tools::Tabix::Fetch;

use strict;
use warnings;

use Genome;
use List::Util qw/min/;

# Tabix can't read regions from a file, so we break up regions into
# batches to avoid exessively long command lines
my $REGION_BATCH_SIZE=1000;

class Genome::Model::Tools::Tabix::Fetch {
    is => 'Genome::Model::Tools::Tabix',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'Tabix indexed file to query',
            shell_args_position => 1,
        },
        regions => {
            is => 'Text',
            doc => 'A list of regions to fetch in chr:start-end format (1 based)',
            is_many => 1,
        },
        print_header => {
            is => 'Boolean',
            doc => 'Include file header in output',
            default_value => 0,
        },
    ],
    has_optional_input => [
        output_file => {
            is => 'Text',
            doc => 'The output file (defaults to stdout)',
        },
    ],
};

sub help_brief {
    return "Fetch one or more regions from a tabix indexed file.";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt tabix fetch in.vcf.gz --regions 1:180000-190000 --output-file=out.vcf --print-header
EOS
}

sub execute {
    my $self = shift;
    my $tabix = $self->tabix_path;
    my $input = $self->input_file;
    my @regions = $self->regions;
    my $print_header = $self->print_header ? '-h' : '';
    my $output = $self->output_file | '';

    if ($output) {
        # Just to ensure the file can be written to
        my $fh = Genome::Sys->open_file_for_writing($output);
        $fh->close;
        $output = ">> $output";
    }

    while (@regions) {
        my $len = min(scalar(@regions), $REGION_BATCH_SIZE);
        my $batch = join(" ", splice(@regions, 0, $len));
        my $cmd = "$tabix $print_header $input $batch $output";
        Genome::Sys->shellcmd(cmd => $cmd);
    }

    return 1;
}

1;
