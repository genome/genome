package Genome::Model::Tools::Sam::FilterUnpairedReads;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::Model::Tools::Sam::FilterUnpairedReads {
    is => 'Genome::Model::Tools::Sam',
    has => [
        input_file => {
            is => 'Path',
            doc => 'Input BAM file',
        },
        output_file => {
            is => 'Text',
            doc => 'Output BAM File',
        },
    ],
};

sub help_detail {
    return <<HELP
Removes unpaired reads in a paired BAM file
HELP
}

sub execute {
    my $self = shift;
    my $input_path = $self->input_file;
    my $output_path = $self->output_file;
    my $input_fh = IO::File->new("samtools view -h $input_path | ");
    my $output_fh = IO::File->new("| samtools view -S -b -o $output_path - ");

    my $prev_line = "";
    while (my $line = <$input_fh>){
        if ($line =~ /^@/){
            print $output_fh $line;
            next;
        }
        my ($name) = split("\t", $line);
        my ($prev_name) = split("\t", $prev_line);
        unless ($prev_name and $prev_name eq $name){
            $prev_line = $line;
            next;
        }
        print $output_fh $prev_line, $line;
        $prev_line = '';
    }

    $input_fh->close;
    $output_fh->close;
    return 1;
}
