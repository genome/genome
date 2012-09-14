package Genome::Model::Tools::Sam::FilterPairedReads;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::Model::Tools::Sam::FilterPairedReads {
    is => 'Genome::Model::Tools::Sam',
    has => [
        input_file => {
            is => 'Path',
            doc => 'Name sorted input BAM file',
        },
        output_file => {
            is => 'Text',
            doc => 'Output BAM File',
        },
    ],
};

sub help_detail {
    return <<HELP
Removes paired reads in a namesorted mixed paired and unpaired BAM file
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

        $DB::single = 1;
        my ($prev_name) = split("\t", $prev_line);
        my ($name) = split("\t", $line);

        if($prev_line eq ""){
            $prev_line = $line;
            next;
        }

        if($prev_name eq $name){
            $prev_line = "";
            next;
        }
        print $output_fh $prev_line; 
        $prev_line = $line;
    }

    $input_fh->close;
    $output_fh->close;
    return 1;
}

1;
