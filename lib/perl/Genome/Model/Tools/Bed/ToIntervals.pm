package Genome::Model::Tools::Bed::ToIntervals;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::ToIntervals {
    is => ['Command'],
    has_input => [
        bed_file => {
            is => 'Text',
            doc => 'The path to a BED format file.',
        },
        seqdict_file => {
            is => 'Text',
            doc => 'The sequence dictionary for a reference genome of which the annotation is based on.',
        },
        interval_file => {
            is => 'Text',
            doc => 'The output Interval format file.',
            is_output => 1
        },
    ],
};

sub help_brief {
    "Convert a BED file into a Interval file.",
}

sub help_synopsis {
my $self = shift;
    return <<"EOS"
gmt bed to-interval...
EOS

};

sub help_detail {
    'This command will take a BED format file and convert it to a Interval format file.
'
}

sub execute {
    my $self = shift;

    my $bed_fh = Genome::Sys->open_file_for_reading($self->bed_file);
    unless ($bed_fh) {
        die('Failed to load BED reader: '. $self->bed_file);
    }

    my @interval_headers = qw/chr start end strand name/;
    my $tmp_file = Genome::Sys->create_temp_file_path();
    my $interval_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $tmp_file,
        headers => \@interval_headers,
        separator => "\t",
        print_headers => 0,
        ignore_extra_columns => 1,
    );
    while (my $line = $bed_fh->getline) {
        chomp($line);
        if ($line =~ /^track/) { next; }
        my ($chr,$start,$end,$name,$score,$strand) = split("\t", $line);
        my %data = (
            name => $name || ($chr .':'. $start .'-'. $end),
            chr => $chr,
            start => $start,
            end => $end,
            strand => $strand || '+',
        );
        $interval_writer->write_one(\%data);
    }
    $interval_writer->output->close;

    my $cmd = 'cat '. $self->seqdict_file .' '. $tmp_file .' > '. $self->interval_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->seqdict_file,$tmp_file],
        output_files => [$self->interval_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
