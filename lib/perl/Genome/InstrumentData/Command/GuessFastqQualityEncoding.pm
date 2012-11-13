package Genome::InstrumentData::Command::GuessFastqQualityEncoding;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::GuessFastqQualityEncoding {
    is => 'Command::V2',
    has => [
        'input_file' => {   is => 'String',
                            shell_args_position => 1,
                            doc => 'Fastq file to read from' },
    ],
    has_optional => [
        'records' =>    {   is => 'Integer',
                            doc => 'Stop after processing this many records',
                            default_value => 10000 },
        'quiet'     => {    is => 'Boolean',
                            doc => 'Do not print anything to stdout',
                            default_value => 0 },
    ],
    has_output => [
        encoding => { is => 'ARRAY', doc => 'The best guess as to what the quality values are encoded in' },
    ],
    doc => 'Guess the quality string encoding for a fastq file',
};

sub help_detail {
    <<EOF
Given a FASTQ file, the command reads the quality values and attempts to
guess whether they are encoded in Sanger, Solexa, Illumina1.3, Illumina1.5
or Illumina1.8 format.  See http://en.wikipedia.org/wiki/FASTQ_format#Encoding
for discussion on the different encodings.

It will run until it has read in the number of records given in the 'records'
option, or until it receives an interrupt (Control-C).
EOF
}

# Class for records read from the file, and the data source
class Genome::FastqRecord {
    id_by => [
        'path'      => { is => 'String', column_name => '__FILE__' },
        'record_no' => { is => 'Integer', column_name => '$.' },
    ],
    has => ['seq_id','sequence', 'quality_string' ],
    data_source => {
        is => 'UR::DataSource::Filesystem',
        handle_class => 'FastQReader',   # FastQReader->getline() is at the end of this module
        path => '$path',
        columns => ['seq_id','sequence',undef, 'quality_string'],
        delimiter => "\n",
        record_separator => "\n",
    },
};

# Possible encoding formats.  Each sub returns false if the quality code is not
# a valid value
my %formats = (
    Sanger => sub {
        my($quals) = @_;
        for (@$quals) {
            if ($_ < 33 or $_ > 73) {
                return 0;
            }
        }
        return 1;
    },

    Solexa => sub {
        my($quals) = @_;
        for (@$quals) {
            if ($_ < 59 or $_ > 104) {
                return 0;
            }
        }
        return 1;
    },

    'Illumina1.3' => sub {
        my($quals) = @_;
        for (@$quals) {
            if ($_ < 64 or $_ > 104) {
                return 0;
            }
        }
        return 1;
    },

    'Illumina1.5' => sub {
        my($quals) = @_;
        for (@$quals) {
            if ($_ < 66 or $_ > 104) {
                return 0;
            }
        }
        return 1;
    },

    'Illumina1.8' => sub {
        my($quals) = @_;
        for (@$quals) {
            if ($_ < 33 or $_ > 74) {
                return 0;
            }
        }
        return 1;
    },
);

sub execute {
    my $self = shift;

    my $iter = Genome::FastqRecord->create_iterator(path => $self->input_file);
    my $count;
    my $max_read = $self->records();

    local $SIG{'INT'} = sub { $count = $max_read + 1; };

    while($count++ < $max_read and my $read = $iter->next()) {
        print "read $count\r" unless ($self->quiet);
        my $qual = $read->quality_string;
        my @raw_qual = map { ord } split(//,$qual);

        foreach my $format (keys %formats) {
            unless ($formats{$format}->(\@raw_qual)) {
                delete $formats{$format};
            }
        }
    }

    $self->encoding([ keys %formats ]);
    unless ($self->quiet) {
        print "\nPossible format(s): ", join(', ', keys %formats),"\n";
    }
        
    1;
}

package FastQReader;
use base 'IO::File';

sub getline {
    my $self = shift;

    local $/ = "\n";
    my $string = '';
    # Read 4 lines in a record
    for (1 .. 4) {
        no warnings 'uninitialized';
        $string .= <$self>;
    }
    return $string;
}

1;
