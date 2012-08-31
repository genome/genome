package Genome::Model::Tools::RepeatMasker::CompareTables;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::RepeatMasker::CompareTables {
    is => 'Genome::Model::Tools::RepeatMasker::TableI',
    has => [
        input_tables => {
            is => 'Text',
            doc => 'A space separated list of RepeatMasker style tables',
        },
        _total_count => {
            is => 'Integer',
            is_optional => 1,
        },
        _total_bp => {
            is => 'Integer',
            is_optional => 1,
        },
    ],
};

sub help_detail {
    'This tools is used to compare the output from multiple runs of RepeatMasker or the gmt bio-samtools repeat-content tool.  Each input table is normalized with the other input tables to generate one file that can be imported into Excel for comparsion.';
}

sub execute {
    my $self = shift;

    my $input_tables = $self->input_tables;
    my @input_tables;
    if (ref($input_tables)) {
        @input_tables = @{$input_tables};
    } else {
        @input_tables = split(' ', $input_tables);
    }
    my %samples;
    for my $table ( @input_tables ) {
        my $sample_name = basename($table);
        my $table_fh = IO::File->new($table,'r');
        unless ($table_fh) {
            die('Failed to open file to read '. $table);
        }
        my $family;
        while ( my $line = $table_fh->getline ) {
            chomp($line);
            my @entry = split("\t",$line);
            if (scalar(@entry) == 2) {
                if ($entry[0] eq 'sequences:') {
                    $samples{$sample_name}{sequences} = $entry[1];
                } elsif ($entry[0] eq 'total length:') {
                    $samples{$sample_name}{'base_pair'} = $entry[1];
                } elsif ($entry[0] eq 'aligned:') {
                    $samples{$sample_name}{'aligned'} = $entry[1];
                } elsif ($entry[0] eq 'repeat aligned:') {
                    $samples{$sample_name}{'repeat_aligned'} = $entry[1];
                }
                # new repeat family
            } elsif (scalar(@entry) == 4) {
                $family = $entry[0];
                $family =~ s/://;
                $samples{$sample_name}{$family}{base_pair} += $entry[2];
                $samples{$sample_name}{masked} += $entry[2];
            } elsif (scalar(@entry) == 5) {
                unless ($family) { die; }
                my $class = $entry[1];
                $class =~ s/://;
                $samples{$sample_name}{$family}{$class}{base_pair} += $entry[3];
            }
        }
        $table_fh->close;
    }
    $self->print_samples_summary_from_hash_ref(\%samples);
    return 1;
}

1;
