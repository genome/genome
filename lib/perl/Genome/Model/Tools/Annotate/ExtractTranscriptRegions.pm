package Genome::Model::Tools::Annotate::ExtractTranscriptRegions;

use Genome;
use strict;
use warnings;

class Genome::Model::Tools::Annotate::ExtractTranscriptRegions {
    is => "Command::V2",
    doc => "Generate a bed file describing transcript regions",
    has_input => [
        transcripts_file => {
            is => "Text",
            doc => "GTF format file describing transcripts",
        },
        output_file => {
            is => "Text",
            doc => "Output bed file",
        },
        feature_type => {
            is => "Text",
            default_value => "CDS",
            doc => "The type of features to pay attention to",
        }
    ],
};

sub _parse_transcript_attr {
    my $attr_str = shift;

    # turn key1 "value1"; key2 "value2"; ... 
    # into hash { key1 => value1, key2 => value2, ... }
    my %attrs = map { /^ *([^ ]*) "([^"]*)"/; ($1, $2)} split(";", $attr_str);
    return %attrs;
}

sub execute {
    my $self = shift;
    my $reader = Genome::Utility::IO::GffReader->create(
        input => $self->transcripts_file,
    );
    my $writer = Genome::Sys->open_file_for_writing($self->output_file);
    my $feature_type = $self->feature_type;

    my %transcripts;
    my $line_num = 0;
    while (my $entry = $reader->next_with_attributes_hash_ref) {
        next if $entry->{type} ne $feature_type;

        my $attrs = delete $entry->{attributes_hash_ref};
        my $tid = $attrs->{transcript_id};
        next unless defined $tid;

        if (!defined $transcripts{$tid}{chr}) {
            $transcripts{$tid}{chr} = $entry->{chr};
        }
        elsif ($transcripts{$tid}{chr} ne $entry->{chr}) {
            die $self->error_message("Sequence name mismatch at line $line_num");
        }
        push(@{$transcripts{$tid}{regions}}, [$entry->{start}, $entry->{end}]);
    }

    for my $tid (keys %transcripts) {
        $writer->print("# $tid\n");
        my $chr = $transcripts{$tid}{chr};
        for my $region (sort {$a->[0] <=> $b->[0]} @{$transcripts{$tid}{regions}}) {
            $writer->print(join("\t", $chr, @$region) . "\n");
        }
    }
    return 1;
}
