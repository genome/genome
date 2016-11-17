package Genome::Model::Tools::Annotate::Sv::RepeatMasker;

use strict;
use warnings;
use Genome;
use List::Util qw(max min);

class Genome::Model::Tools::Annotate::Sv::RepeatMasker {
    is => 'Genome::Model::Tools::Annotate::Sv::Base',
    doc => "Annotate SVs whose breakpoints fall in repeat-masked regions",
    has_input => [
        annotation_file => {
            is => 'String',
            doc => 'Path to annotation file',
            example_values => ["repeat_masker.tsv"],
        },
        length_to_repeat => {
            type => 'Integer',
            doc => 'Look for repeat annotation within +- flanking bp of breakpoint',
            default_value => 200,
        },
        overlap_repeat_size => {
            type => 'Integer',
            doc => 'Overlapped size of +- flanking breakpoint and masked repeat',
            default_value => 5,
        },
        masked_repeat_size => {
            type => 'Integer',
            doc => 'size of masked repeat used for overlap_repeat_size',
            default_value => 100,
        },
    ],
};

sub help_detail {
    return "Determine whether the SV breakpoints fall into any repeat-masked regions"
}

sub process_breakpoint_list {
    my $self = shift;
    my $breakpoints_list = shift;

    my %output;
    my $annotation = $self->read_ucsc_annotation($self->annotation_file);

    $self->annotate_two_breakpoint_overlaps($breakpoints_list, $annotation, "repeatmasker");

    foreach my $chr (keys %{$breakpoints_list}) {
        foreach my $item (@{$breakpoints_list->{$chr}}) {
            my $key = $self->get_key_from_item($item);
            $output{$key} = $item->{repeatmasker};
            if (defined $item->{repeatmasker}) {
                $output{$key} = [$item->{repeatmasker}];
            }
            else {
                $output{$key} = ["N/A"];
            }
        }
    }
    return \%output;
}

sub column_names {
    return ("repeatmasker");
}

sub annotate_two_breakpoint_overlaps {
    my ($self, $breakpoint_list, $annotations, $tag) = @_;

    foreach my $chr (keys %$breakpoint_list) {
        my $chr_breakpoint_list = $breakpoint_list->{$chr};
        foreach my $breakpoint (@$chr_breakpoint_list) {
            my $chr1 = $breakpoint->{chrA};
            my $pos1 = $breakpoint->{bpA};
            my $chr2 = $breakpoint->{chrB};
            my $pos2 = $breakpoint->{bpB};
            my ($length1, $class1) = $self->annotate_breakpoint_overlap($chr1, $pos1, $annotations);
            my ($length2, $class2) = $self->annotate_breakpoint_overlap($chr2, $pos2, $annotations);

            my $repeat_annot = "-";
            my $rm_size = $self->masked_repeat_size;

            if ($length1 > $rm_size or $length2 > $rm_size) {
                $repeat_annot = sprintf "Repeat:%s-%s", $class1 || 'NA', $class2 || 'NA';
            }
            $breakpoint->{$tag} = $repeat_annot;
        }
    }
    return 1;
}

sub annotate_breakpoint_overlap {
    #given one breakpoint, report which intervals in a set of annotations overlap that breakpoint
    my ($self, $chr, $pos, $annotations) = @_;

    my $rp_length = $self->length_to_repeat;
    my $rp_overlap = $self->overlap_repeat_size;
    
    my $start = $pos - $rp_length;
    my $stop = $pos + $rp_length;
    
    #sort annotation intervals by start position
    my $chr_annotations = $annotations->{$chr};
    my %ends;
    foreach my $end (keys %{$chr_annotations}) {
        foreach my $start (keys %{$chr_annotations->{$end}}) {
            $ends{$start} = $end;
        }
    }
    my @sorted_by_start = sort {$a <=> $b} (keys %ends);
    #find all intervals that start before start and stop after start
    my %interval_records;
    while ($sorted_by_start[0] <= $start) {
        if ($ends{$sorted_by_start[0]} < $start) {
            shift @sorted_by_start;
            next;
        }
    #for each of these intervals, deterimine the "tightest" combination of endpoints
        my $start_last = max($sorted_by_start[0], $start);
        my $stop_last = min($ends{$sorted_by_start[0]}, $stop);
    #ensure that the tight interval is within rp_overlap of overlapping the actual position.
        if ($start_last - $rp_overlap <= $pos and $pos <= $stop_last + $rp_overlap) {
            my $overlap_length = $stop_last-$start_last+1;
            my $name = $chr_annotations->{$ends{$sorted_by_start[0]}}->{$sorted_by_start[0]}->[0]->{name};
            if (defined($interval_records{$name}) and $interval_records{$name} > $overlap_length) {
                next;
            }
            #if so, save that as a possible interval, including the name of the annotation and the length of the interval
            else {
                $interval_records{$name} = $overlap_length;
            }
        }
        shift @sorted_by_start;
    }
    #sort annotation intervals by stop position
    #find all intervals that start before stop and stop after stop
    my @sorted_by_end = sort {$b <=> $a} (keys %{$chr_annotations});
    while ($sorted_by_end[0] >= $stop) {
        my @per_end_starts = keys %{$chr_annotations->{$sorted_by_end[0]}};
        foreach my $per_end_start (@per_end_starts) {
            if ($per_end_start > $stop) {
                shift @sorted_by_end;
                next;
            }
            #for each of these intervals, deterimine the "tightest" combination of endpoints
            my $start_last = max($per_end_start, $start);
            my $stop_last = min($sorted_by_end[0], $stop);
            #ensure that the tight interval is within rp_overlap of overlapping the actual position.
            if ($start_last - $rp_overlap <= $pos and $pos <= $stop_last + $rp_overlap) {
                my $overlap_length = $stop_last-$start_last+1;
                my $name = @{$chr_annotations->{$sorted_by_end[0]}->{$per_end_start}}[0]->{name};
                if (defined($interval_records{$name}) and $interval_records{$name} > $overlap_length) {
                    next;
                }
                else {
                    #if so, save that as a possible interval, including the name of the annotation and the length of the interval
                    $interval_records{$name} = $overlap_length;
                }
            }
        }
        shift @sorted_by_end;
    }

    #Find the longest interval and report its length and name.
    my $maxClass;
    my $maxLength = 0;
    foreach my $key (keys %interval_records) {
        if (defined $interval_records{$key} and $interval_records{$key} > $maxLength) {
            $maxClass = $key;
            $maxLength = $interval_records{$key};
        }
    }
    return ($maxLength, $maxClass);
}

1;

