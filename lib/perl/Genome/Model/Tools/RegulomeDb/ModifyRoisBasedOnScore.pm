package Genome::Model::Tools::RegulomeDb::ModifyRoisBasedOnScore;

use strict;
use warnings;
use Genome;
use WWW::Mechanize;

class Genome::Model::Tools::RegulomeDb::ModifyRoisBasedOnScore {
    is => 'Genome::Model::Tools::RegulomeDb',
    has => [
        roi_list => {
            is => 'String',
            doc => 'Path to roi list in bed format (0-based)',
        },
        output_file => {
            is => 'String',
            doc => 'Path of output file',
        },
        valid_scores => {
            is => 'String',
            is_many => 1,
            doc => 'Filter ROIs, keeping only these scores',
            valid_values => [qw(1 2 3 4 5 6 1a 1b 1c 1d 1e 1f 2a 2b 2c 3a 3b)],
        },
    ],
};

sub execute {
    my $self = shift;

    my $roi_string = Genome::Sys->read_file($self->roi_list);
    my @rois = split(/\n/, $roi_string);
    my @modified_rois;
    while (@rois) {
        my @part = splice(@rois, 0, 1000);
        my $part_string = join("\n", @part);
        my $expanded_rois = $self->expand_rois($part_string);
        my $output = $self->fetch_large_annotation(
            'bed', $expanded_rois 
        );
        my @lines = split(/\n/, $output);
        for my $roi (@part) {
            my $modified_roi = $self->process_roi($roi, @lines);
            if ($modified_roi) {
                push @modified_rois, $modified_roi;
            }
        }
    }
    Genome::Sys->write_file($self->output_file, join("\n", @modified_rois));
    return 1;
}

sub expand_rois {
    my $self = shift;
    my $roi_string = shift;
    my @expanded_rois;
    my @rois = split(/\n/, $roi_string);
    for my $roi (@rois) {
        my @fields = split(/\t/, $roi);
        my @extra_fields = @fields[3 .. $#fields];
        for (my $i = $fields[1]; $i < $fields[2]; $i++) {
            push @expanded_rois, join("\t", $fields[0], $i, $i+1, @extra_fields);
        }
    }

    return join("\n", @expanded_rois);
}

sub process_roi {
    my $self = shift;
    my $roi = shift;
    my @lines = @_;
    my @modified_roi;
    my $current_start;
    my $current_stop;
    my $expanded_rois = $self->expand_rois($roi);
    my @positions = split(/\n/, $expanded_rois);
    for my $pos (@positions) {
        my @pos_fields = split(/\t/, $pos);
        my $pos_only = join("\t", @pos_fields[0..2]);
        my @relevant_lines = grep {/^chr$pos_only/} @lines;
        my $num_lines = scalar @relevant_lines;
        unless($num_lines == 1) {
            $self->error_message("There should be exactly one score line per position");
            die $self->error_message;
        }
        if ($self->good_score($relevant_lines[0])) {
            if ($current_start and $current_stop) {
                $current_stop = $pos;
            }
            else {
                $current_start = $pos;
                $current_stop = $pos;
            }
        }
        else {
            if ($current_start and $current_stop) {
                push @modified_roi, $self->combine_rois($current_start, $current_stop);
                $current_start = undef;
                $current_stop = undef;
            }
        }
    }
    if ($current_start and $current_stop) {
        push @modified_roi, $self->combine_rois($current_start, $current_stop);
    }
    return join("\n", @modified_roi);
}

sub good_score {
    my $self = shift;
    my $line = shift;
    my @fields = split(/\t/, $line);
    my @score = split(/;/, $fields[3]);
    return grep {$score[1] =~ /$_/} $self->valid_scores;
}

sub combine_rois {
    my $self = shift;
    my $roi1 = shift;
    my $roi2 = shift;
    my @fields1 = split(/\t/, $roi1);
    my @fields2 = split(/\t/, $roi2);
    unless ($fields1[0] eq $fields2[0]) {
        $self->error_message("Can't combine ROI on different chromosomes");
        die $self->error_message;
    }
    $fields1[2] = $fields2[2];
    my $combined_roi = join("\t", @fields1);
    return $combined_roi;
}
1;

