package Genome::Model::Tools::BioSamtools::StatsSummary;

use strict;
use warnings;

use Genome;
use Statistics::Descriptive;

my %sort_order = (
    targets => 1,
    minimum_depth => 2,
    touched => 3,
    pc_touched => 4,
    target_base_pair => 5,
    covered_base_pair => 6,
    pc_target_space_covered => 7,
    mean_breadth => 8,
    stdev_breadth => 9,
    median_breadth => 10,
    targets_eighty_pc_breadth => 11,
    target_base_pair_eighty_pc_breadth => 12,
    covered_base_pair_eighty_pc_breadth => 13,
    pc_targets_eighty_pc_breadth => 14,
    pc_target_space_covered_eighty_pc_breadth => 15,
    mean_depth => 16,
    stdev_depth => 17,
    depth_quartile_3 => 18,
    median_depth => 19,
    depth_quartile_1 => 20,
    gaps => 21,
    mean_gap_length => 22,
    stdev_gap_length => 23,
    median_gap_length => 24,
);

class Genome::Model::Tools::BioSamtools::StatsSummary {
    is => ['Command'],
    has_input => [
        stats_file => {
            is => 'Text',
            doc => 'A STATS file output from refcov',
        },
        output_directory => {
            doc => 'When run in parallel, this directory will contain all output files. Do not define if output_file is defined.',
            is_optional => 1
        },
    ],
    has_output => [
        output_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'The output file to write summary stats',
        },
    ],
};

sub help_detail {
    'This command takes the STATS file from refcov and generates summary statistics based on all regions/targets.  This command will run on 32 or 64 bit architechture'
}

sub execute {
    my $self = shift;
    if ($self->output_directory) {
        unless (-d $self->output_directory){
            unless (Genome::Sys->create_directory($self->output_directory)) {
                die('Failed to create output directory '. $self->output_directory);
            }
        }
    }
    unless ($self->output_file) {
        my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->stats_file,qw/.tsv/);
        unless (defined($suffix)) {
            die('Failed to recognize stats_file '. $self->stats_file .' without a tsv suffix');
        }
        $self->output_file($self->output_directory .'/'. $basename .'.txt');
    }
    my $stats_fh = Genome::Sys->open_file_for_reading($self->stats_file);
    unless ($stats_fh) {
        die('Failed to read file stats file: '. $self->stats_file);
    }

    my $min_depth;
    my %depth_stats;
    my @headers;
    while (my $line = $stats_fh->getline) {
        chomp($line);
        if ($line =~ /^##/) { next; }
        my @entry = split("\t",$line);
        my $id = $entry[0];
        my $min_depth = $entry[12];
        unless ($min_depth =~ /^\d+$/) {
            #SKIP TEXT HEADER
            next;
        }
        unless (defined($depth_stats{$min_depth}{minimum_depth})) {
            $depth_stats{$min_depth}{minimum_depth} = $min_depth;
        }
        push @{$depth_stats{$min_depth}{breadth}}, $entry[1];
        if ($entry[1] >= 80) {
            $depth_stats{$min_depth}{target_base_pair_eighty_pc_breadth} += $entry[2];
            $depth_stats{$min_depth}{covered_base_pair_eighty_pc_breadth} += $entry[3];
            $depth_stats{$min_depth}{targets_eighty_pc_breadth}++;
        }
        $depth_stats{$min_depth}{target_base_pair} += $entry[2];
        if ($entry[3]) {
            $depth_stats{$min_depth}{covered_base_pair} += $entry[3];
            $depth_stats{$min_depth}{touched}++;
        }
        push @{$depth_stats{$min_depth}{depth}}, $entry[5];
        if ($entry[8]) {
            $depth_stats{$min_depth}{gaps} += $entry[8];
            push @{$depth_stats{$min_depth}{gap_length}}, $entry[9];
        }
        $depth_stats{$min_depth}{targets}++;
    }
    $stats_fh->close;
    for my $min_depth (sort {$a <=> $b} keys %depth_stats) {
        unless (defined($depth_stats{$min_depth}{covered_base_pair})) {
            $depth_stats{$min_depth}{covered_base_pair} = 0;
        }
        unless (defined($depth_stats{$min_depth}{targets_eighty_pc_breadth})) {
            $depth_stats{$min_depth}{targets_eighty_pc_breadth} = 0;
        }
        unless (defined($depth_stats{$min_depth}{touched})) {
            $depth_stats{$min_depth}{touched} = 0;
        }
        unless (defined($depth_stats{$min_depth}{gaps})) {
            $depth_stats{$min_depth}{gaps} = 0;
        }
        unless (defined($depth_stats{$min_depth}{target_base_pair_eighty_pc_breadth})) {
            $depth_stats{$min_depth}{target_base_pair_eighty_pc_breadth} = 0;
        }
        unless (defined($depth_stats{$min_depth}{covered_base_pair_eighty_pc_breadth})) {
            $depth_stats{$min_depth}{covered_base_pair_eighty_pc_breadth} = 0;
        }
        $depth_stats{$min_depth}{pc_touched} = sprintf( "%.03f", ( ( $depth_stats{$min_depth}{touched} / $depth_stats{$min_depth}{targets} ) * 100 ) );
        $depth_stats{$min_depth}{pc_target_space_covered} = sprintf("%.03f",(($depth_stats{$min_depth}{covered_base_pair}/$depth_stats{$min_depth}{target_base_pair})*100));
        $depth_stats{$min_depth}{pc_targets_eighty_pc_breadth} = sprintf("%.03f",(( ($depth_stats{$min_depth}{targets_eighty_pc_breadth} || 0) /$depth_stats{$min_depth}{targets})*100));
        $depth_stats{$min_depth}{pc_target_space_covered_eighty_pc_breadth} = sprintf("%.03f",(( ($depth_stats{$min_depth}{covered_base_pair_eighty_pc_breadth} || 0) /$depth_stats{$min_depth}{target_base_pair})*100));
        
        my $breadth_stat = Statistics::Descriptive::Full->new();
        my @breadth = delete($depth_stats{$min_depth}{breadth});
        $breadth_stat->add_data(@breadth);
        $depth_stats{$min_depth}{mean_breadth} = sprintf("%.03f",$breadth_stat->mean);
        $depth_stats{$min_depth}{stdev_breadth} = sprintf("%.10f",$breadth_stat->standard_deviation);
        $depth_stats{$min_depth}{median_breadth} = sprintf("%.03f",$breadth_stat->median);
        my $depth_stat = Statistics::Descriptive::Full->new();
        my @depth = delete($depth_stats{$min_depth}{depth});
        $depth_stat->add_data(@depth);
        $depth_stats{$min_depth}{mean_depth} = sprintf("%.03f",$depth_stat->mean);
        $depth_stats{$min_depth}{stdev_depth} = sprintf("%.10f",$depth_stat->standard_deviation);
        $depth_stats{$min_depth}{depth_quartile_3} = sprintf("%.03f",$depth_stat->percentile(75));
        $depth_stats{$min_depth}{median_depth} = sprintf("%.03f",$depth_stat->median);
        $depth_stats{$min_depth}{depth_quartile_1} = sprintf("%.03f",$depth_stat->percentile(25));
        my $gaps_stat = Statistics::Descriptive::Full->new();
        my @gaps = delete($depth_stats{$min_depth}{gap_length});
        $gaps_stat->add_data(@gaps);
        $depth_stats{$min_depth}{mean_gap_length} = sprintf("%.03f",$gaps_stat->mean);
        $depth_stats{$min_depth}{stdev_gap_length} = sprintf("%.10f",$gaps_stat->standard_deviation);
        $depth_stats{$min_depth}{median_gap_length} = sprintf("%.03f",$gaps_stat->median);
        unless (@headers) {
            @headers = sort hash_sort_order (keys %{$depth_stats{$min_depth}});
        }
    }
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@headers,
        output => $self->output_file,
    );
    unless ($writer) {
        my $error_message = Genome::Utility::IO::SeparatedValueWriter->error_message;
        $error_message .= 'Failed to open writer handle for output file '. $self->output_file ." with headers:\n\t". join("\n\t",@headers);
        $self->error_message($error_message);
        die($self->error_message);
    }
    for my $min_depth (sort {$a <=> $b } keys %depth_stats) {
        unless ($writer->write_one($depth_stats{$min_depth})) {
            my $error_message = $writer->error_message;
            $writer->output->close;
            unless (unlink($self->output_file)) {
                $error_message .= 'Failed to remove output file '. $self->output_file ."\n";
            }
            $self->error_message($error_message);
            die($self->error_message);
        }
    }
    $writer->output->close;
    return 1;
}

sub hash_sort_order {
    unless ($sort_order{$a}) {
        die('sort order not defined for '. $a);
    }
    unless ($sort_order{$b}) {
        die('sort order not defined for '. $b);
    }
    $sort_order{$a} <=> $sort_order{$b};
}

1;
