package Genome::Model::Tools::Bedpe::EvaluateBedpe;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Bedpe::EvaluateBedpe {
    is => 'Command::V2',
    has_input => [
        bedpe => {
            is => 'File',
            doc => 'bedpe file to be evaluated',
        },
        gold_bedpe => {
            is => 'Path',
            doc => 'bedpe file with gold standard breakpoints'
        },
        bedtools_version => {
            is => 'Text',
        },
        slop => {
            is => 'Integer',
            doc => 'The amount of slop (in b.p.). to be added to one set of breakpoints',
        },
        min_hit_support => {
            is => 'Integer',
            doc => 'Minimum # of hits from gold_bedpe required to report an sv as a true positive.  Note, if != 1, the derived stats are probably not valid',
            default_value => 1,
        },
        true_positive_file => {
            is => 'Path',
            is_optional => 1,
            doc => 'Output the true positive hits to this file',
        },
    ],
    has_transient_optional_output => [
        rawstats => {
            is => "HASH",
            doc => "The raw stats generated during primary execution",
        },
    ],
};

sub execute {
    my $self = shift;
    $self->rawstats({});
    $self->rawstats->{true_positive} = $self->_get_true_positives();
    $self->rawstats->{false_negative} = $self->_get_false_negatives();
    $self->rawstats->{total_unique_calls} = $self->_unique_sv_count($self->bedpe);
    $self->rawstats->{total_unique_gold_calls} = $self->_unique_sv_count($self->gold_bedpe);
    $self->rawstats->{false_positive} = $self->rawstats->{total_unique_calls} - $self->rawstats->{true_positive};
    $self->_set_derivative_stats;
    $self->print_stats;
    return 1;
}

sub common_params {
    my $self = shift;
    return (
        slop => $self->slop,
        slop_strand => "+-",
        ignore_strand => 1,
        require_different_names => 0,
        use_version => $self->bedtools_version,
    );
}

sub _get_pair_to_pair_output {
    my ($self, $file_a, $file_b, $type) = @_;

    my $output_file = Genome::Sys->create_temp_file_path;

    Genome::Model::Tools::BedTools::PairToPair->execute(
        $self->common_params,
        output_file => $output_file,
        input_file_a => $file_a,
        input_file_b => $file_b,
        intersection_type => $type,
    );
    return $output_file;
}

sub _get_false_negatives {
    my $self = shift;
    my $output_file = $self->_get_pair_to_pair_output($self->gold_bedpe, $self->bedpe, 'notboth');
    return $self->_unique_sv_count($output_file);
}

sub _get_true_positives {
    my $self = shift;
    my $output_file = $self->_get_pair_to_pair_output($self->bedpe, $self->gold_bedpe, 'both');
    return $self->_sv_with_min_support_count($output_file, $self->bedpe, $self->true_positive_file);
}

sub _sv_with_min_support_count {
    my ($self, $file, $original_file, $output) = @_;
    my $f = Genome::Sys->open_file_for_reading($file);
    my %results;

    my $id_index = _count_cols($original_file) + 6;

    my $out;
    if ($output) {
        $out = Genome::Sys->open_file_for_writing($output);
    }
    while (my $line = $f->getline) {
        chomp $line;
        my @f = split("\t", $line);
        my $sv_name = $f[6];
        my $read_id = $f[$id_index];
        unless (defined $sv_name) {
            die $self->error_message("Sv name not defined on line:\n$line");
        }
        unless (defined $read_id) {
            die $self->error_message("read id not defined on line:\n$line");
        }
        ++$results{$sv_name}{count} unless exists $results{$sv_name}{reads}{$read_id};
        $results{$sv_name}{reads}{$read_id} = 1;
        $results{$sv_name}{sv} = \@f;
    }
    $f->close();

    my $in = Genome::Sys->open_file_for_reading($original_file);
    my %seen;
    my $count = 0;
    while (my $line = $in->getline) {
        chomp $line;
        my @f = split("\t", $line);
        my $sv_name = $f[6];
        unless (defined $sv_name) {
            die $self->error_message("Sv name not defined on line:\n$line");
        }
        next if exists $seen{$sv_name};
        $seen{$sv_name} = 1;

        next unless exists $results{$sv_name} && ($results{$sv_name}{count} >= $self->min_hit_support);
        $count++;
        if ($out) {
            print $out "$line\n";
        }
    }
    if ($out) {
        $out->close;
    }
    return $count;
}

sub _count_cols {
    my $path = shift;
    my $f = Genome::Sys->open_file_for_reading($path);
    my $line = $f->getline;
    $f->close();
    my @f = split("\t", $line);
    return scalar @f;
}

sub _unique_sv_count {
    my ($self, $file) = @_;
    my $count = `cut -f 1-6 $file | sort -u | wc -l`;
    chomp $count;
    return $count;
}

sub _set_derivative_stats {
    my $self = shift;
    my $tp = $self->rawstats->{true_positive};
    my $fp = $self->rawstats->{false_positive};
    my $fn = $self->rawstats->{false_negative};

    $self->rawstats->{ppv} = $tp/($tp + $fp);
    $self->rawstats->{sensitivity} = $tp/($tp + $fn);
    $self->rawstats->{f1} = 2*$tp/(2*$tp + $fp + $fn);
}

sub print_stats {
    my $self = shift;

    for my $stat (sort keys %{$self->rawstats}) {
        printf("%s\t%s\n", $stat, $self->rawstats->{$stat});
    }
}
1;

