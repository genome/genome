package Genome::Model::Tools::EdgeR::Base;

use Genome;
use Carp qw/confess/;

use strict;
use warnings;

class Genome::Model::Tools::EdgeR::Base {
    is => "Command::V2",
    has => [
        counts_file => {
            is => "Text",
            doc => "Tab delimited file of read counts per structure",
        },

        groups => {
            is => "Text",
            doc => "Comma separated list of conditions corresponding to each " .
                "column in counts_file (e.g, normal,tumor,...)",
        },

        output_file => {
            is => "Text",
            doc => "Output file path",
        },

        p_value => {
            is => "Number",
            doc => "The p-value used for significance testing (0 < p < 1)",
            default_value => 0.05,
        },
    ],
};

sub _validate_params {
    my $self = shift;

    if ($self->p_value <= 0 or $self->p_value >= 1) {
        confess "--p-value parameter must satisfy 0 < p < 1";
    }

    my @groups = split(",", $self->groups);
    my %group_sizes;
    for my $g (@groups) {
        return if ++$group_sizes{$g} > 1;
    }

    confess sprintf("At least one group must contain replication! Groups " .
        "assignments were: %s\n",
        $self->groups);
}

sub execute {
    my $self = shift;

    $self->_validate_params;

    my $cmd = $self->construct_r_command;

    return Genome::Sys->shellcmd(
        cmd => $cmd
        );
}

sub construct_r_command {
    my $self = shift;

    $self->error_message("R command must be specified in the child class");

    return 0;
}

sub generate_counts_file {
    my ($path, $n_genes, $n_normal, $n_tumor) = @_;

    my $fh = Genome::Sys->open_file_for_writing($path);

    srand(1234);
    # genrate some baseline counts
    my @counts = map { 10 + int(rand(15)) } 1..$n_genes;
    my @columns;

    # normal samples
    for my $nid (1..$n_normal) {
        # add a little noise to the counts
        my @new_counts = map {int($_ + rand(2) - 1)} @counts;
        push(@columns, \@new_counts);
    }

    # tumor samples
    for my $tid (1..$n_tumor) {
        # add a little noise to the counts
        my @new_counts = map {int($_ + rand(2) - 1)} @counts;
        # but also greatly increase the expression of gene #0
        $new_counts[0] += 100;
        push(@columns, \@new_counts);
    }

    my @column_names = (
        "Gene",
        (map {"Normal$_"} 1..$n_normal),
        (map {"Tumor$_"} 1..$n_tumor)
        );

    $fh->write(join("\t", @column_names) . "\n");

    for my $i (0..$n_genes - 1) {
        my @fields = ("GENE$i", map {$_->[$i]} @columns);
        $fh->write(join("\t", @fields) . "\n");
    }
    $fh->close;
}



1;
