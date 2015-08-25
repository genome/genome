package Genome::Model::Tools::Bedpe::EvaluateBedpes;

use strict;
use warnings;
use Genome;
use JSON qw(to_json);

class Genome::Model::Tools::Bedpe::EvaluateBedpes {
    is => 'Command::V2',
    has_input => [
        bedtools_version => {
            is => 'Text',
        },
        config_file => {
            is => 'Path',
            doc => 'A tsv file which contains the following columns: bedpe, gold_bedpe, slop (along with any other columns of metadata you wish to include)',
        },
        output_json => {
            is => "Path",
            doc => "Output json file",
            is_output => 1,
        },
    ],
};

sub execute {
    my $self = shift;
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->config_file,
        separator => "\t",
    );
    my @output;
    while (my $line = $reader->next()) {
        my $true_positive_file = Genome::Sys->create_temp_file_path;
        my $stats = $self->_run_one(
            $line->{bedpe},
            $line->{gold_bedpe},
            $line->{slop},
            $line->{min_hit_support},
            $true_positive_file,
        );
        for my $key (keys %$stats) {
            die "Duplicate key $key: this would overwrite the column provided in the config"
                if exists $line->{$key};
            $line->{$key} = $stats->{$key};
        }

        if ($line->{include_tps}) {
            $line->{true_positives} = [_read_tps($true_positive_file)];
        }
        delete $line->{bedpe};
        delete $line->{gold_bedpe};
        delete $line->{include_tps};
        push @output, $line;
    }
    Genome::Sys->write_file($self->output_json,
        to_json(\@output, {canonical => 1, pretty => 1}));
    return 1;
}

sub _read_tps {
    my $file = shift;
    my @tps = Genome::Sys->read_file($file);
    chomp @tps;
    return @tps;
}

sub _run_one {
    my ($self, $bedpe, $gold_bedpe, $slop, $min_hit_support, $true_positive_file) = @_;
    my $cmd = Genome::Model::Tools::Bedpe::EvaluateBedpe->create(
        bedpe => $bedpe,
        gold_bedpe => $gold_bedpe,
        slop => $slop,
        bedtools_version => $self->bedtools_version,
        min_hit_support => $min_hit_support || 1,
        true_positive_file => $true_positive_file,
    );
    $cmd->execute;
    return $cmd->rawstats;
}

1;

