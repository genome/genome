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
            doc => 'A tsv file which contains the following columns: bedpe, gold_bedpe, slop (along with any other columns of metadata you wish to include',
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
        my $stats = $self->_run_one(
            $line->{bedpe},
            $line->{gold_bedpe},
            $line->{slop},
        );
        for my $key (keys %$stats) {
            $line->{$key} = $stats->{$key};
        }
        delete $line->{bedpe};
        delete $line->{gold_bedpe};
        push @output, $line;
    }
    Genome::Sys->write_file($self->output_json,
        to_json(\@output, {canonical => 1, pretty => 1}));
    return 1;
}

sub _run_one {
    my ($self, $bedpe, $gold_bedpe, $slop) = @_;
    my $cmd = Genome::Model::Tools::Bedpe::EvaluateBedpe->create(
        bedpe => $bedpe,
        gold_bedpe => $gold_bedpe,
        slop => $slop,
        bedtools_version => $self->bedtools_version,
    );
    $cmd->execute;
    return $cmd->rawstats;
}

1;

