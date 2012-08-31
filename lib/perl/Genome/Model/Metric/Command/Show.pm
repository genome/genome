package Genome::Model::Metric::Command::Show;

use strict;
use warnings;

use Genome;

class Genome::Model::Metric::Command::Show {
    is => 'Command::V2',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            is_input => 1,
            is_optional => 1,
            require_user_verify => 0,
            doc => 'Models to show metrics. Uses last complete build. Resolved via string from the command line.',
        },
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            is_input => 1,
            is_optional => 1,
            require_user_verify => 0,
            doc => 'Builds to show metrics. Resolved via string from the command line.',
        },
        build_attributes => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            is_input => 1,
            doc => 'Build attributes to list before the metrics. Default is model_name.'
        },
        separator => {
            is => 'Text',
            is_optional => 1,
            is_input => 1,
            default_value => 'TAB',
            doc => 'Separator. Use "TAB" for tab separator.',
        }
    ],
    doc => 'show metrics by model',
};

sub help_detail {
    return 'This will show metrics for models last completed build. They are organized by model. Whereas the lister organizes by metric.'
}

sub execute {
    my $self = shift;

    my %builds;
    for my $model ( $self->models ) {
        my $build = $model->last_complete_build;
        next if not $build;
        $builds{ $build->id } = $build;
    }

    for my $build ( $self->builds ) {
        $builds{ $build->id } = $build;
    }

    if ( not %builds ) {
        $self->error_message('No builds given!');
        return;
    }
    $self->status_message('Builds: '.keys(%builds));

    $self->build_attributes([qw/ model_name /]) if not $self->build_attributes;
    my @build_attributes = $self->build_attributes;

    my @metrics;
    for my $build ( map { $builds{$_} } sort { $a <=> $b } keys %builds ) {
        my @build_metrics = $build->metrics;
        push @metrics, { 
            map({ $_, $build->$_ } @build_attributes),
            map({ $_->name, $_->value } @build_metrics), 
        };
    }

    my %names = map { $_ => 1 } map { keys %{$_} } @metrics;
    delete @names{ @build_attributes };
    my @names = sort { $a cmp $b } keys %names;
    unshift @names, @build_attributes;
    my @headers = @names;
    @headers = map { s/\s+/_/g; $_; } @headers;

    my $sep = $self->separator;
    $sep = "\t" if $sep eq "TAB";
    print join($sep, @headers)."\n";
    for my $metric ( sort { $a->{$build_attributes[0]} cmp $b->{$build_attributes[0]} } @metrics ) {
        print join($sep, map { ( defined $metric->{$_} ? $metric->{$_} : 0 ) } @names)."\n";
    }

    return 1;
}

1;

