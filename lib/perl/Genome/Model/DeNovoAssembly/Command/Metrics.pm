package Genome::Model::DeNovoAssembly::Command::Metrics;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::DeNovoAssembly::Command::Metrics {
    is => 'Command::V2',
    has => [
	    build => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
            doc => 'Build to trun metrics.',
        },
    ],
    has_optional => [
        update_the_db => {
            is => 'Boolean',
            doc => 'Update (or replace) the metrics for the build to the database.',
        },
    ],
};

sub help_brief {
    'Run metrics for a build, optionally saving to the db',
}

sub help_detail {
    return;
}

sub execute {
    my $self = shift;
    $self->debug_message('De Novo metrics...');

    $self->debug_message('Build: '.$self->build->id);
    $self->debug_message('Assembler: '.$self->build->processing_profile->assembler_name);
    my $metrics_class = $self->build->processing_profile->tools_base_class.'::Metrics';
    my $major_contig_length = ( $self->build->processing_profile->name =~ /PGA/ ? 300 : 500 );
    $self->debug_message('Assembly directory: '.$self->build->data_directory);
    $self->debug_message('Major contig length: '.$major_contig_length);
    my $metrics_tool = $metrics_class->create(
        assembly_directory => $self->build->data_directory,
        major_contig_length => $major_contig_length,
    );
    if ( not $metrics_tool ) {
        $self->error_message('Failed to create metrics tool: '.$metrics_class);
        return;
    }
    $metrics_tool->dump_status_messages(1);
    unless( $metrics_tool->execute ) {
        $self->error_message("Failed to execute stats");
        return;
    }

    if ( $self->update_the_db ) {
        $self->debug_message('Update db...');
        my $save = $self->_update_db($metrics_tool->_metrics);
        return if not $save;
    }

    $self->debug_message('Done');
    return 1;
}

sub _update_db {
    my ($self, $metrics) = @_;

    Carp::confess('No metrics to save!') if not $metrics;

    my %build_metrics = map { $_->name => $_ } $self->build->metrics;
    for my $name ( $self->build->metric_names ) {
        my $new_value = $metrics->$name;
        next if not $new_value;
        my $build_metric = $build_metrics{$name};
        if ( $build_metric ) {
            my $current_value = $build_metric->value;
            next if $current_value eq $new_value;
            $build_metric->delete;
        }
        $self->build->add_metric(
            name => $name,
            value => $new_value,
        );
    }

    return 1;
}

1;

