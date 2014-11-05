package Genome::Config::AnalysisProject::Command::View;

use strict;
use warnings FATAL => 'all';

use Genome;
use Genome::Utility::Text qw(justify);

use List::Util qw(max sum);
use YAML::Syck qw();


class Genome::Config::AnalysisProject::Command::View {
    is => [
        'Genome::Command::Viewer',
        'Genome::Command::ColorMixin',
    ],

    has_input => [
       analysis_project  => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis projects to view',
            shell_args_position => 1,
        },
    ],

    has_optional_input => [
        COLUMN_WIDTH => {
            is => 'Number',
            default_value => 40,
        },
        config => {
            is => 'Text',
            valid_values => ['verbose', 'terse', 'quiet'],
            default_value => 'terse',
            doc => 'How to configuration items',
        },
        disabled_configs => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Whether to display disabled config items',
        },
        instrument_data => {
            is => 'Boolean',
            default_value => 1,
            doc => 'Display instrument data summary',
        },
        models => {
            is => 'Boolean',
            default_value => 1,
            doc => 'Display model status summary',
        },
        timeline => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Display timeline events',
        },
        fast_model_summary => {
            is => 'Boolean',
            default_value => 0,
            doc => "Use fast model summary that does not due 'Build Needed' check",
        },
    ],
};


my %STATUS_COLORS = (
    # shared values
    failed => "red",
    new => "white",

    # values for AnalysisProjects
    approved => "blue",
    archived => "magenta",
    completed => "green",
    hold => "red",
    inprogress => "cyan",
    pending => "yellow",

    # values for InstrumentData
    processed => "green",
    skipped => "magenta",

    # values for model status
    buildless => "yellow",
    buildneeded => "blue",
    running => "cyan",
    scheduled => "white",
    succeeded => "green",
    unstartable => "magenta",

    # values for config items
    active => "green",
    disabled => "magenta",
);
my $MAX_STATUS_WIDTH = max(map {length($_)} keys %STATUS_COLORS) + 1;


sub write_report {
    my ($self, $width, $handle) = @_;

    $self->_write_heading($handle);
    $self->_write_instrument_data_summary($handle);
    $self->_write_model_summary($handle);
    $self->_write_config_items($handle);

    $self->_write_timeline($handle);

    1;
}


sub _write_heading {
    my ($self, $handle) = @_;

    $self->_write_section_heading($handle, 'Analysis Project');
    for my $line ($self->_get_heading_lines) {
        $self->_write_pairs_line($handle, @$line);
    }

    if (lc $self->analysis_project->status eq 'hold') {
        print $handle "\n    ",
            $self->_color('** Analysis Project on HOLD **', 'red'), "\n";
    } elsif (lc $self->analysis_project->status eq 'pending') {
        print $handle "\n    ",
            $self->_color('** Analysis Project is PENDING **', 'red'), "\n";
    }

    print $handle "\n";
}


sub _write_section_heading {
    my ($self, $handle, $text) = @_;

    print $handle $self->_color_heading($text), "\n";
}


sub _get_heading_lines {
    my $self = shift;

    my $ap = $self->analysis_project;
    return (
        ["ID", $ap->id, "Name", $ap->name],
        ["Run as", $ap->run_as, "Model Group", $ap->_model_group_id],
        ["Created", $ap->created_at, "Updated", $ap->updated_at],
        ["Created by", $ap->created_by],
    );
}


sub _write_pairs_line {
    my ($self, $handle, $l_label, $l_value, $r_label, $r_value) = @_;

    if ($r_label and $r_value) {
        print $handle justify($self->_color_pair($l_label, $l_value), 'left',
            $self->COLUMN_WIDTH), " ", $self->_color_pair($r_label, $r_value), "\n";

    } else {
        print $handle $self->_color_pair($l_label, $l_value), "\n";
    }
}


sub _write_instrument_data_summary {
    my ($self, $handle) = @_;

    if ($self->_should_write_instrument_data) {
        $self->_write_section_heading($handle, 'Instrument Data');

        my $summary = $self->_get_instrument_data_summary;
        $self->_write_summary_data($handle, $summary);
    }
}


sub _write_summary_data {
    my ($self, $handle, $summary) = @_;

    my @sorted_classes = sort keys %$summary;
    for my $class (@sorted_classes) {
        print $handle $self->_color($class, 'white'), "\n";

        $self->_write_status_table($handle, $summary->{$class});
        print $handle "\n";
    }
}


sub _write_status_table {
    my ($self, $handle, $statuses) = @_;

    my @sorted_statuses = sort keys %$statuses;
    for my $status (@sorted_statuses) {
        $self->_write_status_table_row($handle, $self->_color_status($status),
            $statuses->{$status});
    }
    $self->_write_status_table_row($handle, $self->_color('Total', 'blue'),
        sum(values %$statuses));
}


sub _write_status_table_row {
    my ($self, $handle, $left, $right) = @_;

    print $handle justify($left, 'right', $MAX_STATUS_WIDTH + 4), ' ',
        $right, "\n";
}

sub _should_write_instrument_data {
    my $self = shift;

    unless ($self->instrument_data) {
        return;
    } else {
        my @bridges = $self->analysis_project->analysis_project_bridges;
        return scalar(@bridges);
    }
}


sub _get_instrument_data_summary {
    my $self = shift;

    return $self->_summary_data_from_query($self->_instrument_data_query);
}


sub _summary_data_from_query {
    my ($self, $query) = @_;

    my $dbh = Genome::Model->__meta__->data_source->get_default_handle;

    my $query_object = $dbh->prepare($query);
    $query_object->execute();

    my $summary = {};
    while (my $row = $query_object->fetchrow_arrayref()) {
        my ($class, $status, $count) = @$row;
        $summary->{$class}->{$status} = $count;
    }
    $query_object->finish();

    return $summary;
}


my $INST_DATA_QUERY_TEMPLATE = <<EOS;
SELECT id_class, id_status, count(*) FROM (
    SELECT
        instrument_data.subclass_name as id_class,
        bridge.status as id_status
    FROM instrument.data instrument_data LEFT JOIN
        config.instrument_data_analysis_project_bridge bridge
        ON instrument_data.id = bridge.instrument_data_id
    WHERE
        bridge.analysis_project_id = '%s'
) counts
GROUP BY id_class, id_status;
EOS
sub _instrument_data_query {
    my $self = shift;
    return sprintf($INST_DATA_QUERY_TEMPLATE, $self->analysis_project->id);
}


sub _write_model_summary {
    my ($self, $handle) = @_;

    if ($self->_should_write_models) {
        $self->_write_section_heading($handle, 'Models');

        my $summary = $self->fast_model_summary
                    ? $self->_get_fast_model_summary
                    : $self->_get_model_summary;
        $self->_write_summary_data($handle, $summary);
    }
}


sub _should_write_models {
    my $self = shift;

    unless ($self->models) {
        return;
    } else {
        my @model_bridges = $self->analysis_project->model_bridges;
        return scalar(@model_bridges);
    }
}


sub _get_model_summary {
    my $self = shift;

    my $summary = {};
    my $model_iterator = Genome::Model->create_iterator(
        'analysis_projects.id' => $self->analysis_project->id,
    );
    while (my $model = $model_iterator->next) {
        $summary->{$model->class}->{$model->status}++;
    }
    return $summary;
}


sub _get_fast_model_summary {
    my $self = shift;

    return $self->_summary_data_from_query($self->_model_summary_query);
}


my $MODEL_QUERY_TEMPLATE = <<EOS;
SELECT model_class, model_status, count(*) FROM (
    SELECT
        model.genome_model_id as model_id,
        model.subclass_name AS model_class,
        CASE
            WHEN model.build_requested THEN 'Build Requested'
            WHEN build.status IS NULL THEN 'Buildless'
            ELSE build.status
        END AS model_status
    FROM model.model model LEFT JOIN (
        SELECT ROW_NUMBER() OVER (PARTITION BY model_id ORDER BY created_at DESC) AS r, model_id, status
        FROM model.build
        WHERE status != 'Abandoned'
        AND build.model_id IN (%s)
    ) build ON model.genome_model_id = build.model_id
    WHERE model.genome_model_id IN (%s)
    AND (build.r = 1 OR build.r IS NULL)
) counts
GROUP BY model_class, model_status;
EOS

sub _model_summary_query {
    my $self = shift;

    my @model_ids = $self->_model_ids;
    my $joined_model_ids = join(', ', map {"'$_'"} @model_ids);

    return sprintf($MODEL_QUERY_TEMPLATE, $joined_model_ids, $joined_model_ids);
}


sub _model_ids {
    my $self = shift;

    return map {$_->model_id} $self->analysis_project->model_bridges;
}


sub _write_timeline {
    my ($self, $handle) = @_;

    if ($self->timeline) {
        $self->_write_section_heading($handle, 'Timeline Events');

        my @events = Genome::Timeline::Event::AnalysisProject->get(
            analysis_project => $self->analysis_project);
        for my $event (@events) {
            $self->_write_timeline_event($handle, $event);
        }
    }
}


sub _write_timeline_event {
    my ($self, $handle, $event) = @_;

    print $handle $self->_color_pair($event->updated_at,
        sprintf("%s (%s): %s\n", $event->name, $event->created_by,
            $event->reason));
}


sub _color_status {
    my ($self, $text) = @_;
    return $self->_colorize_text_by_map($text, $text, %STATUS_COLORS);
}


sub _write_config_items {
    my ($self, $handle) = @_;
    my @config_items = $self->analysis_project->config_items;
    if (@config_items) {
        if ($self->config eq 'terse') {
            $self->_write_terse_config_items($handle);

        } elsif ($self->config eq 'verbose') {
            $self->_write_verbose_config_items($handle);
        }
    }
}


sub _write_terse_config_items {
    my ($self, $handle) = @_;

    $self->_write_section_heading($handle, 'Configuration Items');

    for my $config_item ($self->_config_items) {
        $self->_write_config_item_heading($handle, $config_item);
    }
}


sub _config_items {
    my $self = shift;

    if ($self->disabled_configs) {
        return $self->analysis_project->config_items(
            '-order_by' => 'created_at');
    } else {
        return $self->analysis_project->config_items(
            'status' => 'active', '-order_by' => 'created_at');
    }
}

sub _write_config_item_heading {
    my ($self, $handle, $config_item) = @_;

    $self->_write_menu_item($handle, $config_item->analysis_menu_item);
    $self->_write_config_item_properties($handle, $config_item);

    print $handle "\n";
}


sub _write_menu_item {
    my ($self, $handle, $menu_item) = @_;
    if ($menu_item) {
        printf $handle "%s (%s)",
            $self->_color($menu_item->name, 'white'),
            $menu_item->id;

        if ($menu_item->description) {
            print $handle ": ", $menu_item->description;
        }

        print $handle "\n";

    } else {
        print $handle $self->_color('Custom configuration item', 'white'), "\n";
    }
}


sub _write_config_item_properties {
    my ($self, $handle, $config_item) = @_;

    for my $line ($self->_get_config_item_lines($config_item)) {
        print $handle '    ';
        $self->_write_pairs_line($handle, @$line);
    }

}

sub _get_config_item_lines {
    my ($self, $config_item) = @_;

    return (
        ['ID', $config_item->id,
            'Concrete', $self->_get_concrete_string($config_item)],
        ['Created by', $config_item->created_by,
            'Status', $self->_color_status(
                $config_item->status)],
        ['Created', $config_item->created_at,
            'Updated', $config_item->updated_at],
    );
}


sub _get_concrete_string {
    my ($self, $config_item) = @_;

    if ($config_item->is_concrete) {
        return $self->_color('Yes', 'green');

    } else {
        return $self->_color('No', 'yellow');
    }
}


sub _write_verbose_config_items {
    my ($self, $handle) = @_;

    $self->_write_section_heading($handle, 'Configuration Items');

    for my $config_item ($self->_config_items) {
        $self->_write_config_item_heading($handle, $config_item);
        $self->_write_config_item_body($handle, $config_item);
    }
}


sub _write_config_item_body {
    my ($self, $handle, $config_item) = @_;

    my $config_data = Genome::Config::Parser->parse($config_item->file_path);
    my $yaml = YAML::Syck::Dump($config_data);

    $yaml =~ s/^/    /gm;

    print $handle $yaml, "\n";
}


1;
