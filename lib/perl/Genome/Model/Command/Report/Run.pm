package Genome::Model::Command::Report::Run;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Report::Run {
    is => 'Genome::Model::Command::BaseDeprecated',
    has => [
        report_name => {
            is => 'Text',
            doc => "the name of the report",
        },
        build   => {
            is => 'Genome::Model::Build',
            id_by => 'build_id',
            doc => "the specific build of a genome model on which to run",
        },
        build_id => {
            doc => 'the id for the build on which to run',
        },
    ],
    has_optional => [
        save_as => {
            is => 'Text',
            doc => 'instead of saving on the build, save to the specified directory (for testing)',
        },
        test => {
            is => 'Boolean',
            doc => 'runs the report in test mode, if supported',
        },
    ],
};

sub help_brief {
    "generate a report for a given build of a model"
}

sub help_synopsis {
    return <<EOS
genome model report run --build-id 12345 --report-name "dbsnp concordance"

genome model report run -b 12345 -r "dbsnp concordance"

genome model report run -b 12345 -r "dbsnp concordance" --save-as /tmp/mytest
EOS
}

sub help_detail {
    return <<EOS
This launcher generates (or re-generates) a report for some build of a genome model.

By default, the report is automatically linked to the build it reports on.  The
report is then viewable without re-running.  If the optional "save as" parameter is
specified, the report will be saved to disk, and not linked to the model build.

EOS
}

sub sub_command_sort_position { 1 }

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless $self;


    # TODO: move this up
    if (my $build = $self->build) {
        if ($self->model_name and $self->build->model->name ne $self->model_name) {
            $self->error_message("Redundant conflicting parameters.  Specified model name "
                . $self->model_name
                . " does not match the name of build " . $build->id
                . " which is " . $build->model->name
            );
            $self->delete;
            return;
        }
    }
    elsif ($self->model_name) {
        my @models = Genome::Model->get(name => $self->model_name);
        unless (@models) {
            $self->warning_message("No models have exact name " . $self->model_name);
            @models = Genome::Model->get("name like" => $self->model_name . "%");
            unless (@models) {
                $self->error_message("No models have a name beginning with " . $self->model_name);
                $self->delete;
                return;
            }
            if (@models > 1) {
                $self->error_message("Found multiple models with names like " . $self->model_name . ":\n\t"
                    . join("\n\t", map { $_->name } @models)
                    . "\n"
                );
                $self->delete;
                return;
            }
            $self->status_message("Found model " . $models[0]->id . " (" . $models[0]->name . ")");
            my $build = $models[0]->last_complete_build;
            unless ($build) {
                $self->error_message("No complete build for model!");
                $self->delete;
                return;
            }
            $self->status_message("Found build " . $build->id . " from " . $build->date_completed);
            $self->build($build);
        }
    }
    else {
        $self->error_message("A build must be specified either directly, or indirectly via a model name!");
        $self->delete;
        return;
    }

    return $self;
}

sub execute {
    my $self = shift;


    my $build = $self->build;
    my $report_name = $self->report_name;

    my $report_specification_class = $self->_resolve_valid_report_class_for_build_and_name($build,$report_name);

    unless ($report_specification_class) {
        $self->error_message("Unknown report '$report_name'!");
        return;
    }

    my $report_specification = $report_specification_class->create(build_id => $build->id);
    unless ($report_specification) {
        $self->error_message("Failed to initialize report!" . $report_specification_class->error_message);
        return;
    }

    if ($self->test) {
        if ($report_specification->can("test")) {
            $report_specification->test(1);
        }
        else {
            die "There is no test flag on the $report_specification_class report module!";
        }
    }

    my $report_output = $report_specification->generate_report;
    unless ($report_output) {
        $self->error_message("Error generating report!" . $report_specification->error_message);
        $report_specification->delete;
        return;
    }

    my $saved = 0;
    my $path;
    eval {
        if ($path = $self->save_as) {
            $self->status_message("Overriding save location to $path...");
            if (-d $path) {
                my $subdir = $path . '/' . $report_output->name_to_subdirectory($report_output->name);
                if (-e $subdir) {
                    $self->status_message("Sub-directory $subdir exists!   Moving it out of the way...");
                    my $n = 1;
                    while ($n < 10 and -e $subdir . '.' . $n) {
                        $n++;
                    }
                    Genome::Sys->rename($subdir, "$subdir.$n");
                    if (-e $subdir) {
                        die "failed to move old report dir $subdir to $subdir.$n!: $!";
                    }
                }
            }
            else {
                $self->status_message("creating directory $path...");
                unless (mkdir $path) {
                    die "failed to make directory $path!: $!";
                }
            }
            if ($report_output->save($path)) {
                $self->status_message("Saved report to override directory: $path");
                $saved = 1;
            }
            else {
                $self->error_message("Error saving report!: " . $report_output->error_message());
            }
        }
        else {
            if ($build->add_report($report_output)) {
                $self->status_message("Saved report to build directory " . $build->data_directory);
                $saved = 1;
            }
            else {
                $self->error_message("Error saving report to build!: " . $build->error_message());
            }
        }
    };

    # dump this before handling failed saves b/c if the later fails, we'll be glad this at least
    # made it to the terminal...
    #if (my $data = $report_output->get_data) {
    #    print $data->{txt},"\n";
    #}

    if ($saved) {
        # mail whoever ran the tools so they don't have to mail themselves
        my $me = Genome::Config->user_email();
        $self->status_message("Sending you the report (to $me)");
        my $r = Genome::Model::Command::Report::Mail->execute(
            model_id => $build->model_id,
            build_id => $self->build_id,
            report_name => $report_output->name,
            (
                $self->save_as
                    ? (directory => $self->save_as)
                    : ()
            ),
            to => $me,
        );
        unless ($r and $r->result) {
            $self->error_message("Failed to send email!");
        }
    }
    else {
        $self->error_message("Error saving report!");
        $self->error_message($@) if $@;
        my $build_id = $build->id;
        my $path = "/tmp/gm-build-errors-$build_id." . time() . '.' . $$;
        $self->warning_message("Saving to $path on " . Sys::Hostname::hostname() . ' instead for troubleshooting.');
        mkdir $path;
        unless ($report_output->save($path)) {
            $self->error_message("Error saving to temp! $!");
        }
        $self->warning_message("Save complete.");
        $report_output->delete;
        $report_specification->delete;
        return;
    }

    return 1;
}

# TODO: move onto the build or report as a method
sub _resolve_valid_report_class_for_build_and_name{
    my $self = shift;

    my $build = shift;
    my $report_name = shift;

    my $report_class_suffix =
        join('',
            map { ucfirst($_) }
            split(" ",$report_name)
        );

    my @model_classes = (
        grep { $_->isa('Genome::Model') }
        map { $_, $_->inheritance }
        $build->model->__meta__->class_name
    );

    my $report_class_meta;
    for my $model_class_name (@model_classes) {
        $report_class_meta = UR::Object::Type->get($model_class_name . '::Report::' . $report_class_suffix);
        last if $report_class_meta;
    }

    unless ($report_class_meta) {
        $self->error_message("No reports named '$report_name' are available for models of type "
            . join(", ", @model_classes));
        return;
    }

    return $report_class_meta->class_name;
}

1;
