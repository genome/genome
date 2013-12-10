package Genome::Model::Command::Report::Mail;

use strict;
use warnings;
use Mail::Sender;
use Genome;

class Genome::Model::Command::Report::Mail {
    is => 'Genome::Model::Command::BaseDeprecated',
    has => [
        report_name => {
            is => 'Text',
            doc => "the name of the report to mail",
        },
        build   => {
            is => 'Genome::Model::Build',
            id_by => 'build_id',
            doc => "the specific build of a genome model to mail a report for",
        },
        build_id => {
            doc => 'the id for the build on which to report',
        },
        to => {
            is => 'Text',
            doc => 'the recipient list to mail the report to',
        },

    ],
    has_optional => [
        directory => {
            is => 'Text',
            doc => 'the path of report directory to mail (needed only if the report was saved to a non-default location)',
        }
    ],
};

sub help_brief {
    "mail a report for a given build of a model"
}

sub help_synopsis {
    return <<EOS
genome model report mail --build-id 12345 --report-name "Summary" --to user1\@example.com,user2\@example.com

genome model report run -b 12345 -r "DbSnp" --to user\@example.com

genome model report run -b 12345 -r "GoldSnp" --to user\@example.com --directory /homedirs/username/reports

EOS
}

sub help_detail {
    return <<EOS
This launcher mails a report for some build of a genome model.

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

$DB::single = $DB::stopper;

    my $build = $self->build;
    my $model = $build->model;
    my $report_name = $self->report_name;

    my $report_path = $self->directory;
    unless ( defined($report_path) ) {
        $report_path = $build->resolve_reports_directory;
    }

    unless (-d $report_path) {
        $self->error_message("Failed to find report directory $report_path!");
        return;
    }

    my @mail_parts;

    my $report_subdir = $self->report_name;
    $report_subdir =~ s/ /_/g;

    unless (-d "$report_path/$report_subdir") {
        my @others = glob("$report_path/*");
        $self->error_message(
            "Failed to find subdirectory $report_subdir under $report_path!\n"
            . "Found:\n\t" . join("\n\t",@others) . "\n"
        );
        return;
    }

    my $html_file = $report_path."/".$report_subdir."/report.html";
    my $txt_file = $report_path."/".$report_subdir."/report.txt";

    #Need at least one part to send!
    my $file_count = 0;
    if (-e $html_file) {
       #print("\nFound html file.\n");
       $file_count++;
    }
    if (-e $txt_file) {
       #print("\nFound txt file.\n");
    $file_count++;
    }

    if ( $file_count == 0 ) {
        my @others = glob("$report_path/$report_subdir/*");
        $self->error_message(
            "Failed to find either text or html reports at the expected paths under $report_path/$report_subdir!"
            . "\nFound:\n\t" . join("\n\t",@others) . "\n"
        . "Expected HTML report at: $html_file\n"
        . "Expected Text report at: $txt_file"
        );
        return;
    }

    $report_name =~ s/_/ /g;
    my $subject = 'Genome Model '.$model->id.' "'.$model->name.'" '.$report_name.' Report for Build '.$build->id;

    $self->status_message("Sending email...");
    if ($self->send_mail($subject,$html_file,$txt_file)) {
        $self->status_message("Email sent.");
        return 1;
    }

    $self->error_message("Email not sent!");
    return;
}

sub send_mail {

    my $self  = shift;
    my $subject = shift;
    my $msg_html = shift;
    my $msg_txt = shift;

    my $recipients = $self->to;

    eval {

        my $sender = Mail::Sender->new(
                    {
                            smtp => 'gscsmtp.wustl.edu',
                            to => $recipients,
                            from => 'apipe@genome.wustl.edu',
                            subject => $subject,
                            multipart => 'related',
                            on_error => 'die',
                    }
            );

        $sender->OpenMultipart;

        $sender->Part({ctype => 'multipart/alternative'});

        if (-e $msg_txt) {
            my $msg_txt_contents = get_contents($msg_txt);
            $sender->Part({ctype => 'text/plain', disposition => 'NONE', msg => $msg_txt_contents})
        }
        if (-e $msg_html) {
            my $msg_html_contents = get_contents($msg_html);
            $sender->Part({ctype => 'text/html', disposition => 'NONE', msg => $msg_html_contents})
        }

        $sender->EndPart("multipart/alternative");

        $sender->Close;
    };

    if ($@) {
        $self->error_message("Error sending mail!: $@");
        return;
    }

    return 1;
}

sub get_contents {
   my $in_file_name = shift;
   my $in = new IO::File($in_file_name, "r");
   #print ("\nGetting contents for : $in_file_name\n");
   my $ret = "";
   while (<$in>) {
      $ret.= $_;
   }
   $in->close();
   return $ret;
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
        grep { $_->isa('Genome::Model') and $_ ne 'Genome::Model' }
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
