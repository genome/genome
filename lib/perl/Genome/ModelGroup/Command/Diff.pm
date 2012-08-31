package Genome::ModelGroup::Command::Diff;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::Diff {
    is => 'Genome::Command::Base',
    has_input => [
        from => {
            is => 'Genome::ModelGroup',
            id_by => 'from_id',
            shell_args_position => 1,
            doc => 'the first group to examine', 
        },
        to => {
            is => 'Genome::ModelGroup',
            id_by => 'to_id',
            shell_args_position => 2,
            doc => 'the second group to examine', 
        }
    ],
    has_optional_input => {
        data_sets => {
            shell_args_position => 3,
            is => 'Text',
            is_many => 1,
            doc => 'files, directories or metrics to compare on each build'
        },
    },
    doc => 'make a new model group from another, varying properties on the model as specified' 
};

sub help_synopsis {
    return <<EOS
genome model-group compare group1 group2 variants

EOS
}


sub execute {
    my $self = shift;
    my $from = $self->from;
    my $to = $self->to;
    my @data_sets = $self->data_sets;

    $self->status_message("\nselecting builds...\n");

    my @from_models = sort { $a->subject_id <=> $b->subject_id || $a->id <=> $b->id } $from->models;
    my @to_models   = sort { $a->subject_id <=> $b->subject_id || $a->id <=> $b->id } $to->models;
    my @cache_all_builds = Genome::Model::Build->get(model_id => [map { $_->id } ( @from_models,@to_models) ]);

    #print scalar(@from_models),"\n@from_models\n\n";
    #print scalar(@to_models),"\n@to_models\n\n";
    #print scalar(@cache_all_builds),"\n@cache_all_builds\n\n";

    my %selected_build;
    for my $m (@from_models, @to_models) {
        my @b = sort { $b->id <=> $a->id } $m->builds;
        $self->status_message("\n" . $m->__display_name__ . ":");
        for my $b (@b) {
            my $dir = $b->data_directory;
            #system "echo $dir; ls $dir\n";
            my @missing;
            for my $data_set (@data_sets) {
                unless (-e "$dir/$data_set") {
                    $self->status_message(" $dir/$data_set is missing...");
                    push @missing, "$dir/$data_set";        
                }
            }
            unless (@missing) {
                $self->status_message( " > selecting " . $b->__display_name__ . " for model " . $m->__display_name__);
                $selected_build{$m->id} = $b;
                last;
            }
        }
        unless ($selected_build{$m->id}) {
            $self->status_message(" > no builds with appropriate data sets for model " . $m->__display_name__ . "!!!!");
        }
    }

    unless (@data_sets) {
        my %all;
        for my $build (values %selected_build) {
            my $dir = $build->data_directory;
            my @ds = map { s/$dir\///; $_ => 1 } glob("$dir/*");
            %all = (%all, @ds);
        }
        @data_sets = sort keys %all;
        $self->status_message("default to these data sets: @data_sets");
    }

    $self->status_message("\ncomparing...\n\n");

    while (@from_models) {
        my $from_model = shift @from_models;
        my $to_model = shift @to_models;
        
        my $from_build = $selected_build{$from_model->id};
        my $to_build = $selected_build{$to_model->id};

        if (!$from_build and !$to_build) {
            $self->status_message("# no build for either model " . $from_model->__display_name__ . " or " . $to_model->__display_name__);
            next;
        }
        elsif (! $from_build) {
            $self->status_message("# no build for model " . $from_model->__display_name__);
            next;
        }
        elsif (! $to_build) {
            $self->status_message("# no build for model " . $to_model->__display_name__);
            next;
        }

        my $from_dir = $from_build->data_directory;
        my $to_dir = $to_build->data_directory;

        for my $data_set (@data_sets) {
            my $from_path = $from_dir . '/' . $data_set;
            my $to_path = $to_dir . '/' . $data_set;
            
            if (not -e $from_path and not -e $to_path) {
                $self->status_message("MISSING BOTH: $from_path and $to_path");
                next;
            }
            elsif (not -e $from_path) {
                $self->status_message("MISSING: $from_path");
                next;
            }
            elsif (not -e $to_path) {
                $self->status_message("MISSING: $to_path");
                next;
            }

            my $cmd;
            if (-d $from_path) {
                $cmd = "diff -r --brief $from_path $to_path";
            }
            else {
                $cmd = "sdiff -s $from_path $to_path";
            }

            my @diff = `$cmd 2>&1`;
            if (@diff) {
                my $cnt = scalar(@diff);
                if (-d $from_path ) {
                    $self->status_message("DIFF ($cnt): $cmd\n");
                }
                else {
                    my @lines;
                    my $from_lines = (-e $from_path ? scalar(@lines = IO::File->new($from_path)->getlines) : 0);
                    my $to_lines = (-e $to_path ? scalar(@lines = IO::File->new($to_path)->getlines) : 0);
                    $self->status_message("DIFF ($from_lines vs $to_lines with $cnt differences): $cmd\n");
                }
            }
            else {
                $self->status_message("SAME: $cmd\n");
            }
        }
    }

    return 1;
}

1;
