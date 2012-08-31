package Genome::Model::Command::Export::Pcgp;

use strict;
use warnings;
use Genome;

class Genome::Model::Command::Export::Pcgp { 
    is => 'Genome::Command::Base',
    has => [
        paths => { shell_args_position => 1, is => 'FilePath', is_many => 1, doc => 'path(s) to BAMs' }, 
    ],
    doc => 'prepare a directory of symlinks to BAMs for upload to St. Jude'
};

sub help_synopsis {
    return <<EOS
genome model export pcgp /gscmnt/ams1114/info/model_data/2862280028/build104502554/alignments/104502554_merged_rmdup.bam

genome model export pcgp 104502554 104502513 104654215 104566131 

EOS
}

sub help_detail {
    return <<'EOS'
Take one or more bam file paths, or build ids, and create a "pcgp-upload" directory with symlinks to those files.

The name of the symlink will be $PATIENT-$TUMORNORMAL, where $PATIENT is the patient common name, 
and $TUMORNORMAL is the sample common name.

Currently, actually running the uploader is a manual process.

EOS
}

sub execute {
    my $self = shift;

    my @paths_in = $self->paths;

    my @paths;
    my @builds;
    my $errors;
    for my $path_in (@paths_in) {
        my $build_id;
        if ($path_in =~ /^(\d+)$/) {
            $build_id = $1;
        }
        elsif ($path_in =~ /build(\d+)/) {
            $build_id = $1;
        }
        elsif ($path_in =~ /^(\d+)_merged_rmdup.bam/) {
            $build_id = $1;
        }
        else {
            $self->error_message("no build ID found for path $path_in!");
            $errors++;
            next;
        }
        my $build = Genome::Model::Build->get($build_id);
        unless ($build) {
            $self->error_message("build $build_id not found?!");
            $errors++;
            next;
        }
        my $data_directory = $build->data_directory;
        unless (-d $data_directory)  {
            $self->error_message("no data directory " . $data_directory . " found?");
            $errors++;
            next;
        }
        my $bam_path = $data_directory . '/alignments/' . $build_id . '_merged_rmdup.bam';
        unless (-e $bam_path) {
            $self->error_message("missing $bam_path????");
            $errors++;
            next;
        }
        $self->status_message("build $build_id has BAM path $bam_path");
        push @paths, $bam_path;
        push @builds, $build;
    }
    
    if ($errors) {
        $self->status_message("exiting due to errors");
        return;
    }

    my @names;
    for my $build (@builds) {
        my $model = $build->model;
        my $sample = $model->subject;
        unless ($sample) {
            $self->error_message("no sample for build " . $build->__display_name__);
            $errors++;
            next;
        }
        unless ($sample->common_name) {
            $self->error_message("no common name for sample " . $sample->__display_name__); 
            $errors++;
            next;
        }
        my $patient = $sample->patient;
        unless ($patient) {
            $self->error_message("no patient for sample " . $sample->__display_name__);
            $errors++;
            next;
        }
        unless ($patient->common_name) {
            $self->error_message("no common name for patient " . $patient->__display_name__);
            $errors++;
            next;
        }
        my $name = $patient->common_name . "_" . $sample->common_name;
        $self->status_message("build " . $build->id . " gets name $name");
        push @names, $name;
    }

    if ($errors) {
        $self->status_message("exiting due to errors");
        return;
    }

    if (-d "pcgp-upload") {
        $self->status_message("found pcgp-upload directory");
    }
    else {
        Genome::Sys->create_directory("pcgp-upload");
        $self->status_message("created directory pcgp-upload");
    }

    my @links;
    while (@names) {
        my $name = shift @names;
        my $path = shift @paths;
        my $link = "pcgp-upload/$name";
        if (-l $link) {
            if (my $alt = readlink($link) eq $path) {
                $self->warning_message("SYMLINK EXISTS AND IS CORRECT: $path $link");
                next;
            }
            else {
                $self->warning_message("SYMLINK CONFLICT: $link points to " . $alt . " but we want it to point to $path!!!!");
                $errors++;
            }
        }
        push @links, [$path,$link];
    }

    for my $link_params (@links) {
        my ($path,$link) = @$link_params;
        $self->status_message("symlink: $link $path");
        Genome::Sys->create_symlink($path, $link);
    }

    return 1;
}

1;

