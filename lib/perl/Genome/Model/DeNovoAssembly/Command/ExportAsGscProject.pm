package Genome::Model::DeNovoAssembly::Command::ExportAsGscProject;

use strict;
use warnings;

use File::Basename;
use Params::Validate qw/ :types validate_pos /;

class Genome::Model::DeNovoAssembly::Command::ExportAsGscProject {
    is => 'Command::V2',
    has => {
        directory => {
            is => 'Text',
            doc => 'Directory to put the project.',
        },
        project => {
            is => 'Genome::Project',
            doc => 'Work Order',
        },
    },
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    if ( not -d $self->directory ) {
        return UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ directory /],
            desc => 'Directory does not exist: '.$self->directory,
        );
    }

    return;
}

my %supported_assemblers = (
    'velvet one-button' => {
        subdirs => [qw/ edit_dir chromat_dir phd_dir/],
        edit_dir => 'edit_dir',
    },
    'newbler de-novo-assemble' => {
        subdirs => [qw/ edit_dir chromat_dir phd_dir phdball_dir /],
        edit_dir => File::Spec->join('consed', 'edit_dir'),
    },
);
sub supported_assemblers { keys %supported_assemblers }
sub is_model_supported {
    my ($self, $model) = validate_pos(@_, {isa => __PACKAGE__}, {type => OBJECT});

    my $assembler = $model->processing_profile->assembler_name;
    return if ! List::MoreUtils::any { $assembler eq $_ } supported_assemblers();

    my @inst_data = $model->instrument_data;
    return if $assembler =~ /newbler/i && @inst_data > 1;
    return 1;
    #return if ! exists $subjects_to_copy{ $model->subject->name };
}

sub subdirs_for_assembler {
    my ($self, $assembler) = validate_pos(@_, {isa => __PACKAGE__}, {type => SCALAR});
    return if not $supported_assemblers{$_[1]};
    @{$supported_assemblers{$_[1]}->{subdirs}};
}

sub assemblers_edit_dir {
    my ($self, $assembler) = validate_pos(@_, {isa => __PACKAGE__}, {type => SCALAR});
    return if not $supported_assemblers{$_[1]};
    $supported_assemblers{$_[1]}->{edit_dir};
}

sub additional_files_for_assembler {
    my ($self, $build) = validate_pos(@_, {isa => __PACKAGE__}, {isa => 'Genome::Model::Build'});

    my @addl_files;
    if( $build->model->processing_profile->assembler_name eq 'newbler de-novo-assemble' ) {
        for my $file ( glob( File::Spec->join($build->data_directory, "*input.fastq")) ) {
            push @addl_files, [ $file, 'edit_dir' ];
        }
        push @addl_files, [ File::Spec->join($build->data_directory, 'consed', 'phdball_dir', 'phd.ball.1'), File::Spec->join('phdball_dir') ];
    }

    return @addl_files;
}

sub execute {
    my $self = shift;

    my $succeeded_builds = $self->_resolve_builds;
    for my $build ( @$succeeded_builds ) {
        $self->status_message("Export: %s %s", $build->model->__display_name__, $build->data_directory);
        $self->_export_build($build);
    }

    return 1;
}

sub _resolve_builds {
    my $self = shift;

    my $project = $self->project;
    $self->status_message('Resolving builds for project: %s', $project->__display_name__);

    my @project_parts = $project->parts(entity_class_name => 'Genome::Model::DeNovoAssembly');
    $self->fatal_message("No de novo models associated with %s", $project->__display_name__) unless @project_parts;
    $self->status_message('Associated Models: %s', scalar(@project_parts));

    my @models = map { $_->entity } @project_parts;
    $self->fatal_message("Models associated with %s do not exist!", $project->__display_name__) unless @models;
    $self->status_message('Existing Models: %s ', scalar(@models));

    my @succeeded_builds;
    for my $model ( @models ) {
        next if ! $self->is_model_supported($model);
        push @succeeded_builds, $model->last_succeeded_build;
    }

    $self->status_message('Succeeded builds: %s', scalar(@succeeded_builds));
    $self->fatal_message('No succeeded builds found!') if not @succeeded_builds;

    return \@succeeded_builds;
}

sub _export_build {
    my ($self, $build) = @_;

    my $build_dir = $build->data_directory;
    $self->fatal_message("Build data directory does not exist! $build_dir") unless -d $build_dir;

    # project dir
    my $copy_dir = $self->directory;
    my $output_dir = File::Spec->join($copy_dir, $build->model->subject->name);
    $output_dir =~ s/\s+/_/g;
    Genome::Sys->create_directory($output_dir) unless -d $output_dir;

    # create build id empty file
    my $build_id_file = File::Spec->join($output_dir, "build_id=".$build->id);
    Genome::Sys->touch($build_id_file) if not -e $build_id_file;

    # create subdirs
    my $assembler = $build->model->processing_profile->assembler_name;
    my @sub_dirs = $self->subdirs_for_assembler($assembler);
    $self->fatal_message("Can't determine which subdirs to create for assembler! $assembler") unless @sub_dirs;
    for my $dir_name( @sub_dirs ) {
        my $sub_dir = File::Spec->join($output_dir, $dir_name);
        Genome::Sys->create_directory($sub_dir) unless -d $sub_dir;
    }

    # edit_dir files to copy
    my $builds_edit_dir = $self->assemblers_edit_dir($assembler);
    $self->fatal_message("Can't determine locations of builds edit_dir for assembler! $assembler") unless $builds_edit_dir;
    for my $file( glob( File::Spec->join($build_dir, $builds_edit_dir, '*') ) ) {
        my $file_name = basename( $file );
        my $to = File::Spec->join($output_dir, 'edit_dir', $file_name);
        $self->status_message("skipping, $file already exists") and next if -s $to;
        $self->status_message("Copying $file to $to");
        Genome::Sys->copy_file( $file, $to );
    }

    for my $files ( $self->additional_files_for_assembler($build) ) {
        my ($target_path, $destination_subdir) = @$files;
        $self->fatal_message("Link target does not exist! $target_path") unless -s $target_path;
        my $destination_name = File::Basename::basename($target_path);
        my $link_path = File::Spec->join($output_dir, $destination_subdir, $destination_name);
        $self->status_message("Linking $link_path to $target_path");
        Genome::Sys->create_symlink($target_path, $link_path) if not -l $link_path;
    }

    return 1;
}

1;

