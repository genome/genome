package Genome::Model::DeNovoAssembly::Command::ExportAsGscProject;

use strict;
use warnings;

use File::Basename;

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
    'VELVET ONE-BUTTON' => {
        subdirs => [qw/ edit_dir chromat_dir phd_dir/],
        edit_dir => 'edit_dir',
    },
    'NEWBLER DE-NOVO-ASSEMBLE' => {
        subdirs => [qw/ edit_dir chromat_dir phd_dir phdball_dir /],
        edit_dir => File::Spec->join('consed', 'edit_dir'),
    },
);
sub supported_assemblers { keys %supported_assemblers }

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

    my @succeeded_builds = map {
        $_->last_succeeded_build;
    } grep {
        #next if $ARGV[2] && ! exists $subjects_to_copy{ $model->subject->name };
        ! model_is_newbler_assembly_with_multiple_inst_data($_);
    } @models;

    $self->status_message('Succeeded builds: %s', scalar(@succeeded_builds));
    $self->fatal_message('No succeeded builds found!') if not @succeeded_builds;

    return \@succeeded_builds;
}

sub _export_build {
    my ($self, $build) = @_;
    return 1;

    my $build_dir = $build->data_directory;
    die "No data directory found, $build_dir\n" unless -d $build_dir;

    # project dir
    my $copy_dir = $self->directory;
    my $output_dir = $copy_dir.'/'.$build->model->subject->name;
    $output_dir =~ s/\s+/_/g;
    Genome::Sys->create_directory( $output_dir ) unless -d $output_dir;

    # create build id empty file
    my $build_id_file = $output_dir."/build_id=".$build->id;
    Genome::Sys->touch( $build_id_file ) if not -e $build_id_file;

    my $assembler = $build->model->processing_profile->assembler_name;
    # create subdirs
    my @sub_dirs = subdirs_for_assembler( $assembler );
    die "Can't determine which subdirs to create for assembler, $assembler\n" unless @sub_dirs;
    for my $dir_name( @sub_dirs ) {
        my $sub_dir = $output_dir."/$dir_name";
        Genome::Sys->create_directory( $sub_dir ) unless -d $sub_dir;
    }

    # edit_dir files to copy
    my $builds_edit_dir = assemblers_edit_dir( $assembler );
    die "Can't determine locations of builds edit_dir for assembler, $assembler\n" unless $builds_edit_dir;
    for my $file( glob("$build_dir/$builds_edit_dir/*") ) {
        my $file_name = basename( $file );
        my $to = $output_dir.'/edit_dir/'.$file_name;
        print "skipping, $file already exists\n" and next if -s $to;
        print "Copying $file to $to\n";
        Genome::Sys->copy_file( $file, $to );
    }

    my @addl_files = additional_files_for_assembler( $assembler, $build );
    if( @addl_files ) {
        while( @addl_files ) {
            my( $from, $to ) = split(/\s+/, shift @addl_files);
            $from = $build_dir.'/'.$from;
            $to   = $output_dir.'/'.$to;
            die "Can't find file, $from to link\n" unless -s $from;
            print "Linking $from to $to\n";
            Genome::Sys->create_symlink( $from, $to ) if not -l $to;
        }
    }

    return 1;
}

sub model_is_newbler_assembly_with_multiple_inst_data {
    my $model = shift;
    my @inst_data = $model->instrument_data;
    my $assembler = $model->processing_profile->assembler_name;
    return 1 if $assembler =~ /newbler/i && @inst_data > 1;
    return;
}

sub additional_files_for_assembler {
    my $assembler = shift;
    my $build = shift;
    my @addl_files;
    if( uc $assembler eq 'NEWBLER DE-NOVO-ASSEMBLE' ) {
        my @velvet_input_fastqs = velvet_input_fastq_files( $build );
        for( @velvet_input_fastqs ) {
            my $file_name = basename( $_ );
            push @addl_files, $file_name." edit_dir/$file_name";
        }
        push @addl_files, 'consed/phdball_dir/phd.ball.1 phdball_dir/phd.ball.1';
    }
    return @addl_files;
}

sub velvet_input_fastq_files {
    my $build = shift;
    my @input_fastqs = glob( $build->data_directory."/*input.fastq" );
    return @input_fastqs;
}

sub subdirs_for_assembler {
    my $assembler = shift;
    return (qw/ edit_dir chromat_dir phd_dir/) if uc $assembler eq 'VELVET ONE-BUTTON';
    return (qw/ edit_dir chromat_dir phd_dir phdball_dir /) if uc $assembler eq 'NEWBLER DE-NOVO-ASSEMBLE';
    return;
}

sub assemblers_edit_dir {
    my $assembler = shift;
    return 'edit_dir' if uc $assembler eq 'VELVET ONE-BUTTON';
    return 'consed/edit_dir' if uc $assembler eq 'NEWBLER DE-NOVO-ASSEMBLE';
    return;
}

1;

