package Genome::Model::DeNovoAssembly::Command::CopyVelvetBuild;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::DeNovoAssembly::Command::CopyVelvetBuild {
    is => 'Genome::Command::Base',
    has => [
	to => {
	    is => 'Text',
	    doc => 'Where to copy the build to',
	},
    ],
    has_optional => [
	build_id => {
	    is => 'Number',
	    doc => 'Id of build to copy',
	},
	build_directory => {
	    is => 'Text',
	    doc => 'Path of the assembly to copy',
	},
    ],
};

sub help_brief {
    'Command to make a copy of a velvet build to another location for manual assembly manupliation',
}

sub help_detail {
    return <<"EOS"
genome model de-novo-assembly copy-velvet-build --build-id 54958478 --to /gscmnt/111/assembly
genome model de-novo-assembly copy-velvet-build --build-directory /gscmnt/111/2857912274/build105415496 --to /gscmnt/111/assembly
EOS
}

sub execute {
    my $self = shift;

    #dir to copy from
    my $from_dir = $self->_get_build_directory;
    #dir to copy to
    my $to_dir = $self->_resolve_to_directory;

    unless ( $self->_copy_build( $from_dir, $to_dir ) ) {
	$self->error_message("Failed to copy assembly in from dir to to dir");
	return;
    }

    return 1;
}

sub _get_build_directory {
    my $self = shift;
    
    unless ( $self->build_id or $self->build_directory ) {
	$self->error_message("You must supply either a build id or build directory");
	return;
    }

    my $build_directory;

    if ( $self->build_directory ) {
	unless ( -d $self->build_directory ) {
	    $self->error_message("Can't find build directory: ".$self->build_directory);
	    return;
	}
	$build_directory = $self->build_directory;
    }
    elsif ( $self->build_id ) {
	my $build = Genome::Model::Build->get( $self->build_id );
	unless ( $build ) {
	    $self->error_message("Can't get a build for id: ".$self->build_id);
	    return;
	}
	unless ( -d $build->data_directory ) {
	    $self->error_message("Failed to find build directory: ".$build->data_directory);
	}

	$build_directory = $build->data_directory;
    }

    return $build_directory;
}

sub _resolve_to_directory {
    my $self = shift;

    unless ( -d $self->to ) {
	$self->status_message("Did not find directory: ".$self->to." will create it");
	Genome::Sys->create_directory( $self->to );
    }

    unless ( -d $self->to ) {
	$self->error_message("Failed to create directory to copy assembly to: ".$self->to);
	return;
    }
    
    return return $self->to;
}

sub _velvet_files_to_link {
    return (qw/ Graph2 LastGraph Log PreGraph Roadmaps Sequences collated.fastq contigs.fa stats.txt velvet_asm.afg /);
}

sub _copy_build {
    my ($self, $from, $to) = @_;

    foreach my $file ( $self->_velvet_files_to_link ) {
	#check file exists
	unless ( -e $from . "/$file" ) {
	    $self->status_message("Could not find file, $file in build directory, skipping linking of this file");
	    return; #too strict, this files are really not needed??
	}
	#link file
	Genome::Sys->create_symlink( $from . "/$file", $to . "/$file" );
	#check that file has been linked
	unless ( -l $to . "/$file" ) {
	    $self->error_message("Failed to link file, $file to to-directory");
	    return;
	}
    }

    unless ( -d $from."/edit_dir" ) {
	$self->error_message("Failed to find edit_dir in build directory");
	return;
    }

    Genome::Sys->copy_directory( $from . '/edit_dir', $to );

    unless ( -d $to . '/edit_dir' ) {
	$self->error_message("Failed to copy build edit_dir to to-directory");
	return;
    }

    $self->status_message("Completed copying of build to directory: $to");

    return 1;
}


1;
