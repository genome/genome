package Genome::Model::MetagenomicComposition16s::Command::CopyFiles; 

use strict;
use warnings;

use Genome;

require File::Basename;
require File::Copy;

class Genome::Model::MetagenomicComposition16s::Command::CopyFiles {
    is => 'Genome::Model::MetagenomicComposition16s::Command',
    has => [
        file_type => {
            is => 'Text',
            doc => 'Type of file to copy/list.',
            valid_values => [qw/ oriented_fasta oriented_qual chimera_free_fasta chimera_free_qual processed_fasta processed_qual classification chimera_free_classification /],
        },
        destination => {
            is => 'Text',
            is_optional => 1,
            doc => 'The directory to copy the files.',
        },
        force => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'If destination files exist, overwrite.',
        },
        list => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => "List (don't copy) the builds files",
        },
    ],
};

sub sub_command_category { return; }

sub help_brief { 
    return 'List and copy files for MC16s models';
}

sub help_detail {
    return <<HELP;
    Copies files from builds (for each amplicon set) into a destination directory, optionally renaming them as it goes.

    To just see the files, use --list. 
    Use --force to overwrite existing files.

    This command is backward compatible for amplicon assembly builds. The oriented fastas are the same in MC16s and AA. The AA assembly fasta is the MC16s processed fasta. The classification files do not exist for amplicon assembly.

HELP
}

sub execute {
    my $self = shift;

    my $method;
    if ( $self->list ) {
        # list
        $method = '_list';
        $self->debug_message('List mc16s files');
    }
    else {
        # copy
        Genome::Sys->validate_existing_directory( $self->destination )
            or return;
        $method = '_copy';
        $self->debug_message('Copy mc16s files to '.$self->destination);
    }

    my $file_method = $self->file_type.'_file';
    my @builds = $self->_builds;
    return if not @builds;
    $self->debug_message('Found '.@builds.' builds');
    for my $build ( @builds ) {
        # aa backward compatibility
        if ( $build->type_name eq 'metagenomic composition 16s' ) {
            my @amplicon_sets = $build->amplicon_sets;
            unless ( @amplicon_sets ) {
                $self->error_message("No amplicon sets for ".$build->description);
                return;
            }
            for my $amplicon_set ( @amplicon_sets ) {
                my $rv = $self->$method($amplicon_set, $file_method);
                return if not $rv;
            }
        }
        elsif ( $build->type_name eq 'amplicon assembly' ) {
            my $rv = $self->$method($build, $file_method);
            return if not $rv;
        }
        else {
            die "Incompatible build type: ".$build->type_name;
        }
    }

    $self->debug_message('Done');

    return 1;
}

sub _list {
    my ($self, $amplicon_set, $file_method) = @_;

    my $file = $amplicon_set->$file_method;
    return 1 if not -e $file; 

    return print "$file\n";
}

sub _copy {
    my ($self, $amplicon_set, $file_method) = @_;

    my $from = $amplicon_set->$file_method;
    return 1 if not -e $from; 

    my $base_name = File::Basename::basename($from);
    my $dest = $self->destination.'/'.$base_name;

    if ( -e $dest ) {
        if ( $self->force ) {
            unlink $dest;
        }
        else {
            $self->error_message("Can't copy to $dest because it exists. Use --force option to overwrite existing files.");
            return;
        }
    }

    unless ( File::Copy::copy($from, $dest) ) {
        $self->error_message("Can't copy $from\nto $dest\nError: $!");
        return;
    }

    return 1;
}

1;

