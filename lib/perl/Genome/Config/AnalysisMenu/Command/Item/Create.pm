package Genome::Config::AnalysisMenu::Command::Item::Create;

use strict;
use warnings;

use Genome;
use File::Spec;
use File::Basename;

class Genome::Config::AnalysisMenu::Command::Item::Create {
   is => 'Genome::Command::Base',
   has => [
        name => {
            is => 'Text',
            doc => 'name of new AnalysisMenu Item',
        },
        file_path => {
            is => 'Path',
            doc => 'Path to config file',
        },
        description => {
            is => 'Text',
            doc => 'Description of the menu item',
        }
   ],
   doc => 'create Analysis Menu Items',
};

sub execute {
    my $self = shift;
    my $name = $self->name;
    my $file_path = $self->file_path;
    my $description = $self->description;

    unless($self->_validate_file($file_path)){
        die $self->error_message("$file_path is invalid");
    }
    my $allocated_file_path = $self->_copy_file_to_allocation($file_path);
    my $item = Genome::Config::AnalysisMenu::Item->create(name => $name, file_path => $allocated_file_path, description => $self->description);
    $self->status_message("Successfully created AnalysisMenu Item $name for $file_path");
    return 1;
}

sub _validate_file {
    my $self = shift;
    my $file_path = shift;
#TODO: once an actual validator for these files is written, swap this out for a call to that tool
    return -s $file_path;
}

sub _copy_file_to_allocation {
    my $self = shift;
    my $original_file = shift;
    my $destination = File::Spec->join($ENV{GENOME_ANALYSIS_PROJECT_DEFAULTS}, $self->_determine_file_name($original_file));
    unless( Genome::Sys->copy_file($original_file, $destination) ){
        die $self->error_message("Failed to copy $original_file to $destination : $@");
    }
    return $destination;
}

sub _determine_file_name {
    my $self = shift;
    my $file = shift;
    my ($file_name, undef, undef) = File::Basename::fileparse($file);
    return join('_', time, $file_name);
}

1;
