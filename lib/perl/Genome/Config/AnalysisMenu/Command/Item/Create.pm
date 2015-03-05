package Genome::Config::AnalysisMenu::Command::Item::Create;

use strict;
use warnings;

use Genome;
use File::Spec;
use File::Basename;

class Genome::Config::AnalysisMenu::Command::Item::Create {
   is => 'Command::V2',
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
    my $item = Genome::Config::AnalysisMenu::Item->create(name => $name, file_path => $file_path, description => $self->description);
    $self->status_message("Successfully created AnalysisMenu Item $name for $file_path");
    return 1;
}

sub _validate_file {
    my $self = shift;
    my $file_path = shift;
    my ($filename, $dir) = File::Basename::fileparse($file_path);
    $dir =~ s!/$!!;
    unless ($dir eq $ENV{'GENOME_ANALYSIS_PROJECT_DEFAULTS'}){
        return 0;
    }
#TODO: once an actual validator for these files is written, swap this out for a call to that tool
    return -s $file_path;
}

1;
