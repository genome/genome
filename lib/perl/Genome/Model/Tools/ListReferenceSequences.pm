package Genome::Model::Tools::ListReferenceSequences;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::ListReferenceSequences {
    is => 'Command',
    doc => "List available reference-sequence builds.",
    has => [
        dirs => { 
            is => 'boolean',
            default_value => 0,
            doc => "display absolute path to reference-sequence data directory",
        },
    ],
};

sub help_brief {
    'A lister of Reference Sequence builds.'
}

sub execute {
    my $self = shift;
    
        my @build = Genome::Model::Build->get(subclass_name => 'Genome::Model::Build::ImportedReferenceSequence');
        if($self->dirs == 1){
            printf "\n%-15s  %-70s %-25s %-75s \n", 'Build ID', 'Model Name','Version', 'Data Directory';
            print "===========================================================================================================================================\n";
        } else {
            
            printf "\n%-15s  %-70s %-25s\n", 'Build ID', 'Model Name','Version';
            print "==================================================================================================\n";
        }
        for (@build) {
            my $version;
            if(defined($_->version)){ 
                $version = $_->version;
            } else {
                $version = '';
            }
            if($self->dirs == 1){
                printf "%-15s  %-70s %-25s %-75s \n", $_->id, $_->model->name, $version, $_->data_directory;
            } else {
                printf "%-15s  %-70s %-25s\n", $_->id, $_->model->name, $version;
            }
        }    
        if(scalar(@build)==0){
            print " Warning: No Imported Reference Sequence Builds were found...\n";
        }
        print "\n";
}

1;

