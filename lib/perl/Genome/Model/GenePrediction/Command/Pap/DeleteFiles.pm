#$Id$

package Genome::Model::GenePrediction::Command::Pap::DeleteFiles;

use strict;
use warnings;

#use Workflow;

use English;

class Genome::Model::GenePrediction::Command::Pap::DeleteFiles {
    is  => ['Command::V1'],
    has => [
        files => { is => 'ARRAY', doc => 'array of files to delete',
                   is_input => 1, },
    ],
};


sub sub_command_sort_position { 10 }

sub help_brief {
    "Delete a set of files";
}

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}

sub execute {
    
    my $self = shift;

    
    my @files = @{$self->files()};

    foreach my $file (@files) {

    	if (-e $file) {
            unless(unlink($file)) {

                die "failed to unlink '$file': $OS_ERROR";

         	}	
        }

    }

    return 1;

}
 
1;
