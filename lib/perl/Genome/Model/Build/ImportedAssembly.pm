package Genome::Model::Build::ImportedAssembly;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::ImportedAssembly {
    is => 'Genome::Model::Build',
};

sub create {
    my $class = shift;
    
    my $self = $class->SUPER::create(@_)
	or return;
    #CHECK TO MAKE SURE DATA DIRECTORY EXISTS
    unless (-d $self->data_directory) {
	    $self->error_message("Assembly data directory does not exist:\n\t".$self->data_directory);
	    $self->delete;
	    return;
    }
    #CHECK OF ALREADY EXISTING BUILDS WITH SAME ASSEMBLY DIRECTORY
    #THERE SHOULDN'T BE MULTIPLE BUILDS WITH SAME DATA DIRECTORIES
    foreach my $build ($self->model->builds) {
        next if $build->build_id == $self->build_id;
	    if ($build->data_directory eq $self->data_directory) {
	        $self->error_message("A build with data directory: ".$self->data_directory." already exists\n\tBUILD ID: ".
				     $build->build_id."\n\tDIR: ".$build->data_directory);
		$self->delete;
	        return;
	    }
    }

    return $self;
}



