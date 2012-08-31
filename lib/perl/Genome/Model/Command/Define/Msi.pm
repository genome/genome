package Genome::Model::Command::Define::Msi;#manual sequence improvement

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Define::Msi {
    is => 'Genome::Model::Command::Define::HelperDeprecated',
    has => [
        msi_assembly => { #the assembly we are going to perform msi on
            is => 'Genome::Model', 
            id_by => 'input_assembly',
            doc => ''             
        },
        input_assembly => { 
            is => 'Text',
            doc => 'The input assembly that will be improved',
            is_optional => 1,
        },
        
    ],
    has_optional => [
        subject_name =>
        {
            type => 'String',
            is_optional => 1,
            doc => "this parameter is for internal use, any value specified will be over-ridden"
        },
   ],
};

sub help_synopsis {
    return <<"EOS"
genome model define msi
  --input-assembly 54321
EOS
}

sub help_detail {
    return <<"EOS"
This defines a new genome model representing an msi pipeline 
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    return $self;
}

sub execute {
    my $self = shift;


    unless(defined $self->input_assembly) {
        $self->error_message("Could not get a model for input assembly id: " . $self->msi_assembly);
        return;
    }
    
    my $msi_assembly = $self->msi_assembly;
    $self->subject_name($msi_assembly->subject_name);
    $self->subject_type($msi_assembly->subject_type);
     
     # run Genome::Model::Command::Define execute
    my $super = $self->super_can('_execute_body');
    $super->($self,@_);    

    my $model = Genome::Model->get($self->result_model_id);    

    $model->add_from_model(from_model => $self->msi_assembly, role => 'imported_assembly');
        
    return 1; 
}

1;
