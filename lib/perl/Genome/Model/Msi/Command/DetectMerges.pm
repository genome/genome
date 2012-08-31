package Genome::Model::Msi::Command::DetectMerges; 

use strict;
use warnings;

use Genome;

require Workflow::Simple;


class Genome::Model::Msi::Command::DetectMerges {
    is => 'Command',
    has => [
        assembly_build_id => { 
            is => 'Number', 
            shell_args_position => 1,
            is_optional => 1,
            doc => 'the exact build of denovo assembly, imported assembly, or prior msi to examine' 
        },
    ],
    doc => 'identify possible merges in an assembly (usable in the merge command)'
};

sub execute {
    my $self = shift;

    my $assembly_build_id = $self->assembly_build_id;
    unless ($assembly_build_id) {
        # is the user in a build directory?
        # TODO: infer the build from it..
        
        unless ($assembly_build_id) {
            # still not set
            $self->error_message("No assembly build specified, and unable to infer the assembly build from the current directory!");
            return;
        }
    }

    if ($assembly_build_id =~ /\D/) {
        $self->error_message("The specified assembly build id is not a number!: $assembly_build_id");
        return;
    }
    
    my $assembly_build = Genome::Model::Build->get($assembly_build_id);
    unless ($assembly_build) {
        $self->error_message("failed to find a build with id $assembly_build_id!");    
        return;
    }
    
    $self->status_message("Found assembly build " . $assembly_build->__display_name__ . " for model " . $assembly_build->model->__display_name__);

    # get the data directory, run the tool
    my $data_directory = $assembly_build->data_directory;
    my $edit_dir = $data_directory . "/edit_dir";
   
	my @ace_files = `ls  $edit_dir/*.ace*`;    
	chomp @ace_files;

    #filter out singleton acefiles
    @ace_files = grep { !/singleton/ } @ace_files;  
    
    #filter out any non ace files
    
    my @versioned_ace_files;
    @versioned_ace_files = grep { /.+\.ace\.\d+$/ } @ace_files;
    @ace_files = grep { /.+\.ace$/ } @ace_files;
    @ace_files = (@ace_files,@versioned_ace_files);

    my $w = Workflow::Operation->create(
        name => 'detect merges',
        operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Assembly::DetectMerges'),
    );
    
    $w->parallel_by('ace_file');
    
    #$w->log_dir('/gscmnt/936/info/jschindl/MISCEL/wflogs');
    
    $w->validate;
    if(!$w->is_valid)
    {
        $self->error_message("There was an error while validating merge detection workflow.\n") and return;
    }
      
    my $result = Workflow::Simple::run_workflow_lsf(
        $w,
        'ace_file' => \@ace_files,        
    );
    
    unless($result)
    {
        # is this bad?
        foreach my $error (@Workflow::Simple::ERROR) {

            $self->error_message( join("\t", $error->dispatch_identifier(),
                                             $error->name(),
                                             $error->start_time(),
                                             $error->end_time(),
                                             $error->exit_code(),
                                            ) );

            $self->error_message($error->stdout());
            $self->error_message($error->stderr());

        }
        return;

    }
	return 1;  
}
1;
