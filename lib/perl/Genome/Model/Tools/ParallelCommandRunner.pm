package Genome::Model::Tools::ParallelCommandRunner;

use strict;
use warnings;

use Genome;
use Command;
use File::Basename;
use IO::File;
use Genome::Sys;

class Genome::Model::Tools::ParallelCommandRunner {
    is => ['Command'],
    has => [ 
            command_list => { is => 'Array', doc => 'Required List of commands to execute in parallel.' }, 
            log_path => { is => 'String', doc => 'Required path to write log files.' }, 
            log_file => { is => 'String', doc => 'Optional file target to write logs. If provided, all parallel processes will write to this file. If none is provided,each process will write to their own individual log file specified by "[log_path]/parallel_[timestamp].log"', is_optional=>1}, 
           ],
};

sub help_brief {
    "Use workflow to run commands in parallel.";
}

sub help_synopsis {
    return <<"EOS"
    genome-model add-reads postprocess-alignments merge-alignments maq --model-id 5 --ref-seq-id all_sequences
EOS
}

sub help_detail {
    return <<EOS 
This command is usually called as part of the add-reads process
EOS
}

sub create {
    my $self = shift->SUPER::create(@_);
    
    return $self;
}

sub execute {
    	my $self = shift;
        my @commands = @{$self->command_list}; 
        #print("Running commands: ".join(",",@commands)."\n");
        #print('There are '.scalar(@commands).' to run.');
    
        Genome::Sys->create_directory($self->log_path);
        unless (-d $self->log_path) {
	   #print("Can't create directory for logging: ".$self->log_path);
	   $self->error_message("Can't create directory for logging: ".$self->log_path);
           die();
        }

	require Workflow::Simple;

    my $m = Workflow::Model->create(
        name => 'Container',
        input_properties => ['command_list','log_file','log_path'],
        output_properties => ['result','output_file']
    );
        
	my $op = $m->add_operation(
            name => 'Parallel Command.',
            operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::ParallelCommand')
    );

	$op->parallel_by('command_list');

    foreach my $inprop (@{ $m->operation_type->input_properties }) {
        $m->add_link(
            left_operation => $m->get_input_connector,
            left_property => $inprop,
            right_operation => $op,
            right_property => $inprop
        );
    }

    foreach my $outprop (@{ $m->operation_type->output_properties }) {
        $m->add_link(
            left_operation => $op,
            left_property => $outprop,
            right_operation => $m->get_output_connector,
            right_property => $outprop
        );
    }

        my $output = Workflow::Simple::run_workflow_lsf(
            $m,
            'command_list' =>\@commands, 
            'log_file' =>$self->log_file || '', 
            'log_path' =>$self->log_path || '', 
        );
 
        #$self->status_message("Output: ".$output);
        #print("Output: ".$output);

        if (!defined $output) {
            #error occurred.

            foreach my $error (@Workflow::Simple::ERROR) {
                $self->error_message($error->error);
            }
            die 'Commands with errors: ' . scalar @Workflow::Simple::ERROR;
        }

        #while ( my ($key, $value) = each(%$output) ) {
        #        if ( $key eq 'result' ) {	
        #        	#print "Output key: $key \n";
	#		my @result_array = @{$value};
        #        	for my $result_item (@result_array) {
        #            		print "Result: $result_item";
        #        	}
        #	}
	#} 


	#try printing the output
	#for my $output_item (keys %output) {
        #        my $output_string = $output{$output_item};
        #        $self->status_message("Output key: $output_string, Output value: $output_string");
        #        print("Output key: $output_string, Output value: $output_string");
        #}

return 1;

}


1;
