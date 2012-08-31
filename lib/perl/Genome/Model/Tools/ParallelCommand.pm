package Genome::Model::Tools::ParallelCommand;

use strict;
use warnings;

use Command;
use File::Basename;
use IO::File;

class Genome::Model::Tools::ParallelCommand {
    is => ['Command'],
    has_input => [
        command_list => {
            is => 'String',
            doc => 'The list of commands to execute in parallel.' 
        },
        log_path => { is => 'String', doc => 'Required path to write log files.' }, 
        log_file => { is => 'String', doc => 'Optional file target to write logs. If provided, all parallel processes will write to this file. If none is provided,each process will write to their own individual log file specified by "[log_path]/parallel_[process id].log"', is_optional=>1}, 
    ],
    has_output => [
        output_file => { 
            is => 'String', 
            is_optional => 1, 
        }
    ],
    has_param => [
        lsf_resource => {
            default_value => 'select[type==LINUX64]rusage[mem=16000]',
        }
    ],
};

sub create {
    my $self = shift->SUPER::create(@_);
    return $self;
}

sub status_message {
   my $self=shift;
   my $msg=shift;
   #print $msg . "\n";
}

sub execute {
    my $self=shift;
    my $pid = getppid();
    my $log_target;

    if (defined($self->log_file) && $self->log_file ne '' ) {
	$log_target = $self->log_path."/".$self->log_file;
    } else {
	$log_target = $self->log_path."/parallel_$pid.log";
    } 
    #open(STDOUT, ">>$log_target") || die "Can't redirect stdout.";
    #open(STDERR, ">&STDOUT"); 
    $self->status_message("Log target is: ".$log_target);  
    my @list;
    if ( ref($self->command_list) ne 'ARRAY' ) {
       push @list, $self->command_list; 		
    } else {
       @list = @{$self->command_list};   	#the parallelized code will only receive a list of one item. 
    }

    my $now = UR::Time->now;
    my $rv;
    $self->status_message("*** Starting parallel command at $now. ***");
    for my $list_item ( @list  ) {
      	#my $cmd = ${$list_item};    		
        $self->status_message("parallel command: ".$list_item);
       	$rv = `$list_item`;
        #$self->status_message($rv);
       	#$rv = system($list_item);
   	#if($rv) {
       	#	$self->error_message("problem running $list_item");
       	#	return;
    	#}
    }
    $now = UR::Time->now;
    $self->status_message("*** Parallel command completed at $now. ***");
         
    #return 1;
    return $rv;

} #end execute

1;
