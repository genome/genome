package Genome::Sys::Command::Services::TaskRunner;

use warnings;
use strict;
use Genome;
use JSON::XS;
use File::Path;

class Genome::Sys::Command::Services::TaskRunner {
    is => 'Command',
    doc => 'Runner for scheduled tasks',
    has => [
        output_basedir => {
            is => 'String',
            doc => 'Directory to put stdout/stderr from run jobs',    
        },
        restart_file => {
            is => 'String',
            doc => 'File to watch for restart trigger - touch this file to restart daemon',
            is_optional => 1
        },
        _restart_file_mtime => {
            is => 'Integer',
            is_optional => 1,
            is_transient => 1
        }
    ],
};

sub execute {
    my $self = shift;
        
    if ($self->restart_file) {
        $self->_restart_file_mtime(0);
        
        if (-f $self->restart_file) {
            $self->_restart_file_mtime((stat($self->restart_file))[9]);
        }
    }

    while (1) {
        $self->task_loop();
        $self->check_restart_loop() if $self->restart_file;
        sleep 15;
    }

}


sub task_loop {
    my $self = shift;
    
    my $ds = $UR::Context::current->resolve_data_sources_for_class_meta_and_rule(Genome::Task->__meta__);
    my $dbh = $ds->get_default_dbh;

    # get lock on pending jobs and bypass object cache
    my $res = $dbh->selectall_arrayref("SELECT id FROM genome_task WHERE status='submitted' FOR UPDATE");
   
    if (@$res > 0) { 
        my ($id) = $res->[0]->[0];
         
        my $task = Genome::Task->get($id); 
        # before we fork, we'll drop the db connection and lose our lock.  
        # so first, update the status so that nobody else takes it in the meantime.
        $task->status("pending_execute");
        mkpath($self->output_basedir . "/" . $task->id);
        $task->stderr_pathname($self->output_basedir . "/" . $task->id . "/stderr");
        $task->stdout_pathname($self->output_basedir . "/" . $task->id . "/stdout");
        UR::Context->commit;
        
        my $pid = UR::Context::Process->fork();
        if (!$pid) {
           $self->status_message("Forked $$, running $id");
           my ($old_stdout, $old_stderr);
           open ($old_stdout, ">&STDOUT");
           open ($old_stderr, ">&STDERR");

           open (STDERR, ">>" . $task->stderr_pathname);
           open (STDOUT, ">>" . $task->stdout_pathname);
        
           my $ret;
           eval { 
               Genome::Task::Command::Run->create(task=>$task)->execute; 
               $ret = UR::Context->commit;
           };
           # rescue after a failure to sync
           if (!$ret || $@ ) {
               $self->error_message("Failed to execute: $@");
               UR::Context->rollback;
               $task = UR::Context->current->reload('Genome::Task', id=>$id);
               $task->out_of_band_attribute_update(status=>'failed', time_finished=>$UR::Context::current->now);
           }
           exit(0);
        } else {
           $task->unload; 
           waitpid($pid, 0);
        }
    } else {
        # need to commit to throw away the lock from select FOR UPDATE
        UR::Context->commit;
    }
}

sub check_restart_loop() {
    my $self = shift;

    if (-f $self->restart_file) {
        my $mtime = (stat($self->restart_file))[9];
        if ($mtime > $self->_restart_file_mtime) {
            $self->status_message("Restart file updated. Restarting as requested.");
            my $cmd = sprintf("genome sys services task-runner --outut-basedir=%s --restart-file=%s", $self->output_basedir, $self->restart_file);
            exec($cmd);
        }
    }
    
}


1;
