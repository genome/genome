package Genome::Model::Tools::Analysis::AddAndBuild;

use strict;
use warnings;

use Genome;
use Command;
use File::Slurp;
use Text::CSV_XS;
use IPC::Run;

class Genome::Model::Tools::Analysis::AddAndBuild {
    is => 'Command',
    has => [
    model_ids => { 
        type => 'String',
        is_optional => 1,
        doc => "string of space separated model ids to check",
    },
    model_group => {
        type => 'String',
        is_optional => 1,
        doc => "The name of the Model Group to use. Models where is_default is set will be ignored",
    },
    maximum_allowed_error => {
        type => 'Float',
        is_optional => 1,
        doc => "The maximum allowed gerald error rate to assign to a model",
        default => 3.0,
    },
    add => {
        type => 'Boolean',
        is_optional => 1,
        doc => "whether or not to add instrument data to the models or just report what the tool finds (--noadd). Adds by default.",
        default => 1,
    },
    build => {
        type => 'Boolean',
        is_optional => 1,
        doc => "whether or not to build the models once instrument data is added. (implies --add). Off by default.",
        default => 0,
    },

    ]
};


sub execute {
    my $self=shift;
    my @models;
    if($self->model_ids) {
        my @model_ids = split / /, $self->model_ids;
        for my $model_id (@model_ids) {
            my $model = Genome::Model->get($model_id);
            unless(defined($model)) {
                $self->error_message("Unable to find model $model_id");
                return;
            }
            push @models, $model;
        }
    }
    elsif($self->model_group) {
        my $group = Genome::ModelGroup->get(name => $self->model_group);
        unless($group) {
            $self->error_message("Unable to find a model group named " . $self->model_group);
            return;
        }
        push @models, $group->models;
    }
    else {
        $self->error_message("You must provide either model id(s) or a model group name to run this script");
        return;
    }

    if($self->build) {
        $self->add(1);
    }

    #set up csv parser to use for parsing output of instrument data
    my $csv = Text::CSV_XS->new();
    unless($csv) {
        $self->error_message("Couldn't create CSV parser");
        return;
    }

    foreach my $model (@models) {
        my $model_id = $model->id;
        my $model_name = $model->name;
        #Assume that is_default indicates the model is complete/in data freeze
        if($model->is_default) {
            $self->status_message("Skipping $model_name ($model_id) under the assumption that it is complete");
            next;
        }
        #add new data
        $self->status_message("Querying unassigned data for $model_name ($model_id)");

        #TODO need better error handling here
        #Using filt_error_rate_avg since this is properly populated for Standard(fragment) runs and identical to rev_filt_error_rate_avg otherwise
        my @lines = `genome model instrument-data list --model-id $model_id --unassigned --style=csv --show id,flow_cell_id,subset_name,library_name,filt_error_rate_avg,fwd_filt_error_rate_avg --noheaders`; 

        unless(@lines) {
            $self->error_message("Unable to grab unassigned instrument data for $model_name ($model_id) due to error or no available data");
            next;
        }

        #the following ganked from John Osborne's script
        my (@allin, @forward, @reverse);
        foreach my $line (@lines) {
            $csv->parse($line);

            my ($instdata_id, $flowcell, $lane, $libname, $reverr, $fwderr) = $csv->fields();

            #add as all if fwderr is null. Assuming fragment run if this is the case
            if(($reverr < $self->maximum_allowed_error) && ($fwderr eq '<NULL>' || ($fwderr < $self->maximum_allowed_error))) {
                push(@allin, $instdata_id);
            }
            elsif($reverr < $self->maximum_allowed_error) {
                push(@reverse, $instdata_id);
                $self->status_message("Excluding forward read of $flowcell lane $lane with id $instdata_id due to error rate of $fwderr%"); 
            }
            elsif($fwderr ne '<NULL>' && $fwderr < $self->maximum_allowed_error) {
                push(@forward, $instdata_id);
                $self->status_message("Excluding reverse read of $flowcell lane $lane with id $instdata_id due to error rate of $reverr%"); 
            }
            else {
                my $error_report_string = $fwderr eq '<NULL>' ? "($reverr%)" : "($fwderr%, $reverr%)";
                $self->status_message("Excluding $flowcell lane $lane with id $instdata_id due to error rate $error_report_string");
            }
        }

        #report how many lanes of each type were found
        $self->status_message(sprintf("%d all good\n",scalar(@allin)));
        $self->status_message(sprintf("%d forward only\n",scalar(@forward)));
        $self->status_message(sprintf("%d reverse only\n",scalar(@reverse)));

        #skip if no data passed our filter  
        unless(@allin || @forward || @reverse) {
            $self->status_message("No suitable data found for $model_name ($model_id)");
            next;
        }

        #add in the new data
        if($self->add) {    
            my @command = ('genome','model','instrument-data','assign', '--model-id',$model_id);
            my @command_all = (@command, '--instrument-data-ids="'.join(' ',@allin).'"');
            my @command_fwd = (@command, '--instrument-data-ids="'.join(' ',@forward).'"','--filter=forward-only');
            my @command_rev = (@command, '--instrument-data-ids="'.join(' ',@reverse).'"','--filter=reverse-only');

            if(@allin) {
                print join(' ', @command_all), "\n";
                my $cmd1 = join(' ', @command_all);
                system($cmd1);
            }
            if(@forward) {
                print join(' ', @command_fwd), "\n";
                my $cmd2 = join(' ', @command_fwd);
                system($cmd2);
            }
            if(@reverse) {
                print join(' ', @command_rev), "\n";
                my $cmd3 = join(' ', @command_rev);
                system($cmd3);
            }

            
            #build since we added new data
            if($self->build) {
                my $build_command = ("genome model build start $model_id");
                my $rv = system($build_command);
                if($rv == -1) {
                    $self->error_message("Execution of $build_command failed: $!");
                }
                elsif($rv > 0) {
                    $self->error_message(sprintf("Execution of $build_command returned error code %d", $rv>>8));
                }
                else {
                    $self->status_message("$build_command successful");
                }
            }
        }
    }

    return 1;

}


1;

sub help_brief {
    "Iterates through a model group, finds lanes of data with an acceptable error rate and adds them to the models. It will optionally build them."
}

sub help_detail {
    <<'HELP';
This tool take either a list of model ids or a model-group and iterates over each one. Models on which the is_default attribute is set will be assumed to be in "data-freeze" and will not be touched. Adding unassigned instrument data and excluding lanes or parts of lanes that have GERALD calculated error rates above or equal to the set maximum. By default, the maximum error is set to 3. You may also use this tool to merely query out the unassigned instrument data for your model-group by using the --noadd option. You can have it kickoff builds if it finds new data by specifying --build. Logging is pretty extensive. Prepare to be bombarded with information if you have a lot of models.  
HELP
}
