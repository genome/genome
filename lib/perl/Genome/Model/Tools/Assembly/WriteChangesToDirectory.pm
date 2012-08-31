package Genome::Model::Tools::Assembly::WriteChangesToDirectory; 

use strict;
use warnings;

use Genome;
use FindBin;
use Carp::Assert;
use Carp;
use Cwd;
use Utility;
use Workflow::Simple;


class Genome::Model::Tools::Assembly::WriteChangesToDirectory {
    is => 'Command',
    has => [
        cache_dir => {
            is => 'String',
            shell_args_position => 1,
            is_optional => 0,
            is_input => 1,
            doc => 'this is the cache dir that represents the assembly in it\'s current state',
        },
        output_directory => {
            is => 'String',
            shell_args_position => 2,
            is_optional => 1,
            is_input => 1,
            doc => 'this is the output directory that we are writing the assembly out to',
        },
        number =>  {
            is => 'Number',
            shell_args_position => 3,
            is_optional => 1,
            is_input => 1,
            doc => 'this is the number of ace files to write the cached assembly out to, defaults to 1 ace file',
        },
        prefix => {
            is => 'String',
            shell_args_position => 4,
            is_optional => 1,
            is_input => 1,
            doc => 'this is the prefix used to form the base name of the ace file (i.e. a prefix of \'scaffold\' would produce ace files name scaffold0.ace, scaffold1.ace, etc.)',
        },
        index =>  {
            is => 'Number',
            shell_args_position => 5,
            is_optional => 1,
            is_input => 1,
            doc => 'this is a private param that is used by workflow',
        },
    ],
    
    doc => 'writes an updated ace file using an input ace file and the changes stored in the mysql database (this is part of the 4th stage in the Msi pipeline)'
};


sub execute {
    my $self = shift;

    my $index = $self->index;
    my $number = $self->number || 1;    
    my $prefix = $self->prefix||'';
    my $output_directory = $self->output_directory;
    my $cache_dir = $self->cache_dir;
    `mkdir -p $output_directory` if(!-d $output_directory );

    if(defined $index)
    {   
        my $ao = Genome::Model::Tools::Pcap::Ace->new(cache_dir => $cache_dir);
        $ao->write_file(output_file => "$output_directory/$prefix$index.ace", number => $number,index => $index, number => $number);
    }
    else
    {    
        my $w = Workflow::Operation->create(
            name => 'write changes to directory',
            operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Assembly::WriteChangesToDirectory'),
        );

        $w->parallel_by('index');

        #$w->log_dir('/gscmnt/936/info/jschindl/MISCEL/wflogs');

        $w->validate;
        if(!$w->is_valid)
        {
            $self->error_message("There was an error while validating merge detection workflow.\n") and return;
        }
        my @index;
        @index = (0..($number-1));
        my $result = Workflow::Simple::run_workflow_lsf(
            $w,
            index => \@index,
            cache_dir => $cache_dir,
            output_directory => $output_directory,
            number => $number,
            prefix => $prefix,
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
    }

    return 1;
}

1;
