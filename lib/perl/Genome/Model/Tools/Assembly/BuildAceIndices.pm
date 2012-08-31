package Genome::Model::Tools::Assembly::BuildAceIndices; 

use strict;
use warnings;

use Genome;
use IO::File;
use Workflow::Simple;


class Genome::Model::Tools::Assembly::BuildAceIndices {
    is => 'Command',
    has => [
        ace_directory => {
            is => 'String',
            shell_args_position => 1,
            is_optional => 1,
            doc => 'the directory containing the ace files (usually edit_dir), if no ace files are specified, then all files ending in ace are indexed',
        },
        ace_file => {
            is => 'String',
            shell_args_position => 3,
            is_optional => 1,
            is_input => 1,
            doc => 'the ace file(s) that we are indexing, if the ace_directory is specified, then it will be assumed that the ace file is located in the ace_directory unless a full path is specified for each ace file',      
        },
        cache_dir => {
            is => 'String',
            shell_args_position => 2,
            is_optional => 1,
            doc => 'use this option to specify a cache dir that the indices will use, if specific ace files are given, then only those ace files will the specified cache, otherwise all ace files in the ace_directory will share the same cache',
        }
    ],    
    doc => 'Indexes ace files'
};

sub check_multiple_versions
{
    my ($self, $ace_directory) = @_;
    
    #filter out any non ace files
    
    my @ace_files = `ls $ace_directory/*.ace*`;
    my @versioned_ace_files;
    @versioned_ace_files = grep { /.+\.ace\.\d+$/ } @ace_files;
    @ace_files = grep { /.+\.ace$/ } @ace_files;
    @ace_files = (@ace_files,@versioned_ace_files);

    return 1 if(scalar @versioned_ace_files);
    
    return 0;
}

sub execute
{
    my $self = shift;
    my $ace_directory = $self->ace_directory; 
    my $ace_file = $self->ace_file;
    my $cache_dir = $self->cache_dir;
    my @ace_files;
    @ace_files = split /,/, $ace_file if (defined $ace_file);
    if($self->ace_directory||($#ace_files > 1))
    {
        chdir($ace_directory);
        $self->warning_message("Versioned ace files detected.  The merging toolkit only works on files ending in .*.ace, other ace files will be ignored\n") if($self->check_multiple_versions($ace_directory));
        @ace_files = `ls $ace_directory/*.ace`;
    
        $self->error_message( "There are no valid ace files in $ace_directory\n") and return unless (scalar @ace_files);  
        chomp @ace_files;
        
        if($ace_file)#if the user specified a subset of ace files, then we will use those instead
        {
            @ace_files = split /,/,$ace_file;
        }
        my @ace_files_to_index;
        @ace_files_to_index = grep { ! -d $_.'.idx' } @ace_files;
        print "number of ace files is ".scalar(@ace_files_to_index)."\n";
        
        if(scalar(@ace_files_to_index))
        {
            my $w = Workflow::Operation->create(
                name => 'build indices',
                operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Assembly::BuildAceIndices'),
            );

            $w->parallel_by('ace_file');

            $w->log_dir('/gscmnt/936/info/jschindl/MISCEL/wflogs');

            $w->validate;
            if(!$w->is_valid)
            {
                $self->error_message("There was an error while validating parallel merge workflow.\n") and return;
            }


            my $result = Workflow::Simple::run_workflow_lsf(
                $w,
                'ace_file' => \@ace_files_to_index,        

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
    }
    elsif($ace_file && -e $ace_file)
    {
        my $ao = Genome::Model::Tools::Pcap::Ace->new(input_file => $ace_file);
        $self->error_message("There was a problem indexing $ace_file.\n") and return unless(defined $ao);
    }
    
    if($cache_dir)
    {
        $self->error_message("Cache directory already exists, caches can only be created once, please remove the old one and try again.\n") and return if(-d $cache_dir);
        my $ao = Genome::Model::Tools::Pcap::Ace->new(input_files => \@ace_files, cache_dir => $cache_dir);        
    }    

    return 1;
}

1;
