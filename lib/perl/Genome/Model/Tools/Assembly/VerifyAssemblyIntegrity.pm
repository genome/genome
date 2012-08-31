package Genome::Model::Tools::Assembly::VerifyAssemblyIntegrity;

use strict;
use warnings;

use Genome;
use File::Basename;
use IO::String;

use Data::Dumper;
use Workflow::Simple;

class Genome::Model::Tools::Assembly::VerifyAssemblyIntegrity
{
    is => 'Command',
    has => 
    [        
        assembly_dir =>
        {
            type => 'String',
            is_optional => 1,
            doc => "The directory containing the ace files to verify",
            #default_value => '/gscmnt/936/info/jschindl/MISCEL/dir',
        },
        ace_file_extension =>
        {
            type => 'String',
            is_optional => 1,
            doc => "The extension of the ace files, default is 'ace.1', i.e. Trichinella_spiralis-3.0_070112.pcap.scaffold0.ace.1.  This allows the user to verify only a single version of an assembly dir contianing multiple versions."        
        },
        all_versions =>
        {
            type => 'Integer',
            is_optional => 1,
            default_value => 0,
            doc => "If this flag is set to 1, then all files containing the *.ace* extension will be evalulated.  This means that all versions, and any random ace files will be evaluated at the same time.  This is turned off by default, as it probably is not what you want"
        },       
        verbose =>
        {
            type => 'Integer',
            is_optional => 1,
            default_value => 0,
            doc => "If this flag is set to 1, then ace files with errors are listed, if it is set to 2 then status messages are also printed out for valid ace files"
        }
        
    ]
};

sub help_brief {
    "This tool verifies the ace files in a pcap assembly"
}   

sub help_synopsis { 
    return;
}
sub help_detail {
    return <<EOS 
    Verify the integrity of an ace file
EOS
}




############################################################
sub execute { 
    my $self = shift;
    $Storable::eval=1;
    my $verbose = $self->verbose;
    my $dir = $self->assembly_dir || `pwd`;
    chomp $dir;
    $self->error_message("Assembly dir $dir does not exist") and return unless (-d $dir);
    my $ace_file_extension = $self->ace_file_extension || 'ace.1';
    $ace_file_extension = '.ace*' if $self->all_versions;
    my @ace_files = `ls $dir/*$ace_file_extension`;
    chomp @ace_files;
    
    $self->error_message("No ace with extension $ace_file_extension found in assembly dir $dir\n") and return unless (scalar @ace_files);
        
    my $w = Workflow::Operation->create(
        name => 'verify ace files',
        operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Assembly::VerifyAceFileIntegrity'),
    );
    
    $w->parallel_by('ace_file');
    
#    $w->log_dir('/gsc/var/log/verifyassemblyintegrity');
    
    $w->validate;
    $w->is_valid;
        
    my $result = run_workflow_lsf(
        $w,
        'ace_file' => \@ace_files,
    );
    
    if($result)
    {
        my $error = 1;    
        for(my $i = 0; $i <= $#ace_files; $i++) {
        
            if(!$result->{result}->[$i] && $verbose >0 || 
               $result->{result}->[$i]==1 && $verbose >1)
            {
                print "$ace_files[$i] result " . $result->{result}->[$i] . " output: " . $result->{output}->[$i];
            }
            $error = 0 if(!$result->{result}->[$i]);
        }
        $self->error_message("There were errors parsing the ace files.\n") and return if(!$error);
    }
    else 
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
    
    print "Assembly verified successfully\n";    

    return 1;
    

}



1;
