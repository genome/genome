package Genome::Model::Tools::Assembly::VerifyAceFileIntegrity;

use strict;
use warnings;

use Workflow;
use Genome;
use File::Basename;
use IO::String;
use IPC::Open3;

class Genome::Model::Tools::Assembly::VerifyAceFileIntegrity
{
    is => 'Command',
    has => 
    [        
        ace_file =>
        {
            type => 'String',
            is_optional => 0,
            doc => "the name of the ace file to verify",
            is_input => 1,
        },
        verbose =>
        {
            type => 'Integer',
            is_optional => 1,
            default_value => 0,
            doc => "If this flag is set to 1, then consed's error messaging is printed, if it is set to 2, then all of consed's output is printed while parsing",
            is_input => 1,
        },
        output => {
            doc => 'output',
            is_output => 1,  
            is_optional => 1,     
        }
        
    ]
};

sub help_brief {
    "Move Pooled BAC assembly into separate projects"
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
    my $ace_file = $self->ace_file;
    my $verbose = $self->verbose;
    $self->error_message("Ace file:$ace_file does not exist!") and return unless(-e $ace_file);
    print "Verifying $ace_file ..\n" if ($verbose >=0);
    # First, determine what dir we should look in
    my $p = __PACKAGE__ . '.pm';
    $p =~ s/::/\//g;
    my $loaded_dir = dirname($INC{$p});

    my $c_command = "$loaded_dir/acecheck/bin/acecheck";
    my $ld_path_set = "LD_LIBRARY_PATH=$loaded_dir/acecheck/lib:\$LD_LIBRARY_PATH";

    
    $c_command = "$ld_path_set $c_command $ace_file";
    my $pid = open3(undef, \*RDSTDOUT, \*RDSTDERR,$c_command); 
    waitpid( $pid, 0);
    
    my $rc = $? >> 8;
    if($rc)
    {
        if(-e 'core')
        {
            print "core file detected after consed crashed, removing...\n" if ($verbose >0);
            `/bin/rm ./core`;
        }    
    }
    my @out = <RDSTDOUT>;
    my @err = <RDSTDERR>;
    print @out if($verbose&&$verbose>1);
    print @err if($verbose>0 );
    if(grep {$_=~/parsed correctly/} @out)
    {
        $self->output("$ace_file is valid.\n");
        print $self->output if($verbose >=0);
        return 1;
    }
    else
    {
        my $line;
        foreach (@err)
        {
            if(/line \d+/)
            {
                ($line) = $_=~ /.+line (\d+)/;
                last;
            }
        }
        $self->output("$ace_file is invalid at line $line.\n");
        print $self->output if ($verbose >=0);
        return 0;
    }    

}



1;
