
package Genome::Model::Tools::Example2;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Example2 {
    is => 'Command',                       
    has => [
        in  => { shell_args_position => 1 },
        out => { shell_args_position => 2 },
    ], 
};

sub _is_hidden_in_docs { $ENV{GENOME_EXAMPLES} ? () : 1 }

sub sub_command_sort_position { -1 }

sub help_brief {                            # keep this to just a few words <---
    "WRITE A ONE-LINE DESCRIPTION HERE"                 
}

sub help_synopsis {                         # replace the text below with real examples <---
    return <<EOS
gmt example2 IN OUT
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
This is a dummy command.  Copy, paste and modify the module! 
CHANGE THIS BLOCK OF TEXT IN THE MODULE TO CHANGE THE HELP OUTPUT.
EOS
}

sub execute {                               # replace with real execution logic.
    my $self = shift;
    print "Running example command:\n" 
        . "    in is " . (defined $self->in ? $self->in : '<not defined>')
        . "\n" 
        . "    out is " . (defined $self->out ? $self->out : '<not defined>') 
        . "\n";     
    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

1;

