
package Genome::Model::Tools::Example1;     # rename this when you give the module file a different name <--

use strict;
use warnings;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Example1 {
    is => 'Command',                       
    has => [                                # specify the command's single-value properties (parameters) <--- 
        foo      => { is => 'Text',       doc => "some foozy thing" },
        bar      => { is => 'Boolean',    doc => "a flag to turn on and off", is_optional => 1 },
        baz      => { is => 'My::Object', doc => "this won't appear as a command-line option since it's not a primative", is_optional => 1 },
    ], 
    has_many => [                           # specify the command's multi-value properties (parameters) <--- 
        infiles  => { is => 'Text', doc => 'this is a list of values' },
        outfiles => { is => 'Text', doc => 'also a list of values' },
        bare_args => { is => 'Text', shell_args_position => 1, is_optional => 1, doc => 'list of bare arguments' },
    ], 
};

sub _is_hidden_in_docs { $ENV{GENOME_EXAMPLES} ? () : 1 }

sub sub_command_sort_position { -2 }

sub help_brief {                            # keep this to just a few words <---
    "WRITE A ONE-LINE DESCRIPTION HERE"                 
}

sub help_synopsis {                         # replace the text below with real examples <---
    return <<EOS
gmt example1 --bar
gmt example1 --foo=hello
gmt example1 --foo=goodbye --bar
gmt example1 --foo=hello barearg1 barearg2 barearg3
gmt example1 --infiles=a,b,c --outfiles=d,e,f 

EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
This is a dummy command.  Copy, paste and modify the module! 
CHANGE THIS BLOCK OF TEXT IN THE MODULE TO CHANGE THE HELP OUTPUT.
EOS
}

#sub create {                               # rarely implemented.  Initialize things before execute.  Delete unless you use it. <---
#    my $class = shift;
#    my %params = @_;
#    my $self = $class->SUPER::create(%params);
#    # ..do initialization here
#    return $self;
#}

sub execute {                               # replace with real execution logic.
    my $self = shift;
    print "Running example command:\n" 
        . "    foo is " . (defined $self->foo ? $self->foo : '<not defined>')
        . "\n" 
        . "    bar is " . (defined $self->bar ? $self->bar : '<not defined>') 
        . "\n" 
        . "    baz is " . (defined $self->baz ? $self->bar : '<not defined>') 
        . "\n" 
        . "    infiles are :" . join(',',$self->infiles) 
        . "\n" 
        . "    outfiles are :" . join(',',$self->outfiles) 
        . "\n" 
        . "    bare args are: " . ($self->bare_args ? join(",",$self->bare_args)  : '<none>') 
        . "\n";     
    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

1;

