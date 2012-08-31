package Genome::Sys::Command::Completion::Bash;

local $ENV{PERL_ABOVE_QUIET} = 1;
use Genome;
use strict;
use warnings;

class Genome::Sys::Command::Completion::Bash {
    is => 'Genome::Command::Base',
    doc => 'outputs completion code for Bash',
};

sub help_detail {
    my $help_detail;

    $help_detail .= "genome and gmt come with support for command line tab completion in bash. To enable it you must add the output of 'genome sys completion bash' to your ~/.profile, e.g.:\n\n";
    $help_detail .= "\$ genome sys completion bash >> ~/.profile\n\n";
    $help_detail .= "Alternatively, you can use the result of the completion command directly by adding the following line to your ~/.profile:\n\n";
    $help_detail .= "eval \"\`genome sys completion bash 2> /dev/null\`\"";

    return $help_detail;
}

sub execute {
    my $self = shift;
    
    my $text;

    $text .= "\n";
    $text .= "### Begin Bash completion for genome and gmt ###\n";
    $text .= "source /gsc/scripts/opt/genome/current/usr/bin/getopt_complete.sh\n";
    $text .= "complete -F _getopt_complete genome\n";
    $text .= "complete -F _getopt_complete gmt\n";
    $text .= "### End Bash completion for genome and gmt ###\n";
    $text .= "\n";

    print $text;
        
}
