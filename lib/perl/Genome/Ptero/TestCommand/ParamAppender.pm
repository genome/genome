package Genome::Ptero::TestCommand::ParamAppender;

use strict;
use warnings;
use Genome;

class Genome::Ptero::TestCommand::ParamAppender {
    is => "Command::V2",
    has_input => [
        prefix => {
            is => "Text",
        },
        suffix => {
            is => "Text",
        },
    ],
    has_output => [
        output => {
            is => "Text",
        },
    ],
    has => [
        lsf_resource => {
            is => "Text",
            default => "-R 'select[gtmp>1 && mem>3000] span[hosts=1] rusage[gtmp=1 && mem=3000]' -n 1 -M 3000000",
        }
    ],
};

sub shortcut {
    my $self = shift;
    return $self->execute(@_);
}

sub execute {
    my $self = shift;
    my $p = $self->prefix;
    my $s = $self->suffix;
    my $a = "$p$s";
    print "$p + $s = $a!!\n";
    $self->output($a);
    return 1;
};

1;
