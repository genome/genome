package Genome::Ptero::TestCommand::NullCommand;

use strict;
use warnings;
use Genome;

use File::Path qw(make_path);
use File::Basename qw(dirname);

class Genome::Ptero::TestCommand::NullCommand {
    is => "Command::V2",
    has_input => [
        res1 => {
            is => "Text",
            is_optional => 1,
        },
        res2 => {
            is => "Text",
            is_optional => 1,
        },
        catcher => {
            is => "Text",
            is_optional => 1,
        },
        param => {
            is => "Text",
            doc => "A number",
            default_value => "x",
        },
    ],
    has_output => [
        animal => {
            is => "Text",
        },
    ],
    has => [
        lsf_resource => {
            default_value => "-M 200000 -n 4 -R 'rusage[mem=200:gtmp=5]'",
        },
        lsf_queue => {
            default_value => Genome::Config::get('lsf_queue_short'),
        },
    ]
};

sub execute {
    my $self = shift;
    print "param: " . $self->param . "\n";
    print "catcher: " . $self->catcher . "\n";
    $self->animal("frog");
    return 1;
}

1;
