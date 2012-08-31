package Genome::Model::Tools::ContaminationScreen::Hmp;

use strict;
use warnings;

use Genome;    

class Genome::Model::Tools::ContaminationScreen::Hmp {
    is => 'Genome::Model::Tools::ContaminationScreen',
    is_abstract => 1,
    has => [
            filter_list=> {
                             doc => 'file to store list of id\'s to filter from original fasta',
                             is => 'String',
                             is_output => 1,
                             is_optional => 1,
            },
    ],
};

sub help_brief {
    'Tools for human screening.'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools hmp...
EOS
}

sub xhelp_detail {                           
    return <<EOS 
EOS
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return $self;
}

sub execute
{
    die("implement execute in inheriting class");
}

1;
