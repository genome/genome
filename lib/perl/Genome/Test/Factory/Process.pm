package Genome::Test::Factory::Process;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;

our @required_params = qw(subclass_name);

sub generate_obj {
    my $self = shift;

    return Genome::Process->create(@_);
}

sub create_subclass_name {
    return 'Genome::Test::Factory::Process::TestProcess';
}

package Genome::Test::Factory::Process::TestProcess;

class Genome::Test::Factory::Process::TestProcess {
    is => 'Genome::Process',
};

1;
