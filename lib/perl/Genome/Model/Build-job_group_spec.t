use strict;
use warnings;

use Test::More tests => 5;

use Genome::Model::Build qw();
use Genome::Sys qw();

sub _job_group_spec {
    return Genome::Model::Build::_job_group_spec(@_);
};

do {
    is(_job_group_spec({ job_group => undef }), '');
    is(_job_group_spec({ job_group => '' }), '');
    is(_job_group_spec({ job_group => 'foo' }), ' -g foo');
    my $username = getpwuid($<);
    is(_job_group_spec({}), " -g /apipe-build/$username");
    is(_job_group_spec(), " -g /apipe-build/$username");
};
