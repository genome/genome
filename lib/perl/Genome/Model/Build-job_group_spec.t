use strict;
use warnings;

use Test::More tests => 5;

use above 'Genome';
use Genome::Model::Build qw();
use Genome::Sys qw();

for my $subname (qw(_default_job_group _job_group_spec)) {
    my $fullname = 'Genome::Model::Build::' . $subname;
    no strict 'refs';
    *{"main::$subname"} = \&{$fullname};
};

do {
    is(_job_group_spec({ job_group => undef }), '');
    is(_job_group_spec({ job_group => '' }), '');
    is(_job_group_spec({ job_group => 'foo' }), ' -g foo');
    my $username = getpwuid($<);
    is(_job_group_spec({}), ' -g ' . _default_job_group());
    is(_job_group_spec(), ' -g ' . _default_job_group());
};
