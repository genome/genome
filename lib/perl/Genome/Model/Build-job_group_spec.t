use strict;
use warnings;

use Test::More tests => 5;

use above 'Genome';
use Genome::Model::Build qw();
use Genome::Sys qw();

for my $subname (qw(_default_job_group _job_group)) {
    my $fullname = 'Genome::Model::Build::' . $subname;
    no strict 'refs';
    *{"main::$subname"} = \&{$fullname};
};

do {
    ok(!_job_group({ job_group => undef }));
    ok(!_job_group({ job_group => '' }));
    is(_job_group({ job_group => 'foo' }), 'foo');
    my $username = getpwuid($<);
    is(_job_group({}), _default_job_group());
    is(_job_group(), _default_job_group());
};
