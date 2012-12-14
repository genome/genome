#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 7;

use_ok('Genome::Model::Tools::Picard');

my $picard = Genome::Model::Tools::Picard->create();

isa_ok($picard, 'Genome::Model::Tools::Picard');

# should get default versions since we do not specify
my $picard_version = $picard->use_version();
ok(-e $picard->path_for_picard_version($picard_version), "picard version ($picard_version) exists");


my $newest_picard_version = $picard->latest_version();
ok(-e $picard->path_for_picard_version($newest_picard_version), "newest picard version ($newest_picard_version) exists");

## email test
SKIP: {
    skip 'monitor_shellcmd test can be annoying', 3 unless($ENV{'MONITOR_SHELLCMD_TEST'});

    my $test_subcmd;

    ok($test_subcmd = File::Temp->new(SUFFIX => ".pl"), 'opening temp pl');
    $test_subcmd->autoflush(1);
    ok($test_subcmd->print(q|
#!/usr/bin/env genome-perl

use strict;
use warnings;

use IO::Handle;

STDOUT->autoflush(1);

for (0..6) {
    sleep 1;

    if ($_ && $_ % 5 == 0) {
        sleep 7;
    }
    print $_ . "\n";
}
|),'writing temp pl');

    my $picard_cmd = Genome::Model::Tools::Picard->create(
        _monitor_check_interval => 1,
        _monitor_stdout_interval => 5,
        _monitor_mail_to => Genome::Sys->username,
    );

    my $rv;
    ok($rv = $picard_cmd->monitor_shellcmd({
        cmd => 'perl ' . $test_subcmd->filename
    }),'run temp pl');

};
