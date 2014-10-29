use strict;
use warnings;

use Test::More tests => 4;
use Test::Deep qw(cmp_bag);

use Genome;

UR::Object::Type->define(
    class_name => 'Sports::Player',
    has => [
        name => { is => 'Text' },
    ],
    has_optional => [
        team_id => { is => 'Text' },
        team => {
            is => 'Sports::Team',
            id_by => 'team_id',
        },
    ],
);

UR::Object::Type->define(
    class_name => 'Sports::Team',
    has => [
        name => {
            is => 'Text',
        },
    ],
    has_optional => [
        players => {
            is => 'Sports::Player',
            is_many => 1,
            reverse_as => 'team',
        },
    ],
);

UR::Object::Type->define(
    class_name => 'Sports::Team::Command::Copy',
    is => 'Genome::Command::Copy',
);

my $lakers = Sports::Team->create(name => 'Lakers');
my $mj = Sports::Player->create(team_id => $lakers->id, name => 'Magic Johnson');

{
    my $cmd = Sports::Team::Command::Copy->create(
        source => $lakers,
        changes => ['name+=II'],
    );
    ok(!$cmd->execute, 'command failed when trying to += Text');
};


Sports::Team::Command::Copy->dump_status_messages(0);
Sports::Team::Command::Copy->execute(
    source => $lakers,
    changes => ['name.=II'],
);
my $lakersII = Sports::Team->get(name => $lakers->name . 'II');
ok($lakersII, 'copied team with appended name');

cmp_bag([$lakersII->players], [], 'lakersII have no players');
cmp_bag([$lakers->players], [$mj], 'lakers have correct players');
