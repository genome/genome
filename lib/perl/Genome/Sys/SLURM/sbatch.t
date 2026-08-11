use strict;
use warnings;

use above 'Test::More';
use Test::Builder;
use Test::Fatal qw(exception);
use Genome::Sys::SLURM::sbatch qw();

if( Genome::Config::get('job_dispatch_backend') ne 'slurm' ) {
    plan skip_all => 'These tests only work on slurm';
} else {
    plan tests => 11;
}

my $fake_partition = 'fake_partition';
my @partitions = (map { Genome::Config::get($_) } qw(lsf_queue_build_worker lsf_queue_short));

ok( !Genome::Sys::SLURM::sbatch::_valid_slurm_partition($fake_partition),
    qq('$fake_partition' is not a valid partition));

ok( Genome::Sys::SLURM::sbatch::_valid_slurm_partition($partitions[0]),
    qq('$partitions[0]' is a valid partition'));

like( exception { Genome::Sys::SLURM::sbatch::_args(queue => $fake_partition, cmd => 'true') },
    qr/valid SLURM partition/,
    qq('$fake_partition' triggered exception));

is( exception { Genome::Sys::SLURM::sbatch::_args(queue => $partitions[0], cmd => 'true') },
    undef,
    qq('$partitions[0]' did not trigger exception));

do {
    my $_option_mapper = \&Genome::Sys::SLURM::sbatch::_option_mapper;
    no warnings 'redefine';
    local *Genome::Sys::SLURM::sbatch::_option_mapper = sub {
        my $option = $_option_mapper->(@_);
        $option =~ s/^--/-s:/;
        return $option;
    };
    is( Genome::Sys::SLURM::sbatch::_option_mapper('project'),
        '-s:comment',
        q('project' arg maps to '-s:comment' option with _option_mapper overridden));
};

is( Genome::Sys::SLURM::sbatch::_option_mapper('project'),
    '--comment',
    q('project' arg maps to '--comment' option));

my @cases = (
    [
        [
            email => 'nnutter@genome.wustl.edu',
            cmd => 'true',
        ], [qw(--mail-user nnutter@genome.wustl.edu --wrap true)], 'single option',
    ],[
        [
            email => 'nnutter@genome.wustl.edu',
            project => 'HighPriority',
            cmd => 'true',
        ], [qw(--mail-user nnutter@genome.wustl.edu --comment HighPriority --wrap true)], 'multiple options',
    ],[
        [
            hold_job => 0,
            cmd => 'true',
        ], [qw(--wrap true)], 'disabled flag',
    ],[
        [
            hold_job => 1,
            cmd => 'true',
        ], [qw(--hold --wrap true)], 'enabled flag',
    ],
    [
        [
            job_group => '/genome/test',
            user_group => 'compute-test',
            cmd => 'true',
        ], [qw(--account compute-test --wrap true)], 'job group ignored, user group mapped to account',
    ]
);
for my $case (@cases) {
    my @input = @{$case->[0]};
    my $expected = $case->[1];
    my $name = $case->[2];
    my $got = [Genome::Sys::SLURM::sbatch::args_builder(@input)];
    is_deeply($got, $expected, $name);
}
