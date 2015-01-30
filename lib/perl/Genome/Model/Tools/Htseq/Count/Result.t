#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 3;

use File::Spec;
use Sub::Install;

use above "Genome";
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::SoftwareResult::User;

my $pkg = 'Genome::Model::Tools::Htseq::Count::Result';

use_ok($pkg);

my $htseq_run_count = 0;
Sub::Install::reinstall_sub({
    into => $pkg,
    as => '_run_htseq_count',
    code => sub {
        my $self = shift;
        $htseq_run_count++;
        for my $file (qw(gene-counts.tsv transcript-counts.tsv)) {
            my $path = File::Spec->join($self->temp_staging_directory, $file);
            Genome::Sys->write_file($path, '');
        }
        return 1;
    },
});

my @ar;
for(1..3) {
    my $data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(lane => $_);
    push @ar, Genome::InstrumentData::AlignmentResult::PerLaneTophat->__define__(aligner_version => 'test', aligner_params => '', instrument_data => $data);
}

my @params = (
    alignment_results => \@ar,
    app_version => '0.5.4p1',
    result_version => 1,
    limit => 2000,
    minaqual => 1,
    mode => 'intersection-strict',
    users => Genome::Test::Factory::SoftwareResult::User->setup_user_hash(),
);

my $result = $pkg->get_or_create(@params);
isa_ok($result, $pkg, 'generated result');

is($htseq_run_count, 3, 'Htseq Count was run once per alignment result');
