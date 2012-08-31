#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 13;
use File::Compare;
use File::Temp;
use File::Basename;
use Bio::SeqIO;

BEGIN {
        use_ok ('Genome::Model::Tools::Fasta::Orient')
            or die;
}

#< Test if it really works >#
my $tmp_dir = Genome::Sys->create_temp_directory('Genome-Model-Tools-Fasta-Orient');
my $orig_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fasta-Orient';
my %file_path = (
    fasta => '/assembled.fasta',
    qual => '/assembled.fasta.qual',
    s_fasta => '/primers_sense.fasta',
    as_fasta => '/primers_anti_sense.fasta',
);
for my $key (keys %file_path) {
    my $orig_path = $orig_dir . $file_path{$key};
    my $tmp_path = $tmp_dir . $file_path{$key};
    Genome::Sys->create_symlink($orig_path, $tmp_path);
    $file_path{$key} = $tmp_path;
}

my $orient = Genome::Model::Tools::Fasta::Orient->create(
    fasta_file       => $file_path{'fasta'},
    sense_fasta_file => $file_path{'s_fasta'},
    anti_sense_fasta_file => $file_path{'as_fasta'},
);

my $confirmed_fasta_file = $orient->confirmed_fasta_file;
my $confirmed_qual_file  = $orient->confirmed_qual_file;
my $unconfirmed_fasta_file = $orient->unconfirmed_fasta_file;
my $unconfirmed_qual_file  = $orient->unconfirmed_qual_file;

my $confirmed_fasta_file_name = basename $confirmed_fasta_file;
my $confirmed_qual_file_name  = basename $confirmed_qual_file;
my $unconfirmed_fasta_file_name = basename $unconfirmed_fasta_file;
my $unconfirmed_qual_file_name  = basename $unconfirmed_qual_file;

is($confirmed_fasta_file_name,'assembled.confirmed.fasta', "confirmed_fasta_file return ok");
is($confirmed_qual_file_name,'assembled.confirmed.fasta.qual', "confirmed_qual_file return ok");
is($unconfirmed_fasta_file_name,'assembled.unconfirmed.fasta', "unconfirmed_fasta_file return ok");
is($unconfirmed_qual_file_name,'assembled.unconfirmed.fasta.qual', "unconfirmed_qual_file return ok");

ok($orient->execute, "Orient test fine");

is(compare($confirmed_fasta_file, "$orig_dir/expected.confirmed.fasta"), 0, 'Expected and generated confirmed fasta matches');
is(compare($confirmed_qual_file, "$orig_dir/expected.confirmed.fasta.qual"), 0, 'Expected and generated confirmed qual matches');
is(compare($unconfirmed_fasta_file, "$orig_dir/expected.unconfirmed.fasta"), 0, 'Expected and generated unconfirmed fasta matches');
is(compare($unconfirmed_qual_file, "$orig_dir/expected.unconfirmed.fasta.qual"), 0, 'Expected and generated unconfirmed qual matches');

unlink $confirmed_fasta_file;
unlink $confirmed_qual_file;
unlink $unconfirmed_fasta_file;
unlink $unconfirmed_qual_file;

#< Test failing conditions #>
# no primer files
$orient = Genome::Model::Tools::Fasta::Orient->create(
    fasta_file => $file_path{'fasta'},
);
ok(!$orient, "Neither sense nor anti_sense file provided");

# non existing primer files
$orient = Genome::Model::Tools::Fasta::Orient->create(
    fasta_file            => $file_path{'fasta'},
    anti_sense_fasta_file => $tmp_dir.'/no_way_this_exists.fasta',
);
ok(!$orient, "anti_sense file is invalid");

$orient = Genome::Model::Tools::Fasta::Orient->create(
    fasta_file            => $file_path{'fasta'},
    sense_fasta_file => $tmp_dir.'/no_way_this_exists.fasta',
);
ok(!$orient, "anti_sense file is invalid");

exit;

#$HeadURL$
#$Id$
