#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

require File::Temp;
use IPC::Run;
use Test::More;

use_ok('Genome::Model::Tools::Consed::TracesToNav') or die;

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Consed-TracesToNav",
my $file_base = '10_126008345_126010576';
my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $refseq = $tmpdir."/$file_base.c1.refseq.fasta";
symlink(
    "$dir/$file_base.c1.refseq.fasta",
    $refseq,
);
ok (-s $refseq, 'linked refseq');

my $ace = $tmpdir."/10_126008345_126010576.ace.1";
symlink(
    "$dir/$file_base.ace.1",
    $ace,
);
ok(-s $ace, 'linked ace');

my $list = $tmpdir."/Nav.list";
symlink(
    "$dir/Nav.list",
    $list,
);
ok(-s $list, 'linked list');

# Run
my @command = ["gmt" , "consed" , "traces-to-nav" , "--ace" , "$ace" , "--convert-coords" , "$refseq" , "--input-type" , "simple" , "--name-nav" , "$file_base" , "--list" , "$list"];
&ipc_run(@command);

# Verify
my $date = &get_date_tag;
ok (-s $tmpdir."/$file_base.$date.nav", 'created navigator');
ok (-s $tmpdir."/$file_base.$date.csv", 'created spreadsheet');

done_testing();


###

sub ipc_run {

    my (@command) = @_;
    my ($in, $out, $err);
    IPC::Run::run(@command, \$in, \$out, \$err);
    if ($err) {
	print qq(Error $err\n);
    }
    if ($out) {
#	print qq($out\n);
    }
}

sub get_date_tag {
    
    my $time=`date`;
    #my ($handle) = getpwuid($<);
    my $date = `date +%D`;
    (my $new_date) = $date =~ /(\d\d)\/(\d\d)\/(\d\d)/ ;
    my $date_tag = "$3$1$2" ;
    return $date_tag;
}
