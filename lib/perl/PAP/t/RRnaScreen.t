
use strict;
use warnings;
use above;
use Test::More tests => 5;

BEGIN {
    use_ok('PAP::Command');
    use_ok('PAP::Command::RRnaScreen');
}

SKIP: {
    skip "this module will probably be remove; it's not used", 3 unless 0;
my $tempdir = File::Temp::tempdir(
                                  'PAP_RRnaScreen_test_XXXXXXXX',
                                  DIR     => '/tmp',
                                  CLEANUP => 1,
                                 );


my $command = PAP::Command::RRnaScreen->create(
                                           'fasta_file'      => 'data/B_coprocol
a.chunk.fasta', 
                                           'report_save_dir' => $tempdir, 
                                           'blast_report'    => 'data/rrna_screen_blast_report',
                                          );
                                          
isa_ok($command, 'PAP::Command::RRnaScreen');

ok($command->parse_result(),'parse blast result');
my @query_names = @{$command->dead_genes()};
#print $#query_names, "\n"; # should be an array with 5 items...
is($#query_names, 4, 'number of results from parsed output');

}



