#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Data::Dumper;
use File::Temp;
use Test::More tests => 15;

BEGIN {
        use_ok('Genome::Model::Tools::WuBlast::Xdformat');
        use_ok('Genome::Model::Tools::WuBlast::Xdformat::Create');
        use_ok('Genome::Model::Tools::WuBlast::Xdformat::Append');
        use_ok('Genome::Model::Tools::WuBlast::Xdformat::Verify');
}

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $database = $tmp_dir .'/test_db';

my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-WuBlast-Xdformat';
my @fasta_names = (qw/ 11 12 13 /);
my @fasta_files = map { sprintf('%s/%s.fa', $test_data_dir, $_) } @fasta_names;
is(grep({ -s } @fasta_files), 3, 'The three fasta files exist');

# CREATE
my $create_success = Genome::Model::Tools::WuBlast::Xdformat::Create->create(
                                                                     database => $database,
                                                                     fasta_files => \@fasta_files,
                                                                 );
isa_ok($create_success,'Genome::Model::Tools::WuBlast::Xdformat::Create');
ok($create_success->execute,'execute command '. $create_success->command_name);

# The object should never get created because the database already exists.
my $create_fail = Genome::Model::Tools::WuBlast::Xdformat::Create->create(
                                                                          database => $database,
                                                                          fasta_files => \@fasta_files,
                                                                      );
ok($create_fail, 'Create w/ a duplicate db');
ok(!$create_fail->execute, 'Execute failed as exepected');

# APPEND
my $append_success = Genome::Model::Tools::WuBlast::Xdformat::Append->create(
                                                                             database => $database,
                                                                             fasta_files => \@fasta_files,
                                                                         );
isa_ok($append_success,'Genome::Model::Tools::WuBlast::Xdformat::Append');
ok($append_success->execute,'execute command '. $append_success->command_name);

# VERIFY
my $verify_success = Genome::Model::Tools::WuBlast::Xdformat::Verify->create(
                                                                     database => $database,
                                                                 );
isa_ok($verify_success,'Genome::Model::Tools::WuBlast::Xdformat::Verify');
ok($verify_success->execute,'execute command '. $verify_success->command_name);

# Remove a necessary database file and try to verify again(should fail)
unlink($database.'.xnt') || die "Failed to remove database file for test";
my $verify_fail = Genome::Model::Tools::WuBlast::Xdformat::Verify->create(
                                                                  database => $database,
                                                              );
isa_ok($verify_fail,'Genome::Model::Tools::WuBlast::Xdformat::Verify');
ok(!$verify_fail->execute,'expected verify to fail');

exit;

#$HeadURL$
#$Id$
