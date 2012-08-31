#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use Test::More;

use_ok('Genome::InstrumentData::Command::Dacc::Status') or die;

my $status = Genome::InstrumentData::Command::Dacc::Status->create(
    sra_sample_id => 'SRS000000',
    format => 'sff',
);
ok($status, 'create');
$status->dump_status_messages(1);
ok($status->execute, 'execute'); # will fail at sample., but does some testing

done_testing();
exit;

=pod

=head1 Tests

=head1 Disclaimer

 Copyright (C) 2010 Washington University Genome Sequencing Center

 This script is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Eddie Belter <ebelter@watson.wustl.edu>

=cut

