package Genome::Annotation::Detail::TestHelpers;

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Set::Scalar;

use Exporter 'import';

our @EXPORT_OK = qw(test_cmd_and_result_are_in_sync);

sub test_cmd_and_result_are_in_sync {
    my $cmd = shift;

    my $cmd_set = Set::Scalar->new($cmd->input_names);
    my $sr_set = Set::Scalar->new($cmd->software_result->param_names,
        $cmd->software_result->metric_names, $cmd->software_result->input_names);
    is_deeply($cmd_set - $sr_set, Set::Scalar->new(), 'All command inputs are persisted SoftwareResult properties');
}
