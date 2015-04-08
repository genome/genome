#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Regexp::Common;
use Test::Exception;
use Test::More;

use_ok('Genome::Model::Tools::Sx::Validate') or die;

class Genome::Model::Tools::Sx::Tester {
    has => [
        param => { is => 'Number', },
        thread => { is => 'Number', },
    ],
};
sub Genome::Model::Tools::Sx::Tester::__errors__ {
    my $self = shift;

    my @errors;
    my $param = $self->param;
    if ( not $param or $param !~ /^$RE{num}{int}$/ or $param < 0 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ param /],
            desc => 'Param must be a positive integer: '.$param,
        );
    }

    return @errors;
}

my $good_cmd = 'gmt sx';
my $validator = Genome::Model::Tools::Sx::Validate->create(command => $good_cmd);
throws_ok(sub{ $validator->number_of_threads_required; }, qr/No command params set/, 'failed to get number_of_threads_required beofre validating command');
ok($validator->is_valid, 'is_valid');
ok($validator->command_class, 'command_class');
ok($validator->command_params, 'command_params');
is($validator->number_of_threads_required, 1);

$good_cmd = 'gmt sx tester -param 10 -thread 8';
$validator = Genome::Model::Tools::Sx::Validate->create(command => $good_cmd);
ok($validator->is_valid, 'good command validates');
is($validator->number_of_threads_required, 8);

my $bad_cmd = 'gmt tools';
ok(!Genome::Model::Tools::Sx::Validate->create(command => $bad_cmd)->is_valid, 'bad command does not validate');

$bad_cmd = 'gmt sx no-way';
ok(!Genome::Model::Tools::Sx::Validate->create(command => $bad_cmd)->is_valid, 'bad command does not validate');

$bad_cmd = 'gmt sx tester -invalid-param 1';
ok(!Genome::Model::Tools::Sx::Validate->create(command => $bad_cmd)->is_valid, 'bad command does not validate');

$bad_cmd = 'gmt sx tester -param A';
ok(!Genome::Model::Tools::Sx::Validate->create(command => $bad_cmd)->is_valid, 'bad command does not validate');

done_testing();
