#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Regexp::Common;
use Test::More;

use_ok('Genome::Model::Tools::Sx::Validate') or die;

class Genome::Model::Tools::Sx::Tester {
    has => [
        param => { is => 'Number', },
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
ok(Genome::Model::Tools::Sx::Validate->validate_command($good_cmd), 'good command validates');

$good_cmd = 'gmt sx tester -param 10';
ok(Genome::Model::Tools::Sx::Validate->validate_command($good_cmd), 'good command validates');

my $bad_cmd = 'gmt tools';
ok(!Genome::Model::Tools::Sx::Validate->validate_command($bad_cmd), 'bad command does not validate');

$bad_cmd = 'gmt sx no-way';
ok(!Genome::Model::Tools::Sx::Validate->validate_command($bad_cmd), 'bad command does not validate');

$bad_cmd = 'gmt sx tester -invalid-param 1';
ok(!Genome::Model::Tools::Sx::Validate->validate_command($bad_cmd), 'bad command does not validate');

$bad_cmd = 'gmt sx tester -param A';
ok(!Genome::Model::Tools::Sx::Validate->validate_command($bad_cmd), 'bad command does not validate');

done_testing();
