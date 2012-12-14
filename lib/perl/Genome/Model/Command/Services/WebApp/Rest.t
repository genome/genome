#!/usr/bin/env genome-perl

use strict;
use warnings;
use Test::More;
plan tests => 1;

use above 'Genome';

my $restapp = require Genome::Model::Command::Services::WebApp->base_dir . '/Rest.psgi';

ok( $restapp, 'loaded Rest.psgi' );
