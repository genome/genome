#!/usr/bin/env genome-perl;

use strict;
use above 'Workflow::Model';

my $w = Workflow::Model->create_from_xml('Subjects.xml');

my @errors = $w->validate;

die 'Too many problems: ' . join("\n", @errors) unless $w->is_valid;

$w->as_png('/tmp/Subjects.png');

print $w->as_text;
exit;
