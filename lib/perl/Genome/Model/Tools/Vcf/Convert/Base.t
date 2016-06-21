#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;

my $class = 'Genome::Model::Tools::Vcf::Convert::Base';
use_ok($class);

subtest 'query_tcga_barcode' => sub {
    plan tests => 3;

    my $test_class = join('::', $class, 'Test');
    UR::Object::Type->define(
        class_name => $test_class,
        is => $class,
    );

    my $self = $test_class->create();
    $self->dump_error_messages(0);

    {
        my $barcode = ['TCGA-A2-A25B-01A-11R-A169-07'];
        my $content = $self->query_tcga_barcode($barcode, sleep => 1);
        ok($content, qq(received content for barcode: $barcode));
    }

    {
        my $barcode = ['NOT-A-BARCODE'];
        my $content = eval { $self->query_tcga_barcode($barcode, sleep => 1) };
        my $error = $@;
        ok(!$content, qq(did not receive content for barcode: $barcode));
        my $error_re = sprintf($class->query_tcga_barcode_error_template, '.*', '.*');
        like($error, qr/$error_re/, qq(expected exception thrown for barcode: $barcode));
    }
};
