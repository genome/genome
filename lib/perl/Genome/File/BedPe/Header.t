#!/usr/bin/env perl

use above 'Genome';
use Genome::File::BedPe::Entry;
use Test::More;

use strict;
use warnings;

my $pkg = 'Genome::File::BedPe::Header';
use_ok($pkg);

sub make_custom_fields_test_case {
    my %params = @_;
    my $header_pre_lines = delete $params{header_pre_lines};
    my $num_custom_fields = delete $params{num_custom_fields};
    if (%params) {
        die "Unexpected params to test case! what are you doing?"
    }

    my @custom_fields = map {"C$_"} 1..$num_custom_fields;

    my @fields = (@Genome::File::BedPe::Entry::ALL_FIELDS,
        @custom_fields);

    my $last_header_line = '#' . join("\t", @fields);
    my @header_lines = (@$header_pre_lines, $last_header_line);

    return sub {
        my $hdr = new Genome::File::BedPe::Header(\@header_lines);

        ok($hdr, 'Created header');
        is_deeply($hdr->{lines}, \@header_lines, 'Header text is as expected');

        is_deeply($hdr->{custom_fields}, [], 'No custom fields yet');

        ok($hdr->guess_custom_fields, 'We can guess custom fields!');
        is_deeply($hdr->{custom_fields}, \@custom_fields, 'We guessed correctly!');

        for my $i (0..$#custom_fields) {
            my $field_name = $custom_fields[$i];
            is($hdr->custom_field_index($field_name), $i,
                    "Custom field index for $field_name");
        }

        ok(!$hdr->custom_field_index("banana"), 'no custom field "banana"');
    };
}

subtest "Empty header" => sub {
    my $hdr = new Genome::File::BedPe::Header([]);
    ok($hdr, "Created empty header");
    is(ref $hdr->{lines}, "ARRAY", "lines is an arrayref");
    is($#{$hdr->{lines}}, -1, "lines is empty");
};

subtest "Just fields (1 custom)" => make_custom_fields_test_case(
    header_pre_lines => [],
    num_custom_fields => 1);

subtest "Just fields (3 custom)" => make_custom_fields_test_case(
    header_pre_lines => [],
    num_custom_fields => 3);

subtest "Fields with extra nonsense (1 custom)" => make_custom_fields_test_case(
    header_pre_lines => ['#hi', '#there'],
    num_custom_fields => 1);

subtest "Fields with extra nonsense (3 custom)" => make_custom_fields_test_case(
    header_pre_lines => ['#hi', '#there'],
    num_custom_fields => 3);

done_testing();
