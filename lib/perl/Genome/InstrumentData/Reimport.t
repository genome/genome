#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::Exception;
use Test::More;

my $class = 'Genome::InstrumentData::Reimport';
use_ok($class) or die;

my $attribute_label_for_reimported_from = $class->attribute_label_for_reimported_from;
ok($attribute_label_for_reimported_from, 'attribute_label_for_reimported_from');
ok($class->attribute_labels_to_ignore_when_reimporting, 'attribute_labels_to_ignore_when_reimporting');

throws_ok { $class->attributes_for_reimport_from_instrument_data; } qr/^No instrument data given!/;
my $instdata = Genome::InstrumentData::Solexa->__define__(
    run_name => 'RUN',
    subset_name => 'SUBSET_NAME',
    library => Genome::Library->__define__(name => 'LIBRARY'),
);
my %attrs = (
    lane => 1,
    index_sequence => 'AAGGTTCC',
);
for my $attribute_label ( keys %attrs ) {
    $instdata->add_attribute(attribute_label => $attribute_label, attribute_value => $attrs{$attribute_label});
}
my %expected_attrs = %attrs;
$expected_attrs{library_name} = $instdata->library->name;
$expected_attrs{run_name} = $instdata->run_name;
$expected_attrs{subset_name} = $instdata->subset_name;
$expected_attrs{ $class->attribute_label_for_reimported_from } = $instdata->id;
my $reimport = $class->attributes_for_reimport_from_instrument_data($instdata);
is_deeply($reimport, \%expected_attrs, 'attributes_for_reimport');

is_deeply(
    [ $class->headers_for_reimport_attributes($reimport) ],
    [qw/ library_name source_files index_sequence lane /, $attribute_label_for_reimported_from, qw/ run_name subset_name /],
    'headers_for_attributes_for_reimport',
);

my $reimported_from = Genome::InstrumentData::Imported->__define__();
$instdata->add_attribute(attribute_label => $attribute_label_for_reimported_from, attribute_value => $reimported_from->id);
is($class->reimported_from($instdata), $reimported_from, 'reimported_from');

done_testing();
