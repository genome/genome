#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::Exception;
use Test::More;

use above 'Genome';

class Thing::RemovesInputFiles {
    is => 'Command::V2',
    roles => [qw/
        Genome::InstrumentData::Command::Import::WorkFlow::Role::WithWorkingDirectory
        Genome::InstrumentData::Command::Import::WorkFlow::Role::RemovesInputFiles
    /],
    has => {
        property_is_input_file => { is => 'File', is_input => 1, },
        property_is_input_file2 => { is => 'File', is_input => 1, is_optional => 1, },
        property_is_input_not_file => { is => 'Boolean', is_input => 1, },
        property_not_input_is_file => { is => 'File', },
    },
};
sub Thing::RemovesInputFiles::execute { return 1; }

is_deeply(
    [ Thing::RemovesInputFiles->input_file_property_names ],
    [qw/ property_is_input_file property_is_input_file2 /],
    'input_file_property_names',
);
my $tmp_dir = Genome::Sys->create_temp_directory;
my $tmp_dir2 = Genome::Sys->create_temp_directory;
my %files;
my @input_file_properties_to_be_removed = (qw/ property_is_input_file /);
my @input_file_properties_to_be_kept = (qw/ property_is_input_file2 property_is_input_not_file property_not_input_is_file /);
my $cnt = 0;
for my $key ( @input_file_properties_to_be_removed, @input_file_properties_to_be_kept ) {
    $files{$key} = File::Spec->join((++$cnt % 2 == 1 ? $tmp_dir : $tmp_dir2 ), $key);
    print "FILE: $files{$key}\n";
    Genome::Sys->write_file($files{$key}, $key);
}

my $thing = Thing::RemovesInputFiles->execute(working_directory => $tmp_dir, %files);
ok($thing->result, 'execute thing');
is_deeply(
    [ $thing->input_files_to_remove ],
    [ map { glob( sprintf('%s*', $thing->$_) ) } @input_file_properties_to_be_removed ],
    'input_files_to_remove',
);
is(
    grep({ !-e $files{$_} } @input_file_properties_to_be_removed),
    @input_file_properties_to_be_removed,
    'removed correct files',
);
is(
    grep({ -e $files{$_} } @input_file_properties_to_be_kept),
    @input_file_properties_to_be_kept,
    'kept correct files',
);

done_testing();
