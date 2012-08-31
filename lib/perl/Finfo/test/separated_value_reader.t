#!/usr/bin/env genome-perl

use strict;
use warnings;

use Data::Dumper;
use Finfo::SeparatedValueReader;
use IO::String;

use Test::More tests => 29;
use Test::Differences;

my $class = 'Finfo::SeparatedValueReader';
my @expected_headers = (qw/ chromosome_name start stop reference variant reference_type variant_type reference_reads variant_reads maq_score /);

# Tests using comma delimiter
my $csv_reader = Finfo::SeparatedValueReader->new(
    io => _create_csv_io(),
);
isa_ok($csv_reader,$class);
is_deeply($csv_reader->headers, \@expected_headers, 'header parsed correctly for csv');
my @csv_data = $csv_reader->all;
is(scalar(@csv_data), 20, 'got 20 from csv reader');

# Tests using tab delimiter
my $tsv_reader = $class->new(
    io => _create_tsv_io(),
    separator => '\t',
    is_regex => 1,
);
isa_ok($tsv_reader, $class);
is_deeply($tsv_reader->headers, \@expected_headers, 'header parsed correctly for tsv');
my @tsv_data = $tsv_reader->all;
is(scalar(@tsv_data), 20, 'got 20 from tsv reader');

# tests using many space delimiter
my $ssv_reader = $class->new(
    io => _create_ssv_io(),
    separator => '\s+',
    is_regex => 1,
);
isa_ok($ssv_reader,$class);
is_deeply($ssv_reader->headers, \@expected_headers, 'header parsed correctly for ssv');
my @ssv_data = $ssv_reader->all;
is(scalar(@ssv_data), 20, 'got 20 from ssv reader');

# tests using odd delimiter
my $osv_reader = $class->new(
    io => _create_osv_io(),
    separator => '#\w+#',
    is_regex => 1,
);
isa_ok($osv_reader,$class);
is_deeply($osv_reader->headers, \@expected_headers, 'header parsed correctly for osv');
my @osv_data = $osv_reader->all;
is(scalar(@osv_data), 20, 'got 20 from osv reader');

# Test equality of data
is_deeply(\@tsv_data, \@csv_data, 'data produced by tab and comma delimited files');
is_deeply(\@ssv_data, \@csv_data, 'data produced by multi space and comma delimited files');
is_deeply(\@osv_data, \@csv_data, 'data produced by multi space and comma delimited files');

# Test error for uneven headers/values
my $uneven_reader = $class->new(
    io => _create_uneven_io(),
);
is_deeply($uneven_reader->headers, \@expected_headers, 'header parsed correctly for csv');
my $eval;
eval {
    $eval = $uneven_reader->all;
};
is($eval, undef, 'caught error, eval\'d correctly');
is($uneven_reader->line_number, 22, 'stopped on line 22');

# Test null fields 
my $null_reader = $class->new(
    io => _create_null_fields_io(),
);
isa_ok($null_reader,$class);
is_deeply($null_reader->headers, \@expected_headers, 'header parsed correctly for null');
my @null_eval;
eval {
    @null_eval = $null_reader->all;
};
is(defined @null_eval, 1, 'no error, eval\'d correctly');
is(scalar @null_eval, 41, 'got null fields w/o error');
is(scalar keys %{pop @null_eval}, 10, 'null fields line has correct # of columns');

# Test overly long row
my $long_reader = $class->new(
    io => _create_too_long_io(),
);
is_deeply($long_reader->headers, \@expected_headers, 'header parsed correctly for csv');
my $long_eval;
eval {
    $long_eval = $long_reader->all;
};
is($long_eval, undef, 'caught error, eval\'d correctly');
is($long_reader->line_number, 22, 'stopped on line 22');

# Test too short row
my $short_reader = $class->new(
    io => _create_too_short_io(),
);
is_deeply($short_reader->headers, \@expected_headers, 'header parsed correctly for csv');
my $short_eval;
eval {
    $short_eval = $short_reader->all;
};
is($short_eval, undef, 'caught error, eval\'d correctly');
is($short_reader->line_number, 22, 'stopped on line 22');


exit 0;

sub _string {
    return <<"EOS"
22,14431069,14431069,A,C,ref,SNP,49,1,10
22,14431684,14431684,C,A,ref,SNP,56,1,18
22,14431687,14431687,C,A,ref,SNP,59,1,19
22,14433379,14433379,T,G,ref,SNP,71,1,3
22,14433766,14433766,C,T,ref,SNP,73,3,3
22,14434538,14434538,C,G,ref,SNP,62,2,3
22,14434540,14434540,G,T,ref,SNP,69,3,3
22,14434548,14434548,G,T,ref,SNP,93,2,10
22,14435207,14435207,C,A,ref,SNP,71,1,8
22,14436240,14436240,C,T,ref,SNP,60,1,2
22,14436332,14436332,A,C,ref,SNP,70,4,8
22,14436338,14436338,A,C,ref,SNP,61,7,8
22,14436360,14436360,C,A,ref,SNP,61,1,17
22,14436374,14436374,T,C,ref,SNP,61,1,6
22,14436547,14436547,T,G,ref,SNP,120,2,3
22,14436554,14436554,T,A,ref,SNP,112,2,3
22,14436895,14436895,G,C,ref,SNP,46,1,5
22,14437012,14437012,G,C,ref,SNP,67,1,2
22,14437058,14437058,G,T,ref,SNP,83,3,2
22,14437065,14437065,C,T,ref,SNP,78,1,2
EOS
}

sub _create_csv_io {
    my $io = IO::String->new();
    $io->print( join(',', @expected_headers), "\n" );
    $io->print( _string() );
    $io->seek(0, 0);
    return $io;
}

sub _create_tsv_io {
    my $io = IO::String->new();
    $io->print( join("\t", @expected_headers), "\n" );
    my $string = _string();
    $string =~ s/,/\t/g;
    $io->print($string);
    $io->seek(0, 0);
    return $io;
}

sub _create_ssv_io {
    my $io = IO::String->new();
    $io->print( join('          ', @expected_headers), "\n" );
    my $string = _string();
    $string =~ s/,/       /g;
    $io->print($string);
    $io->seek(0, 0);
    return $io;
}

sub _create_osv_io {
    my $io = IO::String->new();
    $io->print( join('#SEP#', @expected_headers), "\n" );
    my $string = _string();
    $string =~ s/,/#SEP#/g;
    $io->print($string);
    $io->seek(0, 0);
    return $io;
}

sub _create_uneven_io {
    my $io = IO::String->new();
    $io->print( join(',', @expected_headers), "\n" );
    $io->print( _string() );
    $io->print("a,b,c\n");
    $io->print( _string() );
    $io->seek(0, 0);
    return $io;
}

sub _create_null_fields_io {
    my $io = IO::String->new();
    $io->print( join(',', @expected_headers), "\n" );
    $io->print( _string() );
    $io->print(",b,,,e,f,,h,i,\n");
    $io->print( _string() );
    $io->seek(0, 0);
    return $io;
}

sub _create_too_long_io {
    my $io = IO::String->new();
    $io->print( join(',', @expected_headers), "\n" );
    $io->print( _string() );
    $io->print("a,b,c,d,e,f,g,h,i,j,k,l,m\n");
    $io->print( _string() );
    $io->seek(0, 0);
    return $io;
}

sub _create_too_short_io {
    my $io = IO::String->new();
    $io->print( join(',', @expected_headers), "\n" );
    $io->print( _string() );
    $io->print("a,b,c,d,e,f,g,h,i\n");
    $io->print( _string() );
    $io->seek(0, 0);
    return $io;
}

#$Header$
#$Id$
