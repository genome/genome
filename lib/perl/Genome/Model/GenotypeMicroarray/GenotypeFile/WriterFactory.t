#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Test::More;

use_ok('Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory') or die;
use_ok('Genome::File::Vcf::Header') or die;
no warnings;
*Genome::File::Vcf::Header::to_string = sub{ return ''; };
use warnings;

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $cnt = 0;
my $output_file;
my $header = bless({}, 'Genome::File::Vcf::Header');
my $writer_params_generator = sub{ 
    my %params = @_;
    if ( exists $params{output} ) {
        $params{output} = $output_file = $tmpdir.'/FILE'.++$cnt;
    }
    return ( 
        header => $header,
        string => join(':', map { join('=', $_, $params{$_}) } keys %params),
    );
};

## TEST ERRORS ##
# Missing required params
ok(!Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer(), 'failed to build writer w/o header');
is(Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->error_message, 'No header given to create writer!', 'correct error');

# Invalid format
ok(!Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer( $writer_params_generator->(format => 'supervcf') ), 'failed to create writer w/ invalid format');

# Dup key
ok(!Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer( $writer_params_generator->(format => 'vcf').':format=vcf'), 'failed to create writer w/ dup key');

### VCF ###
# Default is VCF to STDOUT
my $writer = Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer( $writer_params_generator->() );
isa_ok($writer, 'Genome::File::Vcf::Writer');
is($writer->{name}, '-', 'output is STDOUT');

# Output w/o key is specified
$writer = Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer( $writer_params_generator->(output => '') );
isa_ok($writer, 'Genome::File::Vcf::Writer');
is($writer->{name}, $output_file, "original output is $output_file");

# Output and format specified
$writer = Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer( $writer_params_generator->(output => '', format => 'vcf') );
isa_ok($writer, 'Genome::File::Vcf::Writer');
is($writer->{name}, $output_file, "original output is $output_file");

## CSV ###
# Defaults
$writer = Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer( $writer_params_generator->(output => '', format => 'csv', separator => 'TAB') );
isa_ok($writer, 'Genome::Model::GenotypeMicroarray::GenotypeFile::WriteCsv');
is($writer->get_original_output, $output_file, "original output is $output_file");
is($writer->separator, "\t", 'separator is TAB');
is_deeply($writer->headers, [qw/ chromosome position alleles reference id sample_name log_r_ratio gc_score cnv_value cnv_confidence allele1 allele2 /], 'headers are correct');
ok($writer->print_headers, 'print_headers is true');

# Everything specified
$writer = Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory->build_writer( $writer_params_generator->(output => '', format => 'csv', separator => ',', fields => 'chromosome,allele1', headers => 0, in_place_of_null_value => 'NULL'),
);
isa_ok($writer, 'Genome::Model::GenotypeMicroarray::GenotypeFile::WriteCsv');
is($writer->get_original_output, $output_file, "original output is $output_file");
is($writer->separator, ',', 'separator is comma');
is_deeply($writer->headers, [qw/ chromosome allele1 /], 'headers are correct');
ok(!$writer->print_headers, 'print_headers is false');

done_testing();
