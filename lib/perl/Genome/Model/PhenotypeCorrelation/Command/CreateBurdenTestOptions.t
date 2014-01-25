#!/usr/bin/env perl

use above 'Genome';
use Test::More;
use File::Temp qw/tempdir/;
use File::Slurp qw/read_file/;

use strict;
use warnings;

my $pkg = 'Genome::Model::PhenotypeCorrelation::Command::CreateBurdenTestOptions';
use_ok($pkg);

my $fake_clinical_data = <<EOS
Sample_Name\tT1\tT2\tT3
A\tHigh\t0.344\tX
B\tLow\t-1.256\tX
C\tLow\t-0.347\tY
D\tHigh\t1.434\tY
E\tNA\t0.0162\tZ
EOS
;

my $tmpdir = tempdir(CLEANUP => 1);
my $clinical_data_file = sprintf("%s/clinical_data.txt", $tmpdir);
my $fh = new IO::File($clinical_data_file, "w");
$fh->write($fake_clinical_data);
$fh->close;


my $cmd = $pkg->create(
        mutation_file => "matrix.txt",
        glm_clinical_data_file => $clinical_data_file,
        vep_annotation_file => "anno.vep",
        project_name => "Project",
        output_directory => $tmpdir,
        maf_cutoff => 0.02,
        missing_value_markers => [qw(LOL HAHA -4)],
        num_cores => 8,
        );
my @bare_expected = sort(split("\n", <<EOF
missing.data=c("LOL","HAHA","-4")
genotype.file="matrix.txt"
gfile.delimiter="\\t"
gfile.vid="Project"
gfile.sid="FIRST_ROW"
phenotype.file="$clinical_data_file"
pfile.delimiter="\\t"
pfile.sid="Sample_Name"
anno.file="anno.vep"
afile.delimiter="\\t"
afile.vid="Uploaded_variation"
gene.col="Gene"
vtype.col="Consequence"
vtype.use=c("ALL")
out.dir="$tmpdir"
if (!file.exists(out.dir)==T) dir.create(out.dir)
samplelist.dir="NONE"
covariates="NONE"
maf.cutoff=0.02
num.cores=8
EOF
));

chomp(@bare_expected);

ok($cmd, "created command");
ok($cmd->execute, "executed command");
my $output_file = $cmd->output_file;
ok(-s $output_file, "output file exists and is non-empty");

my @lines =
    sort
    grep {$_ ne ''}
    map {s/ *$//g; $_}
    map {s/#.*//g; $_}
    map {chomp; $_}
    split("\n", read_file($output_file));


is_deeply(\@lines, \@bare_expected, "output is as expected")
    or diag("expected:\n\t" . join("\n\t", @bare_expected) . "\n" .
            "actual:\n\t" . join("\n\t", @lines) . "\n");

done_testing();
