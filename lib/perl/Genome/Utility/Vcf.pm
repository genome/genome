package Genome::Utility::Vcf;

use strict;
use warnings;

use Genome;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(open_vcf_file
                    diff_vcf_file_vs_file
                    diff_vcf_text_vs_text
                    );

# opens up either a .vcf file or a .vcf.gz file for reading.
sub open_vcf_file {
    my ($filename) = @_;

    my $fh;
    my $file_type = Genome::Sys->file_type($filename);
    #NOTE: debian bug #522441 - `file` can report gzip files as any of these....
    if ($file_type eq "gzip" or $file_type eq "Sun" or $file_type eq "Minix") {
        $fh = Genome::Sys->open_gzip_file_for_reading($filename);
    } else {
        $fh = Genome::Sys->open_file_for_reading($filename);
    }
    return $fh;
}

# diffs vcf files ignoring the fileDate header entry.
sub diff_vcf_file_vs_file {
    my ($filename1, $filename2) = @_;

    my $fh1 = open_vcf_file($filename1);
    my $fh2 = open_vcf_file($filename2);

    my $lines1 = [<$fh1>];
    my $lines2 = [<$fh2>];

    return diff_vcf_text_vs_text($lines1, $lines2);
}

sub _strip_fileDate_from_header {
    my @lines = @_;
    return grep {!/fileDate/} @lines;
}

# diffs vcf formatted lines ignoring the fileDate header entry.
sub diff_vcf_text_vs_text {
    my ($lines1, $lines2) = @_;

    my @clean_lines1 = _strip_fileDate_from_header(@{$lines1});
    my @clean_lines2 = _strip_fileDate_from_header(@{$lines2});

    my $t1 = join('', @clean_lines1);
    my $t2 = join('', @clean_lines2);

    return Genome::Sys->diff_text_vs_text($t1, $t2);
}

