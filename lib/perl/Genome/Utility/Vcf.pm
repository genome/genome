package Genome::Utility::Vcf;

use strict;
use warnings;

use Genome;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(open_vcf_file
                    get_vcf_header
                    diff_vcf_file_vs_file
                    diff_vcf_text_vs_text
                    );

# opens up either a .vcf file or a .vcf.gz file for reading.
sub open_vcf_file {
    my ($filename) = @_;

    my $fh;
    if(Genome::Sys->file_is_gzipped($filename)) {
        $fh = Genome::Sys->open_gzip_file_for_reading($filename);
    } else {
        $fh = Genome::Sys->open_file_for_reading($filename);
    }
    return $fh;
}

# returns only the header of the vcf file in one big string.
sub get_vcf_header {
    my ($filename) = @_;

    my $fh = open_vcf_file($filename);
    my $header = "";
    while(my $line = $fh->getline()) {
        if($line =~ m/^#/) {
            $header .= $line;
        } else {
            chomp($header);
            return $header
        }
    }
    Carp::croak("Vcf file ($filename) has no body... only header.");
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

