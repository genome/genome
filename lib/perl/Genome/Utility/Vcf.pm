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
                    convert_string_to_hash
                    convert_hash_to_string
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

# return the vcf header as a string.
# if evaluated in ARRAY context it also returns the open file-handle 
# positioned at the first line that is not a header.
sub get_vcf_header {
    my $filename = @_;
    my $input_fh = open_vcf_file($filename);

    my $header = '';
    my $header_end = 0;
    while (!$header_end) {
        my $line = $input_fh->getline();
        if ($line =~ m/^##/) {
            $header .= $line;
        } elsif ($line =~ m/^#/) {
            $header .= $line;
            $header_end = 1;
        } else {
            Carp::croak("Missed the final header line with the" .
                    " sample list? Last line examined: $line Header so far: " .
                    $header);
        }
    }

    if(wantarray()) {
        return $header, $input_fh;
    } else {
        $input_fh->close();
        return $header;
    }
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

# return the hash representation of the formatted string
sub convert_string_to_hash {
    my ($string, $format) = @_;

    my @keys = split(":", $format);
    my @values = split(":", $string);

    my %result;
    for my $i (0..$#keys) {
        $result{$keys[$i]} = $values[$i];
    }
    return %result;
}

# return the formatted string representation of the sample hash.
sub convert_hash_to_string {
    my ($hash, $format) = @_;

    my @keys = split(":", $format);
    my $result = '';
    my $count = 0;
    for my $key (@keys) {
        $result .= $hash->{$key};
        $result .= ":" if $count < $#keys;
        $count += 1;
    }

    return $result;
}



