package Genome::Utility::Vcf;

use strict;
use warnings;

use Genome;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(open_vcf_file
                    get_vcf_header
                    get_vcf_header_arrayref
                    diff_vcf_file_vs_file
                    diff_vcf_text_vs_text
                    convert_string_to_hash
                    convert_hash_to_string
                    parse_vcf_line
                    deparse_vcf_line
                    get_samples_from_header
                    parse_header_format
                    convert_file_with_alpha_gt_values_to_numeric
                    convert_line_with_alpha_gt_values_to_numeric
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
    my ($filename) = @_;
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

# Return the header as an arrayref rather than a string
sub get_vcf_header_arrayref {
    my ($filename) = @_;

    if(wantarray()) {
        my ($header, $input_fh) = get_vcf_header($filename);
        my @header_array = split("\n", $header);
        return \@header_array, $input_fh;
    } else {
        my $header = get_vcf_header($filename);
        my @header_array = split("\n", $header);
        return \@header_array;
    }
}


# diffs vcf files ignoring the fileDate header entry.
sub diff_vcf_file_vs_file {
    my ($filename1, $filename2, %options) = @_;

    my $fh1 = open_vcf_file($filename1);
    my $fh2 = open_vcf_file($filename2);

    my $lines1 = [<$fh1>];
    my $lines2 = [<$fh2>];

    return diff_vcf_text_vs_text($lines1, $lines2, %options);
}

sub _strip_ignored_patterns {
    my ($lines, $patterns) = @_;
    return grep {
        my $x = $_;
        !grep {$x =~ /$_/i} @$patterns
    } @$lines
}

# diffs vcf formatted lines ignoring the fileDate header entry.
sub diff_vcf_text_vs_text {
    my ($lines1, $lines2, %options) = @_;
    my $ignore_patterns = delete $options{ignore_patterns};
    if (%options) {
        die "Unknown options to diff_vcf_text_vs_test: " . join(", ", keys(%options));
    }

    # we always ignore fileDate
    my %ignore = map {$_ => undef} (@$ignore_patterns, "fileDate");
    @$ignore_patterns = keys %ignore;

    my @clean_lines1 = _strip_ignored_patterns($lines1, $ignore_patterns);
    my @clean_lines2 = _strip_ignored_patterns($lines2, $ignore_patterns);

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

# Given a vcf line, parse it. Return info as a hash, return samples as an arrayref of hashes.
#FIXME probably move this to a base class
sub parse_vcf_line {
    my $line = shift;
    my $smpl_names = shift;
    my @sample_names = @$smpl_names;

    chomp $line;
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info_string, $format, @sample_data) = split "\t", $line;

    my $parsed_line;
    $parsed_line->{"chromosome"} = $chrom;
    $parsed_line->{"position"} = $pos;
    $parsed_line->{"id"} = $id;
    $parsed_line->{"reference"} = $ref;
    $parsed_line->{"alt"} = $alt;
    $parsed_line->{"qual"} = $qual;
    $parsed_line->{"filter"} = $filter;

    # Parse out the info field
    # Example: NS=3;DP=14;AF=0.5;DB;H2 (H2 is a flag)
    if ($info_string eq ".") {
        $parsed_line->{"info"} = $info_string;
    } else {
        my $info;
        my @info_values = split(";", $info_string);
        my @info_keys = ();
        for my $info_value (@info_values) {
            if($info_value =~ /=/) {
                my ($key, $value) = split("=", $info_value);
                unless (defined($value)) {
                    die("Invalid info tag. Has equal sign but no value.");
                }
                $info->{$key} = $value;
                push @info_keys, $key;
            }
            else {
                # This must mean the value is a "flag" in the vcf header, we could check this to be extra careful
                $info->{$info_value} = undef;
                push @info_keys, $info_value;
            }
        }
        $parsed_line->{"info"} = $info;
        $parsed_line->{'_info_tags'} = \@info_keys;
    }

    # Check and parse out the format/sample fields
    # FIXME do I want this to be an array or a hash?
    #my @parsed_samples;
    my $parsed_samples;
    my @format_fields = split ":", $format;
    $parsed_line->{'_format_fields'} = \@format_fields; #have to store this on each line since they can vary across each line
    for my $sample (@sample_data) {
        my $sample_values;
        my $sample_name = shift @sample_names;

        if ($sample eq ".") {
            for my $format_tag (@format_fields) {
                $sample_values->{$format_tag} = ".";
            }
        } else {
            my @sample_fields = split ":", $sample;
            unless (scalar(@format_fields) == scalar(@sample_fields) ) {
                die("Format field ($format) and sample ($sample) field do not have the same number of values");
            }
            for my $format_tag (@format_fields) {
                $sample_values->{$format_tag} = shift @sample_fields;
            }
        }

        $parsed_samples->{$sample_name} = $sample_values;
    }
    #$parsed_line->{"sample"} = \@parsed_samples;
    $parsed_line->{"sample"} = $parsed_samples;

    return $parsed_line;
}

sub deparse_vcf_line {
    my $parsed_line = shift;
    my $smpl_names = shift;
    my %parsed = %$parsed_line;

    #construct info fields
    my $info_fields;
    if($parsed_line->{info} eq '.') {
        $info_fields = ".";
    }
    else {
        for my $tag (@{$parsed_line->{'_info_tags'}}) {
            if($info_fields) {
                $info_fields .= ";";
            }
            $info_fields .= $tag;
            if(defined($parsed_line->{info}->{$tag})) {
                $info_fields .= '=' . $parsed_line->{info}->{$tag};
            }
            delete $parsed_line->{info}->{$tag}
        }
    }

    #do the format fields
    my $format_fields = join(":",@{$parsed_line->{'_format_fields'}});
    my @sample_data;
    for my $sample (@{$parsed_line->{sample}}{@$smpl_names}) {
        #dlarson 4/18/2012 Normally I would just apologize, but eat a slice of my fine hash pie.
        #This spits back the format fields in order and should be slightly more efficient by using the hash slice
        push @sample_data, join(":",@{$sample}{@{$parsed_line->{'_format_fields'}}});
    }
    my $line = join("\t", (@parsed{qw( chromosome position id reference alt qual filter )}, $info_fields, $format_fields, @sample_data));
    return $line . "\n";
}

# Given the header of a vcf, return an array of samples in the final header line
#FIXME probably move this to a base class
sub get_samples_from_header {
    my $header = shift;

    my $sample_line = @$header[-1];
    chomp $sample_line;
    my @fields = split "\t", $sample_line;
    splice(@fields, 0, 9);

    return @fields;
}

# Parse each FORMAT line in the header and store it in a hash
# example: ##FORMAT=<ID=FA,Number=1,Type=Float,Description="Fraction of reads supporting ALT">
sub parse_header_format {
    my $header = shift;

    my $header_format;
    for my $line (@$header) {
        next unless $line =~ /^##FORMAT=/;
        my ($value_string) = $line =~ (m/^##FORMAT=<(.*)>$/);
        unless (defined $value_string) {
            die("Format line malformed? : $line");
        }

        my ($id) = $value_string =~ (m/ID=(.*?),/);
        for my $key ("Number", "Type") {
            my ($value) = $value_string =~ (m/$key=(.*?),/);
            unless (defined $value) {
                die("FORMAT entry malformed? : $key in $value_string");
            }
            $header_format->{$id}->{$key} = $value;
        }
    }

    return $header_format;
}

# This method fixes all instances in a file where GT values are bases rather than numbers which refer to the alt
# Input: Vcf with alpha A/C/G/T GT values
# Output: Vcf with numberic GT values
sub convert_file_with_alpha_gt_values_to_numeric {
    my $input_file = shift;
    my $output_file = shift;

    unless (-s $input_file) {
        print "$input_file does not exist\n";
        return;
    }

    my $ofh = Genome::Sys->open_gzip_file_for_writing($output_file);
    my ($header, $ifh) = get_vcf_header_arrayref($input_file);
    my $header_string = join("\n", @$header);
    $ofh->print($header_string."\n");
    my @sample_names = get_samples_from_header($header);

    while (my $line = $ifh->getline) {
        my $parsed_line = parse_vcf_line($line, \@sample_names);
        my $fixed_line = convert_line_with_alpha_gt_values_to_numeric($parsed_line, \@sample_names);
        $ofh->print($fixed_line);
    }

    return 1;
}

sub convert_line_with_alpha_gt_values_to_numeric {
    my $parsed_line = shift;
    my $sample_names = shift;

    for my $sample_name (@$sample_names) {
        my $old_gt = $parsed_line->{sample}{$sample_name}{"GT"};
        die "No gt found?\n" unless $old_gt;

        # Don't fix it if it is already numeric
        unless ($old_gt =~ m/A|C|T|G/) {
            next;
        }

        my @old_gt_alleles = split "/", $old_gt;
        my @alt_alleles = split ",", $parsed_line->{alt};
        my @new_alleles;

        ALLELE: for my $allele (@old_gt_alleles) {
            if ($allele eq $parsed_line->{"reference"}) {
                push @new_alleles, 0;
            } else {
                for (my $index = 0; $index <= $#alt_alleles; $index++) {
                    my $current_alt = $alt_alleles[$index];
                    if ($current_alt eq $allele) {
                        push @new_alleles, $index+1;
                        next ALLELE;
                    }
                }

                # New alt encountered
                push @alt_alleles, $allele;
                push @new_alleles, scalar(@alt_alleles);
                $parsed_line->{alt} = join ",", @alt_alleles;
            }
        }

        unless (scalar(@new_alleles) == 2) {
            print "Could not get new alleles for line " . Data::Dumper::Dumper $parsed_line . "\n";
            return;
        }
        my $new_gt = join("/", sort(@new_alleles));
        $parsed_line->{sample}{$sample_name}{"GT"} = $new_gt;
    }

    my $deparsed_line = deparse_vcf_line($parsed_line, $sample_names);

    return $deparsed_line;
}
