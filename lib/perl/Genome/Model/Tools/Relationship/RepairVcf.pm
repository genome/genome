package Genome::Model::Tools::Relationship::RepairVcf;

use strict;
use warnings;
use List::Util qw(first);
use Genome::Utility::Vcf qw(convert_string_to_hash
                            convert_hash_to_string);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(fix_alt_and_GT_field);

# return the corrected alt, format, and sample strings
#   deduplicates the ref in alt and correctly updates
#   the GT fields.
sub fix_alt_and_GT_field {
    my $ref = shift;
    my $alt = shift;
    my $format = shift;
    my @samples = @_;

    my @bases = ($ref, split(",", $alt));

    # 1) replace numbers in format with base-letters
    my @hashes;
    for my $sample (@samples) {
        my %hash = convert_string_to_hash($sample, $format);
        my $value = $hash{GT};
        if(defined($value)) {
            my $new_value = '';
            for my $i (0..(length($value)-1)) {
                my $token = substr($value, $i, 1);
                if($token =~ m/\d/) {
                    $new_value .= $bases[$token];
                } else {
                    $new_value .= $token;
                }
            }
            $hash{GT} = $new_value;
        }
        push(@hashes, \%hash);
    }

    # 2) remove ref from alt
    @bases = ($ref);
    my @fixed_alt;
    for my $a (split(",", $alt)) {
        if($ref ne $a) {
            push(@bases, $a);
            push(@fixed_alt, $a);
        }
    }
    my $fixed_alt;
    if(scalar(@fixed_alt)) {
        $fixed_alt = join(",", @fixed_alt);
    } else {
        $fixed_alt = '.';
    }

    # 3) replace base-letters with numbers
    my %pos;
    for my $i (0..$#bases) {
        my $b = $bases[$i];
        $pos{$b} = $i;
    }

    my @fixed_samples;
    for my $hash (@hashes) {
        my %hash = %{$hash};
        my $value = $hash{GT};
        if(defined($value)) {
            my $new_value = '';
            for my $i (0..(length($value)-1)) {
                my $token = substr($value, $i, 1);
                if($token =~ m/[A-Z]/) {
                    $new_value .= $pos{$token};
                } else {
                    $new_value .= $token;
                }
            }
            $hash{GT} = $new_value;
        }

        # ensure GT is first part of the $format
        my @elements = split(":", $format);
        my @new_elements;
        my $GT_index = first {$_ eq 'GT'} @elements;
        if(defined($GT_index)) {
            @new_elements = ("GT");
        }
        for my $element (@elements) {
            if($element ne 'GT') {
                push(@new_elements, $element);
            }
        }
        $format = join(":", @new_elements);

        push(@fixed_samples, convert_hash_to_string(\%hash, $format));
    }
    return ($fixed_alt, $format, @fixed_samples);
}

1;
