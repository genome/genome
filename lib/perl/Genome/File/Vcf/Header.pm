package Genome::File::Vcf::Header;

use Carp qw/confess/;
use Genome;

use strict;
use warnings;

class Genome::File::Vcf::Header {
    is => ['UR::Object'],
    has_transient_optional => [
        lines => {
            is => 'Text',
            doc => 'The raw text lines of the vcf header',
        },
        sample_names => {
            is => 'Text',
            doc => 'The sample names found in the vcf header',
            is_many => 1,
        },
    ],
};

sub create {
    my $class = shift;
    return $class->SUPER::create(@_);
}

sub parse {
    my ($self, @lines) = @_;
    my @header_lines;
    my $line_num = 0;
    while (my $line = shift @lines) {
        ++$line_num;
        if ($line =~ /^##/) {
            push(@header_lines, $line);
        } elsif ($line =~ /^#CHROM/) {
            my @fields = split("\t", $line);
            confess "Invalid vcf header (line $line_num)" if @fields < 7;
            $self->sample_names([@fields[9..$#fields]]) if @fields >= 9;
        }
    }
}

1;
