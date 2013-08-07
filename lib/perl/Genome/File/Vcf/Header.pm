package Genome::File::Vcf::Header;

use Carp qw/confess/;
use Genome;

use strict;
use warnings;

my %VALID_VCF_TYPE_NON_NUMERIC_NUMBERS = (
    'A' => 1, # per alt
    'G' => 1, # per genotype
    '.' => 1, # variable length
);

my %VALID_VCF_TYPE_NAMES = (
    String => 1,
    Integer => 1,
    Float => 1,
    Flag => 1,
    Character => 1,
);

class Genome::File::Vcf::Header {
    is => ['UR::Object'],
    has => [
        lines => {
            is => 'Text',
            doc => 'The raw text lines of the vcf header',
            is_many => 1,
        },
    ],
    has_transient_optional => [
        sample_names => {
            is => 'Text',
            doc => 'The sample names found in the vcf header',
            is_many => 1,
        },
        fileformat => {
            is => 'Text',
            doc => 'The vcf spec version for this header/file',
            is_optional => 1, # set by _parse and friends
        },
        _unparsed_lines => {
            is => 'ARRAY',
            doc => 'An array of lines that were not specifically interpreted',
            default_value => [],
        },
        format_types => {
            is => 'HASH',
            doc => 'Hash of format field type definitions',
            default_value => {},
        },
        info_types => {
            is => 'HASH',
            doc => 'Array of info field type definitions',
            default_value => {},
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    $self->_parse;
    return $self;
}

sub to_string {
    my $self = shift;
    return join("\n", $self->lines);
}

sub _parse {
    my $self = shift;
    my $line_num = 0;
    for my $line ($self->lines) {
        ++$line_num;
        if ($line =~ /^##/) {
            $self->_parse_meta_info($line);
        } elsif ($line =~ /^#CHROM/) {
            my @fields = split("\t", $line);
            confess "Invalid vcf header (line $line_num)" if @fields < 7;
            $self->sample_names([@fields[9..$#fields]]) if @fields >= 9;
        } else {
            confess "Unexpected line in vcf header (line $line_num):\n$line";
        }
    }
}

sub _parse_meta_info {
    my ($self, $line) = @_;
    # if you call me with a line that starts with something
    # other than ##, you are going to have a bad time.
    $line = substr($line, 2);

    my ($key, $value) = split("=", $line, 2);

    # right now we're only parsing lines of the form key=something
    if (!defined $value) {
        push(@{$self->_unparsed_lines}, $line);
        return;
    }

    KEY: {
        $key eq 'fileformat' &&
            do { $self->fileformat($value); last KEY; };

        $key eq 'INFO' && 
            do {
                $self->add_info_str($value);
                last KEY;
            };

        $key eq 'FORMAT' && 
            do {
                $self->add_format_str($value);
                last KEY;
            };

        # Don't know what this is...
        push(@{$self->_unparsed_lines}, $line);
    }
}

# parse types specified in vcf headers like:
# <Key1=Value1,Key2="Value2",...>
sub _parse_vcf_type {
    my $type_str = shift;
    my @fields = $type_str =~
        /^<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="(.*)">$/;
    confess "Malformed vcf type $type_str" unless @fields == 4;

    my $type = {
        id => $fields[0],
        number => $fields[1],
        type => $fields[2],
        description => $fields[3],
    };
    _validate_vcf_type($type);
    return $type;
}

sub _valid_vcf_type_number {
    my $num = shift;
    return $num =~ /[0-9]+/ || exists $VALID_VCF_TYPE_NON_NUMERIC_NUMBERS{$num};
}

sub _validate_vcf_type {
    my $type = shift;
    my $msg;
    $msg = "Vcf type with no ID field: " unless $type->{id};
    $msg = "Vcf type with invalid Number argument" unless _valid_vcf_type_number($type->{number});
    $msg = "Vcf type with invalid Type argument" unless exists $VALID_VCF_TYPE_NAMES{$type->{type}};
    confess "$msg: " . _vcf_type_to_string($type) if $msg;
}

sub _vcf_type_to_string {
    my $type = shift;
    return "<ID=$type->{id},Number=$type->{number},Type=$type->{type},Description=\"$type->{description}\">";
}

sub _index_types {
    my $result = {};
    for my $type (@_) {
        $result->{$type->{id}} = $type;
    }
    return $result;
}

sub add_format_str {
    my ($self, $str) = @_;
    my $type = _parse_vcf_type($str);
    warn "Duplicate entry for FORMAT type $type->{id} in Vcf header" if exists $self->format_types->{$type->{id}};
    $self->format_types->{$type->{id}} = $type;
}

sub add_info_str {
    my ($self, $str) = @_;
    my $type = _parse_vcf_type($str);
    warn "Duplicate entry for INFO type $type->{id} in Vcf header" if exists $self->info_types->{$type->{id}};
    $self->info_types->{$type->{id}} = $type;
}

1;
