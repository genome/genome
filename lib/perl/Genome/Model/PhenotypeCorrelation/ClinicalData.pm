package Genome::Model::PhenotypeCorrelation::ClinicalData;

use Digest::MD5;
use Genome;
use Sort::Naturally qw/nsort ncmp/;
use Data::Dumper;
use Carp qw/confess/;

use strict;
use warnings;

my $CATEGORICAL_THRESHOLD = 0.05;

class Genome::Model::PhenotypeCorrelation::ClinicalData {
    is => "UR::Object",
    doc => "Clinical (phenotype) data for a set of samples",
    has_transient_optional => [
        _attributes => {
            is => "HASH",
            doc => "Contains type information (datatype, quantitative vs categorical, etc.)",
        },
        _sample_indices => {
            is => "HASH",
        },
        _phenotypes => {
            is => "HASH",
        },
    ],
};

sub sample_count {
    my $self = shift;
    return scalar keys %{$self->_sample_indices};
}

sub attributes {
    my ($self, %params) = @_;
    if (defined $params{categorical}) {
        return nsort grep {
                $self->_attributes->{$_}->{categorical} == $params{categorical}
            } keys %{$self->_attributes};
    }
    return nsort keys %{$self->_attributes};
}

sub attribute_types {
    my $self = shift;
    return $self->_attributes;
}

sub _validate_attribute {
    my ($self, $attr) = @_;
    unless (defined $self->_phenotypes->{$attr}) {
        my $valid_attrs = join("\n\t", $self->attributes);
        confess "Unknown attribute name '$attr'. Valid attributes are\n\t$valid_attrs";
    }
}

sub attribute_values {
    my ($self, $attr) = @_;
    $self->_validate_attribute($attr);
    return $self->_phenotypes->{$attr};
}

sub attribute_value_for_sample {
    my ($self, $sample_name, $attr) = @_;
    $self->_validate_attribute($attr);
    my $sample_idx = $self->_sample_indices->{$sample_name};
    return $self->_phenotypes->{$attr}->[$sample_idx];
}

sub from_database {
    my ($class, $nomenclature, @samples) = @_;
    my $self = $class->create;

    @samples = sort { ncmp($a->name, $b->name) } @samples;
    my %sample_indices;
    @sample_indices{map {$_->name} @samples} = 0..$#samples;
    $self->_sample_indices(\%sample_indices);

    $DB::single=1;
    my %phenotypes;
    my %attributes;
    foreach my $sample (@samples) {
        my @attrs = $sample->source->attributes_for_nomenclature($nomenclature);
        my $sample_idx = $sample_indices{$sample->name};

        foreach my $attr (@attrs) {
            my $label = $attr->attribute_label;
            $attributes{$label} = {};

            $phenotypes{$label} = [(undef) x scalar @samples] unless defined $phenotypes{$label};
            $phenotypes{$label}->[$sample_idx] = $attr->attribute_value;
        }
    }
    $self->_phenotypes(\%phenotypes);

    confess "No clinical attributes found!" unless scalar(keys %attributes);

    for my $attr (keys %attributes) {
        $attributes{$attr}->{categorical} = $self->_is_categorical($attr);
    }
    $self->_attributes(\%attributes);

    return $self;
}

sub from_file {
    my ($class, $fh, %params) = @_;
    my $self = $class->create;
    my $missing_string = $params{missing_string} || "NA";

    my $header_line = $fh->getline or confess "Failed to read clinical data header";
    chomp $header_line;
    my @attr_names = split("\t", $header_line);
    shift @attr_names; # drop the leading Sample_name in header

    my @sample_data;
    while (my $line = $fh->getline) {
        chomp $line;
        my @fields = map { $_ eq $missing_string ? undef : $_ } split("\t", $line);
        push(@sample_data, \@fields);
    }

    @sample_data = sort { ncmp($a->[0], $b->[0]) } @sample_data;
    my %sample_indices;
    @sample_indices{map {$_->[0]} @sample_data} = 0..$#sample_data;

    my %phenotypes;
    @phenotypes{@attr_names} = map { [(undef) x scalar(@sample_data)] } @attr_names;
    for my $i (0..$#sample_data) {
        my @values = @{$sample_data[$i]};
        my $sample_name = shift @values;
        for my $j (0..$#attr_names) {
            $phenotypes{$attr_names[$j]}->[$i] = $values[$j];
        }
    }

    $self->_sample_indices(\%sample_indices);
    $self->_phenotypes(\%phenotypes);
    my %attributes = map { $_ => { categorical => $self->_is_categorical($_) } } @attr_names;
    $self->_attributes(\%attributes);

    return $self;
}

sub write_file {
    my ($self, $fh, %params) = @_;
    my $missing_string = $params{missing_string} || "NA";
    my @attr_names;
    if (defined $params{attribute_names}) {
        @attr_names = @{$params{attribute_names}};
        for my $attr (@attr_names) {
            confess "unknown clinical attribute $attr" unless defined $self->_attributes->{$attr};
        }
    } else {
        @attr_names = nsort keys %{$self->_attributes};
    }

    my $md5 = new Digest::MD5;
    my $write = sub {
        my $data = shift;
        $md5->add($data);
        $fh->write($data);
    };

    &$write("Sample_name\t" . join("\t", @attr_names) . "\n");
    for my $sample_name (nsort keys %{$self->_sample_indices}) {
        my $sample_idx = $self->_sample_indices->{$sample_name};
        my @values = map {
             defined $self->_phenotypes->{$_}->[$sample_idx] ?
                $self->_phenotypes->{$_}->[$sample_idx] : $missing_string
            } @attr_names;
        &$write("$sample_name\t" . join("\t", @values) . "\n");
    }

    return $md5->hexdigest;
}

sub _is_categorical {
    my ($self, $attribute) = @_;
    my $values = $self->_phenotypes->{$attribute};
    confess "Unknown attribute '$attribute'" unless defined $values;
    my %uniq = map {$_ => 1} grep {$_} @$values;
    my $n_uniq = scalar keys(%uniq);
    my $ratio = $n_uniq / $self->sample_count;
    my $decision = ($n_uniq <= 10) || ($ratio <= $CATEGORICAL_THRESHOLD);
    return $decision;
}


1;
