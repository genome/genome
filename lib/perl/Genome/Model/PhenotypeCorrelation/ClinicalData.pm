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
        subject_column_header => {
            is => "Text",
            doc => "The name of the subject id column",
        },
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

sub attribute_values_for_samples {
    my ($self, $attr, @sample_names) = @_;
    $self->_validate_attribute($attr);
    my @sample_indices = @{$self->_sample_indices}{@sample_names};
    return @{$self->_phenotypes->{$attr}}[@sample_indices];
}

sub sample_names {
    my $self = shift;
    return keys %{$self->_sample_indices};
}

sub from_database {
    my ($class, $nomenclature, @samples) = @_;
    my $self = $class->create;

    @samples = sort { ncmp($a->name, $b->name) } @samples;
    my %sample_indices;
    @sample_indices{map {$_->name} @samples} = 0..$#samples;
    $self->_sample_indices(\%sample_indices);

    my %phenotypes;
    my %attributes;
    foreach my $sample (@samples) {
        my @attrs = $sample->source->attributes_for_nomenclature($nomenclature);
        my $sample_idx = $sample_indices{$sample->name};

        foreach my $attr (@attrs) {
            my $label = $attr->attribute_label;
            $attributes{$label} = {} unless defined $attributes{$label};
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
    my ($class, $path, %params) = @_;
    return $class->from_filehandle(Genome::Sys->open_file_for_reading($path), %params);
}

sub from_filehandle {
    my ($class, $fh, %params) = @_;
    my $missing_string = $params{missing_string} || "NA";
    my $self = $class->create;

    my $header_line = $fh->getline or confess "Failed to read clinical data header";
    chomp $header_line;
    my @attr_names = split("\t", $header_line);
    my $subject_column_header = shift @attr_names;

    my @sample_data;
    my $line_num = 1;
    while (my $line = $fh->getline) {
        chomp $line;
        ++$line_num;

        my @fields = map {
                my $val = $_;
                # replace NA or '' with undef
                scalar(grep {$val eq $_} ('', $missing_string)) > 0
                    ? undef
                    : $_
            } split("\t", $line, -1);
            # the LIMIT=-1 arg to split here makes sure trailing undef: are detected
            # (the case where the line ends with tab)

        if (scalar @fields - 1 != scalar @attr_names) { # - 1 for Sample_name
            confess "At line $line_num: # of columns does not match header:\n"
                . "HEADER: " . join(",", ("Sample_name", @attr_names)) . "\n"
                . "  LINE: " . join(",", @fields) . "\n";
        }
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
    $self->subject_column_header($subject_column_header);

    return $self;
}

sub to_file {
    my ($self, $path, %params) = @_;
    return $self->to_filehandle(Genome::Sys->open_file_for_writing($path), %params);
}

sub to_filehandle {
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

sub convert_attr_to_factor {
    my ($self, $attribute, %params) = @_;
    my $levels = $params{levels};
    if (!defined $levels || scalar @$levels < 2) {
        confess "Required parameter levels not supplied or contains <2 values.";
    }
    my $values = $self->_phenotypes->{$attribute};
    confess "Unknown attribute '$attribute'" unless defined $values;


    # check if this is already categorical data
    my %invmapping;
    @invmapping{0..$#$levels} = @$levels;
    return unless grep {defined $_ && !exists $invmapping{$_}} @$values;

    # if we need to do conversion, then we had better not allow undef as a level
    # since undef already has a special meaning (NA/missing data)
    confess "Attempted to coerce attribute '$attribute' to a factor with undef values"
        if grep {!defined $_} @$levels;

    # compute mapping of level names => indices
    my %mapping;
    @mapping{@$levels} = 0..$#$levels;

    my %uniq = map {$_ => 0} grep {$_} @$values;
    if (grep {!exists $mapping{$_}} keys %uniq) {
        confess "Unable to coerce attribute '$attribute' to a factor:\n"
            . "Supplied levels:\n\t" . join("\n\t", @$levels) . "\n"
            . "Actual values:\n\t" . join("\n\t", keys %uniq) . "\n"
    }

    $self->_phenotypes->{$attribute} = [ map {
        (defined $_ && defined $mapping{$_})
            ? $mapping{$_}
            : $_
        } @$values ];

    return %mapping;
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
