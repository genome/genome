package Genome::VariantReporting::Reporter::WithHeader;

use strict;
use warnings FATAL => 'all';
use Genome;
use Memoize qw();

class Genome::VariantReporting::Reporter::WithHeader {
    is => 'Genome::VariantReporting::Framework::Component::Reporter::SingleFile',
    is_abstract => 1,
    has => {
        null_character => {
            is => 'Text',
            default => '-'
        },
        delimiter => {
            is => 'Text',
            default => "\t",
        },
    },
    has_transient_optional => [
        _legend_fh => {},
    ],
};

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__;

    my @sample_names = eval{$self->sample_names};
    my %available_fields = $self->available_fields_dict();
    for my $header ($self->headers) {
        my $error_desc;
        if( defined($available_fields{$header}) and $self->header_is_unavailable($header) ) {
            $error_desc = "Header ($header) is provided from an interpreter but marked as unavailable - remove it from unavailable_headers()";
        } elsif ( !defined($available_fields{$header}) and !$self->header_is_unavailable($header) ) {
            $error_desc = "Interpreter field for header ($header) is not defined. Do you need to overwrite available_fields_dict to provide the correct mapping or add it to unavailable_headers()?";
        }
        if (defined $error_desc) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [],
                desc => $error_desc,
            );
        }
    }

    return @errors;
}

sub requires_interpreters_classes {
    my $self = shift;
    my @interpreters;
    for my $interpreter_name ($self->requires_interpreters) {
        push @interpreters, Genome::VariantReporting::Framework::Factory->get_class('interpreters', $interpreter_name);
    }
    return @interpreters;
}

sub initialize {
    my $self = shift;
    my $output_dir = shift;

    $self->SUPER::initialize($output_dir);
    my $legend_fh = Genome::Sys->open_file_for_writing(File::Spec->join($output_dir, 'legend.tsv'));
    $self->_legend_fh($legend_fh);
    $self->write_legend_file();
    $self->print_headers();
}

sub finalize {
    my $self = shift;
    $self->SUPER::finalize(@_);
    $self->_legend_fh->close;
    return
}

sub headers {
    die "abstract";
}

sub print_headers {
    my $self = shift;

    my @headers = $self->headers();
    $self->_output_fh->print(join($self->delimiter, @headers) . "\n");
}

# Default dictionary that maps headers to interpreter fields
# This can be used if the headers are named exactly the same as the
# interpreters' fields. The interpreter fields need to be unique as well.
# This also works if you are only interested in a subset of the interpreters'
# fields as long as the above requirements are met.
# Overwrite in child class if different behavior desired.
sub available_fields_dict {
    my $self = shift;

    my @interpreters = $self->requires_interpreters_classes;
    my %available_fields;
    for my $interpreter (@interpreters) {
        for my $field ($self->available_fields_for_interpreter($interpreter)) {
            if (defined $available_fields{$field}) {
                die $self->error_message("Fields are not unique. Field: %s, Interpreters: %s and %s",
                    $field, $interpreter->name, $available_fields{$field}->{interpreter});
            }
            $available_fields{$field} = {
                interpreter => $interpreter->name,
                field => $field,
            }
        }
    }
    return %available_fields;
}
Memoize::memoize('available_fields_dict');

sub available_fields_for_interpreter {
    my $self = shift;
    my $interpreter = shift;

    return $interpreter->available_fields();
}

# Default report method
# Prints the fields in order of the headers.
# Overwrite in child class if different behavior desired.
sub report {
    my $self = shift;
    my $interpretations = shift;

    my %fields = $self->available_fields_dict();
    for my $allele (keys %{$interpretations->{($self->requires_interpreters)[0]}}) {
        for my $header ($self->headers()) {
            my $interpreter = $fields{$header}->{interpreter};
            my $field = $fields{$header}->{field};

            # If we don't have an interpreter that provides this field, handle it cleanly if the field is known unavailable
            if ($self->header_is_unavailable($header)) {
                $self->_output_fh->print( $self->_format() . "\t");
            } elsif ($interpreter) {
                $self->_output_fh->print($self->_format($interpretations->{$interpreter}->{$allele}->{$field}) . "\t");
            } else {
                # We use $header here because $field will be undefined due to it not being in an interpreter
                die $self->error_message("Field (%s) is not available from any of the interpreters provided", $header);
            }
        }
        $self->_output_fh->print("\n");
    }
}

sub _format {
    my $self = shift;
    my $string = shift;

    if (defined $string ) {
        return $string;
    }
    else {
        return $self->null_character;
    }
}

sub header_is_unavailable {
    my ($self, $header) = @_;

    if (grep {$header eq $_} $self->unavailable_headers) {
        return 1;
    } else {
        return 0;
    }
}

# Override this method if some fields are not available from interpreters but should be in the report as null values
sub unavailable_headers {
    return;
}

sub write_legend_file {
    my $self = shift;

    my %fields = $self->available_fields_dict;
    my $interpreters = $self->interpreters;
    $self->_legend_fh->print("Headers\n");
    for my $header ($self->headers) {
        my $field = $fields{$header}->{field};
        my $interpreter_name = $fields{$header}->{interpreter};
        my $interpreter = $interpreters->{$interpreter_name};
        $self->_legend_fh->print(join($self->delimiter, $header, $interpreter->field_description($field)) . "\n");
    }

    $self->_legend_fh->print("Filters\n");
    my %filters = %{$self->filters || {}};
    while( my ($filter_name, $filter) = each %filters) {
        next if $filter_name eq 'allele-in-genotype';
        $self->_legend_fh->print(join($self->delimiter, $filter_name, $filter->vcf_id, $filter->vcf_description) . "\n");
    }
}

1;
