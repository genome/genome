package Genome::VariantReporting::Report::WithHeader;

use strict;
use warnings FATAL => 'all';
use Genome;
use Set::Scalar;

class Genome::VariantReporting::Report::WithHeader {
    is => 'Genome::VariantReporting::Framework::Component::Report::SingleFile',
    is_abstract => 1,
    has_param => {
        null_character => {
            is => 'Text',
            default => '-',
            doc => 'Character to be used for undefined values in the report',
        },
        delimiter => {
            is => 'Text',
            default => "\t",
            doc => 'Character to be used to separate columns',
        },
        generate_legend_file => {
            is => 'Boolean',
            default => 1,
            doc => 'Boolean flag to control whether or not a legend file will be created',
        }
    },
    has_transient_optional => [
        _legend_fh => {
            is_structural => 1,
        },
    ],
};

sub can_be_merged {
    return 1;
}

sub merge_parameters {
    return {
        sort_columns => [qw(chromosome_name start stop reference variant)],
        contains_header => 1,
    };
}

sub _validation_errors {
    my $self = shift;
    my @errors = $self->SUPER::__errors__;

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

sub legend_file_name {
    my $self = shift;
    return $self->file_name . ".legend";
}

sub legend_path {
    my $self = shift;
    my $path = File::Spec->join($self->output_dir, $self->legend_file_name);

    if (-f $path) {
        return $path;
    } else {
        return;
    }
}

sub initialize {
    my $self = shift;

    my @errors = $self->_validation_errors;
    if (@errors) {
        $self->print_errors(@errors);
        die $self->error_message("Failed to validate");
    }

    $self->SUPER::initialize(@_);
    if ($self->generate_legend_file) {
        my $legend_fh = Genome::Sys->open_file_for_writing(
            File::Spec->join($self->temp_staging_directory,
                $self->legend_file_name));
        $self->_legend_fh($legend_fh);
        $self->write_legend_file();
    }
    $self->print_headers();
}

sub finalize {
    my $self = shift;
    $self->SUPER::finalize(@_);
    if ($self->generate_legend_file) {
        $self->_legend_fh->close;
    }
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

    my @interpreters = values %{$self->interpreters};
    my %available_fields;
    for my $interpreter (@interpreters) {
        for my $field ($interpreter->available_fields) {
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

# Default report method
# Prints the fields in order of the headers.
# Overwrite in child class if different behavior desired.
sub report {
    my $self = shift;
    my $interpretations = shift;

    my %fields = $self->available_fields_dict();
    for my $allele (keys %{$interpretations->{($self->required_interpreters)[0]}}) {
        my @outputs;
        for my $header ($self->headers()) {
            my $interpreter = $fields{$header}->{interpreter};
            my $field = $fields{$header}->{field};

            # If we don't have an interpreter that provides this field, handle it cleanly if the field is known unavailable
            if ($self->header_is_unavailable($header)) {
                push @outputs, $self->_format();
            } elsif ($interpreter) {
                push @outputs, $self->_format($interpretations->{$interpreter}->{$allele}->{$field});
            } else {
                # We use $header here because $field will be undefined due to it not being in an interpreter
                die $self->error_message("Field (%s) is not available from any of the interpreters provided", $header);
            }
        }
        $self->_output_fh->print(join($self->delimiter, @outputs) . "\n");
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
    my $unavailable_headers = Set::Scalar->new($self->unavailable_headers);
    for my $header ($self->headers) {
        # We don't have an interpreter that provides this field
        if ($unavailable_headers->contains($header)) {
            $self->_legend_fh->print(join($self->delimiter, $header, 'undefined') . "\n");
        }
        else {
            my $field = $fields{$header}->{field};
            my $interpreter_name = $fields{$header}->{interpreter};
            my $interpreter = $interpreters->{$interpreter_name};
            $self->_legend_fh->print(join($self->delimiter, $header, $interpreter->field_description($field)) . "\n");
        }
    }

    $self->_legend_fh->print("Filters\n");
    my %filters = %{$self->filters || {}};
    while( my ($filter_name, $filter) = each %filters) {
        next if $self->_skip_list->contains($filter_name);
        $self->_legend_fh->print(join($self->delimiter, $filter_name, $filter->vcf_id, $filter->vcf_description) . "\n");
    }
}

# List of filters that should not be included in the legend file
sub _skip_list {
    return Set::Scalar->new(
        'allele-in-genotype'    #Special filter that is used in all WithHeader reports.
                                #Filters out alternate alleles that are not part of the discovery sample's genotype.
                                #Not useful for analysts to include in legend.
    );
}

1;
