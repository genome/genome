package Genome::VariantReporting::Framework::Component::Interpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Set::Scalar;
use Carp qw(confess);

class Genome::VariantReporting::Framework::Component::Interpreter {
    is => ['Genome::VariantReporting::Framework::Component::WithTranslatedInputs'],
    is_abstract => 1,
};

sub name {
    my $class = shift->class;
    confess "Abstract method 'name' must be defined in subclass '$class'";
}

sub requires_annotations {
    my $class = shift->class;
    confess "Abstract method 'requires_annotations' must be defined in subclass '$class'";
}

sub _interpret_entry {
    my $class = shift->class;
    confess "Abstract method '_interpret_entry' must be defined in subclass '$class'";
}

sub interpret_entry {
    my ($self, $entry, $alt_alleles) = @_;

    my $allele_list   = Set::Scalar->new(@{$entry->{alternate_alleles}});
    my $input_alleles = Set::Scalar->new(@$alt_alleles);

    unless ($input_alleles->is_subset($allele_list)) {
        confess "The input alleles: $input_alleles is not subset of vcf alt_allele list: $allele_list";
    }

    my %interpret_entry = $self->_interpret_entry($entry, $alt_alleles);
    my $output_alleles = Set::Scalar->new(keys %interpret_entry); 

    unless ($output_alleles->is_equal($input_alleles)) {
        confess "The output allele list: $output_alleles is not the same as the input: $input_alleles";
    }
    
    return %interpret_entry;
}

# Set all interpretation fields to null values.
# Override this method if special behavior is needed.
sub null_interpretation {
    my ($self, $alt_alleles) = @_;
    my %return_values;

    for my $allele (@$alt_alleles) {
        $return_values{$allele} = { map { $_ => $self->interpretation_null_character } $self->available_fields };
    }

    return %return_values;
}

sub interpretation_null_character {
    return '.';
}

sub field_description {
    my ($self, $field) = @_;

    my %field_descriptions = $self->field_descriptions;
    return $field_descriptions{$field};
}

sub field_descriptions {
    my $class = shift->class;
    confess "Abstract method 'field_descriptions' must be defined in subclass '$class'. This is a hash of available_fields mapped to their respective description";
}

sub available_fields {
    my $self = shift;
    my %field_descriptions = $self->field_descriptions;
    return keys %field_descriptions;
}

sub vr_doc_sections {
    my $self = shift;
    my @sections = $self->SUPER::vr_doc_sections;
    if ($self->requires_annotations) {
        push @sections,
            {
                header => "REQUIRED EXPERTS",
                items => [$self->requires_annotations],
            };
    }
    my @items;
    my %descriptions = $self->field_descriptions;
    while (my ($name, $description) = each %descriptions) {
        push @items, sprintf
        (
            "  %s\n%s",
            Term::ANSIColor::colored($name, 'bold'),
            Text::Wrap::wrap(
                "    ",
                "    ",
                $description,
            ),
        );
    }
    push @sections,
        {
            header => "AVAILABLE FIELDS",
            items => \@items
        };
    return @sections;
}

1;
