package Genome::VariantReporting::Framework::Component::Interpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Set::Scalar;
use Carp qw(confess);

class Genome::VariantReporting::Framework::Component::Interpreter {
    is => ['Genome::VariantReporting::Framework::Component::Base', 'Genome::VariantReporting::Framework::Component::WithTranslatedInputs'],
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

1;
