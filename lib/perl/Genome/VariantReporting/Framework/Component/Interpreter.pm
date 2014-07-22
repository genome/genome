package Genome::VariantReporting::Framework::Component::Interpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Carp qw(confess);

class Genome::VariantReporting::Framework::Component::Interpreter {
    is => ['Genome::VariantReporting::Framework::Component::Base', 'Genome::VariantReporting::Framework::Component::WithTranslatedInputs'],
    is_abstract => 1,
};

sub name {
    my $self = shift;
    my $class = $self->class;
    confess "Abstract method 'name' must be defined in class '$class'";
}

sub requires_annotations {
    my $self = shift;
    my $class = $self->class;
    confess "Abstract method 'requires_annotations' must be defined in class '$class'";
}

sub interpret_entry {
    my ($self, $entry, $alt_alleles) = @_;

    my @allele_list   = sort @{$entry->{alternate_alleles}};
    my @input_alleles = sort @$alt_alleles;
    my $input_alleles = join ',', @input_alleles;

    for my $input_allele (@input_alleles) {
        unless (grep{$_ eq $input_allele}@allele_list) {
            confess "The input allele $input_allele is not in vcf alt_alleles column";
        }
    }

    my %interpret_entry = $self->process_interpret_entry($entry, $alt_alleles);
    my @output_alleles = sort keys %interpret_entry; 
    my $output_alleles = join ',', @output_alleles;

    unless (@input_alleles ~~ @output_alleles) {
        confess "The output allele list: $output_alleles is not the same as the input: $input_alleles";
    }
    
    return %interpret_entry;
}

1;
