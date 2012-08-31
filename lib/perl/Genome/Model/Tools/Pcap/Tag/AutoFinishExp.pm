package Genome::Model::Tools::Pcap::Tag::AutoFinishExp;

use strict;
use warnings;

use base qw(Genome::Model::Tools::Pcap::Tag);

Genome::Model::Tools::Pcap::Tag::AutoFinishExp->mk_ro_accessors
(qw/
    orientation
    num1 num2 num3 
    chem
    primer_type
    purpose
    fix_cons_errors
    original_cons_errors
    original_single_subclone_bases
    primer
    temp
    id
    exp_id_and_template
    /);

our $VERSION = 0.01;

sub ace
{
    my $self = shift;
    
    return unless defined $self->id;

    my ($ace) = $self->id =~ /^(.+)\.\d+$/;

    return $ace;
}

sub exp_ids_and_templates
{
    my $self = shift;

    return $self->{exp_ids_and_templates} if defined $self->{exp_ids_and_templates};
    
    my %ids_and_templates;
    foreach my $rxn ( split /, /, $self->exp_id_and_template )
    {
        my ($id, $template)  = split / /, $rxn;
        $ids_and_templates{$id} = $template;
    }

    $self->{exp_ids_and_templates} = \%ids_and_templates;
    
    return $self->{exp_ids_and_templates};
}

sub exp_ids
{
    my $self = shift;

    my $ids_and_templates = $self->exp_ids_and_templates;

    return unless defined $ids_and_templates and %$ids_and_templates;
    
    return keys %$ids_and_templates;
}

sub template_for_id
{
    my ($self, $id) = @_;

    die "Need exp id to get template\n" unless defined $id;
    
    my $ids_and_temps = $self->exp_ids_and_templates;

    return unless defined $ids_and_temps and %$ids_and_temps;
    
    return $ids_and_temps->{$id};
}

1;

