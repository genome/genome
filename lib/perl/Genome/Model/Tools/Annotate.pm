package Genome::Model::Tools::Annotate;

use strict;
use warnings;

use Genome;     
use Data::Dumper;

class Genome::Model::Tools::Annotate {
    is => 'Command',
    doc => "annotation tools",
};

sub object_cache_sizes {
    my $class = shift;
    if (@_) {
        my ($low, $high) = @_;
        UR::Context->object_cache_size_lowwater($low);
        UR::Context->object_cache_size_highwater($high);
    }
    else {
        return (UR::Context->object_cache_size_lowwater, UR::Context->object_cache_size_highwater);
    }
}

sub sub_command_sort_position { 15 }

# attributes
sub variant_attributes {
    my $self = shift;
    return (qw/ chromosome_name start stop reference variant /);
}

sub variant_output_attributes {
    my $self = shift;
    return (qw/ type /);
}

sub transcript_attributes {
    my $self = shift;
    return (qw/ gene_name transcript_name transcript_species transcript_source transcript_version strand transcript_status trv_type c_position amino_acid_change ucsc_cons domain all_domains deletion_substructures transcript_error/);
}

sub transcript_report_headers {
    my $self = shift;
    return ($self->variant_attributes, $self->transcript_attributes);
}

# Figures out what the 'type' of this variant should be (snp, dnp, ins, del) based upon
# the start, stop, reference, and variant
# Takes in a variant hash, returns the type
sub infer_variant_type {
    my ($self, $variant) = @_;
    
    if(( (!$variant->{reference})||($variant->{reference} eq '0')||($variant->{reference} eq '-')) &&
        ((!$variant->{variant})||($variant->{variant} eq '0')||($variant->{variant} eq '-'))){
        $self->error_message("Could not determine variant type from variant:");
        $self->error_message(Dumper($variant));
        die;
    }

    # If the start and stop are the same, and ref and variant are defined its a SNP
    if (($variant->{stop} == $variant->{start})&&
        ($variant->{reference} ne '-')&&($variant->{reference} ne '0')&&
        ($variant->{variant} ne '-')&&($variant->{variant} ne '0')) {
        return 'SNP';
    # If start and stop are 1 off, and ref and variant are defined its a DNP
    } elsif (($variant->{stop} - $variant->{start} == 1)&&
             ($variant->{reference})&&($variant->{reference} ne '-')&&($variant->{reference} ne '0')&&
             ($variant->{variant})&&($variant->{variant} ne '-')&&($variant->{variant} ne '0')) {
        return 'DNP';
    # If reference is a dash, we have an insertion
    } elsif (($variant->{reference} eq '-')||($variant->{reference} eq '0')) {
        return 'INS';
    } elsif (($variant->{variant} eq '-')||($variant->{variant} eq '0')) {
        return 'DEL';
    } else {
        $self->error_message("Could not determine variant type from variant:");
        $self->error_message(Dumper($variant));
        die;
    }
}

# Figures out what the 'type' of this variation should be (snp, dnp, ins, del) based upon
# the start, stop, reference, and variation
# Takes in a variation hash, returns the type
sub infer_variation_type {
    my ($self, $variation) = @_;

    # If the start and stop are the same, and ref and variation are defined its a SNP
    if (($variation->stop == $variation->start)&&
        ($variation->reference ne '-')&&($variation->reference ne '0')&&
        ($variation->variant ne '-')&&($variation->variant ne '0')) {
        return 'SNP';
    # If start and stop are 1 off, and ref and variation are defined its a DNP
    } elsif (($variation->stop - $variation->start == 1)&&
             ($variation->reference ne '-')&&($variation->reference ne '0')&&
             ($variation->variant ne '-')&&($variation->variant ne '0')) {
        return 'DNP';
    # If reference is a dash, we have an insertion
    } elsif (($variation->reference eq '-')||($variation->reference eq '0')) {
        return 'INS';
    } elsif (($variation->variant eq '-')||($variation->variant eq '0')) {
        return 'DEL';
    } else {
        $self->error_message("Could not determine variation type from variation:");
        $self->error_message(Dumper($variation));
        die;
    }
}

# Creates output file, deleting old file if necessary
sub _create_file {
    my ($self, $output_file) = @_;
    my $output_fh;

    unlink $output_file if -e $output_file;
    if (-e $output_file) {
        $self->warning_message("found previous output file, removing $output_file");
        unlink $output_file;
        if (-e $output_file) {
            die "failed to remove previous file: $! ($output_file)";
        }
    }
    $output_fh = IO::File->new("> $output_file");
    unless ($output_fh) {
        die "Can't open file ($output_file) for writing: $!";
    }

    return $output_fh;
}

1;

