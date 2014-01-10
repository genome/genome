package Genome::InstrumentData::Gatk::BaseWithKnownSites;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Gatk::BaseWithKnownSites {
    is => 'Genome::InstrumentData::Gatk::Base',
    has_optional_input => [
        known_sites => {# PROVIDES known_sites_vcf
            is => 'Genome::Model::Build::ImportedVariationList',
            is_many => 1,
        },
    ],
};

sub known_sites_vcfs {
    my $self = shift;
    my @vcfs = @{$self->known_sites_indel_vcfs};
    return if not @vcfs;
    push @vcfs, @{$self->known_sites_snv_vcfs};
    return \@vcfs;
}

sub known_sites_indel_vcfs {
    my $self = shift;
    return $self->_known_sites_vcfs('indel');
}

sub known_sites_snv_vcfs {
    my $self = shift;
    return $self->_known_sites_vcfs('snv');
}

sub _known_sites_vcfs {
    my ($self, $type) = @_;

    Carp::confess('No type given to get known sites vcfs!') if not $type;

    if ( not $self->{_known_sites_vcfs} ) {
        my %known_sites_vcfs = $self->link_known_sites_vcfs;
        return if not %known_sites_vcfs;
        $self->{_known_sites_vcfs} = \%known_sites_vcfs;
    }

    return [ @{$self->{_known_sites_vcfs}->{$type}} ];
}

sub link_known_sites_vcfs {
    my $self = shift;

    my %known_sites_vcfs = ( indel => [], snv => [] );
    for my $known_site ( $self->known_sites ) { # all have indels for now...
        # indel
        my $indel_result = $known_site->indel_result;
        if ( $indel_result ) {
            #$self->error_message('No indel result for known site! '.$known_site->__display_name__);
            #return;
            my $link_name = $self->_get_and_link_vcf_from_known_site_result($indel_result, 'indel');
            return if not $link_name;
            push @{$known_sites_vcfs{indel}}, $link_name;
        }

        # snv
        my $snv_result = $known_site->snv_result;
        if ( $snv_result ) {
            my $link_name = $self->_get_and_link_vcf_from_known_site_result($snv_result, 'snv');
            return if not $link_name;
            push @{$known_sites_vcfs{snv}}, $link_name;
        }
    }

    return %known_sites_vcfs;
}

sub _get_and_link_vcf_from_known_site_result {
    my ($self, $result, $type) = @_;
    $self->debug_message("Link VCF from known site result...");

    Carp::confess('No result given to get and link vcf from known site result!') if not $result;
    Carp::confess('No type given to get and link vcf from known site result!') if not $type;

    $self->debug_message('Known site result: '.$result->__display_name__);
    $self->debug_message("Type: $type");

    my $source_known_sites_vcf = $result->path.'/'.$type.'s.hq.vcf';
    if ( not -s $source_known_sites_vcf ) {
        $self->error_message("No $type vcf (${type}s.hq.vcf) in indel result output directory! ".$result->output_dir);
        return;
    }
    $self->debug_message("Target: $source_known_sites_vcf");

    my $link_name = $self->_tmpdir.'/'.$result->id.'.vcf';
    $self->debug_message("Link: $link_name");

    symlink($source_known_sites_vcf, $link_name) if not -l $link_name;

    $self->debug_message("Link VCF from known site result...done");
    return $link_name;
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $link_known_sites_vcfs = $self->link_known_sites_vcfs;
    if ( not $link_known_sites_vcfs ) { # undef on error
        $self->delete;
        return;
    }

    for my $known_sites ( $self->known_sites ) {
        $self->debug_message('Known sites: '.$known_sites->id);
    }

    return $self;
}

1;

