package Genome::Site::TGI::Command::ImportDataForWorkOrder;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Command::ImportDataForWorkOrder {
    is => 'Command::V2',
    has_input => [
        work_order_id => {
            is => 'Text',
            doc => 'The ID of the work order to import',
        },
    ],
};

sub execute {
    my $self = shift;

    my $wo = $self->_resolve_work_order;
    $self->_resolve_samples;

    my @items = $self->_resolve_instrument_data;
    my @existing = Genome::InstrumentData->get(id => [map $_->entity_id, @items]);

    $self->status_message(
        'Work Order has %s instrument data of which %s are already present.',
        scalar(@items),
        scalar(@existing),
    );

    my %found;
    $found{$_->id}++ for @existing;
    my @to_import = grep { !$found{$_->entity_id} } @items;

    my @ii = Genome::Site::TGI::Synchronize::Classes::IndexIllumina->get(id => [map $_->entity_id, @to_import]);
    $self->_import_indexillumina($_) for @ii;
    $self->_import_anp_associations(@ii);

    my @g = Genome::Site::TGI::Synchronize::Classes::Genotyping->get(id => [map $_->entity_id, @to_import]);
    $self->_import_genotyping($_) for @g;
    $self->_import_anp_associations(@g);

    return 1;
}

sub _resolve_work_order {
    my $self = shift;
    my $wo_id = $self->work_order_id;

    my $existing_project = Genome::Project->get($wo_id);
    unless ($existing_project) {
        my $wo = Genome::Site::TGI::Synchronize::Classes::LimsProject->get($wo_id);
        unless ($wo) {
            $self->fatal_message('No work order found for ID %s', $wo_id);
        }

        $existing_project = $wo->create_in_genome;
    }

    return $existing_project;
}

sub _resolve_instrument_data {
    my $self = shift;

    my @existing_items = Genome::ProjectPart->get(
        label => 'instrument_data',
        project_id => $self->work_order_id,
    );

    my @importable_items = Genome::Site::TGI::Synchronize::Classes::LimsProjectInstrumentData->get(
        project_id => $self->work_order_id,
    );

    my %found;
    $found{$_->entity_id}++ for @existing_items;
    my @to_import = grep { !$found{$_->entity_id} } @importable_items;

    for my $item (@to_import) {
        push @existing_items, $item->create_in_genome;
    }

    return @existing_items;
}

sub _resolve_samples {
    my $self = shift;

    my @existing_items = Genome::ProjectPart->get(
        label => 'sample',
        project_id => $self->work_order_id,
    );

    my @importable_items = Genome::Site::TGI::Synchronize::Classes::LimsProjectSample->get(
        project_id => $self->work_order_id,
    );

    my %found;
    $found{$_->entity_id}++ for @existing_items;
    my @to_import = grep { !$found{$_->entity_id} } @importable_items;

    for my $item (@to_import) {
        push @existing_items, $item->create_in_genome;
    }

    my @existing_samples = Genome::Sample->get(id => [map $_->entity_id, @existing_items]);

    my %found_samples;
    $found_samples{$_->id}++ for @existing_samples;
    my @samples_to_import = grep { !$found_samples{$_->entity_id} } @existing_items;
    my @os = Genome::Site::TGI::Synchronize::Classes::OrganismSample->get(id => [map $_->entity_id, @samples_to_import]);
    $self->_import_organismsample($_) for @os;

    return @existing_items;
}

sub _import_indexillumina {
    my $self = shift;
    my $ii = shift;

    my $existing_library = Genome::Library->get($ii->library_id);
    unless ($existing_library) {
        my $ls = Genome::Site::TGI::Synchronize::Classes::LibrarySummary->get($ii->library_id);
        $self->_import_librarysummary($ls);
    }

    $ii->create_in_genome;
}

sub _import_librarysummary {
    my $self = shift;
    my $ls = shift;

    my $existing_sample = Genome::Sample->get($ls->sample_id);
    unless ($existing_sample) {
        my $os = Genome::Site::TGI::Synchronize::Classes::OrganismSample->get($ls->sample_id);
        $self->_import_organismsample($os);
    }

    $ls->create_in_genome;
}

sub _import_organismsample {
    my $self = shift;
    my $os = shift;

    my $existing_individual = Genome::Individual->get($os->source_id);
    unless ($existing_individual) {
        my $oi = Genome::Site::TGI::Synchronize::Classes::OrganismIndividual->get($os->source_id);
        $self->_import_organismindividual($oi);
    }

    $os->create_in_genome;
}

sub _import_organismindividual {
    my $self = shift;
    my $oi = shift;

    my $existing_taxon = Genome::Taxon->get($oi->taxon_id);
    unless ($existing_taxon) {
        my $ot = Genome::Site::TGI::Synchronize::Classes::OrganismTaxon->get($oi->taxon_id);
        $self->_import_organismtaxon($ot);
    }

    $oi->create_in_genome;
}

sub _import_organismtaxon {
    my $self = shift;
    my $ot = shift;

    $ot->create_in_genome;
}

sub _import_anp_associations {
    my $self = shift;
    my @data = @_;

    return 1 unless @data;

    my @existing = Genome::Config::AnalysisProject::InstrumentDataBridge->get(instrument_data_id => [map $_->id, @data]);
    my @potential = Genome::Site::TGI::Synchronize::Classes::InstrumentDataAnalysisProjectBridge->get(instrument_data_id => [map $_->id, @data]);

    my %existing_lookup;
    for (@existing) {
        $existing_lookup{$_->instrument_data_id}{$_->analysis_project_id}++
    }

    for (@potential) {
        unless ($existing_lookup{$_->instrument_data_id}{$_->analysis_project_id}) {
            $_->create_in_genome;
        }
    }

    return 1;
}

sub _import_genotyping {
    my $self = shift;
    my $g = shift;

    my $existing_sample = Genome::Sample->get($g->sample_id);
    unless ($existing_sample) {
        my $os = Genome::Site::TGI::Synchronize::Classes::OrganismSample->get($g->sample_id);
        $self->_import_organismsample($os);
    }

    $g->create_in_genome;
}

1;
