package Genome::Config::Command::ConfigureQueuedInstrumentData;

use strict;
use warnings;

use Genome;
use Genome::Config::Command::ConfigureQueuedInstrumentData; #load real module first

use Sub::Install;

Sub::Install::reinstall_sub ({
    as => '_import_instrument_data',
    code => \&_import_lims_data,

});

sub _import_lims_data {
    my $self = shift;
    my $current_pair = shift;

    my $importer = Genome::Site::TGI::Command::ImportDataFromLims->create(
        instrument_data => $current_pair->instrument_data,
        analysis_project => $current_pair->analysis_project,
    );

    unless($importer->execute) {
        $self->fatal_message('Failed to import data from LIMS.');
    }
}
