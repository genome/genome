#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

my @pses = GSC::PSE->get(
    pse_status => 'inprogress',
    ps_id => 3733,
);
my @status;
for my $pse ( sort { $a->id <=> $b->id } @pses ) {
    my ($instrument_data_id) = $pse->added_param('instrument_data_id');
    my ($instrument_data_type) = $pse->added_param('instrument_data_type');
    my $instrument_data = Genome::InstrumentData->get($instrument_data_id);
    my $tgi_lims_status = eval{ 
        return 'NO_INST_DATA' if not $instrument_data; 
        return $instrument_data->attributes(attribute_label => 'tgi_lims_status')->attribute_value;
    };
    $tgi_lims_status ||= 'NA';
    push @status, join(' ', $pse->id, $pse->date_scheduled, $instrument_data_type, $instrument_data_id, $tgi_lims_status);
}

print join("\n", @status)."\n";

exit;

