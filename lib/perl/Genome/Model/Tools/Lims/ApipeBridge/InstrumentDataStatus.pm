package Genome::Model::Tools::Lims::ApipeBridge::InstrumentDataStatus;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Lims::ApipeBridge::InstrumentDataStatus { 
    is => 'Command::V2',
    has => [
        summary_only => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Only show the summary of the instrument data status.',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            is_optional => 1,
            doc => 'Show the status of these instrument data.',
        },
    ],
};

sub help_brief { return 'Show LIMS-APIPE status of instrument data'; }
sub help_detail { return help_brief(); }

sub execute {
    my $self = shift;

    my ($new_instrument_data, $failed_instrument_data, $qidfgms) = $self->_resolve_instrument_data_and_qidfgms;

    # Map QIDFGMs to instrument data
    my $join = ' ';
    my ($status, %totals, %statuses);
    for my $qidfgm ( @$qidfgms ) {
        my ($instrument_data_id) = $qidfgm->added_param('instrument_data_id');
        my ($instrument_data_type) = $qidfgm->added_param('instrument_data_type');
        $instrument_data_type =~ s/\s+/_/g;
        if ( not $instrument_data_id ) {
            $self->warning_message('No instrument data id for QIDFGM! '.$qidfgm->id);
            next;
        }
        my $instrument_data;
        my $instrument_data_status;
        if ( $instrument_data = delete $new_instrument_data->{$instrument_data_id} ) {
            $totals{synced}++;
            $instrument_data_status = 'new';
        }
        elsif ( $instrument_data = delete $failed_instrument_data->{$instrument_data_id} ) {
            $totals{synced}++;
            $instrument_data_status = 'failed';
        }
        else {
            my $instrument_data = Genome::InstrumentData->get($instrument_data_id);
            if ( not $instrument_data ) { # OK
                $totals{na}++;
                $instrument_data_status = 'na';
            }
        }
        $statuses{ $instrument_data_type.'('.$instrument_data_status.')' }++;
        $statuses{total}++;
        $status .= join(
            ' ', 
            $instrument_data_id, $qidfgm->id, $instrument_data_type, $instrument_data_status,
        )."\n";
    }
    $totals{qidfgm} = @$qidfgms;
    $totals{missing} = keys(%$new_instrument_data) + keys(%$failed_instrument_data);
    $totals{inprogress} = $totals{qidfgm} - $totals{missing};

    print join($join, (qw/ ID QIDFGM TYPE STATUS /))."\n$status" if not $self->summary_only and $status;
    print join(' ', 'STATUS:', map { $_.'='.$statuses{$_} } sort keys %statuses)."\n";
    print join(' ', 'TOTALS:', map { $_.'='.$totals{$_} } sort keys %totals)."\n";

    return 1;
}

sub _resolve_instrument_data_and_qidfgms {
    my $self = shift;

    $instrument_data_params{id} = [ map { $_->id } $self->instrument_data ] if $self->instrument_data;

    # Get 'new' instrument data
    my %new_instrument_data = map { $_->id => $_ } Genome::InstrumentData->get(
        %instrument_data_params,
        'analysis_project_bridges.status' => [qw/ new /],
    );

    # Get 'failed' instrument data
    my %failed_instrument_data = map { $_->id => $_ } Genome::InstrumentData->get(
        %instrument_data_params,
        'analysis_project_bridges.status' => [qw/ failed /],
    );
    # Get the inprogress QIDFGMs mapped to instrument data
    my %qidfgm_params = (
        ps_id => 3733,
        pse_status => 'inprogress',
    );
    $qidfgm_params{pse_id} = [
        map { $_->pse_id }
        GSC::PSEParam->get(
            param_name => 'instrument_data_id',
            param_value => [ map { $_->id } $self->instrument_data ],
        )
    ] if $self->instrument_data;

    my @qidfgms = GSC::PSE->get(
        %qidfgm_params
    );

    return (\%new_instrument_data, \%failed_instrument_data, \@qidfgms);
}

1;

