package Genome::Model::Tools::Lims::ApipeBridge::FixPidfaParamsForBase;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Lims::ApipeBridge::FixPidfaParamsForBase { 
    is => 'Command::V2',
    is_abstract => 1,
    has_optional => [
        pidfa_id => {
            is => 'Integer',
            doc => 'Id of PIDFA to fix params.',
        },
        prior_id => {
            is => 'Integer',
            doc => 'Id of prior PSE to get and fix PIDFA params.',
        },
        qidfgm_id => {
            is => 'Integer',
            doc => 'Id of QIDFGM to get and fix PIDFA params.',
        },
    ],
};

sub pidfa_ps_id { return 3870; }
sub qidfgm_ps_id { return 3733; }

sub help_brief { return 'Fix PIDFA PSE params for '.$_[0]->instrument_data_type; }
sub help_detail { return 'Given a PIDFA or PIDFA\'s prior PSE, this command can fix some broken '.$_[0]->instrument_data_type.' PIDFAs. After fixing, schedule the pse by using "sw --sched ${PIDFA_ID}"'; }

sub _starting_points { 
    my $self = shift;
    my @starting_points = (qw/ prior_id pidfa_id qidfgm_id /);
    push @starting_points, $self->_additional_starting_points if $self->can('_additional_starting_points');
    return @starting_points;
}

sub execute {
    my $self = shift;
    $self->status_message('Fix PIDFA params for '.$self->instrument_data_type.'...');

    my $starting_point_method = $self->_get_init_method;
    return if not $starting_point_method;

    my ($prior, $pidfa, $qidfgm) = $self->$starting_point_method; # qidfgm is optional!
    return if not $pidfa;

    $self->status_message('PIDFA: '.$pidfa->id);
    $self->status_message('QIDFGM: '.( $qidfgm ? $qidfgm->id : 'NA'));
    $self->status_message('Prior PSE id: '.$prior->id);
    $self->status_message('Prior process: '.$prior->process_to);

    if ( not grep { $prior->process_to eq $_ } $self->valid_prior_processes ) {
        $self->error_message('Invalid prior process for '.$self->instrument_data_type.'!');
        next;
    }

    my $sequence_item = $self->_get_sequence_item_from_prior($prior);
    return if not $sequence_item;
    $self->status_message('Instrument data: '.$sequence_item->id);

    my %params_to_fix = (
        instrument_data_id => $sequence_item->id,
        instrument_data_type => $self->instrument_data_type,
    );
    if ( $qidfgm ) {
        $self->status_message('FIX QIDFGM params...');
        my $fixed = $self->_fix_params_for($qidfgm, \%params_to_fix);
        if ( not $fixed ) {
            $self->error_message('Failed to fix params for QIDFGM!');
            return;
        }
    }

    $self->_additional_params_to_fix(\%params_to_fix, $sequence_item);
    $self->status_message('FIX PIDFA params...');
    my $fixed = $self->_fix_params_for($pidfa, \%params_to_fix);
    if ( not $fixed ) {
        $self->error_message('Failed to fix params for PIDFA!');
        return;
    }

    $self->status_message('Done');
    return 1;
}

sub _get_init_method {
    my $self = shift;

    my @starting_points = grep { defined $self->$_ } $self->_starting_points;
    if ( not @starting_points ) {
        $self->error_message('No starting point indicated! Select from '.join(', ', @starting_points));
        return;
    }
    elsif ( @starting_points > 1 ) {
        $self->error_message('Multiple starting points indicated! Plese select only one.');
        return;
    }

    return '_init_with_'.$starting_points[0];
}

sub _init_with_prior_id {
    my $self = shift;

    my $prior = GSC::PSE->get(id => $self->prior_id);
    if ( not $prior ) {
        $self->error_message('Failed to get prior PSE for id!'. $self->prior_id);
        return;
    }

    my $pidfa = $self->_get_future_for_pse($prior, pidfa_ps_id());
    if ( not $pidfa ) {
        $self->error_message('Failed to get PIDFA for prior PSE! '.$prior->id);
        return;
    }

    my $qidfgm = $self->_get_future_for_pse($pidfa, qidfgm_ps_id());

    return ($prior, $pidfa, $qidfgm);
}

sub _init_with_pidfa_id {
    my $self = shift;

    my $pidfa = GSC::PSE->get(id => $self->pidfa_id, ps_id => pidfa_ps_id());
    if ( not $pidfa ) {
        $self->error_message('Failed to get PIDFA for id! '.$self->pidfa_id);
        return;
    }

    my $prior = $self->_get_prior_for_pse($pidfa);
    if ( not $prior ) {
        $self->error_message('Failed to get prior for PIDFA! '.$pidfa->id);
        return;
    }

    my $qidfgm = $self->_get_future_for_pse($pidfa, qidfgm_ps_id());
    
    return ($prior, $pidfa, $qidfgm);
}

sub _init_with_qidfgm_id {
    my $self = shift;

    my $qidfgm = GSC::PSE->get(id => $self->qidfgm_id, ps_id => qidfgm_ps_id());
    if ( not $qidfgm ) {
        $self->error_message('Failed to get QIDFGM for id! '.$self->qidfgm_id);
        return;
    }

    my $pidfa = $self->_get_prior_for_pse($qidfgm, pidfa_ps_id());
    if ( not $pidfa ) {
        $self->error_message('Failed to get PIDFA for QIDFGM! '.$qidfgm->id);
        return;
    }

    my $prior = $self->_get_prior_for_pse($pidfa);
    if ( not $prior ) {
        $self->error_message('Failed to get prior for PIDFA! '.$prior->id);
        return;
    }

    return ($prior, $pidfa, $qidfgm);
}

sub _get_prior_for_pse {
    my ($self, $pse, $ps_id) = @_;

    Carp::confess('No PSE given to get prior!') if not $pse;

    my ($tp_pse) = GSC::TppPSE->get(pse_id => $pse->id);
    if ( not $tp_pse ) {
        $self->error_message('No transfer pattern pse for PSE! '.$pse->id);
        return;
    }

    if ( not $tp_pse->prior_pse_id ) {
        $self->error_message('No prior pse id in transfer pattern for PSE! '.Data::Dumper::Dumper($tp_pse));
        return;
    }

    my %prior_params = ( pse_id => $tp_pse->prior_pse_id );
    $prior_params{ps_id} = $ps_id if $ps_id;
    my $prior = GSC::PSE->get(%prior_params);
    if ( not $prior ) {
        $self->error_message('Failed to get prior PSE for params! '.Data::Dumper::Dumper(\%prior_params));
        return;
    }

    return $prior;
}

sub _get_future_for_pse {
    my ($self, $pse, $ps_id) = @_;

    my ($tp_pse) = GSC::TppPSE->get(pse_id => $pse->id);
    if ( not $tp_pse ) {
        $self->error_message('No transfer pattern pse for prior! '.$pse->id);
        return;
    }

    if ( not $tp_pse->pse_id ) {
        $self->error_message('No PIDFA pse id in transfer pattern for prior!');
        return;
    }

    my %future_params = ( pse_id => $tp_pse->prior_pse_id );
    $future_params{ps_id} = $ps_id if $ps_id;
    my $future = GSC::PSE->get(%future_params);
    if ( not $future ) {
        $self->error_message('Failed to get future PSE for params! '.Data::Dumper::Dumper(\%future_params));
        return;
    }

    return $future;
}

sub _additional_params_to_fix { return; }

sub _fix_params_for {
    my ($self, $pse, $params_to_fix) = @_;

    for my $param_name ( keys %$params_to_fix ) {
        my $new_param_value = $self->_fix_param($pse, $param_name, $params_to_fix->{$param_name});
        return if not $new_param_value;
    }

    return 1;
}

sub _fix_param {
    my ($self, $pidfa, $param_name, $param_value) = @_;

    die 'No PIDFA' if not $pidfa;
    die 'No param name' if not $param_name;
    die 'No param value' if not $param_value;

    my $param_display_name = join(' ', split('_', $param_name));

    # Remove extra params
    my @params = GSC::PSEParam->get(pse_id => $pidfa->id, param_name => $param_name);
    $self->status_message(ucfirst($param_display_name).' params: '.@params);
    for ( my $i = 1; $i < $#params; $i++ ) {
        $params[$i]->delete;
    }

    # Make sure the one left over is correct, if not delete and create a new one
    $self->status_message("Current $param_display_name: ".( @params ? $params[0]->param_value : 'NULL' ));
    if ( not @params or $params[0]->param_value ne $param_value ) {
        $params[0]->delete if @params;
        my $new_param = GSC::PSEParam->create(
            pse_id => $pidfa->pse_id,
            param_name => $param_name,
            param_value => $param_value,
        );
        my $new_param_value = $new_param->param_value;
        if ( $new_param_value ne $param_value ) {
            $self->error_message("Failed to create new param for $param_name!");
            return;
        }
        $self->status_message("New $param_display_name: $new_param_value");
    }

    return 1;
}

1;

