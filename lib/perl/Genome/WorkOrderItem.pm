package Genome::WorkOrderItem;

use strict;
use warnings;

use Genome;

use Carp 'confess';
use Data::Dumper 'Dumper';

class Genome::WorkOrderItem {
    id_by => [
        woi_id => {
            is => 'Integer',
            len => 10,
            column_name => 'WOI_ID',
        },
    ],    
    has => [
        setup_wo_id => {
            is => 'Integer',
            len => 10,
        },
        status => {
            is => 'Text',
            len => 32,
        },
        work_order => {
            is => 'Genome::WorkOrder',
            id_by => 'setup_wo_id',
            doc => 'The work order for this item.',
        },
    ],
    has_optional => [
        barcode => {
            is => 'Text',
            len => 16,
        },
        dna_id => {
            is => 'Text',
            len => 20,
        },
        sample => {
            is => 'Genome::Sample',
            id_by => 'dna_id',
        },
        parent_woi_id => {
            is => 'Integer',
            len => 10,
        },
        pipeline_id => {
            is => 'Integer',
            len => 10,
        },
    ],
    has_many_optional => [
        sequence_products => {
            is => 'Genome::Site::TGI::WorkOrderItemSequenceProduct',
            reverse_as => 'work_order_item',
        },
        sequence_items => {
            is => 'Genome::Site::TGI::SequenceItem',
            via => 'sequence_products',
            to => 'sequence_item',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            via => 'sequence_products',
            to => 'instrument_data',
        },
    ],
    data_source => 'Genome::DataSource::Oltp',
    table_name => 'work_order_item',
};

sub models {
    my $self = shift;

    # Strategy
    # 1 - via instrument data inputs.
    #   a - Only works for solexa and 454
    #
    my @instrument_data = $self->instrument_data();
    if ( @instrument_data ) {
        my @instrument_data_inputs = Genome::Model::Input->get(
            name => 'instrument_data',
            value_id => [ map { $_->id } @instrument_data ],
        );
        return unless @instrument_data_inputs;
        my %model_ids = map { $_->model_id => 1 } @instrument_data_inputs;
        return unless %model_ids;
        return Genome::Model->get(genome_model_id => [ keys %model_ids ]);
    } else {
        $self->warning_message("No sequence products (instrument data) found for work order item (".$self->id.").  Attempting to get models via dna_id.");
    }

    # 2 - Round about, not so accurate way, but works sometimes.
    #      WOI have dna id or barcode.  We only get for barcode until someone
    #      wants to write the logic.  But once the sanger reads are back filled
    #      into the woi seq product table, this should be removed. Ther should be
    #      a dna_id OR barcode
    #
    #  a - get dna_id.  if no dna id, then this doesn't work.  die.
    #  b - get models w/ subject_id == dna_id
    #  
    
    if ( $self->dna_id ) {
        return Genome::Model->get(subject_id => $self->dna_id);
    }

    if ( $self->barcode ) {
        confess "Can't get models for a work order item (".$self->id.") that only has a barcode, and no sequence products or dna_id.";
    }
}

sub _summarize {

    my ($process_to) = @_;
    my @initials;

    if (length($process_to) > 25) {
        my @words = split(/\s+/,$process_to);
        push @initials, uc(substr($_, 0, 1)) for @words; 
        return join('', @initials);
    }

    return $process_to;
}

sub _hash_up {

    my ($events, $pse, $sort_order) = @_;

    my $status     = $pse->pse_status;
    my $process_to = _summarize( $pse->process_to );

    if (not exists $events->{'production'}->{$process_to}->{$status}->{'count'}) {
        $events->{'production'}->{$process_to}->{$status}->{'count'} = 0;
    }
    return {
        count      => $events->{'production'}->{$process_to}->{$status}->{'count'} + 1,
        sort_order => $sort_order,
    };
}

sub event_statuses {

    my ($self) = @_;

    my $events = {};
    for my $sp ($self->sequence_items()) {

        my $status = {};

        # get event statuses from LIMS
        # FIXME this will not work b/c the seq items (Genome::Site::TGI::SequenceItem) here do not have the code to get the creation pse
        #  if adding this functionality, DO NOT add GSC classes, use a view in Genome::Site::TGI.
        if ($sp->sequence_item_type eq 'index illumina') {

            my $creation_pse = $sp->get_creation_event();
            my $setup_pse = $creation_pse->get_first_prior_pse_with_process_to('set up flow cell');
            my $pidfa_pse = $creation_pse->get_first_active_subsequent_pse_with_process_to('prepare instrument data for analysis');
            my $queue_pse = $creation_pse->get_first_active_subsequent_pse_with_process_to('queue instrument data for genome modeling');

            # put the pses in 3 status buckets- "pending", "completed", "failed"
            for my $pse ($creation_pse, $setup_pse, $pidfa_pse, $queue_pse) {

                my $pse_status = $pse->pse_status();
                my @completed = qw(COMPLETED);
                my @failed = qw(FAILED);

                # default status is "pending"
#                $status->{$pse->pse_id} = $pse_status;
                $status->{$pse->pse_id} = 'pending';

                if (grep /^$pse_status$/i, @completed) {
                    $status->{$pse->pse_id} = 'completed';
                } elsif (grep /^$pse_status$/i, @failed) {
                    $status->{$pse->pse_id} = 'failed';
                }
            }

            # generate lane summary
            $events->{'production'}->{_summarize($creation_pse->process_to)}->{$status->{$creation_pse->pse_id}} = _hash_up($events, $creation_pse, 2);

            # set up flow cell
            $events->{'production'}->{_summarize($setup_pse->process_to)}->{$status->{$setup_pse->pse_id}} = _hash_up($events, $setup_pse, 1);
            # PIDFA
            $events->{'production'}->{_summarize($pidfa_pse->process_to)}->{$status->{$pidfa_pse->pse_id}} = _hash_up($events, $pidfa_pse, 3);

            # queue inst data
            $events->{'production'}->{_summarize($queue_pse->process_to)}->{$status->{$queue_pse->pse_id}} = _hash_up($events, $queue_pse, 4);
        }

        # get statuses of latest builds from canonical models
        my $sample = $self->sample();
        my @models = sort { $a->id <=> $b->id } Genome::Model->get(subject_id => $sample->id);
        my $model = $models[0];

        if ($model) {
            my $build = $model->latest_build();
            next if not $build;
            my $build_status = $build->master_event_status();

            # same status buckets as pses- pending, completed, failed
            if ($build_status =~ /SUCCEEDED/i) {
                $build_status = 'completed';
            } elsif ($build_status =~ /FAILED/i) {
                $build_status = 'failed';
            } else {
                # default status is same as PSEs- "pending"
                $build_status = 'pending';
            }

            $events->{'analysis'}->{'builds'}->{$build_status}->{'count'}++;
            $events->{'analysis'}->{'builds'}->{$build_status}->{'sort_order'} = 5;
        }
    }

    return $events;
}





# prepare dna beads
# setup flow cell
   

# 454: demux 454 regions
# illumina: copy sequence files
# 3730: analyze traces





#$HeadURL$
#$Id$
