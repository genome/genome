package Genome::Model::MetagenomicComposition16s::Command::RunStatus; 

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require Carp;

class Genome::Model::MetagenomicComposition16s::Command::RunStatus {
    is => 'Command::V2',
    has => [
        run_name => {
            is_many => 0,
            shell_args_position => 1,
            doc => 'Run name.',
        },
        region_number => {
            shell_args_position => 2,
            doc => 'Region Number.',
        },
    ],
    has_optional => [
        processing_profile => {
            is => 'Genome::ProcessingProfile::MetagenomicComposition16s',
            doc => 'Limit models to this processing profile. Current default profiles are 2571784 (rdp 2.2) and 2752939 (rdp2.5). But any MC16s processing profile can be specified.'
        },
        separator => {
            is => 'Text',
            default_value => ',',
            doc => 'Separator character to use when printing, saving and email attachments. For tab separator, use "tab".',
        },
    ],
};

sub help_brief { 
    return 'Get the status for models for a 454 run';
}

sub help_detail {
    return help_brief();
}

sub execute {
    my $self = shift;
    
    my @instrument_data = Genome::InstrumentData::454->get(
        run_name => $self->run_name,
        region_number => $self->region_number,
    );

    if ( not @instrument_data ) {
        $self->error_message('Cannot find 454 instrument_data for run name ('.$self->run_name.') and region ('.$self->region_number.')');
        return;
    }

    my $processing_profile_id = $self->processing_profile->id;
    my @subjects_to_skip = map { qr/$_/ } (qw/ ^nctrl$ ^n\-cntrl$ /);
    my @rows;
    for my $instrument_data ( @instrument_data ) {
        next if $instrument_data->ignored;
        my $library = $instrument_data->library;
        my $sample = $library->sample;
        next if grep { $sample->name =~ $_ } @subjects_to_skip;
        my @row;
        push @rows, \@row;
        push @row, $sample->name;
        push @row, $instrument_data->id;
        my ($model) = sort { $b->id <=> $a->id } grep { $_->subject_id eq $instrument_data->sample_id } grep { $_->processing_profile_id eq $processing_profile_id } map { $_->model } Genome::Model::Input->get(name => 'instrument_data', value_id => $instrument_data->id,);
        if ( not $model ) {
            push @row, '', '', '', ''. '', '';
            next;
        }
        my $build = $model->current_build;
        if ( not $build ) {
            push @row, $model->id, 'NA', ( $model->build_requested ? 'Requested' : 'None' ), '', '';
            next;
        }
        push @row, map { $build->$_ } (qw/ model_id id status /);
        if ( $build->status eq 'Succeeded' ) {
            push @row, (($build->amplicons_processed_success * 100).'%');
        }
        else {
            push @row, 'NA';
        }
        my $files_string = join ( ' ', grep { -s $_ } $build->oriented_fasta_files );
        push @row, $files_string;
    }

    my $sep = ( $self->separator =~ /^tab$/i ? "\t" : $self->separator );
    print join( "\n", map { join($sep, @$_) } [qw/ sample-name instrument-data-id model-id build-id status process-success oriented-fastas /], @rows)."\n"; 

    return 1;
}

1;

