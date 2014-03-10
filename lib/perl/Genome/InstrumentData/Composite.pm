package Genome::InstrumentData::Composite;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Composite {
    is => 'UR::Object',
    has => [
        inputs => {
            is => 'HASH',
            doc => 'a mapping from keys in the strategy to their values (the source data and reference sequences to use)',
        },
        strategy => {
            is => 'Text',
            doc => 'The instructions of how the inputs are to be aligned and/or filtered',
        },
        force_fragment => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Treat all reads as fragment reads',
        },
        merge_group => {
            is => 'Text',
            default_value => 'sample',
            valid_values => ['sample', 'all'],
            doc => 'When merging, collect instrument data together that share this property',
        },
        _merged_results => {
            is => 'Genome::InstrumentData::AlignedBamResult',
            is_transient => 1,
            is_optional => 1,
            doc => 'Holds the underlying merged results',
            is_many => 1,
        },
    ],
    has_transient_optional => {
        log_directory => {
            is => 'Text',
            doc => 'where to write the workflow logs',
        },
    },
};

#This method should just use the one from Genome::SoftwareResult and then get will return the existing result and create will run the alignment dispatcher
sub get_or_create {
    my $class = shift;

    #my $self = $class->SUPER::create(@_);
    my $self; # = $class->SUPER::create(@_);
    if (!ref($_[0])) {
        # not passed a boolexpr
        my %params = @_;
        my $inputs_from_args;
        $inputs_from_args = delete $params{'inputs'};
        $self = $class->SUPER::create(%params);
        $self->{'inputs'} = $inputs_from_args if $inputs_from_args;
    } else {
        $self = $class->SUPER::create(@_);
    }

    $self->debug_message('Get or create composite instrument data...');

    my $inputs = $self->inputs;
    $inputs->{force_fragment} = $self->force_fragment;
    my $strategy = $self->strategy;

    $self->debug_message('Create composite workflow...');
    my $generator = Genome::InstrumentData::Composite::Workflow->create(
#        inputs => {
#            #%{ $self->inputs },
#            (map { $_ => $inputs->{$_}} keys %$inputs),
#            force_fragment => $self->force_fragment,
#        },
        strategy => $self->strategy,
        merge_group => $self->merge_group,
        log_directory => $self->log_directory,
    );
    $self->debug_message('Create composite workflow...OK');
    $generator->{inputs} = $inputs;

    $self->debug_message('Execute composite workflow...');
    $generator->dump_status_messages(1);
    unless($generator->execute) {
        die $self->error_message('Failed to execute workflow.');
    }
    $self->debug_message('Execute composite workflow...OK');

    $self->debug_message('Get software results...');
    my @result_ids = $generator->_result_ids;
    my @all_results = Genome::SoftwareResult->get(\@result_ids);
    $self->debug_message('Found '.@all_results.' software results');

    #TODO If this is made a result, too, register as a user
    #for my $result (@all_results) {
    #    $result->add_user(label => 'uses', user => $self);
    #}

    my @merged_results = grep($_->isa('Genome::InstrumentData::AlignedBamResult'), @all_results);
    $self->debug_message('Found '.@merged_results.' merged results');
    $self->_merged_results(\@merged_results);

    return $self;
}

sub bam_paths {
    my $self = shift;

    my @results = $self->_merged_results;

    my @bams;
    for my $result (@results) {
        my $bam = $result->bam_file;
        push @bams, $bam;
    }

    return @bams;
}

1;
