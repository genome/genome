package Genome::InstrumentData::Command::AlignReads;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::AlignReads {
    is => ['Command::V2'],
    has_input => [
        name => {
            is => 'Text',
            doc => 'The aligner to use',
        },
        version => {
            is => 'Text',
            doc => 'Version of the aligner to use',
        },
        params => {
            is => 'Text',
            doc => 'Parameters to pass to the aligner',
        },
        instrument_data_id => {
            is => 'Number',
            doc => 'Id of the data to align',
        },
        reference_build_id => {
            is => 'Number',
            doc => 'Id of the reference to which to align',
        },
        picard_version => {
            is => 'Text',
            doc => 'The version of Picard to use when needed by aligners/filters',
        },
        samtools_version => {
            is => 'Text',
            doc => 'The version of Samtools to use when needed by aligners/filters',
        },
    ],
    has_optional_input => [
        annotation_build_id => {
            is => 'Number',
            doc => 'Id of the annotation build to use when aligning',
            is_optional => 1,
        },
        instrument_data_segment_id => {
            is => 'Text',
            doc => 'A specific segment of the data to align',
        },
        instrument_data_segment_type => {
            is => 'Text',
            doc => 'How the data is segmented',
        },
        instrument_data_filter => {
            is => 'Text',
            valid_values => [ 'forward-only', 'reverse-only', undef ],
        },
        force_fragment => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Treat all reads as fragment reads',
        },
        trimmer_name => {
          is => 'Text',
          doc => 'name of the read trimmer to use before alignment is performed'
        },
        trimmer_params => {
          is => 'Text',
          doc => 'params to use with the read trimmmer'
        },
        trimmer_version => {
          is => 'Text',
          doc => 'version of the read trimmer to use'
        },
    ],
    has_transient_optional => [
        _instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'The data to align',
            id_by => 'instrument_data_id',
        },
        _reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference to which to align',
            id_by => 'reference_build_id',
        },
    ],
    has_optional_output => [
        result_id => {
            is => 'Text',
            doc => 'The result generated/found when running the command',
        },
    ],
    has_param => [
        lsf_resource => {
            calculate => q{
                $self->bsub_rusage;
            },
            default_value => &_fallback_lsf_resource, # workflow doesn't support varying this per instance
        },
    ]
};

sub results_class {
    my $self = shift;

    my $aligner_name = $self->name;
    return 'Genome::InstrumentData::AlignmentResult::' .
        Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($aligner_name);
}

sub _fallback_lsf_resource {
    my $tmp_mb = 90000;
    my $mem_mb = 1024 * 14; # increased b/c we have about 16 GB available when 6 jobs run on a 96 Gb server
    my $cpus = 4;

    my $mem_kb = $mem_mb*1024;
    my $tmp_gb = $tmp_mb/1024;

    my $user = getpwuid($<);
    my $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_DEFAULT};

    # TODO: in-house TGI concepts like alignment_pd shouldn't be methods on the generic config module :( -ssmith
    $queue = 'alignment-pd' if (Genome::Config->can('should_use_alignment_pd') and Genome::Config->should_use_alignment_pd);

    #my $host_groups;
    #my $command = qq(bqueues -l $queue | grep ^HOSTS:);
    #$host_groups = qx($command);
    #$host_groups =~ s/\/\s+/\ /;
    #$host_groups =~ s/^HOSTS:\s+//;

    my $select  = "select[ncpus >= $cpus && mem >= $mem_mb && gtmp >= $tmp_gb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb, gtmp=$tmp_gb]";
    my $options = "-M $mem_kb -n $cpus -q $queue";

    my $required_usage = "-R '$select $rusage' $options";

    #check to see if our resource requests are feasible (This uses "maxmem" to check theoretical availability)
    #factor of four is based on current six jobs per host policy this should be revisited later
    #my $select_check = "select[ncpus >= $cpus && maxmem >= " . ($mem_mb * 4) . " && gtmp >= $tmp_gb] span[hosts=1]";

    #$command = qq(bhosts -R '$select_check' $host_groups | grep ^blade);
    #$command = qq(bhosts -R '$select_check' $host_groups );
    #my @selected_blades = qx($command);

    #if (@selected_blades) {
        return $required_usage;
    #} else {
    #    die __PACKAGE__->error_message("Failed to find hosts that meet resource requirements ($required_usage). [Looked with $select_check]");
    #}
}

sub bsub_rusage {
    my $self = shift;
    my $delegate = $self->results_class;
    my $rusage = $delegate->required_rusage(instrument_data => $self->_instrument_data);
    return $rusage;
}

sub params_for_alignment {
    my $self = shift;

    my %params = (
        instrument_data_id => $self->instrument_data_id || undef,
        instrument_data_segment_type => $self->instrument_data_segment_type || undef,
        instrument_data_segment_id => $self->instrument_data_segment_id || undef,

        aligner_name => $self->name || undef,
        aligner_version => $self->version || undef,
        aligner_params => $self->params || undef,

        reference_build_id => $self->reference_build_id || undef,
        annotation_build_id => $self->annotation_build_id || undef,

        force_fragment => $self->force_fragment || undef,
        picard_version => $self->picard_version || undef,
        samtools_version => $self->samtools_version || undef,

        trimmer_name => $self->trimmer_name || undef,
        trimmer_version => $self->trimmer_version || undef,
        trimmer_params => $self->trimmer_params || undef,

        filter_name => $self->instrument_data_filter || undef,

        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );

    return \%params;
}

sub execute {
    my $self = shift;

    return $self->_process_alignments('get_or_create');
}

sub shortcut {
    my $self = shift;

    return $self->_process_alignments('get_with_lock');
}

sub _process_alignments {
    my $self = shift;
    my $mode = shift;

    my $instrument_data = $self->_instrument_data;
    $self->status_message("Finding or generating alignments for " . $instrument_data->__display_name__);

    my $alignment = $self->_fetch_alignment_sets($mode);

    if ($mode eq 'get_or_create') {
        unless ($alignment) {
            $self->error_message("Error finding or generating alignments!");
            return 0;
        }
    } elsif ($mode eq 'get_with_lock') {
        unless ($alignment) {
            return undef;
        }
    } else {
        $self->error_message("process/link alignments mode unknown: $mode");
        die $self->error_message;
    }

    $self->status_message("Verifying...");
    unless ($self->verify_successful_completion($alignment)) {
        $self->error_message("Error verifying completion!");
        return 0;
    }

    $self->_link_alignment_to_inputs($alignment);
    $self->result_id($alignment->id);

    $self->status_message("Complete!");
    return 1;
}

sub _fetch_alignment_sets {
    my $self = shift;
    my $mode = shift;

    my $params = $self->params_for_alignment();
    unless ($params) {
        $self->error_message('Could not get alignment parameters for this instrument data');
        return;
    }
    my $alignment = eval { Genome::InstrumentData::AlignmentResult->$mode(%$params) };
    if($@) {
        $self->error_message($mode . ': ' . $@);
        return;
    }

    return $alignment;
}

sub verify_successful_completion {
    my $self = shift;
    my $alignment = shift;

    unless ($alignment->verify_alignment_data) {
        $self->error_message('Failed to verify alignment data: '.  join("\n",$alignment->error_message));
        return 0;
    }

    return 1;
}

sub _link_alignment_to_inputs { #FIXME move inside alignment result generation
    my $self = shift;
    my $alignment = shift;

    my @results = Genome::SoftwareResult->get([$self->instrument_data_id]);
    for my $result (@results) {
        $result->add_user(user => $alignment, label => 'uses');
    }

    return 1;
}

1;
