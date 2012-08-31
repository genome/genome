package Genome::InstrumentData::AlignmentResult::Merged::CufflinksExpression;

use strict;
use warnings;

use Genome;
use Sys::Hostname;
use File::Path;

class Genome::InstrumentData::AlignmentResult::Merged::CufflinksExpression {
    is => ['Genome::SoftwareResult::Stageable'],
    has_input => [
        alignment_result_id => {
            is => 'Number',
            doc => 'ID of the result for the alignment data upon which to run coverage stats',
        },
    ],
    has_param => [
        cufflinks_params => {
            is => 'Text',
            doc => 'A string of extra params to pass to cufflinks.',
        },
        cufflinks_version => {
            is => 'Text',
            doc => 'The version of cufflinks to use.',
        },
        mask_reference_transcripts => {
            is => 'Text',
            doc => 'The name of an annotation file of known transcripts to mask during expression estimation.',
        },
        annotation_reference_transcripts_mode => {
            is => 'Text',
            doc => 'The mode for cufflinks to utilize the known annotation files.',
            valid_values => ['reference only', 'reference guided', 'de novo'],
        },
    ],
    has_metric => [
        _log_directory => {
            is => 'Text',
            doc => 'Path where workflow logs were written',
        },
        #many other metrics exist--see sub _generate_metrics
    ],
    has => [
        alignment_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged',
            id_by => 'alignment_result_id',
            doc => 'the alignment data upon which to run coverage stats',
        },
    ],
    has_transient_optional => [
        log_directory => {
            is => 'Text',
            doc => 'Path to write logs from running the workflow',
        },
    ],
};

sub resolve_allocation_subdirectory {
    my $self = shift;

    my $hostname = hostname;
    my $user = $ENV{'USER'};
    my $base_dir = sprintf("cufflinks_expression-%s-%s-%s-%s", $hostname, $user, $$, $self->id);

    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'build_merged_alignments','cufflinks_expression',$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    return 'info_genome_models';
}

sub _staging_disk_usage {
    #need the allocation created in advance for this process
    return 5_000_000; #TODO better estimate
}

sub _working_dir_prefix {
    return 'cufflinks-expression';
}

sub _gather_params_for_get_or_create {
    my $class = shift;

    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key};
        }
    }

    my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_input);
    my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_param);

    my %software_result_params = (
        params_id=>$params_bx->id,
        inputs_id=>$inputs_bx->id,
        subclass_name=>$class
    );

    return {
        software_result_params => \%software_result_params,
        subclass => $class,
        inputs=>\%is_input,
        params=>\%is_param,
    };
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    $self->_prepare_staging_directory;

    # Alignment inputs
    my $alignment_result = $self->alignment_result;

    my $log_dir = $self->log_directory;
    unless($log_dir) {
        $log_dir = '' . $self->temp_staging_directory;
    }
    $self->_log_directory($log_dir);
    my $expression_directory = $self->temp_staging_directory;
    unless (Genome::InstrumentData::AlignmentResult::Command::CufflinksExpression->execute(
        alignment_result => $alignment_result,
        reference_build => $alignment_result->reference_build,
        annotation_build => $alignment_result->annotation_build,
        expression_directory => $expression_directory,
        cufflinks_params => $self->cufflinks_params,
        cufflinks_version => $self->cufflinks_version,
        mask_reference_transcripts =>$self->mask_reference_transcripts,
        annotation_reference_transcripts_mode => $self->annotation_reference_transcripts_mode,
    )) {
        return;
    }

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    $self->_generate_metrics();

    return $self;
}

sub _generate_metrics {
    my $self = shift;
    $self->status_message('Currently there are no metrics recorded for '. __PACKAGE__);
    #my $coverage_dir = $self->output_dir;
    #my @metric_files = glob($coverage_dir . "/*_STATS.txt");
    #for my $metric_file (@metric_files) {
    #    my $file_basename = File::Basename::basename($metric_file);
    #    my ($stat_type) = $file_basename =~ m/(.*)_STATS/;
    #    my $metric_reader = Genome::Utility::IO::SeparatedValueReader->create(
    #        input => $metric_file,
    #        separator => "\t",
    #    );
    #    while (my $data = $metric_reader->next) {
    #        for my $key (keys %{$data}) {
    #            my $stat_key = sprintf("%s %s", $stat_type, $key);
    #            $self->add_metric(metric_name => $stat_key, metric_value => $data->{$key});
    #        }
    #    }
    #}
    return 1;
}


1;
