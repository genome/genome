package Genome::Model::ReferenceAlignment::Command::CaptureQc;
use strict;
use warnings;
use Cwd;

class Genome::Model::ReferenceAlignment::Command::CaptureQc {
    is => 'Genome::Command::Base',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            doc => 'use these models and their samples for QC',
            shell_args_position => 1,
        },
    ],
    has_optional => [
        qc_directory => {
            is => 'Text',
            default => cwd(),
            doc => 'path to gmt capture germline-model-group-qc output',
        },
        output_directory => {
            is => 'Text',
            default => cwd(),
            doc => 'Dir to store sample/index/pool summaries',
        },
        debug => {
            is => 'Boolean',
            default => 0,
            doc => 'Print debug messages that could help fix oddities',
        },
        find_qc_models => {
            is => 'Boolean',
            default => 1,
            doc => 'From models given, find and operate on their qc counterparts',
        },
        prefix => {
            is => 'Text',
            default => '',
            doc => 'Prefix generated files with this text',
        },
    ],
    doc => 'Summarize information on models from germline-model-group and germline-model-group-qc',
};

sub help_brief { 'Summarize information on models from germline-model-group and germline-model-group-qc' }

sub help_detail{ help_brief() . "\nExample: genome model reference-alignment capture-qc 10407 -p 10407" }

sub execute {
    my $self = shift;

    $self->output_metrics( $self->get_metrics($self->get_models) );

    return 1;
}

sub get_models {
    my $self = shift;

    my @ref_align_models = grep { $_->isa('Genome::Model::ReferenceAlignment') } $self->models;
    unless (@ref_align_models) {
        die "None of the provided models were reference alignment!";
    }

    my @models;
    if ($self->find_qc_models) {
        my @missing_qc_models;
        for my $ref_align_model (@ref_align_models) {
            my @instrument_data = $ref_align_model->instrument_data;
            my @qc_models = Genome::Model::ReferenceAlignment->default_lane_qc_model_for_instrument_data(@instrument_data);
            if (@qc_models) {
                push @models, @qc_models;
            } else {
                push @models, $ref_align_model;
                push @missing_qc_models, $ref_align_model->id;
            }
        }
        if (@missing_qc_models) {
            $self->status_message('These ' . scalar @missing_qc_models . ' models are missing qc models: ' . join(',', @missing_qc_models));
        }
    } else {
        @models = @ref_align_models;
    }

    return @models;
}

sub metric_names {
    ( #all are wingspan 0
        '%40Depth',
        '%30Depth',
        '%20Depth',
        '%15Depth',
        '%10Depth',
        '%Dup',    # percent duplication
        '%Mapped',
        '%UniqOn', # % of Unique On Target Reads
        '%UniqOff',# % of Unique Off Target Reads
        '%TotalOn',
        '%TotalOff',
        '%Unaligned',
        'TotalReads',
        'SNPsCalled',
        'WithGenotype',
        'MetMinDepth',
        'Reference',
        'RefMatch',
        'RefWasHet',
        'RefWasHom',
        'Variant',
        'VarMatch',
        'HomWasHet',
        'HetWasHom',
        'VarMismatch',
        'VarConcordance',
        'RareHomConcordance',
        'OverallConcordance',
    )
}

sub write_averaged_summary {
    #Accepts builds grouped by a value, and averages all build metrics in each grouping
    my $self = shift;
    my $fh = shift || die;
    my $grouping_value_to_builds = shift || die;
    my $build_to_metrics = shift || die;
    my $grouping_metric = shift || die;

    #Write column headers
    print $fh join ("\t", (
            $grouping_metric,
            $self->metric_names,
        )) . "\n";

    while (my ($grouping_value,$builds) = each %$grouping_value_to_builds) {
        my %sum_value;
        for my $build (@$builds) {
            for my $metric_name ($self->metric_names) {
                $sum_value{$metric_name} += $build_to_metrics->{$build->id}{$metric_name} || 0;
            }
        }
        print $fh join("\t",
            $grouping_value,
            map{ sprintf ("%.2f", $sum_value{$_} / @$builds); } $self->metric_names
        ) . "\n";
    }
}

sub write_full_summary {
    my $self = shift;
    my $fh = shift || die;
    my $build_to_metrics = shift || die;

    #Write column headers
    print $fh join ("\t", (
            'Model',
            'Build',
            'ROI',
            'Sample',
            'Lane',
            'Libraries',
            'Index',
            'Pooled library',
            $self->metric_names,
        )) . "\n";

    for my $build_id (keys %{$build_to_metrics}){
        my $build = Genome::Model::Build->get($build_id);
        my $model = $build->model;

        #Take the first instrument_data's index until a decision is made
        #  on how to handle per-sample data, when per-instrument-data data is unavailable
        my ($index) = map{$_->index_sequence}$model->instrument_data;
        $index = '-' unless defined $index;
        my $pool = $self->_resolve_pooled_sample_name_for_instrument_data($model->instrument_data);
        $pool = '-' unless defined $pool;
        my $libraries = join ' ', map{$_->library->name}$build->instrument_data;

        my $lane = join ' ', map{$_->lane}grep{defined $_->lane}$build->instrument_data;

        print $fh join("\t",
            $model->id,
            $build->id,
            ($build->region_of_interest_set_name || '-'),
            $model->subject->name,
            $lane,
            $libraries,
            $index,
            $pool,
            map{ sprintf ("%.2f", $build_to_metrics->{$build->id}{$_}) } $self->metric_names
        ) . "\n";
    }
}

sub _resolve_pooled_sample_name_for_instrument_data {
    my $self = shift;
    my $instrument_data = shift;

    return unless $instrument_data->can('index_sequence');
    my $index = $instrument_data->index_sequence;
    if($index) {
        my $instrument_data_class = $instrument_data->class;
        my $pooled_subset_name = $instrument_data->subset_name;
        $pooled_subset_name =~ s/${index}$/unknown/;

        my $pooled_instrument_data = $instrument_data_class->get(
                run_name => $instrument_data->run_name,
                subset_name => $pooled_subset_name,
                index_sequence => 'unknown',
                );
        return unless $pooled_instrument_data;

        my $sample = $pooled_instrument_data->sample;
        return unless $sample;

        return $sample->name;
    }

    return;
}

sub output_metrics {
    my $self = shift;
    my ($build_to_metrics, $index_to_builds, $pool_to_builds) = @_;

    print "Summarizing data\n" if $self->debug;
    my $dir = $self->output_directory;
    print "Writing subject_summary\n" if $self->debug;
    my $prefix = $self->prefix;
    $prefix .= '_' if $prefix;
    $self->write_full_summary(
        Genome::Sys->open_file_for_overwriting($dir . "/".$prefix."subject.tsv"),
        $build_to_metrics,
    );
    print "Writing index_summary\n" if $self->debug;
    $self->write_averaged_summary(
        Genome::Sys->open_file_for_overwriting($dir . "/".$prefix."index.tsv"),
        $index_to_builds,
        $build_to_metrics,
        'Index',
    );
    print "Writing pool_summary\n" if $self->debug;
    $self->write_averaged_summary(
        Genome::Sys->open_file_for_overwriting($dir . "/".$prefix."pool.tsv"),
        $pool_to_builds,
        $build_to_metrics,
        'Pool',
    );
}

sub get_metrics {
    my $self = shift;
    my @models = @_;
    my ($build_to_metrics, $index_to_builds, $pool_to_builds);

    for my $model (@models) {
        my $build = $model->last_succeeded_build;
        unless ($build) {
            $self->warning_message("No succeeded build found for model " . $model->__display_name__);
            next;
        }
        print 'Gathering data for build: ' . $build->id . "\n" if $self->debug;

        #Non qc data
        my ($align_stats, $cov_stats);
        eval{
            my $cov_result = $build->coverage_stats_result;
            $align_stats = $cov_result->alignment_summary_hash_ref->{0};
            $cov_stats = $cov_result->coverage_stats_summary_hash_ref->{0};
        };
        if( $@ ){
            eval{
                $align_stats = $build->alignment_summary_hash_ref->{0};
                $cov_stats = $build->coverage_stats_summary_hash_ref->{0};
            };
        }
        $build_to_metrics->{$build->id} = $self->get_metrics_for_non_qc_data(
            $align_stats,
            $cov_stats,
            $build->merged_alignment_result->merged_alignment_bam_flagstat
        );

        #qc data
        $build_to_metrics->{$build->id} = {%{$self->get_metrics_for_qc_data($build)},%{$build_to_metrics->{$build->id}}};


        #Find all the models for each index and pool
        #Take the first instrument_data's index until a decision is made
        #  on how to handle per-sample data, when per-instrument-data data is unavailable
        #  although this won't matter once we start using QC models only
        my $index = eval{(map{$_->index_sequence}$model->instrument_data)[0]};
        $index = '-' unless defined $index;
        my $pool = eval{$self->_resolve_pooled_sample_name_for_instrument_data($model->instrument_data)};
        $pool = '-' unless defined $pool;

        #Build reference of each index and pool to all of their builds to summarize data per each of these
        push @{$index_to_builds->{$index}}, $build;
        push @{$pool_to_builds->{$pool}}, $build;
    }

    return ($build_to_metrics, $index_to_builds, $pool_to_builds);
}

sub get_metrics_for_non_qc_data {
    my $self = shift;
    my $align_stats = shift;
    my $cov_stats = shift;
    my $flagstat = shift;

    my %metric_to_value;

    my $unique_on_target = $align_stats->{unique_target_aligned_bp} || 0;
    my $duplicate_on_target = $align_stats->{duplicate_target_aligned_bp} || 0;
    my $unique_off_target = $align_stats->{unique_off_target_aligned_bp} || 0;
    my $duplicate_off_target = $align_stats->{duplicate_off_target_aligned_bp} || 0;
    my $unaligned = $align_stats->{total_unaligned_bp} || 0;

    my $total = $unique_on_target+$duplicate_on_target+$unique_off_target+$duplicate_off_target+$unaligned;
    my ($percent_unique_on_target, $percent_unique_off_target, $percent_duplicate_off_target, $percent_duplicate_on_target, $percent_unaligned) = (0,0,0,0,0);
    eval{
        $percent_unique_on_target = $unique_on_target/$total*100;
        $percent_duplicate_on_target = $duplicate_on_target/$total*100;
        $percent_unique_off_target = $unique_off_target/$total*100;
        $percent_duplicate_off_target = $duplicate_off_target/$total*100;
        $percent_unaligned = $unaligned/$total*100;
    };

    #Default coverage depth percentage is zero - there may not be depth information, especially for 30x and 40x
    $metric_to_value{'%40Depth'} = $cov_stats->{40}{pc_target_space_covered} || 0;
    $metric_to_value{'%30Depth'} = $cov_stats->{30}{pc_target_space_covered} || 0;
    $metric_to_value{'%20Depth'} = $cov_stats->{20}{pc_target_space_covered} || 0;
    $metric_to_value{'%15Depth'} = $cov_stats->{15}{pc_target_space_covered} || 0;
    $metric_to_value{'%10Depth'} = $cov_stats->{10}{pc_target_space_covered} || 0;

    $metric_to_value{'%Dup'} = $align_stats->{percent_duplicates} || 0;

    $metric_to_value{'%Mapped'} = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat)->{'reads_mapped_percentage'} || 0;

    $metric_to_value{'%UniqOn'} = $percent_unique_on_target || 0;
    $metric_to_value{'%UniqOff'} = $percent_unique_off_target || 0;

    $metric_to_value{'%TotalOn'} = $percent_duplicate_on_target + $percent_unique_on_target || 0;
    $metric_to_value{'%TotalOff'} = $percent_duplicate_off_target + $percent_unique_off_target || 0;
    $metric_to_value{'%Unaligned'} = $percent_unaligned || 0;
    $metric_to_value{'TotalReads'} = $total || 0;


    if ($self->debug) {
        while (my ($metric,$value) = each %metric_to_value) {
            print "Metric for $metric : " . $value . "\n";
        }
    }

    return \%metric_to_value;
}

sub get_metrics_for_qc_data {
    my $self = shift;
    my $build = shift;
    my $qc_file;
    my ($compare_snp_result) = Genome::Model::Tools::Analysis::LaneQc::CompareSnpsResult->get('users.user_id' => $build->id);

    if($compare_snp_result){
        $qc_file = $compare_snp_result->output_file;
    } else {
        my $qc_dir = $self->qc_directory or return _empty_qc_data();
        my $dir_name = $build->subject_name;
        ($qc_file) = `ls $qc_dir/$dir_name/*.qc 2>/dev/null`;
        return _empty_qc_data() unless $qc_file;
        chomp $qc_file;
    }

    print "QC File for build ( " . $build->id . " ) : $qc_file\n" if $self->debug;

    my $fh = Genome::Sys->open_file_for_reading($qc_file);

    my @values;
    my $line_number = 1;
    for my $line (<$fh>){
        if(2 == $line_number++){
            $line =~ s/%//g; #Percent symbols break converting this to a number
            @values = split /\s+/, $line;
            last;
        }
    }

    my %metric_to_value;
    $metric_to_value{SNPsCalled } = $values[1] || 0;
    $metric_to_value{WithGenotype } = $values[2] || 0;
    $metric_to_value{MetMinDepth } = $values[3] || 0;
    $metric_to_value{Reference } = $values[4] || 0;
    $metric_to_value{RefMatch } = $values[5] || 0;
    $metric_to_value{RefWasHet } = $values[6] || 0;
    $metric_to_value{RefWasHom } = $values[7] || 0;
    $metric_to_value{Variant } = $values[8] || 0;
    $metric_to_value{VarMatch } = $values[9] || 0;
    $metric_to_value{HomWasHet } = $values[10] || 0;
    $metric_to_value{HetWasHom } = $values[11] || 0;
    $metric_to_value{VarMismatch } = $values[12] || 0;
    $metric_to_value{VarConcordance } = $values[13] || 0;
    $metric_to_value{RareHomConcordance } = $values[14] || 0;
    $metric_to_value{OverallConcordance } = $values[15] || 0;

    if ($self->debug) {
        while (my ($metric,$value) = each %metric_to_value) {
            print "Metric for $metric : " . $value . "\n";
        }
    }

    return \%metric_to_value;
}

sub _empty_qc_data {
    return {
        SNPsCalled => 0,
        WithGenotype => 0,
        MetMinDepth => 0,
        Reference => 0,
        RefMatch => 0,
        RefWasHet => 0,
        RefWasHom => 0,
        Variant => 0,
        VarMatch => 0,
        HomWasHet => 0,
        HetWasHom => 0,
        VarMismatch => 0,
        VarConcordance => 0,
        RareHomConcordance => 0,
        OverallConcordance => 0,
    };
}
