package Genome::Model::Tools::Capture::GermlineModelGroup2;
use strict;
use warnings;
use Cwd;

class Genome::Model::Tools::Capture::GermlineModelGroup2 {
    is => 'Genome::Command::Base',
    has_input => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            doc => 'use these models and their samples for QC',
        },
        qc_directory => {
            is => 'Text',
            is_optional => 1,
            default => cwd(),
            doc => 'path to gmt capture germline-model-group-qc output',
        },
        output_directory => {
            is => 'Text',
            is_optional => 1,
            default => cwd(),
            doc => 'Dir to store sample/index/pool summaries',
        },
        overwrite_cached_metrics => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Replace all current model metrics associated with this script by generating new ones now.',
        },
    ],
    doc => 'Summarize information on a model group from germline-model-group and germline-model-group-qc',
};

sub help_brief {
    'Summarize information on a model group from germline-model-group and germline-model-group-qc'
}

sub help_detail{
    help_brief() . "\nExample: gmt capture somatic-model-group2 10407 -qc_dir ./germline-model-group-qc_output"
}

sub help_synopsis {
    help_detail()
}

sub execute {
    my $self = shift;

    my $subject_summary_file = Genome::Sys->open_file_for_overwriting($self->output_directory . "/subject_summary.csv");
    my $index_summary_file = Genome::Sys->open_file_for_overwriting($self->output_directory . "/index_summary.csv");
    my $pool_summary_file = Genome::Sys->open_file_for_overwriting($self->output_directory . "/pool_summary.csv");

    if ($self->overwrite_cached_metrics){
        print "Deleting old metrics\n";
        $self->delete_model_metrics();
    }
    unless ($self->non_qc_metrics_exist){
        print "Generating non-qc model metrics\n";
        $self->create_model_metrics_for_non_qc_data();
    }
    unless ($self->qc_metrics_exist){
        print "Generating qc model metrics\n";
        $self->create_model_metrics_for_qc_data();
    }
    print "Summarizing data\n";

    my (%index_to_builds, %pool_to_builds);

    for my $model ($self->models){
        my $build = $model->last_succeeded_build or next;
        next if $model->subject->name =~ /Pooled_Library/;

        #Take the first instrument_data's index until a decision is made
        #  on how to handle per-sample data, when per-instrument-data data is unavailable
        my $index = (map{$_->index_sequence}$model->instrument_data)[0];
        my $pool = Genome::Model::Command::Services::AssignQueuedInstrumentData->_resolve_pooled_sample_name_for_instrument_data((),$model->instrument_data);

        push @{$index_to_builds{$index}}, $build;
        push @{$pool_to_builds{$pool}}, $build;
    }

    $self->subject_summary($subject_summary_file);
    $self->index_summary($index_summary_file, \%index_to_builds);
    $self->pool_summary($pool_summary_file, \%pool_to_builds);
    return 1;
}

sub common_headers {
    (
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

sub write_subject_headers {
    my $self = shift;
    my $fh = shift || die;
    print $fh join ("\t", (
            'Model',
            'Build',
            'Sample',
            'Libraries',
            'Index',
            'Pooled library',
            $self->common_headers,
        )) . "\n";
}

sub write_index_headers {
    my $self = shift;
    my $fh = shift || die;
    print $fh join ("\t", (
            'Index',
            $self->common_headers,
        )) . "\n";
}

sub write_pool_headers {
    my $self = shift;
    my $fh = shift || die;
    print $fh join ("\t", (
            'Pool',
            $self->common_headers,
        )) . "\n";
}

sub index_summary {
    my $self = shift;
    my $fh = shift || die;
    my $index_to_builds = shift || die;
    $self->write_index_headers($fh);

    while (my ($index,$builds) = each %$index_to_builds) {
        my (
            $depth_40,
            $depth_30,
            $depth_20,
            $depth_15,
            $depth_10,
            $duplication,
            $mapped,
            $uniq_on_target,
            $uniq_off_target,
            $total_on_target,
            $total_off_target,
            $unaligned,
            $snps_called,
            $with_genotype,
            $met_min_depth,
            $reference,
            $ref_match,
            $ref_was_het,
            $ref_was_hom,
            $variant,
            $var_match,
            $hom_was_het,
            $het_was_hom,
            $var_mismatch,
            $var_concordance,
            $rare_hom_concordance,
            $overall_concordance,
        ) = (0)x20;
        for my $build (@$builds) {
            my %metric = map{$_->name,$_->value}$build->metrics;
            $depth_40             += $metric{'wingspan_0_40_coverage_depth'} || 0;
            $depth_30             += $metric{'wingspan_0_30_coverage_depth'} || 0;
            $depth_20             += $metric{'wingspan_0_20_coverage_depth'} || 0;
            $depth_15             += $metric{'wingspan_0_15_coverage_depth'} || 0;
            $depth_10             += $metric{'wingspan_0_10_coverage_depth'} || 0;
            $duplication          += $metric{'wingspan_0_percent_duplication'} || 0;
            $mapped               += $metric{'wingspan_0_percent_mapped'} || 0;
            $uniq_on_target       += $metric{'wingspan_0_percent_unique_on_target'} || 0;
            $uniq_off_target      += $metric{'wingspan_0_percent_unique_off_target'} || 0;
            $total_on_target      += $metric{'wingspan_0_percent_total_on_target'} || 0;
            $total_off_target     += $metric{'wingspan_0_percent_total_off_target'} || 0;
            $unaligned            += $metric{'wingspan_0_percent_unaligned'} || 0;
            $snps_called          += $metric{'snps_called'} || 0;
            $with_genotype        += $metric{'with_genotype'} || 0;
            $met_min_depth        += $metric{'met_min_depth'} || 0;
            $reference            += $metric{'reference'} || 0;
            $ref_match            += $metric{'ref_match'} || 0;
            $ref_was_het          += $metric{'ref_was_het'} || 0;
            $ref_was_hom          += $metric{'ref_was_hom'} || 0;
            $variant              += $metric{'variant'} || 0;
            $var_match            += $metric{'var_match'} || 0;
            $hom_was_het          += $metric{'hom_was_het'} || 0;
            $het_was_hom          += $metric{'het_was_hom'} || 0;
            $var_mismatch         += $metric{'var_mismatch'} || 0;
            $var_concordance      += $metric{'var_concordance'} || 0;
            $rare_hom_concordance += $metric{'rare_hom_concordance'} || 0;
            $overall_concordance  += $metric{'overall_concordance'} || 0;
        }
        print $fh join("\t", (
                $index,
                sprintf ("%.2f",$depth_40/@$builds),
                sprintf ("%.2f",$depth_30/@$builds),
                sprintf ("%.2f",$depth_20/@$builds),
                sprintf ("%.2f",$depth_15/@$builds),
                sprintf ("%.2f",$depth_10/@$builds),
                sprintf ("%.2f%%", $duplication/@$builds),
                sprintf ("%.2f%%", $mapped/@$builds),
                sprintf ("%.2f%%", $uniq_on_target/@$builds),
                sprintf ("%.2f%%", $uniq_off_target/@$builds),
                sprintf ("%.2f%%", $total_on_target/@$builds),
                sprintf ("%.2f%%", $total_off_target/@$builds),
                sprintf ("%.2f%%", $unaligned/@$builds),
                sprintf ("%.2f",$snps_called/@$builds),
                sprintf ("%.2f",$with_genotype/@$builds),
                sprintf ("%.2f",$met_min_depth/@$builds),
                sprintf ("%.2f",$reference/@$builds),
                sprintf ("%.2f",$ref_match/@$builds),
                sprintf ("%.2f",$ref_was_het/@$builds),
                sprintf ("%.2f",$ref_was_hom/@$builds),
                sprintf ("%.2f",$variant/@$builds),
                sprintf ("%.2f",$var_match/@$builds),
                sprintf ("%.2f",$hom_was_het/@$builds),
                sprintf ("%.2f",$het_was_hom/@$builds),
                sprintf ("%.2f",$var_mismatch/@$builds),
                sprintf ("%.2f",$var_concordance/@$builds),
                sprintf ("%.2f",$rare_hom_concordance/@$builds),
                sprintf ("%.2f",$overall_concordance/@$builds),
            )) . "\n";
    }
}

sub pool_summary {
    my $self = shift;
    my $fh = shift || die;
    my $pool_to_builds = shift || die;
    $self->write_pool_headers($fh);

    while (my ($pool,$builds) = each %$pool_to_builds) {
        my (
            $depth_40,
            $depth_30,
            $depth_20,
            $depth_15,
            $depth_10,
            $duplication,
            $mapped,
            $uniq_on_target,
            $uniq_off_target,
            $total_on_target,
            $total_off_target,
            $unaligned,
            $snps_called,
            $with_genotype,
            $met_min_depth,
            $reference,
            $ref_match,
            $ref_was_het,
            $ref_was_hom,
            $variant,
            $var_match,
            $hom_was_het,
            $het_was_hom,
            $var_mismatch,
            $var_concordance,
            $rare_hom_concordance,
            $overall_concordance,
        ) = (0)x20;
        for my $build (@$builds) {
            my %metric = map{$_->name,$_->value}$build->metrics;
            $depth_40             += $metric{'wingspan_0_40_coverage_depth'} || 0;
            $depth_30             += $metric{'wingspan_0_30_coverage_depth'} || 0;
            $depth_20             += $metric{'wingspan_0_20_coverage_depth'} || 0;
            $depth_15             += $metric{'wingspan_0_15_coverage_depth'} || 0;
            $depth_10             += $metric{'wingspan_0_10_coverage_depth'} || 0;
            $duplication          += $metric{'wingspan_0_percent_duplication' || 0};
            $mapped               += $metric{'wingspan_0_percent_mapped'} || 0;
            $uniq_on_target       += $metric{'wingspan_0_percent_unique_on_target'} || 0;
            $uniq_off_target      += $metric{'wingspan_0_percent_unique_off_target'} || 0;
            $total_on_target      += $metric{'wingspan_0_percent_total_on_target'} || 0;
            $total_off_target     += $metric{'wingspan_0_percent_total_off_target'} || 0;
            $unaligned            += $metric{'wingspan_0_percent_unaligned'} || 0;
            $snps_called          += $metric{'snps_called'} || 0;
            $with_genotype        += $metric{'with_genotype'} || 0;
            $met_min_depth        += $metric{'met_min_depth'} || 0;
            $reference            += $metric{'reference'} || 0;
            $ref_match            += $metric{'ref_match'} || 0;
            $ref_was_het          += $metric{'ref_was_het'} || 0;
            $ref_was_hom          += $metric{'ref_was_hom'} || 0;
            $variant              += $metric{'variant'} || 0;
            $var_match            += $metric{'var_match'} || 0;
            $hom_was_het          += $metric{'hom_was_het'} || 0;
            $het_was_hom          += $metric{'het_was_hom'} || 0;
            $var_mismatch         += $metric{'var_mismatch'} || 0;
            $var_concordance      += $metric{'var_concordance'} || 0;
            $rare_hom_concordance += $metric{'rare_hom_concordance'} || 0;
            $overall_concordance  += $metric{'overall_concordance'} || 0;
        }
        print $fh join("\t", (
                $pool,
                sprintf ("%.2f",$depth_40/@$builds),
                sprintf ("%.2f",$depth_30/@$builds),
                sprintf ("%.2f",$depth_20/@$builds),
                sprintf ("%.2f",$depth_15/@$builds),
                sprintf ("%.2f",$depth_10/@$builds),
                sprintf ("%.2f%%", $duplication/@$builds),
                sprintf ("%.2f%%", $mapped/@$builds),
                sprintf ("%.2f%%", $uniq_on_target/@$builds),
                sprintf ("%.2f%%", $uniq_off_target/@$builds),
                sprintf ("%.2f%%", $total_on_target/@$builds),
                sprintf ("%.2f%%", $total_off_target/@$builds),
                sprintf ("%.2f%%", $unaligned/@$builds),
                sprintf ("%.2f",$snps_called/@$builds),
                sprintf ("%.2f",$with_genotype/@$builds),
                sprintf ("%.2f",$met_min_depth/@$builds),
                sprintf ("%.2f",$reference/@$builds),
                sprintf ("%.2f",$ref_match/@$builds),
                sprintf ("%.2f",$ref_was_het/@$builds),
                sprintf ("%.2f",$ref_was_hom/@$builds),
                sprintf ("%.2f",$variant/@$builds),
                sprintf ("%.2f",$var_match/@$builds),
                sprintf ("%.2f",$hom_was_het/@$builds),
                sprintf ("%.2f",$het_was_hom/@$builds),
                sprintf ("%.2f",$var_mismatch/@$builds),
                sprintf ("%.2f",$var_concordance/@$builds),
                sprintf ("%.2f",$rare_hom_concordance/@$builds),
                sprintf ("%.2f",$overall_concordance/@$builds),
            )) . "\n";
    }
}

sub subject_summary {
    my $self = shift;
    my $fh = shift || die;
    $self->write_subject_headers($fh);

    for my $model ($self->models){
        next if $model->subject->name =~ /Pooled_Library/;
        my $build = $model->last_succeeded_build || next;

        my %metric = map{$_->name,$_->value}$build->metrics;

        #Take the first instrument_data's index until a decision is made
        #  on how to handle per-sample data, when per-instrument-data data is unavailable
        my $index = (map{$_->index_sequence}$model->instrument_data)[0];
        my $pool = Genome::Model::Command::Services::AssignQueuedInstrumentData->_resolve_pooled_sample_name_for_instrument_data((),$model->instrument_data);
        my $libraries = $build->instrument_data->library->name;

        #add model->instrument_data->lane
        print $fh join("\t", (
                $model->id,
                $build->id,
                $model->subject->name,
                $libraries,
                $index,
                $pool,
                $metric{'wingspan_0_40_coverage_depth'} || '-',
                $metric{'wingspan_0_30_coverage_depth'} || '-',
                $metric{'wingspan_0_20_coverage_depth'} || '-',
                $metric{'wingspan_0_15_coverage_depth'} || '-',
                $metric{'wingspan_0_10_coverage_depth'} || '-',
                $metric{'wingspan_0_percent_duplication'} ? sprintf ("%.1f%%", $metric{'wingspan_0_percent_duplication'}) : '-',
                sprintf ("%.2f%%", $metric{'wingspan_0_percent_mapped'}) || '-',
                sprintf ("%.2f%%", $metric{'wingspan_0_percent_unique_on_target'}) || '-',
                sprintf ("%.2f%%", $metric{'wingspan_0_percent_unique_off_target'}) || '-',
                sprintf ("%.2f%%", $metric{'wingspan_0_percent_total_on_target'}) || '-',
                sprintf ("%.2f%%", $metric{'wingspan_0_percent_total_off_target'}) || '-',
                sprintf ("%.2f%%", $metric{'wingspan_0_percent_unaligned'}) || '-',
                $metric{'snps_called'} || '-',
                $metric{'with_genotype'} || '-',
                $metric{'met_min_depth'} || '-',
                $metric{'reference'} || '-',
                $metric{'ref_match'} || '-',
                $metric{'ref_was_het'} || '-',
                $metric{'ref_was_hom'} || '-',
                $metric{'variant'} || '-',
                $metric{'var_match'} || '-',
                $metric{'hom_was_het'} || '-',
                $metric{'het_was_hom'} || '-',
                $metric{'var_mismatch'} || '-',
                $metric{'var_concordance'} || '-',
                $metric{'rare_hom_concordance'} || '-',
                $metric{'overall_concordance'} || '-',
            )) . "\n";
    }
}

sub create_model_metrics_for_non_qc_data {
    my $self = shift;
    for my $model ($self->models){
        next if grep {$_->index_sequence eq 'unknown'} $model->instrument_data;
        my $build = $model->last_succeeded_build or next;
        my %metrics = map{$_->name,$_->value}$build->metrics;

        my $align_stats;
        my @result_users = $build->result_users('software_result.subclass_name' => 'Genome::InstrumentData::AlignmentResult::Merged::CoverageStats');
        if ($build->alignment_summary_hash_ref) {
            $align_stats = $build->alignment_summary_hash_ref->{0};
        } elsif (@result_users) {
            #There will be multiple coverage stats for somatic validation; one per bam: tumor and normal (aka control)
            #Somatic validation isn't an expected candidate for this tool yet
            $align_stats = $result_users[0]->software_result->alignment_summary_hash_ref->{0}
        } else {
            print "Build " . $build->id . " with sample " . $model->subject_name . " has no alignment stats.\n";
            next;
        }

        my $unique_on_target = $metrics{wingspan_0_unique_target_aligned_bp} || $align_stats->{unique_target_aligned_bp};
        my $duplicate_on_target = $metrics{wingspan_0_duplicate_target_aligned_bp} || $align_stats->{duplicate_target_aligned_bp};
        my $unique_off_target = $metrics{wingspan_0_unique_off_target_aligned_bp} || $align_stats->{unique_off_target_aligned_bp};
        my $duplicate_off_target = $metrics{wingspan_0_duplicate_off_target_aligned_bp} || $align_stats->{duplicate_off_target_aligned_bp};
        my $unaligned = $metrics{wingspan_0_total_unaligned_bp} || $align_stats->{total_unaligned_bp};

        my $total = $unique_on_target+$duplicate_on_target+$unique_off_target+$duplicate_off_target+$unaligned;
        my $percent_unique_on_target = $unique_on_target/$total*100;
        my $percent_duplicate_on_target = $duplicate_on_target/$total*100;
        my $percent_unique_off_target = $unique_off_target/$total*100;
        my $percent_duplicate_off_target = $duplicate_off_target/$total*100;
        my $percent_unaligned = $unaligned/$total*100;
        my $percent_total_on_target = $percent_duplicate_on_target + $percent_unique_on_target;
        my $percent_total_off_target = $percent_duplicate_off_target + $percent_unique_off_target;

        my $stats = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($build->whole_rmdup_bam_flagstat_file);
        $build->add_metric(name => 'wingspan_0_percent_mapped', value => $stats->{reads_mapped_percentage});

        my $cov_stats = $build->coverage_stats_summary_hash_ref->{0};
        #Default coverage depth percentage is zero - there may not be depth information, especially for 30x and 40x
        $build->add_metric(name => 'wingspan_0_40_coverage_depth', value => ($cov_stats->{40}{pc_target_space_covered} || 0));
        $build->add_metric(name => 'wingspan_0_30_coverage_depth', value => ($cov_stats->{30}{pc_target_space_covered} || 0));
        $build->add_metric(name => 'wingspan_0_20_coverage_depth', value => ($cov_stats->{20}{pc_target_space_covered} || 0));
        $build->add_metric(name => 'wingspan_0_15_coverage_depth', value => ($cov_stats->{15}{pc_target_space_covered} || 0));
        $build->add_metric(name => 'wingspan_0_10_coverage_depth', value => ($cov_stats->{10}{pc_target_space_covered} || 0));
        if (my ($dup) = map {$metrics{$_}} grep {$_ =~ /PERCENT_DUPLICATION/} keys %metrics){
            $build->add_metric(name => 'wingspan_0_percent_duplication', value => $dup * 100);
        }
        else{
            $build->add_metric(name => 'wingspan_0_percent_duplication', value => 0);
        }
        $build->add_metric(name => 'wingspan_0_percent_unique_on_target', value => $percent_unique_on_target);
        $build->add_metric(name => 'wingspan_0_percent_unique_off_target', value => $percent_unique_off_target);
        $build->add_metric(name => 'wingspan_0_percent_total_on_target', value => $percent_total_on_target);
        $build->add_metric(name => 'wingspan_0_percent_total_off_target', value => $percent_total_off_target);
        $build->add_metric(name => 'wingspan_0_percent_unaligned', value => $percent_unaligned);
    }
}

sub create_model_metrics_for_qc_data {
    my $self = shift;
    my $qc_dir = $self->qc_directory;

    #Find each of the directories, named after the samples, with qc files in them
    my @dir_names = `find $qc_dir -maxdepth 1 -type d -printf %P\\\\n`;

    for my $dir_name (@dir_names){
        chomp $dir_name;
        my $sample = $dir_name;
        my ($qc_file) = `ls $qc_dir/$dir_name/*.qc 2>/dev/null`;
        next unless $qc_file;
        chomp $qc_file;

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

        my $model = Genome::Model->get(id => [map{$_->id}$self->models], subject_name => $sample) || next;
        my $build = $model->last_succeeded_build;
        my %metrics = map{$_->name,$_->value}$build->metrics;

        $build->add_metric(name => 'snps_called', value => $values[1] || 0);
        $build->add_metric(name => 'with_genotype', value => $values[2] || 0);
        $build->add_metric(name => 'met_min_depth', value => $values[3] || 0);
        $build->add_metric(name => 'reference', value => $values[4] || 0);
        $build->add_metric(name => 'ref_match', value => $values[5] || 0);
        $build->add_metric(name => 'ref_was_het', value => $values[6] || 0);
        $build->add_metric(name => 'ref_was_hom', value => $values[7] || 0);
        $build->add_metric(name => 'variant', value => $values[8] || 0);
        $build->add_metric(name => 'var_match', value => $values[9] || 0);
        $build->add_metric(name => 'hom_was_het', value => $values[10] || 0);
        $build->add_metric(name => 'het_was_hom', value => $values[11] || 0);
        $build->add_metric(name => 'var_mismatch', value => $values[12] || 0);
        $build->add_metric(name => 'var_concordance', value => $values[13] || 0);
        $build->add_metric(name => 'rare_hom_concordance', value => $values[14] || 0);
        $build->add_metric(name => 'overall_concordance', value => $values[15] || 0);
    }
}

sub qc_metrics_exist {
    my $self = shift;
    for my $model ($self->models){
        next if grep {$_->index_sequence eq 'unknown'} $model->instrument_data;
        my $build = $model->last_succeeded_build or next;
        my %metrics = map{$_->name,$_->value}$build->metrics;
        return 0 unless(
            defined($metrics{'snps_called'}) and
            defined($metrics{'with_genotype'}) and
            defined($metrics{'met_min_depth'}) and
            defined($metrics{'reference'}) and
            defined($metrics{'ref_match'}) and
            defined($metrics{'ref_was_het'}) and
            defined($metrics{'ref_was_hom'}) and
            defined($metrics{'variant'}) and
            defined($metrics{'var_match'}) and
            defined($metrics{'hom_was_het'}) and
            defined($metrics{'het_was_hom'}) and
            defined($metrics{'var_mismatch'}) and
            defined($metrics{'var_concordance'}) and
            defined($metrics{'rare_hom_concordance'}) and
            defined($metrics{'overall_concordance'})
        );
    }
    return 1;
}

sub non_qc_metrics_exist {
    my $self = shift;
    for my $model ($self->models){
        next if grep {$_->index_sequence eq 'unknown'} $model->instrument_data;
        my $build = $model->last_succeeded_build or next;
        my %metrics = map{$_->name,$_->value}$build->metrics;
        return 0 unless(
            defined($metrics{'wingspan_0_40_coverage_depth'}) and
            defined($metrics{'wingspan_0_30_coverage_depth'}) and
            defined($metrics{'wingspan_0_20_coverage_depth'}) and
            defined($metrics{'wingspan_0_15_coverage_depth'}) and
            defined($metrics{'wingspan_0_10_coverage_depth'}) and
            defined($metrics{'wingspan_0_percent_duplication'}) and
            defined($metrics{'wingspan_0_percent_mapped'}) and
            defined($metrics{'wingspan_0_percent_unique_on_target'}) and
            defined($metrics{'wingspan_0_percent_unique_off_target'}) and
            defined($metrics{'wingspan_0_percent_total_on_target'}) and
            defined($metrics{'wingspan_0_percent_total_off_target'}) and
            defined($metrics{'wingspan_0_percent_unaligned'})
        );
    }
    return 1;
}

sub delete_model_metrics {
    my $self = shift;
    for my $model ($self->models){
        next if grep {$_->index_sequence eq 'unknown'} $model->instrument_data;
        my $build = $model->last_succeeded_build or next;
        $build->metrics(name => 'wingspan_0_40_coverage_depth')->delete if $build->metrics(name => 'wingspan_0_40_coverage_depth');
        $build->metrics(name => 'wingspan_0_30_coverage_depth')->delete if $build->metrics(name => 'wingspan_0_30_coverage_depth');
        $build->metrics(name => 'wingspan_0_20_coverage_depth')->delete if $build->metrics(name => 'wingspan_0_20_coverage_depth');
        $build->metrics(name => 'wingspan_0_15_coverage_depth')->delete if $build->metrics(name => 'wingspan_0_15_coverage_depth');
        $build->metrics(name => 'wingspan_0_10_coverage_depth')->delete if $build->metrics(name => 'wingspan_0_10_coverage_depth');
        $build->metrics(name => 'wingspan_0_percent_duplication')->delete if $build->metrics(name => 'wingspan_0_percent_duplication');
        $build->metrics(name => 'wingspan_0_percent_mapped')->delete if $build->metrics(name => 'wingspan_0_percent_mapped');
        $build->metrics(name => 'wingspan_0_percent_unique_on_target')->delete if $build->metrics(name => 'wingspan_0_percent_unique_on_target');
        $build->metrics(name => 'wingspan_0_percent_unique_off_target')->delete if $build->metrics(name => 'wingspan_0_percent_unique_off_target');
        $build->metrics(name => 'wingspan_0_percent_total_on_target')->delete if $build->metrics(name => 'wingspan_0_percent_total_on_target');
        $build->metrics(name => 'wingspan_0_percent_total_off_target')->delete if $build->metrics(name => 'wingspan_0_percent_total_off_target');
        $build->metrics(name => 'wingspan_0_percent_unaligned')->delete if $build->metrics(name => 'wingspan_0_percent_unaligned');
        $build->metrics(name => 'snps_called')->delete if $build->metrics(name => 'snps_called');
        $build->metrics(name => 'with_genotype')->delete if $build->metrics(name => 'with_genotype');
        $build->metrics(name => 'met_min_depth')->delete if $build->metrics(name => 'met_min_depth');
        $build->metrics(name => 'reference')->delete if $build->metrics(name => 'reference');
        $build->metrics(name => 'ref_match')->delete if $build->metrics(name => 'ref_match');
        $build->metrics(name => 'ref_was_het')->delete if $build->metrics(name => 'ref_was_het');
        $build->metrics(name => 'ref_was_hom')->delete if $build->metrics(name => 'ref_was_hom');
        $build->metrics(name => 'variant')->delete if $build->metrics(name => 'variant');
        $build->metrics(name => 'var_match')->delete if $build->metrics(name => 'var_match');
        $build->metrics(name => 'hom_was_het')->delete if $build->metrics(name => 'hom_was_het');
        $build->metrics(name => 'het_was_hom')->delete if $build->metrics(name => 'het_was_hom');
        $build->metrics(name => 'var_mismatch')->delete if $build->metrics(name => 'var_mismatch');
        $build->metrics(name => 'var_concordance')->delete if $build->metrics(name => 'var_concordance');
        $build->metrics(name => 'rare_hom_concordance')->delete if $build->metrics(name => 'rare_hom_concordance');
        $build->metrics(name => 'overall_concordance')->delete if $build->metrics(name => 'overall_concordance');
    }
}
