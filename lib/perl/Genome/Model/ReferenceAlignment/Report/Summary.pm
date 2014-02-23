#:boberkfe it would be nice if the other reports stored summary data as build metrics so that this
#:boberkfe could generate without having to parse out the individual files in the report paths.

package Genome::Model::ReferenceAlignment::Report::Summary;

use strict;
use warnings;

use Genome;

use IO::String;
use Template;

my $base_template_path = __PACKAGE__->_base_path_for_templates;

class Genome::Model::ReferenceAlignment::Report::Summary {
    is => 'Genome::Model::Report',
    has => [
        report_templates => {
            is => 'String',
            is_many => 1,
            default_value => [
                 "$base_template_path.html.tt2",
                 "$base_template_path.txt.tt2"
            ],
            doc => 'The paths of template(s) to use to format the report.  (In .tt2 format)',
        },
        name => {
            default_value => 'Summary',
        },
        description => {
            default_value => "Link to summary report will go here",
        },
    ],
};

# TODO: move up into base class
sub _base_path_for_templates
{
    my $module = __PACKAGE__;
    $module =~ s/::/\//g;
    $module .= '.pm';
    my $module_path = $INC{$module};
    unless ($module_path) {
        die "Module " . __PACKAGE__ . " failed to find its own path!  Checked for $module in \%INC...";
    }
    return $module_path;
}

sub _add_to_report_xml
{
    my $self = shift;
    my $template = shift;

    my @templates = $self->report_templates;
    unless (@templates) {
        die "No report templates assigned!  Cannot generate any content."
    }

    #my $data = { description => $self->generate_report_brief };
    my $data = {};

    for my $template (@templates) {
        my $content = $self->generate_report_detail($template);
        my ($format,$key);
        if ($content =~ /\<\s*HTML/i) {
            $format = 'HTML';
            $key = 'html';
        }
        else {
            $format = 'text';
            $key = 'txt';
        }
        if (exists $data->{$key}) {
            die "Multiple templates return content in $format format.  This is not supported, sadly."
                . "  Error processing $template";
        }
        $data->{$key} = $content;
    };
    return $data;
}

sub generate_report_brief
{
    my $self=shift;
    return "Link to summary report will go here";
}

sub generate_report_detail
{
    my $self = shift;
    my $template = shift;
    unless ($template) {
        die "please specify which template to use for this report!";
    }

    my $model = $self->model;
    my $build = $self->build;

    $self->debug_message("Running report summary for build ".$build->id.".");
    my $body = IO::String->new();
    die $! unless $body;
    my $summary = $self->get_summary_information($template);
    $body->print($summary);
    $body->seek(0, 0);
    return join('', $body->getlines);
}

sub get_summary_information
{
    my $self = shift;
    my $template = shift;
    unless ($template) {
        die "please specify which template to use for this report!";
    }

    my $build = $self->build;
    my $model = $build->model;

    my $content;

    #################################
    my $na = "Not Available";

    my $haploid_coverage=$na;

    my $total_unfiltered_snps=$na;
    my $total_filtered_snps=$na;

    #my $unfiltered_dbsnp_positions=$na;
    #my $filtered_dbsnp_positions=$na;

    my $unfiltered_dbsnp_concordance=$na;
    my $filtered_dbsnp_concordance=$na;

    my $report_dir = $build->resolve_reports_directory;

    my $mapcheck_report_file = $report_dir."/Mapcheck/report.html";
    my $goldsnp_report_file = $report_dir."/Gold_SNP_Concordance/report.html";
    my $dbsnp_report_file = $build->dbsnp_file_filtered;
    my $input_base_count_report_file = $report_dir . "/Input_Base_Count/report.html";

    ##match mapcheck report
    my $fh = new IO::File($mapcheck_report_file, "r");
    if ($fh) {
        my $mapcheck_contents = get_contents($fh);
        if ( ($mapcheck_contents =~ m/Average depth across all non-gap regions: (\S+)/g ) || ($mapcheck_contents =~ m/\nAverage Coverage:(\S+)/g ) ) {
            $haploid_coverage=$1 if defined($1);
            if ($haploid_coverage) {
                $build->set_metric( 'haploid_coverage', $haploid_coverage );
            }
        }
        $fh->close();
    } else {
        $self->debug_message("Could not locate RefSeqMaq report at $mapcheck_report_file!");
    }

    ##match goldsnp report
    my %gold_snp_metrics = $self->format_gold_snp_metrics($build);
    if (!%gold_snp_metrics and $build->gold_snp_build) {
        $self->debug_message("Could not generate gold snp metrics, regenerating gold snp concordance and trying again");
        Genome::Model::ReferenceAlignment::Command::CreateMetrics::GoldSnpConcordance->execute(build => $build);
        %gold_snp_metrics = $self->format_gold_snp_metrics($build);
    }

    my $unfiltered_diploid_het_coverage_actual_number = $gold_snp_metrics{unfiltered}{het_hits} || $na;
    my $unfiltered_diploid_het_coverage_percent = $gold_snp_metrics{unfiltered}{het_percent} || $na;
    my $unfiltered_diploid_hom_coverage_actual_number = $gold_snp_metrics{unfiltered}{hom_hits} || $na;
    my $unfiltered_diploid_hom_coverage_percent = $gold_snp_metrics{unfiltered}{hom_percent} || $na;

    my $filtered_diploid_het_coverage_actual_number = $gold_snp_metrics{filtered}{het_hits} || $na;
    my $filtered_diploid_het_coverage_percent = $gold_snp_metrics{filtered}{het_percent} || $na;
    my $filtered_diploid_hom_coverage_actual_number = $gold_snp_metrics{filtered}{hom_hits} || $na;
    my $filtered_diploid_hom_coverage_percent = $gold_snp_metrics{filtered}{hom_percent} || $na;

    ##match dbsnp report
    my $dbsnp_filtered_report_file = $build->dbsnp_file_filtered;
    my $dbsnp_unfiltered_report_file = $build->dbsnp_file_unfiltered;
    unless (-e $dbsnp_filtered_report_file and -e $dbsnp_filtered_report_file) {
        $self->debug_message("Did not find dbSNP concordance files, rerunning dbSNP concordance");
        Genome::Model::ReferenceAlignment::Command::CreateMetrics::DbSnpConcordance->execute(build => $build);
    }

    if (-e $dbsnp_filtered_report_file and -e $dbsnp_unfiltered_report_file) {
        my $dbsnp_filtered_data = Genome::Model::Tools::Joinx::SnvConcordanceByQuality::parse_results_file($dbsnp_filtered_report_file);
        my $dbsnp_unfiltered_data = Genome::Model::Tools::Joinx::SnvConcordanceByQuality::parse_results_file($dbsnp_unfiltered_report_file);

        # get unfiltered data
        if (exists $dbsnp_unfiltered_data->{total_snvs}) {
            $total_unfiltered_snps = $dbsnp_unfiltered_data->{total_snvs};
        }
        else {
            $self->debug_message("Could not extract total unfiltered SNPs from $dbsnp_unfiltered_report_file!");
        }

        if (exists $dbsnp_unfiltered_data->{total_concordance}) {
            $unfiltered_dbsnp_concordance = $dbsnp_unfiltered_data->{total_concordance};
        }
        else {
            $self->debug_message("Could not extract unfiltered concordance from $dbsnp_unfiltered_report_file!");
        }

        # get filtered data
        if (exists $dbsnp_filtered_data->{total_snvs}) {
            $total_filtered_snps = $dbsnp_filtered_data->{total_snvs};
        }
        else {
            $self->debug_message("Could not extract total filtered SNPs from $dbsnp_filtered_report_file!");
        }

        if (exists $dbsnp_filtered_data->{total_concordance}) {
            $filtered_dbsnp_concordance = $dbsnp_filtered_data->{total_concordance};
        } 
        else {
            $self->debug_message("Could not extract filtered concordance from $dbsnp_filtered_report_file!");
        }
    }
    else {
        $self->debug_message("dbSNP filtered report: $dbsnp_filtered_report_file is not available") unless -e $dbsnp_filtered_report_file;
        $self->debug_message("dbSNP unfiltered report: $dbsnp_unfiltered_report_file is not available") unless -e $dbsnp_unfiltered_report_file;
    }

    my @inputs = $build->model->instrument_data_inputs;
    my $total_bases = 0;
    for my $input (@inputs) {
        my $inst_data = $input->value;
        if ($inst_data->can('total_bases_read'))  {
            $total_bases += $inst_data->total_bases_read($input->filter_desc);
        }
    }
    my $total_gigabases = sprintf("%.03f", $total_bases/1000000000);

    if ($model->read_trimmer_name and $model->read_trimmer_name =~ /^trimq2/) {
        my ($total_ct, $total_trim_ct) = $build->calculate_input_base_counts_after_trimq2;
        if ($total_ct and $total_trim_ct) {
            my $gb       = sprintf("%.03f", $total_ct/1000000000);
            my $trim_gb  = sprintf("%.03f", $total_trim_ct/1000000000);
            $total_gigabases = "$trim_gb/$gb";
        }
        else {
            $self->warning_message("Failed to get input base counts after trimq2");
        }
    }

    # summarize the instrument data
    my %library_lane_counts;

    my @inst_data = map { $_->value } @inputs;
    unless ($model->read_aligner_name =~ /Imported$/i) {
        my %library_lanes;
        for my $i (@inst_data) {
            my $library_name = $i->library_name;
            my $a = $library_lanes{$library_name} ||= [];
            push @$a, $i->run_name . "/" . $i->subset_name
        }
        for my $library_name (keys %library_lanes) {
            $library_lane_counts{$library_name} = scalar(@{ $library_lanes{$library_name} })
        }
    }

    # sample variables
    my $sample;
    if($model->subject_type eq 'sample_name' or $model->subject_type eq 'genomic_dna') {
        $sample = $model->subject;
    } elsif ($model->subject_type eq 'library_name') {
        my $library = $model->subject;
        if($library) {
            $sample = $library->sample;
        }
    }

    my ($extraction_label,$tissue_label,$extraction_name,$extraction_id,$extraction_desc,$extraction_type) = ($na,$na,$na,$na,$na,$na);
    if ($sample) {
        $tissue_label = $sample->tissue_label || $na;

        $extraction_label = $sample->extraction_label || $na;
        $extraction_name  = $sample->name || $na;
        $extraction_id    = $sample->id || $na;
        #$extraction_desc  = $sample->description;
        $extraction_type  = $sample->sample_type || $na;
    }
    else {
        $self->warning_message("No sample found for " . $model->subject_name);
    }

    # patient variables
    my $source;
    my ($source_upn,$source_desc) = ($na,$na);
    $source = $sample->source if $sample;
    if ($source) {
        $source_upn = $source->name || $na;
        $source_desc = $source->description || $na;
        #$source_gender = $source->gender;
    }
    else {
        $self->warning_message("No source individual/population found for sample!");
    }

    my $taxon;
    my $taxon_id;
    my $species = $na;
    my $species_latin_name = $na;
    
    if($sample) {
        $taxon = $sample->taxon;

        if ($taxon) {
            $taxon_id = $taxon->taxon_id;
            $species = $taxon->species_name;
            $species_latin_name = $taxon->species_latin_name;
        } else {
            $self->warning_message("No taxon found for sample!");
        }
    }

    my $ref_seq_dir = $self->model->reference_sequence_build->data_directory;

    # processing profile
    my $pp = $model->processing_profile;

    #my @filtered_files = $build->get_variant_bed_file('snvs.hq');
    my @filtered_files = $build->filtered_snvs_bed;
    my $filtered_snp_calls = `wc -l @filtered_files | tail -n 1`;
    $filtered_snp_calls =~ s/\s\S+\s*$//i;
    $filtered_snp_calls =~ s/\s//g;

    my @lq_files = $build->get_variant_bed_file('snvs.lq');
    unless (@lq_files) {
        @lq_files = $build->get_variant_bed_file('snps_all_sequences');
    }
    my $lq_snp_calls = `wc -l @lq_files | tail -n 1`;
    $lq_snp_calls =~ s/\s\S+\s*$//i;
    $lq_snp_calls =~ s/\s//g;
    my $unfiltered_snp_calls = $filtered_snp_calls + $lq_snp_calls;

    my $snp_chromosomes = $self->model->reference_sequence_build->description;
    my $snp_caller = $self->model->snv_detection_strategy;

    my @stat = stat($filtered_files[-1]);
    my $time = POSIX::strftime("%Y-%m-%d %H:%M:%S", localtime($stat[10]));

    my $model_name = $model->name;
    my $build_id = $build->id;
    my $data_directory = $build->data_directory . "/";

    my @vars = (
        model_id                                      => $model->id,
        model_name                                    => $model->name,

        patient_upn                                   => $source_upn,

        taxon_id                                      => $taxon_id,
        species                                       => $species,
        species_latin_name                            => $species_latin_name
    );

    # ehvatum TODO: When ReferecePlaceholder is deleted, remove this if statement and always use ref_seq_build_id,
    # ref_seq_prefix, ref_seq_subject_name, and, if set, ref_seq_version.
    if(defined($model->reference_sequence_build))
    {
        my $refSeqDesc = $self->model->reference_sequence_build->prefix . '-' . $self->model->reference_sequence_build->subject_name;
        if(defined($self->model->reference_sequence_build->version))
        {
            $refSeqDesc .= '-' . $self->model->reference_sequence_build->version;
        }
        push @vars, (ref_seq_desc                     => $refSeqDesc);
        push @vars, (ref_seq_build_id                 => $self->model->reference_sequence_build->build_id);
    }
    else
    {
        push @vars, (ref_seq_name                     => $self->model->reference_sequence_build->name);
    }

    push @vars, (
        ref_seq_dir                                   => $ref_seq_dir,

        tissue_sample_label                           => $tissue_label,

        extraction_label                              => $extraction_label,
        extraction_type                               => $extraction_type,
        extraction_name                               => $extraction_name,
        extraction_id                                 => $extraction_id,
        extraction_desc                               => $extraction_desc,

        processing_profile_type                       => $pp->type_name,
        processing_profile_name                       => $pp->name,
        processing_profile_id                         => $pp->id,

        build_id                                      => $build->id,
        build_date                                    => $time,
        data_directory                                => $data_directory,

        total_number_of_lanes                         => scalar(@inst_data),
        total_gigabases                               => $total_gigabases,
        libraries                                     => [ sort keys %library_lane_counts ],
        lanes_by_library                              => \%library_lane_counts,

        haploid_coverage                              => $haploid_coverage,

        unfiltered_snp_calls                          => commify($unfiltered_snp_calls),
        filtered_snp_calls                            => commify($filtered_snp_calls),

        snp_chromosomes                               => $snp_chromosomes,
        snp_caller                                    => $snp_caller,

        total_filtered_snps                           => commify($total_filtered_snps),
        total_unfiltered_snps                         => commify($total_unfiltered_snps),

        #unfiltered_dnsbp_positions                    => commify($unfiltered_dbsnp_positions),
        #filtered_dnsbp_positions                      => commify($filtered_dbsnp_positions),

        unfiltered_dbsnp_concordance                  => $unfiltered_dbsnp_concordance,
        filtered_dbsnp_concordance                    => $filtered_dbsnp_concordance,

        unfiltered_diploid_het_coverage_actual_number => commify($unfiltered_diploid_het_coverage_actual_number),
        unfiltered_diploid_het_coverage_percent       => $unfiltered_diploid_het_coverage_percent,
        unfiltered_diploid_hom_coverage_actual_number => commify($unfiltered_diploid_hom_coverage_actual_number),
        unfiltered_diploid_hom_coverage_percent       => $unfiltered_diploid_hom_coverage_percent,

        filtered_diploid_het_coverage_actual_number   => commify($filtered_diploid_het_coverage_actual_number),
        filtered_diploid_het_coverage_percent         => $filtered_diploid_het_coverage_percent,
        filtered_diploid_hom_coverage_actual_number   => commify($filtered_diploid_hom_coverage_actual_number),
        filtered_diploid_hom_coverage_percent         => $filtered_diploid_hom_coverage_percent,

        view_url => $ENV{GENOME_SYS_SERVICES_WEB_VIEW_URL},
    );

    ##################################

    my $tt = Template->new({
         ABSOLUTE => 1,
    }) || die "$Template::ERROR\n";

    my $varstest = {
        name     => 'Mickey',
        debt     => '3 riffs and a solo',
        deadline => 'the next chorus',
        files_url => $ENV{GENOME_SYS_SERVICES_FILES_URL},
    };

    $self->debug_message("processing template $template");

    my $rv = $tt->process($template, { @vars }, \$content) || die $tt->error(), "\n";
    if ($rv != 1) {
   	    die "Bad return value from template processing for summary report generation: $rv ";
    }
    unless ($content) {
        die "No content returned from template processing!";
    }

    return $content;
}

sub get_contents {
   my $in = shift;
   my $ret = "";
   while (<$in>) {
      $ret.= $_;
   }
   return $ret;
}

sub commify {
	local $_  = shift;
	1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
	return $_;
}

sub format_gold_snp_metrics {
    my ($self, $build) = @_;

    my $filtered_file = $build->gold_snp_report_file_filtered;
    my $unfiltered_file = $build->gold_snp_report_file_unfiltered;
    if (! -s $filtered_file or ! -s $unfiltered_file) {
        # ... no gold snp data!
        return;
    }

    my $gold_snp_raw = {
        unfiltered => Genome::Model::Tools::Joinx::SnvConcordance::parse_results_file($unfiltered_file),
        filtered => Genome::Model::Tools::Joinx::SnvConcordance::parse_results_file($filtered_file),
    };

    my %metrics;
    for my $type (qw/filtered unfiltered/) {
        my $het_data = $gold_snp_raw->{$type}{'heterozygous snv'};
        my $het_hits = $het_data->{hits}{match}{'heterozygous (1 alleles) snv'}{count};
        my $het_percent = ($het_data->{total} == 0) ? 0 : # don't divide by zero
                    sprintf("%.02f", 100.0 * $het_hits / $het_data->{total});

        my $hom_data = $gold_snp_raw->{$type}{'homozygous snv'};
        my $hom_hits = $hom_data->{hits}{match}{'homozygous snv'}{count};
        my $hom_percent = ($hom_data->{total} == 0) ? 0 : # don't divide by zero
                    sprintf("%.02f", 100.0 * $hom_hits / $hom_data->{total});

        $metrics{$type} = {
            het_hits => $het_hits,
            het_percent => $het_percent,
            hom_hits => $hom_hits,
            hom_percent => $hom_percent,
        };
    }

    return %metrics;
}
1;
