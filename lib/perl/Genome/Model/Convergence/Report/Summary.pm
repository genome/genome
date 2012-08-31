package Genome::Model::Convergence::Report::Summary;

use strict;
use warnings;

use Genome;
use Genome::Info::BamFlagstat;

class Genome::Model::Convergence::Report::Summary{
    is => 'Genome::Model::Report',
};

sub description {
    my $self = shift();

    return 'summary of harmonic convergence model';
}

sub _add_to_report_xml {
    my $self = shift();

    my $build = $self->build;

    my $doc = $self->_xml;

    #List models/builds in this model
    my $members_node = $doc->createElement('members');

    for my $member ($build->members) {
        my $member_node = $members_node->addChild( $doc->createElement('member') );

        $member_node->addChild( $doc->createAttribute('model-id', $member->model->id) );
        $member_node->addChild( $doc->createAttribute('name', $member->model->name) );
        $member_node->addChild( $doc->createAttribute('build-id', $member->id) );
        $member_node->addChild( $doc->createAttribute('completed', $member->date_completed));
    }
    $self->_main_node->addChild($members_node);

    #various reference alignment metrics
    my $metrics_node = $doc->createElement('metrics');

    my @all_subbuilds = $build->all_subbuilds_closure;

    for my $subbuild (@all_subbuilds) {
        next unless $subbuild->type_name eq 'reference alignment';

        my $build_node = $metrics_node->addChild( $doc->createElement('build') );

        my $data = $self->extract_refalign_report_data($subbuild);

        $build_node->addChild( $doc->createAttribute('build-id', $subbuild->id) );
        $build_node->addChild( $doc->createAttribute('name', $subbuild->model->name) );
        $build_node->addChild( $doc->createAttribute('lanes', $data->{lane_count}));
        $build_node->addChild( $doc->createAttribute('haploid-coverage', $subbuild->get_metric('haploid_coverage')) );
        $build_node->addChild( $doc->createAttribute('instrument-data-total-kb', $subbuild->get_metric('instrument data total kb')) );

        $build_node->addChild( $doc->createAttribute('unfiltered-snp-calls', $data->{unfiltered_snp_calls} || '-'));
        $build_node->addChild( $doc->createAttribute('unfiltered-dbsnp-concordance', $data->{unfiltered_dbsnp_concordance} || '-'));
        $build_node->addChild( $doc->createAttribute('filtered-snp-calls', $data->{filtered_snp_calls} || '-'));
        $build_node->addChild( $doc->createAttribute('filtered-dbsnp-concordance', $data->{filtered_dbsnp_concordance} || '-'));

        $build_node->addChild( $doc->createAttribute('unfiltered-diploid-heterozygous-percentage', $data->{unfiltered_diploid_heterozygous_percentage} || '-'));
        $build_node->addChild( $doc->createAttribute('filtered-diploid-heterozygous-percentage', $data->{filtered_diploid_heterozygous_percentage} || '-'));

        my $bam_flag_file  = $subbuild->whole_rmdup_bam_flagstat_file;
        if(Genome::Sys->check_for_path_existence($bam_flag_file)) {
            my $flagstat_stats = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($bam_flag_file);

            for my $key (keys %$flagstat_stats) {
                $build_node->addChild( $doc->createAttribute('flagstat-' . $key, $flagstat_stats->{$key}) );
            }
        }
    }
    $self->_main_node->addChild($metrics_node);

    #File linecounts (== number of SNP calls, etc.)
    my $somatic_stats_node = $doc->createElement('somatic-stats');

    for my $subbuild (@all_subbuilds) {
        next unless $subbuild->type_name eq 'somatic';

        my $build_node = $somatic_stats_node->addChild( $doc->createElement('build') );

        $build_node->addChild( $doc->createAttribute('build-id', $subbuild->id) );
        $build_node->addChild( $doc->createAttribute('name', $subbuild->model->name) );

        my $somatic_data = $self->extract_somatic_report_data($subbuild);

        #TODO Add an SV section? Probably better to break down Breakdancer into SV types
        $build_node->addChild( $somatic_data );
    }
    $self->_main_node->addChild($somatic_stats_node);


    return 1;
}

sub extract_refalign_report_data {
    my $self = shift;
    my $build = shift;

    my $command = Genome::Model::Tools::Convergence::GatherRefalignSummaryData->create(
        build_id => $build->id,
    );

    unless($command->execute) {
        return;
    }

    my $data = $command->data;

    return unless $build->id eq $data->{build_id};

    return $data;
}

sub extract_somatic_report_data {
    my $self = shift;
    my $build = shift;

    my $report = $build->get_report('File_Summary');

    unless($report) {
        my $generate_command = Genome::Model::Somatic::Report::FileSummary->create(
            build_id => $build->id,
        );

        $report = $generate_command->generate_report;

        unless($report) {
            $self->error_message('Failed to generate file summary report.');
            return;
        }

        $build->add_report($report);
    }

    my @found_nodes = $report->xml->findnodes('//files');

    unless(scalar @found_nodes == 1) {
        $self->error_message('Found ' . (scalar @found_nodes) . ' matches instead of 1!');
        return;
    }

    return $found_nodes[0];
}

1;
