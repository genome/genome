package Genome::Model::Event::Build::ReferenceAlignment::RunReports;

#REVIEW fdu 11/18/2009
#Missing ReferenceCoverage report for cDNA or RNA. Either ask jwalker
#to implement it or drop those codes.

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::ReferenceAlignment::RunReports {
    is => [ 'Genome::Model::Event' ],
};

my %REPORT_TYPES = (
    Mapcheck           => 'Mapcheck',
    DbSnpConcordance   => 'dbSNP_Concordance',
    GoldSnpConcordance => 'Gold_SNP_Concordance',
    InputBaseCounts    => 'Input_Base_Counts'
);


sub execute {
    my $self  = shift;
    my $model = $self->model;
    my $build = $self->build;

    if ($model->read_aligner_name =~ /^Imported$/i) {
        $self->debug_message("This build uses Imported as alinger so skip InputBaseCounts");
        delete $REPORT_TYPES{InputBaseCounts};
    }

    unless (defined $model->dbsnp_build) {
        $self->debug_message("No dbsnp_build defined for model, skipping dbsnp concordance report.");
        delete $REPORT_TYPES{DbSnpConcordance};
    } 

    unless ($self->validate_gold_snp_path) {
        $self->debug_message("No valid gold_snp_path for the build, skip GoldSnpConcordance report");
        delete $REPORT_TYPES{GoldSnpConcordance};
    }
    
    unless ( ($model->dna_type eq 'cdna' || $model->dna_type eq 'rna') && $model->reference_sequence_name =~ /^XStrans_adapt_smallRNA_ribo/ ) {
        delete $REPORT_TYPES{ReferenceCoverage};
    }
    
    unless ($self->create_directory($build->resolve_reports_directory) ) {
	    die('Could not create reports directory at: '. $build->resolve_reports_directory);
    }

    for my $report_type (keys %REPORT_TYPES) {
        my $report_class = 'Genome::Model::ReferenceAlignment::Report::' . $report_type;
        my $report_name = 'unknown';
        $self->debug_message("Starting $report_type report.");

        my $report_def = $report_class->create(build_id => $build->id);
        unless ($report_def) {
            $self->error_message("Error creating report $report_class!: " . $report_class->error_message());
            die($self->error_message);
        }
        $report_name = $report_def->name;
        $self->debug_message("Defined report with name $report_name");
        
        my $report = $report_def->generate_report;
        unless ($report) {
            $self->error_message("Error generating $report_name ($report_class)!: " . $report_class->error_message());
            $report_def->delete;
            die($self->error_message);
        } 
        else {
            $self->debug_message("Successfully generated report: $report_name");
        }
        $self->debug_message("About to add report: $report_name to build: ".$self->build->id);
        
        if ($build->add_report($report)) {
            $self->debug_message('Saved report: '.$report);
        } 
        else {
            $self->error_message('Error saving '.$report.'. Error: '. $self->build->error_message);
            die($self->error_message);
        }
    }

    ##############################################
    #Summary Report

    $self->debug_message('Starting report summary.');
    my $r = Genome::Model::ReferenceAlignment::Report::Summary->create( build_id => $build->id );

    my @templates = $r->report_templates;
    $self->debug_message("Using report templates: ". join(",",@templates));

    my $summary_report = $r->generate_report;

    unless ($build->add_report($summary_report)) {
        $self->error_message("Failed to save report: ". $summary_report->name .' to '. $build->resolve_reports_directory);
        return;
    }
    $self->debug_message('Report summary complete.');

    return $self->verify_successful_completion;
}


sub validate_gold_snp_path {
    my $self = shift;

    my $gold_snp_path = $self->build->gold_snp_path;
    unless ($gold_snp_path and -s $gold_snp_path) {
        $self->debug_message('No gold_snp_path provided for the build or it is empty');
        return;
    }

    my $head    = `head -1 $gold_snp_path`;
    my @columns = split /\s+/, $head;
    
    unless (@columns and @columns == 9) {
        $self->debug_message("Gold snp file: $gold_snp_path is not 9-column format");
        return;
    }
    return 1;
}
    

sub verify_successful_completion {
    my $self = shift;
    my $model = $self->model;
    my $build = $self->build;
    unless ($build) {
        $self->error_message('Failed verify_successful_completion of RunReports step. Build is undefined.');
        return;
    }

    my $report_dir = $build->resolve_reports_directory;
    
    my $gold_snp_path = $self->build->gold_snp_path;
    delete $REPORT_TYPES{GoldSnpConcordance} unless $gold_snp_path and -s $gold_snp_path;

    delete $REPORT_TYPES{InputBaseCounts} if $model->read_aligner_name =~ /^Imported$/i;

    unless ( ($model->dna_type eq 'cdna' || $model->dna_type eq 'rna') && $model->reference_sequence_name =~ /^XStrans_adapt_smallRNA_ribo/ ) {
        delete $REPORT_TYPES{ReferenceCoverage};
    }
        
    my @sub_dirs = (values %REPORT_TYPES, 'Summary');  
   
    for my $sub_directory (@sub_dirs) {
        unless (-d $report_dir .'/'. $sub_directory) {
            $self->error_message('Failed verify_successful_completeion of RunReports step.  Failed to find directory: '. $report_dir .'/'. $sub_directory);
            return;
        }
    }
    return 1;
}

1;

#$HeadURL$
#$Id$
