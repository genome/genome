package Genome::Model::ReferenceAlignment::Command::GenotypeConcordance;

use strict;
use warnings;

use Genome;

my %STATES_TO_IGNORE = (
    'MISSING' => 1,
    'NO_CALL' => 1,
    'VC_FILTERED' => 1,
    'LOW_DP' => 1,
    'LOW_GQ' => 1,
);

class Genome::Model::ReferenceAlignment::Command::GenotypeConcordance {
    is => 'Genome::Command::Base',
    doc => 'Determine concordance between a genotype microarray build and a Reference Alignment build',
    has => [
        models => {
            is => 'Genome::Model::ReferenceAlignment',
            is_many => 1,
            doc => 'use these models to find genotype microarray builds and report concordance',
            shell_args_position => 1,
        },
        output_file => {
            is => 'Text',
            doc => 'The basename of the output location.'
        },
        minimum_depth => {
            is => 'Integer',
            doc => 'The minimum depth required to consider concordant.',
            default_value => 4,
        }
    ],
};

sub help_detail {
    return "Generate tab-delimited (*.tsv) files with genotype concordance metrics between the genotype microarray build and a reference alignment build.";
}

sub execute {
    my $self = shift;
    my @headers = qw/model_name model_id build_id match total concordance/;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->output_file,
        separator => "\t",
        headers => \@headers,
    );
    for my $call_model ($self->models) {
        my %data;
        my @instrument_data = $call_model->instrument_data;
        unless (scalar(@instrument_data) eq 1) {
            die('Expecting Lane-QC model with a single piece of instrument data.  Found '. scalar(@instrument_data) .' instrument data associated with model '. $call_model->id);
        }

        my $ref_seq_build = $call_model->reference_sequence_build;
        my $seq_dict = $ref_seq_build->get_sequence_dictionary('sam',$ref_seq_build->species_name,'1.122');

        my $genotype_model = $call_model->genotype_microarray;
        unless ($genotype_model) {
            die('Failed to find a genotype microarray model for reference alignment model '. $call_model->id);
        }
        my $genotype_build = $genotype_model->last_succeeded_build;
        unless ($genotype_build) {
            die('Failed to find a succeeded build for genotype microarray model '. $genotype_model->id);
        }
        
        my $truth_vcf = $genotype_build->original_genotype_vcf_file_path;
        my $truth_sample = $genotype_model->subject->name;

        my $sorted_truth_vcf = Genome::Sys->create_temp_file_path('sorted_truth.vcf');
        unless (Genome::Model::Tools::Picard::SortVcf->execute(
            input_vcf => $truth_vcf,
            output_vcf => $sorted_truth_vcf,
            sequence_dictionary => $seq_dict,
        )) {
            die('Failed to sort genotype microarray VCF: '. $truth_vcf);
        }


        my $call_build = $call_model->last_succeeded_build;
        my $call_vcf = $call_build->data_directory .'/variants/snvs.detailed.vcf.gz';
        my $call_sample = $call_model->subject->name;

        my $sorted_call_vcf = Genome::Sys->create_temp_file_path('sorted_call.vcf');
        unless (Genome::Model::Tools::Picard::SortVcf->execute(
            input_vcf => $call_vcf,
            output_vcf => $sorted_call_vcf,
            sequence_dictionary => $seq_dict,
        )) {
            die('Failed to sort reference alignment VCF: '. $call_vcf);
        }
        my $model_output = Genome::Sys->create_temp_file_path($call_model->id);
        my $gc_cmd = Genome::Model::Tools::Picard::GenotypeConcordance->create(
            truth_vcf => $sorted_truth_vcf,
            call_vcf => $sorted_call_vcf,
            output => $model_output,
            truth_sample => $truth_sample,
            call_sample => $call_sample,
            min_dp => $self->minimum_depth,
        );
        unless ($gc_cmd) {
            die('Failed to create the command to run GenotypeConcordance!');
        }
        unless ($gc_cmd->execute) {
            die('Failed to execute GenotypeConcordance!');
        }
        my $detailed_metrics_file = $model_output .'.detailed_metrics.txt';
        my $hash_ref = Genome::Model::Tools::Picard::GenotypeConcordance->parse_file_into_metrics_hashref($detailed_metrics_file);
        my $snps = $hash_ref->{'SNP'}->{$truth_sample}->{$call_sample};
        for my $truth_state (keys %{$snps}) {
            if ($self->_ignore_state($truth_state)) { next; }
            for my $call_state (keys %{$snps->{$truth_state}} ) {
                if ($self->_ignore_state($call_state)) { next; }
                if ($truth_state eq $call_state) {
                    $data{match} += $snps->{$truth_state}->{$call_state};
                }
                $data{total} += $snps->{$truth_state}->{$call_state};
            }
        }
        $data{model_name} = $call_model->name;
        $data{model_id} = $call_model->id;
        $data{build_id} = $call_build->id;
        $data{concordance} = sprintf("%.04f",  (($data{match} / $data{total}) * 100) );
        $writer->write_one(\%data);
    }
    $writer->output->close;
}

sub _ignore_state {
    my $self = shift;
    my $state = shift;
    if ( defined($STATES_TO_IGNORE{$state}) ) {
        return $STATES_TO_IGNORE{$state};
    }
    return 0;
}



1;

