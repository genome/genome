package Genome::Model::ReferenceAlignment::Command::GenotypeMicroarrayConcordance;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceAlignment::Command::GenotypeMicroarrayConcordance {
    is => 'Command::V2',
    doc => 'Determine concordance between a genotype microarray VCF and a Reference Alignment build',
    has => [
        build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            doc => 'use this build to find genotype microarray builds and report concordance for all instrument data',
            shell_args_position => 1,
        },
        minimum_depth => {
            is => 'Integer',
            doc => 'The minimum depth required to consider concordant.',
            example_values => ['4'],
        },
        picard_version => {
            is => 'Text',
            doc => 'The version of picard to use.',
            example_values => ['1.123'],
        },
        dbsnp_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'The dbsnp build used to generate the microarray VCFs.',
            example_values => ['141765601'],
        },
    ],
    has_optional =>[
        _seqdict => {},
    ],
};

sub help_detail {
    return "Generate picard genotype concordance metrics between the genotype microarray VCF and a reference alignment build.";
}

sub execute {
    my $self = shift;

    my $ref_seq_build = $self->build->reference_sequence_build;
    my $seqdict = $ref_seq_build->get_sequence_dictionary('sam',$ref_seq_build->species_name,$self->picard_version);
    $self->_seqdict($seqdict);

    my @data;
    my @instrument_data = $self->build->instrument_data;
    unless (@instrument_data) {
        $self->error_message('Found no instrument data assigned to build: '. $self->build->display_name);
        die($self->error_message);
    }
    for my $instrument_data (@instrument_data) {
        my $qc_build = $instrument_data->lane_qc_build;
        unless ($qc_build) {
            $self->error_message('Failed to find lane qc build for instrument data : '. $instrument_data->display_name);
            die($self->error_message);
        }
        
        my $lane_qc_vcf = $self->resolve_lane_qc_vcf($qc_build);
        my ($microarray_vcf,$genotype_sample) = $self->resolve_genotype_microarray_vcf_and_sample($qc_build);
        
        my $intersect_vcf = Genome::Sys->create_temp_file_path($genotype_sample->name .'_x_'. $instrument_data->id .'.vcf');
        my $intersect_cmd = Genome::Model::Tools::BedTools::Intersect->create(
            input_file_a => $lane_qc_vcf,
            input_file_a_format => 'bed',
            input_file_b => $microarray_vcf,
            output_file => $intersect_vcf,
            header => 1,
        );
        unless ($intersect_cmd) {
            $self->error_message('Failed to create bedtools intersect!');
            die($self->error_message);
        }
        unless ($intersect_cmd->execute) {
            $self->error_message('Failed to execute bedtools intersect!');
            die($self->error_message);
        }

        # TODO: For exome limit to ROI from input model
        my $output = Genome::Sys->create_temp_file_path($instrument_data->id);
        my $gc_cmd = Genome::Model::Tools::Picard::GenotypeConcordance->create(
            truth_vcf => $microarray_vcf,
            call_vcf => $intersect_vcf,
            output => $output,
            truth_sample => $genotype_sample->name,
            call_sample => $qc_build->model->subject->name,
            min_dp => $self->minimum_depth,
            use_version => $self->picard_version,
        );
        unless ($gc_cmd) {
            $self->error_message('Failed to create GenotypeConcordance!');
            die($self->error_message);
        }
        unless ($gc_cmd->execute) {
            $self->error_message('Failed to execute GenotypeConcordance!');
            die($self->error_message);
        }

        my %data = (
            lane_qc_build_id => $qc_build->id,
            instrument_data_id => $instrument_data->id,
            sample_name => $qc_build->model->subject->name,
            microarray_sample_name => $genotype_sample->name,
            flow_cell_id => $instrument_data->flow_cell_id,
            lane => $instrument_data->lane,
            clusters => $instrument_data->clusters,
        );

        my $gc_summary_metrics_file = $output .'.genotype_concordance_summary_metrics';
        if (-e $gc_summary_metrics_file) {
            my $gc_summary_metrics_hash_ref = Genome::Model::Tools::Picard->parse_file_into_metrics_hashref($gc_summary_metrics_file);
            my $gc_summary_snp_metrics = $gc_summary_metrics_hash_ref->{'VARIANT_TYPE-SNP'};
            for my $key ($self->picard_metrics) {
                $data{$key} = $gc_summary_snp_metrics->{$key};
            }
        } else {
            die('Failed to find GenotypeConcordance summary file: '. $gc_summary_metrics_file);
        }
        push @data, \%data;
    } # for instrument data qc build

    my @headers = $self->summary_headers;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@headers,
        in_place_of_null_value => 'na',
    );
    for my $data (@data) {
        $writer->write_one($data);
    }
    $writer->output->close;
    return 1;
}

sub picard_metrics {
    my $class = shift;
    return qw/
                 HET_SENSITIVITY
                 HET_PPV
                 HOMVAR_SENSITIVITY
                 HOMVAR_PPV
                 VAR_SENSITIVITY
                 VAR_PPV
             /;
}

sub summary_headers {
    my $class = shift;
    my @headers = qw/
                        instrument_data_id
                        flow_cell_id
                        lane
                        clusters
                        sample_name
                        microarray_sample_name
                        lane_qc_build_id
                    /;
    push @headers, $class->picard_metrics;
    return @headers;
}

sub resolve_lane_qc_vcf {
    my $self = shift;
    my $qc_build = shift;
    
    my $build_vcf = $qc_build->get_detailed_snvs_vcf;
    unless (-e $build_vcf) {
        $self->error_message('Unable to find the snvs VCF for build : '. $qc_build->display_name);
        die($self->error_message);
    }
    
    # Sort the build VCF which adds the contig info to the VCF header (required by Picard GenotypeConcordance)
    # TODO : Add something here or to SortVcf that ensures the file is complete
    my $sorted_build_vcf = Genome::Sys->create_temp_file_path($qc_build->id .'_sorted.vcf');
    my $sort_cmd = Genome::Model::Tools::Picard::SortVcf->create(
        input_vcf => $build_vcf,
        output_vcf => $sorted_build_vcf,
        sequence_dictionary => $self->_seqdict,
        use_version => $self->picard_version,
    );
    unless ($sort_cmd) {
        $self->error_message('Failed to create sort VCF command!');
        die($self->error_message);
    }
    unless ($sort_cmd->execute) {
        $self->error_message('Failed to sort reference alignment VCF: '. $build_vcf);
        die($self->error_message);
    }
    return $sorted_build_vcf;
}

sub resolve_genotype_microarray_vcf_and_sample {
    my $self = shift;
    my $qc_build = shift;

    my $genotype_sample = $qc_build->model->subject;
    my $microarray_build = $qc_build->genotype_microarray_build;
    
    # get the VCF from an existing genotype microarray build if it exists already
    my $microarray_vcf;
    if ($microarray_build) {
        $microarray_vcf = $microarray_build->original_genotype_vcf_file_path;
        if (-e $microarray_vcf) {
            $genotype_sample = $microarray_build->model->subject;
        }
    }

    # there is no existing microarray build or the VCF does not exist (ie. old genotype build)
    unless (-e $microarray_vcf) {
        $self->debug_message('Get or create genotype VCF result for sample: '. $genotype_sample->display_name);
        
        my $vcf_result = Genome::InstrumentData::Microarray::Result::Vcf->get_or_create(
            sample => $genotype_sample,
            known_sites_build => $self->dbsnp_build,
        );
        $microarray_vcf = $vcf_result->vcf_path;
        unless (-e $microarray_vcf) {
            $self->error_message('Failed to get or create microarray VCF for sample: '. $genotype_sample->display_name);
            die($self->error_message);
        }
    }
    
    my $sorted_microarray_vcf = Genome::Sys->create_temp_file_path($genotype_sample->name .'_microarray_sorted.vcf');
    my $sort_cmd = Genome::Model::Tools::Picard::SortVcf->create(
        input_vcf => $microarray_vcf,
        output_vcf => $sorted_microarray_vcf,
        sequence_dictionary => $self->_seqdict,
        use_version => $self->picard_version,
    );
    unless ($sort_cmd) {
        $self->error_message('Failed to create sort VCF command!');
        die($self->error_message);
    }
    unless ($sort_cmd->execute) {
        $self->error_message('Failed to sort genotype microarray VCF: '. $microarray_vcf);
        die($self->error_message);
    }
    return ($sorted_microarray_vcf, $genotype_sample);
}

1;

