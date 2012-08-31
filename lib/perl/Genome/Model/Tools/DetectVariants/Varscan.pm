
package Genome::Model::Tools::DetectVariants::Varscan;

use strict;
use warnings;

use FileHandle;

use Genome;

class Genome::Model::Tools::DetectVariants::Varscan {
    is => 'Genome::Model::Tools::DetectVariants',
    has => [
        reference_sequence_input => {
            default => "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa",
        },
    ],
    has_optional => [
        snv_params => {
            default => '--min-var-freq 0.10 --p-value 0.10 --somatic-p-value 0.01',
        },
        indel_params => {
            default => '--min-var-freq 0.10 --p-value 0.10 --somatic-p-value 0.01',
        },
        detect_snvs => {
            default => '1',
        },
        detect_indels => {
            default => '1',
        },
    ],
   
    has_optional => [
        detect_snvs => {
            default => 1,
        },
        detect_indels => {
            default => 1,
        },
    ],

    has_param => [
        lsf_resource => {
            default => "-R 'select[model!=Opteron250 && type==LINUX64] span[hosts=1] rusage[mem=16000]' -M 1610612736",
        }
    ],
    has_constant_optional => [
        sv_params=>{},
    ],
};

sub help_brief {
    "Use Varscan for variant detection.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants var-scan --aligned_reads_input input.bam --reference_sequence_input reference.fa --output-directory ~/example/
EOS
}

sub help_detail {
    return <<EOS 
This tool runs Varscan for detection of SNPs and/or indels.
EOS
}

sub _detect_variants {
    my $self = shift;

    ## Get required parameters ##
    my $output_snp = $self->_snv_staging_output;
    my $output_snp_filtered = $self->_filtered_snv_staging_output;
    my $output_indel = $self->_indel_staging_output;
    my $output_indel_filtered = $self->_filtered_indel_staging_output;

    ## Get Varscan parameters ##
    my $snv_params = $self->snv_params || "";
    my $indel_params = $self->indel_params || "";
    my $result;
    if ( ($self->detect_snvs && $self->detect_indels) && ($snv_params eq $indel_params) ) {
        $result = $self->_run_varscan($output_snp, $output_snp_filtered, $output_indel, $output_indel_filtered, $snv_params);
    } else {
        # Run twice, since we have different parameters. Detect snps and throw away indels, then detect indels and throw away snps
        if ($self->detect_snvs && $self->detect_indels) {
            $self->status_message("Snp and indel params are different. Executing Varscan twice: once each for snps and indels with their respective parameters");
        }
        my ($temp_fh, $temp_name) = Genome::Sys->create_temp_file();
        my ($filtered_temp_fh, $filtered_temp_name) = Genome::Sys->create_filtered_temp_file();

        if ($self->detect_snvs) {
            $result = $self->_run_varscan($output_snp, $output_snp_filtered, $temp_name, $filtered_temp_name, $snv_params);
        }
        if ($self->detect_indels) {
            if($self->detect_snvs and not $result) {
                $self->status_message('Varscan did not report success for snv detection. Skipping indel detection.')
            } else {
                $result = $self->_run_varscan($temp_name, $filtered_temp_name, $output_indel, $output_indel_filtered, $indel_params);
            }
        }
    }

    return $result;
}

sub _run_varscan {
    my $self = shift;
    my ($output_snp, $output_snp_filtered, $output_indel, $output_indel_filtered, $varscan_params) = @_;

    my $reference = $self->reference_sequence_input;
    my $bam_file = $self->aligned_reads_input;

    my $varscan = Genome::Model::Tools::Varscan::Germline->create(
        version => $self->version,
        bam_file => $bam_file,
        reference => $reference,
        output_snp => $output_snp,
        output_snp_filtered => $output_snp_filtered,
        output_indel => $output_indel,
        output_indel_filtered => $output_indel_filtered,
        varscan_params => $varscan_params,
    );

    unless($varscan->execute()) {
        $self->error_message('Failed to execute Varscan: ' . $varscan->error_message);
        return;
    }

    return 1;
}

sub generate_metrics {
    my $self = shift;

    my $metrics = {};
    
    if($self->detect_snvs) {
        my $snp_count      = 0;
        
        my $snv_output = $self->_snv_staging_output;
        my $snv_fh = Genome::Sys->open_file_for_reading($snv_output);
        while (my $row = $snv_fh->getline) {
            $snp_count++;
        }
        $metrics->{'total_snp_count'} = $snp_count;
    }

    if($self->detect_indels) {
        my $indel_count    = 0;
        
        my $indel_output = $self->_indel_staging_output;
        my $indel_fh = Genome::Sys->open_file_for_reading($indel_output);
        while (my $row = $indel_fh->getline) {
            $indel_count++;
        }
        $metrics->{'total indel count'} = $indel_count;
    }

    return $metrics;
}

1;
