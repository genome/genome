package Genome::Model::Tools::DetectVariants::Samtools;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::Model::Tools::DetectVariants::Samtools {
    is => ['Genome::Model::Tools::DetectVariants'],
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
            default => "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>1000 && mem>16000] span[hosts=1] rusage[tmp=1000:mem=16000]' -M 1610612736",
        }
    ],
    has_constant_optional => [
        sv_params=>{},
        detect_svs=>{},
    ],
    has => [
        _genotype_detail_base_name => {
            is => 'Text',
            default_value => 'report_input_all_sequences',
            is_input => 1,
        },
        genotype_detail_output => {
            calculate_from => ['_genotype_detail_base_name', 'output_directory'],
            calculate => q{ join("/", $output_directory, $_genotype_detail_base_name); },
            is_output => 1,
        },
        _genotype_detail_staging_output => {
            calculate_from => ['_temp_staging_directory', '_genotype_detail_base_name'],
            calculate => q{ join("/", $_temp_staging_directory, $_genotype_detail_base_name); },
        },
    ]
};

sub help_brief {
    "Use samtools for variant detection.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants samtools --version r453 --aligned_reads_input input.bam --reference_sequence_input reference.fa --working-directory ~/example/
EOS
}

sub help_detail {
    return <<EOS 
This tool runs samtools for detection of SNVs and/or indels.
EOS
}

sub _detect_variants {
    my $self = shift;

    my $snv_params = $self->snv_params || "";
    my $indel_params = $self->indel_params || "";
    my $result;
    if ( ($self->detect_snvs && $self->detect_indels) && ($snv_params eq $indel_params) ) {
        $result = $self->_run_samtools($self->_snv_staging_output, $self->_filtered_snv_staging_output, $self->_indel_staging_output, $self->_filtered_indel_staging_output, $snv_params);
    } else {
        # Run twice, since we have different parameters. Detect snvs and throw away indels, then detect indels and throw away snvs
        if ($self->detect_snvs && $self->detect_indels) {
            $self->status_message("Snv and indel params are different. Executing Samtools twice: once each for snvs and indels with their respective parameters");
        }
        my ($temp_fh, $temp_name) = Genome::Sys->create_temp_file();
        my ($filtered_temp_fh, $filtered_temp_name) = Genome::Sys->create_temp_file();

        if ($self->detect_snvs) {
            $result = $self->_run_samtools($self->_snv_staging_output, $self->_filtered_snv_staging_output, $temp_name, $filtered_temp_name, $snv_params);
        }
        if ($self->detect_indels) {
            if($self->detect_snvs and not $result) {
                $self->status_message('Samtools did not report success for snv detection. Skipping indel detection.')
            } else {
                $result = $self->_run_samtools($temp_name, $filtered_temp_name, $self->_indel_staging_output, $self->_filtered_indel_staging_output, $indel_params);
            }
        }
    }

    return $result;
}

sub _run_samtools {
    my $self = shift;
    my ($snv_output_file, $filtered_snv_output_file, $indel_output_file, $filtered_indel_file, $parameters) = @_;
    
    my $ref_seq_file = $self->reference_sequence_input;
    my $bam_file = $self->aligned_reads_input;
    my $sam_pathname = Genome::Model::Tools::Sam->path_for_samtools_version($self->version);

    # Remove the result files from any previous run
    unlink($snv_output_file, $filtered_snv_output_file, $indel_output_file, $filtered_indel_file);

    
    #two %s are switch to indicate snvs or indels and output file name

##############################################################
#find a way to detect the version

   my $samtools_cmd = "$sam_pathname pileup -c $parameters -f $ref_seq_file %s $bam_file > %s";
    #Originally "-S" was used as SNP calling. In r320wu1 version, "-v" is used to replace "-S" but with 
    #double indel lines embedded, this need sanitized
    #$rv = system "$samtools_cmd -S $bam_file > $snv_output_file"; 
    
    
    my $snv_cmd = sprintf($samtools_cmd, '-v', $snv_output_file);
    
    my $rv = Genome::Sys->shellcmd(
        cmd => $snv_cmd,
        input_files => [$bam_file, $ref_seq_file],
        output_files => [$snv_output_file],
        allow_zero_size_output_files => 1,
    );
    unless($rv) {
        $self->error_message("Running samtools SNP failed.\nCommand: $snv_cmd");
        return;
    }

    if (-e $snv_output_file and not -s $snv_output_file) {
        $self->warning_message("No SNVs detected.");
    }
    else {
        my $snp_sanitizer = Genome::Model::Tools::Sam::SnpSanitizer->create(snp_file => $snv_output_file);
        $rv = $snp_sanitizer->execute;
        unless($rv and $rv == 1) {
            $self->error_message("Running samtools snp-sanitizer failed with exit code $rv");
            return;
        }
    }




    my $indel_cmd = sprintf($samtools_cmd, '-i', $indel_output_file);
    $rv = Genome::Sys->shellcmd(
        cmd => $indel_cmd,
        input_files => [$bam_file, $ref_seq_file],
        output_files => [$indel_output_file],
        allow_zero_size_output_files => 1,
    );
    unless($rv) {
        $self->error_message("Running samtools indel failed.\nCommand: $indel_cmd");
        return;
    }

    if (-e $indel_output_file and not -s $indel_output_file) {
        $self->warning_message("No indels detected.");
    }

    #for capture models we need to limit the snvs and indels to within the defined target regions
    if ($self->capture_set_input) {
        my $bed_file = $self->capture_set_input;
        for my $var_file ($snv_output_file, $indel_output_file) {
            unless (-s $var_file) {
                $self->warning_message("Skip limiting $var_file to target regions because it is empty.");
                next;
            }
            my $tmp_limited_file = $var_file .'_limited';
            my $no_limit_file = $var_file .'.no_bed_limit';
            unless (Genome::Model::Tools::Sam::LimitVariants->execute(
                variants_file => $var_file,
                bed_file => $bed_file,
                output_file => $tmp_limited_file,
            )) {
                $self->error_message('Failed to limit samtools variants '. $var_file .' to within capture target regions '. $bed_file);
                die($self->error_message);
            }
            unless (move($var_file,$no_limit_file)) {
                $self->error_message('Failed to move all variants from '. $var_file .' to '. $no_limit_file);
                die($self->error_message);
            }
            unless (move($tmp_limited_file,$var_file)) {
                $self->error_message('Failed to move limited variants from '. $tmp_limited_file .' to '. $var_file);
                die($self->error_message);
            }
        }
    }


    if (-s $indel_output_file) {
        my %indel_filter_params = ( indel_file => $indel_output_file, out_file => $filtered_indel_file );
        # for capture data we do not know the proper ceiling for depth
        if ($self->capture_set_input) {
            $indel_filter_params{max_read_depth} = 1000000;
        }
        my $indel_filter = Genome::Model::Tools::Sam::IndelFilter->create(%indel_filter_params);
        unless($indel_filter->execute) {
            $self->error_message("Running sam indel-filter failed.");
            return;
        }
    }
    else {
        Genome::Sys->write_file($filtered_indel_file);
    }

    if (-s $snv_output_file) {
        my $snp_filter = Genome::Model::Tools::Sam::SnpFilter->create(
            snp_file   => $snv_output_file,
            out_file   => $filtered_snv_output_file,
            indel_file => $filtered_indel_file,
        );
        unless($snp_filter->execute) {
            $self->error_message("Running sam snp-filter failed.");
            return;
        }
    }
    else {
        Genome::Sys->write_file($filtered_snv_output_file);
    }

    $rv = $self->generate_genotype_detail_file($snv_output_file);
    unless($rv) {
        $self->error_message('Generating genotype detail file errored out');
        die($self->error_message);
    }
    
    return $self->verify_successful_completion($snv_output_file, $filtered_snv_output_file, $indel_output_file);
}




sub verify_successful_completion {
    my $self = shift;

    for my $file (@_) {
        unless (-e $file) {
            $self->error_message("$file was not successfully created. Failure in verify_successful_completion.");
            return;
        }
    }
    
    return 1;
}

sub generate_genotype_detail_file {
    my $self = shift; 
    my $snv_output_file = shift;
    
    unless(-f $snv_output_file) { # and -s $snv_output_file) {
        $self->error_message("SNV output File: $snv_output_file is invalid.");
        die($self->error_message);
    }

    if (not -s $snv_output_file) {
        $self->warning_message("No report input file generated for SNVs because no SNVs were detected.");
        return 1;
    }
    
    my $report_input_file = $self->_genotype_detail_staging_output;

    my $snp_gd = Genome::Model::Tools::Snp::GenotypeDetail->create(
        snp_file   => $snv_output_file,
        out_file   => $report_input_file,
        snp_format => 'sam',
    );
    
    return $snp_gd->execute;
}

sub generate_metrics {
    my $self = shift;

    my $metrics = {};
    
    if($self->detect_snvs) {
        my $snp_count      = 0;
        my $snp_count_good = 0;
        
        my $snv_output = $self->_snv_staging_output;
        my $snv_fh = Genome::Sys->open_file_for_reading($snv_output);
        while (my $row = $snv_fh->getline) {
            $snp_count++;
            my @columns = split /\s+/, $row;
            $snp_count_good++ if $columns[4] >= 15 and $columns[7] > 2;
        }
        $metrics->{'total_snp_count'} = $snp_count;
        $metrics->{'confident_snp_count'} = $snp_count_good;
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
