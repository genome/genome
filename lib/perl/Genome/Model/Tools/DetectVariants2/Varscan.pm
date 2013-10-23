package Genome::Model::Tools::DetectVariants2::Varscan;

use strict;
use warnings;

use FileHandle;

use Genome;

class Genome::Model::Tools::DetectVariants2::Varscan {
    is => ['Genome::Model::Tools::DetectVariants2::VarscanBase'],
    has => [
        params => {
            default => '--min-var-freq 0.10 --p-value 0.10 --somatic-p-value 0.01',
        },
    ],
    has_param => [
        lsf_resource => {
            default => "-R 'select[ncpus>=2] span[hosts=1] rusage[mem=16000]' -M 1610612736 -n 2",
        }
    ],
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants2 varscan --aligned_reads_input input.bam --reference_sequence_input reference.fa --output-directory ~/example/

gmt detect-variants2 varscan \
    --version 2.2.6
    --reference-build NCBI-human-build36 \ 
    --output-directory /gscuser/ssmith/od \ 
    --aligned-reads-input /gscmnt/gc7001/info/build_merged_alignments/merged-alignment-blade13-4-10.gsc.wustl.edu-rlong-14103-116553088/116553088.bam \
    --aligned-reads-sample Indel_Validation_H_IJ-NA19238-NA19238
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
    my $output_snp = $self->_temp_staging_directory."/snvs.hq";
    my $output_indel = $self->_temp_staging_directory."/indels.hq";

    unless ($self->version) {
        die $self->error_message("A version of Varscan must be specified");
    }

    # Grab the map_quality param and pass it separately
    my $params = $self->params;
    my $map_quality;
    if ($params =~ m/map-quality/) {
        ($map_quality) = ($params =~ m/--map-quality\s*(\d+)/);
        $params =~ s/--map-quality\s*(\d+)\s*//;
    }

    my ($samtools_params, $varscan_params) = $self->_split_params($params);
    my ($samtools_version, $use_baq, $other_params) = $self->_process_samtools_params($samtools_params);

    my %optional_samtools_params;
    $optional_samtools_params{samtools_version} = $samtools_version if $samtools_version;
    $optional_samtools_params{samtools_use_baq} = $use_baq if $use_baq;
    $optional_samtools_params{samtools_params} = $other_params if $other_params;

    my $varscan = Genome::Model::Tools::Varscan::Germline->create(
        bam_file => $self->aligned_reads_input,
        reference => $self->reference_sequence_input,
        output_snp => $output_snp,
        output_indel => $output_indel,
        varscan_params => $params,
        map_quality => $map_quality,
        no_headers => 1,
        version => $self->version,
        %optional_samtools_params,
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

sub parse_line_for_bed_intersection {
    my $class = shift;
    my $line = shift;

    unless ($line) {
        die $class->error_message("No line provided to parse_line_for_bed_intersection");
    }

    my ($chromosome, $position, $_reference, $consensus) = split "\t",  $line;

    if ($consensus =~ /\-|\+/) {
        return $class->_parse_indel_for_bed_intersection($line);
    } else {
        return $class->_parse_snv_for_bed_intersection($line);
    }
}

sub _parse_indel_for_bed_intersection {
    my $class = shift;
    my $line = shift;

    my ($chromosome, $position, $_reference, $consensus, @extra) = split "\t",  $line;
    
    my @variants;
    my @indels = Genome::Model::Tools::Bed::Convert::Indel::VarscanToBed->convert_indel($line);

    for my $indel (@indels) {
        my ($reference, $variant, $start, $stop) = @$indel;
        if (defined $chromosome && defined $position && defined $reference && defined $variant) {
            push @variants, [$chromosome, $stop, $reference, $variant];
        }
    }

    unless(@variants){
        die $class->error_message("Could not get chromosome, position, reference, or variant for line: $line");
    }

    return @variants;
}

sub _parse_snv_for_bed_intersection {
    my $class = shift;
    my $line = shift;

    my ($chromosome, $position, $reference, $consensus, @extra) = split("\t", $line);

    return [$chromosome, $position, $reference, $consensus];
}

1;
