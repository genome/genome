package Genome::Model::Tools::DetectVariants::Maq;

use strict;
use warnings;

use Genome;

use File::Path;
use File::Temp;
use IO::File;

class Genome::Model::Tools::DetectVariants::Maq {
    is => ['Genome::Model::Tools::DetectVariants'],
    has => [
        _pileup_base_name => {
            is => 'Text',
            default_value => 'pileup_all_sequences',
            is_input => 1,
        },
        _pileup_staging_output => {
            calculate_from => ['_temp_staging_directory', '_pileup_base_name'],
            calculate      => q{ join('/', $_temp_staging_directory, $_pileup_base_name); },
        },
        pileup_output => {
            calculate_from => ['output_directory', '_pileup_base_name'],
            calculate      => q{ join('/', $output_directory, $_pileup_base_name); },
            is_output => 1,
        },
        _consensus_directory_base_name => {
            is => 'Text',
            default_value => 'consensus',
            is_input => 1,
        },
        _consensus_staging_directory => {
            calculate_from => ['_temp_staging_directory', '_consensus_directory_base_name'],
            calculate      => q{ join('/', $_temp_staging_directory, $_consensus_directory_base_name); },
        },
        consensus_directory => {
            calculate_from => ['output_directory', '_consensus_directory_base_name'],
            calculate      => q{ join('/', $output_directory, $_consensus_directory_base_name); },
            is_output => 1,
        },
        _consensus_output_base_name => {
            is => 'Text',
            default_value => 'all_sequences.cns',
            is_input => 1, 
        },
        _consensus_staging_output => {
            calculate_from => ['_consensus_staging_directory', '_consensus_output_base_name'],
            calculate      => q{ join('/', $_consensus_staging_directory, $_consensus_output_base_name); },
        },
        consensus_output => {
            calculate_from => ['consensus_directory', '_consensus_output_base_name'],
            calculate      => q{ join('/', $consensus_directory, $_consensus_output_base_name); },
            is_output => 1,
        },
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
            calculate_from => ['output_directory', '_genotype_detail_base_name'],
            calculate => q{ join("/", $output_directory, $_genotype_detail_base_name); },
        },
        _indelpe_base_name => {
            is => 'Text',
            default_value => 'indelpe.out',
            is_input => 1,
        },
        _indelpe_staging_output => {
            calculate_from => ['_temp_staging_directory', '_indelpe_base_name'],
            calculate      => q{ join('/', $_temp_staging_directory, $_indelpe_base_name); },
        },
        indelpe_output => {
            calculate_from => ['output_directory', '_indelpe_base_name'],
            calculate      => q{ join('/', $output_directory, $_indelpe_base_name); },
            is_output => 1,
        },
        _sorted_indelpe_base_name => {
            is => 'Text',
            default_value => 'indelpe.sorted.out',
            is_input => 1,
        },
        _sorted_indelpe_staging_output => {
            calculate_from => ['_temp_staging_directory', '_sorted_indelpe_base_name'],
            calculate      => q{ join('/', $_temp_staging_directory, $_sorted_indelpe_base_name); },
        },
        sorted_indelpe_output => {
            calculate_from => ['output_directory', '_sorted_indelpe_base_name'],
            calculate      => q{ join('/', $output_directory, $_sorted_indelpe_base_name); },
            is_output => 1,
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
    has_constant_optional => [
        sv_params => { },
        detect_svs => { },
    ],
    has_param => [
        lsf_resource => {
            default => "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>1000 && mem>16000] span[hosts=1] rusage[tmp=1000:mem=16000]' -M 1610612736",
        }
    ],
};

sub help_brief {
    "Use maq for variant detection.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants maq --version 0.7.1 --aligned_reads_input input.map --reference_sequence_input reference.bfa --working-directory ~/example/
EOS
}

sub help_detail {
    return <<EOS 
This tool runs maq for detection of SNVs and/or indels.
EOS
}

sub _detect_variants {
    my $self = shift;

    my $snv_params = $self->snv_params || "";
    my $indel_params = $self->indel_params || "";
    my $genotyper_params;

    # make sure the params are the same, or we are only detecting one type of variant
    if ( ($self->detect_snvs && $self->detect_indels) && ($snv_params ne $indel_params) ) {
        $self->status_message("Snv and indel params are different. This is not supported, as these parameters only affect the update genotype step");
        die;
    }
    if ($self->detect_snvs) {
        $genotyper_params = $snv_params;
    } else {
        $genotyper_params = $indel_params;
    }

    my $snv_output = $self->_snv_staging_output;
    my $filtered_snv_output = $self->_filtered_snv_staging_output;
    my $indel_output = $self->_indel_staging_output;

    my $result;
    if ($self->detect_snvs && $self->detect_indels) {
        $result = $self->_run_maq($snv_output, $filtered_snv_output, $indel_output, $genotyper_params);
    } else {
        # Run just snvs or indels if we dont want both. Throw away the other type of variant
        my ($temp_fh, $temp_name) = Genome::Sys->create_temp_file();
        my ($filtered_temp_fh, $filtered_temp_name) = Genome::Sys->create_temp_file();

        if ($self->detect_snvs) {
            $result = $self->_run_maq($snv_output, $filtered_snv_output, $temp_name, $genotyper_params);
        }
        if ($self->detect_indels) {
            $result = $self->_run_maq($temp_name, $filtered_temp_name, $indel_output, $genotyper_params);
        }
    }

    return $result;
}

sub _run_maq {
    my ($self, $snv_output, $filtered_snv_output, $indel_output, $genotyper_params) = @_;

    my $maq_pathname    = Genome::Model::Tools::Maq->path_for_maq_version($self->version);
    my $maq_pl_pathname = Genome::Model::Tools::Maq->proper_maq_pl_pathname($self->version);

    # ensure the reference sequence exists.
    my $reference_sequence = $self->reference_sequence_input;

    $self->update_genotype($genotyper_params);

    my $assembly_output = $self->_consensus_staging_output;
    unless ( Genome::Sys->check_for_path_existence($assembly_output) ) {
        $self->error_message("Assembly output file $assembly_output does not exist");
        return;
    }

    my $pileup_output = $self->pileup_output;
    my $accumulated_alignments = $self->aligned_reads_input;
    my $indelpe_output           = $self->_indelpe_staging_output;
    my $sorted_indelpe_output    = $self->_sorted_indelpe_staging_output;
    
    # Remove the result files from any previous run
    unlink ($snv_output, $filtered_snv_output, $indel_output, $pileup_output, $indelpe_output, $sorted_indelpe_output);

    my $cmd = "$maq_pathname cns2snp $assembly_output > $snv_output";
    unless (Genome::Sys->shellcmd(cmd => $cmd, input_files => [$assembly_output]) ) {
        $self->error_message("cns2snp.\ncmd: $cmd");
        return;
    }

    $cmd = "$maq_pathname indelsoa $reference_sequence $accumulated_alignments > $indel_output";
    unless (Genome::Sys->shellcmd(cmd => $cmd, input_files => [$reference_sequence, $accumulated_alignments]) ) {
        $self->error_message("indelsoa.\ncmd: $cmd");
        return;
    }

    my $filter = 'perl -nae '."'".'print if $F[2] =~ /^(\*|\+)$/'."'";
    $cmd = "$maq_pathname indelpe $reference_sequence $accumulated_alignments | $filter > $indelpe_output";
    unless (Genome::Sys->shellcmd(cmd => $cmd, input_files => [$reference_sequence, $accumulated_alignments]) ) {
        $self->error_message("indelpe.\ncmd: $cmd");
        return;
    }

    my $rv = Genome::Model::Tools::Snp::Sort->execute(
        snp_file    => $indelpe_output,
        output_file => $sorted_indelpe_output,
    );
    unless ($rv) {
        $self->error_message('Failed to run gmt snp sort');
        return;
    }
        
    my $indel_param;
    if (-s $sorted_indelpe_output) {
        $indel_param = "-F '$sorted_indelpe_output'";
    }
    else {
        $self->warning_message('Omitting indelpe data from the SNPfilter results because no indels were found');
        $indel_param = '';
    }

    $cmd = "$maq_pl_pathname SNPfilter $indel_param $snv_output > $filtered_snv_output";
    unless (Genome::Sys->shellcmd(cmd => $cmd, input_files => [$snv_output]) ) {
        $self->error_message("SNPfilter.\ncmd: $cmd");
        return;
    }
    
    # Running pileup requires some parsing of the snv file
    my ($tmp_fh, $temp_filename) = Genome::Sys->create_temp_file;
    my $snv_fh = Genome::Sys->open_file_for_reading($snv_output);
    unless ($snv_fh) {
        $self->error_message("Can't open snv output file for reading: $!");
        return;
    }
    while(<$snv_fh>) {
        chomp;
        my ($id, $start, $ref_sequence, $iub_sequence, $quality_score,
            $depth, $avg_hits, $high_quality, $unknown) = split("\t");
        $tmp_fh->print("$id\t$start\n");
    }
    $tmp_fh->close();
    $snv_fh->close();

    $cmd = sprintf(
        "$maq_pathname pileup -v -l %s %s %s > %s",
        $temp_filename,
        $reference_sequence,
        $accumulated_alignments,
        $pileup_output
    );
    unless (Genome::Sys->shellcmd(cmd => $cmd, input_files => [$temp_filename, $reference_sequence, $accumulated_alignments]) ) {
        $self->error_message("pileup.\ncmd: $cmd");
        return;
    }

    unless ($self->generate_genotype_detail_file($snv_output, $pileup_output)) {
        $self->error_message('Generating genotype detail file errored out');
        return;
    }

    return $self->verify_successful_completion($snv_output, $pileup_output, $filtered_snv_output, $indel_output);
}

#This is being overwritten just so the indelpe output can be used instead of the indelsoa output
#(There may be more elegant approaches than this...)
sub _generate_standard_files {
    my $self = shift;
    
    my $detector = 'Maq';
    my $module_base = 'Genome::Model::Tools::Bed::Convert';
    
    my $retval = 1;
    
    if($self->detect_snvs) {
        my $snv_module = join('::', $module_base, 'Snv', $detector . 'ToBed'); 
        
        for my $variant_file ($self->_snv_staging_output, $self->_filtered_snv_staging_output) {
            if(Genome::Sys->check_for_path_existence($variant_file)) {
                $retval &&= $self->_run_converter($snv_module, $variant_file);
            }  
        }
    }
    
    if($self->detect_indels) {
        my $snv_module = join('::', $module_base, 'Indel', $detector . 'ToBed'); 
        
        for my $variant_file ($self->_indelpe_staging_output, $self->_sorted_indelpe_staging_output) {
            if(Genome::Sys->check_for_path_existence($variant_file)) {
                $retval &&= $self->_run_converter($snv_module, $variant_file);
            }  
        }
    }
    
    return $retval;
}


sub verify_successful_completion {
    my $self = shift;

    for my $file (@_) {
        unless (-e $file) {
            $self->error_message("File $file doesn't exist or has no data");
            return;
        }
    }

    return 1;
}

sub generate_genotype_detail_file {
    my $self  = shift;
    my ($snv_output, $pileup_output) = @_;
    
    my $report_input_file = $self->_genotype_detail_staging_output;

    for my $file ($snv_output, $pileup_output) {
        unless (-s $file) {
            $self->error_message("File $file dosen't exist or has no data");
            return;
        }
    }

    unlink $report_input_file if -e $report_input_file;
    
    my $snp_gd = Genome::Model::Tools::Snp::GenotypeDetail->create(
        snp_file   => $snv_output,
        out_file   => $report_input_file,
        snp_format => 'maq',
        maq_pileup_file => $pileup_output,
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
            my ($r,$p,$a1,$a2,$q,$c) = split /\s+/, $row;
            $snp_count_good++ if $q >= 15 and $c > 2;
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

sub update_genotype {
    my $self = shift;
    my $genotyper_params = shift;

    $DB::single = $DB::stopper;

    my $maq_pathname = Genome::Model::Tools::Maq->path_for_maq_version($self->version);
    my $consensus_dir = $self->_consensus_staging_directory;
    unless (-d $consensus_dir) {
        unless (Genome::Sys->create_directory($consensus_dir)) {
            $self->error_message("Failed to create consensus directory $consensus_dir:  $!");
            return;
        }
    }

    my $consensus_file = $self->_consensus_staging_output;
    my $ref_seq_file = $self->reference_sequence_input;
    my $accumulated_alignments_file = $self->aligned_reads_input;

    my $cmd = $maq_pathname .' assemble '. $genotyper_params.' '. $consensus_file .' '. $ref_seq_file .' '. $accumulated_alignments_file;
    $self->status_message("\n************* UpdateGenotype cmd: $cmd *************************\n\n");
    Genome::Sys->shellcmd(
                    cmd => $cmd,
                    input_files => [$ref_seq_file,$accumulated_alignments_file],
                    output_files => [$consensus_file],
                );

    return $self->update_genotype_verify_successful_completion;
}

sub update_genotype_verify_successful_completion {
    my $self = shift;

    my $consensus_file = $self->_consensus_staging_output;
    unless (-e $consensus_file && -s $consensus_file > 20) {
        $self->error_message("Consensus file $consensus_file is too small");
        return;
    }
    return 1;
}

1;

