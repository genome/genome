package Genome::Model::Tools::DetectVariants2::Samtools;

use strict ;
use warnings;
use Genome;

class Genome::Model::Tools::DetectVariants2::Samtools {
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    has_param => [
        lsf_resource => {
            default => "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>1000 && mem>16000] span[hosts=1] rusage[tmp=1000:mem=16000]' -M 1610612736",
        }
    ],
    has_optional => [
        _genotype_detail_base_name => {
            is => 'Text',
            default_value => 'report_input_all_sequences',
            is_input => 1,
        },
        _genotype_detail_staging_output => {
            calculate_from => ['_temp_staging_directory', '_genotype_detail_base_name'],
            calculate => q{ join("/", $_temp_staging_directory, $_genotype_detail_base_name); },
        },
        _vcf_staging_output => {
            calculate_from => ['_temp_staging_directory'],
            calculate => q{ join("/", $_temp_staging_directory, 'vars.vcf'); },
        },
        _vcf_sani_staging_output => {
            calculate_from => ['_temp_staging_directory'],
            calculate => q { join("/", $_temp_staging_directory, 'vars.vcf.sanitized'); },
        },
        _header_scratch_output => {
            calculate_from => ['_temp_scratch_directory'],
            calculate => q { join("/", $_temp_scratch_directory, 'vars.vcf.header'); },
        },
        _snv_scratch_output => {
            calculate_from => ['_temp_scratch_directory'],
            calculate => q { join("/", $_temp_scratch_directory, 'snvs.vcf'); },
        },
        _indel_scratch_output => {
            calculate_from => ['_temp_scratch_directory'],
            calculate => q { join("/", $_temp_scratch_directory, 'indels.vcf'); },
        },
    ]
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detec,t-variants2 samtools --version r963 --aligned_reads_input input.bam --reference_sequence_input reference.fa --working-directory ~/example/
EOS
}

sub help_detail {
    return <<EOS 
This tool runs samtools for detection of SNVs and/or indels.
EOS
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    my @versions = Genome::Model::Tools::Sam->available_samtools_versions;
    for my $v (@versions){
        if($v eq $version){
            return 1;
        }
    }
    return 0;
}

sub _detect_variants {
    my $self = shift;

    my $ref_seq_file = $self->reference_sequence_input;
    my $bam_file     = $self->aligned_reads_input;
    my $snv_output_file     = $self->_snv_staging_output;
    my $indel_output_file   = $self->_indel_staging_output;
    my $vcf_output_file     = $self->_vcf_staging_output;
    my $filtered_indel_file = $self->_filtered_indel_staging_output;

    my $bed_file = $self->region_of_interest->merged_bed_file if $self->region_of_interest;
    my ($is_mpileup, $parameters) = $self->_mpileup_or_pileup;

    #two %s are switch to indicate snvs or indels and output file
    #name, snv and indel share the same parameter
    my $samtools_path   = Genome::Model::Tools::Sam->path_for_samtools_version($self->version);
    my $samtools_pu_cmd = "$samtools_path pileup -c $parameters -f $ref_seq_file %s $bam_file > %s";
    my ($samtools_cmd, $snv_cmd, $check_out);

    if ($is_mpileup) { #mpileup params contain the word mpileup
        my $bcftools_path = Genome::Model::Tools::Sam->path_for_bcftools($self->version);
        $parameters .= " -l $bed_file" if $self->region_of_interest;
        $snv_cmd   = "$samtools_path $parameters -f $ref_seq_file $bam_file | $bcftools_path view -Avcg - > $vcf_output_file";
        $check_out = $vcf_output_file;
    } 
    else {
        $snv_cmd   = sprintf($samtools_pu_cmd, '-v', $snv_output_file);
        $check_out = $snv_output_file;
    }

    my $rv = Genome::Sys->shellcmd(
        cmd          => $snv_cmd,
        input_files  => [$bam_file, $ref_seq_file],
        output_files => [$check_out],
        set_pipefail => 1,
        allow_zero_size_output_files => 1,
    );
    unless ($rv) {
        $self->error_message("Running samtools SNP failed.\nCommand: $snv_cmd");
        return;
    }

    #steps needed with pileup that mpileup doesn't use -i only shows lines/consensus with indels
    if ($is_mpileup) {
        $self->warning_message("No vars detected.") if -z $vcf_output_file;
        #want to make sure  the headers from the VCF form aren't used because if they are
        #the test to compare the line counts will fail everytime
        #we also want to just look at the indels in the indels.hq.v2.bed file
        unless ($self->create_snv_indel_output_file) {
            $self->error_message("Failed to create snv and indel files from vcf file: $vcf_output_file");
            return;
        }
    }
    else { #pileup
        $self->warning_message("No snvs detected.") if -z $snv_output_file;
        my $snp_sanitizer = Genome::Model::Tools::Sam::SnpSanitizer->create(snp_file => $snv_output_file);
        $rv = $snp_sanitizer->execute;
        unless ($rv and $rv == 1) {
            $self->error_message("Running samtools snp-sanitizer failed with exit code $rv");
            return;
        }

        my $indel_cmd = sprintf($samtools_pu_cmd, '-i', $indel_output_file);
        $rv = Genome::Sys->shellcmd(
            cmd          => $indel_cmd,
            input_files  => [$bam_file, $ref_seq_file],
            output_files => [$indel_output_file],
            allow_zero_size_output_files => 1,
        );
        unless($rv) {
            $self->error_message("Running samtools indel failed.\nCommand: $indel_cmd");
            return;
        }
        $self->warning_message("No indels detected.") if -z $indel_output_file;
    
        #for capture models we need to limit the snvs and indels to within the defined target regions
        if ($self->region_of_interest) {
            for my $var_file ($snv_output_file, $indel_output_file) {
                unless (-s $var_file) {
                    $self->warning_message("Skip limiting $var_file to target regions because it is empty.");
                    next;
                }
                my $tmp_limited_file = $var_file .'_limited';
                my $no_limit_file    = $var_file .'.no_bed_limit';
                unless (Genome::Model::Tools::Sam::LimitVariants->execute(
                    variants_file => $var_file,
                    bed_file      => $bed_file,
                    output_file   => $tmp_limited_file,
                )) {
                    $self->error_message('Failed to limit samtools variants '. $var_file .' to within capture target regions '. $bed_file);
                    die $self->error_message;
                }
                unless (move($var_file, $no_limit_file)) {
                    $self->error_message('Failed to move all variants from '. $var_file .' to '. $no_limit_file);
                    die $self->error_message;
                }
                unless (move($tmp_limited_file, $var_file)) {
                    $self->error_message('Failed to move limited variants from '. $tmp_limited_file .' to '. $var_file);
                    die $self->error_message;
                }
            }
        }
        if (-s $indel_output_file) {
            my %indel_filter_params = ( 
                indel_file => $indel_output_file, 
                out_file   => $filtered_indel_file, 
            );
            # for capture data we do not know the proper ceiling for depth
            $indel_filter_params{max_read_depth} = 1000000 if $self->region_of_interest;

            my $indel_filter = Genome::Model::Tools::Sam::IndelFilter->create(%indel_filter_params);

            unless ($indel_filter->execute) {
                $self->error_message("Running sam indel-filter failed.");
                return;
            }
        }
        else {
            Genome::Sys->write_file($filtered_indel_file);
        }
        #For now run genotype_detail only on pileup output
        $rv = $self->generate_genotype_detail_file($snv_output_file);
        unless ($rv) {
            $self->error_message('Generating genotype detail file errored out');
            die $self->error_message;
        }
    }
    return $self->verify_successful_completion($snv_output_file, $indel_output_file);
}


#need separate indel and snv by "INDEL", also need sanitize the vcf
#file (removing "N" and "." lines) 
sub create_snv_indel_output_file {
    my $self = shift;
    my $vcf_file    = $self->_vcf_staging_output;
    my $header_file = $self->_header_scratch_output;
    my $sani_file   = $self->_vcf_sani_staging_output;
    my $snv_tmp     = $self->_snv_scratch_output;
    my $indel_tmp   = $self->_indel_scratch_output;
    my $snv_file    = $self->_snv_staging_output;
    my $indel_file  = $self->_indel_staging_output;

    unless (-s $vcf_file) {
        $self->error_message("Invalid vcf file: $vcf_file");
        return;
    }

    my $vcf_fh   = Genome::Sys->open_file_for_reading($vcf_file) or return;
    my $snv_fh   = Genome::Sys->open_file_for_writing($snv_tmp) or return;
    my $indel_fh = Genome::Sys->open_file_for_writing($indel_tmp) or return;
    my $head_fh  = Genome::Sys->open_file_for_writing($header_file) or return;
    my $sani_fh  = Genome::Sys->open_file_for_writing($sani_file) or return;

    while (my $line = $vcf_fh->getline) {
        if ($line =~ /^#/) {
            $head_fh->print($line);
            next;
        }
        my @columns = split /\s+/, $line;
        if ($columns[3] eq 'N' or $columns[4] eq '.') {
            $sani_fh->print($line);
            next;
        }
        if ($columns[7] =~ /INDEL/) {
            $indel_fh->print($line);
        }
        else {
            $snv_fh->print($line);
        }
    }
    map{$_->close}($vcf_fh, $snv_fh, $indel_fh, $head_fh, $sani_fh);

    if (-z $snv_tmp) {
        $self->warning_message("snv output is empty");
        `touch $snv_file`;
    }
    else {
        Genome::Sys->cat(input_files => [$header_file, $snv_tmp],   output_file => $snv_file);
    }
    if (-z $indel_tmp) {
        $self->warning_message("indel output is empty");
        `touch $indel_file`;
    }
    else {
        Genome::Sys->cat(input_files => [$header_file, $indel_tmp], output_file => $indel_file);
    }

    return 1;
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

# decide which of mpileup or pileup will be used
# specify mpileup in snv(indel)_detection_strategy like:
# samtools r963 [mpileup]  or  samtools r963 [mpileup -u]
sub _mpileup_or_pileup {
    my $self    = shift;
    my $param   = $self->params;
    my $version = $self->version;
    my $message = "is not compatible with samtools version $version";

    if ($param =~ /^mpileup/) {
        unless ($self->is_mpileup_compatible) {
            $self->error_message("mpileup $message");
            die;
        }
        return (1, $param);
    }
    else { #for now treat all non mpileup as pileup
        unless ($self->is_pileup_compatible) {
            $self->error_message("pileup $message");
            die;
        }
        return (0, $param);
    }
}

sub is_mpileup_compatible {
    my ($ver_num) = shift->version =~ /r(\d+)/;
    return $ver_num >= 599;   #samtools r599 is the first version to have mpileup
}


sub is_pileup_compatible {
    my ($ver_num) = shift->version =~ /r(\d+)/;
    return $ver_num <= 963;   #samtools r963 is the last version to have pileup
}


sub generate_genotype_detail_file {
    my ($self, $snv_output_file) = @_; 

    unless (-f $snv_output_file) { # and -s $snv_output_file) {
        $self->error_message("SNV output File: $snv_output_file is invalid.");
        die $self->error_message;
    }

    if (not -s $snv_output_file) {
        $self->warning_message("No report input file generated for SNVs because no SNVs were detected.");
        return 1;
    }

    my $report_input_file = $self->_genotype_detail_staging_output;

    my %params = ( 
        snp_file => $snv_output_file,
        out_file => $report_input_file,
    );

    my ($is_mpileup, undef) = $self->_mpileup_or_pileup;

    if ($is_mpileup) {
        $params{snp_format} = 'vcf';
    }
    else {
        $params{snp_format} = 'sam';
    }

    my $snp_gd = Genome::Model::Tools::Snp::GenotypeDetail->create(%params);
    return $snp_gd->execute;
}


#TODO clean all of this up. It is usually/should be based on logic from
#Genome::Model::Tools::Bed::Convert logic in process_source...  this should be
#smarter about using that work ... perhaps process_source should call a method
#that just parses one line, and this method can be replaced by a call to that
#instead
sub parse_line_for_bed_intersection {
    my $class = shift;
    my $line = shift;

    unless ($line) {
        die $class->error_message("No line provided to parse_line_for_bed_intersection");
    }

    # If this is a snv line, we can just rely on the default behavior in the super class to return the right answer
    my ($chr, $pos, $ref, $var) = split "\t",  $line;
    if ( ($ref ne $var) && ($ref ne "*") ) {
        return $class->SUPER::parse_line_for_bed_intersection($line);
    }

    # Otherwise, must be parsing an indel file. Use logic from the bed converter

    my ($chromosome, $position, $star,
        $_calls, $consensus_quality, $_ref_quality, $_mapping_quality, $read_depth,
        $indel_call_1, $indel_call_2, @extra) = split("\t", $line);

    return unless $star eq '*'; #samtools indel format includes reference lines as well

    my @parsed_variants;
    for my $indel ($indel_call_1, $indel_call_2) {
        my ($reference, $variant, $start, $stop);
        next if $indel eq '*'; #Indicates only one indel call...and this isn't it!

        $start = $position;
        if(substr($indel,0,1) eq '+') {
            $reference = '*';
            $variant = substr($indel,1);
            $stop = $start; #Two positions are included-- but an insertion has no "length" so stop and start are the same
        } elsif(substr($indel,0,1) eq '-') {
            $reference = substr($indel,1);
            $variant = '*';
            $stop = $start + length($reference);
        } else {
            $class->warning_message("Unexpected indel format encountered ($indel) on line:\n$line");
            #return; skip wrong indel format line instead of failing for now
            next;
        }
        if (defined $chromosome && defined $stop && defined $reference && defined $variant) {
            push @parsed_variants, [$chromosome, $stop, $reference, $variant];
        }
    }

    unless (@parsed_variants) {
        die $class->error_message("Could not get chromosome, position, reference, or variant for line: $line");
    }

    return @parsed_variants;
}


1;

