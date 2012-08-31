package Genome::Model::Tools::Snp::GenotypeDetail;

use strict;
use warnings;

use Genome;
use Command;
use Genome::Info::IUB;


class Genome::Model::Tools::Snp::GenotypeDetail {
    is  => 'Command',
    has => [
    snp_file => {
        is  => 'String',
        doc => 'The input sam/bam snp file',
    },
    snp_format => {
        is  => 'String',
        doc => 'Must be either sam or maq or vcf',
    },
    ],
    has_optional => [
    maq_pileup_file => {
        is  => 'String',
        doc => 'maq pileup output file used to figure number of reads calling the base, it is required for maq snp_format',
    },
    out_file => {
        is  => 'String',
        doc => 'The output snp genotype detail file. if not given, a output will be created in the same dir as snp-file with .genotype_detail as suffix',
    },
    ],
};


sub help_brief {
    'Get detailed snp genotype info';
}

sub help_detail {
    return <<EOS
    Get detailed SNP genotype info, also used in reference alignment solexa pipeline, 
    step Variant Detection
EOS
}

$ENV{UR_COMMAND_DUMP_STATUS_MESSAGES}=1;


sub execute {
    my $self = shift;
    my $maq_pileup_fh;

    my $snp_format = uc $self->snp_format;
    my $maq_pileup = $self->maq_pileup_file;
    my $snp_file   = $self->snp_file;
    my $out_file   = $self->out_file || $snp_file.'.genotype_detail';

    if ($snp_format =~ /^(MAQ|SAM|VCF)$/) { #VCF file is for the mpileup runs
        if ($snp_format eq 'MAQ') {
            $self->error_message("Maq pileup file: $maq_pileup is invalid") and return 
            unless $maq_pileup and -s $maq_pileup;
            $maq_pileup_fh = Genome::Sys->open_file_for_reading($maq_pileup);
            $self->error_message("Fail to open maq pileup: $maq_pileup for writing") and return
            unless $maq_pileup_fh;
        }
    }
    else {
        $self->error_message("Unrecognized snp format: $snp_format");
        return;
    }

    unless (-s $snp_file) {
        $self->error_message('Can not find valid snp file: '.$snp_file);
        return;
    }

    my $snp_fh = Genome::Sys->open_file_for_reading($snp_file);
    unless ($snp_fh) {
        $self->error_message("Failed to open snp out file $snp_file");
        return;
    }

    my $out_fh = Genome::Sys->open_file_for_writing($out_file);
    unless ($out_fh) {
        $self->error_message("Failed to open genotype detail out file $out_file for writing");
        return;
    }

    while (my $snp_line = $snp_fh->getline) {
        my @columns = split /\s+/, $snp_line;
        my ($chr, $pos, $ref, $cns, $cns_qual, $ref_base_ct, $read_bases);

        if ($snp_format eq 'SAM') {
            ($chr, $pos, $ref, $cns, $cns_qual, $read_bases) = map{$columns[$_]}(0..4, 8);
        }

        elsif ($snp_format eq 'VCF') { #mpileup.
            next if $snp_line=~ /^#/;
            ($chr, $pos, $ref, $cns, $cns_qual)= map{$columns[$_]}(0..3,6);
            $read_bases = '';
        }

        else {
            ($chr, $pos, $ref, $cns, $cns_qual) = map{$columns[$_]}(0..4);
            my $pileup_line = $maq_pileup_fh->getline;
            my ($pu_chr, $pu_pos, $pu_ref, $pu_read_bases) =  split("\t", $pileup_line);
            $read_bases = $pu_read_bases;

            unless ($chr eq $pu_chr and $pos == $pu_pos and $ref eq $pu_ref) {
                $self->error_message("Data is not consistent. Maq SNP file line ".$snp_fh->input_line_number." maq pileup line ".$maq_pileup_fh->input_line_number);
                return:
            }
        } 

        $read_bases  = uc $read_bases;
        $ref_base_ct = $read_bases =~ tr/\.\,//;

        my %base_ct = (
            A => $read_bases =~ tr/A//,
            T => $read_bases =~ tr/T//,
            C => $read_bases =~ tr/C//,
            G => $read_bases =~ tr/G//,
        ); 
        my @cns_bases;
        @cns_bases = split ('',$cns);

        my @alleles =map{Genome::Info::IUB->iub_to_alleles($_)}@cns_bases;
        unless (@alleles) {
            $self->error_message("Failed to find IUB code for base $cns");
            return;
        }
        $alleles[1] = 'X' if @alleles > 2; #Maq called more than 2 possible alleles? The original code specified that we should report the variation as 'X'


        for my $i (0..1) {
            my $base = $alleles[$i];
            next if $base eq $ref;
            $ref_base_ct ||= 0;
            $base_ct{$base} ||= 0;
            $out_fh->print(join("\t", $chr, $pos, $pos, $ref, $base, 'ref', 'SNP', $ref_base_ct, $base_ct{$base}, $cns_qual), "\n");

            last if $alleles[0] eq $alleles[1];
        }

    }
    map{$_->close}($snp_fh, $out_fh);
    $maq_pileup_fh->close if $maq_pileup_fh;

    return 1;
}

1;
