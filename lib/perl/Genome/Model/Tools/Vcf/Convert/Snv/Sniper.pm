package Genome::Model::Tools::Vcf::Convert::Snv::Sniper;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;
use File::Basename;

class Genome::Model::Tools::Vcf::Convert::Snv::Sniper {
    is  => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from sniper output',
    has_transient_optional => [
        vcf_ok => {
            is => 'Boolean',
        },
    ],
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from sniper snv output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the snvs.
HELP
}

sub source {
    my $self = shift;
    return "Sniper";
}

sub initialize_filehandles {
    my $self = shift;

    if ($self->_input_fh || $self->_output_fh) {
        return 1; #Already initialized
    }

    my $input  = $self->input_file;
    my $output = $self->output_file;

    my $line  = `head -1 $input`;
    my @cols  = split /\t/, $line;
    my $value = scalar @cols == 26 ? 1 : 0;
    $self->vcf_ok($value);

    #Sniper specific VCF output. For now modify this file. Hacky way
    if ($self->vcf_ok) {
        $self->debug_message("Assume this sniper run also gets vcf output");
        my $dir = dirname $input;
        my ($raw_vcf) = glob($dir . "/*raw.vcf"); 

        unless (-e $raw_vcf) {
            die $self->error_message("sniper raw.vcf is not available under the input directory");
        }

        my $sorted_vcf = Genome::Sys->create_temp_file_path;

        my $sort_cmd = Genome::Model::Tools::Bed::ChromSort->create(
            input => $raw_vcf,
            output => $sorted_vcf,
        );
        unless ($sort_cmd->execute) {
            $self->error_message("Couldn't sort Sniper vcf");
            return;
        }

        $input = $sorted_vcf;
}

    my $input_fh  = Genome::Sys->open_file_for_reading($input)
        or die "Failed to open input $input for reading\n";
    my $output_fh = Genome::Sys->open_gzip_file_for_writing($output) 
        or die "Failed to open output $output for writing\n";
    
    $self->_input_fh($input_fh);
    $self->_output_fh($output_fh);

    return 1;
}

sub parse_line {
    my ($self, $line) = @_;
    return if $line =~ /^#/; # no vcf header here
    
    if ($self->vcf_ok) { #keep original sniper vcf
        my @columns = split /\t/, $line;
        $columns[6] = 'PASS';
        return join "\t", @columns;
    }

    # TODO snv_qual is not used...
    my ($chr, $pos, $ref, $genotype, $tumor_vaq, $snv_qual, $tumor_mq, $tumor_bq, $tumor_dp, $normal_dp ) = split("\t",$line);
    my $tumor_gq = $tumor_mq; #TODO Should these really be the same?

    #replace ambiguous/IUPAC bases with N in ref
    $ref =~ s/[^ACGTN\-]/N/g;

    my @alt_alleles = Genome::Info::IUB->variant_alleles_for_iub($ref, $genotype);
    my @alleles = Genome::Info::IUB->iub_to_alleles($genotype);
    my $alt = join(",", @alt_alleles);

    #add the ref and alt alleles' positions in the allele array to the GT field
    my $tumor_gt = $self->generate_gt($ref, \@alt_alleles, \@alleles);

    # We do not have access to much of the normal information from somatic output
    my $normal_gt = "./.";
    # genotype quality (consensus quality)
    my $normal_gq = ".";
    # avg mapping quality ref/var
    my $normal_mq = ".";
    # avg mapping quality ref/var
    my $normal_bq = ".";
    # allele depth
    my $normal_ad =  ".";
    my $tumor_ad =  ".";
    # fraction of reads supporting alt
    my $normal_fa =  ".";
    my $tumor_fa =  ".";
    # vaq
    my $normal_vaq = ".";
    # somatic status
    my $normal_ss = ".";
    my $tumor_ss  = 2;

    # Placeholder for later adjustment
    my $dbsnp_id = ".";
    my $qual = "."; # Can also be $tumor_vaq
    my $filter = "PASS";
    my $format = "GT:GQ:DP:BQ:MQ:AD:FA:VAQ:SS";
    my $info = ".";
    my $tumor_sample_string = join (":", ($tumor_gt, $tumor_gq, $tumor_dp, $tumor_bq, $tumor_mq, $tumor_ad, $tumor_fa, $tumor_vaq, $tumor_ss));
    my $normal_sample_string = join (":", ($normal_gt, $normal_gq, $normal_dp, $normal_bq, $normal_mq, $normal_ad, $normal_fa, $normal_vaq, $normal_ss));

    my $vcf_line = join("\t", $chr, $pos, $dbsnp_id, $ref, $alt, $qual, $filter, $info, $format, $normal_sample_string, $tumor_sample_string);

    return $vcf_line;
}


sub get_format_meta {
    my $self = shift;

    # Get all of the base FORMAT lines
    my @tags = $self->SUPER::get_format_meta;
    return @tags unless $self->vcf_ok;

    my $amq = {MetaType => "FORMAT", ID => "AMQ",    Number => ".", Type => "Integer", Description => "Average mapping quality for each allele present in the genotype"};
    my $igt = {MetaType => "FORMAT", ID => "IGT",    Number => 1,   Type => "String",  Description => "Genotype when called independently (only filled if called in joint prior mode)"};
    my $bct = {MetaType => "FORMAT", ID => "BCOUNT", Number => 4,   Type => "Integer", Description => "Occurrence count for each base at this site (A,C,G,T)"};
    my $jgq = {MetaType => "FORMAT", ID => "JGQ" ,   Number => 1,   Type => "Integer", Description => "Joint genotype quality (only filled if called in join prior mode)"};
    my $ssc = {MetaType => "FORMAT", ID => "SSC",    Number => 1,   Type => "Integer", Description => "Somatic score between 0 and 255"};
    
    return (@tags, $amq, $igt, $bct, $jgq, $ssc);
}


