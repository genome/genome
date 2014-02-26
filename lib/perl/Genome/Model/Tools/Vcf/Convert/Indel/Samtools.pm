package Genome::Model::Tools::Vcf::Convert::Indel::Samtools;

use strict;
use warnings;
use Genome;
use POSIX 'strftime';

class Genome::Model::Tools::Vcf::Convert::Indel::Samtools {
    is  => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from Samtools indel output',
    has => [
        is_mpileup => {
            type          => 'Boolean',
            is_transient  => 1,
            default_value => 0,
        },
    ],
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from samtools indel output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the indels.
HELP
}


sub source {
    return 'Samtools';
}

#single sample for now
sub _get_header_columns {
    my $self = shift;
    my @header_columns = ("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",$self->aligned_reads_sample);
    return @header_columns;
}


sub convert_file {
    my $self = shift;

    my $input_file = $self->input_file;
    my $token = `head -1 $input_file`;

    $self->is_mpileup(1) if $token =~ /^#/;
    $self->SUPER::convert_file;

    return 1;
}


sub parse_line {
    my ($self, $lines) = @_;
    
    if ($self->is_mpileup) {
        return if $lines =~ /^#/;
        my @columns = split("\t", $lines);
        my ($DP, $MQ) = $columns[7] =~ /DP=(\d+);\S+MQ=(\d+);/;
        unless (defined $DP and defined $MQ) {
            $self->warning_message("Failed to get DP and MQ for line:\n$lines");
            return;
        }
        my ($GT) = $columns[9] =~ /^(\S+?)\:/;
        $columns[6] = 'PASS';
        $columns[7] = '.';
        $columns[8] = 'GT:DP:MQ';
        $columns[9] = join ':', ($GT, $DP, $MQ);

        my $col_ct   = scalar $self->_get_header_columns;
        my $new_line = join "\t", splice(@columns, 0, $col_ct); #remove unwanted columns in test
        return $new_line;
    }

    my @lines = split "\n", $lines;

    my @first_line  = split "\t", $lines[0];
    my @second_line = split "\t", $lines[1];

    my $ref = $first_line[2];

    my $chr = $second_line[0];
    my $pos = $second_line[1];
    my $leading_base = $ref;
    my $indel_string = $second_line[3];
    my $consensus_quality = $second_line[4];
    my $ref_quality = $second_line[5];
    my $mapping_quality = $second_line[6];
    my $read_depth = $second_line[7];
    my $indel_call_1 = $second_line[8];
    my $indel_call_2 = $second_line[9];
    my $allele_depth_1 = $second_line[10];
    my $allele_depth_2 = $second_line[11];

    my @alt_alleles;
    my $alt_allele;
    my $ref_allele;

    if ($indel_call_1 ne '*') {
        push(@alt_alleles, $indel_call_1);
    }
    if ($indel_call_2 ne '*') {
        push(@alt_alleles, $indel_call_2);
    }
    if (@alt_alleles == 0) {
        die $self->error_message("No indel calls were made on this line: $indel_call_1/$indel_call_2");
    }
    
    # Partition calls into insertions and deletions
    my @insertions;
    my @deletions;
    for my $alt (@alt_alleles) {
        if ($alt =~ /^\+/) {
            push(@insertions, substr($alt, 1));
        }
        elsif ($alt =~ /^\-/) {
            push(@deletions, substr($alt, 1));
        }
        else {
            die $self->error_message("Insertion/deletion type not recognized: ".$alt_alleles[0]);
        }
    }

    @deletions = sort { length($b) <=> length($a) } @deletions;

    my $ref_tail = "";
    if (!@deletions) {
        $ref_allele = $leading_base;
    } else {
        $ref_allele = $leading_base . $deletions[0];
        $ref_tail = $deletions[0];
    }

    my @final_alts;
    for my $ins (@insertions) {
        push(@final_alts, $leading_base . $ins . $ref_tail);
    }

    for my $del (@deletions) {
        my $tail = "";
        if (length($del) < length($ref_tail)) {
            $tail = substr($ref_tail, length($del));
        }
        push(@final_alts, $leading_base . $tail);
    }

    $alt_allele = join(",", @final_alts);


    #TODO this is turned off for now because it interferes with applying filters (bed coordinates will be different once left shifted)
    # ($chr, $pos, $ref_allele, $alt_allele) = $self->normalize_indel_location($chr, $pos, $ref_allele, $alt_allele);
    
    my $GT;
    my @indel_string_split = split(/\//, $indel_string);
    if (@indel_string_split != 2) {
        $self->warning_message("Genotype in unexpected format: $indel_string at chr $chr pos $pos");
        return;
    }
    if ($indel_string eq "*/*") {
        $GT = "0/0";
    }
    elsif ($indel_string =~/\*/){
        $GT = "0/1";
    }
    elsif ($indel_string_split[0] eq $indel_string_split[1]) {
        $GT = "1/1";
    }
    else {
        $GT = "1/2";
    }

    my $DP = $read_depth;
    my $MQ = $mapping_quality;

    my $filter = "PASS";

    my $dbsnp_id = ".";
    my $qual = ".";
    my $info = ".";

    my $format = "GT:DP:MQ";
    my $sample_string = join (":", ($GT,$DP,$MQ));
    my $vcf_line = join("\t", $chr, $pos, $dbsnp_id, $ref_allele, $alt_allele, $qual, $filter,
                        $info, $format, $sample_string);
    return $vcf_line;
}

sub get_record {
    my ($self, $input_fh) = @_;

    if ($self->is_mpileup) { #for mpileup just need return one vcf line
        return $self->SUPER::get_record($input_fh);
    }
    #For samtools indel, we need to get two lines at a time.
    my $lines;
    my $line1 = $input_fh->getline; 
    my $line2;
    my $num_lines = 1;
    while ($line1 && $num_lines < 2) { 
        $line2 = $input_fh->getline;

        #Check to make sure the lines are correctly paired
        if ($line2) {
            my @fields1 = split (/\t/, $line1);
            my @fields2 = split (/\t/, $line2);
            if (($fields1[0] eq $fields2[0]) && ($fields1[1] eq $fields2[1])) {
                $lines = $line1.$line2;
                $num_lines++;
            }
            else {
                $line1 = $line2;
            }
        }
        else { #The file ended, so we couldn't get line2
            return undef;
        }
    }

    return $lines;
}

sub get_format_meta {
    my $self   = shift;
    my $MQ_tag = {MetaType => "FORMAT", ID => "MQ",  Number => 1, Type => "Integer", Description => "Phred style probability score that the variant is novel with respect to the genome's ancestor"};

    return ($self->common_format_meta, $MQ_tag);
}

1;

