package Genome::Model::Tools::Bed::Convert::VcfToBed;

use strict;
use warnings;

use Carp qw/confess/;
use Genome;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Genotype;
use Genome::File::Vcf::Entry;
use Genome::Info::IUB;
use Genome::Utility::Vcf qw(convert_indel_gt_to_bed);

class Genome::Model::Tools::Bed::Convert::VcfToBed {
    is => ['Genome::Model::Tools::Bed::Convert'],
    has => [
        sample_name => {
            is => "Text",
            is_optional => 1,
            doc => "Name of specific sample to extract from VCF.  If not specified, the first sample will be used",
        },
        remove_filtered_calls => {
            is => 'Boolean',
            default => 1,
            doc => "Do not convert calls with a filtered status",
        }
    ],
};

sub process_source {
    my $self = shift;
    my $vcf_reader = new Genome::File::Vcf::Reader($self->source);

    my $sample_index = _sample_index_for_name($vcf_reader, $self->sample_name);

    while(my $entry = $vcf_reader->next) {
        # Skip this entry if there is no data for the sample
        next if (scalar @{$entry->sample_data->[$sample_index]} == 0);

        if ($self->remove_filtered_calls) {
            next if $entry->is_filtered;
            my $ft = $entry->sample_field($sample_index,"FT");
            next if(defined $ft and ($ft ne '.' and $ft ne 'PASS'));
        }

        my @genotype_alleles = $entry->bases_for_sample($sample_index);

        if($entry->has_indel) {
            my ($indel_gts, $indel_shifts) = convert_indel_gt_to_bed($entry->{reference_allele}, @genotype_alleles);
            for my $indel (@$indel_gts) {
                my ($ref, $var) = @$indel;
                my $shift = shift @$indel_shifts;
                my $pos = $entry->{position} + $shift;
                if($ref eq '*') {
                    #insertion
                    $self->write_bed_line($entry->{chrom}, $pos-1, $pos-1, $ref, $var, '-', '-');
                }
                if($var eq '*') {
                    #deletion
                    $self->write_bed_line($entry->{chrom}, $pos-1, $pos-1 + length($ref), $ref, $var, '-', '-');
                }
            }
        }
        else {
            my $bed_gt = $self->_convert_snv_gt_to_bed($entry->{reference_allele}, @genotype_alleles);
            $self->write_bed_line($entry->{chrom}, $entry->{position}-1, $entry->{position}, @$bed_gt, '-', '-');
        }
    }
    return 1;
}

sub _sample_index_for_name {
    my $reader = shift;
    my $name = shift;

    unless ($name) {
        return 0;
    }

    return $reader->header->index_for_sample_name($name);
}

sub _convert_snv_gt_to_bed {
    my ($self, $ref, @alleles) = @_;
    if(@alleles == 1) {
        push  @alleles, $alleles[0];
    }
    my $ret = Genome::Info::IUB->iub_for_alleles(@alleles); #this will return undef if the ploidy is off;
    unless($ret) {
        confess "Unable to convert SNV to IUB code";
    }
    else {
        return [$ref, $ret];
    }
}

#1       929701  929701  0/A     -       -
#1       966294  966294  */TG    5       20
#1       1302523 1302523 */T     2       18
#1       1302526 1302526 */C     2       17
#1       8979303 8979307 AAAT/*  13      17
#1       9288000 9288002 AA/0    -       -
#1       9369360 9369364 GTGT/*  12      6
#1       9936903 9936905 AC/*    11      21


sub help_brief {
    "Tool to convert the first sample in a VCF to TGI BED.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed convert vcf-to-bed...
EOS
}

sub help_detail {                           
    return <<EOS
    This tool takes a VCF file and converts the first sample in it to a TGI variant BED file with no information about the depth and quality
EOS
}


1;
