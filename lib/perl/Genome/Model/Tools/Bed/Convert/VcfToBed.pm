package Genome::Model::Tools::Bed::Convert::VcfToBed;

use strict;
use warnings;

use Carp qw/confess/;
use Genome;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Genotype;
use Genome::File::Vcf::Entry;
use Genome::Info::IUB;

class Genome::Model::Tools::Bed::Convert::VcfToBed {
    is => ['Genome::Model::Tools::Bed::Convert'],
};

sub process_source {
    my $self = shift;
    my $vcf_reader = Genome::File::Vcf::Reader->fhopen($self->_input_fh, $self->source);

    while(my $entry = $vcf_reader->next) {
        #do something intelligent here
        next if $entry->is_filtered;
        #for now, only look at the first sample
        my $ft = $entry->sample_field(0,"FT");
        next if(defined $ft and ($ft ne '.' and $ft ne 'PASS'));

        #ok, we're converting this bad boy!
        my $gt = $entry->genotype_for_sample(0);
        my @entry_alleles = $entry->alleles;
        my @genotype_alleles = map {$entry_alleles[$_]} $gt->get_alleles;

        if($entry->has_indel) {
            my ($indel_gts, $indel_shifts) = $self->_convert_indel_gt_to_bed($entry->{reference_allele}, @genotype_alleles);
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

sub _convert_snv_gt_to_bed {
    my ($self, $ref, @alleles) = @_;
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
sub _convert_indel_gt_to_bed {
    my ($self, $reference_allele, @alleles) = @_;

    #the challenge here is that VCF can represent multi-base calls in both the reference and variant and we clearly do NOT support that in our tgi-bed representation

    my %alleles = map { $_ => 1 } @alleles; #our bed doesn't care about indel genotypes.

    my @bed_alleles;
    my @position_shifts;
    for my $allele (keys %alleles) {
        next if $allele eq $reference_allele;
        my ($ref, $var, $right_shift) = $self->_simplify_indel_allele($reference_allele, $allele);
        unless($var eq q{} || $ref eq q{}) {
            confess "Complex indels cannot be converted to TGI bed\n";
        }
        ($ref, $var) = map { $_ ne q{} ? $_ : '*' } ($ref, $var);

        push @bed_alleles, [$ref,$var];
        push @position_shifts, $right_shift;
    }
    return (\@bed_alleles, \@position_shifts);
}

sub _simplify_indel_allele {
    my ($self, $ref, $var) = @_;
    #these could be padded e.g. G, GT for a T insertion or GCC G for a 2bp deletion
    #they could also be complex e.g. GCCCGT, GCGT for a 2bp deletion
    #they could also represent an insertion deletion event e.g. GCCCGT GCGGGGT; these cannot be represented in genome bed. Throw an error or warn.
    #
    #I think the algorithm should be trim end (no updating of coords required)
    #trim beginning and return number of bases trimmed

    my @ref_array = map { uc } split //, $ref;
    my @var_array = map { uc } split //, $var;

    while(@ref_array and @var_array and $ref_array[-1] eq $var_array[-1]) {
        pop @ref_array;
        pop @var_array;
    }

    my $right_shift = 0;
    while(@ref_array and @var_array and $ref_array[0] eq $var_array[0]) {
        shift @ref_array;
        shift @var_array;
        $right_shift++;
    }

    return (join("",@ref_array), join("",@var_array), $right_shift);
}

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
