package Genome::Model::Tools::Vcf::Convert::Indel::SniperSlow;

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::Tools::Vcf::Convert::Indel::SniperSlow {
    is  => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from Sniper indel output',
};


sub help_synopsis {
    <<'HELP';
    Generate a VCF file from Sniper indel output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the indels.
HELP
}

sub source {
    return 'Sniper';
}


sub parse_line {
    my ($self, $line) = @_;

    my @columns = split /\t/, $line;
    my ($chr, $pos, $call_1, $call_2, $len_1, $len_2, $t_gt, $t_dp, $t_read_1, $t_read_2, $n_gt, $n_dp, $n_read_1, $n_read_2) = map{$columns[$_]}qw(0 1 4 5 6 7 8 12 13 14 21 25 26 27);
    
    my ($t_ss, $n_ss) = _parse_ss($t_gt, $n_gt);
    die $self->error_message("Failed to parse ss for line: $line") unless defined $t_ss and defined $n_ss;

    my $rv = $self->_convert_indel($chr, $pos, $call_1, $len_1, $t_gt, $n_gt, $t_dp, $n_dp, $t_read_1, $n_read_1, $t_ss, $n_ss);
    return $rv if $rv;
    $rv = $self->_convert_indel($chr, $pos, $call_2, $len_2, $t_gt, $n_gt, $t_dp, $n_dp, $t_read_2, $n_read_2, $t_ss, $n_ss);
    return $rv if $rv;

    $self->warning_message("$line can not be parsed to vcf format");
    return;
}


sub _convert_indel {
    my ($self, $chr, $pos, $call, $length, $t_gt, $n_gt, $t_dp, $n_dp, $t_read_ct, $n_read_ct, $t_ss, $n_ss) = @_;
    unless($call) {
        $self->warning_message("No call made");
        return;
    }
    return if $call eq '*'; #Indicates only one indel call...and this isn't it!

    my $ref_seq      = $self->reference_sequence_input;
    my $samtool_path = Genome::Model::Tools::Sam->path_for_samtools_version;
    my $region       = "$chr:$pos-$pos";

    my $ref_base = `$samtool_path faidx $ref_seq $region | grep -v \">\"`;
    chomp $ref_base;

    my $alt_base;

    if ($length > 0) {
        $alt_base = $ref_base . $call;
    }
    elsif ($length < 0) {
        $alt_base  = $ref_base;
        $ref_base .= $call;
    }
    else {
        die $self->error_message("$length is not valid for indel length at $chr : $pos");
    }

    $t_gt = _parse_gt($t_gt);
    $n_gt = _parse_gt($n_gt);

    return join "\t", $chr, $pos, '.', $ref_base, $alt_base, '.', 'PASS', '.', 'GT:DP:AD:BQ:SS', $n_gt.':'.$n_dp.':'.$n_read_ct.':.:'.$n_ss, $t_gt.':'.$t_dp.':'.$t_read_ct.':.:'.$t_ss;
}

#get ss
#*/*   */*   =>  wildtype 
#*/+T  */+T  =>  germline
#*/*   */+T  =>  wildtype
#*/+T  */*   =>  somatic
sub _parse_ss {
    my ($t_gt, $n_gt) = @_;

    return (0, '.') if $t_gt eq '*/*';
    return (1, '.') if $t_gt =~ /[\+\-]/ and $n_gt =~ /[\+\-]/;
    return (2, '.') if $t_gt =~ /[\+\-]/ and $n_gt eq '*/*';
    return (undef, undef);
}


sub _parse_gt {
    my $gt = shift;
    my @alleles = split /\//, $gt;

    my @gt;
    for my $allele (@alleles) {
        if ($allele eq '*') {
            push @gt, 0;
        }
        else {
            push @gt, 1;
        }
    }
    return join '/', sort{$a <=> $b}@gt;
}


1;

