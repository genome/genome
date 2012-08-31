package Genome::Model::Tools::Vcf::Convert::Indel::VarscanSomatic;

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::Tools::Vcf::Convert::Indel::VarscanSomatic {
    is  => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from Varscan somatic indel output',
};


sub help_synopsis {
    <<'HELP';
    Generate a VCF file from Varscan somatic indel output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the indels.
HELP
}

sub source {
    return 'VarscanSomatic';
}

#Varscan somatic indel output columns:
#chrom					chromosome name
#position				position (1-based from the pileup)
#ref					reference allele at this position
#var					variant allele at this position
#normal_reads1			reads supporting reference allele
#normal_reads2			reads supporting variant allele
#normal_var_freq		frequency of variant allele by read count
#normal_gt				genotype call for Normal sample
#tumor_reads1			reads supporting reference allele
#tumor_reads2			reads supporting variant allele
#tumor_var_freq			frequency of variant allele by read count
#tumor_gt				genotype call for Tumor sample
#somatic_status			status of variant (Germline, Somatic, or LOH)	
#variant_p_value		Significance of variant read count vs. baseline error rate
#somatic_p_value		Significance of tumor read count vs. normal read count
#tumor_reads1_plus      Ref-supporting reads from + strand in tumor
#tumor_reads1_minus     Ref-supporting reads from - strand in tumor
#tumor_reads2_plus      Var-supporting reads from + strand in tumor
#tumor_reads2_minus		Var-supporting reads from - strand in tumor

#Examples:
#1	797165	A	+T	15	7	31.82%	*/+T	4	3	42.86%	*/+T	Germline	3.8386754900572427E-4	0.4586245338869044	8	0	3	0
#1	805541	G	-GGGTT	5	0	0%	G	6	2	25%	*/-GGGTT	Somatic	1.0	0.3589743589743588	4	3	0	2
#1	805664	A	+G	51	12	19.05%	A	22	15	40.54%	*/+G	Somatic	1.0	0.018486781242907306	23	14	8	7
#1	805674	C	-TTCA	59	4	6.35%	C	29	12	29.27%	*/-TTCA	Somatic	1.0	0.002064031614857694	25	18	7	5
#1	806275	C	-T	16	10	38.46%	*/-T	10	10	50%	*/-T	Germline	6.71709572022225E-8	0.31453504019953177	11	0	10	0
#1  965319  A   -TG     3       3       50%     */-TG   38      2       5%      A       LOH     1.0     0.011822690285784018    19      21      1       1
#1  965346  C   -TG     3       2       40%     */-TG   37      0       0%      C       LOH     1.0     0.01161440185830484     22      15      0       0
#1	790696	C	+AT	    3	    0	    0%	    C	    3	    29	90.62%	+AT/+AT	Somatic	    1.0	    0.003055767761650106	20	    12	    17	    12
#
#somatic status (SS): 0, 1, 2, 3, or 4 depending on whether relative to normal the variant is wildtype, germline, somatic, LOH or unknown respectively.


sub parse_line {
    my ($self, $line) = @_;

    my @columns = split /\s+/, $line;
    my ($ref, $alt);

    if ($columns[3] =~ /^\+(\S+)/) {
        $alt = $columns[2].$1;
        $ref = $columns[2];
    }
    elsif ($columns[3] =~ /^\-(\S+)/) {
        $ref = $columns[2].$1;
        $alt = $columns[2];
    }
    else {
        die $self->error_message("Invalid indel type for line: $line");
    }

    my $n_gt = $self->_convert_gt($columns[2], $columns[3], $columns[7],  $columns[12]);
    my $t_gt = $self->_convert_gt($columns[2], $columns[3], $columns[11], $columns[12]);

    unless ($n_gt and $t_gt) {
        die $self->error_message("Invalid gt str in line: $line");
    }

    my $n_dp  = $columns[4] + $columns[5];
    my $t_dp  = $columns[8] + $columns[9];
    my $n_dp4 = '.';
    my $t_dp4 = join ',', $columns[15], $columns[16], $columns[17], $columns[18];

    my %ss = (
        WILDTYPE => 0,
        GERMLINE => 1,
        SOMATIC  => 2,
        LOH      => 3,
        UNKNOWN  => 4,
    );

    my $t_ss  = $ss{uc($columns[12])};

    unless (defined $t_ss) {
        die $self->error_message("Failed to get somatic status from line: $line");
    }

    my $n_str = join ':', $n_gt, $n_dp, $n_dp4, '.', '.';
    my $t_str = join ':', $t_gt, $t_dp, $t_dp4, '.', $t_ss;

    my $new_line = join "\t", $columns[0], $columns[1], '.', $ref, $alt, '.', 'PASS', '.', 'GT:DP:DP4:BQ:SS', $n_str, $t_str;
    return $new_line;
}


sub _convert_gt {
    my ($self, $ref, $var, $gt, $ss) = @_;

    if ($gt =~ /^[ATCG]$/) {
        unless ($ref eq $gt) {
            $self->error_message("Genotype: $gt conflict with reference: $ref");
            return;
        }
        return '0/0';
    }

    if ($gt =~ /\//) {
        my @alleles = split /\//, $gt;
        my @new_gt;

        for my $allele (@alleles) {
            if ($allele eq '*') {
                push @new_gt, 0;
            }
            elsif ($allele =~ /[\+\-]/) {
                unless (uc($ss) eq 'UNKNOWN') { #Special case that allele does not exactly matches var
                    unless ($allele eq $var) {
                        #put warning here for now for the conflict, probably a Varscan bug 
                        #2	233744130	G	-CC	6	3	33.33%	*/-CCC	1	10	90.91%	-CC/-CC	LOH	1.0	0.012383900928792655	10	28	2
                        $self->warning_message("Genotype: $gt conflict with variant: $var");                
                    }
                }
                push @new_gt, 1;
            }
            else {
                $self->error_message("Invalid genotype: $gt");
                return;
            }
        }
        return join '/', sort{$a <=> $b} @new_gt;
    }
}

1;

