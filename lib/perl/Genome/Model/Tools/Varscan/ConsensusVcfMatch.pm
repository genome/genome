package Genome::Model::Tools::Varscan::ConsensusVcfMatch;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Varscan::ConsensusVcfMatch {
    is => ['Genome::Command::Base'],
    has_input => [
        indel_vcf => {
            is => 'Text',
            doc => '', #TODO: include some doc
        },
        sample_cns_list_file => {
            is => 'Text',
            doc => '', #TODO: include some doc
        },
        output_file => {
            is => 'Text',
            is_output => '1',
            doc => '', #TODO: include some doc
        },
        output_cns_file => {
            is => 'Text',
            is_output => '1',
            doc => '', #TODO: include some doc
        },
    ],
    has_optional_transient => [
        stats => { is => 'Hash', default_value => {} },
        target_indel_genotypes => { is => 'Hash',default_value => {}  },
        target_indels => { is => 'Hash', default_value => {}  },
        samples => { is => 'Hash' ,default_value => {} },
    ],
};

sub execute {
    my $self = shift;
    $self->_initialize;
    my %target_indels = %{$self->target_indels};

	my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
	my $output_cns_fh = Genome::Sys->open_file_for_writing($self->output_cns_file);
    print $output_fh join("\t", 'chrom', 'chr_start', 'chr_stop', 'ref', 'alt', 'obs_hom1', 'obs_het', 'obs_hom2', 'exp_hom1', 'exp_het', 'exp_hom2', 'chi_sum', 'original_fet', 'varscan_ref', 'varscan_het', 'varscan_hom', 'varscan_hwe'), "\n";

    foreach my $key (keys %target_indels) {
        if($target_indels{$key}) {
            my @targetIndels = split(/\n/, $target_indels{$key});
            
            foreach my $target_indel (@targetIndels) {
                my $num_samples = 0;
                my $num_samples_missing = 0;
                my $num_samples_ref = 0;
                my $num_samples_var = 0;
                my $num_samples_het = 0;
                my $num_samples_hom = 0;
                my $num_samples_snp = 0;
                my $num_samples_other = 0;
                
                foreach my $sample (sort keys %{$self->samples}) {
                    my $indel_key = join("\t", $sample, $target_indel);
                    $num_samples++;
                    if($self->target_indel_genotypes->{$indel_key}) {
                        my ($call, $cns, $reads1, $reads2, $var_freq) = split(/\:/, $self->target_indel_genotypes->{$indel_key});
                        
                        print $output_cns_fh join("\t", $indel_key, $call, $cns, $reads1, $reads2, $var_freq) , "\n";

                        if($call eq "Het" || $call eq "Hom") {
                            $num_samples_var++;
                            $num_samples_het++ if($call eq "Het");
                            $num_samples_hom++ if($call eq "Hom");
                        }
                        elsif($call eq "Ref") {						
                            $num_samples_ref++;
                        }
                        elsif($call eq "SNP") {
                            $num_samples_snp++;
                        }
                        elsif($call eq "Other") {
                            $num_samples_other++;
                        }
                        else {
                            die $self->error_message("Unrecognized genotype " .$self->target_indel_genotypes->{$indel_key} . " for indel $indel_key");
                        }
                    }
                    else {
                        $num_samples_missing++;
                    }
                }
                
                my $vaf = "NA";
                if($num_samples_ref || $num_samples_het || $num_samples_hom) {
                    $vaf = ($num_samples_het + 2 * $num_samples_hom) / ($num_samples_ref * 2 + $num_samples_het * 2 + $num_samples_hom * 2);
                    $vaf = sprintf("%.3f", $vaf);
                    
                    ## Run the hwe calculation ##
                    
                    my $hwe = $self->_get_hwe_p($num_samples_het, $num_samples_ref, $num_samples_hom);

                    print join("\t", $target_indel, $num_samples_ref, $num_samples_het, $num_samples_hom, $vaf, $hwe) . "\n";
                    print $output_fh join("\t", $target_indel, $num_samples, $num_samples_ref, $num_samples_het, $num_samples_hom, $hwe) , "\n";
                }
            }
        }
    }

    $output_fh->close;
    $output_cns_fh->close;

    my %stats = %{$self->stats};
    foreach my $key (sort keys %stats) {
        print $stats{$key} . "\t$key\n";
    }
    return 1;
}

sub _initialize {
    my $self = shift;
    $self->target_indels($self->_load_indels);
    $self->samples($self->_parse_sample_list);
    return 1;
}

sub _load_indels {
    my $self = shift;
    my $indel_file = $self->indel_vcf;
	my %indels = ();
	
	my $fh = Genome::Sys->open_file_for_reading($indel_file);
	
	while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^#/;

        my ($chrom, $position, $id, $ref, $alt, $qual, $filter, $observed, $expected, $chi_sum, $original_p) = split(/\t/, $line);
			
        next if $filter ne 'PASS';

        my @vars = split(/\,/, $alt);
        my $key = join("\t", $chrom, $position);

        foreach my $var (@vars) {				
            my $variant_type = "indel";
            my $chr_start = 0;
            my $chr_stop = 0;
            my ($allele1,  $allele2) = ('', '');
            if(length($var) > length($ref)) {
                $variant_type = "ins";
                ## Insertion ##
                $allele2 = $var;
                $allele2 =~ s/$ref//;
                $allele1 = "-";
                $chr_start = $position;
                $chr_stop = $position + 1;
            }
            else {
                $variant_type = "del";
                ## Deletion ##
                $allele1 = $ref;
                $allele1 =~ s/$var//;
                $allele2 = "-";
                my $indel_size = length($allele1);
                $chr_start = $position + 1;
                $chr_stop = $chr_start + $indel_size - 1;
            }
            
            $observed =~ s/\//\t/g;
            $expected =~ s/\//\t/g;
            my $variant = join("\t", $chrom, $chr_start, $chr_stop, $allele1, $allele2, $observed, $expected, $chi_sum, $original_p);
            $indels{$key} .= "\n" if($indels{$key});
            $indels{$key} .= $variant;
            $self->stats->{'indels_in_vcf'}++;
        }
	}
    $fh->close;
	return \%indels;
}

sub _parse_sample_list {
    my $self = shift;
    my $sample_cns_list_file = $self->sample_cns_list_file;
	my %samples = ();
	my $line_counter = 0;
	my $fh = Genome::Sys->open_file_for_reading($sample_cns_list_file);
	
	while (my $line = <$fh>){
		chomp $line;
		$line_counter++;
	
		my ($sample, $cns_file) = split(/\t/, $line);
		
        $self->status_message("Loading CNS for $line_counter $sample...\n");
        $samples{$sample} = 1;
        my %indels = $self->_parse_cns_file($cns_file);

        foreach my $target_indel (keys %indels){
            my $sample_indel = join("\t", $sample, $target_indel);
            $self->target_indel_genotypes->{$sample_indel} = $indels{$target_indel};
        }			
	}
    
    $fh->close;
	return \%samples;
}

sub _parse_cns_file {
    my $self = shift;
    my $cns_file = shift;
	my %genotypes = ();
	my $line_counter = 0;
	my $fh = Genome::Sys->open_file_for_reading($cns_file);
	
	while (my $line = <$fh>) {
		chomp $line;
		$line_counter++;		
        next if $line_counter <= 1;

        my ($chrom, $position, $ref, $cns, $reads1, $reads2, $var_freq) = split(/\t/, $line);
        
        my $genotype = join(":", $cns, $reads1, $reads2, $var_freq);
        my $key = join("\t", $chrom, $position);
        
        if($self->target_indels->{$key}) {
            my @targetIndels = split(/\n/, $self->target_indels->{$key});
            
            foreach my $target_indel (@targetIndels) {
                my ($target_chrom, $target_chr_start, $target_chr_stop, $target_allele1, $target_allele2) = split(/\t/, $target_indel);
                
                ## If we have a reference call, record wild-type for every possible indel here ##	
                if($cns eq $ref) {
                    $genotype = join(":", "Ref", $cns, $reads1, $reads2, $var_freq);
                    $genotypes{$target_indel} = $genotype;
                }
                elsif(length($cns) == 1 && $cns ne "N") {
                    ## A SNP was called ##
                    $genotype = join(":", "SNP", $cns, $reads1, $reads2, $var_freq);
                    $genotypes{$target_indel} = $genotype;
                }
                elsif($cns =~ '\/') {
                    ## Determine type and size of target indel ##
                    my ($target_indel_type, $target_indel_size, $target_indel_bases) = "";
                    if($target_allele1 eq '-')
                    {
                        $target_indel_type = "ins";
                        $target_indel_bases = $target_allele2;
                        $target_indel_size = length($target_indel_bases);
                    }
                    else
                    {
                        $target_indel_type = "del";
                        $target_indel_bases = $target_allele1;
                        $target_indel_size = length($target_indel_bases);
                    }

                    ## Determine type and size of sample indel ##
                    
                    my ($var1, $var2) = split(/\//, $cns);
                    my $indel_type = '';
                    my $indel_size = '';
                    my $indel_bases = '';
                    my $indel_zygosity = '';

                    if(substr($var2, 0, 1) eq '+') {
                        $indel_type = "ins";
                        $indel_bases = substr($var2, 1, 99);
                        $indel_size = length($indel_bases);
                        $indel_zygosity = "Het";
                        $indel_zygosity = "Hom" if($var1 eq $var2);
                    }
                    elsif(substr($var2, 0, 1) eq '-') {
                        $indel_type = "del";
                        $indel_bases = substr($var2, 1, 99);
                        $indel_size = length($indel_bases);
                        $indel_zygosity = "Het";
                        $indel_zygosity = "Hom" if($var1 eq $var2);
                    }

                    if($indel_type) {
                        ## Same type and size, so record it ##
                        if($target_indel_type eq $indel_type && $target_indel_size == $indel_size) {
                            $genotype = join(":", $indel_zygosity, $cns, $reads1, $reads2, $var_freq);
                            $genotypes{$target_indel} = $genotype;
                        }
                        else {
                            $self->warning_message("Expected $target_indel_type $target_indel_size but got $indel_type $indel_size from $cns\n");
                            ## A different indel was called, so we should report a wildtype for this target indel ##
                            my $genotype = join(":", "Other", $cns, $reads1, $reads2, $var_freq);
                            $genotypes{$target_indel} = $genotype;
                        }
                    }
                    ## else: Not sure what this is, so don't report sample genotype ##
                }
            }
        }
	}
    $fh->close;
	return %genotypes;
}

sub _get_hwe_p {
    my $self = shift;
    my $obs_hets = shift;
    my $obs_hom1 = shift;
    my $obs_hom2 = shift;

    return -1 if($obs_hom1 < 0 || $obs_hom2 < 0 || $obs_hets <0);

    # rare homozygotes
    my $obs_homr;

    # common homozygotes
    my $obs_homc;
    if($obs_hom1 < $obs_hom2) {
        $obs_homr = $obs_hom1;
        $obs_homc = $obs_hom2;
    } else {
        $obs_homr = $obs_hom2;
        $obs_homc = $obs_hom1;
    }

    # number of rare allele copies
    my $rare_copies = 2 * $obs_homr + $obs_hets;

    # total number of genotypes
    my $genotypes = $obs_homr + $obs_homc + $obs_hets;

    return -1 if($genotypes <= 0);

    # Initialize probability array
    my @het_probs;
    for(my $i=0; $i<=$rare_copies; $i++) {
        $het_probs[$i] = 0.0;
    }

    # start at midpoint
    my $mid = int($rare_copies * (2 * $genotypes - $rare_copies) / (2 * $genotypes));

    # check to ensure that midpoint and rare alleles have same parity
    if(($rare_copies & 1) ^ ($mid & 1)) {
        $mid++;
    }

    my $curr_hets = $mid;
    my $curr_homr = ($rare_copies - $mid) / 2;
    my $curr_homc = $genotypes - $curr_hets - $curr_homr;

    $het_probs[$mid] = 1.0;
    my $sum = $het_probs[$mid];
    for($curr_hets = $mid; $curr_hets > 1; $curr_hets -= 2) {
        $het_probs[$curr_hets - 2] = $het_probs[$curr_hets] * $curr_hets * ($curr_hets - 1.0) / (4.0 * ($curr_homr + 1.0) * ($curr_homc + 1.0));
        $sum += $het_probs[$curr_hets - 2];

        # 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
        $curr_homr++;
        $curr_homc++;
    }

    $curr_hets = $mid;
    $curr_homr = ($rare_copies - $mid) / 2;
    $curr_homc = $genotypes - $curr_hets - $curr_homr;
    for($curr_hets = $mid; $curr_hets <= $rare_copies - 2; $curr_hets += 2) {
        $het_probs[$curr_hets + 2] = $het_probs[$curr_hets] * 4.0 * $curr_homr * $curr_homc / (($curr_hets + 2.0) * ($curr_hets + 1.0));
        $sum += $het_probs[$curr_hets + 2];

        # add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
        $curr_homr--;
        $curr_homc--;
    }

    for(my $i=0; $i<=$rare_copies; $i++) {
        $het_probs[$i] /= $sum;
    }

    # Initialise P-value
    my $p_hwe = 0.0;

    # P-value calculation for p_hwe
    for(my $i = 0; $i <= $rare_copies; $i++) {
        if($het_probs[$i] > $het_probs[$obs_hets]) {
            next;
        }
        $p_hwe += $het_probs[$i];
    }

    if($p_hwe > 1) {
        $p_hwe = 1.0;
    }
    return $p_hwe;
}

1;
