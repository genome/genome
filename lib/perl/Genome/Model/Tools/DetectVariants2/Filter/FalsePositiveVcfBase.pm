package Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcfBase;

use warnings;
use strict;

use Genome;
use Workflow;
use Workflow::Simple;
use Carp;
use Data::Dumper;
use Genome::Utility::Vcf ('open_vcf_file', 'parse_vcf_line', 'deparse_vcf_line', 'get_vcf_header', 'get_samples_from_header');


class Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcfBase {
    is => 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositive',
    is_abstract => 1,

    has_input => [
        _filters => {
            is => 'HashRef',
            is_optional => 1,
            doc => 'The filter names and descriptions',
        },
    ],
};

sub filter_name { 'FalsePositiveVcf' };

sub _get_readcount_line {
    my $self = shift;
    my ($readcount_fh,$chr,$pos) = @_;
    while( my $line = $readcount_fh->getline){
        chomp $line;
        my ($rc_chr,$rc_pos) = split "\t", $line;
        if(($chr eq $rc_chr)&&($pos == $rc_pos)){
            return $line;
        }
    }
    return;
}

# This method scans the lines of the readcount file until the matching line is found
sub make_buffered_rc_searcher {
    my $self = shift;
    my $fh = shift;
    unless($fh) {
        croak "Can't create buffered search from undefined file handle\n";
    }
    my $current_line;
    return sub {
        my ($chr, $pos) = @_;

        unless($current_line) {
            if($current_line = $fh->getline) {
                chomp $current_line;
            }
            else {
                return;
            }
        }

        my ($rc_chr,$rc_pos) = split /\t/, $current_line;
        if(($chr eq $rc_chr)&&($pos == $rc_pos)){
            my $temp = $current_line;
            $current_line = undef;
            return $temp;
        }
        else {
            #our line should be coming later
            return;
        }
    }
}

sub add_filter_to_vcf_header {
    my ($self, $parsed_header, $filter_name, $filter_description) = @_;
    my $column_header = pop @$parsed_header;
    my $filter_line = qq{##FILTER=<ID=$filter_name,Description="$filter_description">\n};
    push @$parsed_header, $filter_line, $column_header;
}

sub set_format_field {
    my ($self, $parsed_line, $sample, $format_tag, $format_value, %params) = @_;
    if(!exists($parsed_line->{sample}{$sample}{$format_tag})) {
        push @{$parsed_line->{'_format_fields'}}, $format_tag;
        for my $sample (keys %{$parsed_line->{sample}}) {
            #initialize the new format tag
            $parsed_line->{sample}{$sample}{$format_tag} = ".";
        }
    }
    if( $params{is_filter_fail} && $parsed_line->{sample}{$sample}{$format_tag} eq 'PASS') {
        $parsed_line->{sample}{$sample}{$format_tag} = $format_value;   #should really do type checking...
    }
    elsif( $params{append} && $parsed_line->{sample}{$sample}{$format_tag} ne "." ) {
        $parsed_line->{sample}{$sample}{$format_tag} .= ";$format_value";
    }
    else {
        $parsed_line->{sample}{$sample}{$format_tag} = $format_value;   #should really do type checking...
    }
}

# Find the variant allele for this sample (look at the GT and ALT)
#FIXME This is old hacky logic from the old FalsePositive filter and should probably not be used in a base class
sub get_variant_for_sample {
    my $self = shift;
    my $alt = shift;
    my $gt = shift;
    my $reference = shift;

    unless (defined $alt and defined $gt) {
        die $self->error_message("Either alt or gt are undefined. Alt: $alt, GT:$gt");
    }

    if($gt eq '.') {
        #there is no genotype here and thus no variant!
        return;
    }

    my @alts = split(",",$alt);
    my @gt = split("/", $gt);

    # Make sure the gt is either all numbers (standard) or all alleles (denovo polymutt)
    my $is_numeric = 0;
    my $is_alpha = 0;
    for my $gt_value (@gt) {
        if ($gt_value =~ m/\d/) {
            $is_numeric = 1;
        } elsif ($gt_value =~ /[ACGT]/) {
            $is_alpha = 1;
        } else {
            die $self->error_message("GT values do not appear to be either numeric or A/C/T/G! GT: " . join(",", @$gt));
        }
    }

    unless ($is_numeric xor $is_alpha) {
        die $self->error_message("GT values do not appear to be all numeric or all alpha! GT: " . join(",", @$gt));
    }

    # FIXME GT here might be actual alleles instead of numbers if this comes from polymutt denovo... not ideal but deal with it for now
    if ($is_numeric) {
        @gt = $self->convert_numeric_gt_to_alleles(\@alts, \@gt, $reference); # TODO
    }

    my %nonref_alleles;
    for my $allele (@gt) {
        unless ($allele eq $reference) {
            $nonref_alleles{$allele} = 1;
        }
    }
    my @nonref_alleles = keys %nonref_alleles;

    # If we got no nonref calls, no variants!
    unless (@nonref_alleles) {
        return;
    }

# FIXME these cases are actually ok, because the point of denovo not using GT numbers is that the denovo sites are not segregating sites
# Once we clean up the denovo file, remove all this hacky stuff
=cut
    # If there is only one alt, not much work in determining our variant
    if (scalar (@alts) == 1) {
        unless (scalar @nonref_alleles == 1) {
            die $self->error_message("Found only one alt value but more than one nonref_allele unique value! Alt: $alt GT: $gt");
        }
        unless ($alts[0] eq $nonref_alleles[0]) {
            #die $self->error_message("Single alt and single GT allele found, but they do not match! Alt: $alt GT: $gt");
        }
        return $alt;
    }
=cut

    # If we have only one nonref allele present, return that
    if (scalar @nonref_alleles == 1) {
        return $nonref_alleles[0];
    }

    # If we have more than one nonref allele, break the tie and return a single allele (this is the part that is hacky) #FIXME
    my $variant = $self->prioritize_allele(\@nonref_alleles);

    unless ($variant) {
        die $self->error_message("Failed to get a variant from prioritize_alleles for alleles: " . join(",",@nonref_alleles));
    }

    return $variant;
}

# Convert numeric gt values to A/C/T/G
sub convert_numeric_gt_to_alleles {
    my $self = shift;
    my $alt = shift;
    my $gt = shift;
    my $reference = shift;

    # Put reference in the alt list so that GT 0 (index 0) points to reference, 1 points to the first alt, etc)
    my @alt_including_ref = ($reference, @$alt);

    # Make sure the GT makes sense and also find the unique non-ref alleles
    my %alleles;
    for my $gt_number (@$gt) {
        unless ($gt_number <= scalar(@alt_including_ref) ) {
            die $self->error_message("GT does not make sense for the number of alts present. GT: $gt Alt: $alt");
        }

        my $current_allele = $alt_including_ref[$gt_number];
        $alleles{$current_allele} = 1;
    }

    my @alleles = keys %alleles;

    return @alleles;
}

# Break the tie between more than one variant allele
# FIXME This is arbitrary. A better way would be to test both and pass/fail each allele OR pass if either is ok
sub prioritize_allele {
    my $self = shift;
    my $alleles = shift;

    # Prioritization from old iupac_to_base is:  A > C > G > T ... so we can just alpha sort and everything is fine
    my @sorted_alleles = sort (@$alleles);

    return $sorted_alleles[0];
}

sub generate_filter_names {
    my $self = shift;
    my %filters;
    $filters{'position'} = [sprintf("PB%0.f",$self->min_read_pos*100), "Average position on read less than " . $self->min_read_pos . " or greater than " . (1 - $self->min_read_pos) . " fraction of the read length"];
    $filters{'strand_bias'} = [sprintf("SB%0.f",$self->min_strandedness*100), "Reads supporting the variant have less than " . $self->min_strandedness . " fraction of the reads on one strand, but reference supporting reads are not similarly biased"];
    $filters{'min_var_count'} = [ "MVC".$self->min_var_count, "Less than " . $self->min_var_count . " high quality reads support the variant"];
    $filters{'min_var_freq'} = [ sprintf("MVF%0.f",$self->min_var_freq*100), "Variant allele frequency is less than " . $self->min_var_freq];
    $filters{'mmqs_diff'} = [ sprintf("MMQSD%d",$self->max_mm_qualsum_diff), "Difference in average mismatch quality sum between variant and reference supporting reads is greater than " . $self->max_mm_qualsum_diff];
    $filters{'mq_diff'} = [ sprintf("MQD%d",$self->max_mapqual_diff), "Difference in average mapping quality sum between variant and reference supporting reads is greater than " . $self->max_mapqual_diff];
    $filters{'read_length'} = [ sprintf("RLD%d",$self->max_readlen_diff), "Difference in average clipped read length between variant and reference supporting reads is greater than " . $self->max_readlen_diff];
    $filters{'dist3'} = [ sprintf("DETP%0.f",$self->min_var_dist_3*100), "Average distance of the variant base to the effective 3' end is less than " . $self->min_var_dist_3];
    $filters{'homopolymer'} = [ sprintf("HPMR%d",$self->min_homopolymer), "Variant is flanked by a homopolymer of the same base and of length greater than or equal to " . $self->min_homopolymer];
    $filters{'var_mmqs'} = [ sprintf("MMQS%d",$self->max_var_mm_qualsum), "The average mismatch quality sum of reads supporting the variant is greater than " . $self->max_var_mm_qualsum] if($self->max_var_mm_qualsum);
    $filters{'no_var_readcount'} = [ "NRC", "Unable to grab readcounts for variant allele"];
    $filters{'incomplete_readcount'} = [ "IRC", "Unable to grab any sort of readcount for either the reference or the variant allele"];
    $self->_filters(\%filters);
}

# Print out stats from a hashref
sub print_stats {
    my $self = shift;
    my $stats = shift;
    print $stats->{'num_variants'} . " variants\n";
    print $stats->{'num_MT_sites_autopassed'} . " MT sites were auto-passed\n";
    print $stats->{'num_random_sites_autopassed'} . " chrN_random sites were auto-passed\n" if($stats->{'num_random_sites_autopassed'});
    print $stats->{'num_no_allele'} . " failed to determine variant allele\n";
    print $stats->{'num_no_readcounts'} . " failed to get readcounts for variant allele\n";
    print $stats->{'num_fail_pos'} . " had read position < " , $self->min_read_pos."\n";
    print $stats->{'num_fail_strand'} . " had strandedness < " . $self->min_strandedness . "\n";
    print $stats->{'num_fail_varcount'} . " had var_count < " . $self->min_var_count. "\n";
    print $stats->{'num_fail_varfreq'} . " had var_freq < " . $self->min_var_freq . "\n";

    print $stats->{'num_fail_mmqs'} . " had mismatch qualsum difference > " . $self->max_mm_qualsum_diff . "\n";
    print $stats->{'num_fail_var_mmqs'} . " had variant MMQS > " . $self->max_var_mm_qualsum . "\n" if($stats->{'num_fail_var_mmqs'});
    print $stats->{'num_fail_mapqual'} . " had mapping quality difference > " . $self->max_mapqual_diff . "\n";
    print $stats->{'num_fail_readlen'} . " had read length difference > " . $self->max_readlen_diff . "\n";
    print $stats->{'num_fail_dist3'} . " had var_distance_to_3' < " . $self->min_var_dist_3 . "\n";
    print $stats->{'num_fail_homopolymer'} . " were in a homopolymer of " . $self->min_homopolymer . " or more bases\n";

    print $stats->{'num_pass_filter'} . " passed the strand filter\n";

    return 1;
}

# Given an input vcf and an output filename, print a region_list from all variants in the input vcf
sub print_region_list {
    my $self = shift;
    my $input_file = shift;
    my $output_file = shift;
    my $output_fh = Genome::Sys->open_file_for_writing($output_file);

    # Strip out the header
    my ($input_fh, $header) = $self->parse_vcf_header($input_file);

    ## Print each line to file in order to get readcounts
    $self->status_message("Printing variants to temporary region_list file $output_file...");
    while(my $line = $input_fh->getline) {
        if ($self->should_print_region_line($line)) {
            my ($chr, $pos) = split /\t/, $line;
            print $output_fh "$chr\t$pos\t$pos\n";
        }
    }

    $output_fh->close;
    $input_fh->close;

    return 1;
}

# Given an input vcf file name, return the file handle at the position of the first variant line, and the header lines as an array
#FIXME probably move this to a base class
sub parse_vcf_header {
    my $self = shift;
    my $input_file = shift;
    my $input_fh = $self->open_input_file($input_file);

    my @header;
    my $header_end = 0;
    while (!$header_end) {
        my $line = $input_fh->getline;
        if ($line =~ m/^##/) {
            push @header, $line;
        } elsif ($line =~ m/^#/) {
            push @header, $line;
            $header_end = 1;
        } else {
            die $self->error_message("Missed the final header line with the sample list? Last line examined: $line Header so far: " . join("\n", @header));
        }
    }

    return ($input_fh, \@header);
}

# Info fields this filter requires, override in each filter
sub required_info_fields {
    return;
}

sub get_readcounts_from_vcf_line {
    my ($self, $parsed_vcf_line, $readcount_searcher_by_sample, $sample_name) = @_;

    my $readcount_searcher = $readcount_searcher_by_sample->{$sample_name};
    unless ($readcount_searcher) {
        die $self->error_message("Could not get readcount searcher for sample $sample_name " . Data::Dumper::Dumper $readcount_searcher_by_sample);
    }

    #want to do this BEFORE anything else so we don't screw up the file parsing...
    my $readcounts;
    # FIXME this should get readcount lines until we have the correct chrom and pos matching our input line
    # FIXME no readcount output line is a valid output condition if there are no reads at all covering the position. Maybe bam-readcount should be fixed to not do this...
    $readcounts = $readcount_searcher->($parsed_vcf_line->{chromosome},$parsed_vcf_line->{position});
    unless($readcounts) {
        #no data at this site, set FT to null
        #FIXME This breaks what little encapsulation we have started...
        if(!exists($parsed_vcf_line->{sample}{$sample_name}{$self->filter_sample_format_tag})) {
            $self->set_format_field($parsed_vcf_line, $sample_name,
                $self->filter_sample_format_tag, ".");
        }
        return;
    }

    return $readcounts;
}

# Given a vcf line structure and a sample/readcount file, determine if that sample should be filtered
sub filter_one_sample {
    my ($self, $parsed_vcf_line, $readcount_searcher_by_sample, $stats, $sample_name, $var) = @_;

    my $readcounts = $self->get_readcounts_from_vcf_line($parsed_vcf_line,
        $readcount_searcher_by_sample, $sample_name);
    unless ($readcounts) {
        return;
    }

    my $chrom = $parsed_vcf_line->{chromosome};
    my $pos = $parsed_vcf_line->{position};

    my $min_read_pos = $self->min_read_pos;
    my $max_read_pos = 1 - $min_read_pos;
    my $min_var_freq = $self->min_var_freq;
    my $min_var_count = $self->min_var_count;

    my $min_strandedness = $self->min_strandedness;
    my $max_strandedness = 1 - $min_strandedness;

    my $max_mm_qualsum_diff = $self->max_mm_qualsum_diff;
    my $max_mapqual_diff = $self->max_mapqual_diff;
    my $max_readlen_diff = $self->max_readlen_diff;
    my $min_var_dist_3 = $self->min_var_dist_3;
    my $max_var_mm_qualsum = $self->max_var_mm_qualsum if($self->max_var_mm_qualsum);

    my $original_line = $parsed_vcf_line->{original_line};

    ## Run Readcounts ##
    # TODO var has some special logic from iubpac_to_base that prefers an allele in the case of het variants

    my $ref = uc($parsed_vcf_line->{reference});
    $var = $self->update_variant_for_sample($parsed_vcf_line,
        $sample_name, $var);
    unless(defined($var)) {
        #no variant reads here. Won't filter. Add FT tag of "PASS"
        $self->pass_sample($parsed_vcf_line, $sample_name);
        $stats->{'num_no_allele'}++;
        return;
    }

    # TODO Evaluate if this really helps improve performance. We've already sunk the readcounts at this point...
    ## Skip MT chromosome sites,w hich almost always pass ##
    if($chrom eq "MT" || $chrom eq "chrMT") {
        ## Auto-pass it to increase performance ##
        $stats->{'num_MT_sites_autopassed'}++;
        $self->pass_sample($parsed_vcf_line, $sample_name);
    }


    ## Parse the results for each allele ##
    my $ref_result = $self->read_counts_by_allele($readcounts, $ref);
    my $var_result = $self->read_counts_by_allele($readcounts, $var);

    if($ref_result && $var_result) {
        ## Parse out the bam-readcounts details for each allele. The fields should be: ##
        #num_reads : avg_mapqual : avg_basequal : avg_semq : reads_plus : reads_minus : avg_clip_read_pos : avg_mmqs : reads_q2 : avg_dist_to_q2 : avgRLclipped : avg_eff_3'_dist
        my ($ref_count, $ref_map_qual, $ref_base_qual, $ref_semq, $ref_plus, $ref_minus, $ref_pos, $ref_subs, $ref_mmqs, $ref_q2_reads, $ref_q2_dist, $ref_avg_rl, $ref_dist_3) = split(/\t/, $ref_result);
        my ($var_count, $var_map_qual, $var_base_qual, $var_semq, $var_plus, $var_minus, $var_pos, $var_subs, $var_mmqs, $var_q2_reads, $var_q2_dist, $var_avg_rl, $var_dist_3) = split(/\t/, $var_result);

        my $ref_strandedness = my $var_strandedness = 0.50;
        $ref_dist_3 = 0.5 if(!$ref_dist_3);

        ## Use conservative defaults if we can't get mismatch quality sums ##
        $ref_mmqs = 50 if(!$ref_mmqs);
        $var_mmqs = 0 if(!$var_mmqs);
        my $mismatch_qualsum_diff = $var_mmqs - $ref_mmqs;

        ## Determine map qual diff ##

        my $mapqual_diff = $ref_map_qual - $var_map_qual;


        ## Determine difference in average supporting read length ##

        my $readlen_diff = $ref_avg_rl - $var_avg_rl;


        ## Determine ref strandedness ##

        if(($ref_plus + $ref_minus) > 0) {
            $ref_strandedness = $ref_plus / ($ref_plus + $ref_minus);
            $ref_strandedness = sprintf("%.2f", $ref_strandedness);
        }

        ## Determine var strandedness ##

        if(($var_plus + $var_minus) > 0) {
            $var_strandedness = $var_plus / ($var_plus + $var_minus);
            $var_strandedness = sprintf("%.2f", $var_strandedness);
        }

        my $has_failed = 0;

        if($var_count && ($var_plus + $var_minus)) {
            ## We must obtain variant read counts to proceed ##

            my $var_freq = $var_count / ($ref_count + $var_count);

            ## FAILURE 1: READ POSITION ##
            if(($var_pos < $min_read_pos)) {
                $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{position}->[0]);
                print "$original_line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tReadPos<$min_read_pos\n"if ($self->verbose);
                $stats->{'num_fail_pos'}++;
                $has_failed = 1;
            }

            ## FAILURE 2: Variant is strand-specific but reference is NOT strand-specific ##
            if(($var_strandedness < $min_strandedness || $var_strandedness > $max_strandedness) && ($ref_strandedness >= $min_strandedness && $ref_strandedness <= $max_strandedness)) {
                $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{strand_bias}->[0]);
                ## Print failure to output file if desired ##
                print "$original_line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tStrandedness: Ref=$ref_strandedness Var=$var_strandedness\n"if ($self->verbose);
                $has_failed = 1;
                $stats->{'num_fail_strand'}++;
            }

            ## FAILURE : Variant allele count does not meet minimum ##
            if($var_count < $min_var_count) {
                $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{min_var_count}->[0] );
                print "$original_line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarCount:$var_count\n" if ($self->verbose);
                $has_failed = 1;
                $stats->{'num_fail_varcount'}++;
            }

            ## FAILURE : Variant allele frequency does not meet minimum ##
            if($var_freq < $min_var_freq) {
                $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{min_var_freq}->[0]);
                print "$original_line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarFreq:$var_freq\n" if ($self->verbose);
                $has_failed = 1;
                $stats->{'num_fail_varfreq'}++;
            }

            ## FAILURE 3: Paralog filter for sites where variant allele mismatch-quality-sum is significantly higher than reference allele mmqs
            if($mismatch_qualsum_diff> $max_mm_qualsum_diff) {
                ## Print failure to output file if desired ##
                $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{mmqs_diff}->[0]);
                print "$original_line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tMismatchQualsum:$var_mmqs-$ref_mmqs=$mismatch_qualsum_diff" if ($self->verbose);
                $has_failed = 1;
                $stats->{'num_fail_mmqs'}++;
            }

            ## FAILURE 4: Mapping quality difference exceeds allowable maximum ##
            if($mapqual_diff > $max_mapqual_diff) {
                $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{mq_diff}->[0]);
                print "$original_line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tMapQual:$ref_map_qual-$var_map_qual=$mapqual_diff" if ($self->verbose);
                $has_failed = 1;
                $stats->{'num_fail_mapqual'}++;
            }

            ## FAILURE 5: Read length difference exceeds allowable maximum ##
            if($readlen_diff > $max_readlen_diff) {
                $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{read_length}->[0]);
                print "$original_line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tReadLen:$ref_avg_rl-$var_avg_rl=$readlen_diff" if ($self->verbose);
                $has_failed = 1;
                $stats->{'num_fail_readlen'}++;
            }

            ## FAILURE 5: Read length difference exceeds allowable maximum ##
            if($var_dist_3 < $min_var_dist_3) {
                $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{dist3}->[0]);
                print "$original_line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarDist3:$var_dist_3\n" if ($self->verbose);
                $has_failed = 1;
                $stats->{'num_fail_dist3'}++;
            }

            if($self->fails_homopolymer_check($self->reference_sequence_input, $self->min_homopolymer, $chrom, $pos, $pos, $ref, $var)) {
                $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{homopolymer}->[0]);
                print "$original_line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tHomopolymer\n" if ($self->verbose);
                $has_failed = 1;
                $stats->{'num_fail_homopolymer'}++;
            }

            if($max_var_mm_qualsum && $var_mmqs > $max_var_mm_qualsum) {
                $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{var_mmqs}->[0]);
                print "$original_line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarMMQS: $var_mmqs > $max_var_mm_qualsum\n" if ($self->verbose);
                $has_failed = 1;
                $stats->{'num_fail_var_mmqs'}++;
            }

            ## SUCCESS: Pass Filter ##
            unless($has_failed) {
                $stats->{'num_pass_filter'}++;
                ## Print output, and append strandedness information ##
                $self->pass_sample($parsed_vcf_line, $sample_name);
                print "$original_line\t\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tPASS\n" if($self->verbose);
            }
        } else {
            $stats->{'num_no_readcounts'}++;
            $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{no_var_readcount}->[0]);
            print "$chrom\t$pos\t$pos\t$ref\t$var\tFAIL no reads in $var_result\n" if($self->verbose);
        }
    } else {
        # was unable to get read_counts_by_allele for both the reference and variant... just consider this low quality for lack of a better plan
        # This seemingly only happens when the reference allele is an IUB code, which appears to be rare (1 in a dataset of millions, and only some data sets)
        $stats->{'num_no_readcounts'}++;
        $self->fail_sample($parsed_vcf_line, $sample_name, $self->_filters->{incomplete_readcount}->[0]);
    }
}

# In the line provided, set the sample provided to passed (FT)
sub pass_sample {
    my $self = shift;
    my $parsed_vcf_line = shift;
    my $sample_name = shift;

    $self->status_message("Entering pass sample\n");
    $self->set_format_field($parsed_vcf_line, $sample_name,
        $self->filter_sample_format_tag, "PASS");

    return 1;
}

sub get_all_sample_names_and_bam_paths {
    my $self = shift;
    my @sample_names;
    my @bam_paths;

    if (!$self->aligned_reads_input) {
        for my $alignment_result ($self->alignment_results) {
            my $sample_name = $self->find_sample_name_for_alignment_result($alignment_result);
            push @sample_names, $sample_name;
            push @bam_paths, $alignment_result->merged_alignment_bam_path;
        }
    } else {
        push @bam_paths, $self->aligned_reads_input;
        push @sample_names, $self->aligned_reads_sample;
        if ($self->control_aligned_reads_input) {
            push @bam_paths, $self->control_aligned_reads_input;
            push @sample_names, $self->control_aligned_reads_sample;
        }
    }

    unless (scalar(@sample_names) == scalar(@bam_paths)) {
        die $self->error_message("Did not get the same number of sample names as bam paths.\n" . Data::Dumper::Dumper @sample_names . Data::Dumper::Dumper @bam_paths);
    }

    return (\@sample_names, \@bam_paths);
}

sub add_per_bam_params_to_input {
    my ($self, %inputs) = @_;

    my ($sample_names, $bam_paths) = $self->get_all_sample_names_and_bam_paths;

    for my $sample_name (@$sample_names) {
        my $bam_path = shift @$bam_paths;

        $self->status_message("Running BAM Readcounts for sample $sample_name...");
        my $readcount_file = $self->_temp_staging_directory . "/$sample_name.readcounts";  #this is suboptimal, but I want to wait until someone tells me a better way...multiple options exist
        if (-f $bam_path) {
            $inputs{"bam_${sample_name}"} = $bam_path;
        } else {
            die "merged_alignment_bam_path does not exist: $bam_path";
        }
        $inputs{"readcounts_${sample_name}"} = $readcount_file;
    }

    return %inputs;
}

sub generate_and_run_readcounts_in_parallel {
    my $self = shift;
    my $region_path = shift;

    my %inputs;

    #set up global readcount params
    $inputs{reference_fasta} = $self->reference_sequence_input;
    $inputs{region_list} = $region_path;
    $inputs{minimum_base_quality} = $self->bam_readcount_min_base_quality;
    $inputs{use_version} = $self->bam_readcount_version;

    my ($sample_names) = $self->get_all_sample_names_and_bam_paths;
    %inputs = $self->add_per_bam_params_to_input(%inputs);

    my $workflow = Workflow::Model->create(
        name=> "FalsePositiveVcf parallel readcount file creation",
        input_properties => [
        keys %inputs
        ],
        output_properties => [
        'output',
        ],
    );
    for my $sample (@$sample_names) {
        my $op = $workflow->add_operation(
            name=>"readcount creation for $sample",
            operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Sam::Readcount"),
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"reference_fasta",
            right_operation=>$op,
            right_property=>"reference_fasta",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"region_list",
            right_operation=>$op,
            right_property=>"region_list",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"use_version",
            right_operation=>$op,
            right_property=>"use_version",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"minimum_base_quality",
            right_operation=>$op,
            right_property=>"minimum_base_quality",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"bam_$sample",
            right_operation=>$op,
            right_property=>"bam_file",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"readcounts_$sample",
            right_operation=>$op,
            right_property=>"output_file",
        );
        $workflow->add_link(
            left_operation=>$op,
            left_property=>"output_file",
            right_operation=>$workflow->get_output_connector,
            right_property=>"output",
        );
    }

    my $log_dir = $self->output_directory;
    if(Workflow::Model->parent_workflow_log_dir) {
        $log_dir = Workflow::Model->parent_workflow_log_dir;
    }
    $workflow->log_dir($log_dir);

    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }
    $self->status_message("Now launching readcount generation jobs");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %inputs);
    unless($result) {
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("parallel readcount generation workflow did not return correctly.");
    }

    #all succeeded so open files
    my $readcount_searcher_by_sample;
    for my $sample (@$sample_names) {
        $readcount_searcher_by_sample->{$sample} = $self->make_buffered_rc_searcher(Genome::Sys->open_file_for_reading($inputs{"readcounts_$sample"}));
    }
    return $readcount_searcher_by_sample;
}

# In the line provided, set the sample provided to failed (FT)
sub fail_sample {
    my $self = shift;
    my $parsed_vcf_line = shift;
    my $sample_name = shift;
    my $filter_reason = shift;

    $self->status_message("Entering fail sample\n");
    $self->set_format_field($parsed_vcf_line, $sample_name,
        $self->filter_sample_format_tag, $filter_reason,
        append => 1, is_filter_fail => 1);
    return 1;
}

sub region_path {
    my $self = shift;
    return $self->output_directory."/regions";
}

##########################################################################################
# Capture filter for high-depth, lower-breadth datasets
# Contact: Dan Koboldt (dkoboldt@genome.wustl.edu)
##########################################################################################
sub _filter_variants {
    my $self = shift;

    unless ($self->bam_readcount_version) {
        die $self->error_message("Bam readcount version is not specified");
    }

    ## Determine the strandedness and read position thresholds ##

    ## Initialize counters ##
    my $stats;
    $stats->{'num_variants'}  = $stats->{'num_no_readcounts'} = $stats->{'num_pass_filter'} = $stats->{'num_no_allele'} = 0;
    $stats->{'num_fail_varcount'} = $stats->{'num_fail_varfreq'} = $stats->{'num_fail_strand'} = $stats->{'num_fail_pos'} = $stats->{'num_fail_mmqs'} = $stats->{'num_fail_mapqual'} = $stats->{'num_fail_readlen'} = $stats->{'num_fail_dist3'} = 0;
    $stats->{'num_MT_sites_autopassed'} = $stats->{'num_fail_homopolymer'} = 0;

    #First, need to create a variant list file to use for generating the readcounts.
    # FIXME this will work after the polymuttdenovo filter, but not directly after polymutt due to the separate denovo and standard filenames
    my $input_file = $self->input_file_path;

    $self->print_region_list($input_file, $self->region_path);

    if ($self->should_skip_filter) {
        $self->status_message("No variants found, copying file over and reporting success (this filter was a no-op");
        Genome::Sys->copy_file($input_file, $self->output_file_path);
        return 1;
    }
    ## Run BAM readcounts in batch mode to get read counts for all positions in file ##
    my $readcount_searcher_by_sample = $self->generate_and_run_readcounts_in_parallel($self->region_path);

    # Advance to the first variant line and print out the header
    my ($input_fh, $header) = $self->parse_vcf_header($input_file);
    #check here to see if header has FT format tag
    unless($self->header_already_added($header)) {
        my $col_header = $header->[-1];
        $header->[-1] = $self->filter_status_header;
        push @$header, $col_header;
    }
    #here embed the filter codes. This should be more centralized to be less craptastic
    $self->generate_filter_names;
    for my $filter (values %{$self->_filters}) {
        my $filter_name = $filter->[0];
        unless(grep {/$filter_name/} @$header) {
            $self->add_filter_to_vcf_header($header,@$filter);
        }
    }

    my $output_fh = Genome::Sys->open_gzip_file_for_writing(
        $self->output_file_path);
    $output_fh->print(join("",@$header));
    my @sample_names = get_samples_from_header($header);

    ## Parse the variants file ##
    while(my $line = $input_fh->getline) {
        # FIXME should just pass in a single sample here instead of a whole line. Or a sample joined with a line to make a whole single sample vcf line?
        my $parsed_line = parse_vcf_line($line, \@sample_names);
        $self->filter_one_line($parsed_line, $readcount_searcher_by_sample, $stats, \@sample_names);
        $output_fh->print(deparse_vcf_line($parsed_line,\@sample_names));
    }

    $input_fh->close;
    $output_fh->close;
    $self->print_stats($stats);

    $self->_convert_to_standard_formats($self->output_file_path);

    return 1;
}

#############################################################
# Read_Counts_By_Allele - parse out readcount info for an allele
#
#############################################################

sub fails_homopolymer_check {
    my ($self, $reference, $min_homopolymer, $chrom, $start, $stop, $ref, $var) = @_;

    ## Auto-pass large indels ##

    my $indel_size = length($ref);
    $indel_size = length($var) if(length($var) > $indel_size);

    return(0) if($indel_size > 2);

    ## Build strings of homopolymer bases ##
    my $homoRef = $ref x $min_homopolymer;
    my $homoVar = $var x $min_homopolymer;

    ## Build a query string for the homopolymer check ##

    my $query_string = "";

    $query_string = $chrom . ":" . ($start - $min_homopolymer) . "-" . ($stop + $min_homopolymer);

    my $samtools_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
    my $sequence = `$samtools_path faidx $reference $query_string | grep -v \">\"`;
    chomp($sequence);

    if($sequence) {
        if($sequence =~ $homoVar) { #$sequence =~ $homoRef || {
            print join("\t", $chrom, $start, $stop, $ref, $var, "Homopolymer: $sequence") . "\n" if($self->verbose);
            return($sequence);
        }
    }

    return(0);
}

1;
