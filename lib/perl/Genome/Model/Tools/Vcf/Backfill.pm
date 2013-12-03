package Genome::Model::Tools::Vcf::Backfill;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;
use Genome::Utility::Vcf "open_vcf_file";
use Sort::Naturally;

class Genome::Model::Tools::Vcf::Backfill{
    is => 'Command',
    has_input => [
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => "Output backfilled VCF",
        },
        merged_positions_file => {
            is => 'Text',
            doc => "Merged input file for this sample. Contains all positions and their total possible ALT alleles.",
        },
        pileup_file => {
            is => 'Text',
            doc => "Input pileup file for this sample",
        },
        vcf_file => {
            is => 'Text',
            doc => "Input vcf file for this sample to be backfilled",
        },
        use_bgzip => {
            is => 'Boolean',
            doc => 'Expect pileup (but not VCF) input in bgzip format and bgzips the output',
            default => 0,
        },
    ],
    doc => "Backfill a single sample VCF with reference information from pileup",
};


sub help_synopsis {
    <<'HELP';
Backfill a single sample VCF with reference information from pileup,
HELP
}

sub execute {
    my $self = shift;

    my ($pileup_fh, $output_fh);
    if ($self->use_bgzip) {
        $pileup_fh = Genome::Sys->open_gzip_file_for_reading($self->pileup_file);
        $output_fh = Genome::Sys->open_gzip_file_for_writing($self->output_file);
    } else {
        $pileup_fh = Genome::Sys->open_file_for_reading($self->pileup_file);
        $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    }
    my $vcf_fh = open_vcf_file($self->vcf_file);
    my $merged_fh = open_vcf_file($self->merged_positions_file);

    # Copy the header from the input vcf to the output vcf
    my $header_copied = 0;
    while (!$header_copied) {
        my $header_line = $vcf_fh->getline;
        if($header_line=~ m/^##/){
            $output_fh->print($header_line);
        } elsif( $header_line =~ m/^#CHROM/){
            # Add in a filter line to represent a backfilling error
            my $backfill_filter = '##FILTER=<ID=ForcedGenotype,Description="Forced genotyping for an uncalled variant">';
            $output_fh->print("$backfill_filter\n");
            $output_fh->print($header_line);
            $header_copied = 1;
        }
    }

    # Loop through both files and interleave the lines into the output
    my $pileup_line = $self->get_next_valid_pileup_line($pileup_fh);
    my $vcf_line = $vcf_fh->getline;
    while( $vcf_line && $pileup_line ){
        chomp $vcf_line;
        chomp $pileup_line;
        my $comparison = $self->compare_lines($pileup_line,$vcf_line);

        # The vcf file is ahead. Formulate and print VCF lines for the pileup line and advance the pileup filehandle 
        if($comparison == -1){
            my $merged_alts = $self->get_alts_for_position($merged_fh, $pileup_line);
            my $new_vcf_line = $self->create_vcf_line_from_pileup($pileup_line, $merged_alts);
            if ($new_vcf_line) {
                $output_fh->print("$new_vcf_line\n");
            }
            $pileup_line = $self->get_next_valid_pileup_line($pileup_fh);
        # The pileup file is ahead. Print the VCF line to output and advance the vcf filehandle 
        } elsif ($comparison == 1) {
            my $merged_alts = $self->get_alts_for_position($merged_fh, $vcf_line);
            $output_fh->print("$vcf_line\n");
            $vcf_line = $vcf_fh->getline;
        # They are both at the same position. Add pileup information for alts not already known and print.
        } else {
            my $merged_alts = $self->get_alts_for_position($merged_fh, $pileup_line);
            $vcf_line = $self->add_alt_information_to_vcf_line($vcf_line, $pileup_line, $merged_alts);
            $output_fh->print("$vcf_line\n");
            $pileup_line = $self->get_next_valid_pileup_line($pileup_fh);
            $vcf_line = $vcf_fh->getline;
        }
    }

    # Finish off whichever filehandle is not exhausted (including the single line left from the previous loop)
    while ($vcf_line) {
        chomp $vcf_line;
        my $merged_alts = $self->get_alts_for_position($merged_fh, $vcf_line);
        $output_fh->print("$vcf_line\n");
        $vcf_line = $vcf_fh->getline
    }
    while ($pileup_line) {
        chomp $pileup_line;
        my $merged_alts = $self->get_alts_for_position($merged_fh, $pileup_line);
        my $new_vcf_line = $self->create_vcf_line_from_pileup($pileup_line, $merged_alts);
        $output_fh->print("$new_vcf_line\n");
        $pileup_line = $self->get_next_valid_pileup_line($pileup_fh);
    }

    return 1;
}

# Advance the merged positions filehandle, make sure it is on the same chrom/position as the other line(s), and return the alternate allele string
sub get_alts_for_position {
    my $self = shift;
    my $fh = shift;
    my $current_line = shift;

    my ($current_chrom, $current_pos) = split "\t", $current_line;
    my $chrom = "";
    my $pos = "";
    my ($stop, $alts);

    while ($chrom ne $current_chrom || $pos != $current_pos) {
        my $line = $fh->getline;
        unless ($line) {
            die $self->error_message("Could not find a line in the merged positions/alts file to match the current line.\nLine:$current_line");
        }
        chomp $line;
        ($chrom, $pos, $stop, $alts) = split "\t", $line;
    }

    return $alts;
}

# Get the next pileup line that is not an indel
sub get_next_valid_pileup_line {
    my $self = shift;
    my $fh = shift;

    while (my $pileup_line = $fh->getline) {
        chomp $pileup_line;
        my ($chr, $pos, $ref, $genotype, $gq, $vaq, $mq, $dp, $read_bases, $base_qualities) = split("\t", $pileup_line);

        if ($ref ne "*") {
            return $pileup_line;
        }
    }

    return;
}

# Given a pileup line, return a vcf line containing information about the position
sub create_vcf_line_from_pileup {
    my $self = shift;
    my $pileup_line = shift;
    my $alt_string = shift;

    chomp $pileup_line;
    my ($chr, $pos, $ref, $genotype, $gq, $vaq, $mq, $dp, $read_bases, $base_qualities) = split("\t", $pileup_line);

    # Ignore insertions and deletions
    if ($ref eq "*") {
        return;
    }

    my @tags = qw(GT GQ DP BQ MQ AD);
    my %sample;

    # If the genotype differs from the reference, force genotype since no variant was called here
    if ($genotype eq $ref) {
        $sample{'GT'}= "0/0";
    } elsif ($genotype eq "N") {
        $sample{'GT'}= ".";
    } else {
        my @new_alts = Genome::Info::IUB->variant_alleles_for_iub($ref, $genotype);
        my @previous_alts = split ",", $alt_string;

        # Add new alts to existing alt string, if none are duplicates, without disrupting the order
        for my $new_alt (@new_alts) {
            my $found = 0;
            for my $old_alt (@previous_alts) {
                if ($new_alt eq $old_alt) {
                    $found = 1;
                }
            }
            unless ($found) {
                push (@previous_alts, $new_alt);
            }
        }

        # Calculate new GT
        my @genotype_alleles = Genome::Info::IUB->iub_to_alleles($genotype);
        my @gt;
        for my $genotype_allele (@genotype_alleles) {
            if ($genotype_allele eq $ref) {
                push @gt, 0;
            } else {
                my $alt_number = 1;
                for my $alt (@previous_alts) {
                    if ($genotype_allele eq $alt) {
                        push @gt, $alt_number;
                        last;
                    }
                    $alt_number++;
                }
            }
        }
        
        $sample{'GT'}= join( "/", sort(@gt) );
        $sample{'FT'}= "ForcedGenotype";
        push @tags, "FT";
        $alt_string = join ",", @previous_alts;
    }

    my ($ad_string, $bq_string) = $self->calculate_stats_for_alts($ref, $alt_string, $read_bases, $base_qualities);

    $sample{'GQ'}= $gq;
    $sample{'BQ'}= $bq_string;
    $sample{'MQ'}= $mq;
    $sample{'DP'}= $dp;
    $sample{'AD'}= $ad_string;

    my $format_tags = join ":", @tags;
    my @output;
    for my $tag (@tags) {
        my $value = $sample{$tag};
        unless (defined $value) {
            $value = "."
        }
        push @output, $value;
    }
    my $sample_string = join ":", @output;

    my $id = ".";
    my $qual = ".";
    my $info = ".";
    my $filter = "."; 
    my $new_vcf_line = join "\t", ($chr, $pos, $id, $ref, $alt_string, $qual, $filter, $info, $format_tags, $sample_string); 

    return $new_vcf_line;
}

# Given a vcf line and a pileup line 
# Return a vcf line that merges in the pileup information for AD and BQ that were not previously present in the line
sub add_alt_information_to_vcf_line {
    my $self = shift;
    my $vcf_line = shift;
    my $pileup_line = shift;
    my $new_alt = shift;

    # Gather values from both inputs
    my ($chrom, $pos, $id, $ref, $old_alt, $qual, $filter, $info, $format, $sample) = split "\t", $vcf_line;
    my ($pileup_chr, $pileup_pos, $pileup_ref, $pileup_genotype, $pileup_gq, $pileup_vaq, $pileup_mq, $pileup_dp, $read_bases, $base_qualities) = split("\t", $pileup_line);
    my @format_fields = split ":", $format;
    my @sample_fields = split ":", $sample;
    unless (scalar(@format_fields) == scalar(@sample_fields) ) {
        die $self->error_message("Format field (format) and sample ($sample) field do not have the same number of values");
    }
    my %sample_values;
    for my $format_tag (@format_fields) {
        $sample_values{$format_tag} = shift @sample_fields;
    }
    unless ($sample_values{"GT"}) {
        die $self->error_message("Could not find GT values in the line $vcf_line");
    }

    # compare the vcf gt to the vcf alt string, find out what alts we have information for and which we do not
    my @ad_values = split ",", $sample_values{"AD"};
    my @bq_values = split ",", $sample_values{"BQ"};
    my @old_alt_values = split ",", $old_alt;
    my @new_alt_values = split ",", $new_alt;

    # Find a unique list of GT values... but exclude 0. If the original GT is null we have no pre-existing information to consider
    my (%ad_per_alt, %bq_per_alt);

    # Shift the first ad and bq value off. These are for the reference.
    my $ref_ad = shift @ad_values;
    my $ref_bq = shift @bq_values;

    unless ($sample_values{"GT"} =~ /^(\.|\.\/\.)$/) {  # old ., new ./.
        my @gt_values = sort (split "/", $sample_values{"GT"});

        # This assumption may change at some point but for now we assume every VCF will have two values in the GT
        unless (scalar (@gt_values) >= 2) {
            die $self->error_message("This GT does not appear to have at least two values: " . $sample_values{"GT"});
        }

        my %unique_gt_values;
        for my $gt_value (@gt_values) {
            unless ($gt_value eq 0) {
                $unique_gt_values{$gt_value} = 1;
            }
        }
        @gt_values = keys %unique_gt_values;

        # Make sure the number of values we have for everything checks out
        unless ( (scalar(@ad_values) == scalar (@bq_values)) && (scalar(@ad_values) == scalar (@gt_values)) ){
            my $error_message = "For this vcf line: $vcf_line\n" . 
                "Have differing numbers of values for AD " . $sample_values{"AD"} . " and BQ " . $sample_values{"BQ"}. " and GT " . $sample_values{"GT"} . "\n" . 
                "There should be one AD and BQ value for each alt PLUS the reference. If this is not true, rebuilding your model group will fix this.";
            die $self->error_message($error_message);
        }

        # If we have information for all alts present, return the line as is
        if (scalar(@new_alt_values) == scalar(@ad_values) && scalar(@new_alt_values) == scalar(@bq_values) ) {
            #$self->status_message("Information is already present for all alts in $vcf_line");
            return $vcf_line;
        }

        # Gather AD and BQ values already present in the line
        for my $gt (@gt_values ) {
            my $alt_represented = $old_alt_values[$gt-1]; # GT will be a 1 based index of the alt allele values, since 0 is reference
            unless ($alt_represented) {
                die $self->error_message("Could not find an alt to represent gt value $gt");
            }
            # These should be sorted by GT number already, so shifting should be safe
            $ad_per_alt{$alt_represented} = shift @ad_values;
            $bq_per_alt{$alt_represented} = shift @bq_values;
        }
    }

    # Generate AD and BQ information for all alts for which we have no information
    my ($ad_string, $bq_string) = $self->calculate_stats_for_alts($ref, $new_alt, $read_bases, $base_qualities);
    my @backfilled_ad_values = split ",", $ad_string;
    my @backfilled_bq_values = split ",", $bq_string;

    # Combine AD and BQ information already in the line with information from pileup so we have info for every alt
    for my $alt (@new_alt_values) {
        my $backfilled_ad = shift @backfilled_ad_values;
        my $backfilled_bq = shift @backfilled_bq_values;
        # Set AD and BQ from the pileup information if not previously present in the vcf
        unless ($ad_per_alt{$alt} && $bq_per_alt{$alt}) {
            $ad_per_alt{$alt} = $backfilled_ad;
            $bq_per_alt{$alt} = $backfilled_bq;
        }
    }
    # Join BQ and AD in the new alt order
    my (@new_ad_values, @new_bq_values);
    # Refs first...
    push @new_ad_values, $ref_ad;
    push @new_bq_values, $ref_bq;
    for my $alt (@new_alt_values) {
        push @new_ad_values, $ad_per_alt{$alt};
        push @new_bq_values, $bq_per_alt{$alt};
    }
    $sample_values{"AD"} = join(",", @new_ad_values);
    $sample_values{"BQ"} = join(",", @new_bq_values);

    # Get the new GT value since the alt string may be different
    $sample_values{"GT"} = $self->regenerate_gt($ref, $old_alt, $sample_values{"GT"}, $new_alt);

    # Reconstruct the sample field with new information
    my @final_sample_values;
    for my $format_tag (@format_fields) {
        push @final_sample_values, $sample_values{$format_tag};
    }
    my $final_sample_string = join ":", @final_sample_values;

    # Join the line up and return it, using the new alt string and new sample string
    my $final_vcf_line = join "\t", ($chrom, $pos, $id, $ref, $new_alt, $qual, $filter, $info, $format, $final_sample_string);
    return $final_vcf_line;
}

# Given a comma separated list of alts and the read_bases and base_quality strings from pileup
# Generate and return comma separated strings for allele depth and base quality for the vcf
# Now we also include the reference base AD and BQ values.
sub calculate_stats_for_alts {
    my $self = shift;
    my $ref = shift;
    my $alt_string = shift;
    my $read_bases = shift;
    my $base_qualities = shift;

    # Parse the pileup and quality strings so that they have the same length and can be mapped to one another
    if ($read_bases =~ m/[\$\^\+-]/) {
        $read_bases =~ s/\^.//g; #removing the start of the read segment mark
        $read_bases =~ s/\$//g; #removing end of the read segment mark
        while ($read_bases =~ m/[\+-]{1}(\d+)/g) {
            my $indel_len = $1;
            $read_bases =~ s/[\+-]{1}$indel_len.{$indel_len}//; # remove indel info from read base field
        }
    }
    if ( length($read_bases) != length($base_qualities) ) {
        die $self->error_message("After processing, read base string and base quality string do not have identical lengths: $read_bases $base_qualities");
    }

    # Count the number of times the variant occurs in the pileup string (AD) and its quality from the quality string (BQ)
    my @alts = split ",", $alt_string;
    my %ad;
    # Iinitialize values
    for my $alt ($ref, @alts) {
        $ad{$alt} = 0;
    }
    my %bq_total;
    my %bq;
    my (@bq_values, @ad_values);
    for my $alt ($ref, @alts) {
        my @bases = split("", $read_bases);
        my @qualities = split("", $base_qualities);
        for (my $index = 0; $index < scalar(@bases); $index++) {
            my $base = $bases[$index];
            if (lc $base eq lc $alt) { 
                #http://samtools.sourceforge.net/pileup.shtml base quality is the same as mapping quality
                $ad{$alt}++;
                $bq_total{$alt} += ord($qualities[$index]) - 33;
            }
        }

        # Get an average of the quality for BQ
        if ($ad{$alt} == 0) {
            $bq{$alt} = 0;
        } else {
            $bq{$alt} = int($bq_total{$alt} / $ad{$alt});
        }
        push @ad_values, $ad{$alt};
        push @bq_values, $bq{$alt};
    }

    my $ad_string = join(",", @ad_values );
    my $bq_string = join(",", @bq_values );

    return ($ad_string, $bq_string);
}

#return -1 if $chr_a,$pos_a represents a lower position than $chr_b,$pos_b, 0 if they are the same, and 1 if b is lower
sub compare_lines {
    my $self = shift;
    my ($line_a, $line_b) = @_;
    my ($chr_a, $pos_a) = split "\t", $line_a;
    my ($chr_b, $pos_b) = split "\t", $line_b;

    if(($chr_a eq $chr_b) && ($pos_a == $pos_b)){
        return 0;
    }
    if($chr_a eq $chr_b){
        return ($pos_a < $pos_b) ? -1 : 1;
    }
    return ($self->chr_cmp($chr_a,$chr_b)) ? 1 : -1;
}

# return 0 if $chr_a is lower than $chr_b, 1 otherwise
sub chr_cmp {
    my $self = shift;
    my ($chr_a, $chr_b) = @_;
    my @chroms = ($chr_a,$chr_b);
    my @answer = nsort @chroms;
    return ($answer[0] eq $chr_a) ? 0 : 1;
}

# This method is called when you are joining multiple samples together with potentially different ALT strings.
# It will figure out what the new GT string should be based upon the old GT string, reference, new and old ALT values
# old and new ALT are expected to be comma separated strings (as in VCF files).
# old_gt is expected to be the / separated list (0/1) found in VCF files.
sub regenerate_gt {
    my ($self, $reference, $old_alt, $old_gt, $new_alt) = @_;

    # Keep it null if it was null before
    if ($old_gt eq ".") {
        return $old_gt;
    }

    my @old_alt = split(",", $old_alt);
    my @old_gt = split("/", $old_gt);

    # This assumption may change at some point but for now we assume every VCF will have two values in the GT
    unless (scalar (@old_gt) >= 2) {
        die $self->error_message("This GT does not appear to have at least two values: $old_gt");
    }

    # Translate the old GT string into a list of alleles that sample contained
    my @alleles;
    for my $genotype_number (@old_gt) {
        my $allele;
        if ($genotype_number == 0) {
            $allele = $reference;
        } else {
            $allele = $old_alt[$genotype_number - 1]; # Genotype number will be 1 based in regards to the alt string
            unless (defined $allele) {
                die $self->error_message("Could not match genotype number $genotype_number to any allele from the ALT field $old_alt");
            }
        }

        push(@alleles, $allele);
    }

    # Now that we have the list of alternate alleles from the original line/ALT ... calculate what the new GT string should be based upon the new ALT
    my @new_alt_alleles = split(",", $new_alt);

    return $self->generate_gt($reference, \@new_alt_alleles, \@alleles);
}

# Generates the "GT" field. A 0 indicates matching reference. Any other number indicates matching that variant in the available "alt" alleles.
# I.E. REF: A ALT: C,T ... a A/C call in the GT field would be: 0/1. A C,T call in the GT field would be: 1/2
# alt alleles is an arrayref of  the alleles from the "ALT" column, all calls for this position that don't match the reference.
# genotype alleles is an arrayref of the alleles called at this position for this sample, including those that match the reference
sub generate_gt {
    my ($self, $reference, $alt_alleles, $genotype_alleles) = @_;

    my @gt_string;
    for my $genotype_allele (@$genotype_alleles) {
        my $allele_number;
        if ($genotype_allele eq $reference) {
            $allele_number = 0;
        } else {
            # Find the index of the alt allele that matches this genotype allele, add 1 to offset 0 based index
            for (my $i = 0; $i < scalar @$alt_alleles; $i++) {
                if ($genotype_allele eq @$alt_alleles[$i]) {
                    $allele_number = $i + 1; # Genotype index starts at 1
                }
            }
        }
        unless (defined $allele_number) {
            die $self->error_message("Could not match genotype allele $genotype_allele to any allele from the ALT field");
        }

        push(@gt_string, $allele_number);
    }

    # the GT field is sorted out of convention... you'll see 0/1 but not 1/0
    return join("/", sort(@gt_string));
}


1;
