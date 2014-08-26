package Genome::Model::Tools::Vcf::VcfFilter;
##########################################################################
# Given a VCF file and a list of things to filter, removes them or
# marks them appropriately in the VCF file
#
#
#       AUTHOR:         Chris Miller (cmiller@genome.wustl.edu)
#
#       CREATED:        05/04/2011 by CAM
#       MODIFIED:
#
#       NOTES:
#
###########################################################################
use strict;
use warnings;
use Genome;
use Genome::Utility::Vcf "open_vcf_file";

class Genome::Model::Tools::Vcf::VcfFilter {
    is => 'Command',
    has => [
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => "filtered VCF file",
            is_optional => 0,
        },
        vcf_file => {
            is => 'Text',
            is_input => 1,
            doc => "mutations in Vcf format",
            is_optional => 0,
        },
        filter_file => {
            is => 'Text',
            doc => "files containing SNVs to be filtered (assumes first col CHR, second col POS",
            is_optional => 1,
            default => "",
        },
        filter_keep => {
            is => 'Boolean',
            doc => "the filter file contains variants that *passed* filters (all *other* SNVs will be marked invalid). If false, the opposite is assumed - the file contains variants that did not pass filters",
            is_optional => 1,
            default => 0,
        },
        filter_name => {
            is => 'Text',
            doc => "name to add to the FILTER field for variants newly marked filtered",
            is_optional => 1,
            default => "",
        },
        filter_description => {
            is => 'Text',
            doc => "description of the FILTER for the header",
            is_optional => 1,
            default => "",
        },
        remove_filtered_lines => {
            is => 'Boolean',
            is_optional => 1,
            default => 0 ,
            doc => 'remove the filtered lines, as opposed to marking them as non-passing in the VCF',
        },
        bed_input => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'filter file is in bed (0-based format). Default false (expects 1-based coordinates)',
            default=>0,
        },
        variant_type => {
            is => 'Text',
            is_optional => 1,
            doc => 'apply filters only to variants of this type (usually "SNP" or "INDEL")',
        },
        pass_when_alts_disagree => {
            is => 'Boolean',
            default => 1,
            doc => 'If there are two alts on a line and one is filtered and one is not, set the filter status to PASS. Make this false to instead mark as filtered in this case',
        },
    ],
};


sub help_brief {                            # keep this to just a few words <---
    "apply filter labels to a VCF file"
}


sub help_synopsis {
    <<'HELP';
    apply filter labels to a VCF file
HELP
}

sub help_detail {                  # this is what the user will see with the longer version of help. <---
    <<'HELP';

    Takes a VCF and a file containing filtered output, then annotates the VCF with the filter names, and adds the filter info  to the header
HELP
}




################################################################################################
# Execute - the main program logic
################################################################################################

sub execute {
    my $self = shift;
    my $output_file  = $self->output_file;
    my $vcf_file     = $self->vcf_file;
    my $filter_file  = $self->filter_file;
    my $filter_keep  = $self->filter_keep;
    my $filter_name  = $self->filter_name;
    my $variant_type = $self->variant_type;
    my $filter_description    = $self->filter_description;
    my $remove_filtered_lines = $self->remove_filtered_lines;

    my $bgzip_out = ($output_file =~ m/gz$/) ? 1 : 0;

    # first, read the filter file and store the locations
    my %filter;
    my $filter_fh = Genome::Sys->open_file_for_reading($filter_file);
    while( my $line = $filter_fh->getline ) {
        #skip header lines
        chomp $line;
        next if $line =~ /^#/;

        my @fields = split("\t", $line);
        if ($self->bed_input) {
            unless ($fields[3] =~ /\*|0|\-/) { #in bed, this would be an insertion or deletion and we should not naively ++ the start
                $fields[1]+=1;
            }
        }
        my ($ref, $var) = split "/", $fields[3];

        # If this is an insertion, we need the length of it to uniquely identify the variant between bed and vcf
        my $key;
        if ($ref eq "*" or $ref eq "-" or $ref eq "0") {
            $key = join(":", ($fields[0], $fields[1], $fields[2], length($var)));
        } 
        else {
            $key = join(":", ($fields[0], $fields[1], $fields[2]));
        }
        #add this filter to the hash
        $filter{$key} = 0;
    }
    $filter_fh->close;

    #open the output file
    my $outfile = ($bgzip_out) ? Genome::Sys->open_gzip_file_for_writing($output_file) : Genome::Sys->open_file_for_writing($output_file);

    #read the vcf
    my $inFh = open_vcf_file($vcf_file);

    my $found_pass_line    = 0;
    my $found_format_lines = 0;
    my $done_with_header   = 0;
    my $found_ft_header    = 0;

    #if this is a header line
    while(!$done_with_header) {
        my $line = $inFh->getline;
        chomp $line;
        if ($line =~ /^##/) {
            if ($line =~/^##FILTER=<ID=PASS/){
                $found_pass_line = 1;
            }

            # if this is the first FORMAT line, drop our
            # filter headers into the VCF here
            if ($line =~ /^##FORMAT/ && $found_format_lines == 0) {
                unless ($found_pass_line){
                    print $outfile "##FILTER=<ID=PASS,Description=\"Passed all filters\">" . "\n";
                }
                print $outfile "##FILTER=<ID=" . $filter_name . ",Description=\"" . $filter_description . "\">" . "\n";
                $found_format_lines = 1;
            }
            print $outfile $line . "\n";

            if ($line =~ /^##FORMAT=<ID=FT,/) {
                $found_ft_header = 1;
            }
        } 
        elsif ($line =~ /^#CHROM/) {
            $done_with_header = 1;
            unless ($found_ft_header) {
                print $outfile '##FORMAT=<ID=FT,Number=1,Type=String,Description="Sample genotype filter">' . "\n";
            }
            print $outfile $line . "\n";
        } 
        else {
            die $self->error_message("Failed to find the final header line");
        }
    }

    while( my $line = $inFh->getline ) {
        my $remove_line = 0;
        chomp $line;
        my @fields = split("\t",$line);

        # if this is not of the correct variant type, skip it
        if (defined($variant_type) && !($fields[7] =~ /VT=$variant_type/)) {
            print $outfile $line . "\n";
        } 
        else {
            my $final_filter_value = "";

            # only check the filters if this snv hasn't already been filtered
            # (is passing). If it has been filtered, then accept the prior filter
            # and move on
            my $filter_value;
            if (($fields[6] eq "") || ($fields[6] eq "PASS") || ($fields[6] eq ".")) {

                my $ref = $fields[3];
                my $alt_field = $fields[4];
                my @alts = split ",", $alt_field;

                my @filter_values;
                for my $alt (@alts) {
                    # Calculate the stop position that should be in the bed file for this variant
                    # If this is an insertion, we need the length of it to uniquely identify the variant between bed and vcf
                    # convert_indel_string is used to reset the indel start position. Some vcf like the ones from samtools mpileup
                    # will have some common bases shown in both ref and alt columns.  
                    # 1	83433754	.	ATTT	ATTTT
                    # 1	107286849	.	GATAT	GAT
                    my $key;
                    if (length($alt) > length($ref)) { # Insertion
                        my (undef, $start) = Genome::Model::Tools::Bed::Convert::Indel::SamtoolsToBed->convert_indel_string($fields[1], $ref, $alt); 
                        my $insertion_length = length($alt) - length($ref);
                        my $stop = $start;
                        $key = join(":", ($fields[0], $start, $stop, $insertion_length));
                    } 
                    elsif (length($alt) < length($ref)) {# Deletion
                        my (undef, $start) = Genome::Model::Tools::Bed::Convert::Indel::SamtoolsToBed->convert_indel_string($fields[1], $ref, $alt);
                        my $stop = $start + (length($ref) - length($alt));
                        $key = join(":", ($fields[0], $start, $stop));
                    } 
                    else {# SNV
                        my $stop = $fields[1]; 
                        $key = join(":", ($fields[0], $fields[1], $stop));
                    }

                    if ($filter_keep){
                        if (exists($filter{$key})){
                            $filter_value = "PASS";
                        } 
                        else {
                            if($remove_filtered_lines){
                                $remove_line = 1;
                            } 
                            else {
                                $filter_value = $filter_name;
                            }
                        }
                    } 
                    else {
                        if (exists($filter{$key})){
                            if($remove_filtered_lines){
                                $remove_line = 1;
                            } 
                            else {
                                $filter_value = $filter_name;
                            }
                        } 
                        else {
                            $filter_value = "PASS";
                        }
                    }
                    push @filter_values, $filter_value;
                }

                # Find out if there is any disagreement between alts as far as filter status
                my $at_least_one_passes = 0;
                my $at_least_one_fails  = 0;
                for my $filter_value (@filter_values) {
                    if ($filter_value eq "PASS") {
                        $at_least_one_passes = 1;
                    } 
                    else {
                        $at_least_one_fails = 1;
                    }
                }

                # Decide if the variant as a whole should be PASS or fail based upon flags
                if ($self->pass_when_alts_disagree) {
                    if ($at_least_one_passes) {
                        $final_filter_value = "PASS";
                    } 
                    else {
                        $final_filter_value = $filter_name;
                    }
                } 
                else {
                    if ($at_least_one_fails) {
                        $final_filter_value = $filter_name;
                    } 
                    else {
                        $final_filter_value = "PASS";
                    }
                }
                $fields[6] = $final_filter_value;
            } 
            else {
                $final_filter_value = $fields[6];
            }

            # Add a FT field with the same information as the filter field
            my $format_field = $fields[8];
            my @format_keys  = split ":", $format_field;
            # Find the FT field if it exists
            my $ft_index;
            for (my $i = 0; $i <= $#format_keys; $i++) {
                if ($format_keys[$i] eq "FT") {
                    $ft_index = $i;
                }
            }

            # If FT was not previously present, add it to the format field
            unless ($ft_index) {
                push @format_keys, "FT";
                $fields[8] = join(":", @format_keys);
                $ft_index = $#format_keys;
            }

            # For each sample present in the file, either replace the old FT value if one was present or insert a new one
            for (my $sample_index = 9; $sample_index <= $#fields; $sample_index++) {
                my $sample_field = $fields[$sample_index];
                my @sample_fields = split ":", $sample_field;

                $sample_fields[$ft_index] = $final_filter_value;

                for (0 .. $#sample_fields) {
                    $sample_fields[$_] = '.' unless defined $sample_fields[$_];
                }

                $fields[$sample_index] = join(":", @sample_fields);
            }

            #output the line
            unless($remove_line){
                print $outfile join("\t", @fields) . "\n";
            }
        }
    }
    $outfile->close;
    $inFh->close;

    return 1;
}
