package Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcfDenovo;

use warnings;
use strict;

use Genome;
use Workflow;
use Workflow::Simple;
use Carp;
use Data::Dumper;
use Genome::Utility::Vcf ('parse_vcf_line', 'deparse_vcf_line', 'get_samples_from_header');

class Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcfDenovo {
    is => 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcfBase',
    doc => "This module uses detailed readcount information from bam-readcounts to filter likely false positives",
};

sub output_file_path {
    my $self = shift;
    return $self->_temp_staging_directory . "/snvs.vcf.gz";
}

sub input_file_path {
    my $self = shift;
    return $self->input_directory . "/snvs.vcf.gz";
}

sub should_skip_filter {
    my $self = shift;

    return -z $self->region_path;
}

sub header_already_added {
    my ($self, $header) = @_;
    return grep { $_ =~/FORMAT=<ID=DNFT,/ } @$header;
}

sub filter_status_header {
    return qq{##FORMAT=<ID=DNFT,Number=1,Type=String,Description="Denovo Filter Status">\n};
}

sub _convert_to_standard_formats {
}

sub should_print_region_line {
    my ($self, $line) = @_;

    return $line =~ m/DA=/;
}

sub open_input_file {
    my ($self, $input_file) = @_;

    return Genome::Sys->open_gzip_file_for_reading($input_file);
}

sub set_info_field {
    my ($self, $parsed_line, $info_tag, $info_value) = @_;
    if(!exists($parsed_line->{info}{$info_tag})) {
        push @{$parsed_line->{'_info_tags'}}, $info_tag;
    }
    $parsed_line->{info}{$info_tag} = $info_value;   #should really do type checking...
}


#FIXME probably move this to a base class
# Format fields this filter requires, override in each filter
sub required_format_fields {
#    return qw(DNGT);
}

# Given a parsed vcf line structure, filter each sample on the line
sub filter_one_line {
    my $self = shift;
    my $parsed_vcf_line = shift;
    my $readcount_searcher_by_sample = shift;
    my $stats = shift;
    my $samples = shift;

    # FIXME this will be only correct for the number of lines we have
    $stats->{'num_variants'}++;

    # FIXME run this for each sample in the line that is not "." and has a non ref GT
    my $denovo_allele = $parsed_vcf_line->{info}{"DA"};
    if(!defined($denovo_allele)) {
        return;
    }
    my $alt = $parsed_vcf_line->{alt};
    my @alts = split ",", $alt;
    my ($dn_allele) = $self->convert_numeric_gt_to_alleles(\@alts, [$denovo_allele], $parsed_vcf_line->{reference});


    my $denovo_found=0;
    for my $sample_name (@$samples) {
        my $gt = $parsed_vcf_line->{sample}{$sample_name}{DNGT};
        my @gt = split("/", $gt);
        my @alleles  = $self->convert_numeric_gt_to_alleles(\@alts, \@gt, $parsed_vcf_line->{reference});
        for my $allele(@alleles) {
            if($allele eq $dn_allele) {
                $self->filter_one_sample($parsed_vcf_line, $readcount_searcher_by_sample, $stats, $sample_name, $dn_allele);
                if ($denovo_found) {
                    $self->warning_message("I've already filtered a denovo allele, something might be wrong or we have more than one family with a denovo allele.");
                }
                $denovo_found=1;
            }
        }
    }
    return 1;
}

sub filter_sample_format_tag {
    return 'DNFT';
}

sub update_variant_for_sample {
    my ($self, $parsed_vcf_line, $sample_name, $var) = @_;

    return $var;
}

#override the default scratch directories in order to allow for network available temp dirs
sub _create_temp_directories {
    my $self = shift;
    $self->_temp_staging_directory($self->output_directory);
    $self->_temp_scratch_directory($self->output_directory);
    return 1;
}

sub _promote_staged_data {
    my $self = shift;
    return 1;
}

1;
