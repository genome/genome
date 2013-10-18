package Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcf;

use warnings;
use strict;

use Genome;
use Workflow;
use Workflow::Simple;
use Carp;
use Data::Dumper;
use Genome::File::Vcf::Reader;
use Genome::Utility::Vcf ('open_vcf_file', 'parse_vcf_line', 'deparse_vcf_line', 'get_vcf_header', 'get_samples_from_header');

class Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcf {
    is => 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcfBase',
    doc => "This module uses detailed readcount information from bam-readcounts to filter likely false positives",
};

sub output_file_path {
    my $self = shift;
    return $self->_temp_staging_directory . "/snvs.raw.gz";
}

sub input_file_path {
    my $self = shift;

    # FIXME Things need to be cleaned up here...
    # We try to use the incoming snvs.hq file if it exists because of the refalign/normal pipelines
    # But currently, phenotype correlation only produces vcfs and nothing else... so fall back to using the vcf if snvs.hq doesnt exist.
    my $vcf_file = $self->input_directory . "/snvs.vcf.gz";
    my $hq_file = $self->input_directory . "/snvs.hq";
    if (-s $hq_file ) {
        return $hq_file;
    } elsif (-s $vcf_file) {
        return $vcf_file;
    } else {
        die $self->error_message("Could not find $vcf_file or $hq_file");
    }
}

sub should_skip_filter {
    return;
}

sub header_already_added {
    my ($self, $header) = @_;
    return grep { $_ =~/FORMAT=<ID=FT,/ } @$header;
}

sub filter_status_header {
    return qq{##FORMAT=<ID=FT,Number=1,Type=String,Description="Per Sample Filter Status">\n};
}

# Create snvs.hq snvs.lq and bed files
sub _convert_to_standard_formats {
    my ($self, $raw_file) = @_;


    # When multiple samples are present only VCF really makes sense. So don't convert to bed.
    my @alignment_results = $self->alignment_results;
    if ( @alignment_results and ( scalar(@alignment_results >= 2) )  ) {
        my $vcf_file = $self->_temp_staging_directory . "/snvs.vcf.gz";
        Genome::Sys->copy_file($raw_file, $vcf_file);
        unlink $raw_file;
        return 1;
    }

    my $reader = Genome::File::Vcf::Reader->new($raw_file);
    my $header = $reader->header;

    # Put header and only PASS variants in HQ file
    my $pass_snv_output_file = $self->_temp_staging_directory . "/snvs.hq";
    my $pass_fh = Genome::Sys->open_file_for_writing($pass_snv_output_file);
    $pass_fh->print($header->to_string);

    # Put header and only non PASS variants in LQ file
    my $fail_snv_output_file = $self->_temp_staging_directory . "/snvs.lq";
    my $fail_fh = Genome::Sys->open_file_for_writing($fail_snv_output_file);
    $fail_fh->print($header->to_string);

    while(my $entry = $reader->next) {
        my $pass = 1;
        for (my $sample_index=0; $sample_index < scalar($header->sample_names); $sample_index++) {
            my $ft_value = $entry->sample_field($sample_index, 'FT');
            unless ($ft_value eq 'PASS' or $ft_value eq '.') {
                $pass = 0;
            }
        }

        if ($pass) {
            $pass_fh->print($entry->to_string);
        } else {
            $fail_fh->print($entry->to_string);
        }
    }

    $pass_fh->close;
    $fail_fh->close;

    # Convert both files to bed
    my $convert = Genome::Model::Tools::Bed::Convert::Snv::SamtoolsToBed->create(
        source => $pass_snv_output_file,
        output => $self->_temp_staging_directory . "/snvs.hq.bed",
    );
    unless($convert->execute){
        $self->error_message("Failed to convert filter output to bed.");
        die $self->error_message;
    }

    my $convert_lq = Genome::Model::Tools::Bed::Convert::Snv::SamtoolsToBed->create(
        source => $fail_snv_output_file,
        output => $self->_temp_staging_directory . "/snvs.lq.bed",
    );
    unless($convert_lq->execute){
        $self->error_message("Failed to convert failed-filter output to bed.");
        die $self->error_message;
    }

    return 1;
}

sub should_print_region_line {
    return 1;
}

sub open_input_file {
    my ($self, $input_file) = @_;

    return open_vcf_file($input_file);
}

#FIXME probably move this to a base class
# Format fields this filter requires, override in each filter
sub required_format_fields {
    return qw(GT);
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
    for my $sample_name (@$samples) {
        $self->filter_one_sample($parsed_vcf_line, $readcount_searcher_by_sample, $stats, $sample_name);
    }

    return 1;
}

sub filter_sample_format_tag {
    return 'FT';
}

sub get_readcounts_from_vcf_line {
    my ($self, $parsed_vcf_line, $readcount_searcher_by_sample, $sample_name) = @_;

    my $readcount_searcher = $readcount_searcher_by_sample->{$sample_name};

    # This is a special case for when we are not considering two parents. #FIXME This is super hacky and should be improved somehow.
    if ((not defined $readcount_searcher) and (lc($sample_name) =~ m/fake_parent/i)) {
        $self->set_format_field($parsed_vcf_line, $sample_name,
            $self->filter_sample_format_tag, ".");
        return;
    }

    return $self->SUPER::get_readcounts_from_vcf_line($parsed_vcf_line,
        $readcount_searcher_by_sample, $sample_name);
}

sub update_variant_for_sample {
    my ($self, $parsed_vcf_line, $sample_name, $var) = @_;

    return $self->get_variant_for_sample($parsed_vcf_line->{alt},
        $parsed_vcf_line->{sample}->{$sample_name}->{GT},
        $parsed_vcf_line->{reference});
}

#override the default scratch directories in order to allow for network available temp dirs
sub _create_temp_directories {
    my $self = shift;
    local %ENV = %ENV;
    $ENV{TMPDIR} = $self->output_directory;
    return $self->SUPER::_create_temp_directories(@_);
}

sub _promote_staged_data {
    my $self = shift;
    my $output_dir = $self->SUPER::_promote_staged_data(@_);

    #remove the staging directory before the reallocation occurs so that it doesn't get double counted.
    my $staging_dir = $self->_temp_staging_directory;
    Genome::Sys->remove_directory_tree($staging_dir);

    return $output_dir;
}

sub _check_native_file_counts {
    my $self = shift;
    return $self->_check_native_file_counts_vcf;
}

1;
