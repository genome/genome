package Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcfDenovo;

use warnings;
use strict;

use Genome;
use Workflow;
use Workflow::Simple;
use Carp;
use Data::Dumper;
use Genome::Utility::Vcf ('parse_vcf_line', 'deparse_vcf_line');

class Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcfDenovo {
    is => 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositiveVcfBase',
    doc => "This module uses detailed readcount information from bam-readcounts to filter likely false positives",
};

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
    my $input_file = $self->input_directory . "/snvs.vcf.gz";
    ## Build temp file for positions where readcounts are needed ##
    #my $region_path = $self->_temp_scratch_directory."/regions";
    my $region_path = $self->output_directory."/regions";

    $self->print_region_list($input_file, $region_path);

    my $output_file = $self->_temp_staging_directory . "/snvs.vcf.gz";
    unless(-s $region_path) { #no denovo alleles in this file
        $self->status_message("No denovo alleles found, copying file over and reporting success (this filter was a no-op");
        `cp $input_file $output_file`;
        return 1;  ### pass the file along and report successful
    }
    my $output_fh = Genome::Sys->open_gzip_file_for_writing($output_file);
    ## Run BAM readcounts in batch mode to get read counts for all positions in file ##
    my $readcount_searcher_by_sample = $self->generate_and_run_readcounts_in_parallel($region_path);

    # Advance to the first variant line and print out the header
    my ($input_fh, $header) = $self->parse_vcf_header($input_file);
    #check here to see if header has FT format tag
    unless(grep { $_ =~/FORMAT=<ID=DNFT,/ } @$header) {
        my $col_header = $header->[-1];
        $header->[-1] = qq{##FORMAT=<ID=DNFT,Number=1,Type=String,Description="Denovo Filter Status">\n};
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
    $output_fh->print(join("",@$header));
    my @sample_names = $self->get_samples_from_header($header);

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

    return 1;
}

sub should_print_region_line {
    my ($self, $line) = @_;

    return $line =~ m/DA=/;
}

sub open_input_file {
    my ($self, $input_file) = @_;

    return Genome::Sys->open_gzip_file_for_reading($input_file);
}

# Given the header of a vcf, return an array of samples in the final header line
#FIXME probably move this to a base class
sub get_samples_from_header {
    my $self = shift;
    my $header = shift;

    my $sample_line = @$header[-1];
    chomp $sample_line;
    my @fields = split "\t", $sample_line;
    splice(@fields, 0, 9);

    return @fields;
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

# In the line provided, set the sample provided to failed (FT)
sub fail_sample {
    my $self = shift;
    my $parsed_vcf_line = shift;
    my $sample_name = shift;
    my $filter_reason = shift;

    $self->status_message("Entering fail sample\n");
    # XXX @gsanders Should this have the optiosn from FPVcf?
    $self->set_format_field($parsed_vcf_line, $sample_name,
        $self->filter_sample_format_tag, $filter_reason);
    return 1;
}

sub generate_and_run_readcounts_in_parallel {
    my $self = shift;
    my $region_path = shift;

    my %inputs;
    my (@outputs, @inputs);

    #set up global readcount params
    $inputs{reference_fasta} = $self->reference_sequence_input;
    $inputs{region_list} = $region_path;
    $inputs{minimum_base_quality} = $self->bam_readcount_min_base_quality;
    $inputs{use_version} = $self->bam_readcount_version;

    my @sample_names;
    for my $alignment_result ($self->alignment_results) {
        my $sample_name = $self->find_sample_name_for_alignment_result($alignment_result);
        push @sample_names, $sample_name;
        $self->status_message("Running BAM Readcounts for sample $sample_name...");

        my $readcount_file = $self->_temp_staging_directory . "/$sample_name.readcounts";  #this is suboptimal, but I want to wait until someone tells me a better way...multiple options exist
        push @outputs, $readcount_file;
        my $bam_path = $alignment_result->merged_alignment_bam_path;
        if (-f $bam_path) {
            $inputs{"bam_${sample_name}"} = $bam_path;
        } else {
            die "merged_alignment_bam_path does not exist: $bam_path";
        }
        $inputs{"readcounts_${sample_name}"} = $readcount_file;
    }

    my $workflow = Workflow::Model->create(
        name=> "FalsePositiveVcf parallel readcount file creation",
        input_properties => [
        keys %inputs
        ],
        output_properties => [
        'output',
        ],
    );
    for my $sample (@sample_names) {
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
    for my $sample (@sample_names) {
        $readcount_searcher_by_sample->{$sample} = $self->make_buffered_rc_searcher(Genome::Sys->open_file_for_reading($inputs{"readcounts_$sample"}));
    }
    return $readcount_searcher_by_sample;
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
