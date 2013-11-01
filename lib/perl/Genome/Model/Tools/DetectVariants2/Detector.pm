package Genome::Model::Tools::DetectVariants2::Detector;

use strict;
use warnings;

use File::Copy;
use File::Basename;
use Genome;

class Genome::Model::Tools::DetectVariants2::Detector {
    is => ['Genome::Model::Tools::DetectVariants2::Base'],
    is_abstract => 1,
    has_optional => [
        params => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => 'The full parameter list coming in from the dispatcher. It is one string before being parsed.',
        },
        version => {
            is => 'String',
            is_input => 1,
            doc => 'The version of the detector to run',
        },
        region_of_interest => {
            is => 'Genome::FeatureList',
            doc => '',
            id_by => 'region_of_interest_id',
        },
        region_of_interest_id => {
            is => 'Text',
            doc => 'FeatureList for the region of interest (if present, only variants in the set will be reported)',
            is_input => 1,
        },
        _snv_base_name => {
            is => 'Text',
            default_value => 'snvs.hq',
            is_input => 1,
        },
        snv_output => {
            calculate_from => ['_snv_base_name', 'output_directory'],
            calculate => q{ join("/", $output_directory, $_snv_base_name); },
            doc => "Where the SNV output should be once all work has been done",
            is_output => 1,
        },
        snv_bed_output => {
            calculate_from => ['_snv_base_name', 'output_directory'],
            calculate => q{ join("/", $output_directory, $_snv_base_name) . ".bed"; },
            doc => "Where the SNV output which has been converted to .bed format should be once all work has been done",
            is_output => 1,
        },
        _snv_staging_output => {
            calculate_from => ['_temp_staging_directory', '_snv_base_name'],
            calculate => q{ join("/", $_temp_staging_directory, $_snv_base_name); },
            doc => 'Where the SNV output should be generated (It will be copied to the snv_output in _promote_staged_data().)',
        },
        _indel_base_name => {
            is => 'Text',
            default_value => 'indels.hq',
            is_input => 1,
        },
        indel_output => {
            calculate_from => ['_indel_base_name', 'output_directory'],
            calculate => q{ join("/", $output_directory, $_indel_base_name); },
            is_output => 1,
        },
        indel_bed_output => {
            calculate_from => ['_indel_base_name', 'output_directory'],
            calculate => q{ join("/", $output_directory, $_indel_base_name) . ".bed"; },
            is_output => 1,
        },
        _indel_staging_output => {
            calculate_from => ['_temp_staging_directory', '_indel_base_name'],
            calculate => q{ join("/", $_temp_staging_directory, $_indel_base_name); },
        },
        _sv_base_name => {
            is => 'Text',
            default_value => 'svs.hq',
            is_input => 1,
        },
        sv_output => {
            calculate_from => ['_sv_base_name', 'output_directory'],
            calculate => q{ join("/", $output_directory, $_sv_base_name); },
            is_output => 1,
        },
        _sv_staging_output => {
            calculate_from => ['_temp_staging_directory', '_sv_base_name'],
            calculate => q{ join("/", $_temp_staging_directory, $_sv_base_name); },
        },
        _filtered_indel_base_name => {
            is => 'Text',
            default_value => 'indels_all_sequences.filtered',
            is_input => 1,
        },
        filtered_indel_output => {
            calculate_from => ['_filtered_indel_base_name', 'output_directory'],
            calculate => q{ join("/", $output_directory, $_filtered_indel_base_name); },
            is_output => 1,
        },
        filtered_indel_bed_output => {
            calculate_from => ['_filtered_indel_base_name', 'output_directory'],
            calculate => q{ join("/", $output_directory, $_filtered_indel_base_name) . ".bed"; },
            is_output => 1,
        },
        _filtered_indel_staging_output => {
            calculate_from => ['_temp_staging_directory', '_filtered_indel_base_name'],
            calculate => q{ join("/", $_temp_staging_directory, $_filtered_indel_base_name); },
        },
    ],
    has_optional_transient => [
        _result => {
            is => 'UR::Object',
            doc => 'SoftwareResult for the run of this detector',
            id_by => "_result_id",
            id_class_by => '_result_class',
            is_output => 1,
        },
        _vcf_result => {
            is => 'UR::Object',
            doc => 'SoftwareResult for the vcf output of this detector',
            id_by => "_vcf_result_id",
            id_class_by => '_vcf_result_class',
            is_output => 1,
        },
        _result_class => {
            is => 'Text',
            is_output => 1,
        },
        _result_id => {
            is => 'Number',
            is_output => 1,
        },
        _vcf_result_class => {
            is => 'Text',
            is_output => 1,
        },
        _vcf_result_id => {
            is => 'Number',
            is_output => 1,
        },
        _previous_output_directory => {
            is => 'Text',
            is_input => 1,
            doc => 'If this detector is the same as another for a different variant type and is the second run, the prior run will be referenced here',
        },
    ],
    has_constant => [
        #These can't be turned off--just pass no detector name to skip
        detect_snvs => { value => 1 },
        detect_indels => { value => 1 },
        detect_svs => { value => 1 },
    ],
    has_param => [
        lsf_queue => {
            default => 'apipe',
        },
    ],
    doc => 'This is the base class for all detector classes',
};

sub help_brief {
    my $self = shift;
    my $class = ref($self) || $self;
    my ($name) = ($class =~ /Genome::Model::Tools::DetectVariants2::(.*)/);
    my @words = map { lc($_) } split(/(?=[A-Z])/,$name);
    return "directly run the @words variant detector";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
This is just an abstract base class for variant detector modules.
EOS
}

sub help_detail {
    return <<EOS
This is just an abstract base class for variant detector modules.
EOS
}

sub _supports_cross_sample_detection {
    my ($class, $version, $vtype, $params) = @_;
    return 0; # not by default
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    for my $input ('aligned_reads_input', 'control_aligned_reads_input') {
        if($self->$input) {
            my $canonical_path = Cwd::abs_path($self->$input);
            unless($canonical_path) {
                die $self->error_message('Failed to resolve real path to ' . $input);
            }

            $self->$input($canonical_path);
        }
    }

    return $self;
}

sub shortcut {
    my $self = shift;

    $self->_resolve_output_directory;

    $self->status_message("Attempting to shortcut detector result");
    unless($self->shortcut_detector){
        $self->status_message("Could not shortcut detector result.");
        return;
    }

    if($self->_try_vcf){
        $self->status_message("Attempting to shortcut vcf result");
        unless($self->shortcut_vcf){
            $self->status_message("Could not shortcut vcf result.");
            return;
        }
    }

    return 1;
}

sub _try_vcf {
    my $self = shift;
    my @types;
    for ("snvs","indels"){
        my $func = "detect_".$_;
        if($self->$func){
            push @types,$_;
        }
    }

    my $try_vcf=undef;

    for (@types){
        if(Genome::Model::Tools::DetectVariants2::Result::Vcf->conversion_class_name($self->class,$_)){
            return 1;
        }
    }
    return 0;
}

sub shortcut_detector {
    my $self = shift;
    my ($params) = $self->params_for_detector_result;
    $self->status_message("Params for shortcut_detector: " . Data::Dumper::Dumper $params);
    my $result = Genome::Model::Tools::DetectVariants2::Result->get_with_lock(%$params);
    unless($result) {
        $self->status_message('No existing result found.');
        return;
    }

    $self->_result($result);
    $self->status_message('Using existing result ' . $result->__display_name__);
    $self->_link_output_directory_to_result;

    return 1;
}

sub shortcut_vcf {
    my $self = shift;
    my ($params) = $self->params_for_vcf_result;
    $self->status_message("Params for shortcut_vcf: " . Data::Dumper::Dumper $params);
    my $result = Genome::Model::Tools::DetectVariants2::Result::Vcf::Detector->get_with_lock(%$params);
    unless($result) {
        $self->status_message('No existing result found.');
        return;
    }

    $self->_vcf_result($result);
    $self->status_message('Using existing result ' . $result->__display_name__);
    $self->_link_vcf_output_directory_to_result;

    return 1;
}

sub _resolve_output_directory {
    my $self = shift;
    #Subclasses override this
    return 1;
}


sub execute {
    my $self = shift;

    $self->_resolve_output_directory;
    $self->_summon_detector_result;

    if($self->_try_vcf){
        $self->_summon_vcf_result;
    }

    return 1;
}

sub _summon_vcf_result {
    my $self = shift;

    my ($params) = $self->params_for_vcf_result;
    my $result = Genome::Model::Tools::DetectVariants2::Result::Vcf::Detector->get_or_create(%$params); #, _instance => $self);

    unless($result) {
        die $self->error_message('Failed to create generate vcf result!');
    }

    $self->_vcf_result($result);
    $self->status_message('Generated vcf result.');
    $self->_link_vcf_output_directory_to_result;

    return 1;
}

sub _summon_detector_result {
    my $self = shift;

    my ($params) = $self->params_for_detector_result;
    my $result = Genome::Model::Tools::DetectVariants2::Result->get_or_create(%$params, _instance => $self);
    unless($result) {
        die $self->error_message('Failed to create generate detector result!');
    }

    $self->_result($result);
    $self->status_message('Generated detector result.');
    unless(-e $self->output_directory){
        $self->_link_output_directory_to_result;
    }

    return 1;
}

sub _generate_result {
    my $self = shift;

    unless($self->_verify_inputs) {
        die $self->error_message('Failed to verify inputs.');
    }

    unless($self->_create_directories) {
        die $self->error_message('Failed to create directories.');
    }

    #disconnect database before long-running commands
    Genome::DataSource::GMSchema->disconnect_default_dbh if Genome::DataSource::GMSchema->has_default_dbh;

    unless($self->_detect_variants) {
        die $self->error_message('Failed in main execution logic.');
    }

    unless($self->_sort_detector_output){
        die $self->error_message('Failed in _sort_detector_output');
    }

    unless($self->_generate_standard_files) {
        die $self->error_message('Failed to generate standard files from detector-specific files');
    }

    unless($self->_promote_staged_data) {
        die $self->error_message('Failed to promote staged data.');
    }

    return 1;
}

sub _sort_detector_output {
    my ($self, $skip) = @_;

    my @detector_files = glob($self->_temp_staging_directory."/*.hq");

    for my $detector_file (@detector_files){
        my $detector_unsorted_output = $self->_temp_scratch_directory . "/" . basename($detector_file) . ".unsorted";

        unless(rename($detector_file, $detector_unsorted_output)) {
            my $m = sprintf(q(Failed to move '%s' to '%s' for sorting: %s), $detector_file, $detector_unsorted_output, $!);
            $self->error_message($m);
            return;
        }

        my $sort_cmd = Genome::Model::Tools::Bed::ChromSort->create(
            input => $detector_unsorted_output,
            output => $detector_file,
            skip_lines => $skip,
        );

        unless ($sort_cmd->execute()) {
            $self->error_message("Failed to sort detector file " . $detector_unsorted_output);
            return;
        }
    }

    return 1;
}


sub params_for_detector_result {
    my $self = shift;

    my @alignment_results = $self->alignment_results;
    my @control_alignment_results = $self->control_alignment_results;

    my %params = (
        detector_name => $self->class,
        detector_params => $self->params,
        detector_version => $self->version,

        # old
        aligned_reads => $self->aligned_reads_input,
        control_aligned_reads => $self->control_aligned_reads_input,

        #new
        alignment_results => (@alignment_results? [map { $_->id } @alignment_results] : undef),
        control_alignment_results => (@control_alignment_results? [map { $_->id } @control_alignment_results] : undef),
        pedigree_file_path => $self->pedigree_file_path,
        roi_list => $self->roi_list,
        roi_wingspan => $self->roi_wingspan,

        reference_build_id => $self->reference_build_id,
        region_of_interest_id => $self->region_of_interest_id,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
        chromosome_list => undef,
    );

    return \%params;
}

sub params_for_vcf_result {
    my $self = shift;
    my $vcf_version = Genome::Model::Tools::Vcf->get_vcf_version;
    my %params = (
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
        input_id => $self->_result->id,
        vcf_version => $vcf_version,
        aligned_reads_sample => $self->aligned_reads_sample,
    );
    $params{control_aligned_reads_sample} = $self->control_aligned_reads_sample if defined $self->control_aligned_reads_sample;

    return \%params;
}

sub _link_output_directory_to_result {
    my $self = shift;

    my $result = $self->_result;
    return unless $result;

    unless(-e $self->output_directory) {
        Genome::Sys->create_symlink($result->output_dir, $self->output_directory);
    }

    return 1;
}

sub _link_vcf_output_directory_to_result {
    my $self = shift;
    $self->status_message("Linking in vcfs from vcf_result");

    my $result = $self->_vcf_result;
    return unless $result;
    my @vcfs = glob($result->output_dir."/*.vcf.gz");
    my $output_directory = $self->output_directory;
    for my $vcf (@vcfs){
        my $target = $output_directory . "/" . basename($vcf);
        $self->status_message("Attempting to link : " .$vcf."  to  ". $target);
        if(-l $target) {
            if (readlink($target) eq $vcf) {
                $self->status_message("Already found a vcf linked in here, and it already has the correct target. Continuing.");
                next;
            } else {
                $self->status_message("Already found a vcf linked in here, unlinking that for you.");
                unless(unlink($target)){
                    die $self->error_message("Failed to unlink a link to a vcf at: ".$target);
                }
            }
        } elsif(-e $target) {
            die $self->error_message("Found something in place of the vcf symlink.");
        }
        # Symlink both the vcf and the tabix
        Genome::Sys->create_symlink($vcf, $target);
        Genome::Sys->create_symlink("$vcf.tbi", "$target.tbi");
    }
    return 1;
}

# Given a line of output from this detector, parse and return the chromosome, position, reference, and variant
# The position must be converted to the same position that a bed would consider the STOP position
# This is used for intersecting the detector specific file with the bed version
# Override this method in each detector if the format varies from this
#TODO clean all of this up. It is usually/should be based on logic from Genome::Model::Tools::Bed::Convert logic in process_source...
# this should be smarter about using that work ... perhaps process_source should call a method that just parses one line, and this method can be replaced by a call to that instead
sub parse_line_for_bed_intersection {
    my $class = shift;
    my $line = shift;

    unless ($line) {
        die $class->error_message("No line provided to parse_line_for_bed_intersection");
    }

    my ($chromosome, $position, $reference, $variant) = split "\t",  $line;

    unless (defined $chromosome && defined $position && defined $reference && defined $variant) {
        die $class->error_message("Could not get chromosome, position, reference, or variant for line: $line");
    }

    return [$chromosome, $position, $reference, $variant];
}

1;
