package Genome::InstrumentData::AlignmentResult::Merged::RefCov;

use strict;
use warnings;

use Genome;
use Sys::Hostname;
use File::Path;


class Genome::InstrumentData::AlignmentResult::Merged::RefCov{
    is => ['Genome::SoftwareResult::Stageable'],
    has_input => [
        alignment_result_id => {
            is => 'Number',
            doc => 'ID of the result for the alignment data upon which to run refcov',
        },
        region_of_interest_set_id => {
            is => 'Text',
            doc => 'ID of the feature list containing the regions over which to gather refcov',
        },
    ],
    has_param => [
        minimum_depths => {
            is => 'Text',
            doc => 'comma-separated list of minimum depths at which to evaluate coverage',
        },
        use_short_roi_names => {
            is => 'Boolean',
            doc => 'Whether or not to shorten the names in the BED file for processing',
        },
        roi_track_name => {
            is => 'Text',
            doc => 'For multi-tracked ROI use this named track',
        },
	maximum_depth => {
            is => 'Integer',
            doc => 'The maximum depth to evaluate by samtools pileup.',
            is_optional => 1,
        },
	alignment_count => {
            is => 'Boolean',
            doc => 'Calculate the number of alignments that overlap region.',
            is_optional => 1,
        },
        print_min_max => {
            is => 'Boolean',
            doc => 'Print the minimum and maximum depth of coverage.',
            is_optional => 1,
        },
        output_directory => {
            doc => 'When run in parallel, this directory will contain all output and intermediate STATS files. Sub-directories will be made for wingspan and min_depth_filter params. Do not define if stats_file is defined.',
            is_optional => 1,
        },
        print_headers => {
            is => 'Boolean',
            doc => 'Print a header describing the output including column headers.',
            is_optional => 1,
        },
        embed_bed => {
            is		  => 'Boolean',
            doc		  => 'Include associated BED information (REF:START-STOP) per target ID when reporting. Returns 0-based BED coords.',
            is_optional	  => 1,
        },
        normalize_with_formula => {
            is		  => 'String',
            doc		  => 'Runs normalization on AVE DEPTH and MIN & MAX DEPTHs (opt) given a Perl compatible formula, replacing $X with depth value per ROI. EXAMPLE: --normalize-with-formula=\'log( ( $X / 250_000_000 ) * 1_000_000 ) / log( 2 )\'',
            is_optional	  => 1,
        },

    ],
    has_metric => [
        _log_directory => {
            is => 'Text',
            doc => 'Path where workflow logs were written',
        },
        #many other metrics exist--see sub _generate_metrics
    ],
    has => [
        alignment_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged',
            id_by => 'alignment_result_id',
            doc => 'the alignment data upon which to run refcov',
        },
        region_of_interest_set => {
            is => 'Genome::FeatureList',
            id_by => 'region_of_interest_set_id',
            doc => 'regions over which to gather refcov',
        },
    ],
    has_transient_optional => [
        log_directory => {
            is => 'Text',
            doc => 'Path to write logs from running the workflow',
        },
    ],
};

sub resolve_allocation_subdirectory {
    my $self = shift;

    my $hostname = hostname;
    my $user = $ENV{'USER'};
    my $base_dir = sprintf("refcov-%s-%s-%s-%s", $hostname, $user, $$, $self->id);

    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'build_merged_alignments','refcov_stats',$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub _staging_disk_usage {
    #need the allocation created in advance for this process
    return 5_000_000; #TODO better estimate
}

sub _working_dir_prefix {
    return 'refcov-stats';
}

sub _prepare_staging_directory {
    my $self = shift;

    return $self->temp_staging_directory if ($self->temp_staging_directory);

    unless($self->output_dir) {
        $self->_prepare_output_directory;
    }

    #Stage to network disk because of inner workflow
    my $staging_tempdir = File::Temp->newdir(
        $self->_working_dir_prefix . '-staging-XXXXX',
        DIR     => $self->output_dir,
        CLEANUP => 1,
    );

    $self->temp_staging_directory($staging_tempdir);
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    $self->_prepare_staging_directory;

    my $bed_file = $self->_dump_bed_file;
    my $bam_file = $self->alignment_result->merged_alignment_bam_path;

    die $self->error_message("Bed File ($bed_file) is missing") unless -s $bed_file;
    die $self->error_message("Bam File ($bam_file) is missing") unless -s $bam_file;
    
	
    my $log_dir = $self->log_directory;
    unless($log_dir) {
        $log_dir = '' . $self->temp_staging_directory;
    }
    $self->_log_directory($log_dir);

    my %refcov_params = (
        output_directory => '' . $self->temp_staging_directory,
        roi_file_path => $bed_file,
        alignment_file_path => $bam_file,
        min_depth_filter => $self->minimum_depths,
	embed_bed => $self->embed_bed,
	alignment_count =>$self->alignment_count,
	print_min_max => $self->print_min_max,
	print_headers => $self->print_headers,
	roi_file_format => 'bed',
	alignment_file_format => 'bam',
   );

    unless($] > 5.010) {
        #need to shell out to a newer perl #TODO remove this once 5.10 transition complete
        my $cmd = 'genome-perl5.10 -S gmt ref-cov standard ';
        while (my ($key, $value) = (each %refcov_params)) {
            $key =~ s/_/-/g;
            $cmd .= " --$key=$value";
        }

        Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files => [$bed_file, $bam_file],
        );
    } else {
        my $cmd = Genome::Model::Tools::RefCov::Standard->create(%refcov_params);
        unless($cmd->execute) {
            die('Failed to run RefCov');
        }
    }

    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    $self->_generate_metrics;

    #quick sanity check--all wingspans should run on the same number of total reads
    return $self;
}

sub _dump_bed_file {
    my $self = shift;

    my $roi_set = $self->region_of_interest_set;
    return unless $roi_set;

    my $alt_reference;
    my $reference = $self->alignment_result->reference_build;
    unless($reference->is_compatible_with($roi_set->reference)) {
        $alt_reference = $reference;
    }
    my $use_short_names = $self->use_short_roi_names;

    my $bed_file_path = $self->temp_staging_directory .'/'. $roi_set->id .'.bed';
    unless (-e $bed_file_path) {
        my %dump_params = (
            feature_list => $roi_set,
            output_path => $bed_file_path,
            alternate_reference => $alt_reference,
            short_name => $use_short_names,
        );
        if ($self->roi_track_name) {
            $dump_params{track_name} = $self->roi_track_name;
        }
        my $dump_command = Genome::FeatureList::Command::DumpMergedList->create(%dump_params);
        unless ($dump_command->execute) {
            die('Failed to print bed file to path '. $bed_file_path);
        }
    }

    return $bed_file_path;
}


sub stats_file {
    my $self = shift;
    my $refcov_stats_directory = $self->output_dir;
    my @stats_files = glob($refcov_stats_directory.'/*STATS.tsv');
    unless (@stats_files) {
        return;
    }
    unless (scalar(@stats_files) == 1) {
        die("Found multiple stats files:\n". join("\n",@stats_files));
}

}


sub _generate_metrics
{
return 1;
}

1;
