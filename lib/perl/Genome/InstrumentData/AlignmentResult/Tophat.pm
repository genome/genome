package Genome::InstrumentData::AlignmentResult::Tophat;

use strict;
use warnings;

use File::Basename;
use Sys::Hostname;
use File::stat;
use File::Path 'rmtree';

use Genome;

class Genome::InstrumentData::AlignmentResult::Tophat {
    is => 'Genome::SoftwareResult',
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            calculate => q{
                return Genome::InstrumentData->get([$self->instrument_data_id]);
            }
        },
        reference_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            id_by => 'reference_build_id',
        },
        reference_name => {
            via => 'reference_build',
            to => 'name',
            is_mutable => 0,
            is_optional => 1
        },
        _disk_allocation => {
            is => 'Genome::Disk::Allocation',
            is_optional => 1,
            is_many => 1,
            reverse_as => 'owner'
        },
    ],

    has_input => [
        reference_build_id => {
            is => 'Number',
            doc => 'the database id of the reference to use for this alignment',
        },
        instrument_data_id => {
            is => 'Number',
            doc => 'the local database ids of the instrument data (reads) for this merged alignment',
            is_many => 1,
        },
    ],
    has_param => [
        (map
            {$_->property_name => { is => $_->data_type, doc => $_->doc, is_optional => $_->is_optional} }
            Genome::InstrumentData::AlignmentResult->__meta__->_legacy_properties(via => 'params')
        ),
        bowtie_version => {
            is => "Text",
            is_optional => 1,
            doc => "the version of bowtie that tophat should run internally"
        }
    ],

    has_constant => [
        aligner_name => { value => 'tophat', is_param=>1 },
    ],
    has_transient_optional => [
        temp_staging_directory  => {
            is => 'Text',
            doc => 'A directory to use for staging the alignment data while working.',
        },
        temp_scratch_directory  => {
            is => 'Text',
            doc => 'A directory for working files not intended to be kept',
        },
    ],
    has_calculated => [
         aligner_output_file => {
            calculate_from => ['output_dir'],
            calculate => q|
                return $output_dir .'/tophat.aligner_output';
            |,
        },
        sam_file => {
            calculate_from => ['output_dir'],
            calculate => q|
                return $output_dir .'/accepted_hits.sam';
            |,
        },
        bam_file => {
            calculate_from => ['output_dir'],
            calculate => q|
                return $output_dir .'/accepted_hits.bam';
            |,
        },
        coverage_file => {
            calculate_from => ['output_dir'],
            calculate => q|
                return $output_dir .'/coverage.wig';
            |,
        },
        junctions_file => {
            calculate_from => ['output_dir'],
            calculate => q|
                return $output_dir .'/junctions.bed';
            |,
        },
    ],
};

sub required_arch_os { 'x86_64' }

sub required_rusage {
     "-R 'select[model!=Opteron250 && type==LINUX64 && mem>16000 && tmp>150000] span[hosts=1] rusage[tmp=150000, mem=16000]' -M 16000000 -n 4";
}

sub create {
    my $class = shift;

    #This will do some locking and the like for us.
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    my $rv = eval {
        $self->status_message('Preparing directories...');
        $self->_prepare_output_directory(); #This gets a disk allocation
        my @tmp_dirs = $self->_prepare_working_directories(); #need to keep these in scope while in use

        $self->status_message('Collecting FASTQs...');
        my ($left_fastqs, $right_fastqs, $unaligned_bams) = $self->_gather_input_fastqs;

        $self->status_message('Running Tophat...');
        $self->_run_aligner($left_fastqs, $right_fastqs);

        $self->status_message('Merging and calculating stats...');
        $self->_merge_and_calculate_stats($unaligned_bams);

        $self->_promote_validated_data;
        @tmp_dirs = (); #clear out temp directories

        return 1;
    };
    if(my $error = $@) {
        $self->_cleanup;
        die $error;
    } elsif ($rv ne 1) {
        $self->error_message('Unexpected return value: ' . $rv);
        $self->_cleanup;
        die $self->error_message;
    }

    $self->status_message("Resizing the disk allocation...");
    if ($self->_disk_allocation) {
        unless ($self->_disk_allocation->reallocate) {
            $self->warning_message("Failed to reallocate disk allocation: " . $self->_disk_allocation->id);
        }
    }

    $self->status_message('All processes completed.');

    return $self;
}

sub _cleanup {
    my $self = shift;

    return unless $self->_disk_allocation;

    $self->status_message('Now deleting allocation with owner_id = ' . $self->id);
    my $allocation = $self->_disk_allocation;
    if ($allocation) {
        my $path = $allocation->absolute_path;
        unless (rmtree($path)) {
            $self->error_message("could not rmtree $path");
            return;
       }
       $allocation->deallocate;
    }
}

sub _prepare_output_directory {
    my $self = shift;

    return $self->output_dir if $self->output_dir;

    my $subdir = $self->resolve_alignment_subdirectory;
    unless ($subdir) {
        $self->error_message("failed to resolve subdirectory for instrument data.  cannot proceed.");
        die $self->error_message;
    }

    my $allocation = $self->_disk_allocation;

    unless($allocation) {
        my %allocation_parameters = (
            disk_group_name => $ENV{GENOME_DISK_GROUP_MODELS},
            allocation_path => $subdir,
            owner_class_name => $self->class,
            owner_id => $self->id,
            kilobytes_requested => $self->estimated_kb_usage(),
        );

        $allocation = Genome::Disk::Allocation->allocate(%allocation_parameters);
    }

    my $output_dir = $allocation->absolute_path;
    unless (-d $output_dir) {
        $self->error_message("Allocation path $output_dir doesn't exist!");
        die $self->error_message;
    }

    $self->output_dir($output_dir);

    return $output_dir;
}

sub resolve_alignment_subdirectory {
    my $self = shift;

    my $hostname = hostname;
    my $user = $ENV{'USER'};
    my $base_dir = sprintf("alignment-%s-%s-%s-%s", $hostname, $user, $$, $self->id);
    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'build_merged_alignments', $base_dir);
    return $directory;
}

sub _prepare_working_directories {
    my $self = shift;

    return $self->temp_staging_directory if $self->temp_staging_directory;

    my $output_dir = $self->output_dir;

    #file sizes are so large that /tmp/ would be exhausted--stage files to the allocation itself instead
    my $staging_tempdir = File::Temp->newdir(
        "staging-XXXXX",
        DIR     => $output_dir,
        CLEANUP => 1,
    );

    my $scratch_tempdir = File::Temp->newdir(
        "scratch-XXXXX",
        DIR     => $output_dir,
        CLEANUP => 1,
    );

    # fix permissions on this temp dir so others can clean it up later if need be
    chmod(0775,$staging_tempdir);
    chmod(0775,$scratch_tempdir);

    $self->temp_staging_directory($staging_tempdir->dirname);
    $self->temp_scratch_directory($scratch_tempdir->dirname);

    return ($staging_tempdir, $scratch_tempdir);
}


sub _gather_input_fastqs {
    my $self = shift;

    my @instrument_data = $self->instrument_data;
    my @left_fastqs = ();
    my @right_fastqs = ();
    my @unaligned_bams = ();
    for my $instrument_data (@instrument_data) {
        my $directory = join('/', $self->temp_scratch_directory, $instrument_data->id);
        unless(-d $directory) {
            Genome::Sys->create_directory($directory);
        }
        my @files = $instrument_data->dump_trimmed_fastq_files(
            trimmer_name => $self->trimmer_name,
            trimmer_version => $self->trimmer_version,
            trimmer_params => $self->trimmer_params,
            discard_fragments => 0,
            directory => $directory,
        );
        if( scalar( @files ) == 3 ) {
            $self->warning_message( "Discarding fragmented reads since Tophat cannot support a mix of single-end and paired-end reads" );
        }
        elsif( scalar( @files ) != 1 && scalar( @files ) != 2 ) {
            $self->error_message("Unexpected number of FASTQ files (got ".scalar(@files)."):\n\t".join("\n\t", @files));
            die $self->error_message;
        }
        push @left_fastqs, $files[0];
        push @right_fastqs, $files[1] if( scalar(@files) == 2 );

        #needed for post-alignment statistics
        my $unaligned_bam = $directory .'/s_'. $instrument_data->subset_name .'_sequence.bam';

        if($self->trimmer_name or not $instrument_data->bam_path) {
            my %params = (
                fastq => $files[0],
                output => $unaligned_bam,
                quality_format => 'Standard',
                sample_name => $instrument_data->sample_name,
                library_name => $instrument_data->library_name,
                log_file => $directory .'/s_'. $instrument_data->subset_name .'_sequence.log',
                platform => 'illumina',
                platform_unit => $instrument_data->flow_cell_id .'.'. $instrument_data->subset_name,
                read_group_name => $instrument_data->id,
                sort_order => 'queryname',
                use_version => $self->picard_version,
                maximum_memory => 12,
                maximum_permgen_memory => 256,
                max_records_in_ram => 3000000,
            );
            # If these are paired-end reads, then don't forget the second fastq file
            $params{fastq2} = $files[1] if( scalar(@files) == 2 );

            unless (Genome::Model::Tools::Picard::FastqToSam->execute( %params )) {
                die $self->error_message('Failed to create per lane, unaligned BAM file: '
                .$self->temp_scratch_directory .'/s_'. $instrument_data->subset_name .'_sequence.bam');
            }
        }
        else {
            my $existing_bam = $instrument_data->bam_path;
            Genome::Sys->create_symlink($existing_bam, $unaligned_bam);
        }

        push @unaligned_bams, $unaligned_bam;
    }

    return (\@left_fastqs, \@right_fastqs, \@unaligned_bams);
}

sub _run_aligner {
    my $self = shift;
    my $left_fastq_files = shift;
    my $right_fastq_files = shift;

    my $reference_build = $self->reference_build;
    my $reference_path = $reference_build->full_consensus_path('bowtie');

    my $read_1_fastq_list = join(',',@$left_fastq_files);
    my $read_2_fastq_list = join(',',@$right_fastq_files);

    my ($sum_insert_sizes, $sum_insert_size_std_dev, $sum_read_length, $reads);
    for my $instrument_data ($self->instrument_data) {
        my $median_insert_size   = $instrument_data->resolve_median_insert_size;
        my $sd_above_insert_size = $instrument_data->resolve_sd_insert_size;
        # Use the number of reads to somewhat normalize the averages we will calculate later
        # This is not the best approach, any ideas?
        my $clusters = int($instrument_data->read_count / 2);
        if ($median_insert_size && $sd_above_insert_size) {
            $sum_insert_sizes += ($median_insert_size * $clusters);
            $sum_insert_size_std_dev += ($sd_above_insert_size * $clusters);
        } else {
            # These seem like reasonable default values given most libraries are 300-350bp
            $sum_insert_sizes += (300 * $clusters);
            $sum_insert_size_std_dev += (20 * $clusters);
        }
        # TODO: This could be skewed if Read 2 is conconcatanated or trimmed reads are used
        $sum_read_length += ($instrument_data->read_length * $clusters);
        $reads += $clusters;
    }

    unless ($reads) {
        die $self->error_message('Failed to calculate the number of reads across all lanes');
    }
    unless ($sum_insert_sizes) {
        die $self->error_message('Failed to calculate the sum of insert sizes across all lanes');
    }
    unless ($sum_insert_size_std_dev) {
        die $self->error_message('Failed to calculate the sum of insert size standard deviation across all lanes');
    }
    unless ($sum_read_length) {
        die $self->error_message('Failed to calculate the sum of read lengths across all lanes');
    }
    my $avg_read_length = int($sum_read_length / $reads);

    # The inner-insert size should be the predicted external-insert size(300) minus the read lengths(2x100=200). Example: 300-200=100
    # If there is an insert size, it's paired-end mutliply read_length by 2
    my $insert_size = int( $sum_insert_sizes / $reads ) - ($avg_read_length * 2);

    unless (defined($insert_size)) {
        die $self->error_message('Failed to get insert size with '. $reads .' reads and a sum insert size of '. $sum_insert_sizes);
    }
    # TODO: averaging the standard deviations does not seem statisticly sound
    my $insert_size_std_dev = int( $sum_insert_size_std_dev / $reads );
    unless (defined($insert_size_std_dev)) {
        die $self->error_message('Failed to get insert size with '. $reads .' $reads and a sum insert size standard deviation of '. $sum_insert_size_std_dev);
    }

    my %params = (
        reference_path => $reference_path,
        read_1_fastq_list => $read_1_fastq_list,
        insert_size => $insert_size,
        insert_std_dev => $insert_size_std_dev,
        aligner_params => $self->aligner_params,
        alignment_directory => $self->temp_staging_directory,
        use_version => $self->aligner_version,
        bowtie_version => $self->bowtie_version,
    );

    # If these are paired-end reads, then don't forget the second fastq file
    $params{read_2_fastq_list} = $read_2_fastq_list if( scalar( @$right_fastq_files ) > 0 );

    my $tophat_cmd = Genome::Model::Tools::Tophat::AlignReads->create(%params);

    eval {
        unless($tophat_cmd->execute) {
            die 'Execute did not return a true value.';
        }
    };
    if($@) {
        my $error = $@ || '_run_aligner failed.';
        $self->error_message('Failed to execute tophat command.');

        #Try to record a copy of the aligner logs before they get blasted in cleanup
        my $aligner_log = $self->temp_staging_directory . '/tophat.aligner_output';
        if(-e $aligner_log) {
            my $log_text = Genome::Sys->read_file($aligner_log);
            $self->status_message("Aligner log:\n" . $log_text);
        }

        die($error);
    }

    return 1;
}

sub _merge_and_calculate_stats {
    my $self = shift;
    my $unaligned_bams = shift;

    my $tmp_all_reads_bam_file = $self->temp_scratch_directory . '/all_fastq_reads.bam';
    unless (Genome::Model::Tools::Picard::MergeSamFiles->execute(
        input_files => $unaligned_bams,
        output_file => $tmp_all_reads_bam_file,
        maximum_memory => 12,
        maximum_permgen_memory => 256,
        sort_order => 'queryname',
        use_version => $self->picard_version,
    )) {
        die $self->error_message('Failed to merge unaligned BAM files!');
    }

    # queryname sort the aligned BAM file
    my $tmp_aligned_bam_file = $self->temp_scratch_directory . '/accepted_hits_queryname_sort.bam';
    unless (Genome::Model::Tools::Picard::SortSam->execute(
        sort_order => 'queryname',
        input_file => $self->temp_staging_directory .'/accepted_hits.bam',
        output_file => $tmp_aligned_bam_file,
        max_records_in_ram => 3000000,
        maximum_memory => 12,
        maximum_permgen_memory => 256,
        temp_directory => $self->temp_scratch_directory,
        use_version => $self->picard_version,
    )) {
        die $self->error_message('Failed to queryname sort the aligned BAM file!');
    }

    # Find unaligned reads and merge with aligned while calculating basic alignment metrics
    my $tmp_unaligned_bam_file = $self->temp_scratch_directory . '/unaligned_reads.bam';
    my $unaligned_bam_file = $self->temp_staging_directory . '/unaligned_reads.bam';
    my $alignment_stats_file = $self->temp_staging_directory . '/alignment_stats.txt';
    my $cmd = "genome-perl5.10 -S gmt bio-samtools tophat-alignment-stats --aligned-bam-file=$tmp_aligned_bam_file --all-reads-bam-file=$tmp_all_reads_bam_file --unaligned-bam-file=$tmp_unaligned_bam_file --alignment-stats-file=$alignment_stats_file";
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$tmp_aligned_bam_file,$tmp_all_reads_bam_file],
        output_files => [$tmp_unaligned_bam_file,$alignment_stats_file],
    );
    unlink($tmp_all_reads_bam_file);
    unlink($tmp_aligned_bam_file);

    # Scrub the eland alignment information from the unaligned reads bam
    my @eland_attributes = qw(H0 H1 H2 XD SM AS);
    my $reset_command = Genome::Model::Tools::Picard::RevertSam->create(
        remove_alignment_information => 1,
        remove_duplicate_information => 1,
        restore_original_qualities => 1,
        input_file => $tmp_unaligned_bam_file,
        output_file => $unaligned_bam_file,
        use_version => $self->picard_version,
        attribute_to_clear => \@eland_attributes,
    );

    unless ($reset_command->execute == 1) {
        die $self->error_message("Failed to execute gmt picard revert-sam on the unaligned bam file");
    }

    return 1;
}

sub _promote_validated_data {
    my $self = shift;

    my $staging_dir = $self->temp_staging_directory;
    my $output_dir  = $self->output_dir;

    $self->status_message("Now de-staging data from $staging_dir into $output_dir");

    for my $staged_file (glob("$staging_dir/*")) {
        my $destination = $staged_file;
        $destination =~ s/$staging_dir/$output_dir/;
        rename($staged_file, $destination);
    }

    chmod 02775, $output_dir;
    for my $subdir (grep { -d $_  } glob("$output_dir/*")) {
        chmod 02775, $subdir;
    }

    # Make everything in here read-only
    for my $file (grep { -f $_  } glob("$output_dir/*")) {
        chmod 0444, $file;
    }

    $self->status_message("Files in $output_dir: \n" . join "\n", glob($output_dir . "/*"));

    return $output_dir;
}

sub estimated_kb_usage {
    my $self = shift;
    my @instrument_data = $self->instrument_data;
    #50GB per lane
    #TODO Take into account the size of each instrument data rather than assuming a constant size
    my $estimate_from_instrument_data = scalar(@instrument_data) * 52428800;

    return $estimate_from_instrument_data;
}

1;
