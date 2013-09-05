package Genome::Model::MetagenomicShotgun::Build::ExtractFrom;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun::Build::ExtractFrom {
    is => 'Command::V2',
    has_input => [
        input_build => {
            is => 'Genome::Model::Build::MetagenomicShotgun',
            is_many => 1,
            doc => 'The MetaShot build to work with.',
        },
        sub_model_label => { 
            is => 'Text',
            valid_values => [ Genome::Model::MetagenomicShotgun->sub_model_labels ],
        },
        type => {
            is => 'Text',
            value => [qw/ aligned unaligned /],
            doc => 'Type of reads to extract.',
        },
    ],
    has_output => [
        build => {
            is => 'Genome::Model::Build::MetagenomicShotgun',
            calculate_from => ['input_build'],
            calculate => sub{ return $_[0]; },
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'Extracted instrument data.'
        },
    ],
};

sub execute {
    my $self = shift;

    my $sub_build = $self->build->model->last_complete_build_for_sub_model($self->sub_model_label);
    return if not $sub_build;

    my $extraction_type = $self->type;

    if ( not defined $sub_build ){
        $self->status_message("No previous build provided, skipping $extraction_type data extraction");
        return;
    }
    $self->status_message("Extracting $extraction_type reads from ".$sub_build->__display_name__);
    my @assignments = $sub_build->instrument_data_inputs;
    my %assignments_and_extracted_instrument_data;
    for my $assignment (@assignments) {
        my @alignment_results = $sub_build->alignment_results_for_instrument_data($assignment->value);
        if (@alignment_results > 1) {
            die $self->error_message( "multiple alignment_results found for instrument data assignment: " . $assignment->__display_name__);
        }
        if (@alignment_results == 0) {
            die $self->error_message( "no alignment_results found for instrument data assignment: " . $assignment->__display_name__);
        }
        $self->status_message("processing instrument data assignment ".$assignment->__display_name__." for unaligned reads import");

        my $alignment_result = $alignment_results[0];
        my @extracted_instrument_data_for_alignment_result = $self->_extract_data_from_alignment_result(
            $alignment_result, $extraction_type, $self->build->model->filter_duplicates,
        );

        $assignments_and_extracted_instrument_data{ $assignment->id } = \@extracted_instrument_data_for_alignment_result;
        #push @extracted_instrument_data, @extracted_instrument_data_for_alignment_result;
    }

    unless (keys(%assignments_and_extracted_instrument_data) == @assignments) {
        #unless (@extracted_instrument_data == @assignments) {
        die $self->error_message("The count of extracted instrument data sets does not match screened instrument data assignments.");
    }
    $self->instrument_data([map { @$_ } values %assignments_and_extracted_instrument_data]);

    return 1;
}

sub _extract_data_from_alignment_result{
    my ($self, $alignment, $extraction_type, $filter_duplicates) = @_;

    my $instrument_data = $alignment->instrument_data;
    my $lane = $instrument_data->lane;
    my $instrument_data_id = $instrument_data->id;

    my $aligner_name = $alignment->aligner_name; #we need to do special handling for rtg mapx alignments, which are protein space and do not produce bam files

    my $dir = $alignment->output_dir;
    my $alignment_output;
    if ($aligner_name eq 'rtg mapx'){
        unless ($extraction_type eq 'aligned'){
            die $self->error_message("metagenomic shotgun pipeline should only extract aligned reads from mapx output");
        }
        $alignment_output = $dir . '/alignments.txt';
    }else{
        $alignment_output = $dir . '/all_sequences.bam';
    }
    unless (-e $alignment_output) {
        die $self->error_message("Failed to find expected alignment output file $alignment_output\n");
    }

    # come up with the import "original_data_path", used to find existing data, and when uploading new data

    my $tmp_dir = "/tmp/extracted_reads";
    $tmp_dir .= "/".$alignment->id;
    my @instrument_data;
    my $subdir = $extraction_type;
    $subdir =~ s/\s/_/g;
    if ($filter_duplicates){
        $subdir .= "/deduplicated";
    }
    $self->status_message("Preparing imported instrument data for import path $tmp_dir/$subdir");

    my $forward_basename = "s_$lane" . "_1_sequence.txt";
    my $reverse_basename = "s_$lane" . "_2_sequence.txt";
    my $fragment_basename = "s_$lane" . "_sequence.txt";

    my $expected_data_path0 = "$tmp_dir/$subdir/$fragment_basename";
    my $expected_data_path1 = "$tmp_dir/$subdir/$forward_basename";
    my $expected_data_path2 = "$tmp_dir/$subdir/$reverse_basename";

    my $expected_se_path = $expected_data_path0;
    my $expected_pe_path = $expected_data_path1 . ',' . $expected_data_path2;

    # get any pre-existing imported instrument data for the given se_path and pe_path

    my ($se_lock, $pe_lock);

    $self->status_message("Checking for previously imported extracted reads from: $tmp_dir/$subdir");
    my $se_instdata = Genome::InstrumentData::Imported->get(original_data_path => $expected_se_path);
    if ($se_instdata) {
        $self->status_message("imported instrument data already found for path $expected_se_path, skipping");
    }
    else {
        $se_lock = $self->lock($instrument_data_id, $fragment_basename);
        unless ($se_lock) {
            die $self->error_message("Failed to lock $expected_se_path.");
        }
    }

    my $pe_instdata;
    if ( $instrument_data->is_paired_end ){
        $pe_instdata = Genome::InstrumentData::Imported->get(original_data_path => $expected_pe_path);
        if ($pe_instdata) {
            $self->status_message("imported instrument data already found for path $expected_pe_path, skipping");
        }
        else  {
            $pe_lock = $self->lock($instrument_data_id, "s_${lane}_1-2_sequence.txt");
            unless ($pe_lock) {
                die $self->error_message("Failed to lock $expected_pe_path.");
            }
        }
    }

    if (!$se_lock and !$pe_lock) {
        $self->status_message("skipping read processing since all data is already processed and uploaded");
        return grep { defined $_ } ($se_instdata, $pe_instdata);
    }

    # extract what we did not already find...

    my $working_dir = Genome::Sys->create_temp_directory();
    my $forward_unaligned_data_path     = "$working_dir/$instrument_data_id/$forward_basename";
    my $reverse_unaligned_data_path     = "$working_dir/$instrument_data_id/$reverse_basename";
    my $fragment_unaligned_data_path    = "$working_dir/$instrument_data_id/$fragment_basename";

    my $cmd;
    if ($aligner_name eq 'rtg mapx'){
        $cmd = "genome-perl -S genome model metagenomic-shotgun mapx-alignment-to-aligned-fastq --mapx-alignment $alignment_output --output-directory $working_dir --instrument-data $instrument_data_id";
    }else{
        $cmd = "genome-perl -S gmt sam bam-to-unaligned-fastq --bam-file $alignment_output --output-directory $working_dir --ignore-bitflags"; #add ignore bitflags here because some of the aligners used in this pipeline produce untrustworthy flag information
        if ($extraction_type eq 'aligned'){
            $cmd.=" --print-aligned";
        }
        elsif($extraction_type eq 'unaligned'){
            #default
        }
        else{
            die $self->error_message("Unhandled extraction_type $extraction_type");
        }
        if ($filter_duplicates){
            $cmd.=" --filter-duplicates";
        }
    }

    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        die "Failed to extract unaligned reads: $@";
    }

    my @expected_output_fastqs = ( $instrument_data->is_paired_end )
        ?  ($forward_unaligned_data_path, $reverse_unaligned_data_path, $fragment_unaligned_data_path)
        :  ($fragment_unaligned_data_path);

    my @missing = grep {! -e $_} grep { defined($_) and length($_) } @expected_output_fastqs;
    if (@missing){
        die $self->error_message(join(", ", @missing)." unaligned files missing after bam extraction");
    }
    $self->status_message("Extracted unaligned reads from bam file (@expected_output_fastqs)");

    # upload

    $self->status_message("uploading new instrument data from the post-processed unaligned reads...");
    my @properties_from_prior = qw/
        run_name
        sequencing_platform
        median_insert_size
        sd_above_insert_size
        library_name
        sample_name
    /;
    my @errors;
    my %properties_from_prior;
    for my $property_name (@properties_from_prior) {
        my $value = $instrument_data->$property_name;
        no warnings;
        $self->status_message("Value for $property_name is $value");
        $properties_from_prior{$property_name} = $value;
    }
    $properties_from_prior{subset_name} = $instrument_data->lane;

    my @source_paths = ( $instrument_data->is_paired_end )
        ?  ($fragment_unaligned_data_path, "$forward_unaligned_data_path,$reverse_unaligned_data_path" )
        :  ($fragment_unaligned_data_path);

    for my $source_data_files (@source_paths) {
        $self->status_message("Attempting to upload $source_data_files...");
        my $original_data_path;
        if ($source_data_files =~ /,/){
            $properties_from_prior{is_paired_end} = 1;
            $original_data_path = $expected_pe_path;
        }
        else {
            $properties_from_prior{is_paired_end} = 0;
            $original_data_path = $expected_se_path;
        }
        my %params = (
            %properties_from_prior,
            source_data_files => $source_data_files,
            import_format => 'sanger fastq',
        );
        $self->status_message("importing fastq with the following params:" . Data::Dumper::Dumper(\%params));

        my $command = Genome::InstrumentData::Command::Import::Fastq->create(%params);
        unless ($command) {
            $self->error_message( "Couldn't create command to import unaligned fastq instrument data!");
        };
        my $result = $command->execute();
        unless ($result) {
            die $self->error_message( "Error importing data from $source_data_files! " . Genome::InstrumentData::Command::Import::Fastq->error_message() );
        }
        $self->status_message("committing newly created imported instrument data");

        my $new_instrument_data = Genome::InstrumentData::Imported->get($command->generated_instrument_data_id);
        unless ($new_instrument_data) {
            die $self->error_message( "Failed to find new instrument data $source_data_files!");
        }

        $new_instrument_data->original_data_path($original_data_path);

        UR::Context->commit();

        if ($new_instrument_data->__changes__) {
            die "unsaved changes present on instrument data $new_instrument_data->{id} from $original_data_path!!!";
        }
        if ( $se_lock ) {
            $self->status_message("Attempting to remove lock on $se_lock...");
            unless(Genome::Sys->unlock_resource(resource_lock => $se_lock)) {
                die $self->error_message("Failed to unlock $se_lock.");
            }
            undef($se_lock);
        }
        if ( $pe_lock ) {
            $self->status_message("Attempting to remove lock on $pe_lock...");
            unless(Genome::Sys->unlock_resource(resource_lock => $pe_lock)) {
                die $self->error_message("Failed to unlock $pe_lock.");
            }
            undef($pe_lock);
        }

        push @instrument_data, $new_instrument_data;
    }

    if ( not @instrument_data ) {
        die $self->error_message("Error processing unaligned reads!");
    }

    return @instrument_data;
}

sub lock {
    my $self = shift;
    my @parts = @_;
    my $lock_key = join('_', @parts);
    $self->debug_message("Creating lock on $lock_key...");
    my $resource_lock = File::Spec->join($ENV{GENOME_LOCK_DIR}, $lock_key);
    my $lock = Genome::Sys->lock_resource(
        resource_lock => $resource_lock,
        max_try => 2,
    );
    return $lock;
};

1;

