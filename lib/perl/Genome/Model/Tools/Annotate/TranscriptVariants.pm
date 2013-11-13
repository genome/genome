package Genome::Model::Tools::Annotate::TranscriptVariants;

use strict;
use warnings;

use Genome; 

use Command;
use Data::Dumper;
use IO::File;
use Benchmark;
use DateTime;
use Sys::Hostname;
use Cwd;
use File::Basename;
use File::Copy;
use Genome::Info::IUB;
use Genome::Info::UCSCConservation;

class Genome::Model::Tools::Annotate::TranscriptVariants {
    is => 'Command',
    has => [ 
        variant_file => {
            is => 'FilePath',   
            is_input => 1,
            is_optional => 1,
            doc => "File of variants. Tab separated columns: chromosome_name start stop reference variant.",
        },
        variant_bed_file => {
            is => 'FilePath',   
            is_input => 1,
            is_optional => 1,
            doc => "File of variants in BED format.  ", #TODO: figure out how one specifies a variant file in bed, then list the format
        },
        output_file => {
            is => 'Text',
            is_input => 1,
            is_output=> 1,
            doc => "Store annotation in the specified file. Defaults to STDOUT if no file is supplied.",
            default => "STDOUT",
        },
        _version_subclass_name => {
            is => 'Text',# is_mutable => 0, 
            calculate_from => ['use_version'],
            calculate => q( return $self->_annotator_subclass_for_version($use_version) ),
        },
        available_versions => {
            is_calculated => 1,
            is_many => 1,
            calculate => q( my $path = Genome::Model::Tools::Annotate::TranscriptVariants->base_dir;
                            my @versions = glob($path.'/Version*.pm');
                            @versions = map { m/Version(\w+).pm/; $1 } @versions;
                            return @versions; ),
        },
    ],
    has_optional => [
        use_version => {
            is_input => 1,
            is => 'Text',
            default_value => __PACKAGE__->default_annotator_version,
            doc => 'Annotator version to use',
            valid_values => [0..4],#__PACKAGE__->available_versions,
        },
        # IO Params
        _is_parallel => {
            is => 'Boolean',
            is_input => 1,
            default => 0,
        },
        no_headers => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'Exclude headers in report output',
        },
        extra_columns => {
            is => 'Text',
            is_input => 1,
            is_optional => 1,
            is_input => 1,
            doc => "A comma delimited list of any extra columns that exist after the expected 5 in the input. Use this option if it is desired to preserve additional columns from the input file, which will then appear in output.Preserved columns must be contiguous and in order as they appear in the infile after the mandatory input columns. Any desired naming or number of columns can be specified so long as it does not exceed the actual number of columns in the file."
        },
        # Transcript params
        annotation_filter => {
            is => 'String',
            is_optional => 1,
            is_input => 1,
            default => 'gene',
            doc => 'The type of filtering to use on the annotation results. There are currently 3 valid options:
                    "none" -- This returns all possible transcript annotations for a variant. All transcript status are allowed including "unknown" status.
                    "gene" -- This returns the top transcript annotation per gene. This is the default behavior.
                    "top" -- This returns the top priority annotation for all genes. One variant in, one annotation out.',
        },
        reference_transcripts => {
            is => 'String',
            is_input => 1,
            is_optional => 1, 
            doc => 'provide name/version number of the reference transcripts set you would like to use ("NCBI-human.combined-annotation/0").  Leaving off the version number will grab the latest version for the transcript set, and leaving off this option and build_id will default to using the latest combined annotation transcript set. Use this or --build-id to specify a non-default annoatation db (not both). See full help output for a list of available reference transcripts.'
        },
        data_directory => {
            is => 'String',
            is_optional => 1,
            is_input => 1,
            is_deprecated => 1,
            doc => 'Alternate method to specify imported annotation data used in annotation.  This option allows a directory w/o supporting model and build, not reccomended except for testing purposes',
        },
        build_id =>{
            is => "Number",
            is_optional => 1,
            is_input => 1,
            doc => 'build id for the imported annotation model to grab transcripts to annotate from.  Use this or --reference-transcripts to specify a non-default annotation db (not both)',
        },
        build => {
            is => "Genome::Model::Build",
            id_by => 'build_id',
            is_optional => 1, 
        },
        extra_details => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enabling this flag produces an additional four columns: flank_annotation_distance_to_transcript, intron_annotation_substructure_ordinal, intron_annotation_substructure_size, and intron_annotation_substructure_position',
        },
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
        },
        accept_reference_IUB_codes => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to allow the annotator to accept variants that have IUB codes in the reference position.  The IUB code will be translated to a single base represented by the code.  A will be used whenever possible, followed by C, G, T in that order.  This option is only usable for annotator versions 2 or higher',
        },
        # Performance options
        cache_annotation_data_directory => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            is_input => 1,
            is_deprecated => 1,
            doc => 'enable this flag to cache annotation data locally, useful if annotation is being run repeatedly on a pipeline',
        },
        benchmark => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            is_input => 1,
            doc => 'if set, run times are displayed as status messages after certain steps are completed (annotation of whole chromosomes, caching times, etc)',
        },
        check_variants => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            is_input => 1,
            doc => 'If set, the annotator will check all variant reference sequences against the respective species reference before annotating. Annotation is skipped for those variants with mismatches.',
        },
        get_frame_shift_sequence => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'If set, the annotator will get all coding region sequence after a frame shift',
        },
    ], 
    has_param => [
        lsf_resource => {
            is => 'Text',
            default => "-M 7000000 -R 'select[type==LINUX64 && mem>7000 && tmp>10240] rusage[mem=7000]'",
        },
        lsf_queue => {
            is => 'Text',
            default => 'long',
        },
    ],
};

sub is_sub_command_delegator { 0 }

############################################################

sub help_synopsis { 
    return <<EOS
gmt annotate transcript-variants --variant-file variants.tsv --output-file transcript-changes.tsv
EOS
}

sub default_annotator_version {
    return 4;
}

sub help_detail {
    my $self = shift;

    #Generate the currently available annotation models on the fly
    my @currently_available_models = Genome::Model->get(subclass_name => 'Genome::Model::ImportedAnnotation');
    my $currently_available_builds; 
    foreach my $model (@currently_available_models) {
        next unless $model;
        foreach my $build ($model->builds) {
            if($build and $build->status eq 'Succeeded' and $build->name =~ /ensembl/) {  #probably implicit in the loops, but in case we get undefs in our list
                 $currently_available_builds .= "\t" . $model->name . "/" . $build->version . "\n" if $build->version !~ /old/ and $model->name and $build->version;
            }
        }
    }

    my $available_versions = join (', ', $self->available_versions);
    my $version_docs = join("\n", map { "version $_: " . $self->_annotator_subclass_for_version($_)->__meta__->doc } $self->available_versions);

    return <<EOS 
This launches the variant annotator.  It takes genome sequence variants and outputs transcript variants, 
with details on the gravity of the change to the transcript.

The current version presumes that the variants are human, and that positions are relative to Hs36.  The transcript data 
set is a mix of Ensembl 45 and Genbank transcripts.  Work is in progress to support newer transcript sets, and 
variants from a different reference.

The variant (if it is a SNP) can be an IUB code, in which case every possible variant base will be annotated.

INPUT COLUMNS (TAB SEPARATED)
chromosome_name start stop reference variant

The mutation type will be inferred based upon start, stop, reference, and variant alleles.

Any number of additional columns may be in the input following these columns, but they will be disregarded.

OUTPUT COLUMNS (COMMMA SEPARATED)
chromosome_name start stop reference variant type gene_name transcript_name transcript_species transcript_source transcript_version strand transcript_status trv_type c_position amino_acid_change ucsc_cons domain all_domains deletion_substructures transcript_error

The --use-version parameter is used to choose which version of the trasncript variant annotator to run.  The available versions are: $available_versions
$version_docs

COMMONLY USED REFERENCE TRANSCRIPTS WITH VERSIONS
$currently_available_builds
EOS
}

############################################################


sub _annotator_subclass_for_version {
    return 'Genome::Model::Tools::Annotate::TranscriptVariants::Version' . $_[1];
}

sub variant_attributes {
    return (qw/ chromosome_name start stop reference variant /);
}
sub variant_output_attributes {
    return (qw/ type /);
}
sub transcript_attributes {
    my $self = shift;
    my @attrs = qw( gene_name transcript_name transcript_species transcript_source
                    transcript_version strand transcript_status trv_type c_position
                    amino_acid_change ucsc_cons domain all_domains deletion_substructures
                    transcript_error);
    if ($self->use_version >= 3) {
        push @attrs, qw( default_gene_name gene_name_source ensembl_gene_id);
    }
    if ($self->extra_details) {
        push @attrs, qw( flank_annotation_distance_to_transcript
                         intron_annotation_substructure_ordinal intron_annotation_substructure_size
                         intron_annotation_substructure_position );
    }
    return @attrs;
}


# Figures out what the 'type' of this variant should be (snp, dnp, ins, del) based upon
# the start, stop, reference, and variant
# Takes in a variant hash, returns the type
sub infer_variant_type {
    my ($self, $variant) = @_;

    if(( (!$variant->{reference})||($variant->{reference} eq '0')||($variant->{reference} eq '-')) &&
        ((!$variant->{variant})||($variant->{variant} eq '0')||($variant->{variant} eq '-'))){
        $self->error_message("Could not determine variant type from variant:");
        $self->error_message(Dumper($variant));
        die;
    }

    # If the start and stop are the same, and ref and variant are defined its a SNP
    if (($variant->{stop} == $variant->{start})&&
        ($variant->{reference} ne '-')&&($variant->{reference} ne '0')&&
        ($variant->{variant} ne '-')&&($variant->{variant} ne '0')) {
        return 'SNP';
    # If start and stop are 1 off, and ref and variant are defined its a DNP
    } elsif (($variant->{stop} - $variant->{start} == 1)&&
             ($variant->{reference})&&($variant->{reference} ne '-')&&($variant->{reference} ne '0')&&
             ($variant->{variant})&&($variant->{variant} ne '-')&&($variant->{variant} ne '0')) {
        return 'DNP';
    # If reference is a dash, we have an insertion
    } elsif (($variant->{reference} eq '-')||($variant->{reference} eq '0')) {
        return 'INS';
    } elsif (($variant->{variant} eq '-')||($variant->{variant} eq '0')) {
        return 'DEL';
    } else {
        $self->error_message("Could not determine variant type from variant:");
        $self->error_message(Dumper($variant));
        if(defined($self->variant_bed_file)){
            $self->error_message("Perhaps you meant to use --variant-file instead of --variant-bed-file?");
        }
        die;
    }
}


# This is for version 0 of the annotator; the original one that operated on the transcript window
# It required different args to create than the new one does
sub _create_old_annotator {
    my($self, $annotator_version_subclass) = @_;

    my $full_version = $self->build->version;
    my ($version) = $full_version =~ /^\d+_(\d+)[a-z]/;
    my %ucsc_versions = Genome::Info::UCSCConservation->ucsc_conservation_directories;

    my $annotator = $annotator_version_subclass->create(
                        check_variants => $self->check_variants,
                        get_frame_shift_sequence => $self->get_frame_shift_sequence,
                        ucsc_conservation_directory => $ucsc_versions{$version},
                        annotation_build_version => $self->build->version,
                        build => $self->build,
                     );
    return $annotator;
}

sub _convert_bed_file {
    my $self = shift; 
    my $converted_bed_file = Genome::Sys->create_temp_file_path();
    my $bed_converter_class = $self->_version_subclass_name. "::BedToAnnotation";
    my %bed_converter_params = (snv_file => $self->variant_bed_file, output => $converted_bed_file);
    if ($self->extra_columns and $self->use_version >= 3){
        $bed_converter_params{'extra_columns'} = $self->extra_columns;
    }
    $bed_converter_class->execute(%bed_converter_params) || ($self->error_message("Could not convert BED file to annotator format") and return);
    $self->variant_file($converted_bed_file);
}

# debugging help
my $we_are_done_flag;
my $last_variant_annotated;
END {
    if (defined $we_are_done_flag and ! $we_are_done_flag) {
        print STDERR "\n\nThe last variant we worked on is\n",Data::Dumper::Dumper($last_variant_annotated),"\n\n";
    }
};

sub _validate_parameters {
    my $self = shift;

    unless($self->variant_file xor $self->variant_bed_file){
        $self->error_message("Please specify either a --variant-file or a --variant-bed-file");
        return;
    }

    if ($self->variant_bed_file) {
        unless (-s $self->variant_bed_file) {
            die $self->error_message("Variant bed file has no size, exiting");
        }

        $self->_convert_bed_file;

        unless(-s $self->variant_file){
            die $self->error_message("After converting to annotation format, the variant file has no size, exiting");
        }
    } else {
        unless(-s $self->variant_file){
            $self->error_message("Variant file has no size, exiting");
            die;
        }
    }

    if (defined $self->data_directory) {
        $self->error_message("Due to a recent change to the annotation data file format, allowing " .
            "user-specified data directories has been deprecated. Specifying a data directory containing " .
            "data that does not meet the new format will result in some wonky errors, so this is " .
            "one way to avoid that mess. If you have questions, contact apipe.");
        die;
    }

    #check to see if reference_transcripts set name and build_id given
    unless ($self->build xor $self->reference_transcripts){
        $self->error_message("Please provide EITHER a build id OR a reference transcript set name, but not both");
        return;
    }

    unless($self->build) {
        my $ref = $self->reference_transcripts;
        my ($name, $version) = split(/\//, $ref); # For now, version is ignored since only v2 is usable
                                        # This will need to be changed when other versions are available

        my $model = Genome::Model->get(name => $name);
        unless ($model){
            $self->error_message("couldn't get reference transcripts set for $name");
            return;
        }

        my $build = $model->build_by_version($version);
        unless ($build){
            $self->error_message("couldn't get build from reference transcripts set $name");
            return;
        }
        $self->build($build);
    }

    if ($self->build and $self->cache_annotation_data_directory) {
        $self->cache_annotation_data_directory(0);
        $self->status_message("--cache-annotation-data-directory is currently disabled.  Operating from the annotation data directory instead.");
        $self->_notify_cache_attempt;
    }

    my $annotation_filter = $self->annotation_filter( lc $self->annotation_filter );
    unless ($self->_annotation_method_for_annotation_filter) {
       $self->error_message("Unknown annotation_filter value: " . $annotation_filter);
        return;
    }

    if($self->accept_reference_IUB_codes and $self->use_version < 2) {
        $self->warning_message("accept-reference-IUB-codes is only available in version 2 or better");
    }

    return 1;
}

sub _print_starting_message {
    my $self = shift;

    # Useful information for debugging...
    my ($date, $time) = split(' ',$self->__context__->now());
    my $host = hostname;
    $self->status_message("Executing on host $host on $date at $time");
}

# generate an iterator for the input list of variants
sub _create_variant_reader {
    my $self = shift;

    # preserve additional columns from input if desired
    my @columns = (($self->variant_attributes), $self->get_extra_columns);
    my $variant_svr = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->variant_file,
        headers => \@columns,
        separator => "\t",
        is_regex => 1,
        ignore_extra_columns => 1,
    );
    unless ($variant_svr) {
        $self->error_message("error opening file " . $self->variant_file);
        return;
    }
    return $variant_svr;
}

sub _setup_report_fh {
    my $self = shift;

    # establish the output handle for the transcript variants
    my $output_fh;
    my $output_file = $self->output_file;
    my $temp_output_file;
    if ($self->output_file =~ /STDOUT/i) {
        $output_fh = 'STDOUT';
    }
    else {
        my ($output_file_basename) = File::Basename::fileparse($output_file);
        ($output_fh, $temp_output_file) = File::Temp::tempfile(
            "$output_file_basename-XXXXXX",
            DIR => '/tmp/',
            UNLINK => 1,
            CLEANUP => 1,
        );
        chmod(0664, $temp_output_file);
    }
    $self->_transcript_report_fh($output_fh);
    $self->_temp_output_file($temp_output_file);

    # omit headers as necessary
    $output_fh->print( join("\t", $self->transcript_report_headers), "\n" ) unless $self->no_headers;

    return $output_fh;
}

sub _close_report_fh {
    my $self = shift;

    my $output_fh = $self->_transcript_report_fh;
    $output_fh->close unless $output_fh eq 'STDOUT';
    if ($self->_temp_output_file){
        unless (move($self->_temp_output_file, $self->output_file)) {
            $self->error_message("Failed to mv results at ".$self->_temp_output_file." to final location at ".$self->output_file.": $!");
            return 0;
        }
        $output_fh->close unless $output_fh eq 'STDOUT';
    }
}

sub _create_annotator {
    my $self = shift;

    # Initialize the annotator object
    if (! defined ($self->use_version)) {
        # We're going to use the "default" annotator version for these old processing profiles
        # that don't have a --transcript-variant-annotator-version defined
        $self->use_version($self->default_annotator_version);
    }
    my $annotator_version_subclass = $self->_version_subclass_name;
    my $annotator;
    eval {
        if ($self->use_version == 0) {
            $annotator = $self->_create_old_annotator($annotator_version_subclass);
        } 
        elsif ($self->use_version < 3) {
            my @directories = $self->build->determine_merged_data_directory($self->cache_annotation_data_directory);
            $annotator = $annotator_version_subclass->create(
                data_directory => \@directories,
                check_variants => $self->check_variants,
                get_frame_shift_sequence => $self->get_frame_shift_sequence,
            );
        }
        else {
            my @directories = $self->build->determine_merged_data_directory($self->cache_annotation_data_directory);
            $annotator = $annotator_version_subclass->create(
                data_directory => \@directories,
                check_variants => $self->check_variants,
                get_frame_shift_sequence => $self->get_frame_shift_sequence,
                reference_sequence_id => $self->build->reference_sequence_id,
            );
        }
    };
    unless ($annotator){
        $self->error_message("Couldn't create annotator of class $annotator_version_subclass");
        die;
    }

    return $annotator;
}

sub _annotation_method_for_annotation_filter {
    my $self = shift;

    my %map = (
            gene => 'prioritized_transcripts',  # Top annotation per gene
            top  => 'prioritized_transcript',   # Top annotation between all genes
            none => 'transcripts',              # All transcripts, no filter
        );
    return $map{ $self->annotation_filter };
}

sub _make_reference_translator {
    my $self = shift;

    if($self->accept_reference_IUB_codes and $self->use_version >= 2) {
        #This will transform any IUB codes into a single base as detailed in the doc for the accept_reference_IUB_codes flag
        return sub {
            my $variant = $_[0];
            $variant->{reference} = Genome::Info::IUB->reference_iub_to_base($variant->{reference});
        };
    } else {
        return sub {};  # noop
    }
}

sub _main_annotation_loop {
    my($self, $annotator, $variant_svr) = @_;

    my $processed_variants = 0;
    my $chromosome_name = '';
    my ($annotation_start, $annotation_stop);

    my $annotation_method = $self->_annotation_method_for_annotation_filter;

    my $reference_translator = $self->_make_reference_translator();

    while ( my $variant = $variant_svr->next ) {

        # these are tracked by an END{} block for debugging
        $we_are_done_flag = 0;
        $last_variant_annotated = $variant;

        $variant->{type} = $self->infer_variant_type($variant);
        #make sure both the reference and the variant are in upper case
        $variant->{reference} = uc $variant->{reference};
        $variant->{variant} = uc $variant->{variant};

        # make a note when we begin and when we switch chromosomes
        unless ($variant->{chromosome_name} eq $chromosome_name) {
            if ($annotation_start) {
                $annotation_stop = Benchmark->new;
                my $annotation_time = timediff($annotation_stop, $annotation_start);
                $self->status_message("Annotating chromosome $chromosome_name took " . timestr($annotation_time)) if $self->benchmark;
            }

            $chromosome_name = $variant->{chromosome_name};

            $annotation_start = Benchmark->new;
            if ($self->benchmark) {
                $self->status_message("Annotation start for chromosome $chromosome_name");
            }
        }

        $reference_translator->($variant);

        if ($variant->{type} eq 'SNP') {
            # If we have an IUB code, annotate once per base... doesnt apply to things that arent snps
            my @variant_alleles = Genome::Info::IUB->variant_alleles_for_iub($variant->{reference}, $variant->{variant});

            for my $variant_allele (@variant_alleles) {
                # annotate variant with this allele
                $variant->{variant} = $variant_allele;
                $DB::single = 1;
                my @transcripts = $annotator->$annotation_method(%$variant);
                $self->_print_annotation($variant, \@transcripts);
            }
        } else {
            # get the data and output it
            my @transcripts = $annotator->$annotation_method(%$variant);
            $self->_print_annotation($variant, \@transcripts);
        }
        $processed_variants++;
        $self->status_message("$processed_variants variants processed " . scalar(localtime)) unless ($processed_variants % 10000);
    }

    $we_are_done_flag = 1;

    $annotation_stop = Benchmark->new;
    my $annotation_time = timediff($annotation_stop, $annotation_start);
    $self->status_message("Annotating $chromosome_name took " . timestr($annotation_time) . "\n") if $self->benchmark;

    return $processed_variants || '0 but true';
}

sub execute {
    my $self = shift;

    $self->_validate_parameters || return;

    $self->_print_starting_message();

    if ($self->_is_parallel) {
        $self->output_file($self->variant_file . ".out");
    }

    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    my $variant_svr = $self->_create_variant_reader() || return;

    $self->_setup_report_fh() || return;

    # annotate all of the input variants
    $self->status_message("Annotation start") if $self->benchmark;
    my $annotation_total_start = Benchmark->new;

    my $annotator = $self->_create_annotator();

    $self->status_message("Starting annotation loop at ".scalar(localtime));
    $self->status_message("  with annotator version ".$self->use_version);
    my $annotation_loop_start_time = time();

    Genome::DataSource::GMSchema->disconnect_default_handle if Genome::DataSource::GMSchema->has_default_handle;

    # This is where all the action happens!
    my $processed_variants = $self->_main_annotation_loop($annotator, $variant_svr);

    my $annotation_loop_stop_time = time();

    my $annotation_total_stop = Benchmark->new;
    my $total_time = timediff($annotation_total_stop, $annotation_total_start);
    $self->status_message('Total time to complete: ' . timestr($total_time, 'noc') . "\n\n") if $self->benchmark;

    my $timediff = $annotation_loop_stop_time - $annotation_loop_start_time;
    $timediff ||= 1;  # avoid division by zero below
    my $variants_per_sec = $processed_variants / $timediff;
    $self->status_message("Annotated $processed_variants variants in " . $timediff . " seconds.  "
                          . sprintf("%2.2f", $variants_per_sec) . " variants per second");

    $self->_close_report_fh();

    #clean up the temporary annotation data file 
    if ($self->variant_bed_file and $self->variant_file){
        unlink $self->variant_file || die("Could not remove converted variant file " . $self->variant_file);
    }

    return 1;
}

sub _transcript_report_fh {
    my ($self, $fh) = @_;
    $self->{_transcript_fh} = $fh if $fh;
    return $self->{_transcript_fh};
}

sub _temp_output_file {
    my ($self, $name) = @_;
    $self->{_temp_output_file} = $name if $name;
    return $self->{_temp_output_file};
}

sub _print_annotation {
    my ($self, $snp, $transcripts) = @_;

    # Basic SNP Info 
    my $snp_info_string = join
    (
        "\t", 
        map { $snp->{$_} } ($self->variant_attributes, $self->variant_output_attributes, $self->get_extra_columns),
    );

    # If we have no transcripts, print the original variant with dashes for annotation info
    unless( @$transcripts ) {
        $self->_transcript_report_fh->print
        (
            join
            (
                "\t",                   
                $snp_info_string,
                map({ '-' } $self->transcript_attributes),
            ), 
            "\n",
        );
        return 1;
    }

    # Otherwise, print an annotation line for each transcript we have
    for my $transcript ( @$transcripts )
    {
        $self->_transcript_report_fh->print
        (
            join
            (
                "\t",                   
                $snp_info_string,
                map({ $transcript->{$_} ? $transcript->{$_} : '-' } $self->transcript_attributes),
            ), 
            "\n",
        );
    }
    return 1;
}

sub get_extra_columns {
    my $self = shift;

    my $unparsed_columns = $self->extra_columns;
    return unless $unparsed_columns;

    my @columns = split(",", $unparsed_columns);
    chomp @columns;

    return @columns;
}

sub transcript_report_headers {
    my $self = shift;
    return ($self->variant_attributes, $self->variant_output_attributes, $self->get_extra_columns, $self->transcript_attributes);
}

sub _notify_cache_attempt{
    my $self = shift;
    
    my $message_content = <<END_CONTENT;
Hey Jim,

This is a cache attempt on %s running as PID $$ and user %s.

I told the user I'm not freaking doing it, and am working normally.  Just wanted to give you a heads up.

Your pal,
Genome::Model::Tools::Annotate::TranscriptVariants

END_CONTENT
    require MIME::Lite;
    my $msg = MIME::Lite->new(From    => sprintf('"Genome::Utility::Filesystem" <%s@genome.wustl.edu>', $ENV{'USER'}),
            To      => 'jweible@genome.wustl.edu',
            Subject => 'Attempt to cache annotation data directory',
            Data    => sprintf($message_content, hostname, $ENV{'USER'}),
            );
    $msg->send();
}
1;

=pod

=head1 Name

Genome::Model::Tools::Annotate::TranscriptVariants

=head1 Synopsis

Goes through each variant in a file, retrieving annotation information from Genome::Transcript::VariantAnnotator.

=head1 Usage

 in the shell:

     gmt annotate transcript-variants --variant-file myinput.csv --output-file myoutput.csv

 in Perl:

     $success = Genome::Model::Tools::Annotate::TranscriptVariants->execute(
         variant_file => 'myoutput.csv',
         output_file => 'myoutput.csv',
     );

=head1 Methods

=over

=item variant_file

An input list of variants.  The format is:
 chromosome_name
 start
 stop 
 reference
 variant

The mutation type will be inferred based upon start, stop, reference, and variant alleles.

 Any number of additional columns may be in the input, but they will be disregarded.

=item output_file

The list of transcript changes which would occur as a result of the associated genome sequence changes.

One variant may result in multiple transcript entries if it intersects multiple transcripts.  One 
transcript may occur multiple times in results if multiple variants intersect it.

=item 

=back

=head1 See Also

B<Genome::Transcript::VariantAnnotator>, 

=head1 Disclaimer

Copyright (C) 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

Core Logic:

 B<Xiaoqi Shi> I<xshi@genome.wustl.edu>

Optimization, Testing, Data Management:

 B<Dave Larson> I<dlarson@genome.wustl.edu>
 B<Eddie Belter> I<ebelter@watson.wustl.edu>
 B<Gabriel Sanderson> I<gsanderes@genome.wustl.edu>
 B<Adam Dukes> I<adukes@genome.wustl.edu>
 B<Anthony Brummett> I<abrummet@genome.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Genome/Model/Tools/Annotate/TranscriptVariants.pm $
#$Id: TranscriptVariants.pm 44679 2009-03-16 17:55:52Z adukes $
