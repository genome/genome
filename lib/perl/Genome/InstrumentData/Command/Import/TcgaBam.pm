package Genome::InstrumentData::Command::Import::TcgaBam;

use strict;
use warnings;

use Genome;

use XML::Simple;

use Genome::Sample::Command::Import;
require File::Basename;

class Genome::InstrumentData::Command::Import::TcgaBam {
    is  => 'Command::V2',
    has => [
        original_data_path => {
            is => 'Text',
            doc => 'Original data path of import data file',
        },
        tcga_name => {
            is => 'Text',
            doc => 'TCGA name for imported file (must be found in metadata.xml or specified)',
            is_optional => 1,
        },
        target_region => {
            is => 'Text',
            doc => 'Target region set name (capture) or \'none\' (must be found in metadata.xml or specified)',
            is_optional => 1,
        },
        remove_original_bam => {
            is => 'Boolean',
            doc => 'Remove (DELETE!) the original bam file after importation, without warning',
            default => 0,
            is_optional => 1,
        },
        skip_refalign => {
            is => 'Boolean',
            doc => 'Skip creation of a reference alignment model (e.g. for RNA-seq)',
            default => 0,
            is_optional => 1,
        },
        no_md5 => {
            is => 'Boolean',
            default => 0,
            is_optional => 1,
            doc => 'Completely skip the md5 bam file validation process'
        },
        bam_md5 => {
            is => 'Text',
            is_optional => 1,
            doc => "The md5 hash of the bam file (if not specified it will be looked for in the metadata file or a .md5 file in the bam directory)",
        },
        import_source_name => {
            is => 'Text',
            doc => 'The source name for the imported file, i.e. \'broad\' (must be found in metadata.xml or specified)',
            is_optional => 1,
        },
        description  => {
            is => 'Text',
            doc => 'General description of import data, like which software maq/bwa/bowtie to use to generate this data',
            is_optional => 1,
        },
        read_count  => {
            is => 'Number',
            doc => 'Total read count of import data',
            is_optional => 1,
        },
        base_count  => {
            is => 'Number',
            doc => 'Total base count of import data',
            is_optional => 1,
        },
        metadata_file => {
            is => "Text",
            is_optional => 1,
            doc => "The fullpath of a .xml metadata file which contains (among other things) the md5 of the bam file.  If not specified, a metadata file is looked for in the same directory as the bam file (named metadata.xml).",
        },
        analysis_id => {
            is => 'Text',
            is_optional => 1,
            doc => 'From CGHub (may be found in metadata.xml or specified)',
        },
        aliquot_id => {
            is => 'Text',
            is_optional => 1,
            doc => 'From CGHub (may be found in metadata.xml or specified)',
        },
        participant_id => {
            is => 'Text',
            is_optional => 1,
            doc => 'From CGHub (may be found in metadata.xml or specified)',
        },
        sample_id => {
            is => 'Text',
            is_optional => 1,
            doc => 'From CGHub (may be found in metadata.xml or specified)',
        },
        reference_sequence_build_id => {
            is => 'Text',
            valid_values => [qw/ 101947881 106942997 /],
            doc => 'The id of the reference sequence that the data was aligned against. Versions: 101947881 (NCBI-human-build36), 106942997 (GRCh37-lite-build37).',
        },
        reference_sequence_build => {
            is_optional => 1,
            calculate_from => [qw/ reference_sequence_build_id /],
            calculate => q| return Genome::Model::Build::ImportedReferenceSequence->get($reference_sequence_build_id); |,
        },
        _model => { is_optional => 1, },
        _inst_data => { is_optional => 1, },
        import_instrument_data_id => { via => '_inst_data', to => 'id', is_optional => 1, },
        _allocation => { via => '_inst_data', to => 'allocations', is_optional => 1, },
        _absolute_path => { via => '_allocation', to => 'absolute_path', is_optional => 1, },
        _new_bam => {
            is_optional => 1,
            calculate_from => [qw/ _absolute_path /],
            calculate => q| $_absolute_path.'/all_sequences.bam' |,
        },
        _new_md5 => {
            is_optional => 1,
            calculate_from => [qw/ _new_bam /],
            calculate => q| $_new_bam.'.md5' |,
        },
    ],
    doc => 'Create an instrument-data AND and alignment model for a BAM',
};

sub help_detail {
    return <<HELP;
    This command imports a BAM for a TCGA patient. Workflow:
    * creates an instrument data
    * copies the BAM into the allocated spaace of the instruemt data
    * validates the MD5 of the BAM (optional)
    * creates a model and requests a build (optional)
HELP
}

sub execute {
    my $self = shift;

    # resolve command-line arguments and read in from metadata (if avail)
    $self->_resolve_args;

    # Ref Seq
    if ( not $self->reference_sequence_build ) {
        $self->error_message('No reference sequence build given.');
        return;
    }

    # Validate BAM
    my $bam_ok = $self->_validate_bam;
    return if not $bam_ok;

    # Create inst data
    my $inst_data = $self->_create_imported_instrument_data;
    return if not $inst_data;

    # create attributes out of some of the metadata
    $self->_create_attributes($inst_data);

    # Copy and create md5 at the same time w/ tee
    my $copy = $self->_copy_and_generate_md5;
    if ( not $copy ) {
        $self->_bail;
        return;
    }

    # Validate copied BAM
    my $validate = $self->_validate_copied_bam;
    if ( not $validate ) {
        $self->_bail;
        return;
    }

    # Add stats to the instrument-data taken from flagstat, etc
    unless($self->_add_stats){
        die $self->error_message("Could not complete flagstat operation on imported bam");
    }

    # Rm Original BAM
    if($self->remove_original_bam){
        $self->_remove_original_bam; # no error check
    }

    # Unless otherwise specified, realign the imported BAM with our standard ref-alignment pipe
    if($self->skip_refalign) {
        $self->status_message("Skipping creation of a ref-alignment model.");
    }
    else {
        $self->_create_model_and_request_build; # no error check, prints messages
    }

    $self->status_message("Importation of BAM completed successfully.");
    $self->status_message("Your instrument-data id is ".$self->import_instrument_data_id);

    return 1;
}

sub _resolve_args {
    my ($self) = @_;

    my %metadata = %{$self->_read_in_metadata};

    my @optional_arg_names = qw| aliquot_id analysis_id participant_id sample_id |;
    for my $arg_name (@optional_arg_names){
        $self->$arg_name($self->_resolve_single_arg($arg_name, $metadata{$arg_name}));
    }


    my @required_arg_names = qw| import_source_name tcga_name target_region |;
    for my $arg_name (@required_arg_names){
        $self->$arg_name($self->_resolve_single_arg($arg_name, $metadata{$arg_name}));
        if(not defined $self->$arg_name) {
            die $self->error_message("Required argument ($arg_name) was not passed and couldn't be found in the metadata.");
        }
    }

    # handle bam_md5 carefully
    unless($self->no_md5) {
        $self->bam_md5($self->_resolve_single_arg('bam_md5', $metadata{'bam_md5'}));
        $self->bam_md5($self->_resolve_bam_md5);
        if(not defined $self->bam_md5) {
            die $self->error_message("Required argument (bam_md5) was not passed and couldn't be found in the metadata or an .md5 file in the directory where the bam file is.");
        }
    }
    return 1;
}

sub _resolve_bam_md5{
    my ($self) = @_;

    my $md5;
    if(not $self->bam_md5){
        my $bam = $self->original_data_path;
        my $md5_file = $bam.".md5";
        if(not -s $md5_file) {
            $self->error_message("Did not find md5 file ($md5_file) for bam ($bam)");
            return;
        }
        $md5 = $self->_get_md5_from_file($md5_file);
        return if not $md5;
        $self->bam_md5($md5);
    } else {
        $md5 = $self->bam_md5;
    }
    return $md5
}

# resolve the value of a single argument, falls back to metadata value.
sub _resolve_single_arg {
    my ($self, $arg_name, $meta_value) = @_;

    if(defined $self->$arg_name){
        $self->status_message(sprintf('Argument (%s) was passed in as (%s)', $arg_name, $self->$arg_name));
        return $self->$arg_name;
    }
    if(defined $meta_value){
        $self->status_message(sprintf('Argument (%s) was found in the metadata file in as (%s)', $arg_name, $meta_value));
    }
    return $meta_value;
}

sub _read_in_metadata {
    my ($self) = @_;

    # check for xml or look for one in directory.
    my $no_xml = 1;
    if(not defined $self->metadata_file){
        (undef, my $base_path) = File::Basename::fileparse($self->original_data_path);
        $self->metadata_file($base_path . '/metadata.xml');
        if(-e $self->metadata_file){
            $no_xml = 0;
        }
    }
    my %metadata;
    if(not $no_xml) {
        return $self->_parse_metadata_file($self->metadata_file);
    }
    return {};
}

sub _create_attributes {
    my ($self, $inst_data) = @_;

    my @to_become_attributes = qw| aliquot_id analysis_id participant_id sample_id |;
    for my $attr_name (@to_become_attributes) {
        if( defined $self->$attr_name) {
            my $attribute = $inst_data->add_attribute(
                    attribute_label => $attr_name,
                    attribute_value => $self->$attr_name,
                    nomenclature => 'CGHub',
            );
            if(not $attribute) {
                $self->_bail;
                die $self->error_message(sprintf("Couldn't add attribute '%s' to instrument-data(%s).", $attr_name, $inst_data->id))
            }
        }
    }
}

sub _parse_metadata_file {
    my ($self, $metadata_file) = @_;

    my $md = XMLin($metadata_file);

    my %metadata;
    $metadata{bam_md5}            = $self->_md5_checksum_content($md);
    $metadata{tcga_name}          = $md->{Result}->{legacy_sample_id};
    $metadata{import_source_name} = $md->{Result}->{center_name};
    $metadata{analysis_id}        = $md->{Result}->{analysis_id};
    $metadata{aliquot_id}         = $md->{Result}->{aliquot_id};
    $metadata{participant_id}     = $md->{Result}->{participant_id};
    $metadata{sample_id}          = $md->{Result}->{sample_id};

    my $library_strategy          = $md->{Result}->{library_strategy};
    my $target_region;
    if($library_strategy eq 'WGS') { #whole genome
        $target_region = 'none';
    } elsif($library_strategy eq 'WXS') { #exome
        $target_region = $self->reference_sequence_build_id eq '106942997'?
            'agilent_sureselect_exome_version_2_broad_refseq_cds_only_hs37':
            'agilent sureselect exome version 2 broad refseq cds only';
    } elsif($library_strategy eq 'RNA-Seq') { #RNA
        $target_region = 'none';
    } else {
        $self->warning_message('Unknown library strategy: ' . $library_strategy);
    }
    $metadata{target_region}      = $target_region;

    return \%metadata
}

sub _md5_checksum_content {
    my( $self, $md ) = @_;

    my $content = eval{ $md->{Result}->{files}->{file}->{checksum}->{content}; };
    if( not $content ) {
        my $filename = File::Basename::basename( $self->original_data_path );
        my ($fileinfo) = eval{ grep {$_->{filename} eq $filename} @{$md->{Result}->{files}->{file}} };
        if( not $fileinfo ) {
            $self->error_message("Failed to get files md5 info from metadata file");
            return;
        }
        $content = $fileinfo->{checksum}->{content};
    }
    if( not $content ) {
        $self->error_message("Failed to get checksum content for file");
        return;
    }

    return $content;
}

sub _validate_bam {
    my $self = shift;

    return 1;

    my $bam = $self->original_data_path;
    if ( not -s $bam ) {
        $self->error_message('BAM does not exist: '.$bam);
        return;
    }

    if ( $bam !~ /\.bam$/ ) { # why?
        $self->error_message('BAM does not have extension ".bam": '.$bam);
        return;
    }

    return 1;
}

sub _validate_md5 {
    my $self = shift;

    $self->status_message('Validate BAM MD5...');

    my $bam = $self->original_data_path;
    my $md5_file = $bam . ".md5";
    if ( not -s $md5_file ) {
        $self->error_message("Did not find md5 file ($md5_file) for bam ($bam)");
        return;
    }

    $self->status_message("Getting BAM MD5 from file: $md5_file");
    my $md5_fh = eval{ Genome::Sys->open_file_for_reading($md5_file); };
    if ( not $md5_fh ) {
        $self->error_message("Cannot open BAM MD5 file ($md5_file): $@");
        return;
    }
    my $line = $md5_fh->getline;
    chomp $line;
    my ($md5) = split(" ", $line);
    if ( not $md5 ) {
        $self->error_message('No MD5 in file: '.$md5_file);
        return;
    }
    $self->status_message("Got BAM MD5 from file: $md5");

    $self->status_message("Calculate MD5 for bam: $bam");
    $self->status_message("This may take a bit...");
    my $calculated_md5 = $self->_md5_for_file($bam);
    if ( not $calculated_md5 ) {
        $self->error_message("Failed to calculate md5 for BAM: $bam.");
        return;
    }
    $self->status_message("Calculated MD5: $calculated_md5");

    $self->status_message('Validate MD5...');
    if ( $md5 ne $calculated_md5 ) {
        $self->error_message("Calculated BAM MD5 ($calculated_md5) does not match MD5 from file ($md5)");
        return;
    }
    $self->status_message('Validate MD5...OK');

    return 1;
}

sub _md5_for_file {
    my ($self, $file) = @_;

    if ( not -e $file ) {
        $self->error_message("Cannot get md5 for non existing file: $file");
        return;
    }

    my ($md5) = split(/\s+/, `md5sum $file`);

    if ( not defined $md5 ) {
        $self->error_message("No md5 returned for file: $file");
        return;
    }

    return $md5;
}

sub _create_imported_instrument_data {
    my $self = shift;

    $self->status_message('Create imported instrument data...');

    my $tcga_name = $self->tcga_name;

    # Get or create library
    my $sample_importer = Genome::Sample::Command::Import::Tcga->create(
        name => $tcga_name,
    );
    if ( not $sample_importer ) {
        $self->error_message('Could not create TCGA sample importer to get or create library');
        return;
    }
    $sample_importer->dump_status_messages(1);
    if ( not $sample_importer->execute ) {
        $self->error_message('Could not execute TCGA sample importer to get or create library');
        return;
    }
    my $library = $sample_importer->_library;

    my $description = $self->description || "imported ".$self->import_source_name." bam, tcga name is ".$tcga_name;
    if($self->no_md5){
        $description = $description . ", no md5 file was provided with the import.";
    }

    my %params = (
        original_data_path => $self->original_data_path,
        sequencing_platform => "solexa",
        import_format => "bam",
        reference_sequence_build => $self->reference_sequence_build,
        library => $library,
        description => $description,
    );

    my $target_region;
    unless ($self->target_region eq 'none') {
        $target_region = $self->target_region;
        my @feature_lists = Genome::FeatureList->get(name => $target_region);
        if (not @feature_lists or @feature_lists > 1) {
            $self->error_message("Invalid target region: " . $target_region);
            return;
        }
    }

    if ( $target_region ) {
        $params{target_region_set_name} = $target_region
    }

    my $import_instrument_data = Genome::InstrumentData::Imported->create(%params);

    unless ($import_instrument_data) {
       $self->error_message('Failed to create imported instrument data for '.$self->original_data_path);
       return;
    }
    $self->_inst_data($import_instrument_data);

    my $instrument_data_id = $import_instrument_data->id;
    $self->status_message("Instrument data: $instrument_data_id is imported");

    my $kb_usage = $import_instrument_data->calculate_alignment_estimated_kb_usage;
    unless ($kb_usage) {
        $self->error_message('Cannot calculate kb usage for BAM: '.$self->original_data_path);
        $import_instrument_data->delete;
        return 1;
    }

    my $alloc_path = sprintf('alignment_data/imported/%s', $instrument_data_id);

    my %alloc_params = (
        disk_group_name     => 'info_alignments',
        allocation_path     => $alloc_path,
        kilobytes_requested => $kb_usage,
        owner_class_name    => $import_instrument_data->class,
        owner_id            => $import_instrument_data->id,
    );

    my $disk_alloc = Genome::Disk::Allocation->allocate(%alloc_params);
    unless ($disk_alloc) {
        $self->error_message("Failed to get disk allocation with params:\n". Data::Dumper::Dumper(%alloc_params));
        $import_instrument_data->delete;
        return 1;
    }
    $self->status_message("Alignment allocation created for $instrument_data_id .");

    $self->status_message('Create imported instrument data...OK');

    return $self->_inst_data;
}

sub _copy_and_generate_md5 {
    my $self = shift;

    $self->status_message("Copy BAM and generate MD5");

    my $bam = $self->original_data_path;
    my $new_bam = $self->_new_bam;
    my $new_md5 = $self->_new_md5;
    my $cmd = "tee $new_bam < $bam | md5sum > $new_md5";
    $self->status_message('Cmd: '.$cmd);
    my $eval = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $eval ) {
        $self->status_message('Copy BAM and generate MD5 FAILED: '.$@);
        return;
    }

    $self->status_message("Copy BAM and generate MD5");

    return 1;
}


sub _validate_copied_bam {
    my $self = shift;

    $self->status_message('Validate copied BAM...');

    $self->status_message('Validate size...');
    my $bam = $self->original_data_path;
    my $bam_size = -s $bam;
    my $new_bam = $self->_new_bam;
    my $new_bam_size = -s $new_bam;
    if ( $bam_size != $new_bam_size ) {
        $self->error_message("Copied BAM ($new_bam) size ($new_bam_size) does not match original BAM ($bam) size ($bam_size)");
        return;
    }
    $self->status_message("Validate size OK: $bam_size v. $new_bam_size");

    if ( $self->no_md5 ) {
        $self->status_message('Validate copied BAM...OK');
        return 1;
    }

    my $md5 = $self->bam_md5;

    my $new_md5_file = $self->_new_md5;
    if ( not -s $new_md5_file ) {
        $self->error_message("Did not find md5 file ($new_md5_file) for copied bam ($new_bam)");
        return;
    }
    my $new_md5 = $self->_get_md5_from_file($new_md5_file);
    return if not $new_md5;

    if ( $md5 ne $new_md5 ) {
        $self->error_message("Copied BAM MD5 ($new_md5) does not match MD5 from file ($md5)");
        return;
    }
    $self->status_message('Validate MD5...OK');

    $self->status_message('Validate copied BAM...OK');

    return 1;
}



sub _get_md5_from_file {
    my ($self, $file) = @_;

    $self->status_message("Get MD5 from file: $file");

    my $fh = eval{ Genome::Sys->open_file_for_reading($file); };
    if ( not $fh ) {
        $self->error_message("Cannot open file ($file): $@");
        return;
    }
    my $line = $fh->getline;
    chomp $line;
    my ($md5) = split(" ", $line);
    if ( not $md5 ) {
        $self->error_message('No MD5 in file: '.$file);
        return;
    }

    $self->status_message("MD5 for file ($file): $md5");

    return $md5;
}

sub XX_rsync_bam {
    my $self = shift;

    my $bam = $self->original_data_path;
    my $new_bam = $self->_new_bam;
    $self->status_message("Rsync BAM from $bam to $new_bam");

    my $cmd = "rsync -acv $bam $new_bam";
    $self->status_message('Rsync cmd: '.$cmd);
    my $rsync = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rsync ) {
        $self->error_message('Rsync cmd failed: '.$@);
        $self->status_message('Removing disk allocation: '.$self->_allocation->id);
        unlink $new_bam if -e $new_bam;
        $self->_allocation->deallocate;
        $self->error_message('Removing instrument data: '.$self->_inst_data->id);
        $self->_inst_data->delete;
        return;
    }
    $self->status_message('Rync BAM...OK');

    return 1;
}

sub _bail {
    my $self = shift;

    $self->status_message('Copy BAM and generate MD5 FAILED: '.$@);
    $self->status_message('Removing disk allocation: '.$self->_allocation->id);
    my $new_bam = $self->_new_bam;
    unlink $new_bam if -e $new_bam;
    my $new_md5 = $self->_new_md5;
    unlink $new_md5 if -e $new_md5;
    $self->_allocation->deallocate;
    $self->status_message('Removing instrument data: '.$self->_inst_data->id);
    $self->_inst_data->delete;

    return 1;
}

sub _add_stats {
    my $self = shift;
    my $data = $self->_run_flagstat;
    my $inst_data = $self->_inst_data;
    $inst_data->read_count($data->{'total_reads'});
    $inst_data->fragment_count($data->{'total_reads'}*2);
    $inst_data->read_length($self->_read_length);
    $inst_data->base_count(int($inst_data->read_length) * int($inst_data->fragment_count));
    if($data->{'reads_paired_in_sequencing'} > 0){
        $inst_data->is_paired_end(1);
    }
    else {
        $inst_data->is_paired_end(0);
    }
    $self->_inst_data($inst_data);
    return 1;
}

sub _read_length {
    my $self = shift;
    $self->status_message("Now calculating read_length via gmt sam->read_length");
    my $sam = Genome::Model::Tools::Sam->create;
    my $read_length = $sam->read_length($self->_inst_data->bam_path);
    unless(defined($read_length)){
        die $self->error_message("Was not able to run gmt sam->read_length");
    }
    $self->status_message("Finished calculating read_length. read_length=".$read_length);
    return $read_length;
}

sub _run_flagstat {
    my $self = shift;
    unless(defined($self->_inst_data)){
        die $self->error_message("No instrument data found in self->_inst_data");
    }
    my $flagstat_file = $self->_inst_data->bam_path . ".flagstat";
    my $flagstat_object = Genome::Model::Tools::Sam::Flagstat->create( bam_file => $self->_inst_data->bam_path, output_file => $flagstat_file);
    unless(-s $flagstat_file){
        $self->status_message("Generating flagstat file now...");
        unless($flagstat_object->execute){
            die $self->error_message("Failed to run gmt sam flagstat");
        }
        $self->status_message("Flagstat file created");
    }
    my $flag_data = $flagstat_object->parse_file_into_hashref($flagstat_file);

    return $flag_data;
}

sub _remove_original_bam {
    my $self = shift;

    $self->status_message("Now removing original bam in 10 seconds.");
    for (1..10){
        sleep 1;
        print "slept for ".$_." seconds.\n";
    }
    my $bam_path = $self->original_data_path;
    unless(-s $bam_path){
        $self->error_message("Could not locate file to remove at ".$bam_path."\n");
        die $self->error_message;
    }
    unlink($bam_path);
    if(-s $bam_path){
        $self->error_message("Could not remove file at ".$bam_path."\n");
        $self->error_message("Check file permissions.");
    }else{
        $self->status_message("Original bam file has been removed from ".$bam_path);
    }

    return 1;
}

sub _create_model_and_request_build {
    my $self  = shift;
    my $pp_id = Genome::ProcessingProfile::ReferenceAlignment->default_profile_id;

    $self->status_message('Create model and request build');

    my $pp = Genome::ProcessingProfile::ReferenceAlignment->get($pp_id);
    if ( not $pp ) {
        $self->error_message("Cannot find ref align processing profile for $pp_id to create model");
        return;
    }
    $self->status_message('Processing profile: '.$pp->name);

    my $sample = $self->_inst_data->sample;
    $self->status_message('Sample: '.$sample->name);

    my $refseq = $self->reference_sequence_build;
    $self->status_message('Reference build: '.$refseq->__display_name__);

    my $base_name = $sample->name.'.'.$refseq->version.'.refalign';
    my $name = $base_name;
    my $i = 0;
    while ( my $model = Genome::Model::ReferenceAlignment->get(name => $name) ) {
        $name = $base_name.'-'.++$i;
    }

    my $model = Genome::Model::ReferenceAlignment->create(
        name => $name,
        processing_profile => $pp,
        reference_sequence_build => $refseq,
        subject_id => $sample->id,
        subject_class_name => $sample->class,
        build_requested => 1,
        auto_assign_inst_data => 0,
    );
    if ( not $model ) {
        $self->error_message('Failed to create model');
        return;
    }
    $self->_model($model);
    $self->status_message('Model id: '.$model->id);
    $self->status_message('Model name: '.$model->name);

    my $dbsnp = Genome::Model::ImportedVariationList->dbsnp_build_for_reference($refseq);
    if ( $dbsnp ) {
        $self->status_message('dbSNP build: '.$dbsnp->__display_name__);
        $model->dbsnp_build($dbsnp);
    }

    my $annotation = Genome::Model::ImportedAnnotation->annotation_build_for_reference($refseq);
    if ( $annotation ) {
        $self->status_message('Annotation build: '.$annotation->__display_name__);
        $model->annotation_reference_build($annotation);
    }

    my $add = $model->add_instrument_data( $self->_inst_data );
    if ( not $add ) {
        $self->error_message('Failed to add instrument data to model');
        $model->delete;
        return;
    }

    $self->status_message('Create model and request build...OK');

    return $model;
}

1;

