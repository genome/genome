package Genome::Model::Build::ReferenceSequence;
use strict;
use warnings;
use Genome;
use Genome::Sys::LockProxy qw();
use File::Spec;
use File::Path;
use File::Copy;

require Carp;
use Regexp::Common;
use POSIX;
use Set::Scalar;
use Filesys::Df qw();

class Genome::Model::Build::ReferenceSequence {
    is => 'Genome::Model::Build',
    has => [
        name => {
            via => '__self__',
            to => 'build_name',
        },
        calculated_name => {
            calculate_from => ['model_name','version'],
            calculate => q{
                my $name = "$model_name-build";
                $name .= $version if defined $version;
                $name =~ s/\s/-/g;
                return $name;
            },
        },

        species_name => { via => 'subject', to => 'name' },

        manifest_file_path => {
            is => 'Text',
            calculate_from => ['data_directory'],
            calculate => q(
                if($data_directory){
                    return join('/', $data_directory, 'manifest.tsv');
                }
            ),
        },
        _sequence_filehandles => {
            is => 'Hash',
            is_optional => 1,
            is_transient => 1,
            doc => 'file handle per chromosome for reading sequences so that it does not need to be constantly closed/opened',
        },
        _local_cache_dir_is_verified => { is => 'Boolean', default_value => 0, is_optional => 1, is_transient => 1,},
        fasta_md5 => {
            is => 'Text',
            is_optional => 1,
            is_input => 1,
            doc => 'MD5 of the reference FASTA.  Allows for shortcutting of rederivable references.',
        },

    ],

    has_optional_input => [
        # In order to set this property use Genome::Model::ReferenceSequence::Command::CreateFeatureListInput.
        # In order to obtain the bed file used in feature list creation for segmental duplications go to the UCSC table browser and download the segmental duplications track for the appropriate reference sequence
        # current URL: http://genome.ucsc.edu/cgi-bin/hgTables - track: Segmental Dups - table: genomicSuperDups - output format: BED
        segmental_duplications => {
            is => 'Genome::FeatureList',
        },
        # In order to set this property use Genome::Model::ReferenceSequence::Command::CreateFeatureListInput.
        # In order to obtain the bed file used in feature list creation for segmental duplications go to the UCSC table browser and download the self chain track for the appropriate reference sequence
        # current URL: http://genome.ucsc.edu/cgi-bin/hgTables - track: Self Chain - table: chainSelf - output format: BED
        self_chain => {
            is => 'Genome::FeatureList',
        },
        build_name => {
            is => 'Text',
        },

        allosome_names => {
            is => 'Text',
        },

        derived_from => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'Identifies the parent build from which this one is derived, if any.',
        },
        coordinates_from => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'Used to indicate that this build is on the same coordinate system as another.',
        },
        append_to => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'If specified, the created reference will be logically appended to the one specified by this parameter for aligners that support it.',
        },
        combines => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'If specified, merges several other references into one.',
            is_many => 1,
        },
    ],
    has_many_optional => [
        convertible_to_bridges => {
            is => 'Genome::Model::Build::ReferenceSequence::Converter',
            reverse_as => 'source_reference_build',
            doc => 'converters which convert this reference to other references',
        },
        convertible_to => {
            is => 'Genome::Model::Build::ReferenceSequence',
            via => 'convertible_to_bridges',
            to => 'destination_reference_build',
            doc => 'other references to which this reference can be converted',
        },
        convertible_from_bridges => {
            is => 'Genome::Model::Build::ReferenceSequence::Converter',
            reverse_as => 'destination_reference_build',
            doc => 'converters which convert other references to this reference',
        },
        convertible_from => {
            is => 'Genome::Model::Build::ReferenceSequence',
            via => 'converters_from_bridges',
            to => 'source_reference_build',
            doc => 'other references that can be converted to this reference',
        },
    ],
    has_transient_optional => [
        _seqdict => {
            is => 'HASH',
            doc => 'A hash ref version of the sequence dictionary', 
        },
    ],

    doc => 'a specific version of a reference sequence',
};


sub create {
    my $self = shift;
    my $build = $self->SUPER::create(@_);

    # Let's store the name as an input instead of relying on calculated properties
    $build->name($build->calculated_name);

    if ($build->generate_sequence_uri) {
        $build->sequence_uri($build->external_url);
    }

    if ($build->can('derived_from') && $build->derived_from) {
        $build->coordinates_from($build->derived_from_root);
    }

    # set this for the assembly name as well if there is not one already.
    if (!$build->assembly_name) {
        $build->assembly_name($build->calculated_name);
    }

    $self->status_message("Created reference sequence build with assembly name " . $build->name);

    return $build;
}

sub validate_for_start_methods {
    my $self = shift;
    my @methods = $self->SUPER::validate_for_start_methods;
    push @methods, 'check_derived_from_links';
    return @methods;
}

sub check_derived_from_links {
    my $self = shift;
    my @tags;

    # this will die on circular links
    eval { my $coords = $self->derived_from_root(); };
    if ($@) {
        push @tags, UR::Object::Tag->create(
            type => 'error',
            properties => ['derived_from'],
            desc => $@);
    }

    if (defined $self->derived_from and $self->derived_from->id eq $self->id) {
        push @tags, UR::Object::Tag->create(
            type => 'error',
            properties => ['derived_from'],
            desc => "A build cannot be explicitly derived from itself!");
    }

    return @tags;
}

sub get {
    my $self = shift;
    my @results = $self->SUPER::get(@_);
    return $self->SUPER::get(@_) if @results;

    my @caller = caller(1);
    if($caller[3] && $caller[3] =~ m/Genome::Model::Build::ImportedReferenceSequence::get/){
        return;
    }else{
        return Genome::Model::Build::ImportedReferenceSequence->get(@_);
    }
}

#TODO: we need to figure out the ucsc api to download new ones when needed
#For now, this is hard-coded
sub get_or_create_ucsc_tiering_directory {
    my $self = shift;

    my $dir = File::Spec->join($self->data_directory, 'ucsc_tiering_files');

    if (-d $dir) {
        return $dir;
    } else {
        $self->status_message('UCSC Tiering Directory not currently available for this species: '. $self->species_name);
        return;
    }
}

sub is_derived_from {
    my ($self, $build, $seen) = @_;
    $seen = {} if !defined $seen;
    if (exists $seen->{$self->id}) {
        die "Circular link found in derived_from chain. Current build: " . $self->__display_name__ . ", derived from: " .
            $self->derived_from->__display_name__ . ", seen: " . join(',', keys %{$seen});
    }

    return 1 if $build->id eq $self->id;
    return 0 if !defined $self->derived_from;

    # recurse
    $seen->{$self->id} = 1;
    return $self->derived_from->is_derived_from($build, $seen);
}

sub derived_from_root {
    my ($self) = @_;
    my $from = $self;
    my %seen = ($self->id => 1);
    while (defined $from->derived_from) {
        $from = $from->derived_from;
        if (exists $seen{$from->id}) {
            die "Circular link found in derived_from chain while calculating 'derived_from_root'.".
                " Current build: " . $self->__display_name__ . ", derived from: " .
                $from->derived_from->__display_name__ . ", seen: " . join(',', keys %seen);
        }
        $seen{$from->id} = 1;
    }
    return $from;
}

# check compatibility with another reference sequence build
sub is_compatible_with {
    my ($self, $rsb) = @_;
    return if !defined $rsb;
    my $coords_from = $self->coordinates_from; # $self;
    my $other_coords_from = $rsb->coordinates_from; # $rsb;

    return 1 if $self->id eq $rsb->id;

    if($coords_from and $other_coords_from) {
        return 1 if $coords_from->id eq $other_coords_from->id;
    }

    if($coords_from) {
        return 1 if $coords_from->id eq $rsb->id;
        return 1 if $coords_from->is_compatible_with($rsb);
    }

    if($other_coords_from) {
        return 1 if $self->id eq $other_coords_from->id;
        return 1 if $self->is_compatible_with($other_coords_from);
    }

    return;
}

sub __display_name__ {
    my $self = shift;
    my $txt = $self->name . " (" . $self->id . ")";
    return $txt;
}

sub calculate_estimated_kb_usage {
    my $self = shift;
    for my $i ($self->inputs) {
        my $k = $i->name;
        my $v = $i->value_id;
        $self->status_message("INPUT: $k=$v\n");
    }

    my $fastaSize = -s $self->fasta_file;
    if(defined($fastaSize) && $fastaSize > 0)
    {
        $fastaSize = POSIX::ceil($fastaSize * 3 / 1024);
    }
    else
    {
        $fastaSize = $self->SUPER::calculate_estimated_kb_usage();
    }
    return $fastaSize;
}

sub sequence {
    my ($self, $chromosome, $start, $stop) = @_;

    my $f = $self->get_or_create_sequence_filehandle($chromosome);
    return unless ($f);

    my $seq;
    $f->seek($start - 1, 0);
    $f->read($seq, $stop - $start + 1);

    return $seq;
}

sub get_or_create_sequence_filehandle {
    my ($self, $chromosome) = @_;
    my $filehandles = $self->_sequence_filehandles;

    my $basesFileName = $self->get_bases_file($chromosome);
    return $filehandles->{$chromosome} if ($filehandles->{$chromosome});

    my $fh = IO::File->new($basesFileName);
    unless ($fh) {
        $self->error_message("Failed to open bases file \"$basesFileName\".");
        return;
    }

    $filehandles->{$chromosome} = $fh;
    $self->_sequence_filehandles($filehandles);
    return $fh;
}

sub get_bases_file {
    my $self = shift;
    my ($chromosome) = @_;

    my $bases_dir = join('/', $self->data_directory, 'bases');
    unless(-d $bases_dir) {
        #for backwards-compatibility--old builds stored the .bases files directly in the data directory
        #TODO remove this conditional once snapshots prior to this change are in use and the files have been moved
        #in all older I.R.S. builds
        $bases_dir = $self->data_directory;
    }

    # grab the dir here?
    my $bases_file = $bases_dir . "/" . $chromosome . ".bases";

    return $bases_file;
}

sub primary_consensus_path {
    my ($self, $format) = @_;

    return $self->full_consensus_path($format) unless $self->append_to;

    $format ||= 'bfa';
    my $file = $self->data_directory . '/appended_sequences.'. $format;
    return $file;
}

sub full_consensus_path {
    my ($self, $format) = @_;
    $format ||= 'bfa';

    my $file = $self->data_directory . '/all_sequences.'. $format;
    unless (-e $file){
        $file = $self->data_directory . '/ALL.'. $format;
        unless (-e $file){
            $self->error_message("Failed to find " . $self->data_directory . "/all_sequences.$format");
            return;
        }
    }
    return $file;
}

#This is for samtools faidx output that can be used as ref_list for
#SamToBam conversion
sub full_consensus_sam_index_path {
    my $self        = shift;
    my $sam_version = shift;

    my $data_dir = $self->data_directory;
    my $fa_file  = $self->full_consensus_path('fa');
    my $idx_file = $fa_file.'.fai';

    unless (-e $idx_file) {
        my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($sam_version);
        my $cmd      = $sam_path.' faidx '.$fa_file;

        $self->warning_message("no faidx file at $idx_file!");

        my $lock = Genome::Sys::LockProxy->new(
            resource => 'reference-sequence-' . $self->id . '-faidx',
            scope => 'site',
        )->lock(
            max_try       => 2,
        );
        unless ($lock) {
            $self->error_message("Failed to lock resource: $data_dir");
            return;
        }

        my $rv = Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files  => [$fa_file],
            output_files => [$idx_file],
        );

        unless ($lock->unlock()) {
            $self->error_message("Failed to unlock resource: " . $lock->resource);
            return;
        }
        unless ($rv == 1) {
            $self->error_message("Failed to run samtools faidx on fasta: $fa_file");
            return;
        }
    }
    return $idx_file if -e $idx_file;
    return;
}

sub description {
    my $self = shift;
    my $path = $self->data_directory . '/description';
    unless (-e $path) {
        return 'all';
    }
    my $fh = IO::File->new($path);
    my $desc = $fh->getline;
    chomp $desc;
    return $desc;
}

sub external_url {
    my $self = shift;
    my $url = 'https://genome.wustl.edu/view/genome/model/build/reference-sequence/consensus.fasta?id=' . $self->id;
    $url .= "/".$self->name."/all_sequences.bam";
    return $url;
}

sub sequence_dictionary_path {
    my $self = shift;
    my $file_type = shift;

    return File::Spec->join($self->data_directory, 'seqdict', "seqdict.$file_type");
}

sub get_sequence_dictionary {
    my $self = shift;
    my $file_type = shift;
    my $species = shift;
    my $picard_version = shift;
    my $create_ok = shift;
    $create_ok = 1 unless defined $create_ok;

    my $picard_path = Genome::Model::Tools::Picard->path_for_picard_version($picard_version);

    my $seqdict_dir_path = $self->data_directory.'/seqdict';
    my $path = $self->sequence_dictionary_path($file_type);

    if (-s $path) {
        return $path;
    } elsif (not $create_ok) {
        return;
    }

    $self->warning_message("No seqdict at path $path.  Creating...");

    my %new_seqdict = (
        resource => 'reference-sequence-' . $self->id . '-seqdict',
        scope => 'site',
    );
    my $lock = Genome::Sys::LockProxy->new(%new_seqdict)->lock(
        max_try => 2,
    );

    # if it couldn't get the lock after 2 tries, pop a message and keep trying as much as it takes
    unless ($lock) {
        $self->status_message("Couldn't get a lock after 2 tries, waiting some more...");
        $lock = Genome::Sys::LockProxy->new(%new_seqdict)->lock();
        unless($lock) {
            $self->error_message("Failed to lock resource: $seqdict_dir_path");
            return;
        }
    }

    $self->status_message("Failed to find sequence dictionary file at $path.  Generating one now...");
    my $seqdict_dir = $self->data_directory."/seqdict/";
    my $cd_rv =  Genome::Sys->create_directory($seqdict_dir);
    if ($cd_rv ne $seqdict_dir) {
        $self->error_message("Failed to to create sequence dictionary directory for $path. Quitting");
        return;
    }

    my $append_ref = $self->append_to;
    my $remap = $self->primary_consensus_path('fa.remap');

    if ($append_ref and -s $remap) {
        $self->status_message("Detected a remap file, and we're appending to another build. We'll skip sequence dictionary creation for the remap file and just copy the sequence dictionary from the reference build we're appending to.");
        Genome::Sys->copy_file($append_ref->get_sequence_dictionary($file_type, $species, $picard_version), $path);

        unless ($lock->unlock()) {
            $self->error_message("Failed to unlock resource: " . $lock->resource);
            return;
        }
    } else {
        my $uri = $self->sequence_uri;
        if (!$uri) {
            $self->warning_message("No sequence URI defined on this model!  Using generated default: " . $self->external_url);
            $uri = $self->external_url;
        }
        my $ref_seq = $self->full_consensus_path('fa');
        my $assembly_name = $self->assembly_name;

        # fall back to the build name if the assembly name came up short.
        if (!$assembly_name) {
            $assembly_name = $self->name;
        }

        my $create_seq_dict_cmd = "java -Xmx4g -XX:MaxPermSize=256m -cp $picard_path/CreateSequenceDictionary.jar net.sf.picard.sam.CreateSequenceDictionary R='$ref_seq' O='$path' URI='$uri' species='$species' genome_assembly='$assembly_name' TRUNCATE_NAMES_AT_WHITESPACE=true";

        my $csd_rv = Genome::Sys->shellcmd(cmd=>$create_seq_dict_cmd);

        unless ($lock->unlock()) {
            $self->error_message("Failed to unlock resource: " . $lock->resource);
            return;
        }

        if ($csd_rv ne 1) {
            $self->error_message("Failed to to create sequence dictionary for $path. Quitting");
            return;
        }
    }

    $self->reallocate_disk_allocations;

    return $path;
}

sub buckets {
    my $self = shift;
    my $users = shift;

    my $buckets = Genome::Model::Build::ReferenceSequence::Buckets->get_with_lock(
        reference_sequence_build => $self,
        users => $users || Genome::SoftwareResult::User->user_hash_for_build($self),
    );

    return $buckets;
}

sub get_by_name {
    my ($class, $name) = @_;

    unless ( $name ) {
        Carp::confess('No build name given to get imported reference sequence build');
    }

    # we now record the build name explicitly so we can do faster lookups, so this method should not be needed
    # try to get the build by name, and only continue through heuristic logic if it fails
    my $new = $class->get(name => $name);
    if ($new) {
        return $new;
    }

    # This method is not adequate as spaces are substitued in the model anme and version
    #  when creating the build name. But we'll try.
    my ($model_name, $build_version) = $name =~ /^(.+)-build(.*?)$/;
    if ( not defined $model_name ) {
        $class->status_message("Could not parse out model name and build version from build name: $name");
        return;
    }

    $class->status_message("Getting imported reference sequence builds for model ($model_name) and version ($build_version)");

    my $model = Genome::Model::ImportedReferenceSequence->get(name => $model_name);
    if ( not $model ) {
        # ok - model name may have spaces that were sub'd for dashes
        $class->status_message("No imported reference sequence model with name: $model_name");
        return;
    }

    $class->status_message("Getting builds for imported reference sequence model: ".$model->__display_name__);

    my @builds = $model->builds;
    if ( not @builds ) {
        Carp::confess("No builds for imported reference sequence model: ".$model->__display_name__);
    }

    unless($build_version) {
        my @builds_without_version;
        for my $build (@builds) {
            next if defined $build->version;

            push @builds_without_version, $build;
        }

        unless (scalar @builds_without_version > 0) {
            Carp::confess("No builds found with no version for imported reference sequence model: ".$model->__display_name__);
        }
        if ( @builds_without_version > 1 ) {
            Carp::confess("Multiple builds with no version found for model: ".$model->__display_name__);
        }

        return $builds_without_version[0];
    } else {
        my @builds_with_version;
        for my $build ( @builds ) {
            my $version = $build->version;
            if ( not defined $version or $version ne $build_version ) {
                next;
            }
            push @builds_with_version, $build;
        }
        if ( not @builds_with_version ) {
            Carp::confess("No builds found with version $build_version for imported reference sequence model: ".$model->__display_name__);
        }
        elsif ( @builds_with_version > 1 ) {
            Carp::confess("Multiple builds with version $build_version found for model: ".$model->__display_name__);
        }

        return $builds_with_version[0];
    }
}

sub chromosome_array_ref {
    my $self = shift;
    my %params = @_;

    my $format = delete($params{format});
    unless ($format) { $format = 'sam'; }

    my $species = delete($params{species});
    unless ($species) { $species = $self->species_name; }

    my $picard_version = delete($params{picard_version});
    unless ($picard_version) { $picard_version = '1.36'; }

    my $create_if_necessary = delete($params{create_seqdict});
    unless (defined $create_if_necessary) { $create_if_necessary = 1; }

    my $seq_dict = $self->get_sequence_dictionary($format,$species,$picard_version,$create_if_necessary);
    return unless $seq_dict;

    my $cmd = Genome::Model::Tools::BioSamtools::ListChromosomes->create(
        input_file => $seq_dict,
    );
    $cmd->execute or die $self->error_message('Failed to execute chromosome lister.');
    return $cmd->chromosome_array_ref;
}

sub cached_full_consensus_path {
    my ($self, $format) = @_;

    unless ($self->local_cache_basedir and -d $self->local_cache_basedir) {
        $self->status_message('Using original full consensus path.');

        return $self->full_consensus_path($format);
    }

    $self->status_message('Using cached full consensus path');

    if ( not $format or $format ne 'fa' ) {
        $self->error_message('Unsupported format ('.($format ? $format : 'none').') to get cached full consensus path');
        return;
    }

    my $cache_dir = $self->verify_or_create_local_cache;
    return if not $cache_dir;

    my $file = $cache_dir.'/all_sequences.fa';
    unless (-e $file){
        $self->error_message("Failed to find " . $file);
        return;
    }

    $self->status_message('Cached directory: '.$cache_dir);
    $self->status_message('Cached full consensus path: '.$file);

    return $file;
}

sub local_cache_basedir {
    return Genome::Config::get('fs_local_network_cache');
}

sub local_cache_dir {
    my $self = shift;
    # data_directory is usually (always?) an absolute path so remove leading / if it exists before cating with /
    my $data_directory = $self->data_directory;
    $data_directory =~ s/^\///;
    return $self->local_cache_basedir . "/" . $data_directory;
}

sub local_cache_lock {
    my $self = shift;
    return "LOCK-".$self->id;
}

sub available_kb {
    my ($self, $directory) = @_;

    Carp::confess('No directory to get available kb!') if not $directory;
    Carp::confess("Directory ($directory) does not exist! Cannot get available kb!") if not -d $directory;

    $self->status_message('Get available kb for '.$directory);

    my $df = Filesys::Df::df($directory);

    if ( not defined $df->{bavail} ) {
        $self->error_message('Failed to get kb available from df command');
        return;
    }

    $self->status_message('KB available: '. $df->{bavail});

    return $df->{bavail};
}

sub copy_file {
    my ($self, $file, $dest) = @_;

    Carp::confess('No file to copy!') if not $file;
    $self->status_message('File: '.$file);
    my $sz = -s $file;
    Carp::confess("File ($file) does not exist!") if not $file;
    $self->status_message('Size: '.$sz);

    Carp::confess('No destination file!') if not $dest;
    $self->status_message('Destination file: '.$dest);
    my $dest_sz = -s $dest;
    Carp::confess("Destination file ($dest) already exists!") if defined $dest_sz;
    $self->status_message('Destination size: '.( $dest_sz || 0));

    #Look for a gzipped version of the file. If it is present, use zcat to xfer the compressed file,
    # rather than the uncompressed.
    my $gzipped_file = $file.".gz";
    if(-s $gzipped_file){
        $self->status_message("Found a gzipped version of: $file, preparing to zcat this file into place.");

        #zcat command line
        my $zcat_cmd = "zcat ".$gzipped_file." > ".$dest;
        $self->status_message("Begin zcat of ".$file." at ".localtime()."\n");
        my $cmd = Genome::Sys->shellcmd(
            cmd => $zcat_cmd,
        );
        $self->status_message("Completed zcat of ".$file." at ".localtime()."\n");

        unless($cmd){
            die $self->error_message("Could not complete shellcmd to zcat: $gzipped_file into $dest.");
        }
    }
    #No gzipped version was found, use File::Copy::copy
    else {
        $self->status_message("Copy $file to $dest");
        my $cp = File::Copy::copy($file, $dest);
        if ( not $cp ) {
            $self->error_message('Copy failed! '.$@);
            return;
        }
        $self->status_message("Completed copy of $file.");
    }

    $dest_sz = -s $dest;
    $self->status_message('Destination size: '.( $dest_sz || 0));
    if ( $dest_sz != $sz ) {
        $self->error_message('Copy returned ok, but source/destination files are not the same size!');
        return;
    }

    return 1;
}
#

sub verify_or_create_local_cache {
    my $self = shift;

    $self->status_message('Verify or create local cache');

    my $local_cache_dir = $self->local_cache_dir;
    $self->status_message('Local cache directory: '.$local_cache_dir);
    if(  $self->_local_cache_dir_is_verified ) {
        $self->status_message('Local cache directory has already been verified!');
        return $local_cache_dir;
    }

    # lock
    $self->status_message('Lock local cache directory');
    my $lock_name = $self->local_cache_lock;
    $self->status_message('Lock name: '.$lock_name);
    my $lock = Genome::Sys::LockProxy->new(
        resource => $lock_name,
        scope => 'host',
    )->lock(
        max_try => 20, # 20 x 180 sec each = 1hr
        block_sleep => 180,
    );
    unless ($lock) {
        $self->error_message("Failed to get lock for %s!", $lock->resource);
        return;
    }
    $self->status_message('Lock obtained');

    # create the cache dir, if necessary
    if ( not -d $local_cache_dir ) {
        $self->status_message('Create local cache directory');
        my $mk_path = Genome::Sys->create_directory($local_cache_dir);
    }

    # verify files
    $self->status_message('Verify files in local cache directory');
    my %files_to_copy_to_local_cache = (
        # file_base_name => is_requried
        'all_sequences.fa' => 1,
        'all_sequences.fa.fai' => 1,
        'all_sequences.dict' => 0,
    );
    my $dir = $self->data_directory;
    my @files_to_copy;
    my $kb_required = 0;
    for my $base_name ( keys %files_to_copy_to_local_cache ) {
        my $file = $dir.'/'.$base_name;
        $self->status_message('File: '.$file);
        my $sz = -s $file;
        if ( not $sz ) {
            next if not $files_to_copy_to_local_cache{$base_name}; # not required
            $self->error_message("File ($file) to copy to cache directory does not exist!");
            return;
        }
        $self->status_message('Size: '.$sz);
        my $cache_file = $local_cache_dir.'/'.$base_name;
        my $cache_sz = -s $cache_file;
        $cache_sz ||= 0;
        $self->status_message('Cache file: '.$cache_file);
        $self->status_message('Cache size: '.$cache_sz);
        next if $cache_sz and $cache_sz == $sz;
        $self->status_message('Need to copy file: '.$cache_file);
        unlink $cache_file if $cache_sz > 0; # partial???
        push @files_to_copy, $base_name;
        $kb_required += sprintf('%.0d', $sz / 1024);
    }
    $self->status_message('Verify...OK');

    # copy
    if ( @files_to_copy ) {
        $self->status_message('Copy '.@files_to_copy.' to local cache directory');
        my $kb_available = $self->available_kb($local_cache_dir);
        $self->status_message('KB required: '.$kb_required);
        if ( $kb_required and $kb_required > $kb_available ) {
            $self->error_message("Cannot copy files to cache directory ($local_cache_dir). The kb required ($kb_required) is greater than the kb available ($kb_available).");
            return;
        }
        for my $base_name ( @files_to_copy ) {
            my $file = $dir.'/'.$base_name;
            my $cache_file = $local_cache_dir.'/'.$base_name;
            return if not $self->copy_file($file, $cache_file);
        }
        $self->status_message('Copy...OK');
    }

    $self->_local_cache_dir_is_verified(1);

    # rm lock
    $self->status_message('Remove lock');
    my $rv = eval { $lock->unlock() };
    if ( not $rv ) {
        $self->error_message("Failed to unlock (%s): %. But we don't care.", $lock->resource, $@);
    }

    $self->status_message('Verify or create local cache...OK');

    return $local_cache_dir;
}


sub get_or_create_genome_file {
    my $self = shift;

    my $genome_file = $self->data_directory .'/all_sequences.genome';
    unless (-s $genome_file) {
        my $chromosomes_with_lengths = $self->chromosomes_with_lengths;

        my $genome_fh = Genome::Sys->open_file_for_writing($genome_file);
        unless ($genome_fh) {
            die('Failed to open genome file for writing: '. $genome_file);
        }

        for my $chr (@$chromosomes_with_lengths) {
            $genome_fh->say( join("\t", @$chr) );
        }

        $genome_fh->close;
    }
    return $genome_file;
}

sub chromosomes_with_lengths {
    my $self = shift;

    my @chr_lengths;

    my $chr_data = $self->get_or_create_seqdict_hash_ref;

    for my $chr (keys %{$chr_data}) {
        push @chr_lengths, [$chr, $chr_data->{$chr}{length}];
    }

    return \@chr_lengths;
}

sub get_or_create_seqdict_hash_ref {
    my $self = shift;

    return $self->_seqdict if $self->_seqdict;

    my $seqdict = $self->sequence_dictionary_path('sam');
    unless(-s $seqdict) {
        $self->fatal_message('No sequence dictionary found!');
    }
    my %seqdict_fields = (
        'name' => {
            key => 'SN',
            position => 2,
        },
        'length' => {
            key => 'LN',
            position => 3,
        },
        'uri' => {
            key => 'UR',
            position => 4,
        },
        'assembly_name' => {
            key => 'AS',
            position => 5,
        },
        'md5' => {
            key => 'M5',
            position => 6,
        },
        species_name => {
            key => 'SP',
            position => 7,
        },
    );
    my @field_names = sort {$seqdict_fields{$a}{position} <=> $seqdict_fields{$b}{position}} keys %seqdict_fields;
    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        separator => "\t",
        input => $seqdict,
        headers => ['tag', @field_names],
        allow_extra_columns => 1,
        ignore_lines_starting_with => '(?!@SQ)',
    );

    my %seqdict;
    while (my $line = $parser->next) {
        unless ($line->{tag} eq '@SQ') {
            Carp::confess 'parser error';
        }
        my $chr;
        for my $field_name (@field_names) {
            my $value = $self->_parse_seqdict_field($seqdict_fields{$field_name}{key},$seqdict_fields{$field_name}{position},$line->{$field_name});
            if ($field_name eq 'name') {
                $chr = $value;
            }
            $seqdict{$chr}{$field_name} = $value;
        }
    }

    $self->_seqdict(\%seqdict);

    return $self->_seqdict;
}

sub _parse_seqdict_field {
    my $self = shift;
    my ($expected_key,$expected_pos,$field) = @_;
    my ($observed_key,$value) = split(':',$field,2);
    unless ($expected_key eq $observed_key) {
        $self->fatal_message('Expected %s key as %s field but got %s',$expected_key,$expected_pos,$observed_key);
    }
    return $value;
}

sub is_superset_of {
    my ($self, $other_refbuild) = @_;

    my $my_chromosomes = Set::Scalar->new(@{$self->chromosome_array_ref});
    my $other_chromosomes = Set::Scalar->new(@{$other_refbuild->chromosome_array_ref});
    return $my_chromosomes >= $other_chromosomes;
}

sub contains {
    my ($self, $other_refbuild) = @_;

    return 1 if $self eq $other_refbuild;

    if(($self->coordinates_from || '') eq $other_refbuild or $self->is_derived_from($other_refbuild)) {
        return $self->is_superset_of($other_refbuild);
    }

    return 0;
}

# Given a feature list accessor, try to get it from myself or my ancestors
sub get_feature_list {
    my ($self, $feature_list_accessor, @ancestry_stack) = @_;
    push @ancestry_stack, $self;

    my $feature_list = $self->$feature_list_accessor;
    if (not defined $feature_list) {
        if ($self->derived_from) {
            $self->debug_message("Could not get_feature_list with accessor (%s) on reference sequence build (%s)... looking at the next ancestor...", $feature_list_accessor, $self->name);
            $feature_list = $self->derived_from->get_feature_list($feature_list_accessor, @ancestry_stack);
        } else {
            $self->error_message("Reference sequence (%s) does not have any parent reference sequence", $self->name);
        }
    }

    # When we have finished our recursion, make sure that the feature list's reference is compatible with this reference
    if (scalar(@ancestry_stack) == 1 and $ancestry_stack[0]->id eq $self->id and $feature_list) {
        unless ($self->is_superset_of($feature_list->reference)) {
            # This case presents a problem because if the feature list is a superset of our coordinates, it will be nonsense if applied to us.
            die $self->error_message("Reference sequence (%s) is not a superset of the reference sequence on the feature list (%s)", $self->name, $feature_list->name);
        }
        return $feature_list;
    } elsif (scalar(@ancestry_stack) == 1 and $ancestry_stack[0]->id eq $self->id and not $feature_list) {
        die $self->error_message("Could not get_feature_list with accessor (%s) on reference sequence build (%s) or any of its ancestors", $feature_list_accessor, $self->name);
    } else {
        return $feature_list
    }
}

sub combined_references {
    my $class = shift;
    my @references = @_;

    my $reference_set = Set::Scalar->new(map $_->id, @references);

    my @combined_reference_inputs = Genome::Model::Build::Input->get(name => 'combines', value_id => [$reference_set->members]);
    my @combined_references = Genome::Model::Build->get([map $_->build_id, @combined_reference_inputs]);

    #filter to only those that combine exactly all our references
    my @exact_combined_references = grep { $reference_set == Set::Scalar->new(map $_->id, $_->combines) } @combined_references;

    unless(@exact_combined_references) {
        #see if one of the provided references is a combination of the others.
        for my $reference (@references) {
            my $remaining_reference_set = $reference_set - $reference->id;
            if($remaining_reference_set <= Set::Scalar->new(map $_->id, $reference->combines)) {
                push @exact_combined_references, $reference;
            }
        }
    }

    return @exact_combined_references;
}

sub per_chromosome_fastas {
    my $self = shift;
    my $users = shift;

    my $per_chromosome_fastas = Genome::Model::Build::ReferenceSequence::PerChromosomeFastas->get_or_create(
        reference_sequence_build => $self,
        users => $users || Genome::SoftwareResult::User->user_hash_for_build($self),
    );
    return $per_chromosome_fastas;
}

sub _disk_usage_result_subclass_names {
    my $self = shift;

    return [qw(
        Genome::Model::Build::ReferenceSequence::Buckets
        Genome::Model::Build::ReferenceSequence::AlignerIndex
        Genome::Model::Build::ReferenceSequence::AnnotationIndex
    )];
}

1;
