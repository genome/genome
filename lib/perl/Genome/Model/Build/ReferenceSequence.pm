package Genome::Model::Build::ReferenceSequence;
use strict;
use warnings;
use Genome;
use File::Path;
use File::Copy;

require Carp;
use Regexp::Common;
use POSIX;

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
            doc => 'file handle per chromosome for reading sequences so that it does not need to be constantly closed/opened',
        },
        _local_cache_dir_is_verified => { is => 'Boolean', default_value => 0, is_optional => 1, },

    ],

    has_optional_input => [
        build_name => {
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
        },
    ],

    doc => 'a specific version of a reference sequence, with cordinates suitable for annotation',
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

sub get{
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
    if ($self->id eq "10194788") {
        return "/gscmnt/sata921/info/medseq/make_tier_bed_files/hg18_build36_ucsc_files";
    }
    elsif ($self->id eq "106942997") {
        return "/gscmnt/sata921/info/medseq/make_tier_bed_files/NCBI-human-build37/hg19_files";
    }
    elsif ($self->id eq "107494762") {
        return "/gscmnt/sata921/info/medseq/make_tier_bed_files/NCBI-mouse-build37/mm9";
    }
    else {
        $self->status_message("UCSC Tiering Directory not currently available for this species: ".$self->species_name);
        return;
    }
}

#TODO: figure out how to get this and make it a real result
sub get_or_create_ucsc_conservation_directory {
    my $self = shift;
    if ($self->id eq "10194788") {
        return "/gscmnt/ams1161/info/model_data/2771411739/build113115679/annotation_data/ucsc_conservation/";
    }
    elsif ($self->id eq "106942997") {
        return "/gscmnt/ams1102/info/model_data/2771411739/build106409619/annotation_data/ucsc_conservation/"
    }
    else {
        $self->status_message("UCSC conservation scores are not currently available in the system for this species: ".$self->species_name);
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
#SamToBam convertion
sub full_consensus_sam_index_path {
    my $self        = shift;
    my $sam_version = shift;

    my $data_dir = $self->data_directory;
    my $fa_file  = $self->full_consensus_path('fa');
    my $idx_file = $fa_file.'.fai';

    unless (-e $idx_file) {
        my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($sam_version);
        my $cmd      = $sam_path.' faidx '.$fa_file;
        
        $self->warning_message("no failx file at $idx_file!");

        my $lock = Genome::Sys->lock_resource(
            resource_lock => $data_dir.'/lock_for_faidx',
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

        unless (Genome::Sys->unlock_resource(resource_lock => $lock)) {
            $self->error_message("Failed to unlock resource: $lock");
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

sub get_sequence_dictionary {
    my $self = shift;
    my $file_type = shift;
    my $species = shift;
    my $picard_version = shift;

    my $picard_path = Genome::Model::Tools::Picard->path_for_picard_version($picard_version);

    my $seqdict_dir_path = $self->data_directory.'/seqdict';
    my $path = "$seqdict_dir_path/seqdict.$file_type";

    if (-s $path) {
        return $path;
    } 

    $self->warning_message("No seqdict at path $path.  Creating...");

    #lock seqdict dir here
    my $lock = Genome::Sys->lock_resource(
        resource_lock => $seqdict_dir_path."/lock_for_seqdict-$file_type",
        max_try       => 2,
    );

    # if it couldn't get the lock after 2 tries, pop a message and keep trying as much as it takes
    unless ($lock) {
        $self->status_message("Couldn't get a lock after 2 tries, waiting some more...");
        $lock = Genome::Sys->lock_resource(resource_lock => $seqdict_dir_path."/lock_for_seqdict-$file_type");
        unless($lock) {
            $self->error_message("Failed to lock resource: $seqdict_dir_path");
            return;
        }
    }

    $self->status_message("Failed to find sequence dictionary file at $path.  Generating one now...");
    my $seqdict_dir = $self->data_directory."/seqdict/";
    my $cd_rv =  Genome::Sys->create_directory($seqdict_dir);
    if ($cd_rv ne $seqdict_dir) {
        $self->error_message("Failed to to create sequence dictionary directory for $path. Quiting");
        return;
    }

    my $append_ref = $self->append_to;
    my $remap = $self->primary_consensus_path('fa.remap');

    if ($append_ref and -s $remap) {
        $self->status_message("Detected a remap file, and we're appending to another build. We'll skip sequence dictionary creation for the remap file and just copy the sequence dictionary from the reference build we're appending to.");
        Genome::Sys->copy_file($append_ref->get_sequence_dictionary($file_type, $species, $picard_version), $path);

        unless (Genome::Sys->unlock_resource(resource_lock => $lock)) {
            $self->error_message("Failed to unlock resource: $lock");
            return;
        }
    } else {
        #my $picard_path = "$ENV{GENOME_SW_LEGACY_JAVA}/samtools/picard-tools-1.04/";
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

        unless (Genome::Sys->unlock_resource(resource_lock => $lock)) {
            $self->error_message("Failed to unlock resource: $lock");
            return;
        }

        if ($csd_rv ne 1) {
            $self->error_message("Failed to to create sequence dictionary for $path. Quiting");
            return;
        }
    }

    $self->reallocate;

    return $path;
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

    my $seq_dict = $self->get_sequence_dictionary($format,$species,$picard_version);

    my $cmd = Genome::Model::Tools::BioSamtools::ListChromosomes->create(
        input_file => $seq_dict,
    );
    $cmd->execute or die $self->error_message('Failed to execute chromosome lister.');
    return $cmd->chromosome_array_ref;
}

sub cached_full_consensus_path {
    my ($self, $format) = @_;

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
    return $ENV{GENOME_FS_LOCAL_NETWORK_CACHE};
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
    return $self->local_cache_basedir."/LOCK-".$self->id;
}

#MOVE TO GENOME::SYS#
sub available_kb {
    my ($self, $directory) = @_;

    Carp::confess('No directory to get available kb!') if not $directory;
    Carp::confess("Directory ($directory) does not exist! Cannot get available kb!") if not -d $directory;

    $self->status_message('Get availble kb for '.$directory);

    my $cmd = "df -k $directory |";
    $self->status_message('DF command: '.$cmd);
    my $fh = IO::File->new($cmd);
    if ( not $fh ) {
        $self->error_message('Failed to create df command: '.$!);
        return;
    }
    $fh->getline; # Filesystem           1K-blocks      Used Available Use% Mounted on
    my $line = $fh->getline;
    my @tokens = split(/\s+/, $line);
    $fh->close;

    if ( not defined $tokens[3] ) {
        $self->error_message('Failed to get kb available from df command');
        return;
    }

    $self->status_message('KB available: '.$tokens[3]);

    return $tokens[3];
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
    my $lock = Genome::Sys->lock_resource(
        resource_lock => $lock_name,
        max_try => 20, # 20 x 180 sec each = 1hr
        block_sleep => 180,
    );
    unless ($lock) {
        $self->error_message("Failed to get lock for $lock_name!");
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
    my $rv = eval { Genome::Sys->unlock_resource(resource_lock=>$lock); };
    if ( not $rv ) {
        $self->error_message("Failed to unlock ($lock): $@. But we don't care.");
    }

    $self->status_message('Verify or create local cache...OK');

    return $local_cache_dir;
}


sub get_or_create_genome_file {
    my $self = shift;

    my $genome_file = $self->data_directory .'/all_sequences.genome';
    unless (-s $genome_file) {
        my $seqdict_sam = $self->get_sequence_dictionary('sam');
        my $genome_fh = Genome::Sys->open_file_for_writing($genome_file);
        unless ($genome_fh) {
            die('Failed to open genome file for writing: '. $genome_file);
        }
        my $seqdict_fh = Genome::Sys->open_file_for_reading($seqdict_sam);
        unless ($seqdict_fh) {
            die('Failed to open seqdict file for reading: '. $seqdict_sam);
        }
        while (my $line = $seqdict_fh->getline) {
            chomp($line);
            if ($line =~ /^\@SQ\s+SN:(.+)\s+LN:(\d+)\s*/) {
                print $genome_fh $1 ."\t". $2 ."\n";
            }
        }
        $seqdict_fh->close;
        $genome_fh->close;
    }
    return $genome_file;
}

1;
