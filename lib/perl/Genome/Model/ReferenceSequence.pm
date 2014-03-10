package Genome::Model::ReferenceSequence;
use strict;
use warnings;
use Genome;
use File::Spec;
use File::Temp;

# this ensures that, when a generic UR::Value of one or zero is gotten, 
# other subclasses of UR::Value
UR::Value::Number->get([0,1]);

# UR::Value has multiple sub-classes with the same identity for different instances.
# Until a recent patch to how UR::Value loads, an attempt to get UR::Value::Boolean after a generic UR::Value

class Genome::Model::ReferenceSequence {
    is => 'Genome::ModelDeprecated',
    has => [
        processing_profile => {
            # hard-coded to use the base class to get around the ::Imported vs not name
            is => 'Genome::ProcessingProfile',
            id_by => 'processing_profile_id',
        },
    ],
    has_optional_input => [
        fasta_file => {
            is => 'FilesystemPath',
            doc => "fully qualified fasta filename to copy to all_sequences.fa in the build's data_directory.",
        },
        sequence_uri => {
            is => 'Text',
            doc => "publicly available URI to the sequence file for the fasta",
        },
        
        # ${modelname}-build${buildname}
        build_name => {
            is => 'Text',
            doc => "human meaningful name of this build",
        },

        # these three combine to form the default build name, 
        # along with the species
        prefix => { 
            is => 'Text', 
            doc => 'The prefix/source of the sequence (such as GRC).  May not contain spaces.',
        },
        desc => {
            is => 'Text',
            doc => 'Optional additional qualifier, i.e. "lite" for the GRC-human-lite refseq model.',
        },
        version => {
            is => 'Text',
            doc => 'Identifies the version of the reference sequence.  This string may not contain spaces.',
        },

        assembly_name => {
            is => 'Text',
            is_deprecated => 1,
            doc => "when the assembly name is not the name of the build, set this",
        },
        generate_sequence_uri => {
            is => 'Boolean',
            #is_transient => 1,
            #is_optional => 1,
            #default_value => 0,
        },
        skip_bases_files => {
            is => 'Boolean', 
            default_value => 1,
            doc => 'If specified, individual bases files are not created for each sequence in the fasta',
        },
        is_rederivable => {
            is => 'Boolean',
            default_value => 0,
            doc => 'If true, indicates that the reference is the product of other analysis and could be rederived. It will be stored as a model/build product rather than imported data.'
        },

    ],
    has_optional => [
        species_name => {
            is => 'Text',
            via => 'subject',
            to => 'name',
        },
        _expected_name => {
            is => 'Text',
            calculate_from => ['prefix','subject_name','desc'],
            calculate => 'no warnings; my $v = ($desc ? "$prefix-$subject_name-$desc" : "$prefix-$subject_name"); $v =~ s/ /_/g; $v'
        },
        sample_names => {
            is => 'ARRAY',
            calculate => q{ return; }
        },
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            id_by => 'processing_profile_id',
        },
    ],
    doc => 'a versioned reference sequence, with cordinates suitable for annotation',
};

sub _has_legacy_input_types { 1 };

sub build_by_version {
    my $self = shift;
    my $version = shift;
    if ($version eq '36' or $version eq '37') {
        # this is present only to help developers troubleshoot problems related to poorly maintained refseq data
        # remove it when this data is cleaned-up (ssmith-2013-02-09)
        Carp::cluck("Searching for version $version is likely due "
                  . "to out-of-date code.  Reference versions now "
                  . "end in -lite, or a patch level indicator in most cases");
    }
    my @b = Genome::Model::Build::ImportedReferenceSequence->get(
        'version' => $version,
        'model_id' => $self->genome_model_id
    );
    if (@b > 1) {
        die "Multiple builds for version $version for model " . $self->genome_model_id;
    }
    return $b[0];
}

sub build_needed {
    my $self = shift;

    #These models always have one build per "version" of the reference.
    #The models themselves do not store information relevant to whether a future build is needed.
    return 0;
}

sub _help_detail_for_model_define {
    return <<EOS
  # for a regular, complete fasta
  genome model define imported-reference-sequence \\
    --prefix=GRClite \\
    --species-name=human \\
    --version=37 \\
    --use-default-sequence-uri
    --fasta-file=complete.fa \\

  # for a reference which extends another
  # (a remap file can have the same name as the fasta and be next-to it)
  genome model define imported-reference-sequence \\
    --prefix=GRC \\
    --species-name=human \\
    --version=37-p8-test7 \\
    --use-default-sequence-uri
    --fasta-file=additions.fa \\
    --append-to 106942997

EOS
}

sub _resolve_resource_requirements_for_build {
    return "-R 'rusage[gtmp=10:mem=6000]' -M 6000000";
}

sub _resolve_disk_group_name_for_build {
    return ($_[0]->is_rederivable ? $ENV{GENOME_DISK_GROUP_MODELS} : $ENV{GENOME_DISK_GROUP_REFERENCES});
}

sub _execute_build {
    my ($self, $build) = @_;
    my $model = $build->model;

    my $fasta_size = -s $build->fasta_file;
    unless(-e $build->fasta_file && $fasta_size > 0) {
        $self->error_message("Reference sequence fasta file \"" . $build->fasta_file . "\" is either inaccessible, empty, or non-existent.");
        return;
    }

    my $build_directory = $build->data_directory;
    my $output_directory = File::Temp->newdir(
        "tmp_XXXXX",
        DIR     => $build_directory,
        CLEANUP => 1,
    );
    chmod(0775, $output_directory); #so can be manually cleaned up by others if need be

    unless ($self->_copy_fasta_file($build, $output_directory)) {
        $self->error_message("fasta copy failed.");
        return;
    }

    $self->debug_message('Promoting files to final location.');
    for my $staged_file (glob($output_directory . '/*')) {
        my ($vol, $dir, $file_base) = File::Spec->splitpath($staged_file);
        my $final_file = join('/', $build_directory, $file_base);
        rename($staged_file, $final_file);
    }

    $self->debug_message('Creating sequence dictionaries');
    return unless $self->_create_sequence_dictionary($build);

    # Reallocate to amount of space actually consumed if the build has an associated allocation and that allocation
    # has an absolute path the same as this build's data_path
    $self->debug_message("Reallocating.");
    if (defined($build->disk_allocation)) {
        unless($build->disk_allocation->reallocate) {
            $self->error_message("Reallocation failed.");
            return;
        }
    }

    #create manifest file
    unless ($self->_create_manifest_file($build)){
        $self->error_message("Could not create manifest file");
    }

    $self->debug_message("Done.");
    return 1;
}

sub _copy_fasta_file {
    my ($self, $build, $output_directory) = @_;

    my @fastas;
    my $primary_fasta_path;
    if ($build->can('append_to') && $build->append_to) {
        $primary_fasta_path = File::Spec->catfile($output_directory, 'appended_sequences.fa');
    } else {
        $primary_fasta_path = File::Spec->catfile($output_directory, 'all_sequences.fa');
    }

    $self->debug_message("Copying primary fasta file");

    #If an error occurs here about refusing to write to an existing file, that was most likely on a re-run of the build
    #and the original error can be found earlier in the logs.  To restart, clear files out of the build directory.
    my $in = Bio::SeqIO->new(-file => $build->fasta_file, -format => 'Fasta');
    my $out = Bio::SeqIO->new(-file => ">$primary_fasta_path", -format => 'Fasta');
    while(my $seq = $in->next_seq) {
        $out->write_seq($seq);
    }
    $in->close;
    $out->close;
    unless (-s $primary_fasta_path) {
        $self->error_message('Failed to copy "' . $build->fasta_file . '" to "' . $primary_fasta_path. '.');
        return;
    }
    push(@fastas, $primary_fasta_path);

    my $remap_file = $build->fasta_file . ".remap";
    if (-s $remap_file) {
        $self->debug_message("Remapping file found. Copying $remap_file...");
        my $destination = $primary_fasta_path . ".remap";
        unless (Genome::Sys->copy_file($remap_file, $destination)) {
            $self->error_message("Failed to copy $remap_file to $destination");
            return;
        }
    }

    if ($build->can('append_to') && $build->append_to) {
        $self->debug_message("Copying full fasta file");
        my $full_fasta_path = File::Spec->catfile($output_directory, 'all_sequences.fa');
        my $cmd = Genome::Model::Tools::Fasta::Concat->create(
            input_files => [
            $build->append_to->full_consensus_path('fa'),
            $build->fasta_file,
            ],
            output_file => $full_fasta_path,
        );
        unless ($cmd->execute()) {
            $self->error_message("Failed to concatenate fasta files");
            return;
        }

        push(@fastas, $full_fasta_path);
    }

    $self->debug_message("Doing samtools faidx.");
    my $samtools_path = Genome::Model::Tools::Sam->path_for_samtools_version(); #uses default version if none passed

    for my $fasta (@fastas) {
        my $samtools_cmd = sprintf('%s faidx %s', $samtools_path, $fasta);
        my $rv = Genome::Sys->shellcmd(
            cmd => $samtools_cmd,
            input_files => [$fasta],
        );

        unless($rv) {
            $self->error_message("samtools faidx failed for $fasta.");
            return;
        }
    }

    my $sequence_count = `wc -l < $primary_fasta_path.fai`;
    chomp $sequence_count;

    if ($build->skip_bases_files) {
        $self->debug_message("Skipping creation of $sequence_count bases files for fasta.");
    } elsif ($sequence_count > 1024) {
        $self->debug_message("You have requested creation of bases files, but there are $sequence_count sequences in the primary fasta file. Creating bases files when there are thousands or millions of sequences is potentially dangerous.")
    } else {
        $self->debug_message("Making $sequence_count bases files from fasta.");
        my $rv = $self->_make_bases_files($primary_fasta_path, $output_directory);

        unless($rv) {
            $self->error_message('Making bases files failed.');
            return;
        }
    }

    return 1;
}

sub _create_sequence_dictionary {
    my ($self, $build) = @_;

    my $abs_seqdict_path = $build->get_sequence_dictionary('sam', $build->species_name, '1.46');
    unless (-e $abs_seqdict_path) {
        $self->error_message(sprintf('Failed to create sequence dictionary (with $build->species_name=%s)',
                $build->species_name));
        return;
    }
    my $build_data_directory = $build->data_directory;
    (my $seqdict_path = $abs_seqdict_path) =~ s/$build_data_directory\/*//;

    my $dict_path = File::Spec->join($build->data_directory, 'all_sequences.dict');
    Genome::Sys->create_symlink($seqdict_path, $dict_path);
}

# This is a simplified version of the previous code for
# finding chromosome names in fasta files, and splitting the content out into .bases files
# It is assumed that the sequence names are indicated via the ">" character rather than the ";" character

sub _make_bases_files {
    my $self = shift;
    my ($fa,$output_dir) = @_;

    my $bases_dir = join('/', $output_dir, 'bases');
    Genome::Sys->create_directory($bases_dir);

    my $fafh = Genome::Sys->open_file_for_reading($fa);
    unless($fafh){
        $self->error_message("Could not open file $fa for reading.");
        die $self->error_message;
    }
    my @chroms;
    my $file;
    while(<$fafh>){
        my $line = $_;
        chomp($line);
        #if the line contains a sequence name, check that name
        if($line =~ /^>/){
            my $chr = $';
            my @rest;
            ($chr, @rest) = split " ",$chr;

            if(length join(" ", @rest) >= 1024) { #see bns_restore_core in bntseq.c of bwa
                $self->warning_message('The extra information on the header for sequence ' . $chr . ' is longer than the maximum size handled by bwa.  This may cause failures if this reference is used for alignment.');
            }

            $chr =~s/(\/|\\)/_/g;  # "\" or "/" are not allowed in sequence names
            push @chroms, $chr;
            if(defined($file)) {
                $file->close;
            }

            my $file_name = join('/', $bases_dir, $chr . ".bases");
            $file = Genome::Sys->open_file_for_writing($file_name);
            unless($file){
                $self->error_message("Could not open file " . $file_name . " for reading.");
                die $self->error_message;
            }
            next;
        }
        print $file $line;

    }
    $file->close;
    $fafh->close;
}

sub _create_manifest_file {
    my $self = shift;
    my $build = shift;

    my $manifest_path = $build->manifest_file_path;
    if (-z $manifest_path){
        $self->warning_message('Manifest file already exists!');
        return;
    }

    my $manifest_fh = IO::File->new($manifest_path, 'w');
    unless ($manifest_fh){
        $self->error_message("Could not open manifest file path, exiting");
        die();
    }

    my @files = $self->_list_bases_files($build);
    for my $file (@files){
        $manifest_fh->print($self->_create_manifest_file_line($file), "\n");
    }

    $manifest_fh->close;
    return 1;
}

sub _create_manifest_file_line {
    my $self = shift;
    my $file = shift;

    my $file_size = -s ($file);
    my $md5 = Genome::Sys->md5sum($file);
    return join("\t", $file, $file_size, $md5);
}

sub _list_bases_files {
    my $self = shift;
    my $build = shift;

    my $data_dir = $build->data_directory;
    my $fa = $build->fasta_file;
    my $bases_dir = join('/', $data_dir, 'bases');
    $bases_dir = $data_dir unless -e $bases_dir; #some builds lack a bases directory

    my @bases_files;

    my $sequence_count = 0;
    if (-e "$fa.fai") {
        my $sequence_count = `wc -l < $fa.fai`;
        chomp $sequence_count;
        #$self->debug_message("There were $sequence_count bases in $fa.fai")
    }

    if (!($build->skip_bases_files) and !($sequence_count > 1024) and $fa and -e $fa){
        my $fafh = Genome::Sys->open_file_for_reading($fa);
        unless($fafh){
            $self->error_message("Could not open file $fa for reading.");
            die $self->error_message;
        }
        while(<$fafh>){
            my $line = $_;
            chomp($line);
            #if the line contains a sequence name, check that name
            if($line =~ /^>/){
                my $chr = $';
                ($chr) = split " ",$chr;
                $chr =~s/(\/|\\)/_/g;  # "\" or "/" are not allowed in sequence names
                #$chr=~s/(\.1)//; #handle chrom names that end in .1
                my $path =  join("/", $bases_dir, "$chr.bases");
                unless (-e $path) {
                    warn "no bases file $path found";
                    next;
                }
                push(@bases_files, $path);
            }
        }
    }else{
        @bases_files = glob($bases_dir . "/*.bases")
    }

    return @bases_files;
}

1;

