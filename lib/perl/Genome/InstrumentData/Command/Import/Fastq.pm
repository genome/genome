package Genome::InstrumentData::Command::Import::Fastq;

use strict;
use warnings;

use Genome;
use File::Copy;
use File::Basename;
use Data::Dumper;
use IO::File;

my %properties = (
    source_data_files => {
        is => 'Text',
        doc => 'source data path of import data file',
    },
    library_name => {
        is => 'Text',
        doc => 'library name, used to fetch sample name',
        is_optional => 1,
    },
    sample_name => {
        is => 'Text',
        doc => 'sample name for imported file, like TCGA-06-0188-10B-01D',
        is_optional => 1,
    },
    import_source_name => {
        is => 'Text',
        doc => 'source name for imported file, like Broad Institute',
        is_optional => 1,
    },
    species_name => {
        is => 'Text',
        doc => 'species name for imported file, like human, mouse',
        is_optional => 1,
    },
    import_format => {
        is => 'Text',
        doc => 'format of import data, like bam',
        valid_values => ['sanger fastq','solexa fastq','illumina fastq'],
        default_value => 'sanger fastq',
        is_optional => 1,
    },
    sequencing_platform => {
        is => 'Text',
        doc => 'sequencing platform of import data, like solexa',
        valid_values => ['solexa'],
        default_value => 'solexa',
        is_optional => 1,
    },
    description  => {
        is => 'Text',
        doc => 'general description of import data, like which software maq/bwa/bowtie to used to generate this data',
        is_optional => 1,
    },
    read_count  => {
        is => 'Number',
        doc => 'total read count of import data',
        is_optional => 1,
    },
    base_count  => {
        is => 'Number',
        doc => 'total base count of import data',
        is_optional => 1,
    },
    reference_name  => {
        is => 'Text',
        doc => 'reference name for imported data aligned against, if this is given, an alignment allocation will be created',
        is_optional => 1,
    },
    allocation => {
        is => 'Genome::Disk::Allocation',
        is_optional => 1,
    },
    read_length => {
        is => 'Number',
        doc => '------',
        is_optional => 1,
    },
    fwd_read_length => {
        is => 'Number',
        doc => '------',
        is_optional => 1,
    },
    rev_read_length => {
        is => 'Number',
        doc => '------',
        is_optional => 1,
    },
    fragment_count => {
        is => 'Number',
        doc => '------',
        is_optional => 1,
    },
    run_name => {
        is => 'Text',
        doc => '------',
        is_optional => 1,
    },
    subset_name => {
        is => 'Text',
        doc => 'This is the lane #. If this is not specified, it will be autopopulated. This may contain more info than just lane number, but it should preceed the lane number and be set off by hyphens.',
        is_optional => 1,
    },
    sd_above_insert_size => {
        is => 'Number',
        doc => '------',
        is_optional => 1,
    },
    median_insert_size => {
        is => 'Number',
        doc => '------',
        is_optional => 1,
    },
    is_paired_end => {
        is => 'Number',
        doc => '------',
        is_optional => 1,
    },
    sra_sample_id => {
        is=>'String',
        doc => 'SRA sample name',
        is_optional => 1,
    },
    generated_instrument_data_id=> {
        is=>'Number',
        doc => 'generated sample ID',
        is_optional => 1,
    }
);


class Genome::InstrumentData::Command::Import::Fastq {
    is  => 'Command',
    has => [%properties],
    has_optional => [
        validate_only => {
            is => 'Boolean',
            doc => 'Do not perform actual import, only validate input.',
        },
    ],
};


sub execute {
    my $self = shift;

    if ($self->validate_only) {
        if ($self->_check_fastq_integrity) {
            $self->status_message("Files successfully validated.");
            return 1;
        } else {
            $self->error_message("Unable to validate files.");
            die $self->error_message;
        }
    }

    unless(defined($self->sample_name) || defined($self->library_name)){
        $self->error_message("In order to import a fastq, a library-name or sample-name is required.");
        die $self->error_message;
    }
    my $sample;
    my $library;

    if(defined($self->sample_name) && defined($self->library_name)){
        $sample = Genome::Sample->get(name => $self->sample_name);
        $library = Genome::Library->get(name => $self->library_name);
        unless(defined($sample)){
            $self->error_message("Could not locate a sample by the name of ".$self->sample_name);
            die $self->error_message;
        }
        unless(defined($library)){
            $self->error_message("Could not locate a library by the name of ".$self->library_name);
            die $self->error_message;
        }
        unless($sample->name eq $library->sample_name){
            $self->error_message("The supplied sample-name ".$self->sample_name." and the supplied library ".$self->library_name." do not match.");
            die $self->error_message;
        }
    } elsif (defined($self->library_name)){
        $library = Genome::Library->get(name => $self->library_name);
        unless(defined($library)) {
            $self->error_message("Library name not found.");
            die $self->error_message;
        }
        $sample = Genome::Sample->get(id => $library->sample_id);
        unless (defined($sample)) {
            $self->error_message("Could not retrieve sample from library name");
            die $self->error_message;
        }
        $self->sample_name($sample->name);
    } elsif (defined($self->sample_name)){
        $sample = Genome::Sample->get(name=>$self->sample_name);
        unless(defined($sample)){
            $self->error_message("Could not locate sample with the name ".$self->sample_name);
            die $self->error_message;
        }
        # create one.
        $library = Genome::Library->get(sample_name => $sample->name, name=>$sample->name . '-extlibs');
        unless(defined($library)){
            $library = Genome::Library->create(sample_name => $sample->name, name=>$sample->name . '-extlibs');
        }
        unless (defined $library) {
            $self->error_message("COuld not locate a library associated with the sample-name ".$sample->name . "-extlibs and couldn't create one either.");
            die $self->error_message;
        }
        $self->library_name($library->name);
    } else {
        $self->error_message("Failed to define a sample or library.");
        die $self->error_message;
    }

    # gather together params to create the imported instrument data object
    my %params = ();
    for my $property_name (keys %properties) {
        unless ($properties{$property_name}->{is_optional}) {
            # required
            unless ($self->$property_name) {
                # null
                $self->error_message ("Required property: $property_name is not given");
                return;
            }
        }
        next if $property_name =~ /^(species|reference)_name$/;
        next if $property_name =~ /^source_data_files$/;
        next if $property_name =~ /^allocation$/;
        next if $property_name =~ /^library_name$/;
        next if $property_name =~ /^sample_name$/;
        $params{$property_name} = $self->$property_name if $self->$property_name;
    }

    $params{sequencing_platform} = $self->sequencing_platform;
    $params{import_format} = $self->import_format;
    $params{library_id} = $library->id;
    $params{library_name} = $library->name;
    if(defined($self->allocation)){
        $params{allocations} = $self->allocation;
    }

    $self->_check_fastq_integrity;

    unless(defined($params{read_count})){
        my ($read_count,$fragment_count) = $self->get_read_count_and_fragment_count;
        unless(defined($read_count)){
            $self->error_message("No read count was specified and none could be calculated from the fastqs");
            die $self->error_message;
        }
        unless(defined($self->fragment_count)){
            $self->fragment_count($fragment_count);
            $params{fragment_count} = $fragment_count;
        }
        $self->read_count($read_count);
        $params{read_count} = $read_count;
    }
    unless(defined($params{fragment_count})){
        if($self->is_paired_end){
            $params{fragment_count} = $params{read_count} * 2;
        }
        else {
            $params{fragment_count} = $params{read_count};
        }
    }
    if(not defined($params{subset_name})){
        my $subset_name = $self->get_subset_name;
        unless(defined($subset_name) and $subset_name =~ /[1-8]/){
            $self->error_message("Subset_name must be between 1-8. Found subset_name of \"".$subset_name."\"");
            die $self->error_message;
        }
        $self->subset_name($subset_name);
        $params{subset_name} = $subset_name;
    }
    else {
        my $subset_name = $self->get_subset_name;
        my ($subset_name_to_compare) = $params{subset_name} =~ m/^(\d+)(?:$|[-.])/; #ignore indices or other subsets of lanes
        unless($subset_name eq $subset_name_to_compare){
            $self->error_message("Subset name is incorrectly specified. The specified name must match the first digit after s_ in the filenames.
                                    You specified ".$params{subset_name}." which suggests the files should have " . $subset_name_to_compare . " but they instead have ".$subset_name);
            die $self->error_message;
        }
        $self->status_message("Subset name verified to jibe with file names.");
    }

    my $import_instrument_data = Genome::InstrumentData::Imported->create(%params);
    unless ($import_instrument_data) {
       $self->error_message('Failed to create imported instrument data for '.$self->original_data_path);
       return;
    }

    unless ($import_instrument_data->library_id) {
        Carp::confess("No library on new instrument data?");
    }

    my $instrument_data_id = $import_instrument_data->id;
    $self->status_message("Instrument data record $instrument_data_id has been created.");
    $self->generated_instrument_data_id($instrument_data_id);

    my $ref_name = $self->reference_name;

    my $sources = $self->source_data_files;

    if( $sources =~ s/\/\//\//g) {
        $self->source_data_files($sources);
    }

    my @input_files = split /\,/, $self->source_data_files;
    for (sort(@input_files)) {
        unless( -s $_) {
            $self->error_message("Input file(s) were not found $_");
            die $self->error_message;
        }
    }
    $self->source_data_files(join( ',',sort(@input_files)));

    $self->status_message("About to get a temp allocation");
    my $tmp_tar_file = File::Temp->new("fastq-archive-XXXX",DIR=>"/tmp");
    my $tmp_tar_filename = $tmp_tar_file->filename;

    my $suff = ".txt";
    my $basename;
    my %basenames;
    my @inputs;
    my %filenames;
    for my $file (sort(@input_files)) {
        my ($filename,$path,$suffix) = fileparse($file, $suff);
        $basenames{$path}++;
        $basename = $path;
        my $fastq_name = $filename.$suffix;
        if(not defined($self->is_valid_filename($fastq_name))){
            my $fixed_name = $self->fix_fastq_filename($path,$fastq_name);
            unless($fixed_name){
                $self->error_message("File basename - $fastq_name - did not have the form:
                    \n\t\t\t\t s_[1-8]_sequence.txt or s_[1-8]_[1-2]_sequence.txt\t\tand could not be fixed.\n");
                die $self->error_message;
            }
            $filenames{$fastq_name} = $fixed_name;
        }else{
            $filenames{$fastq_name} = $fastq_name;
        }
    }
    unless(scalar(keys(%basenames))==1) {
        $self->error_message("Found more than one path to imported files.");
        die $self->error_message;
    }
    my $change_names = 0;
    for my $file (sort(keys(%filenames))){
        if( $file ne $filenames{$file}){
            $change_names++;
        }
    }
    my $tar_cmd;
    # if the file names need to be changed, we define different tar comamands
    if(($change_names>0)&& $self->is_paired_end){
        if($change_names==2){
            my ($file1,$file2) = keys(%filenames);
            if ($self->tar_can_multi_transform) {
                $tar_cmd = sprintf("tar cvzfh %s -C %s %s %s --transform=s/%s/%s/g --transform=s/%s/%s/g", $tmp_tar_filename, $basename, $file1, $file2, $file1, $filenames{$file1}, $file2, $filenames{$file2});
            } else {
                $self->error_message("The version of tar that is installed is not able to correct multiple file names. Please change the names to match the proper format, as labeled above.");
                die $self->error_message;
            }
        } elsif ($change_names==1) {
            # if only one needs changing, find which one and apply the single transform in the tar command. This will work with our current 1.19 tar install.
            my ($file1,$file2) = keys(%filenames);
            my ($changed_name,$same_name);
            if($filenames{$file1} eq $file1){
                $changed_name = $file2;
                $same_name = $file1;
            } else {
                $changed_name = $file1;
                $same_name = $file2;
            }
            $tar_cmd = sprintf("tar cvzfh %s -C %s %s %s --transform=s/%s/%s/g",$tmp_tar_filename,$basename, $file1,$file2,$changed_name,$filenames{$changed_name});
        } else {
            $self->error_message("Found more than 2 file names to change, which should not be possible, since there should never be more than 2 files.");
            die $self->error_message;
        }
    } else {
        if($change_names == 1){
            my ($file1) = keys(%filenames);
            $tar_cmd = sprintf("tar cvzfh %s -C %s %s --transform=s/%s/%s/g",$tmp_tar_filename,$basename, $file1,$file1,$filenames{$file1});
        } else {
            $tar_cmd = sprintf("tar cvzfh %s -C %s %s",$tmp_tar_filename,$basename, join " ", keys(%filenames));
        }
    }
    $self->status_message("About to execute tar command, this could take a long time, depending upon the location (across the network?) and size (MB or GB?) of your fastq's.");
    unless(Genome::Sys->shellcmd(cmd=>$tar_cmd)){
        $self->error_message("Tar command failed to complete successfully. The command looked like :   ".$tar_cmd);
        die $self->error_message;
    }

    $import_instrument_data->original_data_path($self->source_data_files);

    my $kb_usage = $import_instrument_data->calculate_alignment_estimated_kb_usage;

    unless ($kb_usage) {
        $self->warning_message('Failed to get estimate kb usage for instrument data '.$instrument_data_id);
        return 1;
    }

    my $alloc_path = sprintf('instrument_data/imported/%s', $instrument_data_id);


    my %alloc_params = (
        disk_group_name     => $ENV{GENOME_DISK_GROUP_ALIGNMENTS},
        allocation_path     => $alloc_path,
        kilobytes_requested => $kb_usage,
        owner_class_name    => $import_instrument_data->class,
        owner_id            => $import_instrument_data->id,
    );


    my $disk_alloc;


    if($self->allocation) {
        $disk_alloc = $self->allocation;
    } else {
        $disk_alloc = Genome::Disk::Allocation->allocate(%alloc_params);
    }
    unless ($disk_alloc) {
        $self->error_message("Failed to get disk allocation with params:\n". Data::Dumper::Dumper(%alloc_params));
        return 1;
    }
    $self->status_message("Disk allocation created for $instrument_data_id ." . $disk_alloc->absolute_path);

    my $real_filename = sprintf("%s/archive.tgz", $disk_alloc->absolute_path);
    $self->status_message("About to calculate the md5sum of the tar'd fastq's. This may take a long time.");
    my $md5 = Genome::Sys->md5sum($tmp_tar_filename);
    $self->status_message("Copying tar'd fastq's into the allocation, this will take some time.");
    unless(copy($tmp_tar_filename, $real_filename)) {
        $self->error_message("Failed to copy to allocated space (copy returned bad value).  Unlinking and deallocating.");
        unlink($real_filename);
        $disk_alloc->deallocate;
        return;
    }
    $self->status_message("About to calculate the md5sum of the tar'd fastq's in their new habitat on the allocation. This may take a long time.");
    my $copy_md5;
    unless($copy_md5 = Genome::Sys->md5sum($real_filename)){
        $self->error_message("Failed to calculate md5sum.");
        die $self->error_message;
    }
    unless($copy_md5 eq $md5) {
        $self->error_message("Failed to copy to allocated space (md5 mismatch).  Unlinking and deallocating.");
        unlink($real_filename);
        $disk_alloc->deallocate;
        return;
    }
    $self->status_message("The md5sum for the copied tar file is: ".$copy_md5);
    eval { $disk_alloc->reallocate; }; #don't want to fail just for this
    if($@) { $self->warning_message($@); }
    $self->status_message("The instrument-data id of your new record is ".$instrument_data_id);
    return 1;

}

sub tar_can_multi_transform {
    my ($tar_version) = qx(tar --version 2>&1 | head -n 1) =~ /([.\d]+)/;
    my $tar_vstring = eval "v$tar_version";
    return ($tar_vstring ge v1.21);
}

sub _check_fastq_integrity {
    my $self = shift;
    my $answer = 0;
    my @filepaths = split ",",$self->source_data_files;
    for my $path (@filepaths){
        unless(-s $path){
            $self->error_message("The file at ".$path." was not found.");
            die $self->error_message;
        }
        $self->_check_quality_scores($path);
    }
    if(@filepaths==1){
        if(defined($self->is_paired_end)){
            unless(not $self->is_paired_end){
                $self->error_message("The parameters specify paired-end data, but only one fastq was supplied. If the data are paired-end but only one set of reads is provided, please specify non-paired end.");
                die $self->error_message;
            }
        }
        unless($self->check_last_read($filepaths[0])){
            $self->error_message("Could not validate the last read (last 4 lines) of the imported fastq file.");
            die $self->error_message;
        }

    } elsif (@filepaths==2) {
        if(defined($self->is_paired_end)){
            unless($self->is_paired_end){
                $self->error_message("InstrumentData parameters specify non-paired-end data, but two fastq's were found.");
                die $self->error_message;
            }
            my ($forward_read_length,$reverse_read_length);
            my ($forward_read_name,$reverse_read_name);
            unless(($forward_read_name,$forward_read_length) = $self->check_last_read($filepaths[0])){
                $self->error_message("Could not validate the last read (last 4 lines) of the imported fastq file.");
                die $self->error_message;
            }
            unless(($reverse_read_name,$reverse_read_length) = $self->check_last_read($filepaths[1])){
                $self->error_message("Could not validate the last read (last 4 lines) of the imported fastq file.");
                die $self->error_message;
            }
            unless($forward_read_name eq $reverse_read_name){
                $self->error_message("Forward and Reverse read names do not match.");
                die $self->error_message;
            }
        } else {
            $self->error_message("Found two fastq files but is_paired_end parameter was not set.");
            die $self->error_message;
        }
    } else {
        $self->error_message("Found ".scalar(@filepaths)." fastq files. This tool only allows one or two (paired-end) at a time.");
        die $self->error_message;
    }
}

sub _check_quality_scores {
    my ($self, $filename) = @_;

    my $lines_to_validate = 200;

    $self->status_message(sprintf(
            "Validating sample quality scores from first $lines_to_validate lines of $filename for import format %s.",
            $self->import_format));

    # grab some sample data from the file
    my $head = `head -$lines_to_validate $filename`;
    my @sample_lines = split("\n", $head);

    # We just want every 4th line
    my @quality_score_lines = @sample_lines[
        grep{0 == ($_ + 1) % 4} 0..$#sample_lines];

    for my $qs_line (@quality_score_lines) {
        unless ($self->_validate_quality_scores($qs_line)) {
            $self->error_message(sprintf(
                    "Couldn't validate quality scores for first $lines_to_validate lines of $filename as %s format.",
                $self->import_format));
            die $self->error_message;
        }
    }
}

sub _validate_quality_scores {
    my ($self, $line) = @_;

    my %allowed_qs_chars = (
        'sanger fastq' => '[!-~]*',
        'solexa fastq' => '[;-~]*',
        'illumina fastq' => '[@-~]*',
    );

    my $filter_chars = $allowed_qs_chars{$self->import_format};
    $line =~ s/$filter_chars//g;
    chomp $line;

    # there should be nothing left after removing valid quality scores
    if (length($line)) {
        $line =~ s/(.)(?=.*?\1)//g; # find unique characters
        $self->error_message("Invalid characters in quality score: $line");
        return;
    }
    return 1;
}

sub check_last_read {
    my $self = shift;
    my $fastq = shift;
    my $tail = `tail -4 $fastq`;
    my @lines = split "\n",$tail;
    $self->status_message("Checking last read of FastQ file ($fastq):\n$tail");
    unless(@lines==4){
        $self->error_message("Didn't get 4 lines from the fastq.");
        return;
    }
    unless($lines[0] =~ /^@/){
        return;
    }
    my ($read_name) = split m![/ ]!,$lines[0];
    my $read_length = length $lines[1];
    if(defined $self->read_length){
        unless($read_length <= $self->read_length){
            $self->error_message("Read-Length was set to ".$self->read_length." however, a read of length ".$read_length." was found in the last read.");
            return;
        }
    }
    unless((length $lines[2]) >= 1){
        $self->error_message("Quality Score name was too short.");
        return;
    }
    unless($read_length == (length $lines[3])){
        $self->error_message("Length of read did not match length of quality string.");
        die $self->error_message;
    }
    return ($read_name,$read_length);
}

sub get_read_count_and_fragment_count {
    my $self = shift;
    my ($line_count,$read_count,$fragment_count);
    my @files = split ",", $self->source_data_files;
    $self->status_message("Now attempting to determine read_count by calling wc on the imported fastq(s). This may take a while if the fastqs are large.");
    my $last_count;  #use this for checking if paired_end files have the same amt of reads
    for my $file (@files){
        my $sub_count = `wc -l $file`;
        ($sub_count) = split " ",$sub_count;
        unless(defined($sub_count)&&($sub_count > 0)){
            $self->error_message("couldn't get a response from wc.");
            return;
        }
        if ($last_count){
            unless ($last_count = $sub_count){
                die $self->error_message("Two provided fastq files have differing line counts and number of reads! (lines: $last_count and $sub_count)");
            }
        }else{
            $last_count=$sub_count;
        }
        $line_count += $sub_count;
    }
    if($line_count % 4){
        $self->error_message("Calculated a number of lines in the fastq file that was not divisible by 4.");
        return;
    }
    $read_count = $line_count / 4;
    if(scalar(@files)==2){
        $read_count = $read_count / 2;
        $fragment_count = $read_count * 2;
    }
    else {
        $fragment_count = $read_count;
    }
    return $read_count,$fragment_count;
}

sub get_subset_name {
    my $self = shift;

    my $subset_name;
    my @files = split(/,/, $self->source_data_files);
    for my $file (@files) {
        my ($filename, $path, $suffix) = fileparse($file, '.txt');
        my ($sname) = $filename =~ /s_(\d)(_\d)?_sequence/;
        if(not defined $subset_name) {
            $subset_name = $sname;
        } elsif ($sname ne $subset_name) {
            die "The subset names didn't match.";
        }
    }

    return $subset_name;
}

sub is_valid_filename {
    my $self = shift;
    my $fastq_name = shift;
    unless(($fastq_name=~m/^s_[1-8]_sequence.txt$/)||($fastq_name=~m/^s_[1-8]_[1-2]_sequence.txt$/)){
        return;
    }
    return 1;
}

sub fix_fastq_filename {
    my $self = shift;
    my ($path,$file) = @_;
    #my $file = shift;
    print "file = ".$file."\n";
    print "path = ".$path."\n";
    my $new_filename;
    unless(defined($self->subset_name)){
        $self->error_message("Attempting to fix fastq filenames, but no subset_name was specified. This is required if the fastq filenames do not contain it.");
        die $self->error_message;
    }
    if($self->is_paired_end){
        my $fh = new IO::File;
        $fh->open($path.$file);
        my $line = $fh->getline;
        $fh->close;
        chomp $line;
        my ($x,$direction) = split "/", $line;
        print "line = ".$line."  and direction =  ".$direction."\n";
        $new_filename = "s_".$self->subset_name."_".$direction."_sequence.txt";
        print "new_filename = ".$new_filename."\n";
    } else {
        $new_filename = "s_".$self->subset_name."_sequence.txt";
    }
    if($self->is_valid_filename($new_filename)){
        return $new_filename;
    } else {
        return;
    }
}

1;
