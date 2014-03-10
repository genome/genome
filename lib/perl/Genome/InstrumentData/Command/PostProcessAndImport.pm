package Genome::InstrumentData::Command::PostProcessAndImport;

use strict;
use warnings;
use Genome;
use File::Basename;

our $UNALIGNED_TEMPDIR = '/tmp'; #'/gscmnt/sata844/info/hmp-mgs-test-temp';

class Genome::InstrumentData::Command::PostProcessAndImport{
    is => "Genome::Command::Base",
    has => {
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'Instrument Data to dust, n-remove and import as new instrument data',
        },
        n_removal_threshold => {
            is => 'Int',
            doc => 'reads with this amount of ns or more will be removed in post-processing, defaults to 0, which indicates no reads with ns will be removed',
            default => 0,
        },
        non_n_base_threshold => {
            is => 'Int',
            doc => 'reads with less than this amount of non-n bases will be removed in post-processing',
            default => 0,
        },
        dust => {
            is => 'Boolean',
            doc => 'dust before n removal if flag is specified',
            default => 0,
        },
        post_processed_instrument_data => {
            is => 'Genome::InstrumentData::Imported',
            doc => 'this parameter holds the post processed, imported instrument data',
            is_output => 1,
            is_many => 1,
            is_optional => 1,
        },
    }
};


sub execute {
    my $self = shift;
    my $instrument_data = $self->instrument_data;
    my $lane = $instrument_data->lane;
    my $instrument_data_id = $instrument_data->id;

    my $tmp_dir = "$UNALIGNED_TEMPDIR/unaligned_reads";
    unless ( -d $tmp_dir or mkdir $tmp_dir) {
        die "Failed to create temp directory $tmp_dir : $!";
    }

    $tmp_dir .= "/$instrument_data_id";

    unless (-d $tmp_dir) {
        mkdir $tmp_dir or die "Failed to create temp directory $tmp_dir : $!";
    }

    # TODO: dust, n-remove and set the sub-dir based on the formula
    # and use a subdir name built from that formula
    my $subdir=$tmp_dir;
    if ($self->dust){
        $subdir.='/dusted';
    }
    unless (-d $subdir or mkdir $subdir) {
        die "Failed to create temp directory $subdir : $!";
    }
    if ($self->n_removal_threshold){
        $subdir .= '/n-remove_'.$self->n_removal_threshold;
    }elsif($self->non_n_base_threshold){
        $subdir .= '/non-n-base-filter_'.$self->non_n_base_threshold;
    }

    unless (-d $subdir or mkdir $subdir) {
        die "Failed to create temp directory $subdir : $!";
    }

   

    $self->status_message("Preparing imported instrument data for import path $subdir");

    # proceed extracting and uploading unaligned reads into $subdir....

    # resolve the paths at which we will place processed instrument data
    # we're currently using these paths to find previous unaligned reads processed the same way

    my $forward_basename = "s_${lane}_1_sequence.txt";
    my $reverse_basename = "s_${lane}_2_sequence.txt";
    my $n_removed_fragment_basename = "s_${lane}_sequence.txt";
    my $fragment_basename = "s_${lane}_sequence.txt";

    my $expected_path;
    my $expected_path1; #for paired end fastq processing
    my $expected_path2;
    my $expected_path_fragment; #for paired end w/ n-removal, potentially

    if ($instrument_data->is_paired_end){
        $expected_path1 = "$subdir/$forward_basename";
        $expected_path2 = "$subdir/$reverse_basename";
        $expected_path_fragment = "$subdir/$n_removed_fragment_basename";
        $expected_path = $expected_path1 . ',' . $expected_path2;
    }else{
        $expected_path = "$subdir/$fragment_basename";
    }

    $self->status_message("Checking for previously post-processed and reimported reads from: $expected_path");
    my $post_processed_inst_data = Genome::InstrumentData::Imported->get(original_data_path => $expected_path);
    if ($post_processed_inst_data) {
        $self->status_message("post processed instrument data already found for path $expected_path");
        $self->status_message("skipping read processing since all data is already processed and uploaded");

        my @return = ($post_processed_inst_data);
        #check if we produced a fragment from pairwise n-removal to return shortcut

        if ($expected_path_fragment and ($self->n_removal_threshold or $self->non_n_base_threshold) ){
            my $n_removed_fragment = Genome::InstrumentData::Imported->get(original_data_path => $expected_path_fragment);
            push @return, $n_removed_fragment if $n_removed_fragment;
        }

        $self->post_processed_instrument_data(\@return);
        return @return;
    }

    my @instrument_data = eval {

	my $import_lock;

        # check for previous unaligned reads
        $import_lock = $self->lock($instrument_data_id, basename($expected_path));
        unless ($import_lock) {
            die $self->error_message("Failed to lock $expected_path.");
        }
        my $upload_path = $expected_path;

        # extract
        $self->status_message("Preparing imported instrument data for import path $expected_path");

        my $fastq_filenames = $instrument_data->resolve_fastq_filenames;
        for (@$fastq_filenames){
            unless (-s $_){
                $self->error_message("expected fastq ($_) extracted from instrument data ".$instrument_data->__display_name__." doesn't have size!");
                system("ls -l $_");
                die $self->error_message;
            }
        }
        if ($instrument_data->is_paired_end){
            unless (@$fastq_filenames == 2){
                $self->error_message("instrument_data ".$instrument_data->display_name." is paired end but doesn't have 2 fastq files!");
                die $self->error_message;    
            }
            my ($forward) = grep {$_ =~ $forward_basename} @$fastq_filenames;
            my ($reverse) = grep {$_ =~ $reverse_basename} @$fastq_filenames;
            unless($forward and -s $forward and $reverse and -s $reverse){
                $self->error_message("couldn't find expected fastq basenames in ".$instrument_data->display_name);
                die $self->error_message;
            }
            my ($processed_fastq1,$processed_fastq2, $processed_fastq_fragment) = $self->_process_unaligned_fastq_pair($forward,$reverse,$expected_path1, $expected_path2, $expected_path_fragment );
            my @missing = grep {! -s $_} ($expected_path1, $expected_path2);
            if (@missing){
                $self->error_message("Expected data paths do not exist after fastq processing: ".join(", ", @missing));
                die($self->error_message);
            }
        }else{
            unless (@$fastq_filenames == 1){
                $self->error_message("instrument_data ".$instrument_data->display_name." is not paired end but doesn't have exactly 1 fastq file!"); 
                die $self->error_message;
            }
            my ($fragment) = grep {$_ =~ $fragment_basename} @$fastq_filenames;
            unless ($fragment and -e $fragment){
                $self->error_message("couldn't find expected fastq basename in ".$instrument_data->display_name);
                die $self->error_message;
            }

            my $processed_fastq = $self->_process_unaligned_fastq($fragment, $expected_path);
            unless (-s $expected_path){
                die $self->error_message("Expected data path does not exist after fastq processing: $expected_path");
            }
        }

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
            my $value;
            if ($property_name eq 'subset_name'){
                $value = $instrument_data->lane;  #subset name can include import source, we just want lane
            }else{
                $value = $instrument_data->$property_name;
            }
            no warnings;
            $self->debug_message("Value for $property_name is $value");
            $properties_from_prior{$property_name} = $value;
        }
        $properties_from_prior{subset_name} = $instrument_data->lane;

        if ($upload_path =~ /,/){  #technically this can go in the @properties_prior_array above, but i'm trying to keep as much in common with processed_unaligned_reads as possible to simplify refactoring
            $properties_from_prior{is_paired_end} = 1;
        }else{
            $properties_from_prior{is_paired_end} = 0;
        }

        my %params = (
            %properties_from_prior,
            source_data_files => $upload_path,
            import_format => 'sanger fastq', #this is here because of sra data, happens to coincide with skip_contamination_screen, but don't think that will always be the case, probably should improve how this is derived
        );
        $self->status_message("importing fastq with the following params:" . Data::Dumper::Dumper(\%params));

        my $command = Genome::InstrumentData::Command::Import::Fastq->create(%params);
        unless ($command) {
            $self->error_message( "Couldn't create command to import unaligned fastq instrument data!");
        };
        my $result = $command->execute();
        unless ($result) {
            die $self->error_message( "Error importing data from $upload_path! " . Genome::InstrumentData::Command::Import::Fastq->error_message() );
        } 

        #now create fragment instrument data for $expected_path_fragment, as a result of pairwise n-removal, if the file exists and has size

        if ($expected_path_fragment and -s $expected_path_fragment){
            $params{source_data_files} = $expected_path_fragment;
            $params{is_paired_end} = 0;
            $self->debug_message("importing fastq with the following params:" . Data::Dumper::Dumper(\%params));

            my $command = Genome::InstrumentData::Command::Import::Fastq->create(%params);
            unless ($command) {
                $self->error_message( "Couldn't create command to import unaligned fastq instrument data!");
            };
            my $result = $command->execute();
            unless ($result) {
                die $self->error_message( "Error importing data from $upload_path! " . Genome::InstrumentData::Command::Import::Fastq->error_message() );
            } 
        }

        $self->status_message("committing newly created imported instrument data");
        $self->debug_message("UR_DBI_NO_COMMIT: ".$ENV{UR_DBI_NO_COMMIT});
        UR::Context->commit(); # warning: most code should NEVER do this in a pipeline

        my @return_inst_data;

        my $new_instrument_data = UR::Context->current->reload("Genome::InstrumentData::Imported",
            original_data_path => $upload_path
        );
        unless ($new_instrument_data) {
            die $self->error_message( "Failed to find new instrument data $upload_path!");
        }
        if ($new_instrument_data->__changes__) {
            die "Unsaved changes present on instrument data $new_instrument_data->{id} from $upload_path!!!";
        }
        push @return_inst_data, $new_instrument_data;

        if ($expected_path_fragment and -s $expected_path_fragment){
            my $new_instrument_data = UR::Context->current->reload("Genome::InstrumentData::Imported",
                original_data_path => $expected_path_fragment,
            );
            unless ($new_instrument_data) {
                die $self->error_message( "Failed to find new instrument data $upload_path!");
            }
            if ($new_instrument_data->__changes__) {
                die "Unsaved changes present on instrument data $new_instrument_data->{id} from $upload_path!!!";
            }
            push @return_inst_data, $new_instrument_data;
        }

        if ($import_lock) {
            unless(Genome::Sys->unlock_resource(resource_lock => $import_lock)) {
                die $self->error_message("Failed to unlock $expected_path.");
            }
        }
        return @return_inst_data;
    };

    # TODO: add directory removal to Genome::Sys
    if ($@) {
        system "/bin/rm -rf '$tmp_dir'";
        die $self->error_message("Error processing unaligned reads! $@");
    }
    
    system "/bin/rm -rf '$tmp_dir'";
    unless (@instrument_data){
        die $self->error_message("no resulting post_processed instrument data!");
    }

    $self->post_processed_instrument_data(\@instrument_data);
    return @instrument_data;
}

sub lock {
    my $self = shift;
    my @parts = @_;
    my $lock_key = join('_', @parts);
    my $resource_lock = File::Spec->join($ENV{GENOME_LOCK_DIR}, $lock_key);
    $self->debug_message("Creating lock on $lock_key...");
    my $lock = Genome::Sys->lock_resource(
        resource_lock => $resource_lock,
        max_try => 2,
    );
    return $lock;
};

sub _process_unaligned_fastq_pair {
    my $self = shift;
    my ($forward, $reverse, $forward_out, $reverse_out, $fragment_out) = @_;
    #run dust on forward and reverse
    my $forward_dusted;
    my $reverse_dusted;

    if ($self->dust){
        $self->debug_message("Dusting fastq pair $forward, $reverse");
        $forward_dusted = "$forward.DUSTED";
        $reverse_dusted = "$reverse.DUSTED";

        $self->dust_fastq($forward, $forward_dusted);
        $self->dust_fastq($reverse, $reverse_dusted);
    }else{
        $self->debug_message("skipping dusting");
        $forward_dusted = $forward;
        $reverse_dusted = $reverse;
    }

    #run pairwise n-removal
    if ($self->n_removal_threshold or $self->non_n_base_threshold){
        $self->debug_message("running remove-n-pairwise on $forward, $reverse");
	my %params= (forward_fastq => $forward_dusted,
            reverse_fastq => $reverse_dusted,
            forward_n_removed_file => $forward_out,
            reverse_n_removed_file => $reverse_out,
            singleton_n_removed_file => $fragment_out
	    );
	if ($self->n_removal_threshold){
	    $params{n_removal_threshold}=$self->n_removal_threshold;
	}elsif ($self->non_n_base_threshold){
	    $params{non_n_base_threshold}=$self->non_n_base_threshold;
	}
        my $cmd = Genome::Model::Tools::Fastq::RemoveNPairwise->create(
            %params
	    );
	
        unless ($cmd){
            die $self->error_message("couldn't create remove-n-pairwise command for $forward_dusted, $reverse_dusted!");
        }
        my $rv = $cmd->execute;
        unless ($rv){
            die $self->error_message("couldn't create remove-n-pairwise command for $forward_dusted, $reverse_dusted!");
        }
        unless(-e $forward_out && -e $reverse_out && -e $fragment_out){
            die $self->error_message("couldn't find all expected output files! $forward_out, $reverse_out, $fragment_out");
        }
        #clean up, maybe make these temp files
        if ($self->dust){
            #only need to do this if we actually dusted
            unlink $forward_dusted;
            unlink $reverse_dusted;
        }

        #return the 3 processed fastq files
        return ($forward_out, $reverse_out, $fragment_out);
    }else{
        $self->status_message("skipping n-removal");
        Genome::Sys->copy_file($forward_dusted, $forward_out);
        Genome::Sys->copy_file($reverse_dusted, $reverse_out);
        if ($self->dust){
            #only need to do this if we actually dusted
            unlink $forward_dusted;
            unlink $reverse_dusted;
        }
        return ($forward_out, $reverse_out);
    }
}

sub _process_unaligned_fastq {
    my $self = shift;
    my ($fastq_file, $output_path) = @_;

    my $dusted_fastq;
    if ($self->dust){
        $dusted_fastq = "$fastq_file.DUSTED";
        $self->dust_fastq($fastq_file, $dusted_fastq);
    }else{
        $self->status_message("skipping dusting $fastq_file");
        $dusted_fastq = $fastq_file;
    }

    if ($self->n_removal_threshold or $self->non_n_base_threshold){
        $self->status_message("Running n-removal on file $fastq_file");
	my %params=( fastq_file => $dusted_fastq,
		     n_removed_file => $output_path
	    ); 
	if ($self->n_removal_threshold){
	    $params{n_removal_threshold}=$self->n_removal_threshold;
	}elsif ($self->non_n_base_threshold){
	    $params{non_n_base_threshold}=$self->non_n_base_threshold;
	}
        my $cmd = Genome::Model::Tools::Fastq::RemoveN->create(%params);
        unless ($cmd){
            die $self->error_message("couldn't create remove-n command for $dusted_fastq");
        }
        my $rv = $cmd->execute;
        unless ($rv){
            die $self->error_message("couldn't execute remove-n command for $dusted_fastq");
        }
    } else {
        $self->status_message("No n-removal params specified, skipping");
        Genome::Sys->copy_file($dusted_fastq, $output_path);
    }
    if ($self->dust){
        unlink $dusted_fastq;
    }
    return $output_path;
}


sub dust_fastq{
    my ($self, $in, $out) = @_;
    my $cmd = Genome::Model::Tools::Fastq::Dust->create(
        fastq_file => $in,
        output_file => $out,
    );
    unless ($cmd){
        die $self->error_message("couldn't create dust command for $in -> $out!");
    }
    my $rv = $cmd->execute;
    unless ($rv){
        die $self->error_message("failed to execute dust command for $in -> $out! rv:$rv");
    }
    unless (-s $out){
        die $self->error_message("expected output file $out doesn't exist or has 0 size!");
    }
    return $out;
}
