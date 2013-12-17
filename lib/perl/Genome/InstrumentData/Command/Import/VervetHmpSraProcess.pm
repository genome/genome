package VervetHmpSraProcess;
#package Genome::InstrumentData::Command::Hmp:: ;

use strict;
use warnings;
use Genome;
use Cwd;
use IO::File;
use File::Basename;
use File::Path;


#### ADD:
# 1) rename files to conform to naming standards (eg. PE: s_1_1_sequence.txt, s_1_2_sequence.txt    FRAGMENT: s_1_sequence.txt) when putting data into allocation spot for final instrument-data id inside the archive (archive.tgz)
# 2) check to see if instrument data ids already exist (less important, but should be done)




#class Genome::InstrumentData::Command::Hmp::
class VervetHmpSraProcess {
    is  => 'Command::V2',
    has => [
	srs_sample_id => { is_optional => 0, doc => 'SRS sample id to extract ... for Vervet this is the sample_name, such as WFAA-1088-0104014095', },
	srr_accessions => { is_optional => 0, doc => 'space separated list of SRR accession ids for the raw SRA data downloads to use.', },
	picard_dir => { is_optional => 1, default_value => "$ENV{GENOME_SW_LEGACY_JAVA}/samtools/picard-tools-1.27", doc => 'full path to directory containing Picard jar files (note: This path must include the updated EstimateLibraryComplexity that handles redundancy removal)', },
	_working_dir => { is_optional=>1, is_transient=>1, },
	organism => { is_optional => 1, default_value => "Chlorocebus aethiops", }
    ],
};

sub execute {
    my $self = shift;

#___This line stops the perl debugger as though I'd set a break point in the GUI
    $DB::single = 1;
    
    my @srrs = split /\s+/, $self->srr_accessions;


#___grab & setup fastq import parameters ... NOTE: These are arguments that are setup based on the SAMPLE, so grabbing this stuff from ANY single one of the incoming SRR ids (or in the case of Vervet, the 'subset_name' which is in the 'description' field).   This module gets run once per SAMPLE (which usually has 2 'subset_names')
    my %import_params;
    my (@inst_data) = Genome::InstrumentData::Imported->get(
	import_format => 'bam',
	sample_name   => $self->srs_sample_id
	);

    $import_params{library_name} = $inst_data[0]->library_name;
    $import_params{sequencing_platform} = 'solexa';
    $import_params{import_format} = 'sanger fastq';
    $import_params{sample_name} = $self->srs_sample_id; #For Vervet, I'm using the 'sample_name' input argument to hold the sample name (such as WFAA-1253-0104014103) ... this field used to hold the EMMES id (HMPZ-0000000000 or etc...) ... the SRS id was in a different field (sra_sample_id)
    $import_params{species_name} = $self->organism, #This value defaults to Vervet
  
    my $tmp_dir      = Genome::Sys->create_temp_directory;

    my $working_dir = $tmp_dir . "/srr_datasets";
    $self->_working_dir($working_dir);
    mkpath($working_dir);
    for my $srr (@srrs) {
        ####my $instrument_data = Genome::InstrumentData::Imported->get(import_format=>'raw sra download', sra_accession=>$srr);
	my (@instrument_data) = Genome::InstrumentData::Imported->get(import_format => 'bam', sample_name   => $self->srs_sample_id);

#_______This section here deals with the fact that I am NOT able to nicely grab an instrument-data id by an SRR id, since I don't have SRR ids...I only have subset_names in the 'description' field
	foreach my $single_inst_data (@instrument_data) {

#___________First I need to make sure I'm ONLY working with the current 'subset_name'
	    my $description = $single_inst_data->description;
	    my ($subset_name) = $description =~ /\S+\s+\(subset_name\:(\d+)\s+\S+/;
	    next unless ($subset_name eq $srr); #This should ensure that I only do the allocation once per 'subset_name'

	    my ($alloc) = $single_inst_data->allocations; 
	    unless(symlink($alloc->absolute_path . "/" . $srr, $working_dir . "/" . $srr)) {
		$self->error_message("Failed to set up symlink from SRA data dir: " . $alloc->absolute_path . " to " . $working_dir . "/" . $srr); #NOTE: In this case, for Vervet, the $srr is a 'subset_name' (since thats what was entered as 'srr_accessions)
		return;
	    }
	}

    } 



#___Prepare for vervet_process_runs.sh run & for finding output for later fastq importing
    my ($fwd_read, $rev_read, $singleton_read);


    $DB::single = 1;
    
    my $path_to_scripts_dir =$self->get_script_path;
    $self->status_message("Scripts are in: $path_to_scripts_dir");
    
    #Find current path to the script 'trimBWAstyle.usingBam.pl'
    my $current_dir = `pwd`;
    
    my $sra_samples = $working_dir . "/sample_mapping.txt";
    my $list_of_srrs = $working_dir . "/srr_listing.txt";
    
    unless (open (SRA_SAMPLE_MAPPING, ">$sra_samples")) {
	$self->error_message("Failed to open SRR/SRS data mapping file");
	return;
    }
    
    unless (open (SRR_LISTING, ">$list_of_srrs")) {
	$self->error_message("Failed to open SRR/SRS data mapping file");
	return;
    }
    
    for (@srrs) {
	print SRA_SAMPLE_MAPPING sprintf("%s\t%s\n", $_, $self->srs_sample_id);
	print SRR_LISTING sprintf("%s\n", $_);
	
    }
    close SRA_SAMPLE_MAPPING;
    close SRR_LISTING;
    
    $DB::single = 1;
    
    #Run BROAD's processing script
    my $cmd;
    my $errfile = $working_dir . "/ReadProcessing." . $self->srs_sample_id . ".err";
    my $outfile = $working_dir . "/ReadProcessing." . $self->srs_sample_id  . ".out";
    my $picard_dir   = $self->picard_dir;
    
    #Set the $PATH env variable in perl
    my $path = $ENV{'PATH'} . ":" . $path_to_scripts_dir;
    
    
    
    
    #Note: I need to set the path to my scripts INSIDE the shell command
    ####$cmd = "cd $working_dir; export PATH=$path; process_runs.sh $list_of_srrs $sra_samples $picard_dir $tmp_dir > $outfile 2> $errfile";
    
#_______THIS LINE IS FOR REAL DATA...UNCOMMENT TO REALLY DO PROCESSING AND UPLOAD TO DACC!!!! jmartin ... 110927
    ####$cmd = "cd $working_dir; export PATH=$path; process_runs.sh $list_of_srrs $sra_samples $picard_dir $tmp_dir $user $pwd";
    
#_______This line is just for testing, the DACC upload is commented out in this version jmartin ... 110923
    ####$cmd = "cd $working_dir; export PATH=$path; process_runs.for_testing.sh $list_of_srrs $sra_samples $picard_dir $tmp_dir $user $pwd";
    
#_______This line is for using 'vervet_process_runs.sh' ... jmartin 111031
    ####$cmd = "cd $working_dir; export PATH=$path; vervet_process_runs.sh $list_of_srrs $sra_samples $picard_dir $tmp_dir";

#_______This version now points to the modified processing engine that does a Picard 'SortSam.jar SORT_ORDER=queryname' before running EstimateLibraryComplexity ... jmartin 111118
    $cmd = "cd $working_dir; export PATH=$path; vervet_process_runs.w_SortSam.sh $list_of_srrs $sra_samples $picard_dir $tmp_dir";


    
    
# BEFORE:
    # $alloc->reallocate(kilobytes_requested => $alloc->kilobytes_requested * 2.2, allow_reallocate_with_move => 1,)
# AFTER:
     # $alloc->reallocate(); # allocates to current size
    
    $self->status_message("CMD=>$cmd<=\n");
    $self->status_message("PWD=>$current_dir<=\n");
    
    $DB::single = 1;

    eval {
        Genome::Sys->shellcmd(
            cmd => $cmd,
            );
    };
    
    if ($@) {
	$DB::single = 1; 
	$self->error_message("Error running broad script: $@\n");
	return;
    }
    
    $DB::single = 1;
    my @reads = glob($working_dir . "/" . $self->srs_sample_id . "/*.trimmed.*.fastq.bz2");
    
    for (@reads) {
	Genome::Sys->shellcmd(cmd=>"bunzip2 $_");
    }
    
    
    ($fwd_read) = glob($working_dir . "/" . $self->srs_sample_id . "/*.trimmed.1.fastq");
    ($rev_read) = glob($working_dir . "/" . $self->srs_sample_id . "/*.trimmed.2.fastq");
    ($singleton_read) = glob($working_dir . "/" . $self->srs_sample_id . "/*.trimmed.singleton.fastq");
    
#________Reinstated to deal with name conformity issues 120113 ... jmartin
    rename($fwd_read, $working_dir . "/s_1_1_sequence.txt");
    $fwd_read = $working_dir . "/s_1_1_sequence.txt";
    
    rename($rev_read, $working_dir . "/s_1_2_sequence.txt");
    $rev_read = $working_dir . "/s_1_2_sequence.txt";
    
    rename($singleton_read, $working_dir . "/s_1_sequence.txt");
    $singleton_read = $working_dir . "/s_1_sequence.txt";
    
    
    $DB::single = 1;


#___Now that the vervet_process_runs.sh is done, import the redundancy-free & quality trimmed fastq data
    my $subset_name_list = join('_',@srrs);

    print "subset_name_list=>$subset_name_list<= ... fwd_read=>$fwd_read<= ... rev_read=>$rev_read<=\n";



#___Grab $library object for use below
    my $library = Genome::Library->get(
	name => $inst_data[0]->library_name
	);



    # make new inst data
    my $pe_original_data_path = $working_dir . "/" . $self->srs_sample_id;
    my $pe_inst_data = Genome::InstrumentData::Imported->create(
	library             => $library,
	import_format       => 'sanger fastq',
	description         => "qualityTrimmed nonRedundant pairedFastq (subset_namesInThisSample:$subset_name_list)", #Made this change to support having variable extrapolated into this object (using single-quotes was just having it store the string '$subset_name_list') ... jmartin 120214
	sequencing_platform => 'solexa',
	is_paired_end       => 1,
	original_data_path  => $pe_original_data_path,
	srs_sample_id         => $self->srs_sample_id
	);

    # make an allocation
    my $fwd_read_size = `du -sk $fwd_read`;
    chomp($fwd_read_size);
    my @fwd_read_size = split(/\s+/,$fwd_read_size);
    my $rev_read_size = `du -sk $rev_read`;
    chomp($rev_read_size);
    my @rev_read_size = split(/\s+/,$rev_read_size);
    my $pe_read_size_in_kb = $fwd_read_size[0] + $rev_read_size[0];
    my $pe_allocation_request_in_kb = $pe_read_size_in_kb * 1.05;
    my $pe_alloc = Genome::Disk::Allocation->create(
	disk_group_name     => $ENV{GENOME_DISK_GROUP_ALIGNMENTS},
	allocation_path     => 'instrument_data/imported/'.$pe_inst_data->id,
	kilobytes_requested => $pe_allocation_request_in_kb,
	owner_class_name    => $pe_inst_data->class,
	owner_id            => $pe_inst_data->id,
    );





#___Changes made for PE ... 120113 ... jmartin
    # tar non zipped fastqs to temp file
    ####my $tmp_tar_file = $pe_alloc->absolute_path.'/TEMP.archive.tgz';
    ####my $pe_import_cmd = "tar cvzfh $tmp_tar_file -C ".$pe_alloc->absolute_path." $fwd_read $rev_read $singleton_read";
    $current_dir = `pwd`;
    chomp($current_dir);
    chdir($pe_alloc->absolute_path);
    print "manual chdir into =>$pe_alloc->absolute_path<=\n";




###################NOTE: Changed the archive & copy block to use symbolic link ... jmartin 120123
    # link in the forward & reverse read files for the current pair for archiving
    my $pe_fwd_link_cmd = "ln -s $fwd_read ./s_1_1_sequence.txt";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $pe_fwd_link_cmd); };
    if (not $rv) {
	unlink("s_1_1_sequence.txt");
	die "ERROR! failed to link in fwd read (s_1_1_sequence.txt) of read pair\n";
    }
    my $pe_rev_link_cmd = "ln -s $rev_read ./s_1_2_sequence.txt";
    $rv = eval{ Genome::Sys->shellcmd(cmd => $pe_rev_link_cmd); };
    if (not $rv) {
	unlink("s_1_2_sequence.txt");
	die "ERROR! failed to link in rev read (s_1_2_sequence.txt) of read pair\n";
    }

    # build tar archive in the PE allocation directory
    ####my $pe_import_cmd = "tar cvzfh archive.tgz --ignore-failed-read $fwd_read $rev_read"; #NOTE: I removed '$singleton_read' from this command...this instrument data should just have PE reads ... also, fwd_read & rev_read should include full paths to those files
    my $pe_import_cmd = "tar cvzfh archive.tgz --ignore-failed-read s_1_1_sequence.txt s_1_2_sequence.txt"; #NOTE: I had to use symbolic links here, with the tar -h argument, to ensure that the archive.tgz file did not include the full paths to these reads
    print "archiving PE output into PE allocation =>$pe_import_cmd<=\n"; #Log this step
    $rv = eval{ Genome::Sys->shellcmd(cmd => $pe_import_cmd); };
    if ( not $rv ) {
	####unlink($tmp_tar_file);
	unlink("archive.tgz");
	die "ERROR! failed to tar paired end fastq";
    }
    # move temp file to archive file name
    ####$rv = eval{ Genome::Sys->shellcmd(cmd => "mv $tmp_tar_file ".$pe_alloc->absolute_path.'/archive.tgz'); };
    ####if ( not $rv ) {
####	die 'failed to rename PE TEMP.archive.tgz to archive.tgz:'.$@;
    ####}

    # unlink the symbolic links from the PE allocation directory
    unlink("s_1_1_sequence.txt");
    unlink("s_1_2_sequence.txt");


    chdir($current_dir);
    print "manual chdir into directory of execution at =>$current_dir<=\n";


    # make new inst data
    my $sing_original_data_path = $working_dir . "/" . $self->srs_sample_id;
    my $sing_inst_data = Genome::InstrumentData::Imported->create(
	library             => $library,
	import_format       => 'sanger fastq',
	description         => "qualityTrimmed nonRedundant singletonFastq (subset_namesInThisSample:$subset_name_list)", #Made this change to support having variable extrapolated into this object (using single-quotes was just having it store the string '$subset_name_list') ... jmartin 120214
	sequencing_platform => 'solexa',
	is_paired_end       => 0,
	original_data_path  => $sing_original_data_path,
        srs_sample_id         => $self->srs_sample_id
	);

    # make an allocation
    my $singleton_read_size = `du -sk $singleton_read`;
    chomp($singleton_read_size);
    my @singleton_read_size = split(/\s+/,$singleton_read_size);
    my $sing_allocation_request_in_kb = $singleton_read_size[0] * 1.05;
    my $sing_alloc = Genome::Disk::Allocation->create(
	disk_group_name     => $ENV{GENOME_DISK_GROUP_ALIGNMENTS},
	allocation_path     => 'instrument_data/imported/'.$sing_inst_data->id,
	kilobytes_requested => $sing_allocation_request_in_kb,
	owner_class_name    => $sing_inst_data->class,
	owner_id            => $sing_inst_data->id,
    );


#___Changes made for Singleton ... 120113 ... jmartin
    # tar non zipped fastqs to temp file
    ####$tmp_tar_file = $sing_alloc->absolute_path.'/TEMP.archive.tgz';
    ####my $sing_import_cmd = "tar cvzfh $tmp_tar_file -C ".$sing_alloc->absolute_path." $singleton_read";
    $current_dir = `pwd`;
    chomp($current_dir);
    chdir($sing_alloc->absolute_path);
    print "manual chdir into =>$sing_alloc->absolute_path<=\n";


###################NOTE: Changed the archive & copy block to use symbolic link ... jmartin 120123
    # link in the forward & reverse read files for the current pair for archiving
    my $sing_link_cmd = "ln -s $singleton_read ./s_1_sequence.txt";
    $rv = eval{ Genome::Sys->shellcmd(cmd => $sing_link_cmd); };
    if (not $rv) {
	unlink("s_1_sequence.txt");
	die "ERROR! failed to link in singleton read (s_1_sequence.txt) read\n";
    }

    # build tar archive in the PE allocation directory
    ####my $sing_import_cmd = "tar cvzfh archive.tgz --ignore-failed-read $singleton_read"; #NOTE: singleton_read should include full path to that file
    my $sing_import_cmd = "tar cvzfh archive.tgz --ignore-failed-read s_1_sequence.txt"; #NOTE: I had to use a symbolic link here, with the tar -h argument, to ensure that the archive.tgz file did not include the full path to this read
    print "archiving Singleton output into Singleton allocation =>$sing_import_cmd<=\n"; #Log this step
    $rv = eval{ Genome::Sys->shellcmd(cmd => $sing_import_cmd); };
    if ( not $rv ) {
	####unlink($tmp_tar_file);
	unlink("archive.tgz");
	die "ERROR! failed to tar singleton fastq";
    }
    # move temp file to archive file name
    ####$rv = eval{ Genome::Sys->shellcmd(cmd => "mv $tmp_tar_file ".$sing_alloc->absolute_path.'/archive.tgz'); };
    ####if ( not $rv ) {
####	die 'failed to rename SINGLETON TEMP.archive.tgz to archive.tgz:'.$@;
    ####}

    # unlink the symbolic link from the Singleton allocation directory
    unlink("s_1_sequence.txt");



    chdir($current_dir);
    print "manual chdir into directory of execution at =>$current_dir<=\n";




#___Copy metrics for PE & Singleton into /metrics subdirectories with the cleaned fastq data
    my $pe_idid   = $pe_inst_data->id;

    unless ($pe_idid) {
	$self->error_message("did not get instrument data id for paired end fastq ... did this really import?");
	return;
    }
    my $pe_path = $pe_inst_data->allocations->absolute_path;

    unless ($self->copy_metrics($pe_path)) {
	$self->error_message("could not copy paired end metrics into paired end allocation absolute_path");
	return;
    }

    my $sing_idid = $sing_inst_data->id;

    unless ($sing_idid) {
	$self->error_message("did not get instrument data id for singleton fastq ... did this really import?");
	return;
    }
    my $sing_path = $sing_inst_data->allocations->absolute_path;
    unless ($self->copy_metrics($sing_path)) {
	$self->error_message("could not copy singleton metrics into singleton allocation absolute_path");
	return;
    }



#### THIS PART REPLACED WITH ABOVE CODE ... swapped out specific import-fastq methods with generic imports ... 111107

#    my $pe_import_cmd = Genome::InstrumentData::Command::Import::Fastq->create(%import_params,
#								     ####subset_name       => 1,
#								     source_data_files => "$fwd_read,$rev_read",
#								     is_paired_end     => 1,
#								     description       => "qualityTrimmed nonRedundant pairedFastq (subset_namesInThisSample:$subset_name_list)"
#								     );
#    
#
#    unless ($pe_import_cmd->execute) {
#	$self->error_message("Failed to import paired end reads");
#	return;
#    }
#    
#    my $se_import_cmd = Genome::InstrumentData::Command::Import::Fastq->create(%import_params,    ####$self->srs_sample_id
#								     ####subset_name       => 1,
#								     source_data_files => "$singleton_read",
#								     is_paired_end     => 0,
#								     description       => "qualityTrimmed nonRedundant singletonFastq (subset_namesInThisSample:$subset_name_list)"
#								     );
#
#
#
#    unless ($se_import_cmd->execute) {
#	$self->error_message("Failed to import singleton reads");
#	return;
#    }
#
#
#
#    for ($pe_import_cmd, $se_import_cmd) {
#        my $idid = $_->generated_instrument_data_id;
#
#        unless ($idid) {
#            $self->error_message("did not get a generated imported instrument data id.  did this really import?");
#            return;
#        }
#
#
#
#        my $iid = Genome::InstrumentData::Imported->get($idid);
#        unless ($iid) {
#            $self->error_message("could not get the imported instrument data for this id.");
#            return;
#        }
#
#
#        my $path = $iid->allocations->absolute_path;
#    
#
#        unless ($self->copy_metrics($path)) {
#            $self->error_message("could not copy metrics into the destination path!");
#            return;
#        }
#
#
#
#    }






    $self->status_message("Imported all resultant reads and metrics! Done!");

    return 1;
}

sub copy_metrics {
    my $self = shift;
    my $imported_data_path = shift;

    my $working_dir;
   
    $working_dir = $self->_working_dir;

    $working_dir = $working_dir . "/" . $self->srs_sample_id;
    
    my $destination = $imported_data_path . "/metrics";
    
    unless (mkpath($destination)) {
        $self->error_message("Failed to make dest path $destination");
        return; 
    }

#___Changed this for Vervet, since there is no '*masked' files to be moved ... jmartin 111121
    ####my @masks = ("$working_dir/*.masked", "$working_dir/*.denovo_duplicates_marked.counts", "$working_dir/*.denovo_duplicates_marked.metrics", "$working_dir/trimBWAstyle.out");
    my @masks = ("$working_dir/*.denovo_duplicates_marked.counts", "$working_dir/*.denovo_duplicates_marked.metrics", "$working_dir/trimBWAstyle.out");

    my $working_dir_view = `ls -l $working_dir/*`;

    for (@masks) {
       my $cmd = sprintf("rsync -rptgovz --copy-links %s %s", $_, $destination);
       Genome::Sys->shellcmd(cmd=>$cmd);
    }

    
    return 1;
}

sub get_script_path {
    my $self = shift;
    my $file   = __PACKAGE__;
    $file =~ s{::}{/}g;
    $file .= ".pm";

    my $path;
    for my $dir (@INC) {
        $path = "$dir/$file";
        last if -r $path;
        $path = undef;
    }

    $path =~ s/\.pm$//;

#___Added because I put the Vervet processing script in with the normal HMP processing script (under .../HmpSraProcess)
    $path =~ s/Vervet//;

    return $path;
}

#End package
1;
