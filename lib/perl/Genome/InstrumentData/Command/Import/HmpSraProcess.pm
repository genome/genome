package Genome::InstrumentData::Command::Import::HmpSraProcess;

use strict;
use warnings;
use Genome;
use Cwd;
use IO::File;
use File::Basename;
use File::Path;

class Genome::InstrumentData::Command::Import::HmpSraProcess {
    is  => 'Command',
    has_input => [
	ascp_user => {
	    is_optional => 0,
	    doc => 'DACC FTP user_name for aspera transfer',
	},
	ascp_pw => {
	    is_optional => 0,
	    doc => 'DACC FTP password for aspera transfer',
	},
	srs_sample_id => {
	    is_optional => 0,
	    doc => 'SRS sample id to extract',
	},
	srr_accessions => {
	    is_optional => 0,
	    doc => 'space separated list of SRR accession ids for the raw SRA data downloads to use.',
	},
	picard_dir => {
	    is_optional => 1,
	    doc => 'full path to directory containing Picard jar files (note: This path must include the updated EstimateLibraryComplexity that handles redundancy removal)',
            ####default_value => Genome::Model::Tools::Picard->path_for_picard_version,
	    default_value => "$ENV{GENOME_SW_LEGACY_JAVA}/samtools/picard-tools-1.27",
	},
        container_dir => {
            is_optional => 1,
            doc => 'skip processing, pull processed SRS directory from this location.',
        },
    ],
   has => [
     _working_dir => {
        is_optional=>1,
        is_transient=>1
     }
    ],
    doc => 'de-duplicate and quality trim Illumina WGS runs downloaded from SRA',
};

sub execute {
    my $self = shift;


#___This line stops the perl debugger as though I'd set a break point in the GUI
    $DB::single = 1;
    
    my @srrs = split /\s+/, $self->srr_accessions;

    my %samples;

    for my $line (@srrs) {
##_______This does not work as of 110927...jmartin
#        my $sample = Genome::Sample->get(sql=>qq/
#            select os.*
#            from organism_sample os
#            join sra_organism_sample sos on sos.organism_sample_id=os.organism_sample_id
#            join sra_experiment ex on ex.sra_sample_id=sos.sra_sample_id
#            join sra_run ru on ru.sra_experiment_id=ex.sra_item_id
#            join sra_item rui on rui.sra_item_id=ru.sra_item_id
#            join sra_accession ruacc on ruacc.alias=rui.alias
#            where ruacc.accession='$line'
#        /);
#_______Modified to use a 'workaround' to deal with 'deprecated Genome gets'...jmartin 110927
        my $dbh = Genome::DataSource::GMSchema->get_default_handle();
        my ($sample) = $dbh->selectrow_array(qq/
            select os.organism_sample_id from organism_sample os
            join sra_organism_sample sos
            on sos.organism_sample_id=os.organism_sample_id
            join sra_experiment ex
            on ex.sra_sample_id=sos.sra_sample_id
            join sra_run ru
            on ru.sra_experiment_id=ex.sra_item_id
            join sra_item rui
            on rui.sra_item_id=ru.sra_item_id
            join sra_accession ruacc
            on ruacc.alias=rui.alias where ruacc.accession='$line'
            /);

        unless ($sample) {
            $self->error_message("Failed to get a sample object from the warehouse for SRR ID $line.");
            return;
        }

#_______The 'name' attribute doesn't seem to exist, but full_name does, and this value seems only used for the sanity check immediately below, so this change should be OK ... jmartin 110929
        ####$samples{$sample->name} = 1;
        $samples{$sample->full_name} = 1;
    }

    unless (scalar keys %samples == 1) {
        $self->error_message("There is more than one sample represented in this SRR list, can't proceed.  Samples are: " . join "\n", keys %samples);
        return;
    }

    # grab import parameters
    my $dbh = Genome::DataSource::GMSchema->get_default_handle();
    my ($fc_id, $lane) = $dbh->selectrow_array(qq/select ii.flow_cell_id, ii.lane
        from organism_sample os
        join sra_organism_sample sos on sos.organism_sample_id=os.organism_sample_id
        join sra_experiment ex on ex.sra_sample_id=sos.sra_sample_id
        join sra_run ru on ru.sra_experiment_id=ex.sra_item_id
        join sra_item rui on rui.sra_item_id=ru.sra_item_id
        join sra_accession ruacc on ruacc.alias=rui.alias
        join index_illumina ii on ii.seq_id=rui.source_entity_id and rui.source_entity_type='index illumina'
        where ruacc.accession='$srrs[0]'/);

    unless (defined $fc_id && defined $lane) {
        $self->error_message("Couldn't recover original flow cell id and lane for this SRR id $srrs[0]");
        return;
    }

    my $original_inst = Genome::InstrumentData::Solexa->get(flow_cell_id=>$fc_id, lane=>$lane);

    my %import_params;
    $import_params{library_name} = $original_inst->library_name;
    $import_params{sample_name} = $original_inst->sample_name;
    $import_params{sequencing_platform} = 'solexa';
    $import_params{import_format} = 'sanger fastq';
    $import_params{sra_sample_id} = $self->srs_sample_id;
  
    my $tmp_dir      = Genome::Sys->create_temp_directory;

    my $working_dir = $tmp_dir . "/srr_datasets";
    $self->_working_dir($working_dir);
    mkpath($working_dir);
    for my $srr (@srrs) {
        my $instrument_data = Genome::InstrumentData::Imported->get(import_format=>'raw sra download', sra_accession=>$srr);
        
        my ($alloc) = $instrument_data->allocations; 
        unless(symlink($alloc->absolute_path . "/" . $srr, $working_dir . "/" . $srr)) {
	    $self->error_message("Failed to set up symlink from SRA data dir: " . $alloc->absolute_path . " to " . $working_dir . "/" . $srr);
	   return;
	}
    } 

    my ($fwd_read, $rev_read, $singleton_read);

    if ($self->container_dir) {
        unless (-d $self->container_dir . "/" . $self->srs_sample_id) {
            $self->error_message("Can't find srs sample id in " . $self->container_dir);
            return;
        }
    
        my ($fwd_read_bz) = glob($self->container_dir . "/" . $self->srs_sample_id . "/*.trimmed.1.fastq.bz2");
        my ($rev_read_bz) = glob($self->container_dir . "/" . $self->srs_sample_id . "/*.trimmed.2.fastq.bz2");
        my ($singleton_read_bz) = glob($self->container_dir . "/" . $self->srs_sample_id . "/*.trimmed.singleton.fastq.bz2");

        $fwd_read = $working_dir . "/s_1_1_sequence.txt";
        Genome::Sys->shellcmd(cmd=>"bzcat $fwd_read_bz > $fwd_read",
                                      output_files=>[$fwd_read]);
        
        $rev_read = $working_dir . "/s_1_2_sequence.txt";
        Genome::Sys->shellcmd(cmd=>"bzcat $rev_read_bz > $rev_read",
                                      output_files=>[$rev_read]);
        
        $singleton_read = $working_dir . "/s_1_sequence.txt";
        Genome::Sys->shellcmd(cmd=>"bzcat $singleton_read_bz > $singleton_read",
                                      output_files=>[$singleton_read]);
        

    } else {

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

        #Get ascp user/pwd
        my $user = $self->ascp_user;
        my $pwd  = $self->ascp_pw;




        #Note: I need to set the path to my scripts INSIDE the shell command
        ####$cmd = "cd $working_dir; export PATH=$path; process_runs.sh $list_of_srrs $sra_samples $picard_dir $tmp_dir > $outfile 2> $errfile";

#_______THIS LINE IS FOR REAL DATA...UNCOMMENT TO REALLY DO PROCESSING AND UPLOAD TO DACC!!!! jmartin ... 110927
        ####$cmd = "cd $working_dir; export PATH=$path; process_runs.sh $list_of_srrs $sra_samples $picard_dir $tmp_dir $user $pwd";

#_______This line is just for testing, the DACC upload is commented out in this version jmartin ... 110923
	$cmd = "cd $working_dir; export PATH=$path; process_runs.for_testing.sh $list_of_srrs $sra_samples $picard_dir $tmp_dir $user $pwd";




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

        rename($fwd_read, $working_dir . "/s_1_1_sequence.txt");
        $fwd_read = $working_dir . "/s_1_1_sequence.txt";
        
        rename($rev_read, $working_dir . "/s_1_2_sequence.txt");
        $rev_read = $working_dir . "/s_1_2_sequence.txt";
        
        rename($singleton_read, $working_dir . "/s_1_sequence.txt");
        $singleton_read = $working_dir . "/s_1_sequence.txt";

        $DB::single = 1;

    }
    
    my $pe_import_cmd = Genome::InstrumentData::Command::Import::Fastq->create(%import_params,
								     subset_name=>1,
								     source_data_files=>"$fwd_read,$rev_read",
								     is_paired_end=>1
								     );
    
    unless ($pe_import_cmd->execute) {
	$self->error_message("Failed to import paired end reads");
	return;
    }
    
    my $se_import_cmd = Genome::InstrumentData::Command::Import::Fastq->create(%import_params,
								     subset_name=>1,
								     source_data_files=>"$singleton_read",
								     is_paired_end=>0
								     );

    unless ($se_import_cmd->execute) {
	$self->error_message("Failed to import singleton reads");
	return;
    }

    for ($pe_import_cmd, $se_import_cmd) {
        my $idid = $_->generated_instrument_data_id;

        unless ($idid) {
            $self->error_message("did not get a generated imported instrument data id.  did this really import?");
            return;
        }

        my $iid = Genome::InstrumentData::Imported->get($idid);
        unless ($iid) {
            $self->error_message("could not get the imported instrument data for this id.");
            return;
        }

        my $path = $iid->allocations->absolute_path;
    
        unless ($self->copy_metrics($path)) {
            $self->error_message("could not copy metrics into the destination path!");
            return;
        }
    }

    $self->status_message("Imported all resultant reads and metrics! Done!");

    return 1;
}

sub copy_metrics {
    my $self = shift;
    my $imported_data_path = shift;

    my $working_dir;
    if ($self->container_dir) {
        $working_dir = $self->container_dir;
    } else {
        $working_dir = $self->_working_dir;
    }

    $working_dir = $working_dir . "/" . $self->srs_sample_id;
    
    my $destination = $imported_data_path . "/metrics";
    
    unless (mkpath($destination)) {
        $self->error_message("Failed to make dest path $destination");
        return; 
    }

    my @masks = ("$working_dir/*.masked", "$working_dir/*.denovo_duplicates_marked.counts", "$working_dir/*.denovo_duplicates_marked.metrics", "$working_dir/trimBWAstyle.out");

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

    return $path;
}

#End package
1;
