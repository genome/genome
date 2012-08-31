package Genome::Model::Tools::Tcga::ConvertUnalignedMaqReadsToSamFiles;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::Tcga::ConvertUnalignedMaqReadsToSamFiles {
    is  => ['Command'],
    has => [
        model_id => {
            is  => 'String',
            is_input => '1',
            doc => 'The model id to process.',
        },
        working_directory => {
            is => 'String',
            is_input => '1',
            doc => 'The working directory where results will be deposited.',
        },
        delete_intermediates => {
            is => 'Integer',
            is_input =>1,
            is_optional =>1,
            default=>0,
        },
        unaligned_sam_file_directory => {
            is => 'String',
            is_output => 1,
            is_optional => 1,
            doc => 'The unaligned sam directory.',
        },
        read_group_directory => {
            is => 'String',
            is_output =>1,
            is_optional => 1,
            doc => 'All of the readgroup entries for the TCGA header.',
        },
        program_group_directory => {
            is => 'String',
            is_output =>1,
            is_optional=>1,
            doc => 'All of the programgroup entries for the TCGA header.',
        },
    ],
};


sub help_brief {
    'Convert Maq map files into the TCGA format.';
}

sub help_detail {
    return <<EOS
    Convert Maq map files into the TCGA format.
EOS
}


sub execute {
    my $self = shift;
    $self->dump_status_messages(1);
    $self->dump_error_messages(1);
    $self->dump_warning_messages(1);

    my $now = UR::Time->now;
    $self->status_message(">>>Starting ConvertUnalignedReadsToSamFiles execute() at $now"); 
    
    my $working_dir = $self->working_directory;
    
    my $unaligned_dir = "$working_dir/unaligned";
    my $header_dir = "$working_dir/header"; 
    my $model_id = $self->model_id;

    my $combined_unaligned = "$working_dir/combined_unaligned.sam";
    my $combined_rg = "$working_dir/combined.rg";
    my $combined_pg = "$working_dir/combined.pg";

    my $unaligned_log_file = "$working_dir/logs/UNALIGNED_EXCLUDED.txt";
    my $ua_log = Genome::Sys->open_file_for_writing($unaligned_log_file);
    print $ua_log "This file contains the instrument data inputs which were expected to have unaligned read files, but for some reason could NOT be found.\n";
    print $ua_log "If this file is empty, then all the unaligned reads files associated with the expected instrument data were found.\n\n";
    print $ua_log "ALIGNMENT SEQ ID\tFLOW CELL\tLANE\n";
    
    my $model = Genome::Model->get($model_id);

    die "Model $model_id is not defined. Quitting." unless defined($model);

    #for paired end maq
    my $readName = $1;
    my $readData = $2;

    my $filler= "\t*\t0\t0\t*\t*\t0\t0\t";
    my $pair1 = "\t69";
    my $pair2 = "\t133";
    my $frag = "\t4";

    my $rg_tag = "\tRG:Z:";
    my $pg_tag = "\tPG:Z:";

    my @instrument_data = $model->instrument_data;
    my $alignment_count = scalar(@instrument_data);
    my $build = $model->last_complete_build;
    unless ($build) {
        die "No last complete build for model " . $model->__display_name__;
    }

    my $count = 0;
    for my $instrument_data (@instrument_data) {
        my $alignment = $build->alignment_results_for_instrument_data($instrument_data);
        my $idid = $instrument_data->id;
        my $flow_cell = $instrument_data->flow_cell_id;
        my $lane = $instrument_data->lane;
        my $unaligned_file = $alignment->unaligned_reads_list_path;

        if (!-s $unaligned_file) {
            $self->status_message("* * * WARNING  $unaligned_file does not exist for $idid.");
            print $ua_log "$idid\t$flow_cell\t$lane\n"; 
        }
        next if (!-s $unaligned_file);       
 
        my $unaligned_sam_file = "$unaligned_dir/$idid.sam";
        if (-s $unaligned_sam_file) {
            $self->status_message("*** $unaligned_sam_file already exists for $idid.  Skipping processing.");
        } 
        next if (-s $unaligned_sam_file);
 
        my $ua_fh = Genome::Sys->open_file_for_writing($unaligned_sam_file);

        #calculate all the parameters and header info
        my $insert_size_for_header;
        if ($instrument_data->median_insert_size) {
            $insert_size_for_header= $instrument_data->median_insert_size;
        } else {
            $insert_size_for_header = 0;
        }
        
        my $description_for_header;
        if ($instrument_data->is_paired_end) {
            $description_for_header = 'paired end';
        } else {
            $description_for_header = 'fragment';
        }

        my $date_run_field = $instrument_data->run_start_date_formatted();
        my $sample_name_field = $instrument_data->sample_name;
        my $library = $instrument_data->library_name;
        my $platform_unit_field = sprintf("%s.%s",$instrument_data->flow_cell_id,$instrument_data->lane);
        
        my $aligner_version = $alignment->aligner_version;

        my $aligner_params;

        my $seed = 0;
        for my $c (split(//,$instrument_data->flow_cell_id || $alignment->instrument_data_id)) {
            $seed += ord($c)
        }
        $seed = $seed % 65536;
        $aligner_params .= " -s $seed ";

        ##### -a param
        my $upper_bound_on_insert_size;
        if ($instrument_data->is_paired_end && !$alignment->force_fragment) {
            my $sd_above = $instrument_data->sd_above_insert_size;
            my $median_insert = $instrument_data->median_insert_size;
            $upper_bound_on_insert_size= ($sd_above * 5) + $median_insert;
            unless($upper_bound_on_insert_size > 0) {
                #$self->status_message("Unable to calculate a valid insert size to run maq with. Using 600 (hax)");
                $upper_bound_on_insert_size= 600;
            }
        $aligner_params .= ' -a '. $upper_bound_on_insert_size;
        }
      
        #adaptor param 
        my $adaptor_file = $instrument_data->resolve_adaptor_file;
        if ($adaptor_file) {
            $aligner_params .= ' -d '. $adaptor_file;
        }
     
        my $cmdline = Genome::Model::Tools::Maq->path_for_maq_version($alignment->aligner_version) . " map $aligner_params "; 


        $self->status_message("\n**********************************************************************************************************\n"); 
        $self->status_message("*** unaligned: $unaligned_file\n");
        $self->status_message("*** RG: $idid\n");
        $self->status_message ("*** LB: $library\n");
        $self->status_message ("*** PI: $insert_size_for_header\n");
        $self->status_message ("*** DS: $description_for_header\n");
        $self->status_message ("*** SM: $sample_name_field\n");
        $self->status_message ("*** PU: $platform_unit_field\n");
        $self->status_message ("*** DR: $date_run_field\n");
        $self->status_message ("*** AV: $aligner_version\n");
        $self->status_message ("*** AP: $aligner_params\n");
        $self->status_message ("***CMD: $cmdline\n");
       
        my $rg_string = "\@RG\tID:$idid\tPL:illumina\tPU:$platform_unit_field\tLB:$library\tPI:$insert_size_for_header\tDS:$description_for_header\tDT:$date_run_field\tSM:$sample_name_field\tCN:WUGSC\n";

        #program group
        my $pg_string ="\@PG\tID:$idid\tVN:$aligner_version\tCL:$cmdline\n";

        my $rg_file = "$header_dir/$idid.rg";
        my $pg_file = "$header_dir/$idid.pg";

        if (!-s $rg_file) {
            my $rg_fh = Genome::Sys->open_file_for_writing($rg_file);
            print $rg_fh $rg_string; 
            $rg_fh->close;
        }
    
        if (!-s $pg_file) { 
            my $pg_fh = Genome::Sys->open_file_for_writing($pg_file);
            print $pg_fh $pg_string; 
            $pg_fh->close;
        }


        my $fh = Genome::Sys->open_file_for_reading("$unaligned_file");
        if ($instrument_data->is_paired_end) {
            while (my $line1 = $fh->getline) {
                my $line2 = $fh->getline; 
                $line1 =~ m/(^\S.*)\t99\t(.*$)/;
                my $readName1 = $1;
                my $readData1 = $2;
                $line2 =~ m/(^\S.*)\t99\t(.*$)/;
                my $readName2 = $1;
                my $readData2 = $2;

                print $ua_fh $readName1.$pair1.$filler.$readData1.$rg_tag.$idid.$pg_tag.$idid."\n";
                print $ua_fh $readName2.$pair2.$filler.$readData2.$rg_tag.$idid.$pg_tag.$idid."\n";

                $count++;
            }
        } else {
            while (my $line1 = $fh->getline) {
                $line1 =~ m/(^\S.*)\t99\t(.*$)/;
                my $readName1 = $1;
                my $readData1 = $2;
                print $ua_fh $readName1.$frag.$filler.$readData1.$rg_tag.$idid."\n";
                $count++;
            }
        }


        #if ($count == 2) {
        #    die "Dieing for testing purposes...";
        #}

 
    }

    $self->unaligned_sam_file_directory($unaligned_dir);
    $self->read_group_directory($header_dir);
    $self->program_group_directory($header_dir);

    print "\nDone with unaligned map to sam conversion.\n";

    $now = UR::Time->now;
    #$self->dump_status_messages(1); 
    $self->status_message("<<<Completed ConvertUnalignedReadsToSamFiles execute() at $now"); 
    $ua_log->close;
    return 1;
}
1;
