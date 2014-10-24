package Genome::InstrumentData::AlignmentResult::Imported;

use strict;
use warnings;
use File::Basename;

use Genome;

class Genome::InstrumentData::AlignmentResult::Imported {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'Imported', is_param=>1 },
    ],
    has_transient_optional => [
         _bwa_sam_cmd => { is=>'Text' }
    ]
};

sub required_arch_os { 'x86_64' }

sub required_rusage { 
    "-R 'select[mem>10000] span[hosts=1] rusage[mem=10000]' -M 10000000";
}

# No need to extract fastqs & run aligner

sub collect_inputs {
    return 1;
}

sub run_aligner {
    return 1;
}

sub prepare_scratch_sam_file {
    return 1;
}

#Overload the creation of BAM's in AlignmentResult

sub create_BAM_in_staging_directory {
    my $self = shift;
    my $tmp_dir = $self->temp_staging_directory;
    my $instrument_data_id = $self->instrument_data_id;
    my $instrument_data;
    
    if(defined($instrument_data_id)) {
        $self->status_message("Found instrument-data.\n");
        print $self->status_message;
        
        $instrument_data = Genome::InstrumentData::Imported->get(id=>$instrument_data_id);
        unless(defined($instrument_data)) {
            $self->error_message(" Could not locate an instrument data record with an ID of ".$instrument_data_id."\n");
            die $self->error_message;
        }
        my $bam_output_path = $tmp_dir."/all_sequences.bam";
        unless(-e $bam_output_path) {
            $self->status_message("No all_sequences.bam file found at ".$bam_output_path." attempting to create one.");

            #If the source of the imported BAM is broad, apply the debroadifyBam tool to it
            if(($instrument_data->import_source_name =~ /broad/i)||($instrument_data->import_source_name =~ /bcm/i)) {
                $self->status_message("Import source is ".$instrument_data->import_source_name.", call Debroadify tool and convert the broadiferous bam to GSC style.\n");
                my $debroadify_cmd = Genome::Model::Tools::Sam::Debroadify->create(
                    input_file  => $instrument_data->data_directory . "/all_sequences.bam",
                    output_file => $bam_output_path,
                    reference_file  => $self->reference_build->full_consensus_path('fa'),
                );
                unless ($debroadify_cmd->execute) {
                    die $self->error_message("Debroadify failed to complete.");
                }
                unless(-e $bam_output_path) {
                    die $self->error_message("Could not find all_sequences.bam after running the Debroadify tool.");
                }   
                $self->status_message("Successfully created an all_sequences.bam file.");
            } else {     #otherwise simply symlink the bam into the staging directory 
                $self->status_message("Attempting to symlink file from import to alignment scratch dir.\n");
                unless(Genome::Sys->create_symlink($instrument_data->data_directory."/all_sequences.bam",$bam_output_path)) {
                    $self->error_message("Failed to symlink BAM from instrument-data allocation to alignment scratch dir.\n");
                    die $self->errror_message;
                }
                unless(-e $bam_output_path) {  
                    $self->error_message("Could not verify symlink of BAM to staging directory.\n");
                    die $self->error_message;
                }
            }
        }
    } else {
       $self->error_message("Alignment has no instrument data.\n");
        die $self->error_message;
    }
    return 1; 
}

sub _filter_samxe_output {
    my ($self, $sam_cmd, $sam_file_name) = @_;

    my $sam_run_output_fh = IO::File->new( $sam_cmd . "|" );
    $self->status_message("Running $sam_cmd");
    if ( !$sam_run_output_fh ) {
            $self->error_message("Error running $sam_cmd $!");
            return;
    }


    my $sam_map_output_fh = IO::File->new(">>$sam_file_name");
    if ( !$sam_map_output_fh ) {
            $self->error_message("Error opening sam file for writing $!");
            return;
    }
    $self->status_message("Opened $sam_file_name");
    
    while (<$sam_run_output_fh>) {
            #write out the aligned map, excluding the default header- all lines starting with @.
            my $first_char = substr($_,0,1);
                if ($first_char ne '@') {
                $sam_map_output_fh->print($_);
            }
    }
    $sam_map_output_fh->close;
    return 1;
}


sub _verify_bwa_aln_did_happen {
    my $self = shift;
    my %p = @_;

    unless (-e $p{sai_file} && -s $p{sai_file}) {
        $self->error_message("Expected SAI file is $p{sai_file} nonexistent or zero length.");
        return;
    }
    
    unless ($self->_inspect_log_file(log_file=>$p{log_file},
                                     log_regex=>'(\d+) sequences have been processed')) {
        
        $self->error_message("Expected to see 'X sequences have been processed' in the log file where 'X' must be a nonzero number.");
        return 0;
    }
    
    return 1;
}

sub _inspect_log_file {
    my $self = shift;
    my %p = @_;

    my $aligner_output_fh = IO::File->new($p{log_file});
    unless ($aligner_output_fh) {
        $self->error_message("Can't open expected log file to verify completion " . $p{log_file} . "$!"
        );
        return;
    }
    
    my $check_nonzero = 0;
    
    my $log_regex = $p{log_regex};
    if ($log_regex =~ m/\(\\d\+\)/) {
        
        $check_nonzero = 1;
    }

    while (<$aligner_output_fh>) {
        if (m/$log_regex/) {
            $aligner_output_fh->close();
            if ( !$check_nonzero || $1 > 0 ) {
                return 1;
            }
            return;
        }
    }

    return;
}

sub decomposed_aligner_params {
    my $self = shift;
    my $params = $self->aligner_params || ":::";
    
    my @spar = split /\:/, $params;
    
    
    return ('bwa_aln_params' => $spar[0], 'bwa_samse_params' => $spar[1], 'bwa_sampe_params' => $spar[2]);
}

sub aligner_params_for_sam_header {
    my $self = shift;
    
    my %params = $self->decomposed_aligner_params;
    my $aln_params = $params{bwa_aln_params} || "";
    
    my $sam_cmd = $self->_bwa_sam_cmd || "";

    return "bwa aln $aln_params; $sam_cmd ";
}

sub fillmd_for_sam {
    return 1;
}

sub _compute_alignment_metrics {
    return 1;
}

sub _check_read_count { 
    return 1; 
} 
