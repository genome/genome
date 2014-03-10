package Genome::InstrumentData::AlignmentResult::Clc;

use strict;
use warnings;
use File::Basename;
use IO::File;
use File::Path qw(make_path);
use File::Copy;

use Genome;

class Genome::InstrumentData::AlignmentResult::Clc {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'clc', is_param=>1 },
    ],
    has_transient_optional => [
         _Clc_sam_cmd => { is=>'Text' }
    ]
};

sub required_arch_os { 'x86_64' }

sub required_rusage {
    my $class = shift;
    my %p = @_;
    my $instrument_data = delete $p{instrument_data};

    my $estimated_usage_mb = 90000;
    if (defined $instrument_data && $instrument_data->can("calculate_aligned_estimated_kb_usage")) {
        my $kb_usage = $instrument_data->calculate_alignment_estimated_kb_usage;
        $estimated_usage_mb = int(($kb_usage * 5) / 1024)+100;
    }
    my $mem_thousand = 14000;
    my $mem_million = 15000000;
    if (defined $instrument_data && $instrument_data->is_paired_end){
        $mem_thousand = 45000;
        $mem_million = 47000000;
    }
        
    return "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>" . $estimated_usage_mb . " && mem>$mem_thousand] span[hosts=1] rusage[mem=$mem_thousand]' -M $mem_million -n 4 -q hmp -m hmp";
}

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;

    $self->copy_license_file_to_home_directory();


    my $tmp_dir = $self->temp_scratch_directory;

    my $tmp_cas_file = $tmp_dir . "/" . "aligned.cas"; 
    my $tmp_sam_file = $tmp_dir . "/" . "aligned.sam";
    my $staging_sam_file = $tmp_dir . "/" . "all_sequences.sam";
    
    # get refseq info
    my $reference_build = $self->reference_build;
    my $ref_path = $self->reference_build->full_consensus_path('fa');

    my $clc_path = '/gsc/bin/clc_ref_assemble_long';
    my $castosam_path = '/gscmnt/233/analysis/sequence_analysis/species_independant/jmartin/CLC/server_based/castosam';
    my $aligner_params = $self->aligner_params;



#___Fix headers in reference DB to keep only the first, white-space delimited word as the name
    my $ref_db_fh = new IO::File $ref_path;
    my $repaired_db_name = $tmp_dir . "/repaired_all_sequences.fa";
    my $repaired_db_name_fh = new IO::File ">>$repaired_db_name";
    while (<$ref_db_fh>) {
	chomp;
	my $line = $_;
	if ($line =~ /^>/) {
	    my @line = split(/\s+/,$line);
	    print $repaired_db_name_fh "$line[0]\n";
	} else {
	    print $repaired_db_name_fh "$line\n";
	}
    }
    $ref_db_fh->close;
    $repaired_db_name_fh->close;


#___Figure out whether the input is paired end or fragment (can't do both)
    my $align_cmd; #This is assuming the user will enter ONLY the additional params....such as -l & -s.   For now I will just hardcode in the min & max insert size ranges
    if (@input_pathnames == 2) {
	####$align_cmd = "$clc_path -o $tmp_cas_file -d $ref_path -q -p fb ss 180 250 $aligner_params -i @input_pathnames --cpus 4";
	$align_cmd = "$clc_path -o $tmp_cas_file -d $repaired_db_name -q -p fb ss 180 250 $aligner_params -i @input_pathnames --cpus 4";
    } elsif (@input_pathnames == 1) {
	####$align_cmd = "$clc_path -o $tmp_cas_file -d $ref_path -q @input_pathnames $aligner_params --cpus 4";
	$align_cmd = "$clc_path -o $tmp_cas_file -d $repaired_db_name -q @input_pathnames $aligner_params --cpus 4";
    }

    #my $align_cmd = $clc_path . sprintf(' -o %s -d %s -q %s -p fb ss %s %s -i %s %s %s',$tmp_cas_file,$ref_path,$frag_input,$insert_size_range_min,$insert_size_range_max,$pe1_input,$pe2_input,$other_params);



#___Run alignment
    Genome::Sys->shellcmd(
	cmd             => $align_cmd,
	####input_files     => [ $ref_path, @input_pathnames ],
	input_files     => [ $repaired_db_name, @input_pathnames ],
	output_files    => [ $tmp_cas_file ],
	skip_if_output_is_present => 0,
        );


#___Convert cas alignment output into sam file
    my $convert_cmd = $castosam_path . sprintf(' -a %s -o %s -f 33',$tmp_cas_file,$tmp_sam_file);

    Genome::Sys->shellcmd(
	cmd             => $convert_cmd,
	input_files     => [ $tmp_cas_file ],
	output_files    => [ $tmp_sam_file ],
	skip_if_output_is_present => 0,
        );

    die "Failed to process sam command line, error_message is ".$self->error_message unless $self->_filter_sam_output($tmp_sam_file, $staging_sam_file); #MAKE SURE THIS IS APPENDING INTO FINAL STAGING FILE (not overwriting!)

    return 1;
}

sub copy_license_file_to_home_directory {
    my $self = shift;
    $self->status_message("copying latest clc license to user directory");
    ####this method grabs the latest clc license and copies it to the users home directory
    my $allocation_id = 'f392a8fcbbf94b0f865bd172d51ae096'; #This is good through November 2011
    my $allocation = Genome::Disk::Allocation->get($allocation_id);
    die $self->error_message("no allocation for clc license") unless $allocation;
    my $license_filename = "CLC_Assembly_Cell_JCVI-Singh-HMPConsortium.lic";
    my $allocated_license = $allocation->absolute_path."/$license_filename";
    die $self->error_message("no clc license file in allocation") unless -e $allocated_license;
    my $user = getpwuid($<);
    my $local_path = "/gscuser/$user/.clcbio/licenses/";
    make_path($local_path);
    die $self->error_message("unable to copy clc license to users home directory") unless copy($allocated_license, $local_path);
    return 1;
}


sub _filter_sam_output {
    my ($self, $cur_sam_file, $all_sequences_sam_file) = @_;

#    my $sam_run_output_fh = IO::File->new( $sam_cmd . "|" );
#    $self->debug_message("Running $sam_cmd");
#    if ( !$sam_run_output_fh ) {
#            $self->error_message("Error running $sam_cmd $!");
#            return;
#    }

    my $cur_sam_fh = IO::File->new( $cur_sam_file );
    if ( !$cur_sam_fh ) {
        $self->error_message("Error opening current sam file for reading $!");
        return;
    }
    $self->debug_message("Opened $cur_sam_file");

    my $all_seq_fh = IO::File->new(">>$all_sequences_sam_file");
    if ( !$all_seq_fh ) {
        $self->error_message("Error opening all seq sam file for writing $!");
        return;
    }
    $self->debug_message("Opened $all_sequences_sam_file");
    
    while (<$cur_sam_fh>) {
        #write out the aligned map, excluding the default header- all lines starting with @.
        $all_seq_fh->print($_) unless $_ =~ /^@/;
    }

    $cur_sam_fh->close;
    $all_seq_fh->close;
    return 1;
}

sub aligner_params_for_sam_header {
    my $self = shift;
    
    my $params = $self->aligner_params;
    
    return "clc $params";

}

sub fillmd_for_sam {
#___Tells AlignmentResult.pm to run samtools fillmd
    return 0;
}

sub prepare_reference_sequence_index {
    my $class = shift;

    $class->debug_message("CLC doesn't need any index made, doing nothing.");

}



1;
