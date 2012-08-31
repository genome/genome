package Genome::InstrumentData::AlignmentResult::Mblastx;

use strict;
use warnings;
use File::Basename;
use File::Path;
use File::Copy;
use Genome;

my %VERSIONS = (
    '09242010' => "/gscmnt/sata895/research/mmitreva/SOFTWARE/MCW_09242010",
    '12072010' => "/gscmnt/sata895/research/mmitreva/SOFTWARE/MCW_12072010",
    '04272011' => "/gscmnt/sata895/research/mmitreva/SOFTWARE/MCW_04272011",
);

class Genome::InstrumentData::AlignmentResult::Mblastx {
    is => 'Genome::InstrumentData::AlignmentResult',

    has_constant => [ 
        aligner_name => { 
            value => 'mblastx', 
            is_param => 1, 
        }, 
    ],
    has => [
        _max_read_id_seen => { 
            default_value => 0,
            is_optional => 1, 
        },
        _file_input_option => { 
            default_value => 'fasta', 
            is_optional => 1, 
        },
    ],
};

sub required_arch_os { 'x86_64' }

sub required_rusage {
"-R 'select[model!=Opteron250 && type==LINUX64 && tmp>90000 && mem>32000] span[hosts=1] rusage[tmp=90000, mem=32000]' -M 32000000 -n 8 -m hmp -q hmp";
}

sub _decomposed_aligner_params {
	my $self = shift;

	my $aligner_params = ( $self->aligner_params || '' );

	my $cpu_count = $self->_available_cpu_count;
	$aligner_params .= " -T $cpu_count";

	return ( 'mblastx_aligner_params' => $aligner_params );
}

sub aligner_params_for_sam_header {
    my $self = shift;
    
    my %params = $self->_decomposed_aligner_params;
    my $aln_params = $params{mblastx_aligner_params};
    my $mblastx = $self->_mblastx_path;
    
    return "$mblastx $aln_params";
}

sub _mblastx_path {
    my $self = shift;
    my $version = $self->aligner_version;
    return $VERSIONS{$version}."/mblastx";
}

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;
    if ( @input_pathnames == 1 ) {
        $self->status_message("_run_aligner called in single-ended mode.");
    } 
    elsif ( @input_pathnames == 2 ) {
        $self->status_message("_run_aligner called in paired-end mode.  We don't actually do paired alignment with MBlastx though; running two passes.");
    }
    else {
        die $self->error_message( "_run_aligner called with ".scalar(@input_pathnames)." files.  It should only get 1 or 2!" );
    }

    # get refseq info

    my $reference_build = $self->reference_build;
    my $reference_name = $reference_build->prefix;
    my $reference_mblastx_path = $reference_build->data_directory . '/mblastx';

    my $scratch_directory = $self->temp_scratch_directory;
    my $staging_directory = $self->temp_staging_directory;

    my @mblastx_input_fastas;

    for my $i ( 0 ... $#input_pathnames ) {
        my $input_pathname = $input_pathnames[$i];
        my $chunk_path = $scratch_directory . "/chunks/chunk-from-" . $i;
        unless ( mkpath($chunk_path) ) {
            die $self->error_message("couldn't create a place to chunk the data in $chunk_path");
        }
        my $input_fasta = File::Temp::tempnam( $scratch_directory, "input-XXX" ). ".fasta";


        my $fastq_to_fasta =  Genome::Model::Tools::Fastq::ToFasta->create(
                                fastq_file => $input_pathname,
                                fasta_file => $input_fasta);
        
        unless ($fastq_to_fasta->execute && -s $input_fasta){
                die  $self->error_message("Failed to convert fastq $input_pathname to fasta");
        }
        push @mblastx_input_fastas, $input_fasta;

    }

    for my $input_fasta (@mblastx_input_fastas) {
        my $output_file =  $self->temp_scratch_directory."/".basename($input_fasta)."_vs_$reference_name"."_mblastx.out";

        #STEP 2 - run mblastx aligner
        my %aligner_params = $self->_decomposed_aligner_params;

        my $mblastx = $self->_mblastx_path;
        my $mblastx_aligner_params = (
                defined $aligner_params{'mblastx_aligner_params'}
                ? $aligner_params{'mblastx_aligner_params'}
                : ""
        );
        my $cmd = sprintf( '%s -q %s -o %s %s',
                $mblastx, $input_fasta, $output_file, $mblastx_aligner_params );

        local $ENV{MBLASTX_DATADIR} = $reference_mblastx_path;
        $self->status_message("mblastx data dir variable set to $ENV{MBLASTX_DATADIR}");
        
        Genome::Sys->shellcmd(
                cmd                       => $cmd,
                input_files               => [ $input_fasta ],
                output_files              => [ $output_file ],
                skip_if_output_is_present => 0,
        );
        rename($output_file,$self->temp_staging_directory."/".basename($output_file));
    }

    return 1;
}

sub _compute_alignment_metrics {
	return 1;
}

sub prepare_scratch_sam_file {
        return 1;
}

sub create_BAM_in_staging_directory {
	return 1;
}

sub postprocess_bam_file {
	return 1;
}

sub _prepare_reference_sequences {
    my $self = shift;
    my $reference_build = $self->reference_build;

    my $dir = $reference_build->data_directory . '/mblastx';
    if ( -e $dir ) {
        $self->status_message("Found reference data at: $dir");
        return 1;
    }
    $self->status_message("No reference data found at: $dir");
    mkpath($dir);
    my $ref_basename = File::Basename::fileparse( $reference_build->full_consensus_path('fa') );
    my $reference_fasta_path = sprintf( "%s/%s", $reference_build->data_directory, $ref_basename );

    unless ( -e $reference_fasta_path ) {
        $self->error_message("Alignment reference path $reference_fasta_path does not exist");
            die $self->error_message;
    }

    # generate a reference data set at $dir...
    my $mhashgen = $VERSIONS{$self->aligner_version}."/mhashgen";

    unless($self->mhashgen_format){
        die $self->error_message("Database format option (--mhashgen-format=K/N/A) not provided.");
    }

    my $cmd = sprintf( '%s -s %s -T %s', $mhashgen, $reference_fasta_path, $self->mhashgen_format );

    local $ENV{MBLASTX_DATADIR} = $dir;

    #todo figure out how to copy the BLOSUM matrix

    Genome::Sys->shellcmd(
            cmd                       => $cmd,
            input_files               => [$reference_fasta_path],
            skip_if_output_is_present => 0,
    );

    $self->status_message("Reference data generation complete at: $dir");
    my $dat_file = $dir."/BLOSUM62_6_26.dat";
    my $dat_file_source = $VERSIONS{$self->aligner_version} . "/BLOSUM62_6_26.dat";
    unless(-e $dat_file){
        unless(copy($dat_file_source,$dat_file)){
            die $self->error_message("Failed to copy dat file into destination directory: ".$dir);
        }
    }
    
    
    return 1;
}

