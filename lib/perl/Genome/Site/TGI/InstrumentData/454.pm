package Genome::Site::TGI::InstrumentData::454;

use strict;
use warnings;

use Genome;

require Carp;

class Genome::Site::TGI::InstrumentData::454 {
    is  => 'Genome::Site::TGI::InstrumentData',
    table_name => <<'EOS'
        (
            select 
                to_char(case when ri.index_sequence is null then ri.region_id else ri.seq_id end) id,
                '454' sequencing_platform,
                r.region_id genome_model_run_id, --legacy
                BEADS_LOADED,
                COPIES_PER_BEAD,          
                FC_ID,                    
                INCOMING_DNA_NAME,        
                KEY_PASS_WELLS,           
                ri.library_id, --r.LIBRARY_ID,               
                lib.full_name library_name, -- r.LIBRARY_NAME,             
                PAIRED_END,               
                PREDICTED_RECOVERY_BEADS, 
                r.REGION_ID,                
                REGION_NUMBER,            
                RESEARCH_PROJECT,         
                RUN_NAME,                 
                lib.SAMPLE_ID,                
                s.full_name SAMPLE_NAME,              
                SAMPLE_SET,               
                SS_ID,                    
                SUPERNATANT_BEADS,        
                TOTAL_KEY_PASS,           
                TOTAL_RAW_WELLS,
                NUM_BASES,
                NUM_READS,
                INDEX_SEQUENCE
            from run_region_454 r
            join region_index_454 ri on ri.region_id = r.region_id
            join library_summary lib on lib.library_id = ri.library_id
            join organism_sample s on s.organism_sample_id = lib.sample_id
        ) x454_detail
EOS
    ,
    has_constant => [
        sequencing_platform => { value => '454' },
    ],    
    has_optional => [
        _fasta_file => {
                        is => 'String',
                        is_transient => 1,
                        is_mutable => 1,
                  },
        _qual_file => {
                       is => 'String',
                       is_transient => 1,
                       is_mutable => 1,
                   },
        #< Run Region 454 from DW Attrs >#
        run_region_454     => {
            doc => '454 Run Region from LIMS.',
            is => 'GSC::RunRegion454',
            calculate => q| GSC::RunRegion454->get($region_id); |,
            calculate_from => ['region_id']
        },
        region_index_454     => {
            doc => 'Region Index 454 from LIMS.',
            is => 'GSC::RegionIndex454',
            calculate => q| GSC::RegionIndex454->get($id); |,
            calculate_from => ['id']
        },
        region_id           => { },
        region_number       => { },
        total_reads         => { column_name => "NUM_READS" },
        total_bases_read    => { column_name => "NUM_BASES" },
        is_paired_end       => { column_name => "PAIRED_END" },
        index_sequence      => { },

        sample_source       => { via => 'sample', to => 'source', },
        sample_source_name  => { via => 'sample_source', to => 'name' },

        # indirect via the sample source, but we let the sample manage that
        # since we sometimes don't know the source, it also tracks taxon directly
        taxon               => { via => 'sample', to => 'taxon', is => 'Genome::Site::TGI::Taxon' },
        species_name        => { via => 'taxon' },
    ],
};

sub full_path {
    Carp::confess("Full path is not valid for 454 instrument data");
}

sub _default_full_path {
    my $self = shift;
    return sprintf('%s/%s/%s', $self->_data_base_path, $self->run_name, $self->region_id);
}

sub calculate_alignment_estimated_kb_usage {
    my $self = shift;
    return 500000;
}

sub is_external {
    return;
}

#< Fastq, Fasta, Qual ... >#
sub dump_sanger_fastq_files {
    my $self = shift;
    
    my %params = @_;
    
    unless (-s $self->sff_file) {
        $self->error_message(sprintf("SFF file (%s) doesn't exist for 454 instrument data %s", $self->sff_file, $self->id));
        die $self->error_message;
    }
    
    my $dump_directory = delete $params{'directory'} || Genome::Sys->base_temp_directory();
    
    my $output_file = sprintf("%s/%s-output.fastq", $dump_directory, $self->id);
    
    my $cmd = Genome::Model::Tools::454::Sff2Fastq->create(sff_file => $self->sff_file,
                                                           fastq_file => $output_file);
    
    unless ($cmd->execute) {
        $self->error_message("Sff2Fastq failed while dumping fastq file for instrument data " . $self->id);
        die $self->error_message;
    }
    
    unless (-s $output_file) {
        $self->error_message("Sff2Fastq claims it worked, but the output file was gone or empty length while dumping fastq file for instrument data "
                             . $self->id . " expected output file was $output_file");
        die $self->error_message;
    }
    
    return ($output_file);
}

sub dump_fasta_file {
    my ($self, %params) = @_;
    $params{type} = 'fasta';
    return $self->_run_sffinfo(%params);
}

sub dump_qual_file {
    my ($self, %params) = @_;
    $params{type} = 'qual';
    return $self->_run_sffinfo(%params);
}

sub _run_sffinfo {
    my ($self, %params) = @_;

    # Type 
    my $type = delete $params{type};
    my %types_params = (
        fasta => '-s',
        qual => '-q',
    );
    unless ( defined $type and grep { $type eq $_ } keys %types_params ) { # should not happen
        Carp::confess("No or invalid type (".($type || '').") to run sff info.");
    }

    # Verify 64 bit
    unless ( Genome::Config->arch_os =~ /x86_64/ ) {
        Carp::confess(
            $self->error_message('Dumping $type file must be run on 64 bit machine.')
        );
    }
    
    # SFF
    my $sff_file = $self->sff_file;
    unless ( -s $sff_file ) {
        Carp::confess(
            $self->error_message(
                "SFF file ($sff_file) doesn't exist for 454 instrument data (".$self->id.")"
            )
        );
    }

    # File
    my $directory = delete $params{'directory'} 
        || Genome::Sys->base_temp_directory();
    my $file = sprintf("%s/%s.%s", $directory, $self->id, $type);
    unlink $file if -e $file;
    
    # SFF Info
    my $sffinfo = Genome::Model::Tools::454::Sffinfo->create(
        sff_file => $sff_file,
        output_file => $file,
        params => $types_params{$type},
    );
    unless ( $sffinfo ) {
        Carp::confess(
            $self->error_message("Can't create SFF Info command.")
        );
    }
    unless ( $sffinfo->execute ) {
        Carp::confess(
            $self->error_message("SFF Info command failed to dump $type file for instrument data (".$self->id.")")
        );
    }

    # Verify
    unless ( -s $file ) {
        Carp::confess(
            $self->error_message("SFF info executed, but a fasta was not produced for instrument data (".$self->id.")")
        );
    }

    return $file;
}

#< Run Info >#
sub run_identifier {
my $self = shift;

my $ar_454 = $self->run_region_454->get_analysis_run_454;

my $pse = GSC::PSE->get($ar_454->pse_id);
my $loadpse = $pse->get_load_pse;
my $barcode = $loadpse->picotiter_plate;

return $barcode->barcode->barcode;
}

sub run_start_date_formatted {
    my $self = shift;

    my ($y, $m, $d) = $self->run_name =~ m/R_(\d{4})_(\d{2})_(\d{2})/;

    my $dt_format = UR::Time->config('datetime');
    UR::Time->config(datetime=>'%Y-%m-%d');
    my $dt = UR::Time->numbers_to_datetime(0, 0, 0, $d, $m, $y);
    UR::Time->config(datetime=>$dt_format);

    return $dt; 
}
#<>#

1;

