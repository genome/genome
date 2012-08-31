package Genome::Model::Tools::Bed::Convert::Indel::PindelToBed;

use warnings;
use strict;

use Genome;
use Workflow;

class Genome::Model::Tools::Bed::Convert::Indel::PindelToBed {
    is => ['Genome::Model::Tools::Bed::Convert::Indel'],
    has => [
        use_old_pindel => {
            type => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => 'Run on pindel 0.2 or 0.1',
        },
        include_normal => {
            type => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Include events which have some or all normal support alongside events with only tumor support',
        },
        include_bp_ranges => {
            type => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Include pindels calculated bp_range for the location of the indel.',
        },
    ],
    has_transient_optional => [
        _big_output_fh => {
            is => 'IO::File',
            doc => 'Filehandle for the output of large events',
        },
        big_output_file => {
            is => 'String',
            doc => 'File containing large output',
        },
    ],
    has_param => [
        # Make workflow choose 64 bit blades, this is needed for samtools faidx
        lsf_resource => {
            default_value => 'rusage[mem=4000] select[type==LINUX64] span[hosts=1] -M 4000000',
        },
        lsf_queue => {
            default_value => 'long'
        }, 
    ],
};

sub help_brief {
    "Transforms a file from pindel output format to bed format";
}

sub help_synopsis {
    return <<"EOS"
gmt bed convert indel pindel-to-bam --input-file pindel.outfile --output-file pindel.adapted 
EOS
}

sub help_detail {                           
    return <<EOS 
Transforms a file from pindel output format to bed format
EOS
}

sub initialize_filehandles {
    my $self = shift;
    if($self->_big_output_fh) {
        return 1; #Already initialized
    }
    
    my $big_output = $self->output . ".big_deletions";
    $self->big_output_file($big_output); 

    return $self->SUPER::initialize_filehandles(@_);
}

sub close_filehandles {
    my $self = shift;

    # close big deletions fh
    my $big_output_fh = $self->_big_output_fh;
    close($big_output_fh) if $big_output_fh;

    return $self->SUPER::close_filehandles(@_);
}

sub process_source { 
    my $self = shift;
    $self->_input_fh->close;
    my $temp_bed = Genome::Sys->create_temp_file_path;
    my $ppr_cmd = Genome::Model::Tools::Pindel::ProcessPindelReads->create(
                    input_file => $self->source,
                    output_file => $temp_bed,
                    reference_build_id => $self->reference_build_id,
                    mode => 'to_bed', 
                    big_output_file => $self->big_output_file, );
    unless( $ppr_cmd->execute ){
        die $self->error_message( "Failed to run gmt pindel process-pindel-reads" );
    }
                   
    my $bed_fh = Genome::Sys->open_file_for_reading($temp_bed);
    while( my $line = $bed_fh->getline){
        chomp $line;
        $self->write_bed_line(split /\t/, $line);
    }
    $bed_fh->close;
    my $ifh = Genome::Sys->open_file_for_reading($self->source);
    $self->_input_fh($ifh);
    return 1;
}

1;
