package Genome::InstrumentData::AlignmentResult::Nucmer;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::InstrumentData::AlignmentResult::Nucmer{
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'nucmer', is_param => 1 },
    ],
};

sub required_arch_os { 'x86_64' }

sub required_rusage {
    # TODO - not sure yet how much mem to require
    return "-q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT} -R 'select[type==LINUX64 && mem>4000] rusage[mem=4000] span[hosts=1]' -M 4000000";
}

sub _run_aligner {
    my $self = shift;
    my @input_path_names = @_;

    # merge fastq files into single fasta file
    my $input_fasta = $self->temp_scratch_directory.'/input.fasta';
    my $writer = Genome::Model::Tools::Sx::PhredWriter->create(
        file => $input_fasta,
    );
    for my $file ( @input_path_names ) {
        # TODO - read in as paired if paired
        my $reader = Genome::Model::Tools::Sx::FastqReader->create(
            file => $file,
        );
        while ( my $seq = $reader->read ) {
            my $id = $seq->{id};
            $id =~ s/\/\d$//;
            $seq->{id} = $id;
            $writer->write( $seq );
        }
    }

    # reference to align against
    my $reference_fasta = $self->reference_build->full_consensus_path('fa');

    # nucmer is subpackage of mummer
    my $mummer = Genome::Model::Tools::Mummer->create( use_version => $self->aligner_version );
    my $version_nucmer = $mummer->path_for_version.'/nucmer';

    # RUN ALIGNER
    # USAGE: nucmer  [options]  <Reference>  <Query>
    my $aligner_params = $self->aligner_params || '';
    $aligner_params .= ' --prefix '.$self->temp_scratch_directory.'/OUT';

    my $align_cmd = $version_nucmer.' '.$aligner_params.' '.$reference_fasta.' '.$input_fasta;
    $self->status_message( "Align command: $align_cmd\n" );
    Genome::Sys->shellcmd( cmd => $align_cmd );

    # check align output file
    my $align_out_file = $self->temp_scratch_directory.'/OUT.delta';
    if ( not -s $align_out_file ) {
        $self->error_message("NO nucmer align out file found: $align_out_file");
        return;
    }

    # convert output to sam format
    my $sam_file = $self->temp_scratch_directory.'/all_sequences.sam';
    if ( not -s $sam_file ) {
        $self->status_message("Failed to find sam output file: $sam_file");
        return;
    }
    my $sam_fh = IO::File->new(">> $sam_file"); 

    my $sc_fh = Genome::Sys->open_file_for_reading( $align_out_file );
    my %align_info;
    while ( my $line = $sc_fh->getline ) {
        my @tmp = split( /\s+/, $line );
        chomp $line;
        if ( $line =~ /^>/ ) {
            # > ref_name query_name ref_length query_length
            $line =~ s/^>//;
            my @tmp = split(/\s+/, $line);
            $align_info{reference_name} = $tmp[0];
            $align_info{query_name} = $tmp[1];
            
            # next line contains align coordinates
            # 1 429 15413 15841 0 0 0
            my $coordinates = $sc_fh->getline;
            chomp $coordinates;
            @tmp = split( /\s+/, $coordinates );
            if ( not scalar @tmp == 7 ) {
                $self->error_message("Expected line with 7 lines but got: $coordinates\n");
                return;
            }
            # validate line??
            # [0] = ref start
            # [1] = ref end
            # [2] = query start
            # [3] = query end
            $align_info{reference_align_length} = abs( $tmp[1] - $tmp[0] );
            $align_info{query_align_length} = abs( $tmp[3] - $tmp[2] );
            $align_info{query_reverse_aligned} = ( $tmp[2] > $tmp[3] ) ? 1 : 0;
            $align_info{query_left_align_pos} = ( $tmp[2] > $tmp[3] ) ? $tmp[3] : $tmp[2];
        }
        elsif ( $line eq '0' ) {
            # 0 - end of alignment description
            my $sam_string = $self->_nucmer_to_sam_format(%align_info);
            $sam_fh->print ( "$sam_string\n" );
            delete $align_info{indels};
        }
        elsif ( $line =~ /^-\d+$/ or $line =~ /^\d+$/ ) { # $line ne '0'
            # neg number = deletion if ref seq
            # number = insertion in ref sequence
            push @{$align_info{indels}}, $line;
        }
    }

    $sam_fh->close;

    return 1;
}

sub postprocess_bam_file {
    my $self = shift;
    # currently not 1:1 input/output ratio so _verify_bam will fail
    # create bam flagstat
    $self->status_message('Creating flagstat for bam file');
    unless( $self->_create_bam_flagstat ) {
        $self->error_message('Failed to create bam flagstat file');
        die $self->error_message;
    }
    $self->status_message('Indexing bam file');
    unless($self->_create_bam_index) {
        $self->error_message('Fail to create bam md5');
        die $self->error_message;
    }


    return 1;
}

sub aligner_params_for_sam_header {
    # TODO - not sure what this is
    return 'Nucmer test';
}

sub fillmd_for_sam {
    # TODO - not sure about this yet
    return 0;
}

sub _nucmer_to_sam_format {
    my ( $self, %h ) = @_;

    my @saminfo;
    # SAM out
    # [0]  = query name
    # [1]  = bitwise flag
    # [2]  = target contig, * if not aligned to anyting
    # [3]  = read align position, 0 if not aligned
    # [4]  = MAPping qual
    # [5]  = CIGAR string
    # [6]  = mate ref seq '=' if both mates hit the same thing, 0 if no mate
    # [7]  = mate align position, 0 if not aligned
    # [8]  = inferred insert size ( + or - value depending on 1st or 2nd read )
    # [9]  = sequence
    # [10] = qual

    # cigar string
    my $cigar_string = $h{query_align_length}.'M';
    if ( exists $h{indels} ) {
        $cigar_string = $self->_cigar_string(\@{$h{indels}}, $h{query_align_length});
    }

    # flag
    # 8  = read align forward, mate not mapped
    # 24 = read align reverse, mate not mapped
    my $flag = ( $h{query_reverse_aligned} == 1 ) ? 8 : 24;

    $saminfo[0] = $h{query_name};
    $saminfo[1] = $flag;                       # TODO
    $saminfo[2] = $h{reference_name};
    $saminfo[3] = $h{query_left_align_pos};
    $saminfo[4] = 0;
    $saminfo[5] = $cigar_string;
    $saminfo[6] = '*'; # no mate
    $saminfo[7] = 0;   # no mate
    $saminfo[8] = 0;   # no mate
    $saminfo[9] = '*';
    $saminfo[10] = '*';
    
    return join( "\t", map {$_} @saminfo );
}

sub _cigar_string {
    my $self = shift;
    my $indels = shift;
    my $query_align_length = shift;

    #store each Ins or Del operations
    my %h;
    my $op = 1;
    my $query_pos_offset = 0;
    foreach ( @$indels ) {
        my $type = ( $_ > 0 ) ? 'D' : 'I';
        # I = insertion in query
        # D = deletion in query
        # track position after each op
        $query_pos_offset += abs( $_ );
        # if insertion remove ++'ed offset
        $query_pos_offset -- if $type eq 'D';
        if ( abs($_) == 1 ) {
            # make 1D1D or 1I1I to 2D or 2I
            $h{$op - 1}{base_count}++;
        }
        else {
            $h{$op}{base_count}++;
            $h{$op}{type} = $type;
            $h{$op}{base_offset} = abs( $_ );
            ++$op;
        }
    }
    # write cigar string left to right
    my $cigar_string;
    foreach ( sort {$a<=>$b} keys %h ) {
        # add matching operation to left of indel
        $cigar_string .= ( $h{$_}{base_offset} - 1 ) . 'M';
        # append indel operation
        $cigar_string .= $h{$_}{base_count} . $h{$_}{type};
    }

    # append remaing matches to the right end
    my $end_match = $query_align_length - $query_pos_offset;
    if ( not $end_match > 0 ) {
        $self->error_message('Query offset exceeded alignment length');
        # cigar string not available
        return '*';
    }

    $cigar_string .= $end_match . 'M';
    
    return $cigar_string;
}

1;
