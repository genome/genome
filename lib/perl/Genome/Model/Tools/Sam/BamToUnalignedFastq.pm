package Genome::Model::Tools::Sam::BamToUnalignedFastq;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sam::BamToUnalignedFastq {
    is => ['Genome::Model::Tools::Sam'],
    has_input => [
        bam_file => {
            is => 'Text',
            doc => 'A BAM format file of alignment data.'
        },
        output_directory => {
            is => 'Text',
            doc => 'A directory to output s_*_*_sequence.txt files.  Two files for unmapped pairs and one file for unmapped fragments or unmapped mates whose mate-pair is mapped.',
        },
        print_aligned => {
            is => 'Boolean',
            default => 0,
            doc => 'If set, this tool will print only the aligned reads, instead of the unaligned reads.',
        },
        ignore_bitflags => {
            is => 'Boolean',
            default => 0,
            doc => 'If set, this tool will ignore the bitflag column and instead use the RNAME column of the bam to determine if a read has been mapped',
        },
        filter_duplicates => {
            is => 'Boolean',
            default => 0,
            doc => 'If set, this tool will not print reads that are marked as duplicates in the bam file, this happens regardless of the ignore-bitflags flag',
        },
    ],
};

sub execute {
    my $self = shift;

    unless (-d $self->output_directory) {
        unless (Genome::Sys->create_directory($self->output_directory)) {
            die('Failed to create output directory '. $self->output_directory);
        }
    }

    #Parse through header and create filehandles for all lanes
    my $header_fh = IO::File->new('samtools view -H '.$self->bam_file.'|');
    my %fhs;
    while (my $line =$header_fh->getline) {
        if ($line =~ /^\@RG/){
            unless ($line =~ /^\@RG\s+ID:(\d+).+PU:\S+\.(\d+).+DS:(paired end|fragment)/) {
                die('Failed to match line '. $line .' with regex!');
            }
            my $instrument_data_id = $1;
            my $lane = $2;
            my $run_type = $3;
            my $instrument_data_directory = $self->output_directory .'/'. $instrument_data_id;
            unless (Genome::Sys->create_directory($instrument_data_directory)) {
                die('Failed to create instrument data directory '. $instrument_data_directory);
            }
            my $fragment_file = $instrument_data_directory .'/s_'. $lane .'_sequence.txt';
            my $fragment_fh = Genome::Sys->open_file_for_writing($fragment_file);
            unless ($fragment_fh) {
                die('Failed to create fragment filehandle '. $fragment_file);
            }
            $fhs{$instrument_data_id}{'fragment'} = $fragment_fh;
            if ($run_type eq 'paired end') {
                for my $end (1 .. 2) {
                    my $end_file = $instrument_data_directory .'/s_'. $lane .'_'.$end .'_sequence.txt';
                    my $end_fh = Genome::Sys->open_file_for_writing($end_file);
                    unless ($end_fh) {
                        die('Failed to create read end '. $end .' filehandle '. $end_file);
                    }
                    my $key = 'read_'. $end;
                    $fhs{$instrument_data_id}{$key} = $end_fh;
                }
            }
        }
    }
    
    my %read_pairs;
    my $bam_fh = IO::File->new('samtools view '.$self->bam_file.'|');
    if ($self->ignore_bitflags){
        #print everything out to a single file, sort, then print to correct fh
        $self->debug_message("extracting (un)aligned reads based on reference name");
        my $temp_out = Genome::Sys->create_temp_file_path();
        my $tfh = IO::File->new("> $temp_out");
        while ( my $align = $bam_fh->getline) {
            my @cols = split(/\s+/, $align);
            my $mapped;
            my $flag = $cols[1];
            my $rname = $cols[2];
            if ($rname eq '*'){
                $mapped = 0;
            }
            else{
                $mapped = 1;
            }

            if ($mapped == $self->print_aligned){
                unless ($self->filter_duplicates and ($flag & 1024)){
                    $tfh->print($align);
                }
            }
        }
        $tfh->close;
        my $temp_out_sorted = Genome::Sys->create_temp_file_path();
        $self->debug_message("sorting (un)aligned reads");
        system("sort $temp_out > $temp_out_sorted");
        my $sfh = IO::File->new($temp_out_sorted);
        my $cache;
        $self->debug_message("placing (un)aligned reads in paired-end or fragment fastq files");
        while (my $align = $sfh->getline){
            if ( ! $cache ){
                $cache = $align;
                next;
            }
            my ($cache_read) = split(/\s+/, $cache);
            my $cache_base = $cache_read;
            $cache_base =~ s/\/.*$//;
            my ($current_read) = split(/\s+/, $align);
            my $current_base = $current_read;
            $current_base =~ s/\/.*$//;
            if ($cache_read eq $current_read or $cache_base eq $current_base){
                print_align_to_fh($cache,\%fhs,'read_1');
                print_align_to_fh($align,\%fhs,'read_2');
                $cache = undef;
            }else{
                print_align_to_fh($cache,\%fhs,'fragment');
                $cache = $align;
            }
        }
        if ($cache){
            print_align_to_fh($cache, \%fhs, 'fragment');
        }
    }
    else{
        while (my $align = $bam_fh->getline) {
            my @cols = split(/\s+/, $align);
            my ($mapped, $type);
            my $flag = $cols[1];
            if ($flag & 1) { # Is this read part of a pair?
                if ($flag & 64)  {
                    $type = 'read_1';
                } elsif ($flag & 128) {
                    $type = 'read_2';
                } else {
                    die('Read pair info lost for alignment of read '. $align->qname .' from BAM file '. $self->bam_file);
                }

                # if both halves of the mate do not have the same mapping status, this is a fragment
                unless ( ($flag & 4) *2 == ($flag & 8) ){ 
                    $type = 'fragment';
                }

                if ($flag & 4) { # Is this part of the pair unmapped?
                    $mapped = 0;
                } else {
                    $mapped = 1;
                }
                # Else the read is not a part of a pair, and treated as a fragment
            } else {
                # Fragment Read
                $type = 'fragment';
                if ($flag & 4) { # Is this read unmapped?
                    $mapped = 0;
                } else {
                    $mapped = 1;
                }
            }

            # If the read is mapped and we're printing mapped things, or if the read is unmapped and we're printing unaligned things... print it
            if ( $mapped == $self->print_aligned ) {
                unless ($self->filter_duplicates and ($flag & 1024)){
                    print_align_to_fh($align,\%fhs,$type);
                }
            }
        }
    }
    return 1;
}


sub print_align_to_fh {
    my $align = shift;
    my $fhs = shift;
    my $type = shift;

    my @cols = split(/\s+/,$align);

    #print $align."\n";
    my ($instrument_data_id) = $align =~/RG:Z:(\d+)/;
    #print $instrument_data_id."\n";
    my $fh = $fhs->{$instrument_data_id}->{$type};
    unless ($fh) {
        die('Failed to get back fragment filehandle using instrument data id '. $instrument_data_id .' and read type '. $type);
    }
    my $name = $cols[0];
    # TODO: verify orientation with original fastq
    my $seq = $cols[9];
    # TODO: verify orientation and quality conversion with original fastq
    my $qual = $cols[10];
    my $flag = $cols[1];
    print $fh '@' ."$name\n$seq\n+\n$qual\n";
}


1;
