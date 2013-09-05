package Genome::Model::Tools::Newbler;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use Carp 'confess';

use Bio::SeqIO;

class Genome::Model::Tools::Newbler {
    is => 'Command',
    has => [],
};

sub help_detail {
    return <<EOS
    Tools to work with newbler assembler
EOS
}

sub path_to_version_run_assembly {
    my $self = shift;

    my $app_bin;
    for my $bin( $self->possible_app_bin_names ) {
        $app_bin = $ENV{GENOME_SW} . '/454/' . $self->version . "/$bin";
        last if -d $app_bin;
    }

    my $assembler = $app_bin.'/runAssembly';
    unless ( -x $assembler ) {
        $self->error_message( "Invalid version: ".$self->version.' or versions runAssembly is not executable' );
        return;
    }
    return $assembler;
}

sub possible_app_bin_names {
    return (qw/
bin
applicationBin
applicationsBin
/);
}

#< input fastq files >#
sub input_fastq_files {
    my $self = shift;
    my @files = glob( $self->assembly_directory."/*-input.fastq" );
    unless ( @files ) {
        Carp::confess(
            $self->error_message( "No input fastq files found for assembly")
        ); #shouldn't happen but ..
    }
    return @files;
}

#< newbler output files >#
sub newb_ace_file {
    return $_[0]->assembly_directory.'/consed/edit_dir/454Contigs.ace.1';
}

#TODO - rename these wit newb*
sub scaffolds_agp_file {
    return $_[0]->assembly_directory.'/454Scaffolds.txt';
}

sub all_contigs_fasta_file {
    return $_[0]->assembly_directory.'/454AllContigs.fna';
}

sub all_contigs_qual_file {
    return $_[0]->assembly_directory.'/454AllContigs.qual';
}

sub newb_read_status_file {
    return $_[0]->assembly_directory.'/454ReadStatus.txt';
}

sub newb_metrics_file {
    return $_[0]->assembly_directory.'/454NewblerMetrics.txt';
}

#< post assemble output files/dirs >#
sub create_consed_edit_dir {
    my $self = shift;

    my $dir = $self->consed_edit_dir;
    unless ( -d $dir ) {
        Genome::Sys->create_directory($dir);
    }

    return 1;
}

sub consed_edit_dir {
    return $_[0]->assembly_directory.'/consed/edit_dir';
}

sub pcap_scaffold_ace_file {
    return $_[0]->consed_edit_dir.'/Pcap.454Contigs.ace';
}

sub contigs_bases_file {
    return $_[0]->consed_edit_dir.'/contigs.bases';
}

sub contigs_quals_file {
    return $_[0]->consed_edit_dir.'/contigs.quals';
}

sub gap_file {
    return $_[0]->consed_edit_dir.'/gap.txt';
}

sub read_info_file {
    return $_[0]->consed_edit_dir.'/readinfo.txt';
}

sub reads_placed_file {
    return $_[0]->consed_edit_dir.'/reads.placed';
}

sub reads_unplaced_file {
    return $_[0]->consed_edit_dir.'/reads.unplaced';
}

sub reads_unplaced_fasta_file {
    return $_[0]->consed_edit_dir.'/reads.unplaced.fasta';
}

sub supercontigs_fasta_file {
    return $_[0]->consed_edit_dir.'/supercontigs.fasta';
}

sub supercontigs_agp_file {
    return $_[0]->consed_edit_dir.'/supercontigs.agp';
}

sub stats_file {
    return $_[0]->consed_edit_dir.'/stats.txt';
}

#< create assembly sub dirs >#
sub create_consed_dir {
    my $self = shift;

    unless ( -d $self->assembly_directory.'/consed' ) {
        Genome::Sys->create_directory( $self->assembly_directory.'/consed' );
    }
    for my $subdir ( qw/ edit_dir phd_dir chromat_dir phdball_dir / ) {
        unless ( -d $self->assembly_directory."/consed/$subdir" ) {
            Genome::Sys->create_directory( $self->assembly_directory."/consed/$subdir" );
        }
    }
    return 1;
}

#filter out min_contig length
sub get_scaffolding_info { #TODO - reaname this get_valid_scaffolds
    my $self = shift;

    if ( -s $self->scaffolds_agp_file ) {
        return $self->create_scaffolded_contig_info;
    }
    else {
        return $self->create_unscaffolded_contig_info;
    }
    return;
}

#< create scaffolds info >#
sub create_scaffolded_contig_info {
    my $self = shift;

    #create hash of contig info
    my $scaffolds = {};
    my $fh = Genome::Sys->open_file_for_reading( $self->scaffolds_agp_file );
    while ( my $line = $fh->getline ) {
        my @tmp = split( /\s+/, $line );
        #contig describing line
        if ( $tmp[5] =~ /contig\d+/ ) {
            $scaffolds->{$tmp[5]}->{contig_length} = $tmp[7];
            $scaffolds->{$tmp[5]}->{supercontig} = $tmp[0];
            $scaffolds->{$tmp[5]}->{contig_name} = $tmp[5];
            $self->{prev_contig} = $tmp[5];
        }
        #gap describing line .. does not always exist
        if ( $tmp[6] =~ /fragment/ ) { #gap describing line
            my $prev_contig = $self->{prev_contig};
            $scaffolds->{$prev_contig}->{gap_length} = $tmp[5];
        }
    }
    $fh->close;

    #fill in missing gap sizes with default where gap describing line was missing
    my $default_gap_size = ( $self->can('default_gap_size') ) ? $self->default_gap_size : 1;
    for my $contig ( keys %$scaffolds ) {
        $scaffolds->{$contig}->{gap_length} = $default_gap_size
            unless exists $scaffolds->{$contig}->{gap_length};
    }

    #remove contigs less than min_contig_length & update gap size
    for my $contig ( sort keys %$scaffolds ) {
        my $current_scaffold = $scaffolds->{$contig}->{supercontig};
        my $gap_length = $scaffolds->{$contig}->{gap_length};
        my $contig_length = $scaffolds->{$contig}->{contig_length};

        $self->{PREV_SCAFFOLD} = $scaffolds->{$contig} unless defined $self->{PREV_SCAFFOLD};

        if ( $contig_length < $self->min_contig_length ) {
            my $add_to_prev_gap_size = $contig_length + $gap_length;
            $self->{PREV_SCAFFOLD}->{gap_length} += $add_to_prev_gap_size;
            delete $scaffolds->{$contig};            
        }
        else {
            $self->{PREV_SCAFFOLD} = $scaffolds->{$contig};
        }
    }

    #rename remaining contigs to pcap format
    my $pcap_contig = 1;
    my $pcap_supercontig = 0;
    for my $contig ( sort keys %$scaffolds ) {
        my $scaffold_name = $scaffolds->{$contig}->{supercontig};
        $self->{CURR_SCAFFOLD_NAME} = $scaffold_name unless $self->{CURR_SCAFFOLD_NAME};
        if ( not $self->{CURR_SCAFFOLD_NAME} eq $scaffold_name ) {
            $pcap_supercontig++;
            $pcap_contig = 1;
        }
        my $pcap_name = 'Contig'.$pcap_supercontig.'.'.$pcap_contig;
        $scaffolds->{$contig}->{pcap_name} = $pcap_name;
        $pcap_contig++;
        $self->{CURR_SCAFFOLD_NAME} = $scaffold_name;
    }

    #clean up
    delete $self->{CURR_SCAFFOLD_NAME};
    delete $self->{supercontig};
    delete $self->{pcap_name};
    delete $self->{PREV_SCAFFOLD};

    return $scaffolds;
}

sub create_unscaffolded_contig_info {
    my $self = shift;

    my $contigs = {};

    unless( -s $self->all_contigs_fasta_file ) {
        $self->error_message("Failed to find newbler all contigs file: ".$self->all_contigs_fasta_file);
        return;
    }

    my $newb_supercontig = 0;
    my $pcap_supercontig = 0;
    my $io = Bio::SeqIO->new( -format => 'fasta', -file => $self->all_contigs_fasta_file );
    while ( my $seq = $io->next_seq ) {
        next unless length $seq->seq >= $self->min_contig_length;
        my $length = length $seq->seq;
        my $newbler_supercontig_name = $self->_derive_newbler_supercontig_name( ++$newb_supercontig );
        $contigs->{$seq->primary_id}->{supercontig} = $newbler_supercontig_name;
        $contigs->{$seq->primary_id}->{contig_length} = length $seq->seq;
        $contigs->{$seq->primary_id}->{pcap_name} = 'Contig'.$pcap_supercontig++.'.1';
        $contigs->{$seq->primary_id}->{contig_name} = $seq->primary_id;
        $contigs->{$seq->primary_id}->{gap_length} = $self->default_gap_size;
    }

    return $contigs;
}

sub _derive_newbler_supercontig_name {
    my ( $self, $number ) = @_;

    if ( length $number == 1 ) {
        return 'scaffold0000'.$number;
    } elsif ( length $number == 2 ) {
        return 'scaffold000'.$number;
    } elsif ( length $number == 3 ) {
        return 'scaffold00'.$number;
    } elsif ( length $number == 4 ) {
        return 'scaffold0'.$number;
    } else {
        return 'scaffold'.$number;
    }
    return;
}

1;
