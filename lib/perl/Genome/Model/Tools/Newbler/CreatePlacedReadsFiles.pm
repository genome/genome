package Genome::Model::Tools::Newbler::CreatePlacedReadsFiles;

use strict;
use warnings;

use Genome;
use Bio::SeqIO;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::Newbler::CreatePlacedReadsFiles {
    is => 'Genome::Model::Tools::Newbler',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Newbler assembly directory',
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to consider',
        },
        default_gap_size => {
            is => 'Number',
            doc => 'Gap size to assign when newbler does not assign one',
            is_optional => 1,
            default_value => 10,
        }
    ],
};

sub help_brief {
    'Tool to create pcap style reads placed and unplaced files for newbler assemblies';
}

sub help_detail {
    return <<"EOS"
gmt newbler create-placed-reads-files --assembly-directory /gscmnt/111/newbler_assembly --min-contig-length 200
gmt newbler create-placed-reads-files --assembly-directory /gscmnt/111/newbler_assembly --min-contig-length 200 --default-gap-size 25
EOS
}

sub execute {
    my $self = shift;

    unless ( -d $self->consed_edit_dir ) {
        $self->create_consed_dir;
    }

    my $positions;
    unless( $positions = $self->_contigs_supercontig_positions ) {
        $self->error_message("Failed to determine contig supercontig positions");
        return;
    }

    unless( $self->_write_placed_reads_files( $positions )) {
        $self->error_message( "Failed to write reads.placed and readinfo.txt files" );
        return;
    }

    return 1;
}

sub _write_placed_reads_files {
    my ( $self, $positions ) = @_;

    my $in_fh = Genome::Sys->open_file_for_reading( $self->newb_read_status_file );
    unlink $self->read_info_file;
    $in_fh->getline; #header
    my $ri_fh = Genome::Sys->open_file_for_writing( $self->read_info_file );
    unlink $self->reads_placed_file;
    my $rp_fh = Genome::Sys->open_file_for_writing( $self->reads_placed_file );
    while ( my $line = $in_fh->getline ) {
        my @tmp = split( /\s+/, $line );
        #$tmp[0] = read name
        #$tmp[1] = status, assembled or not
        #$tmp[2] = contig name
        #$tmp[3] = 3' position
        #$tmp[4] = 3' complimented or uncomp
        #$tmp[5] = contig name
        #$tmp[6] = 5' position
        #$tmp[7] = 5' complimented or uncomp
        next unless $tmp[1] eq 'Assembled' or $tmp[1] eq 'PartiallyAssembled'; #otherwise reads has not been assembled
        next unless exists $positions->{$tmp[2]}; #otherwise contig was filtered out
        next unless $tmp[2] eq $tmp[5]; #read is split shouldn't happen

        #print readinfo.txt file info
        my $u_or_c = ( $tmp[4] eq '+' ) ? 'U' : 'C';
        my $read_contig_start_pos = ( $tmp[3] > $tmp[6] ) ? $tmp[6] : $tmp[3];
        my $read_length = abs( $tmp[3] - $tmp[6]) + 1;
        my $pcap_contig_name = $positions->{$tmp[2]}->{pcap_contig_name};

        my $ri_string = "$tmp[0] $pcap_contig_name $u_or_c $read_contig_start_pos $read_length\n";
        $ri_fh->print( $ri_string );

        #print reads.placed file info
        $u_or_c = ( $u_or_c eq 'U' ) ? 1 : 0;
        my $contig_supercontig_start_pos = $positions->{$tmp[2]}->{supercontig_start_position};
        my $read_supercontig_start_pos = $contig_supercontig_start_pos + $read_contig_start_pos - 1;
        my ($pcap_supercontig_name) = $pcap_contig_name =~ /(Contig\d+)\.\d+/;

        my $rp_string = "$tmp[0] 1 $read_length $u_or_c $pcap_contig_name $pcap_supercontig_name $read_contig_start_pos $read_supercontig_start_pos\n";
        $rp_fh->print( $rp_string );
    }
    $in_fh->close;
    $ri_fh->close;
    $rp_fh->close;
}

sub _contigs_supercontig_positions {
    my $self = shift;
    #filter contigs my min length
    my $scafs;
    unless ( $scafs = $self->get_scaffolding_info ) {
        $self->error_message("Failed to get scaffolding info");
        return;
    }
    my $pos = {};
    my $current_supercontig;
    my $supercontig_position = 1;
    #create a hash of contig's supercontig position
    for my $ctg ( sort keys %$scafs ) {
        my $pcap_contig_name = $scafs->{$ctg}->{pcap_name};
        $pos->{$ctg}->{pcap_contig_name} = $pcap_contig_name;
        my ( $pcap_supercontig_name ) = $pcap_contig_name =~ /(Contig\d+)\.\d+/;
        if ( not $current_supercontig ) {
            $current_supercontig = $pcap_supercontig_name;
        }
        #re-set supercontig position if start of new scaffold
        if ( $current_supercontig ne $pcap_supercontig_name ) {
            $supercontig_position = 1;
            $current_supercontig = $pcap_supercontig_name;
        }
        $pos->{$ctg}->{supercontig_start_position} = $supercontig_position;
        $supercontig_position += $scafs->{$ctg}->{contig_length};
        $supercontig_position += $scafs->{$ctg}->{gap_length};
    }
    
    return $pos;
}

1;
