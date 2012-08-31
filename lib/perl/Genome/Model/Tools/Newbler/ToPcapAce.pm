package Genome::Model::Tools::Newbler::ToPcapAce;

use strict;
use warnings;

use Genome;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::Newbler::ToPcapAce {
    is => 'Genome::Model::Tools::Newbler',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Newbler assembly directory',
        },   
        default_gap_size => {
            is => 'Number',
            doc => 'Gap size to assign when newbler does not assign one',
            is_optional => 1,
            default_value => 10,
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to include in post',
        },
    ],
};

sub help_brief {
    'Tool to convert newber generated ace file to pcap scaffolded ace file';
}

sub help_detail {
    return <<"EOS"
gmt newbler to-pcap-ace --assembly-directory /gscmnt/111/assembly/newbler_e_coli
gmt newbler to-pcap-ace --assembly-directory /gscmnt/111/assembly/newbler_e_coli --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    if ( not $self->_validate_inputs ) {
        $self->error_message( "Failed to validate inputs" );
    }

    my $scaffolds;
    unless ( $scaffolds = $self->get_scaffolding_info ) {
        $self->error_message( "Failed to parse newbler scaffold file" );
        return;
    }

    unless( $self->_write_pcap_style_ace( $scaffolds ) ) {
        $self->error_message( "Failed to write new scaffolded pcap ace file" );
        return;
    }

    unless ( $self->_write_gap_txt_file( $scaffolds ) ) {
        $self->error_message( "Failed to create gap.txt file" );
        return;
    }

    return 1;
}

sub _write_gap_txt_file {
    my ( $self, $scaffolds ) = @_;
    my $gap_file = $self->gap_file;
    unlink $gap_file;
    my $fh = Genome::Sys->open_file_for_writing( $gap_file );
    for my $newb_contig_name ( sort keys %{$scaffolds} ) {
        my $pcap_contig_name = $scaffolds->{$newb_contig_name}->{pcap_name};
        my $gap_size = ( $scaffolds->{$newb_contig_name}->{gap_size} ) ?
            $scaffolds->{$newb_contig_name}->{gap_size} : $self->default_gap_size;
        $fh->print( $pcap_contig_name.' '.$gap_size."\n" );
    }
    $fh->close;

    return 1;
}

sub _validate_inputs {
    my $self = shift;

    #assembly directory
    if ( not -d $self->assembly_directory ) {
        $self->error_message( "Can not find assembly directory: ".$self->assembly_directory );
        return;
    }
    #input newbler ace file
    if ( not -s $self->newb_ace_file ) {
        $self->error_message( "Failed to find newbler ace file or file is zero size: ".$self->newbler_ace_file );
        return;
    }
    #output pcap style ace file
    if ( -s $self->pcap_scaffold_ace_file ) {
        $self->status_message( "Removing existing pcap ace out file: ".$self->pcap_scaffold_ace_file );
    }

    return 1;
}

sub _write_pcap_style_ace {
    my ( $self, $scaffolds ) = @_;

    $self->status_message( "Found scaffold agp file or file was supplied, generating scaffolded pcap ace file" );

    my $fh = Genome::Sys->open_file_for_reading( $self->newb_ace_file );
    unlink $self->pcap_scaffold_ace_file.'.int';
    my $fh_int = Genome::Sys->open_file_for_writing( $self->pcap_scaffold_ace_file.'.int' );
    my $print_setting = 0;
    my $contig_count = 0;
    my $read_count = 0;
    while ( my $line = $fh->getline ) {
        if ( $line =~ /^CO\s+/ ) {
            my ( $newb_ctg_name ) = $line =~ /^CO\s+(\S+)/;
            my $rest_of_line = "$'";
            
            if ( exists $scaffolds->{$newb_ctg_name} ) {
                my $pcap_name = $scaffolds->{$newb_ctg_name}->{pcap_name};
                $fh_int->print( "CO $pcap_name $rest_of_line" );
                $contig_count++;
                $print_setting = 1;
                next;
            }
            else {
                $print_setting = 0;
            }
        }
        $fh_int->print( $line ) if $print_setting == 1;
        $read_count++ if $line =~ /^RD\s+/;
    }
    $fh->close;
    $fh_int->close;
    
    #append AS line that describes number of contigs and reads to new ace
    my $as_line =  "AS  $contig_count $read_count\n\n";
    unlink $self->pcap_scaffold_ace_file;
    my $fh_out = Genome::Sys->open_file_for_writing( $self->pcap_scaffold_ace_file );
    $fh_out->print( $as_line );
    my $fh_in = Genome::Sys->open_file_for_reading( $self->pcap_scaffold_ace_file.'.int' );
    while ( my $line = $fh_in->getline ) {
        $fh_out->print( $line );
    }
    $fh_out->close;
    $fh_in->close;

    unlink $self->pcap_scaffold_ace_file.'.int';

    return 1;
}

1;
