package Genome::Model::Tools::Assembly::Ace::RenameDuplicateReads;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Assembly::Ace::RenameDuplicateReads {
    is => 'Genome::Model::Tools::Assembly::Ace',
    has => [
        ace_in => {
            is => 'Text',
            doc => 'Input ace file',
        },
        ace_out => {
            is => 'Text',
            doc => 'Output ace file name',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Tool to rename duplicated reads in ace file .. duplicated reads are renamed with -1, -2, etc'
}

sub help_detail {
    return <<EOS
gmt assembly ace rename-duplicated-reads --ace-in /gscmnt/assembly/e_coli/contigs.ace.0
gmt assembly ace rename-duplicated-reads --ace-in /gscmnt/assembly/e_coli/contigs.ace.0 --ace-out contigs.ace.0.fixed
EOS
}

sub execute {
    my $self = shift;

    #check ace file
    unless( -s $self->ace_in ) {
        $self->error_message("Could not find ace input or input file is zero size: ".$self->ace_in );
        return;
    }

    #get reads duplicated in ace file
    if ( not $self->_duplicate_reads ) {
        $self->debug_message("Did not find any duplicate reads in ace file: ".$self->ace_in );
        return 1;
    }
 
    #report counts and read names
    $self->_print_status_message;

    #re-write ace file
    my $ace_out = ( $self->ace_out ) ? $self->ace_out : $self->ace_in.'.duplicate_reads_renamed';
    unlink $ace_out;
    my $fh_out = Genome::Sys->open_file_for_writing( $ace_out );
    my $fh = Genome::Sys->open_file_for_reading( $self->ace_in );

    #re-write each line that uses read name
    while ( my $line = $fh->getline ) {
        chomp $line;
        if ( $line =~ /^AF\s+/ or $line =~ /^BS\s+/ or $line =~ /^RD\s+/ or $line =~ /^DS\s+/ ) {
            my ($line_name) = $line =~ /^(\w\w)\s+/;
            my $method = '_update_'.$line_name.'_line';
            my $updated_line = $self->$method( $line );
            $fh_out->print( $updated_line );
        } else {
            $fh_out->print( $line."\n" );
        }
    }
    $fh_out->close;
    $fh->close;
    return 1;
}

sub _update_DS_line {
    my $self = shift;
    my $line = shift;

    my ( $read_name ) = $line =~ /(\S+)\.phd\./;
    my $right_side_of_line = "$'";
    $read_name =~ s/-\d+$//;
    my ( $chromat_name ) = $line =~ /CHROMAT_FILE:\s+(\S+)/; #chromat name and read names diff for newbler
    $chromat_name =~ s/-\d+$//;
    my $new_read_name = $read_name;
    my $new_chromat_name = $chromat_name;
    if ( exists $self->{_DUP_READS}{$read_name} ) {
        $new_read_name = $self->_check_and_update_dup_reads_status( $read_name );
        if ( not $chromat_name =~ /sff:/ ) {
            $new_chromat_name = ( $self->{_DUP_READS}{$read_name}{current_iteration} > 1 ) ?
                $new_chromat_name.'-'.( $self->{_DUP_READS}{$read_name}{current_iteration} - 1 ) :
                $new_chromat_name;
        }
    }

    return 'DS CHROMAT_FILE: '.$new_chromat_name.' PHD_FILE: '.$new_read_name.'.phd.'.$right_side_of_line."\n";
}

sub _update_RD_line {
    my $self = shift;
    my $line = shift;

    my ( $read_name) = $line =~ /RD\s+(\S+)/;
    my $right_side_of_line = "$'";
    $read_name =~ s/-\d+$//;
    my $new_read_name = ( exists $self->{_DUP_READS}{$read_name} ) ?
        $self->_check_and_update_dup_reads_status( $read_name ) :
        $read_name;

    return "RD $new_read_name".$right_side_of_line."\n";
}

sub _update_AF_line {
    my $self = shift;
    my $line = shift;

    my ( $read_name ) = $line =~ /AF\s+(\S+)/;
    my $right_side_of_line = "$'";
    $read_name =~ s/-\d+$//; #clip off -# if read already duplicated
    my $new_read_name = $read_name;
    if ( exists $self->{_DUP_READS}{$read_name} ) {
        $self->{_DUP_READS}{$read_name}{current_iteration}++;
        $new_read_name = $self->_check_and_update_dup_reads_status( $read_name );
    }

    return "AF $new_read_name"."$right_side_of_line\n";
}

sub _update_BS_line {
    my $self = shift;
    my $line = shift;

    my ( $read_name ) = $line =~ /(\S+)$/;
    my $left_side_of_line = "$`";
    $read_name =~ s/-\d+$//;
    my $new_read_name = ( exists $self->{_DUP_READS}{$read_name} ) ?
        $self->_check_and_update_dup_reads_status( $read_name ) :
        $read_name;

    return $left_side_of_line.$new_read_name."\n";
}

sub _check_and_update_dup_reads_status {
    my $self = shift;
    my $read_name = shift;
    
    my $new_read_name = ( $self->{_DUP_READS}{$read_name}{current_iteration} > 1 ) ?
        $read_name.'-'.( $self->{_DUP_READS}{$read_name}{current_iteration} - 1 ) : $read_name;

    return $new_read_name;
}

sub _duplicate_reads {
    my $self = shift;

    my %duplicate_reads;
    my %all_reads;

    my $c = 0; my $tc = 0;
    my $fh = Genome::Sys->open_file_for_reading( $self->ace_in );

    my $contig_name;

    while ( my $line = $fh->getline ) {
        if ( $line =~ /^CO\s+/ ) {
            ( $contig_name ) = $line =~ /^CO\s+(\S+)/;
        }
        if ( my ( $read_name ) = $line =~ /^AF\s+(\S+)/ ) {
            if ( exists $all_reads{$read_name} ) {
                $self->{_DUP_READS}{$read_name}{times_used}++;
                $self->{_DUP_READS}{$read_name}{contig_name} = $contig_name;
            } else {
                $all_reads{$read_name} = 1;
            }
            #print status message
            ++$tc;
            if ( ++$c == 10000 ) {
                $self->debug_message("Checked $tc reads for duplicates");
                $c = 0;
            }
        }
    }

    $fh->close;

    return \%duplicate_reads;
}

sub _print_status_message {
    $_[0]->status_message( 'Found '. ( scalar keys %{$_[0]->{_DUP_READS}} )." duplicate reads\n" );
    $_[0]->status_message( map { "\t$_\n" } keys %{$_[0]->{_DUP_READS}} );
}

1;

