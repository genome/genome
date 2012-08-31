package Genome::Model::Tools::Assembly::SplitScaffold;

use strict;
use warnings;

use Genome;

use Cwd;
use Data::Dumper;
use Sort::Naturally;

class Genome::Model::Tools::Assembly::SplitScaffold {
    is => 'Command',
    has => [
        ace_file => {
            type => 'Text',
            doc => "This is the input ace file"
        }, 
	split_contigs => {
	    type => 'Text',
	    is_many => 1,
	    doc => 'Multiple scaffold contigs to split scaffolds by',
	},
	out_file_name => {
	    type => 'Text',
	    is_optional => 1,
	    doc => "This is the name of the output file",
	},
    ],
};

sub help_brief {
    'Tool to split pcap style ace file scaffolds';
}

sub help_synopsis { 
    return;
}

sub help_detail {
    return <<EOS 
    gmt split-scaffold --ace-file Pcap.Contigs.ace --split-contigs Contig0.8 --out-file-name Pcap.Contigs.ace.split
    gmt split-scaffold --ace-file Pcap.Contigs.ace --split-contigs Contig0.8,Contig2.5
EOS
}

sub execute {
    my $self = shift;

    unless ( $self->split_contigs ) {
	$self->error_message("You must supply contig(s) to split scaffolds by");
	return;
    }

    unless ( -s $self->ace_file ) {
	$self->error_message("Failed to find input ace file: ".$self->ace_file);
	return;
    }

    unless ( $self->_split_scaffolds ) {
	$self->error_message("Failed to successfully split scaffolds");
	return;
    }

    return 1;
}

sub _split_scaffolds {
    my $self = shift;

    my $rename_contigs = $self->_get_contigs_to_rename;

    my $fh = Genome::Sys->open_file_for_reading( $self->ace_file );

    my $ace_out_file = ( $self->out_file_name ) ? $self->out_file_name : $self->ace_file.'.scaffolds_split';

    unlink $ace_out_file if -e $ace_out_file;
    my $fh_out = Genome::Sys->open_file_for_writing( $ace_out_file );

    unlink $ace_out_file.'.LOG';# if -e $ace_out_file.'LOG';
    my $log_fh = Genome::Sys->open_file_for_writing( $ace_out_file.'.LOG' );

    while ( my $line = $fh->getline ) {
	if ( $line =~ /^CO\s+/ ) {
	    chomp $line;
	    my ( $contig_number ) = $line =~ /CO\s+Contig(\d+\.\d+)\s+/;
	    my $rest_of_line = "$'";

	    if ( exists $rename_contigs->{$contig_number} ) {
		$fh_out->print( 'CO Contig'. $rename_contigs->{$contig_number}.' '.$rest_of_line."\n" );
		$log_fh->print( 'Contig'.$contig_number.' split to Contig'.$rename_contigs->{$contig_number}."\n" );
		$self->status_message( 'Contig'.$contig_number.' split to Contig'.$rename_contigs->{$contig_number}."\n" );
	    }
	    else {
		$fh_out->print( $line."\n" );
	    }
	}
	else {
	    $fh_out->print( $line );
	}
    }

    $fh->close;
    $fh_out->close;
    $log_fh->close;

    return 1;
}

sub _get_contigs_to_rename {
    my $self = shift;

    my @contigs_to_split_at = $self->_get_scaffold_split_contig_names;

    my $split_points = {};
    for my $contig ( @contigs_to_split_at ) {
	my ( $scaffold_number, $contig_number ) = $contig =~ /Contig(\d+)\.(\d+)/i;
	$split_points->{$scaffold_number}->{$scaffold_number.'.'.$contig_number} = 1;
    }

    my $sorted_contigs_list = $self->_get_ace_contig_names;
    
    #assign first scaffold to current scaffold name space
    my ( $current_scaffold ) = @{$sorted_contigs_list}[0] =~ /(\d+)\.\d+/;
    #scalar to hold updated values when scaf number changes
    my ( $updated_scaffold_number ) = @{$sorted_contigs_list}[0] =~ /(\d+)\.\d+/;
    #what the next split off scaffold number should be .. 1 greater than current largest scaffold number;
    my ( $last_scaffold_number ) = @{$sorted_contigs_list}[-1] =~ /(\d+)\.\d+/;
    #scalar to hold updated value of how much contig number will change
    my $contig_offset = 0;
    #store contigs to rename
    my %rename; 
    
    foreach my $contig (  @{$sorted_contigs_list} ) { #make sure it's sorted
	my ( $scaffold_number, $contig_number ) = $contig =~ /(\d+)\.(\d+)/;

	#when next scaffold reached set contig offset to 0 and reset upate scaffold to current scaffold number
	if ( not $scaffold_number == $current_scaffold ) {
	    $contig_offset = 0;
	    $updated_scaffold_number = $scaffold_number;
	}

	#when split point reached, figure out new contig number
	if ( exists $split_points->{$scaffold_number} ) { #split scaffold
	
	    my $updated_contig_number = $contig_number;
	    
	    #reached split point
	    if ( exists $split_points->{$scaffold_number}->{$contig} ) {
		#reset by subtracting contig number from contig number than add 1, ie, 
		my ( $split_contig_number ) = $contig =~ /\d+\.(\d+)/;
		$contig_offset = $split_contig_number - 1;
		$updated_scaffold_number = ++$last_scaffold_number;
	    }

	    $updated_contig_number -= $contig_offset;

	    my $old = $scaffold_number.'.'.$contig_number;
	    my $new = $updated_scaffold_number.'.'.$updated_contig_number;

	    #print $old.' updated to '.$new."\n";

	    #scaf contig 0 to first cutoff don't change so don't store it
	    $rename{$old} = $new unless $old == $new;
	}

	$current_scaffold = $scaffold_number;
    }

    return \%rename;
}
    
sub _get_scaffold_split_contig_names {
    my $self = shift;
    
    #verify contig names in pcap style
    for my $contig_name ( $self->split_contigs ) {
	unless ( $contig_name =~ /Contig\d+\.\d+/i ) {
	    $self->error_message("Expected pcap style contigs names like Contig4.5 but got: $contig_name");
	    return;
	}
    }
    
    return $self->split_contigs;
}

sub _get_ace_contig_names {
    my $self = shift;

    my $contig_numbers = {};

    my $fh = Genome::Sys->open_file_for_reading( $self->ace_file );

    while (my $line = $fh->getline) {
	if ( $line =~ /^CO\s+/ ) {
	    my ( $scaffold_number, $contig_number ) = $line =~ /CO\s+Contig(\d+)\.(\d+)\s+/;

	    unless ( defined $scaffold_number and defined $contig_number ) {
		$self->error_message("Failed to get contig number from line: $line, expected Contig#.#");
		return;
	    }
            # for sorting .. {$a<=>$b} doesn't seem to work with decimal numbers
	    $contig_numbers->{$scaffold_number}->{$contig_number} = 1; 
	}
    }

    $fh->close;

    unless ( scalar keys %{$contig_numbers} > 0 ) {
	$self->error_message("Failed to get any contig numbers from ace file");
	return;
    }

    my $sorted_list = $self->_sorted_array_ref_of_contig_numbers( $contig_numbers );

    $contig_numbers = undef;

    return $sorted_list;
}

sub _sorted_array_ref_of_contig_numbers {
    my ( $self, $hash ) = @_;
    my @ar;

    foreach my $scaf_number ( sort {$a<=>$b} keys %{$hash} ) {
	foreach my $ctg_number ( sort {$a<=>$b} keys %{$hash->{$scaf_number}} ) {
	    push @ar, $scaf_number.'.'.$ctg_number;
	}
    }

    unless ( scalar @ar > 0 ) {
	$self->error_message("Found no sorted contigs after attempting to sort a list of contig numbers");
	return;
    }

    return \@ar;
}

1;
