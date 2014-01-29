package Genome::Model::Tools::Velvet::Base;

use strict;
use warnings;

use POSIX;
use Genome;
use AMOS::AmosLib;
use Data::Dumper;
use Regexp::Common;

use Bio::SeqIO;

class Genome::Model::Tools::Velvet::Base {
    is  => 'Command::V2',
    is_abstract  => 1,
    has_optional => [
        version => {
            is   => 'String',
            doc  => 'velvet version, must be valid velvet version number like 0.7.22, 0.7.30. It takes installed as default.',
            default => 'installed',
        },
    ],
};

sub resolve_version {
    my $self = shift;

    my ($type) = ref($self) =~ /\:\:(\w+)$/;
    $type = 'velvet'.lc(substr $type, 0, 1);

    my $ver = $self->version;
    $ver = 'velvet_'.$ver unless $ver eq 'installed';
    
    my @uname = POSIX::uname();
    $ver .= '-64' if $uname[4] eq 'x86_64';
    
    my $exec = $ENV{GENOME_SW} . "/velvet/$ver/$type";
    unless (-x $exec) {
        $self->error_message("$exec is not excutable");
        return;
    }

    return $exec;
}

sub load_sequence_seek_positions {
    my ($self, $seq_file) = @_;

    my @seek_positions;
    my $fh = Genome::Sys->open_file_for_reading( $seq_file );

    my $seek_pos = $fh->tell;
    my $io = Bio::SeqIO->new(-format => 'fasta', -fh => $fh);
    while (my $seq = $io->next_seq) {
	my ($read_index) = $seq->desc =~ /(\d+)\s+\d+$/; #numerically ordered .. could just do $c++;
	unless ($read_index) {
            $self->error_message("Failed to get read index number from seq->desc: ".$seq->desc);
            return;
        }
	
	$seek_pos = ( $seek_pos == 0 ) ? $seek_pos : $seek_pos - 1;

        $seek_positions[$read_index] = $seek_pos;
        $seek_pos = $fh->tell;
    }
    $fh->close;

    return \@seek_positions;
}

sub resolve_edit_dir {
    my $self = shift;
    return $self->assembly_directory.'/edit_dir';
}

sub create_edit_dir {
    my $self = shift;

    my $edit_dir = $self->resolve_edit_dir;
    unless ( -d $edit_dir ) {
        Genome::Sys->create_directory($edit_dir);
    }

    return 1;
}

sub _global_default_gap_size {
    return 100;
}

sub global_gap_size {
    my $self = shift;
 
    my $gap_size = eval{ $self->default_gap_size; };
    return $gap_size if $gap_size;

    return $self->_global_default_gap_size;
}

sub get_scaffold_info_from_afg_file {
    my $self = shift;

    # return stored scaf info from file if present
    if ( -s $self->scaffold_info_stor_file ) {
        my $scaf_info = Storable::retrieve $self->scaffold_info_stor_file;
        return $self->filter_contigs_by_params( $scaf_info );
    }

    my %scaf_info;
    my %ctgs_in_sctgs;
    my $assembled_position;
    my $fh = Genome::Sys->open_file_for_reading( $self->velvet_afg_file );
    while ( my $record = getRecord($fh) ) {
        my ( $rec, $fields, $recs ) = parseRecord( $record );
        if ( $rec eq 'CTG' ) {
            my $seq = $fields->{seq};
            $seq =~ s/\n//g;
            my( $sctg, $ctg ) = $fields->{eid} =~ /^(\d+)-(\d+)$/;
            die 'Could not get supercontig and contig numbers from name: '.$fields->{eid} if
                not defined $sctg and not defined $ctg;
            # set contig lengths
            $scaf_info{$sctg.'.'.$ctg}{contig_length} = length $seq;
            # assembled position
            $scaf_info{$sctg.'.'.$ctg}{assembled_position} = ++$assembled_position;
            # velvet contig name
            $scaf_info{$sctg.'.'.$ctg}{velvet_name} = $sctg.'.'.$ctg;
            # track contigs in scaffolds
            push @{$ctgs_in_sctgs{$sctg}}, $ctg;
        }
    }
    $fh->close;
    
    my @method_names = (qw/ gap_length supercontig_start_pos pcap_contig_name
                            next_contig_name contig_order scaffold_contig_names /);

    for my $name ( @method_names ) {
        my $method = '_add_'.$name;
        if ( not $self->$method(\%scaf_info,\%ctgs_in_sctgs) ) {
            $self->debug_message("Failed to run add $name to scaffold infos");
            return;
        }
    }

    Storable::nstore \%scaf_info, $self->scaffold_info_stor_file;
    return $self->filter_contigs_by_params(\%scaf_info);
}

sub _add_gap_length {
    my( $self, $scafs, $list ) = @_;
    for my $sctg ( keys %$list ) {
        my @ctgs = @{$list->{$sctg}};
        pop @ctgs; # no gap for last ctg
        next if not @ctgs; # single contig scaffolds
        for my $ctg ( @ctgs ) {
            my $name = $sctg.'.'.$ctg;
            $scafs->{$name}->{gap_size} = $self->global_gap_size;
        }
    }
    return 1;
}

sub _add_supercontig_start_pos {
    my ( $self, $scafs, $list ) = @_;
    for my $sctg ( keys %$list ) {
        my $current_position = 1;
        for my $ctg ( @{$list->{$sctg}} ) {
            my $name = $sctg.'.'.$ctg;
            $scafs->{$name}->{contig_start_position} = $current_position;
            my $contig_length = $scafs->{$name}->{contig_length};
            my $gap_size = ( exists $scafs->{$name}->{gap_size} ) ?
                $scafs->{$name}->{gap_size} :
                0;
            $current_position += ( $contig_length + $gap_size );
        }
    }
    return 1;
}

sub _add_next_contig_name {
    my ($self, $scafs, $list) = @_;
    for my $sctg ( keys %$list ) {
        for ( 0 .. $#{$list->{$sctg}} ) {
            next if $_ == $#{$list->{$sctg}}; # no next ctg for last ctg
            my $ctg = ${$list->{$sctg}}[$_];
            my $next_ctg = ${$list->{$sctg}}[$_ + 1];
            my $curr_name = $sctg.'.'.$ctg;
            my $next_name = $sctg.'.'.$next_ctg;
            $scafs->{$curr_name}->{next_contig_name} = $next_name;
        }
    }
    return 1;
}

sub _add_contig_order {
    my ( $self, $scafs, $list ) = @_;
    for my $sctg ( keys %$list ) {
        my $order = 0;
        for my $ctg ( @{$list->{$sctg}} ) {
            my $name = $sctg.'.'.$ctg;
            my $contig_count = @{$list->{$sctg}};
            $scafs->{$name}->{contig_order} = ++$order.'/'.$contig_count;
        }
    }
    return 1;
}

sub _add_supercontig_length {
    my ( $self, $scafs, $list ) = @_;
    my %lengths;
    for my $sctg ( keys %$list ) {
        for my $ctg ( @{$list->{$sctg}} ) {
            my $name = $sctg.'.'.$ctg;
            my $length = $scafs->{$name}->{contig_length};
            $lengths{$sctg} += $length;
        }
    }
    for my $name ( keys %$scafs ) {
        my ( $sctg, $ctg ) = $self->_sctg_ctg_name_from_velvet( $name );
        $scafs->{$name}->{supercontig_length} = $lengths{$sctg};
    }
    return 1;
}

sub _add_pcap_contig_name {
    my ( $self, $scafs ) = @_;
    for my $name ( keys %$scafs ) {
        my ( $sctg, $ctg ) = $self->_sctg_ctg_name_from_velvet( $name );
        $scafs->{$name}->{pcap_name} = 'Contig'.($sctg - 1).'.'.($ctg + 1);
    }
    return 1;
}

sub _add_scaffold_contig_names {
    my( $self, $scafs, $list ) = @_;
    for my $name ( keys %$scafs ) {
        my ( $sctg, $ctg ) = $self->_sctg_ctg_name_from_velvet( $name );
        @{$scafs->{$name}->{contigs_in_scaffold}} = @{$list->{$sctg}};
    }
    return 1;
}

sub contig_infos {
    my $self = shift;
    my $scafs = $self->get_scaffold_info_from_afg_file;
    return $self->sort_contigs_in_assembled_order($scafs);
}

sub sort_contigs_in_assembled_order {
    my( $self, $scafs ) = @_;
    my @ordered;
    for ( 1 .. scalar ( keys %$scafs ) ) {
        for my $ctg ( keys %$scafs ) {
            if ( $scafs->{$ctg}->{assembled_position} == $_ ) {
                push @ordered, $scafs->{$ctg};
                last;
            }
        }
    }
    return \@ordered;
}

sub filter_contigs_by_params {
    my ($self, $scafs) = @_;
    $self->add_filterable_properties($scafs); # add error
    # filter my min_length
    for my $name ( keys %$scafs ) {
        if ( $scafs->{$name}->{filtered_supercontig_length} < $self->min_contig_length ) {
            $scafs->{$name}->{is_removed} = 1;
        }
        elsif ( $scafs->{$name}->{contig_length} < $self->min_contig_length ) {
            $scafs->{$name}->{is_removed} = 1;
        }
    }
    # update gap to reflect removed reads info
    for my $name ( keys %$scafs ) {
        # no gap info needed for removed contig
        next if exists $scafs->{$name}->{is_removed};
        # no gap info for last contig in scaffold
        next if not exists $scafs->{$name}->{next_contig_name};
        #my $gap_size = $scafs->{$name}->{gap_size};
        my $gap_size = $self->global_gap_size;
        my $next_ctg_names = $self->get_downstream_contig_names($scafs,$name);
        for my $next_ctg_name ( @$next_ctg_names ) {
            if ( not exists $scafs->{$next_ctg_name}->{is_removed} ) {
                $scafs->{$name}->{processed_gap_size} = $gap_size;
                last;
            }
            #elsif ( exists $scafs->{$next_ctg_name}->{is_removed} ) {
            #    # next contigs is filtered out add next contig length and gap size to current
            #    $gap_size += $scafs->{$next_ctg_name}->{contig_length};
            #    $gap_size += $scafs->{$next_ctg_name}->{gap_size};
            #}
        }
    }
    # set last contig in scaffold
    for my $name ( keys %$scafs ) {
        next if $scafs->{$name}->{is_removed};
        my $is_last_ctg = 1;
        my $next_ctg_names = $self->get_downstream_contig_names($scafs,$name);
        for my $next_ctg ( @$next_ctg_names ) {
            if ( $scafs->{$next_ctg}->{is_last_scaffold_contig} ) {
                # last contig already set for this scaffold
                $is_last_ctg = 0;
                last;
            }
            next if $scafs->{$next_ctg}->{is_removed};
            $is_last_ctg = 0;
        }
        if ( $is_last_ctg ) {
            $scafs->{$name}->{is_last_scaffold_contig} = 1;
        }
    }
    return $scafs;
}

sub get_downstream_contig_names {
    my ( $self, $scafs, $name ) = @_;
    my @names;
    my ( $sctg, $ctg ) = $self->_sctg_ctg_name_from_velvet( $name );
    for my $nth ( 0 .. $#{$scafs->{$name}->{contigs_in_scaffold}} ) {
        if ( ${$scafs->{$name}->{contigs_in_scaffold}}[$nth] == $ctg ) {
            # skip current nth ctg
            ++$nth;
            for ( $nth .. $#{$scafs->{$name}->{contigs_in_scaffold}} ) {
                my $nth_ctg_name = $sctg.'.'.${$scafs->{$name}->{contigs_in_scaffold}}[$_];
                push @names, $nth_ctg_name;
            }
        }
    }
    return \@names;
}

sub add_filterable_properties {
    my ( $self, $scafs ) = @_;
    my $lengths = $self->filtered_supercontig_lengths( $scafs );
    # add filtered scaffold length
    for my $name ( keys %$scafs ) {
        my( $sctg, $ctg ) = $self->_sctg_ctg_name_from_velvet( $name );
        $scafs->{$name}->{filtered_supercontig_length} = $lengths->{$sctg};
    }
    return 1;
}

sub filtered_supercontig_lengths {
    my ( $self, $scaf_info ) = @_;
    my %lengths;
    for my $name ( keys %$scaf_info ) {
        my $length = ( $scaf_info->{$name}->{contig_length} >= $self->min_contig_length ) ?
            $scaf_info->{$name}->{contig_length} :
            0;
        my ( $sctg, $ctg ) = $self->_sctg_ctg_name_from_velvet( $name );
        $lengths{$sctg} += $length;
    }
    return \%lengths;
}

sub _sctg_ctg_name_from_velvet {
    my ( $self, $name ) = @_;
    my ( $sctg, $ctg ) = $name =~ /^(\d+)\.(\d+)$/;
    die "Could not derive supercontig and contig number from velvet contig name: $name\nExpected patter like ##.##" if
        not defined $sctg or not defined $ctg;
    return( $sctg, $ctg );
}

#post assemble standard output files
sub contigs_bases_file {
    my $self = shift;
    return $self->resolve_contigs_bases_file;
}

sub resolve_contigs_bases_file {
    my $self = shift;
    return $self->resolve_edit_dir.'/contigs.bases';
}

sub contigs_quals_file {
    my $self = shift;
    return $self->resolve_contigs_quals_file;
}

sub resolve_contigs_quals_file {
    return $_[0]->resolve_edit_dir.'/contigs.quals';
}

sub contigs_cmt_file {
    return $_[0]->assembly_directory.'/edit_dir/contigs.cmt';
}

sub gap_sizes_file {
    my $self = shift;
    return $self->resolve_gap_sizes_file;
}
sub resolve_gap_sizes_file {
    my $self = shift;
    return $self->resolve_edit_dir.'/gap.txt';
}

sub read_info_file {
    my $self = shift;
    return $self->resolve_read_info_file;
}
sub resolve_read_info_file {
    my $self = shift;
    return $self->resolve_edit_dir.'/readinfo.txt';
}

sub reads_placed_file {
    my $self = shift;
    return $self->resolve_reads_placed_file;
}
sub resolve_reads_placed_file {
    my $self = shift;
    return $self->resolve_edit_dir.'/reads.placed';
}

sub reads_unplaced_file {
    return $_[0]->assembly_directory.'/edit_dir/reads.unplaced';
}

sub reads_unplaced_fasta_file {
    return $_[0]->assembly_directory.'/edit_dir/reads.unplaced.fasta';
}

sub stats_file {
    return $_[0]->assembly_directory.'/edit_dir/stats.txt';
}

sub supercontigs_agp_file {
    return $_[0]->assembly_directory.'/edit_dir/supercontigs.agp';
}

sub supercontigs_fasta_file {
    return $_[0]->assembly_directory.'/edit_dir/supercontigs.fasta';
}

#other files
sub scaffold_info_stor_file {
    return $_[0]->assembly_directory.'/scaffolds.stor';
}

sub read_names_sqlite {
    return $_[0]->assembly_directory.'/velvet_reads.sqlite';
}

sub input_collated_fastq_file {
    my $self = shift;

    my @files = glob( $self->assembly_directory."/*collated.fastq" );

    unless ( @files == 1 ) {
	$self->error_message("Expected 1 *collated.fastq file but got " . scalar @files);
	return;
    }

    return $files[0];
}

sub core_gene_survey_file {
    return $_[0]->assembly_directory.'/edit_dir/core_gene_survey_result';
}

#velvet generated files
sub velvet_afg_file {
    my $self = shift;
    return $self->resolve_afg_file;
}
sub resolve_afg_file {
    my $self = shift;
    return $self->assembly_directory.'/velvet_asm.afg';
}

sub velvet_contigs_fa_file {
    return $_[0]->assembly_directory.'/contigs.fa';
}

sub velvet_sequences_file {
    return $_[0]->assembly_directory.'/Sequences';
}

1;

