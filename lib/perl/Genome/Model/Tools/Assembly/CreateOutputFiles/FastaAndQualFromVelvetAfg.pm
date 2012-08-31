package Genome::Model::Tools::Assembly::CreateOutputFiles::FastaAndQualFromVelvetAfg;

use strict;
use warnings;

use Genome;
use AMOS::AmosLib;
use IO::File;

use Bio::Seq;
use Bio::Seq::Quality;
use Bio::SeqIO;

class Genome::Model::Tools::Assembly::CreateOutputFiles::FastaAndQualFromVelvetAfg {
    is => 'Genome::Model::Tools::Assembly::CreateOutputFiles',
    has => [
        afg_file => {
            is => 'Text',
            doc => 'Velvet afg file to get fasta and qual info from',
            is_optional => 1,
        },
        directory => {
            is => 'Text',
            doc => 'Main assembly directory .. above edit_dir',
        },
    ],
};

sub help_brief {
    'Tool to create contigs.bases and contigs.qual files from velvet afg file';
}

sub help_detail {
    "Tool to create a separate files of fastas and qualities of contigs from velvet afg file";
}

sub execute {
    my $self = shift;

    my $afg_file = ($self->afg_file) ? $self->afg_file : $self->directory.'/velvet_asm.afg';

    unless (-s $afg_file) {
	$self->error_message("Can't find velvet afg file: ".$afg_file);
	return;
    }

    #read in afg file
    my $afg_fh = Genome::Sys->open_file_for_reading($afg_file)
	or return;

    #make edit_dir
    unless (-d $self->directory.'/edit_dir') {
	Genome::Sys->create_directory($self->directory.'/edit_dir');
    }

    #write out contigs.bases file
    my $f_io = Bio::SeqIO->new(-format => 'fasta', -file => ">".$self->contigs_bases_file);

    #write out contigs.quals file
    my $q_io = Bio::SeqIO->new(-format => 'qual', -file => ">".$self->contigs_quals_file);

    while (my $record = getRecord($afg_fh)) {
	my ($rec, $fields, $recs) = parseRecord($record);

	if ($rec eq 'CTG') { #contigs
	    #figure out contig naming
	    my ($sctg_num, $ctg_num) = split('-', $fields->{eid});
	    my $contig_id = 'Contig'.--$sctg_num.'.'.++$ctg_num;

	    #write fasta
	    my $seq = $fields->{seq};
	    $seq =~ s/\n//g; #contig seq is written in multiple lines .. remove end of line
	    my $seq_obj = Bio::Seq->new(-display_id => $contig_id, -seq => $seq);
	    $f_io->write_seq($seq_obj);
	    #write qual
	    my $qual = $fields->{qlt};
	    $qual =~  s/\n//g;
	    my @quals;
	    for my $i (0..length($qual)-1) {
                unless (substr($seq, $i, 1) eq '*') {
                    push @quals, ord(substr($qual, $i, 1)) - ord('0');
                }
            }
	    my $qual_obj = Bio::Seq::Quality->new(-display_id => $contig_id, -seq => $seq, -qual => \@quals);
	    $q_io->write_seq($qual_obj);
	}
    }

    $afg_fh->close;
    return 1;
}

1;
