package Genome::Model::Tools::Lims::AssemblyGetCloneReads;

use strict;
use warnings;

use Genome;
use IO::File;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Seq::Quality;
use Bio::SeqIO::fastq;
use Cwd;


class Genome::Model::Tools::Lims::AssemblyGetCloneReads {
    is => 'Command',
    has => [
	clone => {
	    is => 'Text',
	    doc => 'Clone name to get reads for',
	    is_optional => 1,
	},
	list_of_clones => {
	    is => 'Text',
	    doc => 'List of clones names to get read for',
	    is_optional => 1,
	},
	read_type => {
	    is => 'Text',
	    doc => 'Read types to get',
	    valid_values => ['all', 'end'],
	},
	output_dir => {
	    is => 'Text',
	    doc => 'Directory to output data to',
	    is_mutable => 1,
	    is_optional => 1,
	},
	dump_separately => {
	    is => 'Boolean',
	    doc => 'Option create separate fasta and quals for each clone',
	    is_optional => 1,
	    default_value => 0,
	},
    ],
};

sub help_brief {
    'Tool to dump reads for finishing clones'
}

sub help_synopsis {
    return <<"EOS"
gmt assembly get-clone-reads --list-of-clones /gscmnt/111/my_assembly/list --read-type end --output-dir /gscmnt/111/my_assembly
gmt assembly get-clone-reads --list-of-clones /gscmnt/111/my_assembly/list --read-type end --output-dir /gscmnt/111/my_assembly --dump-separately
gmt assembly get-clone-reads --clone H_GY-34H03 --read-type all
EOS
}

sub help_detail {
    return <<EOS
Tools to dump fasta and qual for a clone or list of clones.  Use --read-type end to dump CLONEEND
sequences and --read-type all to dump all reads.  Use --dump-separately option to create separate
fasta and qual for each clone; otherwise, tool will dump all fasta and qual into
combined_clone.fasta/qual files
EOS
}

sub execute {
    my $self = shift;

    my @clones;
    unless (@clones = $self->_get_clones()) {
	$self->warning_message("Failed to get clones names to dump");
	return;
    }

    my @seq_ids;

    foreach my $clone_name (@clones) {
	$clone_name =~ s/\s+$//;
	$self->debug_message("Getting reads for clone: $clone_name\n");

	#get seq_ids
	my @ids;
	unless (@ids = $self->_get_seq_ids($clone_name)) {
	    $self->debug_message("Failed to find any reads for clone: $clone_name");
	    next;
	}

	$self->debug_message("Found ".scalar @ids." reads");

	@seq_ids = (@seq_ids, @ids);

	#create separate fasta and quals for each clone 
	if ($self->dump_separately) {
	    #create re_ids file for seq_dump input
	    unless ($self->_write_read_ids_to_file($clone_name, @ids)) {
		$self->error_message("Failed to write read ids to file $clone_name.re_ids");
		return;
	    }
	    #usee seq_dump to dump fasta and quals
	    unless ($self->_dump_data ($clone_name)) {
		$self->error_message("Failed to dump fasta and qual for $clone_name");
		return;
	    }
	}
    }
    #dump all fasta and qual into combined fasta and qual files
    unless ($self->dump_separately) {
	#make re_ids file
	unless ($self->_write_read_ids_to_file('combined_clones', @seq_ids)) {
	    $self->error_message("Failed to write read ids to file combined_clones.re_ids");
	    return;
	}
	#dump fasta and quals
	unless ($self->_dump_data ('combined_clones')) {
	    $self->error_message("Failed to dump fasta and qual for combined_clones");
	    return;
	}
    }

    #use gsc::sequence::item to dump fasta/qual individually
    #TODO - seems to convert Ns in seq to dashes .. will look into it
    #unless ($self->_dump_reads_individually($clone_name)) {
	#$self->error_message("Failed to dump read individually");
	#return
    #}

    return 1;
}

sub _get_seq_ids {
    my ($self, $clone_name) = @_;

    my $clone = GSC::Clone->get(clone_name => $clone_name);

    my @clone_ids_ncbi;
    my $dna_ext_name = GSC::DNAExternalName->get(name_type  => 'ncbi clone id',
						 dna_id     => $clone->id,
	);
    my $external_id_ncbi = $dna_ext_name->name if($dna_ext_name);
    push @clone_ids_ncbi, $external_id_ncbi if $external_id_ncbi;
    
    my $agi_id_ncbi = $clone->convert_name(output => 'agi');
    push @clone_ids_ncbi, $agi_id_ncbi if $agi_id_ncbi;

    my $clone_id_ncbi = $clone->convert_name(output => 'ncbi');
    push @clone_ids_ncbi, $clone_id_ncbi if $clone_id_ncbi;
    push @clone_ids_ncbi, $clone_name;
    
    #try clone_id and trace_type_code to get reads
    my %params = (
	'clone_id' => [@clone_ids_ncbi],
	);
    $params{'trace_type_code'} = 'CLONEEND' if $self->read_type eq 'end';

    my @srs = GSC::Sequence::Read->get( %params	);

    #try clone_id and stratege
    unless (@srs) {
	delete $params{'trace_type_code'} if exists $params{'trace_type_code'};
	$params{'strategy'} = 'CLONEEND' if $self->read_type eq 'end';
	my @srs2 = GSC::Sequence::Read->get( %params );
	push @srs,@srs2;
    }

    #try template_id adn trace_type_code
    unless (@srs) {
	delete $params{'clone_id'};
	$params{'template_id'} = [@clone_ids_ncbi];
	delete $params{'strategy'} if exists $params{'strategy'};
	$params{'trace_type_code'} = 'CLONEEND' if $self->read_type eq 'end';
	my @srs2 = GSC::Sequence::Read->get( %params );
	push @srs,@srs2;
    }

    return @srs;
}

sub _write_read_ids_to_file {
    my ($self, $clone_name, @seq_reads) = @_;

    $self->output_dir(cwd()) unless $self->output_dir;

    unlink $self->output_dir.'/'.$clone_name.'.re_ids';
    my $fh = Genome::Sys->open_file_for_writing($self->output_dir.'/'.$clone_name.'.re_ids') ||
	return;
   
    foreach (@seq_reads) {
	next unless ($_->pass_fail_tag eq 'PASS');
	$fh->print ($_->seq_id."\n");
    }

    $fh->close;

    unless (-s $self->output_dir.'/'.$clone_name.'.re_ids') {
	$self->error_message("Failed to write $clone_name.re_ids file or file is blank");
	return;
    }

    return 1;
}

sub _dump_data {
    my ($self, $clone_name) = @_;

    my $re_ids_file = $self->output_dir.'/'.$clone_name.'.re_ids';
    my $fasta_out = $self->output_dir.'/'.$clone_name.'.fasta';
    my $qual_out = $self->output_dir.'/'.$clone_name.'.fasta.qual';

    unlink $fasta_out, $qual_out;
    $self->debug_message("Dumping fasta\n");
    if (system("seq_dump --input-file $re_ids_file --output type=fasta,file=$fasta_out,maskq=1,maskv=1,nocvl=35")) {
	$self->error_message("Failed to dump fasta for $clone_name");
	return;
    }
    $self->debug_message("Dumping quality\n");
    if (system ("seq_dump --input-file $re_ids_file --output type=qual,file=$qual_out")) {
	$self->error_message("Failed to dump qual for $clone_name");
	return;
    }

    return 1;
}

sub _get_clones {
    my $self = shift;
    my @clones;

    if ($self->list_of_clones) {
	my $fh = Genome::Sys->open_file_for_reading($self->list_of_clones) ||
	    return;
	while (my $line = $fh->getline) {
	    next if $line =~ /^\s+$/;
	    chomp $line;
	    unless (grep(/^$line$/, @clones)) {
		push @clones, $line;
	    }
	}
	$fh->close;
    }
    push @clones, $self->clone if $self->clone;

    return @clones;
}

sub _dump_reads_individually {
    my ($self, $clone_name) = @_;

    $self->debug_message("Dumping fasta and qual for $clone_name");

    my $fasta_out = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$self->output_dir.'/'.$clone_name.'.fasta');
    my $qual_out = Bio::SeqIO->new(-format => 'qual', -file => '>'.$self->output_dir.'/'.$clone_name.'.fasta.qual');

    my $re_ids_file = $self->output_dir.'/'.$clone_name.'.re_ids';
    my $fh = Genome::Sys->open_file_for_reading($re_ids_file) ||
	return;
    while (my $re_id = $fh->getline) {
	chomp $re_id;
	my $seq = GSC::Sequence::Item->get($re_id);

	unless ($seq) {
	    $self->debug_message("Failed to get sequence object for re_id: $re_id");
	    next;
	}
	unless ($seq->sequence_base_string ge 10) {
	    $self->debug_message("Skipping, read is under 10 bases: ".$seq->gsc_trace_name);
	    next;
	}

	my $seq_obj = Bio::Seq->new(-display_id => $seq->gsc_trace_name, -seq => $seq->sequence_base_string);
	$fasta_out->write_seq($seq_obj);

	my $qual_obj = Bio::Seq::Quality->new(-display_id => $seq->gsc_trace_name, -seq => $seq->sequence_base_string, -qual => $seq->sequence_quality_string);
	$qual_out->write_seq($qual_obj);
    }
    $fh->close;

    return 1;
}
