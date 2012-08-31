package Genome::Model::Tools::Assembly::UnpackTraces;

use strict;
use warnings;
use Genome;
use XML::Simple;
use Bio::SeqIO;
use Getopt::Long;
use Cwd;
use Data::Dumper;
use File::Path;

class Genome::Model::Tools::Assembly::UnpackTraces {
    is => 'Command',
    has => [
	    trace_file => {
		type => 'String',
		is_optional => 0,
		doc => "Ncbi downloaded tar zipped traces file"
	        },
	    clip_vector => {
		type => 'Boolean',
		is_optional => 1,
		doc => "set program to clip vector",
	        },
	    clip_quality => {
		type => 'Boolean',
		is_optional => 1,
		doc => "set program to clip low quality",
	        },
	    data_out_dir => {
		type => 'String',
		is_optional => 1,
		doc => "directory path to ouput files",
	        },
	    zip_files => {
		type => 'Boolean',
		is_optional => 1,
		doc => "set to zip output files",
	        },
	    save_trace_directory => {
		type => 'Boolean',
		is_optional => 1,
		doc => "delete directory created from archive",
		default_value => 0,
	        },
	    ],
};

sub help_brief {
    'Tool to unpack, optionally clip then create data files for assembling'
}

sub help_detail {
    return <<"EOS"
This module is a tool to unpack fasta and quality data downloaded from
ncbi and create input files for assemblies.  User can optionally clip
vector and or poor quality regions as specified in TRACEINFO.xml file 
EOS
}

sub execute {
    my ($self) = @_;

    my $trace_file = $self->trace_file;

    unless (-s $trace_file) {
	$self->error_message("No such file or file is blank: $trace_file");
	return;
    }

    my (@xml_files) = $self->_get_xml_file_names ($trace_file);

    my (undef, $base_path) = File::Spec->splitpath($xml_files[0]);

    my ($base_name) = $base_path =~ /^(\S+)\//;

#OUTPUT DIR
    my $output_dir = cwd();
    $output_dir = $self->data_out_dir if $self->data_out_dir;

#OUTPUT FILES
    my $fasta_out = $output_dir.'/'.$base_name.'.fasta';
    unlink $fasta_out if -s $fasta_out;
    unlink $fasta_out.'.gz' if -s $fasta_out.'.gz';
    my $qual_out = $output_dir.'/'.$base_name.'.fasta.qual';
    unlink $qual_out if -s $qual_out;
    unlink $qual_out.'.gz' if -s $qual_out.'.gz';

#UNPACK TRACE FILE
    if ( system ("tar xvf $trace_file") ) {
	$self->error_message ("Unable to get file info from $trace_file");
	return;
    }

#PROCESS DATA
    foreach my $xml_file (@xml_files) {
	my $xml = XML::Simple->new( Forcearray => [ qw( book ) ], keyattr => { book => 'isbn' } );
	my $xml_data = $xml->XMLin($xml_file);
	foreach my $read_info (@{$xml_data->{trace}}) {
	    my $fasta_file = $base_path.$read_info->{base_file};
	    my $qual_file = $base_path.$read_info->{qual_file};

	    next unless -s $fasta_file and -s $qual_file;

	#CREATE GC READ NAME
	    my $read_name = $self->_configure_read_name ($read_info);

	#BUILD HEADER
	    my $desc = $self->_build_header_desc ($read_info);

	#FIGURE OUT CLIPPING POSITIONS
	    my $read_length = $read_info->{ncbi_trace_archive}->{basecall_length};

	    my ($clip_left_pos, $clip_right_pos);

	    if ($self->clip_vector and $self->clip_quality) {
		$clip_left_pos = ($read_info->{clip_vector_left} > $read_info->{clip_quality_left}) ?
		    $read_info->{clip_vector_left} : $read_info->{clip_quality_left};
		$clip_right_pos = ($read_info->{clip_vector_right} < $read_info->{clip_quality_right}) ?
		    $read_info->{clip_vector_right} : $read_info->{clip_quality_right};
	    }
	    elsif ($self->clip_vector) {
		$clip_left_pos = $read_info->{clip_vector_left};
		$clip_right_pos = $read_info->{clip_vector_right};
	    }
	    elsif ($self->clip_quality) {
		$clip_left_pos = $read_info->{clip_quality_left};
		$clip_right_pos = $read_info->{clip_quality_right};
	    }
	    else {
		$clip_left_pos = 1;
		$clip_right_pos = $read_info->{ncbi_trace_archive}->{basecall_length};
	    }

	#CLIPPING STARTS BEFORE AND AFTER LEFT AND RIGHT POS
	    $clip_left_pos -= 1;
	    $clip_right_pos += 1;

	#MAKE SURE CLIP LEFT POS DOESN'T EXTEND BEYOND CLIP
	#RIGHT POSITION
	    if ($clip_left_pos > $clip_right_pos) {
		$clip_right_pos = $clip_left_pos + 1;
	    }

	#CLIPPED FASTA AND NEW QUAL OUTPUT FILES
	    my $out_fio = Bio::SeqIO->new(-format => 'fasta', -file => ">>$fasta_out");
	    my $out_qio = Bio::SeqIO->new(-format => 'qual', -file => ">>$qual_out");

	    my $fio = Bio::SeqIO->new(-format => 'fasta', -file => $fasta_file);
	    while (my $seq = $fio->next_seq()) {
		my $fasta = $seq->seq;
		my $pre_process_length = length $fasta;

		if ($self->clip_vector or $self->clip_quality) {
		    $fasta = substr $fasta, $clip_left_pos, $clip_right_pos - $clip_left_pos - 1;
		    my $left_xs = $self->_create_x_string ($clip_left_pos);
		    my $right_xs = $self->_create_x_string ($read_length - $clip_right_pos + 1);
		    $fasta = $left_xs.$fasta.$right_xs;
		}
		my $post_process_length = length $fasta;

		unless ($pre_process_length == $post_process_length) {
		    $self->error_message ("Pre and post clipping fasta lenghts did not match for $read_name");
		    return;
		}

		$seq->seq ($fasta);
		$seq->id ($read_name);
		$seq->desc ($desc);
		$out_fio->write_seq($seq);
	    }
	    my $qio = Bio::SeqIO->new(-format => 'qual', -file => $qual_file);
	    while (my $qual = $qio->next_seq()) {
#	    print map {$_} @{$qual->qual};
		$qual->id ($read_name);
		$qual->desc ($desc);
		$qual->qual ($qual->qual);
		$out_qio->write_seq ($qual);
	    }
	}
    }
    #ZIP OUTPUT FILES FOR PCAP
    if ($self->zip_files){
	if ( system ("gzip $fasta_out $qual_out") ) {
	    $self->error_message ("Warning: Unable to zip $fasta_out and $qual_out files");
	}
    }
    #REMOVE DIRECTORY CREATED FROM ARCHIVE SINCE IT'S NOT NEEDED FOR ASSEMBLIES
    if (! $self->save_trace_directory) {
	rmtree $base_path;
    }
}

sub _build_header_desc {
    my ($self, $readinfo) = @_;
    my $desc;
    my $read_name = $self->_configure_read_name ($readinfo);
    #CHROMAT FIELD
    $desc .= ' CHROMAT_FILE: '.$read_name;
    #PHD FIELD
    $desc .= ' PHD_FILE: '.$read_name.'.phd.1';
    #TEMPLATE
    $desc .= ' TEMPLATE: '.$readinfo->{template_id};
    #DIRECTION
    if (exists $readinfo->{trace_end}) {
	my $direction = ($readinfo->{trace_end} eq 'FORWARD') ? 'fwd' : 'rev';
	$desc .= ' DIRECTION: '.$direction;
    }
    #TIME
    my $time = localtime();
    if (exists $readinfo->{ncbi_trace_archive}->{load_date}) {
	$time = $readinfo->{ncbi_trace_archive}->{load_date};
    }
    $desc .= ' TIME: '.$time;
    #INSERT SIZE;
    if (exists $readinfo->{insert_size}) {
	$desc .= ' INSERT_SIZE: '.$readinfo->{insert_size};
    }
    #TI
    if (exists $readinfo->{ncbi_trace_archive}->{ti}) {
	$desc .= ' TI: '.$readinfo->{ncbi_trace_archive}->{ti};
    }
    return $desc;
}

sub _configure_read_name {
    my ($self, $readinfo) = @_;
    my $read_name = $readinfo->{template_id};
    if (exists $readinfo->{trace_end}) {
	my $ext = ($readinfo->{trace_end} eq 'FORWARD') ? '.b1' : '.g1';
	$read_name = $read_name.$ext;
    }
    return $read_name;
}

sub _create_x_string {
    my ($self, $number) = @_;
    return '' if $number < 1;
    my $string;
    my $c = 0;
    while ($c < $number) {
	$c++;
	$string .= 'X';
    }
    return $string;
}

sub _get_xml_file_names {
    my ($self, $file) = @_;
    my @files = `tar tzf $file`;
    unless (@files) {
	$self->error_message("Unable to extract xml file names from trace file");
	return;
    }
    chomp @files;
    my @xml_files;
    foreach (@files) {
	if ($_ =~ /\.xml$/i) {
	    push @xml_files, $_;
	}
    }
    return @xml_files;
}

1;
