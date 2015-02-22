#!/usr/bin/perl

use strict;

use Getopt::Long;
use FileHandle;
use File::Path qw(make_path);
use Storable qw(nstore_fd);
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::DB::Fasta;

our $VERSION = 2.7;

my $config = {};

my $count_args = scalar @ARGV;

GetOptions(
	$config,
	'input|i|gtf|g=s',
	'fasta|f=s',
	'species|s=s',
	'db_version|d=i',
	'cache_region_size=i',
	'dir=s',
	'help',
) or die "ERROR: Failed to parse command line options\n";

if(defined($config->{help}) || !$count_args) {
	usage();
	exit(0);
}

# check for errors
die "ERROR: No species specified\n" unless defined($config->{species});
die "ERROR: No whitespace allowed in species name\n" if $config->{species} =~ /\s/;
die "ERROR: No DB version specified\n" unless defined($config->{db_version});
die "ERROR: No FASTA file/directory specified\n" unless defined($config->{fasta});

# set defaults
$config->{cache_region_size} ||= 1000000;
$config->{write_cache} = 1;
$config->{dir} ||= $ENV{HOME}.'/.vep/';
$config->{compress} ||= 'zcat';

$config->{dir} .= $config->{species}.'/'.$config->{db_version};

if(defined($config->{fasta})) {
	die "ERROR: Specified FASTA file/directory not found" unless -e $config->{fasta};
	
	debug("Checking/creating FASTA index");
	$config->{fasta_db} = Bio::DB::Fasta->new($config->{fasta});
}

# create a coord system
$config->{coord_system} = Bio::EnsEMBL::CoordSystem->new(
	-NAME => 'chromosome',
	-RANK => 1,
);

my @fields = qw(seqname source feature start end score strand frame attributes comments);
my $line_num = 0;
$config->{dbID} = 1;
my ($prev_tr_id, $prev_chr, $by_region);
my $in_file_handle = new FileHandle;

if(defined($config->{input})) {
	
	# check defined input file exists
	die("ERROR: Could not find input file ", $config->{input}, "\n") unless -e $config->{input};
	
	if($config->{input} =~ /\.gz$/){
		$in_file_handle->open($config->{compress}." ". $config->{input} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
	}
	else {
		$in_file_handle->open( $config->{input} ) or die("ERROR: Could not read from input file ", $config->{input}, "\n");
	}
}

# no file specified - try to read data off command line
else {
	$in_file_handle = 'STDIN';
	debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
}


while(<$in_file_handle>) {
	chomp;
	
	my @split = split /\t/, $_;
	
	my $data;
	
	# parse split data into hash
	for my $i(0..$#split) {
		$data->{$fields[$i]} = $split[$i];
	}
	
	debug("Processing chromosome ".$data->{seqname}) if !defined($prev_chr) || $data->{seqname} ne $prev_chr;
	
	# parse attributes
	if(defined($data->{attributes})) {
		$data->{attributes} =~ s/^\s+//g;
		
		my %attribs;
		
		foreach my $pair(split /;\s*/, $data->{attributes}) {
			my ($key, $value) = split /\s+/, $pair;
			next unless defined($key) and defined($value);
			
			# remove quote marks
			$value =~ s/\"//g;
			
			$attribs{$key} = $value;
		}
		
		$data->{attributes} = \%attribs;
	}
	
	my $tr_id = $data->{attributes}->{transcript_id};
	
	my $ref = parse_data($config, $data);
	
	if(defined($prev_tr_id) && $prev_tr_id ne $tr_id) {
		
		my $prev_tr = $config->{transcripts}->{$prev_tr_id};
		
		# add to by_region hash
		push @{$by_region->{$prev_chr || $data->{seqname}}->{get_region($config, $prev_tr)}}, $prev_tr;
		
		# dump if into new region or new chrom
		if(defined($prev_chr) && $prev_chr ne $data->{seqname}) {
			export_data($config, $prev_chr, $by_region);
		}
	}
	
	$prev_tr_id = $tr_id;
	$prev_chr = $data->{seqname};
}

my $prev_tr = $config->{transcripts}->{$prev_tr_id};

# add to by_region hash
push @{$by_region->{$prev_chr}->{get_region($config, $prev_tr)}}, $prev_tr;

# dump remaining transcripts
export_data($config, $prev_chr, $by_region);

debug("All done!");

sub parse_data {
	my ($config, $data) = @_;
	
	return unless defined($data);
	
	# check defined feature type
	return unless defined($data->{feature}) && $data->{feature} =~ /exon|cds/i;#|(start|stop)_codon/i;
	
	# run data fix
	fix_data($config, $data);
	
	# create transcript if not already done
	$config->{transcripts}->{$data->{attributes}->{transcript_id}} ||= create_transcript($config, $data);
	
	my $method = 'parse_'.lc($data->{feature});
    my $method_ref = \&$method;
	
	#print "running $method\n";
    
    return &$method_ref($config, $data);
}

sub fix_data {
	my ($config, $data) = @_;
	
	# fix strand
	$data->{strand} = $data->{strand} =~ /\-/ ? -1 : 1;
}

# creates a new transcript object
sub create_transcript {
	my ($config, $data) = @_;
	
	my $transcript = new Bio::EnsEMBL::Transcript(
		-STABLE_ID => $data->{attributes}->{transcript_id},
		-BIOTYPE   => $data->{source},
		-SLICE     => get_slice($config, $data->{seqname}),
		-VERSION   => 1,
		-dbID      => $config->{dbID}++,
	);
	
	$transcript->{_gene_stable_id} = $data->{attributes}->{gene_id};
	
	return $transcript;
}

# creates a new exon object
sub parse_exon {
	my ($config, $data) = @_;
	
	my $exon = new Bio::EnsEMBL::Exon(
		-START  => $data->{start},
		-END    => $data->{end},
		-STRAND => $data->{strand},
		-SLICE  => get_slice($config, $data->{seqname}),
		-PHASE  => -1,		# default phase to -1
	);
	
	# hidden number field
	$exon->{_number} = $data->{attributes}->{exon_number};
	
	# get sequence
	if(defined($config->{fasta_db})) {
		my $seq = $config->{fasta_db}->seq($data->{seqname}, $data->{start} => $data->{end});
		reverse_comp(\$seq) if $data->{strand} < 0;
		$exon->{_seq_cache} = $seq;
	}
	
	# add it to the transcript
	my $tr = $config->{transcripts}->{$data->{attributes}->{transcript_id}};
	
	$tr->add_Exon($exon);
	
	# update/create the translation
	if($tr->biotype eq 'protein_coding') {
		$tr->{translation} ||= new Bio::EnsEMBL::Translation(
			-TRANSCRIPT => $tr,
			-VERSION    => 1,
		);
		
		weaken($tr->{translation}->{transcript});
	}
	
	else {
		$exon->phase(-1);
	}
	
	return $exon;
}

# modifies a transcript with coding start/end info
sub parse_cds {
	my ($config, $data) = @_;
	
	# update the coding_region_start/end
	my $tr = $config->{transcripts}->{$data->{attributes}->{transcript_id}};
	
	return unless $tr->biotype eq 'protein_coding';
	
	# get matching exon
	my ($matched_exon) = grep {$_->{_number} eq $data->{attributes}->{exon_number}} @{$tr->get_all_Exons};
	
	$tr->{translation}->start_Exon($matched_exon) unless defined($tr->{translation}->start_Exon);
	$tr->{translation}->end_Exon($matched_exon);
	
	my ($start_exon, $end_exon) = ($tr->{translation}->start_Exon, $tr->{translation}->end_Exon);
	
	# match start exon?
	if($data->{attributes}->{exon_number} eq $start_exon->{_number}) {
		my $offset;
		
		if($data->{strand} > 0) {
			$offset = ($data->{start} - $start_exon->start) + 1;
		}
		else {
			$offset = ($start_exon->end - $data->{end}) + 1;
		}
		
		$tr->{translation}->start($offset);
		
		# update phase
		$start_exon->phase($data->{frame});
	}
	
	# match end exon?
	if($data->{attributes}->{exon_number} eq $end_exon->{_number}) {
		my $offset;
		
		if($data->{strand} > 0) {
			$offset = ($data->{end} - $end_exon->start) + 1;
		}
		else {
			$offset = ($end_exon->end - $data->{start}) + 1;
		}
		
		$tr->{translation}->end($offset);
		
		# update phase
		$end_exon->phase($data->{frame});
	}
}

sub get_slice {
	my $config = shift;
	my $chr = shift;
	
	if(!defined($config->{slice_cache}->{$chr})) {
		$config->{slice_cache}->{$chr} = Bio::EnsEMBL::Slice->new(
			-COORD_SYSTEM      => $config->{coord_system},
			-START             => 1,
			-END               => $config->{fasta_db}->length($chr),
			-SEQ_REGION_NAME   => $chr,
			-SEQ_REGION_LENGTH => $config->{fasta_db}->length($chr)
		);
	}
	
	return $config->{slice_cache}->{$chr};
}

sub get_region {
	my $config = shift;
	my $obj = shift;
	
	my ($s, $e);
	$s = $config->{cache_region_size} * int($obj->{start} / $config->{cache_region_size});
	$e = $s + $config->{cache_region_size};
	$s += 1;
	
	return $s.'-'.$e;
}

sub get_dump_file_name {
    my $config = shift;
    my $chr    = shift;
    my $region = shift;
    my $type   = shift;
    
    $type ||= 'transcript';
    
    if($type eq 'transcript') {
        $type = '';
    }
    else {
        $type = '_'.$type;
    }
    
    my $dir = $config->{dir}.'/'.$chr;
    my $dump_file = $dir.'/'.$region.$type.(defined($config->{tabix}) ? '_tabix' : '').'.gz';
    
    # make directory if it doesn't exist
    if(!(-e $dir) && defined($config->{write_cache})) {
        make_path($dir);
    }
    
    return $dump_file;
}

sub export_data {
	my $config = shift;
	my $chr = shift;
	my $hash = shift;
	
	foreach my $region(keys %{$hash->{$chr}}) {
		
		my @array = sort {$a->start <=> $b->start || $a->end <=> $b->end} @{$hash->{$chr}->{$region}};
		
		#foreach my $tr(@array) {
		#	check_transcript($config, $tr);
		#}
		
		dump_transcript_cache($config, {$chr => \@array}, $chr, $region);
		
		# remove all the dumped transcripts
		delete $config->{transcripts}->{$_->stable_id} for @array;
		delete $hash->{$chr}->{$region};
	}
}

sub check_transcript {
	my $config = shift;
	my $tr = shift;
	
	my @errors;
	
	push @errors, "Object is not a transcript" unless $tr->isa('Bio::EnsEMBL::Transcript');
	push @errors, "No exons found" unless scalar @{$tr->get_all_Exons};
	push @errors, "Exon missing phase" if grep {not defined $_->phase} @{$tr->get_all_Exons};
	
	#if(scalar @errors) {
	#	print join "\n", @errors;
	#	print "\n";
	#	
	#	
	#	use Data::Dumper;
	#	$Data::Dumper::Maxdepth = 3;
	#	warn Dumper $tr;
	#	
	#	delete $tr->{slice};
	#	
	#	print $tr->translate."\n";
	#	
	#	die "ERROR!\n";
	#}
}

# dumps out transcript cache to file
sub dump_transcript_cache {
    my $config = shift;
    my $tr_cache = shift;
    my $chr = shift;
    my $region = shift;
    
    #debug("Dumping cached transcript data") unless defined($config->{quiet});
    
    my $dump_file = get_dump_file_name($config, $chr, $region, 'transcript');
    
    #debug("Writing to $dump_file") unless defined($config->{quiet});
    $DB::single = 1;
	
    # storable
    open my $fh, "| gzip -9 -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
    nstore_fd($tr_cache, $fh);
    close $fh;
}


# DEBUG AND STATUS METHODS
##########################

# gets time
sub get_time() {
    my @time = localtime(time());

    # increment the month (Jan = 0)
    $time[4]++;

    # add leading zeroes as required
    for my $i(0..4) {
        $time[$i] = "0".$time[$i] if $time[$i] < 10;
    }

    # put the components together in a string
    my $time =
        ($time[5] + 1900)."-".
        $time[4]."-".
        $time[3]." ".
        $time[2].":".
        $time[1].":".
        $time[0];

    return $time;
}

# prints debug output with time
sub debug {
    my $text = (@_ ? (join "", @_) : "No message");
    my $time = get_time;
    
    print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}



sub usage {
    my $usage =<<END;
#----------------------------#
# GTF to VEP cache converter #
#----------------------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

This script creates a VEP cache from a GTF file containing transcript/exon
definitions and a FASTA file containing the reference sequence for the same
species.

Usage: perl gtf2vep.pl [arguments]

Options
=======

-h | --help               Display this message and quit
-i | --input [file]       GTF files (may be gzipped)
-f | --fasta [file]       FASTA file or directory containing FASTA files
-s | --species [species]  Species name
-d | --db_version [n]     Database version - must match version of API in use
--dir [dir]               Root directory for cache (default = '\$HOME/.vep/')

END

    print $usage;
}
