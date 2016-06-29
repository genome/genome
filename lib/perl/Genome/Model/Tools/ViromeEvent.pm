package Genome::Model::Tools::ViromeEvent;

use strict;
use warnings;

use Genome;
use IO::File;
use File::Basename;

use Bio::SeqIO;
use Bio::SearchIO;

class Genome::Model::Tools::ViromeEvent{
    is => 'Command',
    is_abstract => 1,
    has => [
	dir => {
	    doc => 'directory of inputs',
	    is => 'String',
	    is_input => 1,
	    is_optional => 1,
	},            
	logfile => {
	    is => 'String',
	    doc => 'output file for monitoring progress of pipeline',
	    is_input => 1,
	},
        human_db => {
            is => 'String',
            doc => 'human blast db',
            is_optional => 1,
            is_input => 1,
        },
        nt_db => {
            is => 'String',
            doc => 'nt blast db',
            is_optional => 1,
            is_input => 1,
        },
        virus_db => {
            is => 'String',
            doc => 'virus blast db',
            is_optional => 1,
            is_input => 1,
        },
        taxonomy_db => {
            is => 'String',
            doc => 'taxonomy db',
            is_optional => 1,
            is_input => 1,
        },
    ],
};

sub help_brief {
    "skeleton for gzhao's virome script"
}

sub help_detail {
    'wrapper for script sequence to be utilized by workflow';
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return $self;
}

sub execute {
    die("abstract");
}

sub log_event {
    my ($self,$message) = @_;
    my ($event) = $self->class =~ /ViromeEvent::(\S+)$/;
    my $fh = Genome::Sys->open_file_for_appending( $self->logfile );
    my $time = localtime(time);
    $fh->printf("%-30s%-45s%10s\n", $time, $event, $message);
    $fh->close();
}

my %file_extensions = (
    'hg_blast' => {
	blast_dir_ext => 'fa.cdhit_out.masked.goodSeq_HGblast',
        pooled_out_file_ext => 'fa.cdhit_out.masked.goodSeq',
	prev_pooled_file_ext => 'fa.cdhit_out.masked.goodSeq',
    },
    'blastn' => {
	blast_dir_ext       => 'HGfiltered_BLASTN',                  # dir to run blastn in
	pooled_out_file_ext => 'HGfiltered.fa',                      # file ext for prev blast filtered reads
	prev_blast_dir_ext  => 'fa.cdhit_out.masked.goodSeq_HGblast',# prev blast dir to get filtered reads from
	split_file_ext      => 'HGfiltered.fa_file',                 # file ext for files to run blastn on
    },
    'blastx_nt' => {
	blast_dir_ext => 'BNFiltered_TBLASTX_nt', #<--- same
	pooled_out_file_ext => 'BNfiltered.fa',
	prev_blast_dir_ext => 'HGfiltered_BLASTN',
	split_file_ext => 'BNFiltered.fa_file',
    },
    'blastx_viral' => {
	blast_dir_ext => 'TBXNTFiltered_TBLASTX_ViralGenome', #<--- same
	pooled_out_file_ext => 'TBXNTfiltered.fa',
	prev_blast_dir_ext => 'BNFiltered_TBLASTX_nt',
	split_file_ext => 'TBXNTFiltered.fa_file',
    },
);

sub pool_and_split_sequence {
    my ( $self, $stage, $read_limit ) = @_;

    my $dir = $self->dir;
    my $sample_name = basename($dir);

    #create a blast dir
    my $blast_dir = $dir.'/'.$sample_name.'.'.$file_extensions{$stage}{blast_dir_ext};
    Genome::Sys->create_directory( $blast_dir ) unless -d $blast_dir;
  
    #define a file to pool blast filtered files from previous stage
    my $pooled_file = $dir.'/'.$sample_name.'.'.$file_extensions{$stage}{pooled_out_file_ext};
    my $out = Bio::SeqIO->new(-format => 'fasta', -file => ">$pooled_file");

    #find previous stage blast dir
    my $prev_blast_dir = $dir.'/'.$sample_name.'.'.$file_extensions{$stage}{prev_blast_dir_ext};

    unless (-d $prev_blast_dir) {
	$self->log_event("Failed to fine previous blast stage dir for stage: $stage, sample name: $sample_name.  Expected $prev_blast_dir");
	return;
    }

    #check blast results from previous stage
    my @prev_bl_files = glob("$prev_blast_dir/*fa");
    if (@prev_bl_files == 0) {
	$self->log_event("No further reads available for $stage for $sample_name");
	return 1;
    }

    #find blast filtered files .. if blast input files are present and blast runs
    #output file should get created even if empty
    my $glob_file_ext = $file_extensions{$stage}{pooled_out_file_ext};
    my @filtered_files = glob("$prev_blast_dir/*$glob_file_ext");

    unless (scalar @filtered_files > 0) {
	$self->log_event("Failed to find any $stage filtered data for $sample_name");
	return;
    }

    #pool blast filtered files into a single output file
    foreach my $file (@filtered_files) {
	my $in = Bio::SeqIO->new(-format => 'fasta', -file => $file);
	while (my $seq = $in->next_seq) {
	    $out->write_seq($seq);
	}
    }
    unless (-e $pooled_file) {
	$self->log_event("Failed to create pooled file of $stage filtered reads for sample: $sample_name");
	return;
    }

    if ( not -s $pooled_file ) {
        $self->log_event("No reads available for further processing .. all reads appear to have been filtered out in stage $stage for sample $sample_name");
        return 1;
    }

    #split files
    my $c = 0; my $n = 0; my $limit = $read_limit;
    my $in = Bio::SeqIO->new(-format => 'fasta', -file => $pooled_file);
    my $path = $blast_dir.'/'.$sample_name.'.'.$file_extensions{$stage}{split_file_ext};
    my $split_file = $path.$n.'.fa';

    my $split_out = Bio::SeqIO->new(-format => 'fasta', -file => ">$split_file");
    while (my $seq = $in->next_seq) {
	$c++;
	$split_out->write_seq($seq);
	if ($c == $limit) {
	    $c = 0;
	    my $split_file = $path.++$n.'.fa';
	    $split_out = Bio::SeqIO->new(-format => 'fasta', -file => ">$split_file");
	}
    }

    $self->log_event("Pooled data to run $stage completed for $sample_name");

    #clean up empty files which can be made if total # reads if a multiple of read limit
    for my $fa_file ( glob( "$path*" ) ) {
        unlink $fa_file if not -s $fa_file;
    }

    return 1;
}

sub get_files_for_blast {
    my ( $self, $stage ) = @_;
    
    my $dir = $self->dir;
    my $sample_name = basename ($dir);

    #check to make sure blast dir exists
    my $blast_dir = $dir.'/'.$sample_name.'.'.$file_extensions{$stage}{blast_dir_ext};
    unless ( -d $blast_dir ) {
	$self->log_event("Failed to find $stage blast directory for sample name: $sample_name");
	return; #die
    }

    my @files_for_blast;
    #exclude fa files not for blasting
    for my $file ( glob("$blast_dir/$sample_name*fa") ) {
	next if $file =~ /filtered\.fa$/; #skip blast out files
	next if $file =~ /hits\.fa$/;     #skip parsed file
	push @files_for_blast, $file;
    }
    # all reads filtered out in prev stage
    if ( not @files_for_blast ) {
        $self->log_event("No reads available to blast at stage, $stage, for sample_name, $sample_name");
        #place holder .. otherwise workflow will not proceed to next stage
        $self->files_for_blast( ['no_files_to_blast'] );
	return 1;
    }
    $self->files_for_blast( \@files_for_blast );

    $self->log_event("Completed for sample: $sample_name");

    return 1;
}

sub _blast_params_for_stage {
    my( $self, $stage ) = @_;
    my %p = (
        hg_blast => {
            blast_db => $self->human_db,
            blast_cmd => 'blastall -p blastn -e 1e-8 -I T -b 2',
            out_file_ext => 'HGblast.out',
        },
        blast_n => {
            blast_db => $self->nt_db,
            blast_cmd => 'blastall -p blastn -e 1e-8 -I T',
            out_file_ext => 'blastn.out',
        },
        blastx_nt => {
            blast_db => $self->nt_db,
            blast_cmd => 'blastall -p tblastx -e 1e-2 -I T',
            out_file_ext => 'tblastx.out',
        },
        blastx_viral => {
            blast_db => $self->virus_db,
            blast_cmd => 'blastall -p tblastx -e 0.1 -I T',
            out_file_ext => 'tblastx_ViralGenome.out',
        }
    );
    return if not exists $p{$stage};
    return $p{$stage};
}

sub run_blast_for_stage {
    my ( $self, $stage ) = @_;

    my $blast_params = $self->_blast_params_for_stage( $stage );
    if ( not $blast_params ) {
        $self->log_event("Failed to get blast params for stage: $stage");
        return;
    }

    if ( $self->file_to_run eq 'no_files_to_blast' ) {
        #just place holding to move to next stage
        $self->log_event("No reads available to blast at stage: $stage");
        return 1;
    }

    my $input_file = $self->file_to_run;
    my $input_file_name = File::Basename::basename( $input_file );

    my $blast_done_file_ext = $blast_params->{out_file_ext};
    my $blast_out_file = $input_file;
    $blast_out_file =~ s/fa$/$blast_done_file_ext/;
    
    if (-s $blast_out_file) {
	my $tail = `tail -n 50 $blast_out_file`;
	if ($tail =~ /Matrix/) {
	    $self->log_event("$stage already ran for $input_file_name");
	    return 1;
	}
    }

    #$self->log_event( "Running $stage for $input_file_name" );

    my $blast_cmd = $blast_params->{blast_cmd};
    if ( not $blast_cmd ) {
        $self->log_event("Failed to get blast command for stage: $stage");
        return;
    }
    my $blast_db = $blast_params->{blast_db};
    if ( not $blast_db ) {
        $self->log_event("Failed to get blast db for stage: $stage");
        return;
    }
    my $cmd = $blast_cmd.' -i '.$input_file.' -o '.$blast_out_file.' -d '.$blast_db;
    #$self->log_event("Running blast with command: $cmd"); # for test

    if (system ($cmd)) {
	$self->log_event("$stage failed for $input_file_name, cmd: $cmd");
	return;
    }

    $self->log_event("Blast completed for $input_file_name");

    return 1;
}

sub versioned_taxonomy_dir {
    my $self = shift;

    if( -d $self->taxonomy_db ) {
        return $self->taxonomy_db;
    }
    else {
        my ($version) = $self->taxonomy_db =~ /taxonomy_db_(\S+)$/;
        $self->log_event('Could not derive taxonomy version') and return
            if not $version;
        my $version_dir = "/gscmnt/sata835/info/medseq/virome/taxonomy_$version";
        $self->log_event("Could not get taxonomy dir for version: $version") and return
            if not -d $version_dir;

        return $version_dir
    }
}

sub get_taxids_for_gis {
    my($self, $gis) = @_;

    my %taxids;
    return \%taxids if not %$gis;

    my $taxonomy_dir = $self->versioned_taxonomy_dir;

    my $gi_file = $taxonomy_dir.'/gi_taxid_nucl.dmp';
    $self->log_event('Could not get gi file for taxonomy version $version') and return
        if not -s $gi_file;
        
    my $fh = Genome::Sys->open_file_for_reading($gi_file);
    while ( my $line = $fh->getline ) {
        my @tmp = split( /\s+/, $line);
        #tmp[0] = gi
        #tmp[1] = taxid
        next if not exists $gis->{$tmp[0]};
        $taxids{$tmp[0]} = $tmp[1];
        delete $gis->{$tmp[0]};
        last if not %$gis;
    }
    $fh->close;
    return \%taxids;
}

sub taxonomy_nodes_file {
    my $self = shift;

    my $file = $self->versioned_taxonomy_dir.'/nodes.dmp';

    $self->log_event("Failed to find taxnonomy nodes file: $file") and return
        if not -s $file;

    return $file;
}

sub taxonomy_names_file {
    my $self = shift;

    my $file = $self->versioned_taxonomy_dir.'/names.dmp';

    $self->log_event("Failed to find taxonomy names file: $file") and return
        if not -s $file;

    return $file;
}

sub taxon_db {
    my $self = shift;

    return $self->{taxon_db} if $self->{taxon_db};

    my $tax_dir = Genome::Sys->create_temp_directory();

    my $taxon_db = Bio::DB::Taxonomy->new(
        -source => 'flatfile',
        -directory => $tax_dir,
        -nodesfile => $self->taxonomy_nodes_file,
        -namesfile => $self->taxonomy_names_file,
    );
    return if not $taxon_db;

    $self->{taxon_db} = $taxon_db;

    return $self->{taxon_db};
}

sub get_blast_report {
    my ( $self, %p ) = @_;

    my $file = $p{blast_out_file};
    my $type = $p{blast_type};

    my $report = new Bio::SearchIO(
        -format      => 'blast',
        -file        => $file,
        -report_type => $type,
    );
    
    if ( not $report ) {
        $self->log_event('Failed to create blast report');
        return;
    }
    
    return $report;
}

1;
