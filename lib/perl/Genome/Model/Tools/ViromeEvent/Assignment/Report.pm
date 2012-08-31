package Genome::Model::Tools::ViromeEvent::Assignment::Report;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;
use Switch;
use Bio::SeqIO;
use Bio::SearchIO;
use File::Basename;
use Data::Dumper;

class Genome::Model::Tools::ViromeEvent::Assignment::Report{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    "gzhao's reporting for virome"
}

sub help_detail {
    return <<"EOS"
This script will read corresponding files in the given director and 
generate a report. It will report in each library, for each category,
how many total sequence were assigned to this category, how many were 
assigned by BLASTN, how many were assigned by TBLASTX, the range of 
percent identity. It will also generate four fasta format files which 
contain viral reads from blastn, tblastx, all viral reads and reads
that can not be assigned to any category.
EOS
}

sub execute {
    my $self = shift;
    my $lib_name = basename($self->dir);
    $self->log_event("Assignment Reporting starting for $lib_name");
    #TABLE TO CONVERT TO CORRECT BLAST NAME AND
    #CORRISPONDING FILE EXTENSIONS FOR EACH BLAST
    my $blasts = $self->_blasts_and_file_exts();
    #REPORT WHAT ORGANISMS WHERE HIT IN EACH BLAST RUNS
    my $all_orgs_hit = {};
    my $viral_lineage_hits = {};
    my @unassigned_reads;
    my @all_virus_lineages;
    #BN_HG  = BlastN HumanGenomic    #BN     = BlastN NT
    #TBX_nt = BlastX NT              #TBX_VG = BlastX Viral
    foreach my $blast_type (qw/ BN_HG BN TBX_nt TBX_VG /) {
	my $blast_dir = $self->dir.'/'.$lib_name . $blasts->{$blast_type}->{'dir_ext'};
	unless (-d $blast_dir) {
	    $self->log_event("Failed to find blast dir: ".basename($blast_dir));
	    return;
	}
	my @parse_files = glob("$blast_dir/*parsed");
	unless (@parse_files) {
	    $self->log_event("NO blast parse files to read in ".basename($blast_dir));
	    next;
	}
	foreach my $file (@parse_files) {
	    my $hits = $self->get_organisms_hit($file);

	    foreach my $hit (keys %$hits) {
		next if $hit eq 'virus_lineage';
		next if $hit eq 'unassigned_reads';
		$all_orgs_hit->{$blast_type}->{$hit} += $hits->{$hit};
		$all_orgs_hit->{'total'}->{$hit} += $hits->{$hit};
	    }
	    #KEEP A LIST OF READS TYPED BY LINEAGE FOR EACH BLAST EVENT
	    if (exists $hits->{'virus_lineage'}) {
		foreach my $vl (keys %{$hits->{'virus_lineage'}}) {
		    push @all_virus_lineages, $vl unless grep (/^$vl$/, @all_virus_lineages);
		    foreach my $read_name (keys %{$hits->{'virus_lineage'}->{$vl}}) {
			#READ_NAME => 'E-VALUE'
			$viral_lineage_hits->{$blast_type}->{$vl}->{$read_name} = $hits->{'virus_lineage'}->{$vl}->{$read_name};
		    }
		}
	    }
	    if (exists $hits->{'unassigned_reads'}) {
		push @unassigned_reads, map {$_} @{$hits->{'unassigned_reads'}};
	    }
	}
    }

    my $report_file = $self->dir.'/'.$lib_name.'.AssignmentReport';
    my $report_fh = IO::File->new("> $report_file") ||
	die "Can not create file handle for ".basename($report_file);
    #PRINT TABLE .. COLUMN NAMES: TOTAL BN_HG BN TBX_NT TBX_VG
    my $header;
    foreach my $col (qw/ category total BN_HG BN TBX_nt TBX_VG /) {
	$header .= sprintf "%12s", $col;
    }
    $report_fh->print($self->dir."\n".$header."\n");
    my $string;
    foreach my $org ($self->_valid_organism_category()) {
	$string .= sprintf "%12s", $org; #ORGANISM NAMES
	foreach my $bl (qw/ total BN_HG BN TBX_nt TBX_VG /) {
	    my $hit = (exists $all_orgs_hit->{$bl}->{lc $org}) ? $all_orgs_hit->{$bl}->{lc $org} : 0;
	    $string .= sprintf "%12s", $hit;
	}
	$report_fh->print("$string\n");
	$string = '';
    }
    $report_fh->print("\n###########################################################\n\n");

    my $detailed_read_info = {};
    #PRINT VIRUS LINEAGE INFO .. HAVE TO GO THROUGH EACH BLAST FILE OF EACH BLAST TYPE
    foreach my $bl (qw/ BN TBX_nt TBX_VG /) {
	if (exists $viral_lineage_hits->{$bl}) {
	    #GET A LIST OF BLAST FILES
	    my $blast_dir = $self->dir.'/'.$lib_name . $blasts->{$bl}->{'dir_ext'};
	    my @bl_out_files = glob("$blast_dir/*out");
	    #THERE MUST BE OUT FILES .. 
	    unless (@bl_out_files) {
		$self->log_event("NO blast out files available in ".basename($blast_dir));
		next;
	    }
	    #RE-ARRANGE DATA .. MAKE SURE READ SET CORRISPONDS WITH BLAST OUT FILES
	    my $reads = {};
	    #sorry this is awfully hard to follow
	    foreach my $lineage (keys %{$viral_lineage_hits->{$bl}}) {
		foreach my $read (keys %{$viral_lineage_hits->{$bl}->{$lineage}}) {
		    $reads->{$read} = $viral_lineage_hits->{$bl}->{$lineage}->{$read};# = evalue
		}
	    }
	    my $type = $blasts->{$bl}->{'engine'};
	    my $info = $self->_detailed_virus_info($reads, \@bl_out_files, $type);
	    foreach my $read (keys %$info) {
		$detailed_read_info->{$read} = $info->{$read}; #BAD
	    }
	}
    }
    #print Dumper $detailed_read_info;
    foreach my $lineage (@all_virus_lineages) {
	my $read_count = 0;
	#FIRST PASS TO GET NUMBER OF READS
	foreach my $bl (sort keys %$viral_lineage_hits) {
	    if (exists $viral_lineage_hits->{$bl}->{$lineage}) {
		$read_count += scalar (keys %{$viral_lineage_hits->{$bl}->{$lineage}});
	    }
	}
	#SECOND PASS PRINT DATA
	$report_fh->print("$lineage\ttotal number of reads $read_count\n\n".
	                  "QueryName\tQuerylength\t\tHitName\t\tHitLen\t\t\tHitDesc\t\t\tAlnLen\t%ID\tHitStart\tHitEnd\te\n");
	foreach my $bl (sort keys %$viral_lineage_hits) {
	    if (exists $viral_lineage_hits->{$bl}->{$lineage}) {
		$report_fh->print("reads from ".$blasts->{$bl}->{'engine'}.":\n");
		foreach my $read_name (keys %{$viral_lineage_hits->{$bl}->{$lineage}}) {
		    $self->log_event("Failed to find detailed info for $read_name") unless
			exists $detailed_read_info->{$read_name};
		    $report_fh->print(map {$_."\n"} @{$detailed_read_info->{$read_name}});
		}
	    }
	}
	$report_fh->print("\n###########################################################\n\n");
    }
    $report_fh->close;

    #DEFINE OUTPUT FILES .. WE WANT THESE TO BE THERE EVEN WITH ZERO SIZE
    #VIRAL READS OUT FILE
    my $viral_fasta_out = $self->dir.'/'.$lib_name.'.ViralReads_all.fa';
    my $v_out = Bio::SeqIO->new(-format => 'fasta', -file => ">$viral_fasta_out");
    #UNASSIGNED READS OUT FILE
    my $unassigned_out = $self->dir.'/'.$lib_name.'.unassigned.fa';
    my $unas_out = Bio::SeqIO->new(-format => 'fasta', -file => ">$unassigned_out");

    #CREATE VIRAL READS FASTA AND UNASSIGNED READS FASTA
    my $fasta_file = $self->dir.'/'.$lib_name.'.fa.cdhit_out.masked.goodSeq';
    if (-s $fasta_file) {
	my $in = Bio::SeqIO->new(-format => 'fasta', -file => $fasta_file);
	while (my $seq = $in->next_seq) {
	    my $read_name = $seq->primary_id;
	    if (grep (/^$read_name$/, @unassigned_reads)) {
		$unas_out->write_seq($seq);
	    }
	    if (exists $detailed_read_info->{$read_name}) {
		$v_out->write_seq($seq);
	    }
	}
    }
    elsif (-e $fasta_file) {
	$self->log_event("No data available for analysis in ".basename($lib_name));
	#return 1;
    }
    else {
	$self->log_event("Failed to find repeatMasker goodSeq file");
	return;
    }

    $self->log_event("No viral reads found for $lib_name") unless -s $viral_fasta_out;
    $self->log_event("No unassigned reads found for $lib_name") unless -s $unassigned_out;

    $self->log_event("Assignment Reporting completed for $lib_name");
    return 1;
}

sub get_organisms_hit {
    my $self = shift;
    my $parse_file = shift;
    my $org_hits = {};
    my @valid_orgs = $self->_valid_organism_category();
    my $fh = IO::File->new("< $parse_file") ||
	die "Can not create file handle for $parse_file";
    while (my $line = $fh->getline) {
	next if $line =~ /#/ || $line =~ /^\s+$/ || $line =~ /^QueryName/;
	my @tmp = split (/\t+/, $line);
	#$tmp[0] = read name
	#$tmp[2] = category, ie, organism
	#$tmp[3] = organism lineage
	#$tmp[5] = e-value
	unless (grep (/^$tmp[2]$/i, @valid_orgs)) {
	    $self->log_event("Invalid organism category: $line");
	    return;
	}
	$org_hits->{lc $tmp[2]}++;
	if ($tmp[2] =~ /viruses/i) {
	    #TRACK VIRUS HITS BY LINEAGE $tmp[3]
	    $org_hits->{'virus_lineage'}->{$tmp[3]}->{$tmp[0]} = $tmp[5];
	}
	#ARRAY OF UNASSIGNED READS
	if ($tmp[2] =~ /unassigned/i) {
	    push @{$org_hits->{'unassigned_reads'}}, $tmp[0];
	}
    }
    $fh->close;
    return $org_hits;
}

sub _valid_organism_category {
    return qw/ Bacteria Fungi Homo Mus Phage Viruses other unassigned /;
}

sub _blasts_and_file_exts {
    return {
	'BN_HG' => {
	    'name'    => 'BlastHumanGenomic',
	    'engine'  => 'blastn',
	    'dir_ext' => '.fa.cdhit_out.masked.goodSeq_HGblast',
	},
	'BN' => {
	    'name'    => 'NTBlastN',
	    'engine'  => 'blastn',
	    'dir_ext' => '.HGfiltered_BLASTN',
	},
	'TBX_nt' => {
	    'name'    => 'NTBlastX',
	    'engine'  => 'tblastx',
	    'dir_ext' => '.BNFiltered_TBLASTX_nt',
	},
	'TBX_VG' => {
	    'name'    => 'ViralBlastX',
	    'engine'  => 'tblastx',
	    'dir_ext' => '.TBXNTFiltered_TBLASTX_ViralGenome',
	},
    };
}

sub _detailed_virus_info {
    my $self = shift;
    my $reads = shift;
    my $blast_files = shift;
    my $report_type = shift;
    my $info = {};
    foreach my $file (@$blast_files) {
	my $report = Bio::SearchIO->new(-format => 'blast', -file => $file, -report_type => $report_type);
	while (my $result = $report->next_result) {
	    if (exists $reads->{$result->query_name}) {
		while (my $hit = $result->next_hit) {
		    next unless $hit->significance == $reads->{$result->query_name};
		    my $desc = $result->query_name."\t".$result->query_length."\t".
			    $hit->name."\t".$hit->length."\t".$hit->description."\t";
		    my $hsp = $hit->next_hsp;
		    next unless $hsp; #weird ..next_hsp can return undef hsp ??
		    $desc .= $hsp->length('hit')."\t"; #sometimes this is not defined
		    $desc .= sprintf("%4.1f", $hsp->percent_identity)."\t";
		    $desc .= $hsp->start('hit')."\t".$hsp->end('hit')."\t".$hsp->evalue;
		    push @{$info->{$result->query_name}}, $desc;
		}
	    }
	}
    }
    return $info;
}

1;
