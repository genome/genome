package Genome::Model::Tools::ViromeEvent::BlastX_NT::CheckParseOutput;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;
use File::Temp;
use Bio::SeqIO;
use Bio::SearchIO;
use Data::Dumper;
use File::Basename;

class Genome::Model::Tools::ViromeEvent::BlastX_NT::CheckParseOutput{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    return "gzhao's BlastX NT check parse output";
}

sub help_detail {
    return <<"EOS"
This script will check all .tblastx.parsed file in the .BNfiltered_TBLASTX
subdirectory of the given directory to make sure parsing blastn output 
file is finished for each file.
EOS
}

sub execute {
    my $self = shift;

    my $dir = $self->dir;
    my $sample_name = basename($dir);

    $self->log_event("Checking NT blastX parse for $sample_name");

    my $blast_dir = $dir.'/'.$sample_name.'.BNFiltered_TBLASTX_nt';

    unless (-d $blast_dir) {
	$self->log_event("Failed to find NT blastX dir for $sample_name");
	return;
    }

    my @fa_files = glob("$blast_dir/*fa"); #JUST TO GET ROOT FILE NAME TO DERIVE OTHER FILE NAMES
    unless (scalar @fa_files > 0) {
	if (-s $dir.'/'.$sample_name.'.BNfiltered.fa' > 0) {
	    $self->log_event("Failed to find any NT BlastX out files for $sample_name");
	    return;
	}
	elsif (-e $dir.'/'.$sample_name.'.BNfiltered.fa') {
	    $self->log_event("No NT BlastX out files available for parsing for $sample_name");
	    return 1;
	}
	else {
	    $self->log_event("Failed to find any NT BlastX out files for $sample_name");
	    return;
	}
    }

    foreach my $fa_file (@fa_files) {
	next if $fa_file =~ /TBXNTfiltered\.fa$/; #ALREADY FILTERED FILES
	next if $fa_file =~ /TBXNThits\.fa$/; #FASTA OF BLAST HITS
	my $root_name = $fa_file;
	$root_name =~ s/\.fa$//;

	#BLAST OUTPUT FILE
	my $blast_out_file = $root_name.'.tblastx.out';
	unless (-s $blast_out_file) {
	    $self->log_event("Can not find NT blast out file for ".basename($fa_file));
	    return;
	}
	#BLAST OUT FILE PARSED
	my $blast_parse_file = $root_name.'.tblastx.parsed';
	#BLAST FILTERED FILE
	my $blast_filtered_file = $root_name.'.TBXNTfiltered.fa';

	unless (-s $blast_parse_file) {
	    $self->log_event("Running parse NT blastX for ".basename($blast_out_file));
	    unless ($self->run_parser($blast_out_file)) {
		$self->log_event("Parsing failed for ".basename($blast_out_file));
		return;
	    }
	}
	else {
	    unless ($self->check_completed_parse($blast_parse_file)) {
		unless($self->run_parser($blast_out_file)) {
		    $self->log_event("Parsing failed for ".basename($blast_out_file));
		    return;
		}
	    }
	    $self->log_event("NT blastX parse already completed for ".basename($blast_out_file));
	}
    }

    $self->log_event("Completed checking NT blastX parse for $sample_name");

    return 1;
}

sub check_completed_parse {
    my ($self, $parse_file) = @_;
    my $fh = IO::File->new("< $parse_file") ||
	die "Can not create file handle for $parse_file";
    my $line_count = 0; my $undef_taxon = 0;
    my $saved_seq = 0;  my $total_seq = 0;
    my $parse_completed = 0;
    while (my $line = $fh->getline) {
	$line_count++;
	$undef_taxon++ if $line =~ /undefined\s+taxon/;
	if ($line =~ /#\s+Summary:\s+(\d+)\s+out\s+of\s+(\d+)/) {
	    $saved_seq = $1;
	    $total_seq = $2;
	    $parse_completed = 1; #PARSE COMPLETED RUNNING
	    #FOR NOW DON'T DO ADDITIONAL PHYLOTYPE CHECKS
	    return 1;
	}
    }
    $fh->close;
    
    if ($parse_completed == 0) {
	#PARSING NEVER COMPLETED SO RE-RUN
	return;
    }
    #ANY RUNS BELOW HERE COULD RESULT IN INFINITE LOOP SINCE YOU
    #COULD GET THE SAME BLAST RESULT BACK SO THESE ARE SKIPPED FOR NOW
    my $phylotyped_count = $total_seq - $saved_seq;
    if ($phylotyped_count == 0) {
	return 1;
    }
    if ($phylotyped_count <= $undef_taxon) {
	return;
    }
    if (($line_count -1) == $undef_taxon) {
	return;
    }
    if ($phylotyped_count > $line_count -1) {
	return;
    }
    return 1;
}

sub clean_up_tmp_dir {
    my $path = '/tmp/' . Genome::Sys->username;
    unlink $path . 'nodes', $path . 'parents', $path . 'names2id', $path . 'id2names';
    return 1;
}

sub run_parser {
    my ($self, $blast_out_file) = @_;
    
    my $E_cutoff = 1e-5;

    my @unassigned = (); # query should be kept for further analysis
    my $total_records = 0;

    # create ouput file
    my $parse_out_file = $blast_out_file;
    $parse_out_file =~ s/out$/parsed/;
    my $out_fh = IO::File->new("> $parse_out_file") ||
	die "Can not create file handle for $parse_out_file";

    # get a Taxon from a Bio::DB::Taxonomy object
    my $tax_dir = File::Temp::tempdir (CLEANUP => 1);
    my $dbh_sqlite = DBI->connect("dbi:SQLite:/gscmnt/sata835/info/medseq/virome/taxonomy_db");
    my $dbh = Bio::DB::Taxonomy->new(-source => 'flatfile',
				 -directory=> "$tax_dir",
				 -nodesfile=> '/gscmnt/sata835/info/medseq/virome/taxonomy/nodes.dmp',
				 -namesfile=> '/gscmnt/sata835/info/medseq/virome/taxonomy/names.dmp',);
    my $report = new Bio::SearchIO(-format => 'blast', -file => $blast_out_file, -report_type => 'tblastx');

    unless ($report) {
	$self->log_event("Failed to create Bio SearchIO to parse ".basename($blast_out_file));
	return;
    }

    $out_fh->print("QueryName\tQueryLen\tAssignment\tlineage\tHit\tSignificance\n");

    # Go through BLAST reports one by one      
    while(my $result = $report->next_result) {# next query output
	$total_records++;
	my $haveHit = 0;
	my $have_significant_hit = 0;
	my %assignment = ();
	my $assigned = 0;

	# only take the best hits
	my $hit_count = 0;
	my $determined = 0;

	my $bits_cutoff = 45;
	my $best_bit_value = 0;
	my $highest_bit_value = 0;
	while(my $hit = $result->next_hit) {
	    my $hsp_count = 0;
	    foreach my $hsp ($hit->hsps) {
		$highest_bit_value = $hsp->bits if $hsp->bits > $highest_bit_value;
		$hsp_count++;
	    }

	    my @temp_arr = split(/\|/, $hit->name); # gi|num|database|accessionNum|
	    my $gi = $temp_arr[1];
	    next if $temp_arr[2] eq 'pdb'; # skip data from pdb database

	    $haveHit = 1;
	    $hit_count++;
	    if ($hit_count == 1) {
		$best_bit_value = $highest_bit_value;
	    }
   
	    # check whether the hit should be kept for further analysis
	    if ($hit->significance <= $E_cutoff){ # similar to known, need Phylotyped

		my $have_significant_hit = 1;
		if ($highest_bit_value == $best_bit_value) {

		    # from gi get taxonomy lineage
		    my $sth = $dbh_sqlite->prepare("SELECT * FROM gi_taxid where gi = $gi");
		    $sth->execute();
		    my $ref = $sth->fetchrow_hashref();

		    $sth->finish();
		    my $taxID = $ref->{'taxid'};
		    if ($taxID) { # some gi don't have record in gi_taxid_nucl, this is for situation that has
			my $taxon_obj = $dbh->get_taxon(-taxonid => $taxID);
			if (!(defined $taxon_obj)) {
			    my $description .= "undefined taxon\t".$hit->name."\t".$hit->significance;
			    $assignment{"other"} = $description;
			}
			else {
			    my $tree_function = Bio::Tree::Tree->new();
			    my @lineage = $tree_function->get_lineage_nodes($taxon_obj);
			    # each lineage node is a Bio::Tree::NodeI object

			    if (scalar @lineage) {
				$determined = 1;
				#$self->PhyloType(\@lineage,$hit, $best_e, $dbh_sqlite, $dbh, \%assignment);
				$self->PhyloType(\@lineage,$hit, $dbh_sqlite, $dbh, \%assignment);
			    }
			}
		    }
		    else { # for situations that gi does not have corresponding taxid
			my $desc = $hit->description."\t".$hit->name."\t".$hit->significance;
			$determined = 1;
			$assignment{"other"} = $desc;
		    }
		}
		else { # significant but does not have the same e value as the first best hit
		    last; # skip the rest significant hits 
		}
	    }
	    else { # E value is not significant enough
		if($determined){ # skip the rest hits that are not significant
		    last;
		}
		else {
		    my $desc = "hit not significant\t".$hit->name."\t".$hit->significance;
		    $assignment{"unassigned"} = $desc;
		    last;
		}
	    }
	} # end with all hits

	if (!$haveHit) {
	    $assignment{"unassigned"} = "no hit";
	}

	# consolidate assignment
	# If a query is assigned a real taxon name and "other" for reason like "other sequences
	# artificial sequences", or no taxon id in taxon database it will be reported only as
	# the real taxon name
	# If a query is assigned both Homo and Primates, it will be reported as Homo only

	my $num_assignment = keys %assignment;
	if ($num_assignment > 1) { # have multiple assignment
	    # handle the situation that assigned both "Homo" and "other"
	    my $has_specific = 0;
	    my $has_other = 0;
	    if ((defined $assignment{"Bacteria"}) || 
		(defined $assignment{"Fungi"}) || 
		(defined $assignment{"Homo"}) || 
		(defined $assignment{"Mus"}) || 
		(defined $assignment{"Phage"}) || 
		(defined $assignment{"Viruses"})) {
		$has_specific = 1;
	    }
	    if (defined $assignment{"other"}) {
		$has_other = 1;
	    }
	    if ($has_specific && $has_other) {
		delete $assignment{"other"};
	    }
	}
	# print out assignment
	foreach my $assign (keys %assignment) {
	    if ($assign eq "unassigned") {
		push @unassigned, $result->query_name;
	    }
	    else {
		$out_fh->print($result->query_name, "\t", $result->query_length, "\t", $assign, "\t", $assignment{$assign}, "\n");
	    }
	}
    } # end of report parsing
    $out_fh->print("# Summary: ", scalar @unassigned, " out of $total_records ", (scalar @unassigned)*100/$total_records, "% are unassigned.\n");
    $out_fh->close;
    $dbh_sqlite->disconnect();

    # generate a fasta file that contains all the sequences that do not match
    # to known sequences read in tblastx input sequences
    # read in tblastx input sequences

    my $root_file_name = $blast_out_file;
    $root_file_name =~ s/\.tblastx\.out$//;
    my $orig_fasta = $root_file_name.'.fa';
    my $filter_out_file = $root_file_name.'.TBXNTfiltered.fa';
    my $hits_file = $root_file_name.'.TBXNThits.fa';

    my $in = Bio::SeqIO->new(-format => 'fasta', -file => $orig_fasta);
    my $out = Bio::SeqIO->new(-format => 'fasta', -file => ">$filter_out_file");
    my $hits = Bio::SeqIO->new(-format => 'fasta', -file => ">$hits_file");
    unless ($in and $out) {
	$self->log_event("Failed to create SeqIO to read/write sequence for ".basename($blast_out_file));
	return;
    }
    while (my $seq = $in->next_seq) {
	my $read_name = $seq->primary_id;
	if (grep (/^$read_name$/, @unassigned)) {
	    $out->write_seq($seq);
	}
        else
        {
            $hits->write_seq($seq);
        }
    }

    return 1;
}
		
sub PhyloType {
    my ($self, $lineage_ref, $hit_ref, $dbh_sqlite, $dbh_taxonomy, $assignment_ref) = @_;
    my $description = "";
    my $node_id; 
    my $obj;
    my $name;
    my $assigned = 0;
    
    my $Lineage = "";
    for (my $i = 0; $i <= $#$lineage_ref; $i++) { 
	my $temp_node_id = $lineage_ref->[$i]->id;
	my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid=>$temp_node_id);
	my $temp_name = $temp_obj->scientific_name;
	$Lineage .= $temp_name.";";
    }					

    # check to see if it is a human sequence
    if (scalar @{$lineage_ref} >= 4) {
	$node_id = $lineage_ref->[3]->id;
	$obj = $dbh_taxonomy->get_taxon(-taxonid=>$node_id);
	$name = $obj->scientific_name;

	#looks like if assigned 'Homo' don't do rest of classification??
	if ($name eq "Metazoa") {
	    # make assignment
	    if (not exists $assignment_ref->{'Homo'}) {
		for (my $i = 0; $i <= $#$lineage_ref; $i++) { 
		    my $temp_node_id = $lineage_ref->[$i]->id;
		    my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid=>$temp_node_id);
		    my $temp_name = $temp_obj->scientific_name;
		    
		    if ($temp_name eq "Homo") {
			$description = "Homo\t".$hit_ref->name."\t".$hit_ref->significance;
			$assignment_ref->{"Homo"} = $description;
			$assigned = 1;
			last;
		    }
		}
	    }
	    if (!$assigned and not exists $assignment_ref->{'Mus'}) {
		for (my $i = 0; $i <= $#$lineage_ref; $i++) { 
		    my $temp_node_id = $lineage_ref->[$i]->id;
		    my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid=>$temp_node_id);
		    my $temp_name = $temp_obj->scientific_name;

		    if ($temp_name eq "Mus") {
			$description = "Mus\t".$hit_ref->name."\t".$hit_ref->significance;
			$assignment_ref->{"Mus"} = $description;
			$assigned = 1;
			last;
		    }
		}
	    }
	    if (!$assigned and not exists $assignment_ref->{'other'}) {
		$description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
		$assignment_ref->{"other"} = $description;
		$assigned = 1;
	    }
	}
    }
    
    # check to see if it is bacteria sequence
    if ((scalar @{$lineage_ref} >= 2) && (!$assigned) and not exists $assignment_ref->{'Bacteria'}) {
	$node_id = $lineage_ref->[1]->id;
	$obj = $dbh_taxonomy->get_taxon(-taxonid=>$node_id);
	$name = $obj->scientific_name;

	if ($name eq "Bacteria") {
	    $description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
	    $assignment_ref->{"Bacteria"} = $description;
	    $assigned = 1;
	}
    }

    # check to see if it is a phage virus sequence
    if (!$assigned and not exists $assignment_ref->{'Phage'}) {
	$node_id = $lineage_ref->[0]->id;
	$obj = $dbh_taxonomy->get_taxon(-taxonid=>$node_id);
	$name = $obj->scientific_name;

	if ($name eq "Viruses") {
	    for (my $i = 0; $i <= $#$lineage_ref; $i++) { 
		my $temp_node_id = $lineage_ref->[$i]->id;
		my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid=>$temp_node_id);
		my $temp_name = $temp_obj->scientific_name;
		if (($temp_name eq "Lipothrixviridae")||($temp_name eq "Caudovirales")||($temp_name eq "Corticoviridae")||($temp_name eq "Cystoviridae")||($temp_name eq "Inoviridae")||($temp_name eq "Leviviridae")||($temp_name eq "Microviridae")||($temp_name eq "Tectiviridae")||($temp_name =~ /phage/i)) {
		    $description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
		    $assignment_ref->{"Phage"} = $description;
		    $assigned = 1;
		    last;
		}
	    }
	}
    }
    
    # check to see if it is a virus sequence
    if (!$assigned and not exists $assignment_ref->{'Viruses'}) {
	$node_id = $lineage_ref->[0]->id;
	$obj = $dbh_taxonomy->get_taxon(-taxonid=>$node_id);
	$name = $obj->scientific_name;

	if ($name eq "Viruses") {
	    $description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
	    $assignment_ref->{"Viruses"} = $description;
	    $assigned = 1;
	}
    }
    
    # check to see if it is a fungi sequence
    if ((scalar @{$lineage_ref} >= 4) && (!$assigned) and not exists $assignment_ref->{'Fungi'}) {
	$node_id = $lineage_ref->[3]->id;
	$obj = $dbh_taxonomy->get_taxon(-taxonid=>$node_id);
	$name = $obj->scientific_name;

	if ($name eq "Fungi") {
	    $description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
	    $assignment_ref->{"Fungi"} = $description;
	    $assigned = 1;
	}
    }
    
    if (!$assigned and not exists $assignment_ref->{'other'}) {
	$description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
	$assignment_ref->{"other"} = $description;
	$assigned = 1;
    }
    return $assigned;
}

1;
