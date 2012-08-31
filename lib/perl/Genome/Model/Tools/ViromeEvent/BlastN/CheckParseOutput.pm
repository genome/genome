package Genome::Model::Tools::ViromeEvent::BlastN::CheckParseOutput;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;
use File::Temp;
use File::Basename;
use Bio::SeqIO;
use Bio::SearchIO;

class Genome::Model::Tools::ViromeEvent::BlastN::CheckParseOutput{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    return "gzhao's Blast N check parse output";
}

sub help_detail {
    return <<"EOS"
This script will check all .tblastx.parsed file in the .BNfiltered_TBLASTX
subdirectory of the given directory to make sure parsing blastn output 
file is finished for each file.

perl script <sample dir>
<sample dir> = full path to the folder holding files for a sample library 
	without last "/"
EOS
}

sub execute {
    my $self = shift;

    my $dir = $self->dir;
    my $sample_name = basename($dir);

    $self->log_event("Checking NT blastN parse for $sample_name");

    my $blast_dir = $dir.'/'.$sample_name.'.HGfiltered_BLASTN';
    unless (-d $blast_dir) {
	$self->log_event("Failed to find NT blastN dir for $sample_name");
	return;
    }

    my @fa_files = glob("$blast_dir/*fa");
    unless (scalar @fa_files > 0) {
	if (-s $dir.'/'.$sample_name.'.HGfiltered.fa' > 0) {
	    $self->log_event("Failed to find any NT blastN output files for $sample_name");
	    return;
	}
	elsif (-e $dir.'/'.$sample_name.'.HGfiltered.fa') {
	    $self->log_event("No NT blastN output files available for parsing for $sample_name");
	    return 1;
	}
	else {
	    $self->log_event("Failed to find any NT blastN out files for $sample_name");
	    return;
	}
    }

    foreach my $fa_file (@fa_files) {
	next if $fa_file =~ /BNfiltered\.fa$/; #FILTERED FILES  .. NOT THE ORIG INPUT
	next if $fa_file =~ /BNhits\.fa$/;     #FASTA OF BLAST HITS
	my $root_name = $fa_file;
	$root_name =~ s/\.fa$//;

	#BLAST OUTPUT FILE
	my $blast_out_file = $root_name.'.blastn.out';
	unless (-s $blast_out_file) {
	    $self->log_event("Can not find NT blastN output file for ".basename($fa_file));
	    return;
	}
	#BLAST OUT FILE PARSED
	my $blast_parse_file = $root_name.'.blastn.parsed';
	#BLAST OUT FILTERED OUTPUT FILE
	my $blast_filtered_file = $root_name.'.BNfiltered.fa';

	if (-s $blast_parse_file) {
	    if ($self->check_completed_parse($blast_out_file)) {
		$self->log_event("Parsing already complete for ".basename($blast_out_file));
	    }
	    else {
		if ($self->run_parser($blast_out_file)) {
		    $self->log_event("Parsing completed for ".basename($blast_out_file));
		}
		else {
		    $self->log_event("Parsing failed for ".basename($blast_out_file));
		    return;
		}
	    }
	}
	else {
	    $self->log_event("Running NT blastX parse for ".basename($blast_out_file));
	    if ($self->run_parser($blast_out_file)) {
		$self->log_event("Parsing ran successfully for ".basename($blast_out_file));
	    }
	    else {
		$self->log_event("Parsing failed for ".basename($blast_out_file));
		return;
	    }
	}
    }

    $self->log_event("Completed checking NT blastN parse for $sample_name");

    return 1;
}

# METHODS BELOW STILL NEEDS SOME SIGNIFICANT REFACTORING

sub check_completed_parse {
    my ($self, $parse_file) = @_;
    my $fh = IO::File->new("<$parse_file");
    #DIE STATEMENT SEEM TO GET BURIED IN WORKFLOW??
    $self->log_event("Failed to create file handle for $parse_file") and return
	unless $fh;
    my $line_count = 0; my $undef_taxon = 0; my $parse_completed = 0;
    my $total_seq = 0;  my $saved_seq = 0;
    while (my $line = $fh->getline) {
	$line_count++;
	$undef_taxon++ if $line =~ /undefined\s+taxon/;

	if ($line =~ /\s+Summary:/) {
	    #SKIP THE CHECKS BELOW .. COULD GO INTO INFINITE LOOP
	    return 1;
	    $parse_completed = 1;
	    #$saved_seq = $1;
	    #$total_seq = $2;
	}
    }
    $fh->close;
    if ($parse_completed == 0) {
	#PARSING NEVER COMPLETED
	return 0;
    }

    #CHECKS BELOW ARE NOT DONE TO AVOID POTENTIAL RERUNNING OF PARSE
    #AND GETTING THE SAME RESULT AND RUNNING IT AGAIN
    my $phylotyped_count = $total_seq - $saved_seq;
    if ($phylotyped_count == 0) {
	return 1;
    }
    #RE-RUN PARSING IF MORE THAN HALF RESULTS ARE UNDEFINED TAXON?
    if ($phylotyped_count <= $undef_taxon){
	return 0;
    }
    #RE-RUN BLAST PARSE IF ALL RECORDS ARE UNDEFINED TAXON
    if (($line_count) -1 == $undef_taxon) {
	return 0;
    }
    #???
    # deal with old situation where some read was not recorded because of no 
    # record of gi-taxon record in the database 
    #???
    if ($phylotyped_count > $line_count -1) {
	return 0;
    }
    return 1;
}

sub run_parser {
    my ($self, $blast_out_file) = @_;
    # cutoff value for having a good hit
    my $E_cutoff = 1e-10;

    # create ouput file
    my $parse_out_file = $blast_out_file;
    $parse_out_file =~ s/out$/parsed/;
    my $out_fh = IO::File->new("> $parse_out_file") ||
	die "Can not create file handle for $parse_out_file";
    # get a Taxon from a Bio::DB::Taxonomy object
    my $tax_dir = File::Temp::tempdir (CLEANUP => 1);
    my $dbh_sqlite = DBI->connect("dbi:SQLite:/gscmnt/sata835/info/medseq/virome/taxonomy_db");
    my $dbh = Bio::DB::Taxonomy->new(-source    => 'flatfile',
				     -directory => "$tax_dir",
				     -nodesfile => '/gscmnt/sata835/info/medseq/virome/taxonomy/nodes.dmp',
				     -namesfile => '/gscmnt/sata835/info/medseq/virome/taxonomy/names.dmp',);
    my @keep_for_tblastx = (); # query should be kept for further analysis
    my $total_records = 0;
    my $report = new Bio::SearchIO(-format => 'blast', -file => $blast_out_file, -report_type => 'blastn');
    unless ($report) {
	$self->log_event("Failed to create Bio SearchIO to parse ".basename($blast_out_file));
	return;
    }
    while(my $result = $report->next_result) {# next query output
	$total_records++;
	my $keep_for_tblastx = 1;  my %assignment = ();   my $best_e = 100;  my $hit_count = 0;
	while(my $hit = $result->next_hit) {
	    my @temp_arr = split(/\|/, $hit->name); # gi|num|database|accessionNum|
	    next if $temp_arr[2] eq 'pdb'; # skip data from pdb database
	    $hit_count++;
	    if ($hit_count == 1) {
		$best_e = $hit->significance;
	    }
	    # check whether the hit should be kept
	    if ($best_e <= $E_cutoff) { # similar to known, need Phylotyped
		$keep_for_tblastx = 0;
		if ($hit->significance == $best_e) { # only get best hits
		    # from gi get taxonomy lineage
		    my $sth = $dbh_sqlite->prepare("SELECT * FROM gi_taxid where gi = $temp_arr[1]");
		    $sth->execute();
		    unless ($sth->execute()) {
			$self->log_event("Failed to get taxonomy for gi = $temp_arr[1] in ".basename($blast_out_file));
		    }
		    my $ref = $sth->fetchrow_hashref();
		    $sth->finish();
		    if ($ref->{'taxid'}) { # some gi don't have record in gi_taxid_nucl
			my $taxon_obj = $dbh->get_taxon(-taxonid => $ref->{'taxid'});
			if (!(defined $taxon_obj)) {
			    my $description = "undefined taxon ".$hit->description."\t".$hit->name."\t".$hit->significance;
			    $assignment{"other"} = $description;
			}
			else {
			    my $tree_function = Bio::Tree::Tree->new();
			    my @lineage = $tree_function->get_lineage_nodes($taxon_obj);
			    # each lineage node is a Bio::Tree::NodeI object
			    if (scalar @lineage) {				
				$self->PhyloType(\@lineage,$hit, $best_e, $dbh_sqlite, $dbh, \%assignment);
			    }
			}
		    }	
		    else { # for situations that gi does not have corresponding taxid
			my $desc = $hit->description."\t".$hit->name."\t".$hit->significance;
			$assignment{"other"} = $desc;
		    } 
		}
		else {
		    last;
		}
	    }
	}  #END OF WHILE HIT LOOP
	# consolidate assignment
	# If a query is assigned both Homo and Primates, it will be reported as Homo only
	# If a query is assigned a real taxon name and "other" for reason like"other sequences;
	# artificial sequences", or no taxon id in taxon database it will be reported only as 
	# the real taxon name
	delete $assignment{'other'} if
	    exists $assignment{'other'} && keys %assignment > 1;
	#READS TO KEEP FOR NT BLASTX ANALYSIS
	push @keep_for_tblastx, $result->query_name if
	    $keep_for_tblastx == 1;

	# print out assignment for this query
	foreach my $assign (keys %assignment) {
	    $out_fh->print($result->query_name,"\t",$result->query_length,"\t",$assign,"\t",$assignment{$assign},"\n");
	}
    }#END OF WHILE RESULT LOOP

    $out_fh->print("# Summary: ", scalar @keep_for_tblastx, " out of $total_records ", scalar @keep_for_tblastx/$total_records, " is saved for TBLASTX analysis.\n");
    $out_fh->close;

    $dbh_sqlite->disconnect();

    #CREATE A FASTA FILE OF ALL UNKNOWN SEQUENCES TO RUN NT BLASTX
    my $root_file_name = $blast_out_file;
    $root_file_name =~ s/\.blastn\.out$//;
    my $orig_fasta = $root_file_name.'.fa';
    my $filter_out_file = $root_file_name.'.BNfiltered.fa';
    my $hits_file = $root_file_name.'.BNhits.fa';

    my $in = Bio::SeqIO->new(-format => 'fasta', -file => $orig_fasta);
    my $hits = Bio::SeqIO->new(-format => 'fasta', -file => ">$hits_file");
    my $out = Bio::SeqIO->new(-format => 'fasta', -file => ">$filter_out_file");
    unless ($in and $out) {
	$self->log_event("Failed to create SeqIO to read/write sequence for ".basename($blast_out_file));
	return;
    }
    while (my $seq = $in->next_seq) {
	my $read_name = $seq->primary_id;
	if (grep (/^$read_name$/, @keep_for_tblastx)) {
	    $out->write_seq($seq);
	}
        else
        {
            $hits->write_seq($seq);
        }
    }

    return 1;
}
	
############################################
sub PhyloType {
    my ($self,$lineage_ref, $hit_ref, $best_e, $dbh_sqlite, $dbh_taxonomy, $assignment_ref) = @_;
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
	if ($name eq "Metazoa") {
	    # make assignment
	    for (my $i = 0; $i <= $#$lineage_ref; $i++) { 
		my $temp_node_id = $lineage_ref->[$i]->id;
		my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid=>$temp_node_id);
		my $temp_name = $temp_obj->scientific_name;
		if ($temp_name eq "Homo") {
		    $description .= "Homo\t".$hit_ref->name."\t".$hit_ref->significance;
		    $assignment_ref->{"Homo"} = $description;
		    $assigned = 1;
		    last;
		}
	    }
	    if (!$assigned) {
		for (my $i = 0; $i <= $#$lineage_ref; $i++) { 
		    my $temp_node_id = $lineage_ref->[$i]->id;
		    my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid=>$temp_node_id);
		    my $temp_name = $temp_obj->scientific_name;
		    
		    if ($temp_name eq "Mus") {
			$description .= "Mus\t".$hit_ref->name."\t".$hit_ref->significance;
			$assignment_ref->{"Mus"} = $description;
			$assigned = 1;
			last;
		    }
		}
	    }
	    if (!$assigned) {
		$description .= $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
		$assignment_ref->{"other"} = $description;
		$assigned = 1;
	    }
	}
    }
    
    # check to see if it is bacteria sequence
    if ((scalar @{$lineage_ref} >= 2)&&(!$assigned)) {
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
    if (!$assigned) {
	$node_id = $lineage_ref->[0]->id;
	$obj = $dbh_taxonomy->get_taxon(-taxonid=>$node_id);
	$name = $obj->scientific_name;
	if ($name eq "Viruses") {
	    for (my $i = 0; $i <= $#$lineage_ref; $i++) { 
		my $temp_node_id = $lineage_ref->[$i]->id;
		my $temp_obj = $dbh_taxonomy->get_taxon(-taxonid=>$temp_node_id);
		my $temp_name = $temp_obj->scientific_name;
		$description .= $temp_name.";";
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
    $description = "";
    if (!$assigned) {
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
    if ((scalar @{$lineage_ref} >= 4)&&(!$assigned)) {
	$node_id = $lineage_ref->[3]->id;
	$obj = $dbh_taxonomy->get_taxon(-taxonid=>$node_id);
	$name = $obj->scientific_name;
	if ($name eq "Fungi") {
	    $description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
	    $assignment_ref->{"Fungi"} = $description;
	    $assigned = 1;
	}
    }
    
    # if still not assigned, assigned to "other" category
    if (!$assigned) {
	$description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
	$assignment_ref->{"other"} = $description;
	$assigned = 1;
    }
    
    return $assigned;
}

1;

