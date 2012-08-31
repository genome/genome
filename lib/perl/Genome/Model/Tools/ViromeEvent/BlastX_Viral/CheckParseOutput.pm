package Genome::Model::Tools::ViromeEvent::BlastX_Viral::CheckParseOutput;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use DBI();
use File::Temp;
use File::Basename;

class Genome::Model::Tools::ViromeEvent::BlastX_Viral::CheckParseOutput{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    return "gzhao's Blast X Viral check parse output";
}

sub help_detail {
    return <<"EOS"
This script will check all .tblastx_ViralGenome.parsed file in the 
BNfiltered_TBLASTX_ViralGenome subdirectory of the given directory to make sure 
parsing tblastx output file is finished for each file.
EOS
}

sub execute {
    my $self = shift;

    my $dir = $self->dir;
    my $sample_name = basename($dir);
    
    my $blast_dir = $dir.'/'.$sample_name.'.TBXNTFiltered_TBLASTX_ViralGenome';
    unless (-d $blast_dir) {
	$self->log_event("Failed to find Viral blastX dir for $sample_name");
	return;
    }

    # check files for input fasta and blast filtered fasta files
    my @fa_files = glob("$blast_dir/*fa");
    if ( not @fa_files ) {
        $self->log_event("No files found in viral blastx dir .. probably all reads were filtered out in earlier blast stage");
        return 1;
    } else {
        my $filtered_file = $dir.'/'.$sample_name.'TBXNTFiltered.fa';
        if ( not -s $filtered_file > 0 ) {
            $self->log_event('Viral blast filtered fasta file is empty .. probably all reads were filtered out');
        }
    }

    foreach my $fa_file (@fa_files) {
	my $root_name = $fa_file;
	$root_name =~ s/\.fa$//;
	
	my $blast_out_file = $root_name.'.tblastx_ViralGenome.out';
	unless (-s $blast_out_file) {
	    $self->log_event("Can not find Viral blastX output file for ".basename($fa_file));
	    return;
	}
	my $blast_parse_file = $root_name.'.tblastx_ViralGenome.parsed';

	if (-s $blast_parse_file) { #BLAST PARSE FILE EXISTS
	    $self->log_event("Checking existing Viral blastX parse file: ".basename($blast_parse_file));
	    if ($self->check_completed_parse($blast_out_file)){ #CHECK TO SEE IF PARSE RAN CORRECTLY
		$self->log_event("Viral blastX parse already complete for ".basename($blast_out_file));
	    }
	    else { #RE-RUN PARSING
		$self->log_event("Re-running Viral blastX parse for ".basename($blast_out_file));
		if ($self->run_parse($blast_out_file)) {
		    $self->log_event("Parsing re-ran successfully for ".basename($blast_out_file));
		}
		else {
		    $self->log_event("Parsing failed for ".basename($blast_out_file));
		    return;
		}
	    }
	}
	else { #BLAST PARSE FILE DOES NOT EXIST
	    $self->log_event("Running Viral blastX parse for ".basename($blast_out_file));
	    if ($self->run_parse($blast_out_file)) { #RUN PARSING
		$self->log_event("Parsing ran successfully for ".basename($blast_out_file));
	    }
	    else {
		$self->log_event("Parsing failed for ".basename($blast_out_file));
		return;
	    }
	}
    }

    $self->log_event("Completed checking Viral blastX parse for $sample_name");

    return 1;
}

sub check_completed_parse {
    my ($self, $parse_file) = @_;
    my $fh = IO::File->new("< $parse_file") ||
	die "Can not create file handle to read $parse_file";
    my $line_count = 0;  my $saved_seq = 0;  my $total_seq = 0;
    my $undef_taxon = 0; my $parse_completed = 0;
    while (my $line = $fh->getline) {
	$line_count++;
	$undef_taxon++ if $line =~ /undefined\s+taxon/;
	if ($line =~ /#\s+Summary:\s+(\d+)\s+out\s+of\s+(\d+)/) {
	    $parse_completed = 1;
	    $saved_seq = $1;
	    $total_seq = $2;
	    return 1; #SKIP ADDITIONAL CHECKS AND LEAVE HERE
	}
    }
    $fh->close;

    if ($parse_completed == 0) { #PARSING NEVER COMPLETED - REDO
	return;
    }
    my $phylotyped_count = $total_seq - $saved_seq;
    #???  SKIPPING CHECKS BELOW ???
    # deal with situation where all records showed as undefined taxon
    my $num_phylogypted = $total_seq - $saved_seq; 
    if ( ($num_phylogypted ne 0) && ($num_phylogypted <= $undef_taxon)) { 
	return 0;
    }
    # taxonomy record has to be equal or greater than the number of sequences get 
    # successful phylotyped because some sequence could be assigned multiple taxonomy
    # categories.
    # at this step, every sequence's assignment information should be there even with
    # unassigned. So record number has to be >+ total sequence
    if (($line_count -2 ) < $total_seq ) {
	return 0;
    }
    
    return 1;
}

sub run_parse {
    my ($self, $blast_out_file) = @_;

    my $E_cutoff = 1e-5;

    my $parse_out_file = $blast_out_file;
    $parse_out_file =~ s/out$/parsed/;

    my $out_fh = IO::File->new("> $parse_out_file") ||
	die "Can not create file handle for $parse_out_file";

    my @unassigned = (); # query should be kept for further analysis
    my $total_records = 0;

    # get a Taxon from a Bio::DB::Taxonomy object
    my $dbh_sqlite = DBI->connect("dbi:SQLite:/gscmnt/sata835/info/medseq/virome/taxonomy_db");
    my $tax_dir = File::Temp::tempdir (CLEANUP => 1);
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
        my $best_e = 100;
        my $hit_count = 0;
        my $determined = 0;
        while(my $hit = $result->next_hit) {
	    my @temp_arr = split(/\|/, $hit->name); # gi|num|database|accessionNum|
	    next if $temp_arr[2] eq 'pdb'; # skip data from pdb database
	    $haveHit = 1;
    	    $hit_count++;
	    if ($hit_count == 1) {
	        $best_e = $hit->significance;
	    }

	    if ($hit->significance <= $E_cutoff){ # similar to known, need Phylotyped
	        my $have_significant_hit = 1;
	        if ($hit->significance == $best_e) {
		    # from gi get taxonomy lineage
		    my $sth = $dbh_sqlite->prepare("SELECT * FROM gi_taxid where gi = $temp_arr[1]");
		    $sth->execute();
		    my $ref = $sth->fetchrow_hashref();
		    $sth->finish();
		    if ($ref->{'taxid'}) { # some gi don't have record in gi_taxid_nucl, this is for situation that has
		        my $taxon_obj = $dbh->get_taxon(-taxonid => $ref->{'taxid'});
		        if (!(defined $taxon_obj)) {
			    my $description .= "undefined taxon\t".$hit->name."\t".$hit->significance;
			    $assignment{"Viruses"} = $description;
		        }
		        else {
			    my $tree_function = Bio::Tree::Tree->new();
			    my @lineage = $tree_function->get_lineage_nodes($taxon_obj);
			    # each lineage node is a Bio::Tree::NodeI object
			    if (scalar @lineage) {				
			        $determined = 1;
			        $self->PhyloType(\@lineage,$hit, $best_e, $dbh_sqlite, $dbh, \%assignment);
			    }
		        }
		    }
		    else { # for situations that gi does not have corresponding taxid
		        my $desc = $hit->description."\t".$hit->name."\t".$hit->significance;
		        $determined = 1;
		        $assignment{"Viruses"} = $desc;
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

	# print out assignment
        foreach my $assign (keys %assignment) {
	    if ($assign eq "unassigned") {
	        $out_fh->print($result->query_name, "\t", $result->query_length, "\t", $assign, "\t", $assignment{$assign}, "\n");
	        push @unassigned, $result->query_name;
	    }
	    else {
	        $out_fh->print($result->query_name, "\t", $result->query_length, "\t", $assign, "\t", $assignment{$assign}, "\n");
	    }
        }
    } # end of report parsing
    $out_fh->print("# Summary: ", scalar @unassigned," out of $total_records ", (scalar @unassigned)*100/$total_records, "% is unassigned.\n");
    $out_fh->close;

    $dbh_sqlite->disconnect();

    return 1;
}
		
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
		    $description = "Homo\t".$hit_ref->name."\t".$hit_ref->significance;
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
			$description = "Mus\t".$hit_ref->name."\t".$hit_ref->significance;
			$assignment_ref->{"Mus"} = $description;
			$assigned = 1;
			last;
		    }
		}
	    }
	    if (!$assigned) {
		$description = $Lineage."\t".$hit_ref->name."\t".$hit_ref->significance;
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
