package Genome::Model::Tools::Annotate::TranscriptRegions;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::TranscriptRegions {
    is => 'Command',                       
    has => [ 
	organism => {
	    type  =>  'String',
	    doc   =>  "provide the organism either mouse or human; default is human",
	    is_optional  => 1,
	    default => 'human',
	},
	version => {
	    type  =>  'String',
	    doc   =>  "provide the imported annotation version; default for human is 54_36p_v2 and for mouse is 54_37g_v2",
	    is_optional  => 1,
	},
	   output => {
	    type  =>  'String',
	    doc   =>  "provide a file name to write you transcript information to .txt will be appended to it. Default is to print to stdout.",
	    is_optional  => 1,
	},
	chromosome => {
	    type  =>  'String',
	    doc   =>  "chromosome ie {1,2,...,22,X,Y,M}",
	},
	start => {
	    type  =>  'Number',
	    doc   =>  "Start coordinate",
	},
	stop => {
	    type  =>  'Number',
	    doc   =>  "Stop coordinate;",
	},
    ], 
};


sub help_synopsis {
    return <<EOS

gmt annotate transcript-region -h

EOS
}

sub help_detail {
    return <<EOS 

will provide the transcript substructures in your range of specified coordinates from both the ensembl and then genbank annotation.

EOS
}


sub execute {

    my $self = shift;

    my $chromosome = $self->chromosome;
    my $start = $self->start;
    my $stop = $self->stop;
    my $organism = $self->organism;
    my $version = $self->version;



    unless ($version) {	
	if ($organism eq "mouse") {
	    $version = "54_37g_v2" ;
	} elsif ($organism eq "human") { 
	    $version = "54_36p_v2" ;
	} else { 
	    die "organism is restricted to mouse or human\n";
	}
    }

    my $output = $self->output;

    my ($ncbi_reference) = $version =~ /\_([\d]+)/;
    my $eianame = "NCBI-" . $organism . ".ensembl";
    my $gianame = "NCBI-" . $organism . ".genbank";
    my $build_source = "$organism build $ncbi_reference version $version";
    
    my $ensembl_build = Genome::Model::ImportedAnnotation->get(name => $eianame)->build_by_version($version);
    unless ($ensembl_build) { die qq(Couldn't get ensembl build info for $build_source\n);}

    my ($ensembl_data_directory) = $ensembl_build->determine_data_directory;

    my $genbank_build = Genome::Model::ImportedAnnotation->get(name => $gianame)->build_by_version($version);
    my ($genbank_data_directory) = $genbank_build->determine_data_directory;
    
    my (@et) = Genome::Transcript->get(data_directory => $ensembl_data_directory, reference_build_id => $ensembl_build->reference_sequence_id);
    my (@gt) = Genome::Transcript->get(data_directory => $genbank_data_directory, reference_build_id => $genbank_build->reference_sequence_id);
    
    my @join_array = (@et,@gt);
    

    if ($output) {
	open(OUT,">$output") || die "couldn't open the output file $output\n";
        #print OUT qq(chromosome,transcript_start,transcript_stop,source,organism,version,hugo_gene_name,gene_id,strand,transcript_name,transcript_status\n);
	print OUT qq(chromosome\ttranscript_start\ttranscript_stop\tsource\torganism\tversion\thugo_gene_name\tgene_id\tstrand\ttranscript_name\ttranscript_id\ttranscript_status\ttotal_substructures\n);
	print OUT qq(\tn\tstructure_type\tstructure_start\tstructure_stop\n);
    } else {
	print qq(chromosome\ttranscript_start\ttranscript_stop\tsource\torganism\tversion\thugo_gene_name\tgene_id\tstrand\ttranscript_name\ttranscript_id\ttranscript_status\ttotal_substructures\n);
	print qq(\tn\tstructure_type\tstructure_start\tstructure_stop\n);
    }


    my $transcript_number = 0;

    for my $t (@join_array) {
	
	my $chrom_name = $t->chrom_name;
	
	next unless $chromosome eq $chrom_name;
	
	my $g_id = $t->gene_id;
	my @gid = split(/[\s]+/,$g_id);
	my ($gene_id) = $gid[0];
	
	my $source = $t->source;
	my $transcript_status = $t->transcript_status;
	my $strand = $t->strand;
	my $transcript_name = $t->transcript_name;
	my $transcript_id = $t->transcript_id;
	my $transcript_start = $t->transcript_start;
	my $transcript_stop = $t->transcript_stop;
	
	next unless $start >= $transcript_start && $start <= $transcript_stop ||
	    $start <= $transcript_start && $start >= $transcript_stop ||
	    $stop >= $transcript_start && $stop <= $transcript_stop ||
	    $stop <= $transcript_start && $stop >= $transcript_stop;
	
	my $gene = $t->gene;
	my $hugo_gene_name = $gene->hugo_gene_name;
	unless ($hugo_gene_name) {$hugo_gene_name = "unknown";}
	my @substructures = $t->ordered_sub_structures;
	my $total_substructures = @substructures;
	my $t_n = 0; #substructure counter
	$transcript_number++;

	$self->{transcript}->{$transcript_number}->{transcript_name}=$transcript_name;
	$self->{transcript}->{$transcript_number}->{transcript_start}=$transcript_start;
	$self->{transcript}->{$transcript_number}->{transcript_stop}=$transcript_stop;
	$self->{transcript}->{$transcript_number}->{strand}=$strand;
	$self->{transcript}->{$transcript_number}->{transcript_id}=$transcript_id;
	$self->{transcript}->{$transcript_number}->{hugo_gene_name}=$hugo_gene_name;
	$self->{transcript}->{$transcript_number}->{total_substructures}=$total_substructures;
	$self->{transcript}->{$transcript_number}->{source}=$source;
	$self->{transcript}->{$transcript_number}->{transcript_status}=$transcript_status;

	my $n_ss = $total_substructures - 2; #subtracting out the flanking regions
	$self->{transcript}->{$transcript_number}->{total_substructures}=$n_ss;

	if ($output) {
	    print OUT qq($chromosome\t$transcript_start\t$transcript_stop\t$source\t$organism\t$version\t$hugo_gene_name\t$gene_id\t$strand\t$transcript_name\t$transcript_id\t$transcript_status\t$n_ss\n);
	} else {
	    #print qq($chromosome,$transcript_start,$transcript_stop,$source,$organism,$version,$hugo_gene_name,$gene_id,$strand,$transcript_name,$transcript_status\n);
	    print qq($chromosome\t$transcript_start\t$transcript_stop\t$source\t$organism\t$version\t$hugo_gene_name\t$gene_id\t$strand\t$transcript_name\t$transcript_id\t$transcript_status\t$n_ss\n);
	}

	my $n = 0;
	if (@substructures) {
	    
	    while ($t_n < $total_substructures) {
		my $t_region = $substructures[$t_n];
		$t_n++;

		my $structure_type = $t_region->{structure_type};
		unless ($structure_type eq "flank") {$n++;}

		my $tr_start = $t_region->{structure_start};
		my $tr_stop = $t_region->{structure_stop};
				
		next unless $start >= $tr_start && $stop <= $tr_start ||
		    $start <= $tr_start && $stop >= $tr_start ||
		    $start >= $tr_stop && $stop <= $tr_stop ||
		    $start <= $tr_stop && $stop >= $tr_stop;

		#unless ($structure_type eq "intron" || $structure_type eq "flank") {
		unless ($structure_type eq "flank") {
		    if ($output) {
			#print OUT qq($chromosome,$tr_start,$tr_stop,$source,$organism,$version,$hugo_gene_name,$gene_id,$strand,$transcript_name,$transcript_status,$structure_type\n);
			print OUT qq(\t$n\t$structure_type\t$tr_start\t$tr_stop\n);
		    } else {
			#print qq($chromosome,$tr_start,$tr_stop,$source,$organism,$version,$hugo_gene_name,$gene_id,$strand,$transcript_name,$transcript_status,$structure_type\n);
			print qq(\t$n\t$structure_type\t$tr_start\t$tr_stop\n);
		    }
		    my $structure="$structure_type:$tr_start:$tr_stop";
		    $self->{transcript}->{$transcript_number}->{structure}->{$n}=$structure;
		}
	    }
	}
    }
    return($self);
}

1;

=head1 TITLE

TranscriptRegions

=head1 DESCRIPTION

This script will produce transcript info for defined regions

=head1 Input Options:

chrmosome start stop organism version output

=head1 KNOWN BUGS

Please report bugs to <rmeyer@genome.wustl.edu>

=head1 AUTHOR

Rick Meyer <rmeyer@genome.wustl.edu>

=cut
