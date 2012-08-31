package Genome::Model::Tools::Snp::GetDbsnps;

use strict;
use warnings;
use Genome;
use FileHandle;

use Bio::Seq;

class Genome::Model::Tools::Snp::GetDbsnps {
    is => 'Command',                    
    has => [ # specify the command's properties (parameters) <--- 
	     out => {
		 type  =>  'String',
		 doc   =>  "provide a name for your output file or it will simply be printed as stdout",
		 is_optional  => 1,
	     },
	     chromosome => {
		 type  =>  'String',
		 doc   =>  "chromosome ie {1,2,...,22,X,Y,M}",
		 is_optional  => 1,
	     },
	     start => {
		 type  =>  'Number',
		 doc   =>  "build 36 start coordinate",
		 is_optional  => 1,
	     },
	     stop => {
		 type  =>  'Number',
		 doc   =>  "build 36 stop coordinate;",
		 is_optional  => 1,
	     },
	     ref => {
		 type  =>  'String',
		 doc   =>  "referance allele;",
		 is_optional  => 1,
	     },
	     var => {
		 type  =>  'String',
		 doc   =>  "variant allele;",
		 is_optional  => 1,
	     },
	     list => {
		 type  =>  'string',
		 doc   =>  "a list of positions with ref and variant alleles for match check",
		 is_optional  => 1,
	     },
	     gff => {
		 type  =>  'Boolean',
		 doc   =>  "optional gff file of all dbsnps from your input coverage range with genomic coordiante or if used with refseq-fasta option and you'll get a snp.gff based on the referance coordinate",
		 is_optional  => 1,
	     },
	     refseq_fasta => {
		 type => 'string',
		 doc  => "intended to be used with the gff option to produce a snp.gff based on refseq coordinates output for this gff based on refseq name",
		 is_optional  => 1,
	     },
	     organism => {
		 type  =>  'String',
		 doc   =>  "provide the organism either mouse or human; default is human",
		 is_optional  => 1,
		 default => 'human',
	     },
	]
};


sub help_brief {
    return <<EOS
  This tool was design to retrieve dbsnp\'s for an individual site/range or list of sites and will check for an allele match to your variant if provided or it will simply state which sites on your list/input coords coinside with a dbsnp. It will work for either human build 36 or mouse build 37.
EOS
}

sub help_synopsis {
    return <<EOS

gmt snp get-dbsnps --list

the list is intended to be set up the same as with annotate transcript-variants

if the ref and var columns are omitted, then this script will still return dbsnps found in your input coordinates

using the --out option will allow your to get or list back in a file rather than as standard output to your screen 

chromosome start stop ref and var options will not be used if the list option is

EOS
}

sub help_detail {
    return <<EOS 

gmt snp get-dbsnps --list file

your list should be a tab/space delimited file with five columns
chromosome   start   stop   ref_allele   variant_allele

or

running...
gmt snp get-dbsnps --chromosome 1 --start 202785447 --stop 202785447 --ref A --var T

 will produce
      1 202785447 202785447 A T rs4252743:snp:1:'A/T':dbsnp_match

running...
gmt snp get-dbsnps --chromosome 1 --start 202785447 --stop 202785447 --ref A --var G

 will produce
      1 202785447 202785447 A G rs4252743:snp:1:'A/T':no_match

running...
gmt snp get-dbsnps  --gff --refseq-fasta 1_202785447_202795447.c1.refseq.fasta
 will produce a gff file of all the variation sequence tag entries in the database from said range



when detirmining the validation status of a dbsnp this tool assumes that if the snp was ever entered in the database as being validated that it is still validated.


if a dbsnp coinsides with the coordinates and the reference and variant alleles are provided, this script will check to see if the dbsnp alleles match your variant 

use the gff option to get a gff like file of the dbsnps identified by your input coordinates 

EOS
}


sub execute {

    my $self = shift;
    my $list;

    my $refseq_fasta = $self->refseq_fasta;
    my $gff = $self->gff;
    my $ref_list;
    if ($refseq_fasta && -e $refseq_fasta && $gff) {
	$ref_list = &getDBSNPS_GFF($self);
    }

    my $out = $self->out;
    if ($out) {open(OUT,">$out") || die "\nCould not open $out file\n";}

    my $file = $self->list;
    my $chr = $self->chromosome;
    my $start = $self->start;
    my $stop = $self->stop;

    if ($file) {
	unless (-f $file) {system qq(gmt snp get-dbsnps --help);print qq(Your list was not found.\t\n\tPlease check that your list is in place and try again.\n\n\n);return 0;}
	open(LIST,"$file") || die "\nCould not open $file\n";
	while (<LIST>) {
	    chomp;
	    my $line = $_;
	    my ($chr,$start,$stop,$ref,$var) = split(/[\s]+/,$line);
	    unless ($ref) {$ref = "na";}
	    unless ($var) {$var = "na";}
	    
	    my $alleles = "$ref,$var";
	    $list->{$chr}->{$start}->{$stop}->{$alleles}=1;
	} close(LIST);

    } elsif ($chr && $start && $stop) {
	my $ref = $self->ref;
	my $var = $self->var;
	unless ($ref) {$ref = "na";}
	unless ($var) {$var = "na";}
	my $alleles = "$ref,$var";
	$list->{$chr}->{$start}->{$stop}->{$alleles}=1;
    }
    
    if ($list) {
	$list = &getDBSNPS($self,$list);
	foreach my $chr (sort keys %{$list}) {
	    foreach my $start (sort {$a<=>$b} keys %{$list->{$chr}}) {
		foreach my $stop (sort {$a<=>$b} keys %{$list->{$chr}->{$start}}) {
		    foreach my $alleles (sort keys %{$list->{$chr}->{$start}->{$stop}}) {
			my ($ref,$var) = split(/\,/,$alleles);

			my $dbsnp_info = $list->{$chr}->{$start}->{$stop}->{$alleles};
			if ($dbsnp_info eq "1") {$dbsnp_info = "no_dbsnp_hit";}
			
			if ($out) {
			    unless ($ref eq "na" && $var eq "na") {
				print OUT qq($chr $start $stop $ref $var $dbsnp_info\n); ###add source to dbsnp_info
			    }
			} else {
			    unless ($ref eq "na" && $var eq "na") {
				print qq($chr $start $stop $ref $var $dbsnp_info\n);
			    }
			}
		    }
		}
	    }
	}
	if ($out) {close(OUT);unless(-f $out) {`rm $out`;}}
        return ($list);

    } elsif ($refseq_fasta && $ref_list) {
	$list = $ref_list;
        return ($list);
    } else {
	system qq(gmt snp get-dbsnps --help);return 0;
    }
}

sub getDBSNPS_GFF {
    my ($self) = @_;
    my $refseq_fasta = $self->refseq_fasta;

    my ($chromosome,$start,$stop,$genomic_coord,$orientation,$ref_seq_length);

    my $refseq_file = new FileHandle ($refseq_fasta);
    while (<$refseq_file>) {
	chomp;
	my $line=$_;
	if ($line =~ /\>/) {
	    
	    if ($line=~ /\s+Chr\:(\S+)\,\s+/) {
		$chromosome=$1;
	    }
	    my ($fisrt_coord,$second_coord);
	    if ($line=~ /Ori\s+\(\+\)/) {
		($fisrt_coord,$second_coord)=$line =~ /Coords[\s]+(\d+)\S(\d+)/;
		$orientation="plus";
		$genomic_coord = $fisrt_coord - 1;
	    } elsif ($line=~ /Ori\s+\(\-\)/) {
		($fisrt_coord,$second_coord)=$line =~ /Coords[\s]+(\d+)\S(\d+)/;
		$orientation="minus";
		$genomic_coord = $second_coord + 1;
	    }
	    $ref_seq_length = $second_coord - $fisrt_coord + 1;
	    $start = $fisrt_coord;
	    $stop = $second_coord;
	}
    }

    $self->{refseq_start} = $start;
    $self->{refseq_stop} = $stop;
    $self->{refseq_chr} = $chromosome;
    $self->{refseq_orientation} = $orientation;
    $self->{refseq_genomic_coord} = $genomic_coord;
    $self->{ref_seq_length} = $ref_seq_length;
    my $gff_out = $refseq_fasta;
    $gff_out =~ s/\.c1.refseq.fasta$//;
    $gff_out =~ s/\.refseq.fasta$//;
    $gff_out =~ s/\fasta$//;
    $gff_out = "$gff_out.dbsnp.gff";
    $self->{gff_out}=$gff_out;

    &getDBSNPS($self);

    return $self;
}

sub ref_coord {
    
    #my ($pos,$orientation) = @_;
    my ($self,$pos)  = @_;
    my $orientation = $self->{refseq_orientation};
    my $genomic_coord = $self->{refseq_genomic_coord};

    my $ref_pos;
    if ($orientation eq "plus") {
	$ref_pos = $pos - $genomic_coord;
    } elsif ($orientation eq "minus") {
	$ref_pos = $genomic_coord - $pos;
    }
    return $ref_pos;
}

sub getDBSNPS {

    my ($self,$list) = @_;   
    my $dbsnp;

    my $gff = $self->gff;
    my $organism = $self->organism;
    my $refseq_fasta = $self->refseq_fasta;

    my $g;
    if ($organism eq "human") {
	$g = GSC::Sequence::Genome->get(sequence_item_name => 'NCBI-human-build36');
    } elsif ($organism eq "mouse") {
	$g = GSC::Sequence::Genome->get(sequence_item_name => 'NCBI-mouse-buildC57BL6J');
    } else {
	print qq(Organism choices are restricted to either the default human or mouse.\n);
	exit (1);
    }

    if ($list) {
	my $dbsnp_out;
	if ($gff) {
	    unless ($self->refseq_fasta) {
		my $out = $self->out;
		if ($out) {
		    $dbsnp_out = "$out.dbsnp.gff";
		} else {
		    my $file = $self->list;
		    my $chr = $self->chromosome;
		    my $start = $self->start;
		    my $stop = $self->stop;
		    if ($file) {
			$dbsnp_out = "$file.dbsnp.gff";
		    } elsif ($chr && $start && $stop) {
			$dbsnp_out = "$chr:$start\_$stop.dbsnp.gff";
		    } else {
			$dbsnp_out = "genomic.dbsnp.gff";
		    }
		}
		open(GFF,">$dbsnp_out") || die "\n couldn't open the output dbsnp.gff $dbsnp_out\n";
	    }
	}
	foreach my $chr (sort keys %{$list}) {
	    next unless $chr =~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y)$/;
	    my $c = $g->get_chromosome($chr);
	    foreach my $start (sort {$a<=>$b} keys %{$list->{$chr}}) {
		next unless $start =~ /^[\d]+$/;
		foreach my $stop (sort {$a<=>$b} keys %{$list->{$chr}->{$start}}) {
		    next unless $stop =~ /^[\d]+$/;
		    
		    #for my $pos ($start..$stop) {
		    my @t = $c->get_tags( begin_position => { operator => 'between', value => [$start,$stop] },);
		    #my @t = $c->get_tags( begin_position => {  operator => 'between', value => [$pos,$pos] },);

		    my $refid_source;
		    for my $t (@t) {
			next unless $t->sequence_item_type eq 'variation sequence tag';
			
			my $t_start = $t->begin_position;
			my $t_stop = $t->end_position;
			my $pos = $t_start;
			
			my $variation_type = $t->variation_type;
			my $ref_id = $t->ref_id;
			my $allele_description = $t->allele_description;
			$allele_description =~ s/\'//gi;

			my $validated = $t->is_validated;
			my $seq_length = $t->seq_length;
			my $stag_id = $t->stag_id;
			my $seq_id = $t->seq_id;
			
			my $unzipped_base_string = $t->unzipped_base_string;
			
			my $end = $pos + ($seq_length - 1);
			
			unless ($validated) { $validated = 0; }

			my @sources;
			if ($seq_id) {
			    my @seq_col = GSC::SequenceCollaborator->get(seq_id => $seq_id);
			    if (@seq_col) {
				for my $source (@seq_col) {
				    my $r = $source->role_detail;
				    my $c = $source->collaborator_name;
				    my $e = "$c\_$r";
				    push @sources , $e unless grep (/$e/, @sources);
				    $refid_source->{$ref_id}->{$e}=1;
				}
			    }
			}
			my $source = join "||" , @sources;
						
			my ($vt,$v,$ad,$s);
			if ($dbsnp->{$chr}->{$pos}->{$ref_id}) {
			    ($vt,$v,$ad,) = split(/\:/,$dbsnp->{$chr}->{$pos}->{$ref_id});
			}
			
			if ($validated == 1) {
			    $dbsnp->{$chr}->{$pos}->{$ref_id}="$variation_type\:$validated\:$allele_description";
			} elsif ($vt && $v && $ad ) {
			    if ($v == 1) {$dbsnp->{$chr}->{$pos}->{$ref_id}="$vt\:$v\:$ad";}
			} else {
			    $dbsnp->{$chr}->{$pos}->{$ref_id}="$variation_type\:$validated\:$allele_description";
			}
				
			if ($self->gff) {
			    unless ($self->refseq_fasta) {
				print GFF qq(Chromosome$chr\tDB\t$variation_type\t$t_start\t$t_stop\t.\t+\t.\t$ref_id \; Alleles \"$allele_description\" ; Validation_Status \"$validated\" ; seq_id \"$seq_id\" ; source \"$source\"\n);
				#print GFF qq(Chromosome$chr\tDB\t$variation_type\t$pos\t$end\t.\t+\t.\t$stag_id\t$ref_id \; Alleles $allele_description \; Validation $validated\n);
			    }
			}
		    }
		    
		    foreach my $alleles (sort keys %{$list->{$chr}->{$start}->{$stop}}) {
			my ($ref,$var) = split(/\,/,$alleles);
			foreach my $ref_id (sort keys %{$dbsnp->{$chr}->{$start}}) {
			    my @sources;
			    foreach my $source (sort keys %{$refid_source->{$ref_id}}) {
				push @sources , $source unless grep (/$source/, @sources);
			    }
			    my $combined_source = join "||" , @sources;

			    my ($variation_type,$validated,$allele_description) = split(/\:/,$dbsnp->{$chr}->{$start}->{$ref_id});
			    my $match = &check_match($ref,$var,$allele_description);
			    unless ($match) {$match="no_match";}
			    my $snpfo = $list->{$chr}->{$start}->{$stop}->{$alleles};
			    if ($snpfo eq "1") {
				$list->{$chr}->{$start}->{$stop}->{$alleles}="$ref_id\:$combined_source\:$variation_type\:$validated\:$allele_description\:$match";
			    } else {
				unless ($snpfo =~ /$combined_source\:$variation_type\:$validated\:$allele_description\:$match/) {
				    $list->{$chr}->{$start}->{$stop}->{$alleles}="$snpfo\:\:$ref_id\:$combined_source\:$variation_type\:$validated\:$allele_description\:$match";
				}
			    }
			}
		    }
		}
	    }
	}
	if ($self->gff) { unless ($self->refseq_fasta) { print qq(dbsnp.gff writen to => $dbsnp_out\n); close(GFF); }}
	
	return $list;
	
    } elsif ($gff && $self->refseq_fasta) {
	
	my $start = $self->{refseq_start};
	my $stop = $self->{refseq_stop};
	my $chromosome = $self->{refseq_chr};
	my $gff_out = $self->{gff_out};
	my $root = $gff_out;
	$root =~ s/\.snp\.gff//;
	my @root = split(/\//,$root);
	$root = pop(@root);

	open(GFF,">$gff_out");
	
	my $c = $g->get_chromosome($chromosome);
	my @t = $c->get_tags( begin_position => { operator => 'between', value => [$start,$stop] },);
	for my $t (@t) {
	    next unless $t->sequence_item_type eq 'variation sequence tag';
	    my $ref_start = &ref_coord($self,$t->begin_position);# - $genomic_coord;
	    my $ref_stop = &ref_coord($self,$t->end_position);# - $genomic_coord;
	    
	    my $variation_type = $t->variation_type;
	    my $ref_id = $t->ref_id;
	    my $allele_description = $t->allele_description;

	    $allele_description =~ s/\'//gi;
	    my $validation_status = $t->is_validated;
	    unless ($validation_status) {$validation_status = 0;}
	    my $seq_id = $t->seq_id;
	    my @sources;
	    if ($seq_id) {
		my @seq_col = GSC::SequenceCollaborator->get(seq_id => $seq_id);
		if (@seq_col) {
		    for my $source (@seq_col) {
			my $r = $source->role_detail;
			my $c = $source->collaborator_name;
			my $e = "$c\_$r";
			push @sources , $e unless grep (/$e/, @sources);
		    }
		}
	    }
	    my $source = join "::" , @sources;

	    
	    print GFF qq($root\tDB\t$variation_type\t$ref_start\t$ref_stop\t.\t+\t.\t$ref_id \; Alleles \"$allele_description\" ; Validation_Status \"$validation_status\" ; seq_id \"$seq_id\" ; source \"$source\"\n);
	}
	close (GFF);
	print qq(dbsnp.gff writen to => $gff_out\n);
    } 
}


sub check_match { ### order and orientation of dbsnp allele_description is ambiguous 

    my ($ref,$var,$allele_description) = @_;
    
    $allele_description =~ s/\'//gi;

    my @dbsnp_allele_array = split(/\//,$allele_description);
    my $array_n = @dbsnp_allele_array;


    my ($rm,$vm);
    for my $n (1..$array_n) {
	my $m = $n - 1;
	my $dbsnp_allele = $dbsnp_allele_array[$m];
	if ($ref eq $dbsnp_allele) {
	    $rm = 1;
	} elsif ($var eq $dbsnp_allele) {
	    $vm = 1;
	}
    }

    unless ($rm && $vm) {
	undef($rm);
	undef($vm);
	for my $n (1..$array_n) {
	    my $m = $n - 1;
	    my $dbsnp_allele = $dbsnp_allele_array[$m];
	    my $rev_dbsnp_allele = &reverse_complement_allele ($dbsnp_allele); 
	    if ($ref eq $rev_dbsnp_allele) {
		$rm = 1;
	    } elsif ($var eq $rev_dbsnp_allele) {
		$vm = 1;
	    }
	}
    }
    
    my $dbsnp_match;
    if ($rm && $vm) {
	$dbsnp_match = "dbsnp_match";
    } else {    
	$dbsnp_match = "no_match";
    }
    return ($dbsnp_match);
}	    


sub reverse_complement_allele {
    my ($allele_in) = @_;

    if ($allele_in =~ /[\-\+X]/) {
	return $allele_in;
    } else {
	if ($allele_in =~ /^[ACGT]+$/) {
	    my $seq_1 = new Bio::Seq(-seq => $allele_in);
	    my $revseq_1 = $seq_1->revcom();
	    my $rev1 = $revseq_1->seq;
	    return $rev1;
	} else {
	    return $allele_in;
	}
    }
}


1;
