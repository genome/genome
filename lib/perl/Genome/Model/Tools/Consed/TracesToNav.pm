package Genome::Model::Tools::Consed::TracesToNav;

use strict;
use warnings;
use Genome;
use IPC::Run;

class Genome::Model::Tools::Consed::TracesToNav {
    is => 'Command',                    
    has => [ # specify the command's properties (parameters) <--- 
	     ace              => {
		 type         => 'String',
		 doc          => "the full path to and including the ace file you want to navigate",
	     },
	     list             => {
		 type         => 'String',
		 doc          => "provide a comma delimited list of sample sites to review",
	     },
	     name_nav          => {
		 type         => 'String',
		 doc          => "allows you to provide a descriptive name for the manual review spreadsheet defaults is sites_for_review.date",
		 is_optional  => 1,
	     },
	     convert_coords   => {
		 type         => 'String',
		 doc          => "if the coordinates on your list are genomic, you may use this option by providing the refseq fasta file to get the refseq coordinate to navigate", 
		 is_optional  => 1,
	     },
	     unpaired         => {
		 type         => 'Boolean',
		 doc          => "if your list does not included paired sample info use use this option", 
		 is_optional  => 1,
	     },
	     input_type       => {
		 type         => 'String',
		 doc          => "state simple paired or expanded",

	     },

	     ],
    
    
};


sub help_brief {
    "This tool will make a consed ace navigator from a list of samples and coordinates"     
}

sub help_synopsis {
    return <<EOS
gmt consed traces-to-nav --ace --list

running...

gmt consed traces-to-nav --ace $ENV{GENOME_TEST_INPUTS}Genome-Model-Tools-Consed-TracesToConsed/10_126008345_126010576/edit_dir/10_126008345_126010576.ace.1 --convert-coords $ENV{GENOME_TEST_INPUTS}Genome-Model-Tools-Consed-TracesToConsed/10_126008345_126010576/edit_dir/10_126008345_126010576.c1.refseq.fasta --input-type simple --name-nav test.traces.to.nav --list $ENV{GENOME_TEST_INPUTS}Genome-Model-Tools-Consed-TracesToNav/Nav.list

will produce a navigator ==> test.traces.to.nav.date.nav
and a spreadsheet ==> test.traces.to.nav.date.csv

EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 


Formating Your List
     Your list should be a comma delimited plain text file in one of the following two formats

          if you want to navigate paired samples
	     sample,pos,paired_sample,comment

          if you do not want to navigate paired samples
	     sample,pos,comment

The spreadsheet test.traces.to.nav.date.csv will be a tab delimited file in this format
refpos  sample  reads   note    manual_genotype comments

where manual_genotype and comments are left blank

The comment is what ever was supplied in the list and should not contain any commas


EOS
}


sub execute {                               

    my $self = shift;
    my $ace = $self->ace;
    my $list = $self->list;
    my $name_nav = $self->name_nav;

    my @subdir = split(/\//,$ace);
    my $ace_name = pop(@subdir);
    my $edit_dir = join('/',@subdir);
    chdir $edit_dir;
    
    my $date_tag = &get_date_tag();
    
    
    my ($out);
    if ($name_nav) {
	$out = "$edit_dir/$name_nav.$date_tag";
    } else {
	$out = "$edit_dir/sites_for_review.$date_tag";
    }

    my $refseq = $self->convert_coords;
    my ($refseq_info);
    if ($refseq) {
	($refseq_info) = &parse_ref($self, $refseq);
    } else {
	$refseq_info = 0;
    }

    my ($sites_to_nav) = &get_sites_to_nav($self,$refseq_info); #read in list
    unless ($sites_to_nav) {
	$self->error_message( "\n\ndidn't find sites to nav\n\n");
	return;
    }

    my ($reads_to_nav,$main_contig) = &get_reads_to_nav($self, $ace_name,$sites_to_nav); #find the reads for samples from the list
    unless ($reads_to_nav) {
	$self->error_message( "\n\ndidn't find reads to nav\n\n");
    }


    my ($done) = &make_read_nav($sites_to_nav,$reads_to_nav,$main_contig,$out,$self); #write the navigator and if desired the spreadsheet
    if ($done) {
	return 1;    
    } else {
	$self->error_message( "\n\ndidn't write to $out.nav or $out.csv\n\n");
	return;
    }
}


sub get_sites_to_nav {
    
    my ($self,$refseq_info) = @_;

    my $list = $self->list;
#(file_check);

    my ($sites_to_nav,$samples_pos_info,$tumor_samples,$paired_samples);
    open(LIST,"$list") ||  $self->error_message( "\n\ncouldn't open the list $list\n\n") && return;
    while (<LIST>) {
	chomp;
	
	my ($line) = $_;
	#my ($Gene,$pos,$note,$sample,$sample_gt,$pair,$pair_gt,$somatic_status)= split(/\,/,$line);
	
	my ($sample,$pos,$pair,$note);
	if ($self->input_type eq "simple") {
	    ($sample,$pos,$note) = split(/\,/,$line);
	    $pair = $sample;
	} elsif ($self->input_type eq "paired") {
	    ($sample,$pos,$pair,$note) = split(/\,/,$line);
	} elsif ($self->input_type eq "expanded") {
	    $self->error_message( "\n\ninput_type expanded not yet supported\n\n");
	    return;
	}
	if ($refseq_info) {
	    $pos = &get_refpos($refseq_info,$pos);
	}
	
	unless ($note) {$note="no comment";}
	
	$sites_to_nav->{$pos}->{$sample}->{comment}=$note;
	$sites_to_nav->{$pos}->{$sample}->{primary}=$sample;
	$sites_to_nav->{$pos}->{$sample}->{pair}=$pair;
	$sites_to_nav->{$pos}->{$pair}->{pair}=$sample;
	
    }
    close (LIST);
    return ($sites_to_nav);#,$samples_pos_info,$tumor_samples,$paired_samples);
}





sub get_reads_to_nav {
    

    my ($self, $ace,$sites_to_nav) = @_;

    my $ao = Genome::Model::Tools::Consed::AceReader->create(file => $ace);
    unless ($ao) { $self->error_message( "\n\ndidn't get the ace object\n\n") && return;}
    my ($reads_to_nav,$main_contig);
    while ( my $contig = $ao->next_contig ) {
        if (grep { /\.c1$/ } keys %{ $contig->{reads} }) {
            $main_contig = $contig->{name};
        }
    }

    my %base_pad_count;
    my @hhh;

    my $p;
    my $q;
    my @con_seq;
    
    open (ACE_file, "$ace")  ||  $self->error_message( "\n\ncouldn't open the ace file\n\n") && return;
    my @seq_line = ();
    my @file = <ACE_file>;
    my $file_n = @file;
    close (ACE_file);
    $p = 0;
    
    
    while ($file_n >= $p){
	$q = 1;
	if ($file[$p]) {
	    if ($file[$p] =~ /CO $main_contig\s/) { 
		until ($file[$p + $q] =~ /BQ/){
		    chomp $file[$p + $q];
		    #print ("$file[$p + $q]");
		    @seq_line = split(//, $file[$p + $q]);
		    chomp @seq_line;
		    push @con_seq, @seq_line;
		    
		    $q++;
		}
	    }
	}
	$p++;
    }
    my $pad_count = 0;
    my $base_number = 0;
    foreach my $base (@con_seq){
	#$base_number++;
	if ($base =~ /\*/) { 
	    $pad_count++;
	} else {
	    $base_number++;
	    $base_pad_count{$base_number} = $pad_count;
	}
    }
    
    @hhh = %base_pad_count;
    
    
#	&Count_pads; #count the pads in the consensus
    %base_pad_count = @hhh; # 
    
    
    
    foreach my $pos (sort {$a<=>$b} keys %{$sites_to_nav}) {
	
            $ao->_fh->seek(0, 0);
    while ( my $contig = $ao->next_contig ) {
	    
        my $name = $contig->{name};
	    if ($name eq $main_contig) {
		my $info = $contig->{reads};
		
		foreach my $read (keys %{ $contig->{reads} }) {
		    unless ($read =~ /\.c1$/) {
			
			my ($id);
			
			foreach my $sample (sort keys %{$sites_to_nav->{$pos}}) {
			    
			    if ($read =~ /$sample/) {
				#my $sample = substr($read,$pretty_source_1,$pretty_source_2);
				$id = $sample;
			    } 
			    next if $id;
			}
			
			if ($id) {
			    if ($sites_to_nav->{$pos}->{$id}) {
				
				my $q1 = $info->{$read}->{qual_clip_start};
				my $q2 = $info->{$read}->{qual_clip_end};
				my $position = $info->{$read}->{position};
				
				if ($position) {
				    
				    my $length = length $info->{$read}->{sequence};
				    my $align_clip_start = $info->{$read}->{align_clip_start};
				    my $align_clip_end = $info->{$read}->{align_clip_end};
				    my $position2 = $position + $length - 1;
				    my $p1 = ($position  - $base_pad_count{$position});
				    my $p2 = $p1 + $length;
				    my $padded_pos = $pos + $base_pad_count{$pos};

				    my $pair = $sites_to_nav->{$pos}->{$id}->{pair};
				    my $comment = $sites_to_nav->{$pos}->{$id}->{comment};
				    unless ($comment) {$comment = $sites_to_nav->{$pos}->{$pair}->{comment};}

				    my $diff = "$pos $pos";
				    my $ref_q1 = $position + $align_clip_start + 1;
				    my $bases_clipped_off_end = $length - $align_clip_end;
				    my $ref_q2 = $position2 - $bases_clipped_off_end;
				    
				    if (($padded_pos >= $ref_q1)  && ($padded_pos  <= $ref_q2)) {
					
					print qq($read $p1 $position - $base_pad_count{$position} $length $diff pos\n);
					
					#$nav_select->{$pos}->{$id}->{$read}->{nav} = "BEGIN_REGION\nTYPE: READ\nCONTIG: $name\nREAD: $read\nUNPADDED_CONS_POS: $diff\nCOMMENT: $comment\nEND_REGION\n\n";
					$reads_to_nav->{$pos}->{$id}->{$read}->{nav} = "BEGIN_REGION\nTYPE: READ\nCONTIG: $name\nREAD: $read\nUNPADDED_CONS_POS: $diff\nCOMMENT: $comment\nEND_REGION\n\n";


					my $read_group = $sites_to_nav->{$pos}->{$id}->{nav};
					if ($read_group) {
					    unless ($read_group =~ /$read/) {
						#my $read1 = $sites_to_nav->{$pos}->{$id}->{nav};
						$sites_to_nav->{$pos}->{$id}->{nav} = "$read_group:$read";
					    }
					} else {
					    
					    $sites_to_nav->{$pos}->{$id}->{nav} = $read;
					}
				    } 
				}
			    }
			}
		    }
		}
	    }
	}
    }

    foreach my $pos (sort keys %{$sites_to_nav}) {
	my $navpos;
	foreach my $id (sort keys %{$sites_to_nav->{$pos}}) {
	    my $nav = $sites_to_nav->{$pos}->{$id}->{nav};
	    if ($nav) {
		$navpos = 1;
	    } else {
		$sites_to_nav->{$pos}->{$id}->{nav} = "notinnavigator";
	    }
	} 
	unless ($navpos) {
	    $sites_to_nav->{$pos}->{navconsensus} = "notinnavigator";
	}
    }
    return ($reads_to_nav,$main_contig);
}
####################################################################


sub make_read_nav {
    
    my ($sites_to_nav,$reads_to_nav,$main_contig,$out,$self) = @_;
    
#	$sites_to_nav->{$pos}->{$sample}->{comment}=$note;
#	$sites_to_nav->{$pos}->{$sample}->{pair}=$pair;
    
    open (NAV,">$out.nav") || $self->error_message( "\n\ncouldn't open $out.nav to for writting\n\n") && return;
    open (CSV,">$out.csv") || $self->error_message( "\n\ncouldn't open $out.csv to for writting\n\n") && return;
    if ($self->unpaired) {
	#($sample,$pos,$note) = split(/\,/,$line);
	print CSV qq(refpos\tsample\treads\tnote\tmanual_genotype\tcomments\n);
    } else {
	#($sample,$pos,$pair,$note) = split(/\,/,$line);
	print CSV qq(refpos\tsample\treads\tnote\tmanual_genotype\tsomatic_status\tcomments\n);
    }
    
    print NAV qq(TITLE:\n\n);
    my ($naved,$printed);
    foreach my $pos (sort {$a<=>$b} keys %{$sites_to_nav}) {
	my $navconsensus = $sites_to_nav->{$pos}->{navconsensus};
	
	if ($navconsensus) {

	    my $comment = "no reads in navigator at this position.";
	    my $nav_line =  "BEGIN_REGION\nTYPE: CONSENSUS\nCONTIG: $main_contig\nUNPADDED_CONS_POS: $pos $pos\nCOMMENT: $comment\nEND_REGION\n\n";
	    print NAV qq($nav_line);
	    
	} else {
	    
	    foreach my $sample (sort keys %{$sites_to_nav->{$pos}}) {
		
		my $pair = $sites_to_nav->{$pos}->{$sample}->{pair};
		my $comment = $sites_to_nav->{$pos}->{$sample}->{comment};
		
		my $primary_sample = $sites_to_nav->{$pos}->{$sample}->{primary};
		
		if ($primary_sample) {
		    my @samples;
		    
		    if ($primary_sample eq $pair) {
			@samples=($primary_sample);
		    } else {
			@samples=($primary_sample,$pair);
		    }
		    
		    for my $id (@samples) { 
			
			if ($sites_to_nav->{$pos}->{$id}->{nav} eq "notinnavigator") {

			    unless ($printed->{$pos}->{$id}) {
				$printed->{$pos}->{$id}=1;
				print CSV qq($pos\t$id\tnot in navigator\t$comment\t\n);
			    }

			} else {
			    foreach my $read (sort keys %{$reads_to_nav->{$pos}->{$id}}) {
				my $nav_line = $reads_to_nav->{$pos}->{$id}->{$read}->{nav};
				if ($nav_line) {
				    unless ($naved->{$pos}->{$id}->{$read}) {
					print NAV qq($nav_line);
					$naved->{$pos}->{$id}->{$read}=1;
				    }
				}
			    }
			    unless ($printed->{$pos}->{$id}) {
				my $reads = $sites_to_nav->{$pos}->{$id}->{nav};
				$printed->{$pos}->{$id}=1;
				print CSV qq($pos\t$id\t$reads\t$comment\t\n);
			    }
			}
		    }			
		}
	    }
	}
    }	

    close (NAV);
    close (CSV);
    if (-f "$out.nav" && -f "$out.csv") {
	return 1;
    } else {
	return;
    }

}


####################################################################
####################################################################
sub make_cons_nav {
    
    my ($self, $out,$main_contig,$sites_to_nav) = @_;
    open (NAV3,">$out.consensus.nav") || $self->error_message( "\n\ncouldn't open $out.consensus.nav to for writting\n\n") && return;

    print NAV3 qq(TITLE:\n\n);

    my $cons;    
    foreach my $pos (sort {$a<=>$b} keys %{$sites_to_nav}) {
	foreach my $sample (sort keys %{$sites_to_nav->{$pos}}) {
	    unless( $cons->{$pos} ) {
		my $comment = $sites_to_nav->{$pos}->{$sample}->{comment};
		if ($comment) {
		    $cons->{$pos} = 1;
		    my $comment = $sites_to_nav->{$pos}->{$sample}->{comment};
		    my $nav_line =  "BEGIN_REGION\nTYPE: CONSENSUS\nCONTIG: $main_contig\nUNPADDED_CONS_POS: $pos $pos\nCOMMENT: $comment\nEND_REGION\n\n";
		    print NAV3 qq($nav_line);
		}
	    }
	}
    }
    close(NAV3);
    if (-f "$out.consensus.nav") {
	return 1;
    } else {
	return;
    }
}


sub get_newest {
    my ($file) = @_;
    $file  = `ls $file`;
    chomp($file);
    my @file_array = ();
    @file_array = split(/\n/, $file);
    
    if($file_array[1])
    {
	@file_array = sort byDateDesc @file_array;
    }	
    
    $file = $file_array[0];
    #my @f_array = split(/\//,$file);
    #$file = pop(@f_array);
    
    return $file;
}


sub byDateDesc
{
    my @fileStats = stat $a;
    my $mtime_a = $fileStats[9];
    
    @fileStats = stat $b;
    my $mtime_b = $fileStats[9];	
    
    $mtime_b <=> $mtime_a;
}



sub get_date_tag {
    
    my $time=`date`;
    #my ($handle) = getpwuid($<);
    my $date = `date +%D`;
    (my $new_date) = $date =~ /(\d\d)\/(\d\d)\/(\d\d)/ ;
    my $date_tag = "$3$1$2" ;
    return $date_tag;
}

sub parse_ref {
    
    
    my ($self, $refseq) = @_;
    my $orientation;
    my $genomic_coord;
    my $refseq_info;
    open(REF,"$refseq")  || $self->error_message( "\n\ncouldn't open $refseq\n\n") && return;
    while (<REF>) {
	chomp; 
	my $line = $_;
	
	if ($line =~ /\>/) {
	    #my ($roi) = $line =~ /^\S+\.(\S+)\.c1\.refseq\.fasta\:\>/;
	    #my ($gene) = $line =~ /GeneName:(\S+)\,/;
	    #my ($chromosome) = $line =~ /Chr\:([\S]+)\,/;
	    
	    if ($line=~ /Ori\s+\(\+\)/) {
		my ($fisrt_coord)=$line =~ /Coords\s+(\d+)\S\d+/;
		$orientation="plus";
		$genomic_coord = $fisrt_coord - 1;
	    } elsif ($line=~ /Ori\s+\(\-\)/) {
		my ($second_coord)=$line =~ /Coords\s+\d+\S(\d+)/;
		$orientation="minus";
		$genomic_coord = $second_coord + 1;
	    }
	}
    }
    close(REF);
    $refseq_info->{orientation}=$orientation;
    $refseq_info->{genomic_coord}=$genomic_coord;
    return ($refseq_info);
    
}

sub get_refpos {

    my ($refseq_info,$genpos) = @_;
    my $orientation = $refseq_info->{orientation};
    my $genomic_coord = $refseq_info->{genomic_coord};
    
    my $ref_pos;
    
    if ($orientation eq "plus") {
	$ref_pos = $genpos - $genomic_coord;
	
    } elsif ($orientation eq "minus") {
	$ref_pos = $genomic_coord - $genpos;
    }
    return $ref_pos;
}

1;

