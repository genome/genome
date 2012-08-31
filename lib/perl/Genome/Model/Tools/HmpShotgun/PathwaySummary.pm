package Genome::Model::Tools::HmpShotgun::PathwaySummary;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

use Data::Dumper;

#gmt hmp-shotgun pathway-summary --ko-pathway-file=Genome/Model/Tools/HmpShotgun/sahar/pathway_ko_map.txt --read-ko-file=Genome/Model/Tools/HmpShotgun/sahar/ko_orig.out --pathway-tree-file=Genome/Model/Tools/HmpShotgun/sahar/pathways.txt --output-file=summary.out

class Genome::Model::Tools::HmpShotgun::PathwaySummary {
    is  => ['Command'],
    has => [
        ko_pathway_file => {
            is  => 'String',
            is_input => '1',
            doc => 'File that maps the kos to pathways.',
        },
        pathway_tree_file => {
            is => 'String',
            is_input => '1',
            doc => 'Pathway tree.',
        },
        read_ko_file => {
            is  => 'String',
            is_input => '1',
            doc => 'The read ko file.',
        },
        output_file => {
            is => 'String',
            is_output =>1,
            is_optional => 1,
            default=>0,
        },
    ],
};


sub help_brief {
    'Handle the unaligned reads in a novel fashion.';
}

sub help_detail {
    return <<EOS
    Handle the unaligned reads in a novel fashion.
EOS
}


sub execute {
    my $self = shift;
    $self->dump_status_messages(1);
    $self->dump_error_messages(1);
    $self->dump_warning_messages(1);

    
my $tree;
my $path=new FileHandle($self->pathway_tree_file);
my $main_pathway;
my $sub_pathway;
while(<$path>){
    chomp;
    my $line=$_;
    next if ($line =~ /^\s+$/);
    if ($line =~ /^\d+\.\s+/){
	$main_pathway=join(" ",$line);
    }
    elsif ($line =~ /^\d+\.\d+\s+/){
	$sub_pathway=$line;
    }else{
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;
	push(@{$tree->{$main_pathway}->{$sub_pathway}},$line);
    } 
}


#Read KO map
my $pathway_ko_data;
my $ko_pathway_data;
my $map=new FileHandle($self->ko_pathway_file);
while (<$map>) {
    chomp;
    my @line=split(/\t/,$_);
    my $konum=shift(@line);
    foreach my $k (@line){
	$k=~s/^\s+//;
	$ko_pathway_data->{$konum}->{$k}=1;
	$pathway_ko_data->{$k}->{$konum}=1;
    }
}
  
#Read parsed blast
my $ko_group_data;
my $pb=new FileHandle($self->read_ko_file);
my $header;
my $ko;

while (<$pb>) {
	chomp;
	my $line = $_;
	my @line = split(/\t/,$line);
	my $header = shift @line;
	foreach my $ko_item (@line) {
		$ko_group_data->{$ko_item}->{$header}=1;
		print "\n***** associating: $ko_item with $header\n";
	}
	
}


#Sahar original code 
#while (<$pb>) {
#    chomp;
#    my $line=$_;
#    next if ($line =~ /^\</);
#    next if ($line =~ /^\(/);
#    my @line=split(/\s+/,$line);
#    if ($line[0] =~ /^====/){
#	$ko=0;
#	$header=$line[1];
#    }
#    else{
#	if ($ko == 0){
#	    while ($line =~ /(K\d{5})\s+/){
#		$ko_group_data->{$1}->{$header}=1;
#		print "\n*************Associating: $header \t $1\n";
#		$line = "$'";
#		$ko=1;
#	      }
#	}
#    }
#}

print Dumper($ko_group_data);

my %data;
my $output=new FileHandle("> ".$self->output_file);
#Count pathway hits
foreach my $maincat (sort {$a <=> $b} keys%{$tree}){
    #print $output "$maincat\n";
    my %maincat_members_hash;
    my %maincat_ko_hash;
    foreach my $subcat (sort {$a <=> $b} keys%{$tree->{$maincat}}){
	#print $output "$subcat\n";
	my %subcat_members_hash;
	my %subcat_ko_hash;
	foreach my $path (@{$tree->{$maincat}->{$subcat}}){
	    my $num_kos=0;
	    my $ko_list="";
	    my $all_ce_members="";
	    my %ce_members_hash;
	    foreach my $ko (keys%{$pathway_ko_data->{$path}}){
		if ($ko_group_data->{$ko}){
		    $num_kos++;
		    $ko_list.=" $ko ";
		    $subcat_ko_hash{$ko}=1;
		    $maincat_ko_hash{$ko}=1;
		    foreach my $c (keys%{$ko_group_data->{$ko}}){
			$all_ce_members.=" $c ";
			$ce_members_hash{$c}=1;
			$subcat_members_hash{$c}=1;
			$maincat_members_hash{$c}=1;
		    }
		}
	    }
	    if ($num_kos == 0){
		#print $output "$path\t0\t0\-\t-\n";
		#$data{$path}="$path\t0\t0\t-\t-";
		$data{$path}="$path\t0\t0";
		next;

	    }
	    $all_ce_members =~ s/\s+$//;
	    $all_ce_members =~ s/^\s+//;
	    my $ce_print_line="";
	    foreach my $c (keys%ce_members_hash){
		$ce_print_line.="$c ";
	    }
	    my $cnt=keys%ce_members_hash;
	    print $output "$path\t$cnt\t$num_kos\t$ce_print_line\t$ko_list\n";
	    #$data{$path}="$path\t$cnt\t$num_kos\t$ce_print_line\t$ko_list";
	    $data{$path}="$path\t$cnt\t$num_kos";
	}
	 my $subcat_print_line="";
	 my $subcat_ko_line="";
	 foreach my $c (keys%subcat_members_hash){
	     $subcat_print_line.="$c ";
	 }
	 my $subcat_cnt=keys%subcat_members_hash;
	 foreach my $c (keys%subcat_ko_hash){
	     $subcat_ko_line.="$c ";
	 }
	 my $subcat_num_kos=keys%subcat_ko_hash;
	 #$data{$subcat}="$subcat\t$subcat_cnt\t$subcat_num_kos\t$subcat_print_line\t$subcat_ko_line";
	$data{$subcat}="$subcat\t$subcat_cnt\t$subcat_num_kos";
    }
    my $maincat_print_line="";
    my $maincat_ko_line="";
    foreach my $c (keys%maincat_members_hash){
	$maincat_print_line.="$c ";
    }
    my $maincat_cnt=keys%maincat_members_hash;
    foreach my $c (keys%maincat_ko_hash){
	$maincat_ko_line.="$c ";
    }
    my $maincat_num_kos=keys%maincat_ko_hash;
    $data{$maincat}="$maincat\t$maincat_cnt\t$maincat_num_kos\t$maincat_print_line\t$maincat_ko_line";
    $data{$maincat}="$maincat\t$maincat_cnt\t$maincat_num_kos";
}


foreach my $maincat (sort {$a <=> $b} keys%{$tree}){
    print $output "$data{$maincat}\n";
    foreach my $subcat (sort {$a <=> $b} keys%{$tree->{$maincat}}){
	print $output "$data{$subcat}\n";
	foreach my $path (@{$tree->{$maincat}->{$subcat}}){
	    if ($data{$path}){
		print $output "$data{$path}\n";
	    }
	}
    }
}




    return 1;

}
1;
