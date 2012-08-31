
package Genome::Model::Tools::Pcap::RunStats;

use strict;
use warnings;

use Genome;
use IO::File;
use ChimpaceObjects;
use Bio::SeqIO;
use Bio::Seq::Quality;
use Bio::Seq::SequenceTrace;
use Cwd;
use Data::Dumper;


class Genome::Model::Tools::Pcap::RunStats {
    is => 'Command',
    has => [
	    dir => {
		         type         => 'String',
			 is_optional  => 1,
			 doc          => "directory"
		    },

	    tier_1 => {
		          type        => 'String',
			  is_optional => 1,
			  doc         => "first tier value"
		       },

	    tier_2 => {
		          type        => 'String',
			  is_optional => 1,
			  doc         => "Second tier value"
		       },

	    major_contig_length => {
		                     type        => 'String',
				     is_optional => 1,
				     doc         => "Major contig length"
				    },

	    major_supercontig_length => {
		                          type        => 'String',
					  is_optional => 1,
					  doc         => "Major supercontig length"
					 },

	    output_file => {
		             type        => 'String',
			     is_optional => 1,
			     doc         => "Stats output file name"
			    },
	    ],
};

sub help_brief
{
    "runs pcap stats";
}

sub help_synopsis
{
    return <<EOS
gmt pcap run-stats
EOS
}

sub help_detail
{
    return <<EOS
gmt pcap run-stats --tier_1 <value> --tier_2 <value_2>
If est genome size is 1,000,000 bases, tier_1 is 500,000, tier_2 700,000
EOS
}

sub execute
{
    my ($self) = shift;

    unless ($self->validate)
    {
	$self->error_message("Failed to validate");
	return;
    }

    #will parse this out better
    unless ($self->run_contig_stats)
    {
	$self->error_message("Failed to run stats");
	return;
    }
    return 1;
}

sub validate
{
    my ($self) = shift;

    # validate directory #

    my $project_dir = cwd ();

    $project_dir = $self->dir if $self->dir;

    $self->error_message ("you must be in edit_dir") and return
	unless $project_dir =~ /edit_dir$/;

    $self->{project_dir} = $project_dir;

    #check for files necessary to run stats

    my @dir_files = glob ("$project_dir/*");
    
    $self->error_message("No input fasta files found") and return
	unless @{$self->{input_fasta_files}} = grep (/fasta\.gz$/, @dir_files);

    $self->error_message ("No input qual files found") and return
	unless @{$self->{input_qual_files}} = grep (/fasta\.qual\.gz$/, @dir_files);
    
#    $self->error_message ("No pcap unused output file") and return
#	unless @{$self->{reads_unused_files}} = grep (/unused\d+$/, @dir_files);

    $self->error_message ("No scaffold ace files") and return
	unless @{$self->{scaffold_ace}} = grep (/\.scaffold\d+\.ace$/, @dir_files);

    #There should only be one scaffold ace file
    $self->error_message ("Multiple scaffold ace files") and return
	if @{$self->{scaffold_ace}} > 1;

    $self->error_message ("No pcap results files") and return
	unless @{$self->{results_file}} = grep (/\.results$/, @dir_files);
    #There should only be one results file
    $self->error_message ("Multiple results files") and return
	if scalar @{$self->{results_file}} > 1;
		    
    $self->error_message ("No contigs.bases file") and return
	unless $self->{contigs_bases_file} = $project_dir.'/contigs.bases';

    $self->error_message ("No contigs.quals file") and return
	unless $self->{contigs_quals_file} = $project_dir.'/contigs.quals';

    $self->error_message ("No reads.placed file") and return
	unless $self->{reads_placed_file} = $project_dir.'/reads.placed';

    $self->error_message ("No reads.unplaced file") and return
	unless $self->{reads_unplaced_file} = $project_dir.'/reads.unplaced';


    #validate tier values

    if ($self->tier_1 or $self->tier_2)
    {
	$self->error_message("must define both tier_1 and tier_2 values") and return
	    unless $self->tier_1 and $self->tier_2;

	$self->error_message("tier_1 value must be less than tier_2 value") and return
	    unless $self->tier_2 > $self->tier_1;

	$self->{tier_one} = $self->tier_1;
	$self->{tier_two} = $self->tier_2 - $self->tier_1;
    }
    else
    {
	#estimate tier values based on contigs.bases file
	my $contig_bases_file = $self->{project_dir}.'/contigs.bases';

	$self->error_message("You must have a contigs.bases file") and return
	    unless my $size = -s $contig_bases_file;

	#tier values estimated is tier_1 50% of contigs.bases file size, %20 for tier_two

	$self->{tier_one} = int ($size * 0.50);
	$self->{tier_two} = int ($size * 0.20);
    }

    #validate major contig and supercontig sizes

    $self->{major_ctg_length} = 500;
    $self->{major_sctg_length} = 500;

    $self->{major_ctg_length} = $self->major_contig_length
	if $self->major_contig_length;

    $self->{major_sctg_length} = $self->major_supercontig_length
	if $self->major_supercontig_length;

    return 1;
}


#clean this up later
sub run_contig_stats
{
    my ($self) = shift;

    my $dir = $self->{project_dir};
    
    my $GC_num = 0;
    my $AT_num = 0;
    my $NX_num = 0;
    my $total_ctg_length = 0;

    my $stats;

    #DON'T NEED THIS
    opendir (D, $dir) or die "Cannot open $dir directory.\n";
    my @files = readdir (D);

    my $fh = IO::File->new("< $self->{contigs_bases_file}");
    while (my $line = $fh->getline)
    {
	next if ($line =~ /^>/);
	chomp $line;
	my @bases = split(//, $line);
	foreach my $b (@bases)
	{
	    $total_ctg_length++;
	    if($b =~ /[g,c]/i) { $GC_num++; }
	    elsif($b =~ /[a,t]/i) { $AT_num++; }
	    else { $NX_num++; }
	}
    }
    $fh->close;

    my $total_input_reads;
    my $total_input_read_bases;
    my $total_Q20_bases;
    my $ave_Q20_bases_per_read;
    my $ave_input_read_length;
    my $unplaced_reads = 0;
    my $total_prefin_reads = 0;
    my $unplaced_prefin_reads = 0;
    my $reads_in_scaf = 0;

    foreach my $file ( @{$self->{input_qual_files}} )
    {
	my $fh = IO::File->new ("zcat $file |");
	while (my $line = $fh->getline)
	{
	    chomp $line;
	    if($line =~ /^>/) 
	    {
		$total_input_reads++;
		my ($read_name) = $line =~ /^>(\S+)/;
		$total_prefin_reads++ if $read_name =~ /_t/;
	    }
	    else
	    {
		my @tmp = split(' ', $line);
		foreach my $qual (@tmp)
		{
		    $total_input_read_bases++;
		    $total_Q20_bases++ if $qual >= 20;
		}
	    }
	}
	$fh->close;
    }

    my @unplaced_reasons = qw/ unused chimera short small repeat /;

    my $unplaced = {};
    my $reads_unplaced_fh = IO::File->new("< $self->{reads_unplaced_file}") || die
	"Can not open reads.unplaced file";
    while (my $line = $reads_unplaced_fh->getline) {
	next if $line =~ /^\s+$/;
	my @tmp = split (/\s+/, $line);
	print "Warning: Incorrect format for line in reads.unplaced file: $line\n" and
	    next unless scalar @tmp == 3;
	my $read_name = $tmp[1];
	my $reason = $tmp[2];
	if (grep (/^$reason$/, @unplaced_reasons)) {
	    $unplaced->{$reason}++;
	}
	else {
	    $unplaced->{unknown}++;
	}
	$unplaced->{total_reads}++;
	
	if ($read_name =~ /_t/) {
	    $unplaced->{prefinish_reads}++;
	}
    }

    $reads_unplaced_fh->close;

    my $unplaced_unused = 0;
    $unplaced_unused = $unplaced->{unused} if
	exists $unplaced->{unused};
    my $unplaced_chimera = 0;
    $unplaced_chimera = $unplaced->{chimera} if
	exists $unplaced->{chimera};
    my $unplaced_small = 0;
    $unplaced_small = $unplaced->{small} if
	exists $unplaced->{small};
    my $unplaced_short = 0;
    $unplaced_short = $unplaced->{short} if
	exists $unplaced->{short};
    my $unplaced_repeat = 0;
    $unplaced_repeat = $unplaced->{repeat} if
	exists $unplaced->{repeat};
    my $unplaced_prefinish_reads = 0;
    $unplaced_prefinish_reads = $unplaced->{prefinish_reads} if
	exists $unplaced->{prefinish_reads};
    my $all_unplaced_reads = 0;
    $all_unplaced_reads = $unplaced->{total_reads} if
	exists $unplaced->{total_reads};
    
#    foreach my $file ( @{$self->{reads_unused_files}} )
#    {
#	my $fh = IO::File->new("< $file");

#	while (my $line = $fh->getline)
#	{
#	    next unless ($line =~ /\S+/);
#	    $unplaced_reads++;
#	    $unplaced_prefin_reads++ if $line =~ /_t/;
#	}
#	$fh->close;
#    }


    my $rp_fh = IO::File->new("< $self->{reads_placed_file}");
    while (my $line = $rp_fh->getline)
    {
	next if $line =~ /^\s+$/;
	$reads_in_scaf++;
    }

    $rp_fh->close;

    $ave_Q20_bases_per_read = int($total_Q20_bases / $total_input_reads + 0.5);
    $ave_input_read_length = int($total_input_read_bases / $total_input_reads + 0.5);
    
#    my $placed_reads_num = $total_input_reads - $unplaced_reads;
#    my $reads_in_singleton = $total_input_reads - $unplaced_reads - $reads_in_scaf;

    my $chaff_rate = int( $all_unplaced_reads*10000/$total_input_reads + 0.5 ) / 100;
    my $Q20_redundancy = int( $total_Q20_bases * 10 / $total_ctg_length ) / 10;


    $stats =  "\n*** SIMPLE READ STATS ***\n".
	  "Total input reads: $total_input_reads\n".
	  "Total input bases: $total_input_read_bases bp\n".
	  "Total Q20 bases: $total_Q20_bases bp\n".
	  "Average Q20 bases per read: $ave_Q20_bases_per_read bp\n".
	  "Average read length: $ave_input_read_length bp\n".
	  "Placed reads: $reads_in_scaf\n".
#	  "Placed reads: ", $total_input_reads - $unplaced_reads, "\n".
#	  "  (reads in scaffold: $reads_in_scaf)\n".
#	  "  (reads in singleton: ", $total_input_reads - $unplaced_reads-$reads_in_scaf, ")\n".
#	  "  (reads in singleton: $reads_in_singleton)\n".    
#	  "Unplaced reads: $unplaced_reads\n".
#	  "Chaff rate: ", int( $unplaced_reads*10000/$total_input_reads + 0.5 ) / 100, "%\n".
          "Unplaced reads: $all_unplaced_reads\n".
          "  (Singletons: $unplaced_unused)\n".
          "  (Chimera: $unplaced_chimera)\n".
          "  (Repat: $unplaced_repeat)\n".
          "  (Short: $unplaced_short)\n".
          "  (Small: $unplaced_small)\n".
	  "Chaff rate: $chaff_rate"."%\n".
#	  "Q20 base redundancy: ", int( $total_Q20_bases * 10 / $total_ctg_length ) / 10, "X\n".
	  "Q20 base redundancy: $Q20_redundancy"."X\n".
	  "Total prefin reads input: $total_prefin_reads\n".
	  "Total prefin reads unused: $unplaced_prefin_reads\n\n";
    
    my $tier1 = $self->{tier_one};
    my $tier2 = $self->{tier_two};

    my $major_ctg_length = $self->{major_ctg_length};
    my $major_sctg_length = $self->{major_sctg_length};


    my $q20bases;

    my $total_ctg_num;
    $total_ctg_length = 0;
    my $max_ctg_length = 0;
    my $major_ctg_num = 0;
    my $ave_ctg_length;
    my $N50_ctg_length;
    my $N50_ctg_num;
    my %ctg_length = ();
    my %ctg_q20base = ();
    
    my %sctg_length = ();
    my %sctg_q20base = ();

    my $total_sctg_num;
    my $ave_sctg_length;
    my $max_sctg_length = 0;
    my $major_sctg_num = 0;
    my $N50_sctg_length;
    my $N50_sctg_num;

    my $tier1sum = 0; my $tier1num = 0;
    my $tier2sum = 0; my $tier2num = 0;
    my $tier3sum = 0; my $tier3num = 0;
    my $largest_tier1 = 0;
    my $largest_tier2 = 0;
    my $largest_tier3 = 0;
    my $cummulative_length = 0;
    my $not_reached_yet = 1;

    my $t1_ctg_base = 0;  my $t2_ctg_base = 0;  my $t3_ctg_base = 0;
    my $t1_q20base = 0;   my $t2_q20base = 0;   my $t3_q20base = 0;

    my $t1_sctg_base = 0; my $t2_sctg_base = 0; my $t3_sctg_base = 0;
    my $t1_sq20base = 0;  my $t2_sq20base = 0;  my $t3_sq20base = 0;

    my ($t1_N50_ctg_num, $t1_N50_ctg_length);
    my ($t2_N50_ctg_num, $t2_N50_ctg_length);
    my ($t3_N50_ctg_num, $t3_N50_ctg_length);
    
    my ($t1_N50_sctg_num, $t1_N50_sctg_length);
    my ($t2_N50_sctg_num, $t2_N50_sctg_length);
    my ($t3_N50_sctg_num, $t3_N50_sctg_length);
    
    my $t1_not_reached_yet = 1;
    my $t2_not_reached_yet = 1;
    my $t3_not_reached_yet = 1;
    
    my $major_ctg_bases = 0;
    my $major_ctg_q20_bases = 0;

    my $major_sctg_bases = 0;
    my $major_sctg_q20_bases = 0;

    my $larger_than_1M_scaf = 0;
    my $larger_than_250K_scaf = 0;
    my $larger_than_100K_scaf = 0;
    my $larger_than_10K_scaf = 0;
    my $larger_than_5K_scaf = 0;
    my $larger_than_2K_scaf = 0;
    my $larger_than_0K_scaf = 0;

    my $sctg_name;
    my $ctg_name;

    my $iq_fh = IO::File->new("< $self->{contigs_quals_file}");
    while (my $line = $iq_fh->getline)
    {
	chomp $line;

	if($line =~ /^>/)
	{
	    $total_ctg_num++;
	    ($ctg_name) = $line =~ /Contig(\d+\.\d+)/;
	    ($sctg_name) = $line =~ /Contig(\d+)\.\d+/;
	}
	else
	{
	    my @tmp = split(' ', $line);
	    foreach my $item (@tmp)
	    {
		if($item >=20) 
		{
		    $q20bases++; 
		    $ctg_q20base{$ctg_name}++;
		    $sctg_q20base{$sctg_name}++;
		}
		$total_ctg_length++;
		$ctg_length{$ctg_name}++;
		$sctg_length{$sctg_name}++;
	    }
	}
    }

    $iq_fh->close;

    $ave_ctg_length = int($total_ctg_length / $total_ctg_num + 0.5);
    my $tier3 = $total_ctg_length - $tier1 - $tier2;

    foreach my $c (sort {$ctg_length{$b} <=> $ctg_length{$a}} keys %ctg_length)
    {
	if($ctg_length{$c} > $major_ctg_length)
	{
	    $major_ctg_num++;
	    $major_ctg_bases += $ctg_length{$c};
	    $major_ctg_q20_bases += $ctg_q20base{$c};
	}
    
	$cummulative_length += $ctg_length{$c};
    
	$N50_ctg_num++ if $not_reached_yet;

	$max_ctg_length = $ctg_length{$c} if $ctg_length{$c} > $max_ctg_length;
	
	if($not_reached_yet && ($cummulative_length >= ($total_ctg_length * 0.50)))
	{
	    $N50_ctg_length = $ctg_length{$c};
	    $not_reached_yet = 0;
	}
    
	if($tier1sum < $tier1)
	{
	    $t1_N50_ctg_num++ if $t1_not_reached_yet;

	    $tier1sum += $ctg_length{$c};
	    $tier1num++;
	    $largest_tier1 = $ctg_length{$c} if($largest_tier1 == 0);
	
	    if($t1_not_reached_yet && ($cummulative_length >= ($tier1 * 0.50)))
	    {
		$t1_N50_ctg_length = $ctg_length{$c};
		$t1_not_reached_yet = 0;
	    }
	    $t1_ctg_base += $ctg_length{$c};
	    $t1_q20base += $ctg_q20base{$c};
	
	}
	elsif($tier2sum < $tier2)
	{
	    $t2_N50_ctg_num++ if $t2_not_reached_yet;
	   	    
	    $tier2sum += $ctg_length{$c};
	    $tier2num++;
	    $largest_tier2 = $ctg_length{$c} if($largest_tier2 == 0);
	    
	    if($t2_not_reached_yet && (($cummulative_length - $tier1) >= ($tier2 * 0.50)))
	    {
		$t2_N50_ctg_length = $ctg_length{$c};
		$t2_not_reached_yet = 0;
	    }
	
	    $t2_ctg_base += $ctg_length{$c};
	    $t2_q20base += $ctg_q20base{$c};
	}
	else
	{
	    $t3_N50_ctg_num++ if $t3_not_reached_yet;
	
	    $tier3sum += $ctg_length{$c};
	    $tier3num++;
	    $largest_tier3 = $ctg_length{$c} if($largest_tier3 == 0);
	
	    if($t3_not_reached_yet && (($cummulative_length - $tier1 - $tier2) >= ($tier3 * 0.50)))
	    {
		$t3_N50_ctg_length = $ctg_length{$c};
		$t3_not_reached_yet = 0;
	    }
	
	    $t3_ctg_base += $ctg_length{$c};
	    $t3_q20base += $ctg_q20base{$c};
	}
    }

#need to iterate through this again for some N50 specific stats
    my $N50_not_yet_reached = 1;
    my $N50_cummulative_length = 0;
    my $maj_ctg_N50_ctg_num = 0;
    my $maj_ctg_N50_ctg_length = 0;
    foreach my $c (sort {$ctg_length{$b} <=> $ctg_length{$a}} keys %ctg_length)
    {
	next unless $ctg_length{$c} > $major_ctg_length;
	$N50_cummulative_length += $ctg_length{$c};
	$maj_ctg_N50_ctg_num++ if $N50_not_yet_reached;
	if ( $N50_not_yet_reached and $N50_cummulative_length >= ($major_ctg_bases * 0.50) )
	{
	    $N50_not_yet_reached = 0;
	    $maj_ctg_N50_ctg_length = $ctg_length{$c};
	}
    }

    my $Q20_base_ratio = int($q20bases*1000/$total_ctg_length)/10;

    $stats .= "\n*** Contiguity: Contig ***\n".
	    "Total contig number: $total_ctg_num\n".
	    "Total contig bases: $total_ctg_length bp\n".
	    "Total Q20 bases: $q20bases bp\n".
#	    "Q20 bases %: ", int($q20bases*1000/$total_ctg_length)/10, "%\n".
	    "Q20 bases %: $Q20_base_ratio"."%\n".
	    "Average contig length: $ave_ctg_length bp\n".
	    "Maximum contig length: $max_ctg_length bp\n".
	    "N50 contig length: $N50_ctg_length bp\n".
	    "N50 contig number: $N50_ctg_num\n\n";

    my $maj_ctg_avg_ctg_length = int ($major_ctg_bases/$major_ctg_num);
    my $maj_ctg_q20_ratio = (int ($major_ctg_q20_bases * 1000 / $major_ctg_bases + 0.5) ) /10;

    $stats .= "Major_contig (> $major_ctg_length bp) number: $major_ctg_num\n".
	    "Major_contig total contig bases: $major_ctg_bases bp\n".
#	   "Major_contig avg contig length: ",int ($major_ctg_bases/$major_ctg_num),"\n".
	    "Major_contig avg contig length: $maj_ctg_avg_ctg_length\n".
	    "Major_contig Q20 bases: $major_ctg_q20_bases bp\n".
#	   "Major_contig Q20 base percent: ",(int($major_ctg_q20_bases * 1000 / $major_ctg_bases + 0.5))/10,"%\n".
	    "Major_contig Q20 base percent: $maj_ctg_q20_ratio"."%\n".
	    "Major_contig N50 contig length: $maj_ctg_N50_ctg_length\n".
	    "Major_contig N50 contig number: $maj_ctg_N50_ctg_num\n\n";

    my $ctg_top_tier_avg_lenth = int($tier1sum/$tier1num + 0.5);
    my $ctg_top_tier_q20_base_ratio = (int($t1_q20base * 1000 /$t1_ctg_base))/10;

    $stats .= "Top tier (up to $tier1 bp): \n".
	    "  Contig number: $tier1num\n".
#	   "  Average length: ", int($tier1sum/$tier1num + 0.5)," bp\n".
	    "  Average length: $ctg_top_tier_avg_lenth bp\n".
	    "  Longest length: $largest_tier1 bp\n".
	    "  Contig bases in this tier: $t1_ctg_base bp\n".
	    "  Q20 bases in this tier: $t1_q20base bp\n".
#	   "  Q20 base percentage: ", (int($t1_q20base * 1000 /$t1_ctg_base))/10, "%\n".
	    "  Q20 base percentage: $ctg_top_tier_q20_base_ratio"."%\n".
	    "  Top tier N50 contig length: $t1_N50_ctg_length bp\n".
	    "  Top tier N50 contig number: $t1_N50_ctg_num\n\n";
    
    my $ctg_mid_tier_value = $tier1 + $tier2;
    my $ctg_mid_tier_avg_length = int($tier2sum/$tier2num + 0.5);
    my $ctg_mid_tier_q20_base_ratio = (int($t2_q20base * 1000 /$t2_ctg_base))/10;

#    $stats .= "Middle tier ($tier1 bp -- ", $tier1+$tier2, " bp): \n".
    $stats .= "Middle tier ($tier1 bp .. $ctg_mid_tier_value bp):\n".
	    "  Contig number: $tier2num\n".
#	   "  Average length: ", int($tier2sum/$tier2num + 0.5)," bp\n".
	    "  Average length: $ctg_mid_tier_avg_length bp\n".
	    "  Longest length: $largest_tier2 bp\n".
	    "  Contig bases in this tier: $t2_ctg_base bp\n".
	    "  Q20 bases in this tier: $t2_q20base bp\n".
#	   "  Q20 base percentage: ", (int($t2_q20base * 1000 /$t2_ctg_base))/10, "%\n".
	    "  Q20 base percentage: $ctg_mid_tier_q20_base_ratio"."%\n".
	    "  Middle tier N50 contig length: $t2_N50_ctg_length bp\n".
	    "  Middle tier N50 contig number: $t2_N50_ctg_num\n\n";

    my $ctg_low_tier_value = $tier1 + $tier2;
    my $ctg_low_tier_avg_length = int($tier3sum/$tier3num + 0.5);
    my $ctg_low_tier_q20_base_ratio = (int($t3_q20base * 1000 /$t3_ctg_base))/10;

#    print "Bottom tier (", $tier1+$tier2, " bp -- end): \n".
    $stats .= "Bottom tier ( $ctg_low_tier_value bp -- end ): \n".
            "  Contig number: $tier3num\n".
#	   "  Average length: ", int($tier3sum/$tier3num + 0.5)," bp\n".
	    "  Average length: $ctg_low_tier_avg_length bp\n".
	    "  Longest length: $largest_tier3 bp\n".
	    "  Contig bases in this tier: $t3_ctg_base bp\n".
	    "  Q20 bases in this tier: $t3_q20base bp\n".
#	   "  Q20 base percentage: ", (int($t3_q20base * 1000 /$t3_ctg_base))/10, "%\n".
	    "  Q20 base percentage: $ctg_low_tier_q20_base_ratio"."%\n".
	    "  Bottom tier N50 contig length: $t3_N50_ctg_length bp\n".
	    "  Bottom tier N50 contig number: $t3_N50_ctg_num\n\n";


#3. Contiguity: Supercontig stats
    $tier1sum = 0; $tier1num = 0;
    $tier2sum = 0; $tier2num = 0;
    $tier3sum = 0; $tier3num = 0;
    $largest_tier1 = 0;
    $largest_tier2 = 0;
    $largest_tier3 = 0;
    $cummulative_length = 0;
    $not_reached_yet = 1;
    
    $t1_not_reached_yet = 1;
    $t2_not_reached_yet = 1;
    $t3_not_reached_yet = 1;

    foreach my $sc (sort {$sctg_length{$b} <=> $sctg_length{$a}} keys %sctg_length)
    {
	$total_sctg_num++;
	
	if($sctg_length{$sc} > 1000000) {$larger_than_1M_scaf++;}
	elsif($sctg_length{$sc} > 250000) {$larger_than_250K_scaf++;}
	elsif($sctg_length{$sc} > 100000) {$larger_than_100K_scaf++;}
	elsif($sctg_length{$sc} > 10000) {$larger_than_10K_scaf++;}	
	elsif($sctg_length{$sc} > 5000) {$larger_than_5K_scaf++;}
	elsif($sctg_length{$sc} > 2000) {$larger_than_2K_scaf++;}
	else {$larger_than_0K_scaf++;}
    
	if($sctg_length{$sc} > $major_sctg_length)
	{
	    $major_sctg_bases += $sctg_length{$sc};
	    $major_sctg_q20_bases += $sctg_q20base{$sc};
	}
    
	$cummulative_length += $sctg_length{$sc};
	
	$N50_sctg_num++ if $not_reached_yet;

	$max_sctg_length = $sctg_length{$sc} if $sctg_length{$sc} > $max_sctg_length;

	$major_sctg_num++ if $sctg_length{$sc} >= $major_sctg_length;
	
	if($not_reached_yet && ($cummulative_length >= ($total_ctg_length * 0.50)))
	{
	    $N50_sctg_length = $sctg_length{$sc};
	    $not_reached_yet = 0;
	}
	
	if($tier1sum < $tier1)
	{
	    $t1_N50_sctg_num++ if $t1_not_reached_yet;	
	
	    $tier1sum += $sctg_length{$sc};
	    $tier1num++;
	    $largest_tier1 = $sctg_length{$sc} if $largest_tier1 == 0;
	
	    if($t1_not_reached_yet && (($cummulative_length) >= ($tier1 * 0.50)))
	    {
		$t1_N50_sctg_length = $sctg_length{$sc};
		$t1_not_reached_yet = 0;
	    }
	
	    $t1_sctg_base += $sctg_length{$sc};
	    $t1_sq20base += $sctg_q20base{$sc};
	}
	elsif($tier2sum < $tier2)
	{
	    $t2_N50_sctg_num++ if $t2_not_reached_yet;	
	    
	    $tier2sum += $sctg_length{$sc};
	    $tier2num++;
	    $largest_tier2 = $sctg_length{$sc} if $largest_tier2 == 0;
	
	    if($t2_not_reached_yet && (($cummulative_length - $tier1) >= ($tier2 * 0.50)))
	    {
		$t2_N50_sctg_length = $sctg_length{$sc};
		$t2_not_reached_yet = 0;
	    }
	
	    $t2_sctg_base += $sctg_length{$sc};
	    $t2_sq20base += $sctg_q20base{$sc};
	}
	else
	{
	    $t3_N50_sctg_num++ if $t3_not_reached_yet;
	    
	    $tier3sum += $sctg_length{$sc};
	    $tier3num++;
	    $largest_tier3 = $sctg_length{$sc} if $largest_tier3 == 0;
	    
	    if($t3_not_reached_yet && (($cummulative_length - $tier1 - $tier2) >= ($tier3 * 0.50)))
	    {
		$t3_N50_sctg_length = $sctg_length{$sc};
		$t3_not_reached_yet = 0;
	    }
		
	    $t3_sctg_base += $sctg_length{$sc};
	    $t3_sq20base += $sctg_q20base{$sc};
	}
    }
    $ave_sctg_length = int($total_ctg_length / $total_sctg_num + 0.5);
    
    $N50_not_yet_reached = 1;
    $N50_cummulative_length = 0;
    my $maj_sctg_N50_ctg_num = 0;
    my $maj_sctg_N50_ctg_length = 0;
    
    foreach my $sc (sort {$sctg_length{$b} <=> $sctg_length{$a}} keys %sctg_length)
    {
	next unless $sctg_length{$sc} > $major_sctg_length;
	$N50_cummulative_length += $sctg_length{$sc};
	$maj_sctg_N50_ctg_num++ if $N50_not_yet_reached;
	if ( $N50_not_yet_reached and $N50_cummulative_length >= ($major_sctg_bases * 0.50) )
	{
	    $N50_not_yet_reached = 0;
	    $maj_sctg_N50_ctg_length = $sctg_length{$sc};
	}
    }

    $stats .= "\n*** Contiguity: Supercontig ***\n".
            "Total supercontig number: $total_sctg_num\n".
	    "Average supercontig length: $ave_sctg_length bp\n".
	    "Maximum supercontig length: $max_sctg_length bp\n".
	    "N50 supercontig length: $N50_sctg_length bp\n".
	    "N50 supercontig number: $N50_sctg_num\n\n";

    my $maj_sctg_avg_length = int($major_sctg_bases / $major_sctg_num);
    my $maj_sctg_q20_base_ratio = (int($major_sctg_q20_bases * 1000 / $major_sctg_bases + 0.5))/10;

    $stats .= "Major_supercontig (> $major_sctg_length bp) number: $major_sctg_num\n".
            "Major_supercontig total supercontig bases: $major_sctg_bases bp\n".
#	   "Major_supercontig avg contig length: ",int($major_sctg_bases / $major_sctg_num),"\n".
	    "Major_supercontig avg contig length: $maj_sctg_avg_length\n".
	    "Major_supercontig Q20 bases: $major_sctg_q20_bases bp\n".
#	   "Major_supercontig Q20 base percent: ", (int($major_sctg_q20_bases * 1000 / $major_sctg_bases + 0.5))/10, "%\n".
	    "Major_supercontig Q20 base percent: $maj_sctg_q20_base_ratio"."%\n".
	    "Major_supercontig N50 contig length: $maj_sctg_N50_ctg_length\n".
	    "Major_supercontig N50 contig number: $maj_sctg_N50_ctg_num\n\n";

    $stats .= "Scaffolds > 1M: $larger_than_1M_scaf\n".
	    "Scaffold 250K--1M: $larger_than_250K_scaf\n".
	    "Scaffold 100K--250K: $larger_than_100K_scaf\n".
	    "Scaffold 10--100K: $larger_than_10K_scaf\n".
	    "Scaffold 5--10K: $larger_than_5K_scaf\n".
	    "Scaffold 2--5K: $larger_than_2K_scaf\n".
	    "Scaffold 0--2K: $larger_than_0K_scaf\n\n";

    my $sctg_top_tier_avg_length = int($tier1sum/$tier1num + 0.5);
    my $sctg_top_tier_q20_base_ratio = (int($t1_sq20base * 1000 /$t1_sctg_base))/10;

    $stats .= "Top tier (up to $tier1 bp): \n".
            "  Supercontig number: $tier1num\n".
#	   "  Average length: ", int($tier1sum/$tier1num + 0.5)," bp\n".
	    "  Average length: $sctg_top_tier_avg_length bp\n".
	    "  Longest length: $largest_tier1 bp\n".
	    "  Contig bases in this tier: $t1_sctg_base bp\n".
	    "  Q20 bases in this tier: $t1_sq20base bp\n".
#	   "  Q20 base percentage: ", (int($t1_sq20base * 1000 /$t1_sctg_base))/10, "%\n".
	    "  Q20 base percentage: $sctg_top_tier_q20_base_ratio"."%\n".
	    "  Top tier N50 supercontig length: $t1_N50_sctg_length bp\n".
	    "  Top tier N50 supercontig number: $t1_N50_sctg_num\n\n";

    my $sctg_mid_tier_value = $tier1 + $tier2;
    my $sctg_mid_tier_avg_length = int($tier2sum/$tier2num + 0.5);
    my $sctg_mid_tier_q20_base_ratio = (int($t2_sq20base * 1000 /$t2_sctg_base))/10;

#    print "Middle tier ($tier1 bp -- ", $tier1+$tier2, " bp): \n".
    $stats .= "Middle tier ($tier1 bp -- $sctg_mid_tier_value bp): \n".
            "  Supercontig number: $tier2num\n".
#	    "  Average length: ", int($tier2sum/$tier2num + 0.5)," bp\n".
	    "  Average length: $sctg_mid_tier_avg_length bp\n".
	    "  Longest length: $largest_tier2 bp\n".
	    "  Contig bases in this tier: $t2_sctg_base bp\n".
	    "  Q20 bases in this tier: $t2_sq20base bp\n".
#	    "  Q20 base percentage: ", (int($t2_sq20base * 1000 /$t2_sctg_base))/10, "%\n".
	    "  Q20 base percentage: $sctg_mid_tier_q20_base_ratio"."%\n".
	    "  Middle tier N50 supercontig length: $t2_N50_sctg_length bp\n".
	    "  Middle tier N50 supercontig number: $t2_N50_sctg_num\n\n";

    my $sctg_low_tier_value = $tier1 + $tier2;
    my $sctg_low_tier_avg_length = int($tier3sum/$tier3num + 0.5);
    my $sctg_low_tier_q20_base_ratio = (int($t3_sq20base * 1000 /$t3_sctg_base))/10;

#    print "Bottom tier (", $tier1+$tier2, " bp -- end): \n".
    $stats .= "Bottom tier ( $sctg_low_tier_value bp .. end): \n".
	    "  Supercontig number: $tier3num\n".
#	   "  Average length: ", int($tier3sum/$tier3num + 0.5)," bp\n".
	    "  Average length: $sctg_low_tier_avg_length bp\n".
	    "  Longest length: $largest_tier3 bp\n".
	    "  Contig bases in this tier: $t3_sctg_base bp\n".
	    "  Q20 bases in this tier: $t3_sq20base bp\n".
#	   "  Q20 base percentage: ", (int($t3_sq20base * 1000 /$t3_sctg_base))/10, "%\n".
	    "  Q20 base percentage: $sctg_low_tier_q20_base_ratio"."%\n".
	    "  Bottom tier N50 supercontig length: $t3_N50_sctg_length bp\n".
	    "  Bottom tier N50 supercontig number: $t3_N50_sctg_num\n\n";

#############
#constraints#
#############

    my $results_file = ${$self->{results_file}}[0];
    $stats .= "\n*** Constraints ***\n".
	    `tail -12 $results_file`;

    my $GC_ratio = int(1000 * $GC_num/$total_ctg_length + 0.5) / 10;
    my $AT_ratio = int(1000 * $AT_num/$total_ctg_length + 0.5) / 10;
    my $NX_ratio = int(1000 * $NX_num/$total_ctg_length + 0.5) / 10;


    $stats .= "\n\n*** Genome Contents ***\n".
#           "Total GC count: $GC_num, (".int(1000 * $GC_num/$total_ctg_length + 0.5) / 10 ."%)\n".
	    "Total GC count: $GC_num, (".$GC_ratio."%)\n".
#	    "Total AT count: $AT_num, (".int(1000 * $AT_num/$total_ctg_length + 0.5) / 10 ."%)\n".
	    "Total AT count: $AT_num, (".$AT_ratio."%)\n".
#	    "Total NX count: $NX_num, (".int(1000 * $NX_num/$total_ctg_length + 0.5) / 10 ."%)\n".
	    "Total NX count: $NX_num, (".$NX_ratio."%)\n".
	    "Total: $total_ctg_length\n\n";


#########################
# FOSMID COVERAGE STATS #
#########################

#THIS IS NOT NECESSARY NOW BUT MAY NEED TO BE BROUGHT BACK

# if this is brought back change ARGV to a property with shell_args_position

#   my $fosmid_prefix = $ARGV[2];
#    my $fosmid_prefix;

#    if($fosmid_prefix)
#    {
#	open(R, "> readinfo.txt") or die "cannot create readinfo.txt.\n";
#	foreach my $file (@files)
#	{
#	    next unless ($file =~ /scaffold\d+\.ace/);
#	    $file = $dir.$file;
#	    if($file =~ /gz$/) {open(F, "zcat $file |") or die "Cannot open $file.\n";}
#	    else { open(F, $file) or die "Cannot open $file.\n"; }
#	    my $contig;
#	    my %read_info = ();
#	    while(<F>)
#	    {
#		chomp($_);
#		if ($_ =~ /^CO (\S+) /)
#		{
#		    $contig = $1;
#		    next;
#		}
#		if ($_ =~ /^AF /)
#		{ 
#		    my @t = split(' ', $_);
#		    my $read = $t[1];
#		    $read_info{$read}{UorC} = $t[2];
#		    $read_info{$read}{StartPos} = $t[3];
#		    $read_info{$read}{Ctg} = $contig;
#		}
#		if ($_ =~ /^RD /)
#		{
#		    my @t = split(' ', $_);
#		    my $read = $t[1];
#		    $read_info{$read}{Len} = $t[2];
#		}
#	    }
#	    # print out contents in this file.
#	    foreach my $r (sort keys %read_info)
#	    {
#		print R "$r $read_info{$r}{Ctg} $read_info{$r}{UorC} $read_info{$r}{StartPos} $read_info{$r}{Len}\n";
#	    }
#	}
#	
#	my $used_reads = 0;
#	my $cum_len = 0;
#	print "\n*** Fosmid Coverage ***\n";
#	open(F, "readinfo.txt") or die "cannot open readinfo.txt\n";
#	while(<F>)
#	{
#	    next unless ($_ =~ /^$fosmid_prefix/);
#	    chomp();
#	    my @a = split(/\s+/, $_);
#	    
#	    $used_reads++;
#	    $cum_len += $a[4];
#	}
#	print "total fosmids used: $used_reads\n";
#	print "total length of these fosmids: $cum_len bp\n";
#	my $fos_cov = int(100 * $cum_len * $ave_Q20_bases_per_read /$ave_input_read_length/$total_ctg_length)/100;
#	
#	print "Estimated fosmid coverage (total fosmids Q20 bases over genome): ", $fos_cov, "X\n";
#    }

#####################
# 5X COVERAGE STATS #
#####################

#check to make sure there is only one ace file
#my @aces = glob ("*scaffold*ace");

    my $ace_file = ${$self->{scaffold_ace}}[0];

    my $depth_out_file = $ace_file.'_base_depths';

    #create a new depth file only if one doesn't already exist

    if (! -s $depth_out_file)
    {
	my $depth_fh = IO::File->new(">$depth_out_file");
	my $ace_fh = IO::File->new("<$ace_file");
	my $ace_obj = ChimpaceObjects->new ( -acefilehandle => $ace_fh );
	my %cover_depth;
	my @contig_names = $ace_obj->get_all_contigNames ();
	foreach my $contig (@contig_names)
	{
	    my $ctg_length = $ace_obj->Contig_Length ( -contig => $contig );

	    #CONSIDER ONLY MAJOR CONTIGS
	    next unless $ctg_length > 500;

	    my @all_reads = $ace_obj->ReadsInContig ( -contig => $contig );
	    foreach my $read (@all_reads)
	    {
		my $align_unpadded_clip_start = $ace_obj->getAlignUnpadedClipStart ( -read => $read );
		my $align_unpadded_clip_end = $ace_obj->getAlignUnpadedClipEnd ( -read => $read );

		next if (! defined $align_unpadded_clip_end or ! defined $align_unpadded_clip_start);
		my $start = $align_unpadded_clip_start - 1;
		my $end = $align_unpadded_clip_end - 1;

		#if read has a left end overhang, $start value is negative
		#and this causes the array assignment below to crash
		$start = 0 if $start < 0;

		for my $i ( $start .. $end )
		{
		    $cover_depth{$contig}[$i]++;
		}
	    }
	    my $contig_length = $ace_obj->Contig_Length ( -contig => $contig );
	
	    $depth_fh->print (">$contig\n");

	    for ( my $i = 0; $i < $contig_length; $i++ )
	    {
		$depth_fh->print ("$cover_depth{$contig}[$i] ")
		    if defined $cover_depth{$contig}[$i];
	    
		$depth_fh->print ("\n") if ( ($i % 50) == 49);
	    }
	    
	    $depth_fh->print ("\n");

	    delete $cover_depth{$contig};
	}

	$ace_fh->close;

	$depth_fh->close;
    }
    my $total_read_bases = 0;

    my $one_x_cov = 0;
    my $two_x_cov = 0;
    my $three_x_cov = 0;
    my $four_x_cov = 0;
    my $five_x_cov = 0;

    my $df_fh = IO::File->new("<$depth_out_file");
    while (my $line = $df_fh->getline)
    {
	next if $line =~ /^>/;
	next if $line =~ /^\s+$/;
	chomp $line;
	my @numbers = split (/\s+/, $line);
	$total_read_bases += scalar @numbers;
	foreach my $depth_num (@numbers)
	{
	    $one_x_cov++ if $depth_num >= 1;
	    $two_x_cov++ if $depth_num >= 2;
	    $three_x_cov++ if $depth_num >= 3;
	    $four_x_cov++ if $depth_num >= 4;
	    $five_x_cov++ if $depth_num >= 5;
	}
    }
    $df_fh->close;

    my $five_x_cov_ratio = 0;
    my $four_x_cov_ratio = 0;
    my $three_x_cov_ratio = 0;
    my $two_x_cov_ratio = 0;
    my $one_x_cov_ratio =0;

    $five_x_cov_ratio = $five_x_cov/$total_read_bases if $five_x_cov > 0;
    $four_x_cov_ratio = $four_x_cov/$total_read_bases if $four_x_cov > 0;
    $three_x_cov_ratio = $three_x_cov/$total_read_bases if $three_x_cov > 0;
    $two_x_cov_ratio = $two_x_cov/$total_read_bases if $two_x_cov > 0;
    $one_x_cov_ratio = $one_x_cov/$total_read_bases if $one_x_cov > 0;

    $stats .= "\n*** Read Depth Info ***\n".
	      "Total consensus bases: $total_read_bases\n".

#	  "Depth >= 5: $five_x_cov\t". $five_x_cov/$total_read_bases."\n".
#	  "Depth >= 4: $four_x_cov\t". $four_x_cov/$total_read_bases."\n".
#	  "Depth >= 3: $three_x_cov\t". $three_x_cov/$total_read_bases."\n".
#	  "Depth >= 2: $two_x_cov\t". $two_x_cov/$total_read_bases."\n".
#	  "Depth >= 1: $one_x_cov\t". $one_x_cov/$total_read_bases."\n";
 	  
	      "Depth >= 5: $five_x_cov\t $five_x_cov_ratio\n".
	      "Depth >= 4: $four_x_cov\t $four_x_cov_ratio\n".
	      "Depth >= 3: $three_x_cov\t $three_x_cov_ratio\n".
	      "Depth >= 2: $two_x_cov\t $two_x_cov_ratio\n".
	      "Depth >= 1: $one_x_cov\t $one_x_cov_ratio\n\n";
 
   

#####################
# 5 KB CONTIGS INFO #
#####################

    my $total_ctg_lengths = 0;
    my $five_kb_ctg_lengths = 0;

    my $cb_fh = IO::File->new("< $self->{contigs_bases_file}");
    my $fio = Bio::SeqIO->new(-format => 'fasta', -fh => $cb_fh);

    while (my $f_seq = $fio->next_seq)
    {
	my $ctg_length = length $f_seq->seq;
	$total_ctg_lengths += $ctg_length;
	if ($ctg_length >= 5000)
	{
	    $five_kb_ctg_lengths += $ctg_length;
	}
    }
    $cb_fh->close;

    my $ratio = 0;

    if ($five_kb_ctg_lengths > 0)
    {
	$ratio = int ($five_kb_ctg_lengths / $total_ctg_lengths * 100);
    }

    $stats .= "\n*** 5 Kb and Greater Contigs Info ***\n".
            "Total lengths of all contigs: $total_ctg_lengths\n".
	    "Total lengths of contigs 5 Kb and greater: $five_kb_ctg_lengths\n".
	    "Percentage of genome: $ratio%\n";

    if ($self->output_file)
    {
#	my $stats_file_name = $dir.'/'.$self->output_file;
	my $stats_file = $self->output_file;
	my $out_fh = IO::File->new(">$stats_file");
	$out_fh->print("$stats");
	$out_fh->close;
	return 1;
    }

    print $stats;

    return 1;
}

1;
