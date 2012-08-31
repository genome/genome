package Genome::Model::Tools::Nimblegen::CheckArrayDesign;

use strict;
use warnings;

use Genome;
use IO::File;
use List::Util qw/min max/;

class Genome::Model::Tools::Nimblegen::CheckArrayDesign {
    is => 'Command',
    has => [
    nimblegen_probe_bedfile => {
        type => 'String',
        is_optional => 0,
        doc => 'Zip file provided by nimblegen after designing your array',
    },
    design_files => {
        type => 'Csv',
        is_optional => 1,
        doc => 'A comma-delimited list of files that went into the array design',
    },
    design_file_list => {
        type => 'String',
        is_optional => 1,
        doc => 'Absolute path to a file containing list of all the files (1 per line) that went into the array design. If specified, will override design-files option',
    },
    design_summary_outfile => {
        type => 'String',
        is_optional => 1,
        doc => 'An output summary file giving the number of input sites that came from each array design input file as compared to the number of these respective sites covered in the array design. If this option is not defined, STDOUT will be used.',
    },
    design_file_flank => {
        type => 'Integer',
        is_optional => 1,
        default => '0',
        doc => 'The number of bases to pad each side of events from the design files with when looking for overlap from the probe bedfile',
    },
    probe_bedfile_flank => {
        type => 'Integer',
        is_optional => 1,
        doc => 'The number of bases to pad each side of a nimblegen probe region when looking for overlap with the design regions',
        default => '50',
    },
    covered_perc_cutoff => {
        type => 'Integer',
        is_optional => 1,
        doc => 'The minimum percent [1-100] of bases that must be covered in a single target region by design tilage to be counted as "tiled"',
        default => '80',
    },
    ]
};

sub execute {
    my $self = shift;
    $DB::single=1;

    #parse inputs
    my $probe_bed = $self->nimblegen_probe_bedfile;
    my $summary_file = $self->design_summary_outfile;
    my $design_flank = $self->design_file_flank;
    my $probe_bed_flank = $self->probe_bedfile_flank;
    my $cov_perc_cutoff = $self->covered_perc_cutoff;
    if ($cov_perc_cutoff > 100) {
        $self->error_message("You cannot do better than 100 percent, pardner.");
        return;
    }
    if ($cov_perc_cutoff < 1) {
        $self->error_message("If you only care if the sites are covered 0 percent, then why even run this script?");
        return;
    }

    my @design_files;
    if($self->design_file_list) {
        @design_files = process_design_file_list($self->design_file_list);
    } else {
        @design_files = split(/,/,$self->design_files);
    }


    #put the tiled regions from the probe set into a hash
    my %probes;
    my $probe_fh = new IO::File $probe_bed,"r";
    my $track_found = 0;

    while (my $line = $probe_fh->getline) {
        if (($line =~ /track name=tiled_region description="NimbleGen Tiled Regions"/i) ||
            ($line =~ /track name=hg18_tiled_region description="hg18 NimbleGen Tiled Regions"/ )) {
            $track_found = 1;
            next;
        }
        if ($track_found) {
            my ($chr,$start,$stop) = split /\t/,$line;
            my $modstart = $start - $probe_bed_flank;
            my $modstop = $stop + $probe_bed_flank;
            $probes{$chr}{$modstart} = $modstop;
            #$probes{$chr}{$start} = $stop;
        }
    }
    $probe_fh->close;

    #set up summary filehandle and header
    my $sum_fh;
    if(defined $summary_file) {
        $sum_fh = IO::File->new($summary_file,"w");
        unless($sum_fh) {
            $self->error_message("Unable to open file " . $summary_file . " for writing.");
            return;
        }
    }
    else {
        $sum_fh = IO::File->new_from_fd(fileno(STDOUT),"w");
        unless($sum_fh) {
            $self->error_message("Unable to open STDOUT for writing.");
            return;
        }
    }

    my $header = join("\t","Design_File","#_Target_Sites","#_Tiled_Sites_at_".$cov_perc_cutoff."%","%_Tiled_Sites_at_".$cov_perc_cutoff."%","Total_#_Target_Bases","%_Target_Bases_Covered") . "\n";
    print $sum_fh $header;


    #loop through files, writing outputs as you go
    for my $file (@design_files) {
        my $infh = new IO::File $file,"r";

        #set up file to catch uncovered (untiled) sites
        my $uncovered_file = $file . ".not_tiled";
        my $uncov_fh = new IO::File $uncovered_file,"w";

        #variables to record statistics for summary file
        my $sites = '0';
        my $target_bases = '0';
        my $covered_sites = '0';
        my $covered_bases = '0';
        my $is_SV_file = '0';

        #loop through sites in the file
        while (my $line = $infh->getline) {
            chomp $line;

            #check to see if this is an SV file, with two sides of event to be checked on each line
            if ($line =~ /OUTER_START/) { $is_SV_file = '1'; next; }

            #snp or indel file
            if ( !$is_SV_file ) {

                my ($chr,$start,$stop) = split /\t/,$line;
                if ($chr =~ /MT/) { next; } #ignore MT contig sites (build 36)
                if($chr !~ /^chr/){ $chr = "chr".$chr; } #in case the chromosome number is formatted differently
                $sites++;
                my ($site_target_bases,$site_covered_bases,$site_perc_covered) = count_cov_bases(\%probes,$chr,$start,$stop,$design_flank);
                $target_bases += $site_target_bases;
                $covered_bases += $site_covered_bases;
                #see if site was covered to the appropriate percentage
                if ($site_perc_covered >= $cov_perc_cutoff) {
                    $covered_sites++;
                }
                else { print $uncov_fh "$line\n"; }
            }
            #sv file (two sides to each event)
            else {
                my ($id,$chr1,$start1,$stop1,$chr2,$start2,$stop2) = split /\t/,$line;
                if ($chr1 !~ /^chr/){ $chr1 = "chr".$chr1;}
                if ($chr2 !~ /^chr/){ $chr2 = "chr".$chr2;}
                $sites++;
                #left side
                my ($sv_target_bases_left,$sv_covered_bases_left,$sv_perc_covered_left) = count_cov_bases(\%probes,$chr1,$start1,$stop1,$design_flank);
                #right side
                my ($sv_target_bases_right,$sv_covered_bases_right,$sv_perc_covered_right) = count_cov_bases(\%probes,$chr2,$start2,$stop2,$design_flank);
                #record base stats
                $target_bases += $sv_target_bases_left + $sv_target_bases_right;
                $covered_bases += $sv_covered_bases_left + $sv_covered_bases_right;
                #see if site was covered to the appropriate percentage
                if ($sv_perc_covered_left >= $cov_perc_cutoff && $sv_perc_covered_right >= $cov_perc_cutoff) {
                    $covered_sites++;
                }
                else { print $uncov_fh "$line\n"; }
            }

        }
        $infh->close;
        $uncov_fh->close;

        #add line to summary file for this design file
        #Reminder: HEADER: "Design_File" "#_Target_Sites" "#_Tiled_Sites_at_###%" "%_Tiled_Sites_at_###%" "Total_#_Target_Bases" "%_Target_Bases_Covered"
        my $percent_sites_tiled = 0;
        my $percent_bases_tiled = 0;
        if ($sites > 0) { $percent_sites_tiled = sprintf("%.1f",($covered_sites / $sites * 100)); }
        if ($target_bases > 0) { $percent_bases_tiled = sprintf("%.1f",($covered_bases / $target_bases * 100)); }
        print $sum_fh join("\t",$file,$sites,$covered_sites,$percent_sites_tiled,$target_bases,$percent_bases_tiled) . "\n";
    }
    return 1;
}

sub count_cov_bases {
    my ($probes,$chr,$start,$stop,$design_flank) = @_;
    $start -= $design_flank;
    $stop += $design_flank;
    my $target_bases = $stop - $start + 1;
    my $covered_bases = 0;
    #check all the probes
    for my $probe_start (keys %{$probes->{$chr}}) {
        my $probe_stop = $probes->{$chr}->{$probe_start};
        if ($start <= $probe_stop && $stop >= $probe_start) {
            #some overlap exists here
            my $overlap = min($probe_stop,$stop) - max($probe_start,$start) + 1;
            $covered_bases += $overlap;
            #this next 3 lines of logic is to help remove over-counting if more than one probe covers a target site. It is not perfect, but it is more accurate than not having it.
            if ($probe_start < $start && $probe_stop < $stop) { $start = $probe_stop; }
            if ($probe_stop > $stop && $probe_start > $start) { $stop = $probe_start; }
            if ($covered_bases >= $target_bases) { last; }
        }
    }
    if ($covered_bases >= $target_bases) { $covered_bases = $target_bases; }
    my $perc_covered = '0';
    if ($covered_bases > 0) { $perc_covered = $target_bases / $covered_bases * 100; }
    return ($target_bases,$covered_bases,$perc_covered);
}


sub process_design_file_list {
#return a list of files from a file to process;
#returns reference to array of files

    my $file = shift;
    my $fh = IO::File->new($file,"r");
    my @list=();
    while(my $line = $fh->getline) {
        chomp $line;
        next if($line =~ /^\s+$/ || $line =~ /^\#/ || $line =~ /^$/); #ignore comment and blank lines
        push(@list,$line);
    }
    $fh->close;

    #check to see if files exists and readable
    my @bad_files = grep {!-e $_ } @list;
    my @good_files = grep {-e $_ } @list;
    print STDERR "$_ NOT FOUND\n" for(@bad_files);
    return (@good_files);
}

sub help_brief {
    "Check nimblegen array design. Print uncovered sites."
}

sub help_detail {
    "This script takes in a Nimblegen-designed probe .bed file, and also a comma-delimited list of input files used to make the original probe design sent to nimblegen (or a file which lists these input design files 1 per line), and it does two things: 1) The script checks to see how many of the sites sent to Nimblegen ended up on the probe covered at the percentage level specified by the --coverage-perc-cutoff input parameter, and prints out a summary file listing the filename of the design file, the number of sites in the design file, the number and percentage of sites tiled at the specified level from this design file, the total number of target BASES from the design file, and then the overall number and percentage of bases that actually ended up being tiled. 2) The script also takes the original design file and places a file directly next to it called <design_filename>.not_tiled which contains all of the sites from the original design file that did not end up tiled on the probe .bed file at the appropriate percentage level. Anything not covered at the level specified by --coverage-perc-cutoff will show up in this \".not_tiled\" file, and will not be counted as a tiled site, although whatever bases that ARE tiled from sites ending up in the \".not_tiled\" file will still count towards the overall covered base count."
}

1;
