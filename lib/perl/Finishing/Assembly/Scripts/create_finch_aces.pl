#!/usr/bin/env genome-perl

use strict;
use warnings;
use Data::Dumper;
use Finishing::Assembly::Factory;
use Finishing::Assembly::ContigTools;
use Finishing::Assembly::PhdDB;
use Finishing::Assembly::Phd;
use IO::File;
use IO::Dir;
use Finishing::Assembly::DumpPhd;
use Finfo::Logging;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

my $base_dir = "/gscmnt/temp113/finishing/scratch/finch_merge";
my $assembly_dir = "$base_dir/Assembly_ace_data_for_merges_by_project";
my $input_dir_local ="Single_scaffold_merges";
my $output_dir_local = "Merge_results";
my $result_dir;
my $queue;
my $input_file;
my $bac_name;
my $left_ctg_num;
my $right_ctg_num;
my $left_end_pos;
my $right_end_pos;
my $orientation;
my $on_queue;
my $logfile="$base_dir/run_merges.log";
my $acefile;
my $organism_name='pan_troglodytes';
my $assembly_name='2.1_051011';

my $script_path = "/gscuser/adukes/svn/test_modules/Finishing/Assembly/Scripts/run_finch_merges.pl";

my @args = @ARGV;
GetOptions (
    "queue" => \$queue,
    "input_file=s" => \$input_file,
    "result_dir=s" => \$result_dir,
    "acefile=s" => \$acefile,
    "bac_name=s" => \$bac_name,
    "left_ctg=s"   => \$left_ctg_num,
    "right_ctg=s"  => \$right_ctg_num,
    "left_end_pos=s" => \$left_end_pos,
    "right_end_pos=s" => \$right_end_pos,
    "orientation=s" => \$orientation,
    "organism=s" => \$organism_name,
    "assembly=s" => \$assembly_name,
    "logfile=s" => \$logfile,
    "on-queue" => \$on_queue,
    "input_dir=s" => \$input_dir_local,
    "output_dir=s" => \$output_dir_local,
);

my $input_dir = "$base_dir/$input_dir_local" unless $input_dir_local =~ /^\//;
$input_dir = $input_dir_local unless defined $input_dir;
my $output_dir = "$base_dir/$output_dir_local" unless $output_dir_local =~ /^\//;
$output_dir = $output_dir_local unless defined $output_dir;

if ($queue and !$on_queue){
    exec ("bsub -q long -J RUNMERGES -g /finfo -u adukes perl -I /gscuser/adukes/svn/test_modules $script_path @args --on-queue");
}

@ARGV=();
use GSCApp;
App->init;

if ($queue and $input_file){
    my @jobs;
    my $fh = IO::File->new("< $input_file");
    my $dir = IO::Dir->new($input_dir);
    while (my $line = $fh->getline){
        ($line) = $line =~ /(\S+)/;
        $dir->rewind;
        my $file_path;
        while (defined(my $file = $dir->read)){
            next unless $file =~/$line/ and $file =~/mapbac$/i;
            $file_path = "$input_dir/$file";
            push @jobs, create_job($file_path);
        }
    }
    run_jobs(@jobs);
}

if ($queue and !$input_file)
{
    my @jobs;
    my $dir = IO::Dir->new($input_dir);
    main->fatal_msg("input dir $input_dir doesn't exist") unless defined $dir;
    while(defined(my $file = $dir->read)){
        next unless $file =~ /mapbac$/i;
        my $file_path = "$input_dir/$file";
        main->fatal_msg("Couldn't parse absolute path from $file") unless defined $file_path;
        push @jobs, create_job($file_path);
    }
    run_jobs(@jobs);
}

sub run_jobs{
    my (@jobs) = @_;
    my $scheduler= new PP::JobScheduler(
        job_list => \@jobs,
        day_max => 10,
        night_max => 30,
        refresh_interval => 60,
    );
    $scheduler->start();
    exit;
}

sub create_job{
    my $file_path = shift;
    my $fh = IO::File->new("< $file_path");
    my %ops_hash;
    while (my $line = $fh->getline){
        if ($line =~ /left_ctg|right_ctg|\bbac\b|left_end_pos|right_end_pos|orientation/){
            my ($key, $val) = $line =~ /"(\S+)"\s*=>\s*"?([^"\s,]+)"?/;
            $ops_hash{$key} = $val;
        }
    }
    #my $acefile = $ops_hash{'bac'};
    #my $bac = $acefile;
    #$bac =~ s/\/$//;
    #$bac = basename($acefile);
    #($bac) = $bac =~ /^([^\.]+)\./;
    my $bac = $ops_hash{'bac'};
    my $acefile = Finishing::Assembly::Ace::Dir->new(dir => $ENV{SEQMGR}."/$bac/edit_dir")->recent_acefile;
    my $result_dir = "$output_dir/$bac";
    my $cmd ="perl -I /gscuser/adukes/svn/test_modules $script_path --result_dir $result_dir";
    $cmd .= " --acefile $acefile";
    $cmd .= " --bac_name $bac";
    $cmd .= " --left_ctg ".$ops_hash{'left_ctg'};
    $cmd .= " --right_ctg ".$ops_hash{'right_ctg'};
    $cmd .= " --left_end_pos ".$ops_hash{'left_end_pos'};
    $cmd .= " --right_end_pos ".$ops_hash{'right_end_pos'};
    $cmd .= " --orientation ".$ops_hash{'orientation'};
    main->info_msg($cmd);

    my $pp = undef;

    while (!$pp) {
        $pp = PP->create(
            pp_type => 'lsf',
            q =>'short',
            command => $cmd,
            u => 'adukes@watson.wustl.edu',
            J => "MERGE.$bac",
        );
        if (!$pp){
            warn "Failed to create LSF job for $bac";
            sleep(10);
        }
        else{
            return $pp;
        }
    }
}

#main->fatal_msg("Need input file!") unless defined $input_file;
main->fatal_msg("Need result_dir") unless defined $result_dir;
main->fatal_msg("$result_dir does not exist!") unless -d $result_dir;
#main->fatal_msg("Need left_contig number!") unless defined $left_ctg_num;
#main->fatal_msg("Need right_contig number!") unless defined $right_ctg_num;
#main->fatal_msg("Need left_end_position!") unless defined $left_end_pos;
#main->fatal_msg("Need right_end_position!") unless defined $right_end_pos;
main->fatal_msg("Need orientation!") unless defined $orientation;
main->fatal_msg("Need bac_name!") unless defined $bac_name;
main->fatal_msg("Need acefile!") unless defined $acefile;

main->fatal_msg("Need left or right contig and info to continue!") unless (($left_ctg_num && $left_end_pos) || ($right_ctg_num && $right_end_pos));

main->info_msg("Beginning merge on $bac_name");

if ($orientation == 0){  #####  TODO these need further study, probably a sorting problem
    main->error_msg("Orientation for $bac_name alignment is irregular!");
    exit;
}

my $edit_dir = "$result_dir/edit_dir";
my $phd_dir = "$result_dir/phd_dir";
my $out_acefile = "$edit_dir/$bac_name.merge.ace";
my $split_acefile = "$edit_dir/$bac_name.split.ace";
my $merge_stats_file = "$result_dir/$bac_name.merge_stats";
(mkdir $edit_dir and main->info_msg("$edit_dir created") or main->fatal_msg("couldn't create $edit_dir")) unless -d $edit_dir;
(mkdir $phd_dir and main->info_msg("$phd_dir created") or main->fatal_msg("couldn't create $phd_dir")) unless -d $phd_dir;

#my $dump_phd = Finishing::Assembly::DumpPhd->new(phd_dir => $phd_dir, grab_dir => $bcm_phd_dir);

my $phddb = Finishing::Assembly::PhdDB->new;

my $ct = Finishing::Assembly::ContigTools->new;

#FOR FINCH WE HAVE THE ASSEMBLY CONGTIGS IN ACEFILES
#my $factory = Finishing::Assembly::DBIxFactory->instance;
#$factory->begin_transaction;
#my $factory = Finishing::Assembly::Factory->connect(
#my $organism = $factory->get_organism($organism_name);
#my $assembly = $organism->get_assembly($assembly_name);
#main->fatal_msg("couldn't get organism") unless $assembly;
#`touch temp.ace`;
#my $out_ao = GSC::IO::Assembly::Ace->new(input_file => 'temp.ace');
#main->fatal_msg("couldn't get out ace object") unless $out_ao;

#first we want to perform a split at the sites that fdu indicated
#my ($scaffold_name) = $left_ctg_num =~ /(Contig\d+)\.\d+/;	
#my $scaffold = $assembly->get_scaffold($scaffold_name);
#main->fatal_msg("Couldn't get scaffold!") unless $scaffold;
#my $left_contig = $scaffold->get_contig($left_ctg_num);
#my $right_contig = $scaffold->get_contig($right_ctg_num);

my $assembly_dir_h = IO::Dir->new($assembly_dir);
my $assembly;
my $file;
while ($file = $assembly_dir_h->read and defined $file){
    next unless $file =~ /$bac_name/;
    $assembly = Finishing::Assembly::Factory->connect('ace', "$assembly_dir/$file")->get_assembly;
    last;
}

my $left_contig;
my $right_contig;

if ($left_ctg_num){
    $left_contig = $assembly->get_contig($left_ctg_num);

    unless ($left_end_pos > 0 and $left_end_pos <= $left_contig->length){
        main->fatal_msg("left_end_pos $left_end_pos is out of bounds: ".$left_contig->length);
    }

    $left_end_pos += 1000;
    $left_end_pos = $left_contig->length if $left_end_pos > $left_contig->length;
    #$dump_phd->dump_phd(contig => $left_contig);
}

if ($right_ctg_num){
    $right_contig = $assembly->get_contig($right_ctg_num);
    unless ($right_end_pos > 0 and $right_end_pos <= $right_contig->length){
        main->fatal_msg("right_end_pos $right_end_pos is out of bounds: ".$right_contig->length);
    }

    $right_end_pos -= 1000;
    $right_end_pos = 0 if $right_end_pos < 0;

    #$dump_phd->dump_phd(contig => $right_contig);
}

#####################
#check cut positions against contig boundaries
#####################

#####################
#adjust left_end_pos and right_end_pos to include some overlap area
#####################

my $po = Finishing::Assembly::Phd->new(input_directory => "/home1/watson/seqmgr/$bac_name/phd_dir");
my $output_phd = Finishing::Assembly::Phd->new(input_directory => "$phd_dir");
main->info_msg($phd_dir);
main->fatal_msg("couldn't get phd object") unless $po;
##################
#grab ace contigs#
#my $left_ace_contig_ao = GSC::IO::Assembly::Ace->new(input_file=>$left_contig->acefile);
#my $right_ace_contig_ao = GSC::IO::Assembly::Ace->new(input_file=>$right_contig->acefile);
#my $right_ace_contig = $right_ace_contig_ao->get_contig($right_contig->id);
#my $left_ace_contig = $left_ace_contig_ao->get_contig($left_contig->id);
##################
#next merge the contigs together

my $bac_file = get_bac_file($bac_name, $acefile); #`ls ~seqmgr/$bac_name/edit_dir/$bac_name.fasta.screen.ace.0`;
main->fatal_msg("$bac_file doesn't exist") unless -e $bac_file;

#my $bac_ao = GSC::IO::Assembly::Ace->new(input_file => $bac_file);

my $bac_contig = Finishing::Assembly::Factory->connect('ace', $bac_file)->get_assembly->get_longest_contig;
if ($orientation == -1){
    main->info_msg("Complementing...");
    $bac_contig->complement;
}

main->info_msg("writing out pre-merge acefile");
my $split_ex = Finishing::Assembly::Ace::Exporter->new(file => $split_acefile);
$split_ex->export_contig(contig=>$bac_contig);
if ($left_contig){
    $split_ex->export_contig(contig=>$left_contig, name=>$left_contig->name.".left");
}
if ($right_contig){
    $split_ex->export_contig(contig=>$right_contig, name=>$right_contig->name.".right");
}
$split_ex->close;

my ($left_merge_contig, $right_merge_contig);
if ($left_contig){
    ($left_merge_contig,undef) = $ct->split($left_contig,undef, no_gui=>1, split_position => $left_end_pos, phd_array => [$po,$phddb], allow_no_split => 1);
    main->info_msg('done left');
}
if ($right_contig){
    (undef,$right_merge_contig) = $ct->split($right_contig,undef, no_gui=>1, split_position => $right_end_pos, phd_array => [$po,$phddb], allow_no_split => 1);
    main->info_msg('done right');
}

main->info_msg("writing out pre-merge acefile with split contigs");
my $split_ex_2 = Finishing::Assembly::Ace::Exporter->new(file => $split_acefile.".post_split");
$split_ex_2->export_contig(contig=>$bac_contig);
if ($left_contig){
    $split_ex_2->export_contig(contig=>$left_merge_contig, name=>$left_merge_contig->name.".left");
}
if ($right_contig){
    $split_ex_2->export_contig(contig=>$right_merge_contig, name=>$right_merge_contig->name.".right");
}
$split_ex_2->close;
exit(0);

#$dump_phd->dump_phd(ace_contig => $bac_contig);
unlink $merge_stats_file if -e $merge_stats_file;
`touch $merge_stats_file`;
my $statsfh = IO::File->new(">> $merge_stats_file");
my $merge_contig = $bac_contig;


main->info_msg("merging");
if ($left_merge_contig){
    $merge_contig = $ct->merge($left_merge_contig, $merge_contig,undef,phd_array => [$po,$output_phd, $phddb ], output_phd=>$output_phd, statsfh =>$statsfh);
    main->info_msg("l merge done");
}

if ($right_merge_contig){
    $merge_contig = $ct->merge($merge_contig, $right_merge_contig,undef,phd_array => [$po,$output_phd, $phddb ], output_phd=>$output_phd, statsfh =>$statsfh);
    main->info_msg("l merge done");
}

main->info_msg("merge complete!");
my $exporter = Finishing::Assembly::Ace::Exporter->new(file => $out_acefile);
$exporter->export_contig(contig => $merge_contig, name => "$left_ctg_num".'_x_'."$right_ctg_num");
$exporter->close;

#$out_ao->add_contig($merge_contig);
main->info_msg("contig_added!");

#$out_ao->write_file(output_file => $out_acefile);

sub get_bac_file{
    my ($bac_name, $acefile) = @_;
    return $acefile if -e $acefile;
    return "/home1/watson/seqmgr/$bac_name/edit_dir/$acefile";
}

