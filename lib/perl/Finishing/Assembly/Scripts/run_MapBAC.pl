#!/usr/bin/env genome-perl

use strict;
use warnings;

use IO::File;
use GSCApp;
use Data::Dumper;
use Assembly::Improvement::MapBAC;
use Getopt::Long;
use Finfo::Logging;
my $org = 'pan_troglodytes';
my $base_dir = '/gscmnt/temp113/finishing/scratch/chimp_merge';
my $script_path = "$base_dir/run_MapBAC.pl";
my $output_dir_local = "merge_out";
my $blastdb = "$base_dir/blast_DB/Chimp2_1.contigs.fasta"; #CHIMP_DB
my $bes_prefix = "WUGSC-";
my $logfile_local = "run_MapBAC.out";
my $on_queue;
my $input_file_local;
my $acefile;
my @acefiles;
my $project_name;



GetOptions (
    "organism" => \$org,
    "input_file=s" => \$input_file_local,  #fof of absolute paths to ace files
    "output_dir=s" => \$output_dir_local, #dir where merge info goes
    "acefile=s" => \$acefile,
    "on_queue" => \$on_queue,
    "blastdb=s" => \$blastdb,
    "logfile=s" => \$logfile_local,
    "prefix=s" => \$bes_prefix,
    "project_name=s"=>\$project_name,
);

my $logfile;
if (defined $logfile_local){
    $logfile = "$base_dir/$logfile_local" unless $logfile_local =~ /^\//;
}else{
    $logfile = "/tmp/runMapBAC.out";
}
my $output_dir = "$base_dir/$output_dir_local" unless $output_dir_local =~ /^\//;
$output_dir = $output_dir_local unless defined $output_dir;
my $input_file;
if (defined $input_file_local){
    $input_file = "$base_dir/$input_file_local" unless $input_file_local =~ /^\//;
    $input_file = $input_file_local unless defined $input_file;
}
unless ($on_queue){
    my $cmd = "bsub -q long -u ".'adukes@watson.wustl.edu'." -J ALGNMRG -g /finfo -oo $logfile perl -I /adukes/svn/perl_modules $script_path --input_file $input_file --output_dir $output_dir --blastdb $blastdb --on_queue --prefix $bes_prefix";
    main->info_msg("executing $cmd");
    exec $cmd;
    exit;
}

@ARGV = ();

use GSCApp;
App->init;

if (defined $input_file and -s $input_file) {
    my $fh = IO::File->new($input_file) or main->fatal_msg("can't open $input_file\n");
    my @jobs;
    #statistics
    my @lines;
    my @acefiles_attempted;
    while (my $line = $fh->getline) {
        my ($ace) = $line =~ /^(\S+)\s+/;
        main->info_msg("Creating job from file line $ace");
        push @lines, $ace;
        my $project_name;
        if ($ace =~ /^\/.*ace/){
            ($project_name) = $ace =~ /(\w+)\.((\S+)\.)+?ace/;
        }else{
            $project_name = $ace;
            my @aces = ProjectWorkBench::Model::Ace::Dir->new(dir=>"/home1/watson/seqmgr/$ace/edit_dir")->acefiles;
            ($ace) = grep { $_ =~ /ace\.0$/ } @aces;
            main->error_msg("$line is not an acefile and can't recover from bac name") and next unless $ace =~/ace\.0/;
        }
        my $pp;
        my $cmd = "perl -I /adukes/svn/perl_modules $script_path --acefile $ace --output_dir $output_dir --on_queue --prefix $bes_prefix --blastdb $blastdb --project_name $project_name";
        while (!$pp) {
            $pp = PP->create(
                pp_type => 'lsf',
                q =>'long',
                command => $cmd, 
                u => 'adukes@watson.wustl.edu',
                oo => "$output_dir/$project_name.out"
            );
            if (!$pp){
                warn "Failed to create LSF job for $ace";
                sleep(10);
            }
            else{
                main->info_msg("job created!");
                push @jobs, $pp;
            }
        }
        push @acefiles_attempted, $ace;
    }
    $fh->close;
    foreach (@jobs){
        main->info_msg($_->command);
    }
    my $scheduler= new PP::JobScheduler(
        job_list => \@jobs,
        day_max => 10,
        night_max => 30,
    );
    $scheduler->start();
    main->info_msg(scalar(@lines)." lines in input file, ".scalar @acefiles_attempted." acefiles attempted");
    
    exit;
}
else {
    push @acefiles, $acefile;
}

main->fatal_msg("need acefile!") unless defined $acefile;
main->fatal_msg("need project_name!") unless defined $project_name;
main->fatal_msg("need output_dir!") unless defined $output_dir;

main->info_msg("performing MapBAC analysis of $project_name from file $acefile");
my $acefile_name = $acefiles[0];
my $match = Assembly::Improvement::MapBAC->new(
organism        => $org,
bac_list        => \@acefiles,
blast_db        => $blastdb,
contiguous      => 1
);

#my @aligns = $match->get_bac_alignments($bac);
#my $ends = $match->get_bac_ends;

my $output_file = "$output_dir/$project_name.merge_out";
my $out  = $match->execute;
my $fail = $match->get_failure;

$output_file.=".fail" if $fail;
my $ofh = IO::File->new("> $output_file");
my $outs = Dumper $out;
my $fails = Dumper $fail;
print $outs;
print $fails;
print $ofh $outs if defined $outs;
print $ofh $fails if defined $fails;

