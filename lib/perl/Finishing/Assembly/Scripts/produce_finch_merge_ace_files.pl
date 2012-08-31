#!/usr/bin/env genome-perl

use strict;
use warnings;
use Data::Dumper;
use IO::Dir;
use IO::File;
use File::Copy;
use Finishing::Assembly::Commands::RemoveAndReplaceAce;

my $project = shift @ARGV;

my %multi_merge =( 
    0=>[qw/ C-50H12 TGAC-1H12/ ],
    2=>[qw/ TGAC-38H12 TGAC-41H12 /],
    4=>[qw/ TGAC-4H12 TGAC-47H12 /],
    10=>[qw/ TGAC-46H12 TGAC-22H12 /],
    12=>[qw/ TGAC-35H12 TGAC-5H12 TGAC-13H12 /],
    15=>[qw/ TGAC-17H12 TGAC-24H12 /],
    19=>[qw/ TGAC-9H12 TGAC-31H12 /],
);
my @ordered_inserts = qw(
C-50H12
TGAC-1H12
TGAC-38H12
TGAC-41H12
TGAC-4H12
TGAC-47H12
TGAC-46H12
TGAC-22H12
TGAC-35H12
TGAC-5H12
TGAC-13H12
TGAC-17H12
TGAC-24H12
TGAC-9H12
TGAC-31H12
);

if ($project){
    run_command($project,1);
    exit(0);
}

my $finch_path = '/gscmnt/temp113/finishing/scratch/finch_merge';
my $finch_ace_path = "$finch_path/Assembly_ace_data_for_inserts";
my $finch_merge_path = "$finch_path/Merge_results";

my $merge_dir = IO::Dir->new($finch_merge_path);

while (my $proj = $merge_dir->read){
    next unless defined $proj;
    next unless $proj =~ /^T/;
    next if grep { $proj =~ $_ } @ordered_inserts;
    run_command($proj);
    print "Executed on $proj\n";
}

foreach my $scaf (keys %multi_merge){
    my @order = @{$multi_merge{$scaf}};
    foreach (@order){
        run_command($_); 
        print "Executed on $_, scaffold $scaf";
    }
}

sub run_command{
    my ($proj, $no_queue) = @_;
    my @jobs;
    if ($no_queue){

        my $merge = "$finch_merge_path/$proj/edit_dir/$proj.merge.ace";
        my $mapbac_file = "$finch_merge_path/$proj/$proj.MapBAC";
        my ($scaffold_num, $start, $stop) = parse_mb($mapbac_file);
        my $assembly_ace = get_ace($scaffold_num);
        backup_ace($assembly_ace);

        my $command = Finishing::Assembly::Commands::RemoveAndReplaceAce->new(
            replacement_ace_file => $merge,
            source_ace_file => $assembly_ace,
            scaffold_num => $scaffold_num,
            start => $start,
            stop => $stop,
            output_file => "$assembly_ace.new",
            log_file => "$finch_path/finch_insert.out",
        );

        $command->execute;
        unling $assembly_ace;
        move ("$assembly_ace.new", $assembly_ace);

    }else{
        push @jobs, create_job($proj);
    }
    run_jobs(@jobs) unless $no_queue;
}

sub run_jobs{
    my (@jobs) = @_;
    my $scheduler= new PP::JobScheduler(
        job_list => \@jobs,
        day_max => 20,
        night_max => 30,
        refresh_interval => 60,
    );
    $scheduler->start();
    exit;
}

sub create_job{
    my $proj = shift;
    my $script_path = '/gscuser/adukes/svn/perl_modules/Finishing/Assembly/Scripts/produce_finch_merge_ace_files';
    my $cmd ="perl -I /gscuser/adukes/svn/perl_modules $script_path $proj";
    main->info_msg($cmd);

    my $pp = undef;

    while (!$pp) {
        $pp = PP->create(
            pp_type => 'lsf',
            q =>'short',
            command => $cmd,
            u => 'adukes@watson.wustl.edu',
            oo => "$finch_path/insert_$proj.log",
            J => "INSERT.$proj",
        );
        if (!$pp){
            warn "Failed to create LSF job for $proj";
            sleep(10);
        }
        else{
            return $pp;
        }
    }
}

sub get_ace{
    my $num = shift;
    my $mod_path = $num%200;
    my $file = "$finch_ace_path/Taeniopygia_guttata-3.2_071213.pcap.scaffold$num.ace";
    die "no file $file" unless -e $file;
    return $file;
}

my %backups;
sub backup_ace{
    my $acefile = shift;
    $backups{$acefile}||=0;
    $backups{$acefile}++;
    my $file =  $acefile.".".$backups{$acefile};
    copy($acefile,$file) unless -e $file;
}

sub parse_mb{
    my $file = shift;
    my $fh = IO::File->new("< $file");
    while (my $line = $fh->getline){
        chomp $line;
        my ($scaf, $scaf2, $start, $stop);
        $line =~ /right_ctg/ 
        ?
        ($scaf, $stop) = $line =~ /Contig(\d+)\.(\d+)/
        : 
        $line =~ /left_ctg/
        ?
        ($scaf2, $start) = $line =~ /Contig(\d+)\.(\d+)/
        :
        next;
        $scaf ||= $scaf2;
        $scaf2 ||= $scaf;
        $stop ||= 'end';
        $start ||= 'end';
        die "scafs $scaf $scaf2 start $start stop $stop" unless $scaf == $scaf2 and !($start eq $stop);
        return ($scaf, $start, $stop);
    }
}



=pod

=head1 NAME
ScriptTemplate - template for new perl script

=head1 SYNOPSIS

=head1 DESCRIPTION 

=cut

#$HeadURL$
#$Id$


