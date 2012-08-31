#!/usr/bin/env genome-perl

use strict;

use Getopt::Long;
use Finishing::Assembly::Factory;
use Finishing::Assembly::DBIx::AssemblyImporter;

use PP;
use PP::JobScheduler;
use App::Name;

use IO::File;
use IO::Dir;
use DateTime;

use Finfo::Logging;

my $base_file_path = "/gscmnt/843/finishing/assembly/WholeGenomeImprovement";

my $usage = "\nUSAGE: dbix_schema_assembly_importer.pl -assembly <assembly> -organism <organism> -store-reads -store-tags -store-ace -parallel <acefile1> <acefile2> ...  -or\n dbix_schema_assembly_importer.pl -assembly <assembly> -organism <organism> -store-reads -store-tags -store-ace -parallel -source-dir <acefile directory>  if all the acefiles are in a single directory\n";
$usage .= sprintf('%-18s', "assembly").sprintf('%-5s', "")."name of the assembly\n";
$usage .= sprintf('%-18s', "assembly-file-path").sprintf('%-5s', "")."file path to write contig ace files to(default $base_file_path/<organism>/<assembly>\n";
$usage .= sprintf('%-18s', "organism").sprintf('%-5s', "")."name of the organism\n";
$usage .= sprintf('%-18s', "source-dir").sprintf('%-5s', "")."path to directory with assembly ace files to import\n";
$usage .= sprintf('%-18s', "fof").sprintf('%-5s', "")."fof with paths to acefiles to import\n";
$usage .= sprintf('%-18s', "db").sprintf('%-5s', "")."database, mysql or oracle_dev (default mysql)\n";
$usage .= sprintf('%-18s', "store-reads").sprintf('%-5s', "flag")."store reads in database (default 1)\n";
$usage .= sprintf('%-18s', "store-tags").sprintf('%-5s', "flag")."store tags in database (default 1)\n";
$usage .= sprintf('%-18s', "store-ace").sprintf('%-5s', "flag")."store contig ace files at assembly-file-path (default 1)\n";
$usage .= sprintf('%-18s', "parallel").sprintf('%-5s', "flag")."run in parallel on the queue\n";
$usage .= sprintf('%-18s', "night-max").sprintf('%-5s', "")."max number of simultaneous jobs at night(default 50)\n";
$usage .= sprintf('%-18s', "day-max").sprintf('%-5s', "")."max number of simultaneous jobs during the day(default 30)\n";

my @args = @ARGV; #store these for reconstituting the command-line if using parallel processing
if (@args == 0) {
    main->fatal_msg("Need arguments:$usage");
}

my %db_connections = (
    mysql => [
    'dbi:mysql:cmap_sorting:mysql1',
    'sorting_admin',
    's_admin',
    {AutoCommit => 1}
    ],
    oracle_dev => [
    'dbi:Oracle:dwrac:Assembly',
    'gscuser',
    'g_user',
    {
        AutoCommit => 0,
        LongReadLen => 10000000
    },
    {
        on_connect_do => ["ALTER SESSION SET NLS_DATE_FORMAT = 'YYYY-MM-DD HH24:MI:SS'"],
    }
    ],
);

my (
    $source_dir,
    $organism_name,
    $assembly_name,
    $assembly_file_path,
    $store_ace,
    $store_reads,
    $store_tags,
    $parallel,
    $night_max,
    $day_max,
    $db,
    $queue,
    $on_queue,
    $fof,
);

my $help;

$day_max ||= 5;
$night_max ||= 5;
$db = 'mysql';
$queue ||= '@qblade';
$store_reads ||= 1;
$store_ace ||= 1;
$store_tags ||= 1;


my $result = GetOptions(
    "assembly=s" => \$assembly_name,
    "assembly-file-path" => \$assembly_file_path,
    "organism=s" => \$organism_name,
    'source-dir=s' => \$source_dir,
    'db=s' => \$db,
    'store-reads' => \$store_reads,     ##AceAdaptor params
    'store-tags' => \$store_tags,       ##""
    'store-ace' => \$store_ace,         ##""
    'parallel' => \$parallel,           ##for queue - reconstitute command line
    'q=s'=> \$queue,                    ##'@qblade' by default, don't know any other choices
    'night-max=i' => \$night_max,       ##PP::JobScheduler params
    'day-max=i' => \$day_max,           ##""
    'on-queue' => \$on_queue,
    'help' => \$help,
    'fof=s' => \$fof,
);

print $usage and exit(0) if $help;

unless ($on_queue){
    exec "bsub -q long -J IMPASS -g /finfo -oo /gscuser/adukes/bin/testers/AssemblyImprovement/Output/job_$assembly_name.out -u adukes perl -I /gscuser/adukes/svn/perl_modules -I /gscuser/adukes/svn/perl_modules /gscuser/adukes/svn/perl_modules/Finishing/Assembly/Scripts/assembly_importer.pl @args -on-queue";
}

#build or locate file_path to write contig acefiles to
unless ($assembly_file_path){
    my $org_path = "$base_file_path/$organism_name";
    mkdir $org_path unless -d $org_path;
    $assembly_file_path = "$org_path/$assembly_name";
}
mkdir $assembly_file_path unless -d $assembly_file_path;

#assign default values
main->fatal_msg("Need assembly name! $usage") unless $assembly_name;
main->fatal_msg("Need organism! $usage") unless $organism_name;

my @acefiles;
if (@ARGV){
    foreach (@ARGV){
        if (-e $_){
            push @acefiles, $_;
        }else{
            main->error_msg("Acefile $_ does not exist!");
        }
    }
}elsif ($source_dir){
    main->fatal_msg("Source dir $source_dir does not exist!") unless -d $source_dir;
    my $dir = IO::Dir->new("$source_dir");
    while (defined ($_ = $dir->read) ){
        next if $_ =~ /^\.+$/;
        next if $_ !~ /ace/;
        my $path = "$source_dir/$_";
        push @acefiles, $path;
    }
}elsif($fof){
    my $file = IO::File->new("< $fof");
    while (my $line = $file->getline){
        chomp $line;
        main->fatal_msg("$line isn't an acefile!") unless -e $line;
        push @acefiles, $line;
    }
}else{
    main->fatal_msg("Need a list of acefiles(on command line or using -fof option) to import or the -source-dir containing the acefiles!");
}

main->fatal_msg("No acefiles to import!") unless @acefiles;


#TODO need to create organism and assembly before splitting into parallel to prevent locking over creation
my $db_factory;
if ($db eq 'dw_assembly') {
    $db_factory = Finishing::Assembly::Factory->connect($db);
}
else {
    $db_factory = Finishing::Assembly::Factory->connect('cmap_admin');
}

my ($organism, $assembly);
my $schema = $db_factory->schema;

$schema->txn_do(  #All dbix from here on out
    sub{
        $organism = $schema->create(
            'organism', { name => $organism_name } 
        );
        main->fatal_msg("Can't find or create organism $organism_name") unless $organism;
    }
);

$schema->txn_do(
    sub {
        $assembly = $organism->create_related('assembly', 
            name => $assembly_name,
            organism_id => $organism->id,
            file_path => $assembly_file_path,
        );
    }
);

if ($parallel) {
    my @flags;
    my %seen;

    $seen{$_}=1 foreach @ARGV;

    foreach my $arg (@args) {
        push(@flags, $arg) unless $seen{$arg};  #filter out ace_paths
    }

    @flags = grep {$_ !~ /parallel/} @flags;

    my @jobs;

    my $resource = "\'rusage[db_dw_dev=1, mem=8000]\'";

    my $username = App::Name->real_user_name;

    foreach my $ace_path (@acefiles) {
        my ($acefile) = ($ace_path =~ /\/([^\/]*)$/);
        print "Running $ace_path\n";
        my $pp = undef;
        while (!$pp) {
            $pp = PP->create(
                pp_type => 'lsf',
                q => 'long',
                command => "perl -I /gscuser/adukes/svn/perl_modules -I /gscuser/adukes/svn/perl_modules /gscuser/adukes/svn/perl_modules/Finishing/Assembly/Scripts/assembly_importer.pl @flags $ace_path",
                R => $resource,
                J => "IMPACE",
                u => $username.'@watson.wustl.edu',
                oo => "/gscuser/adukes/bin/testers/AssemblyImprovement/Output/$acefile.out",
            );

            if (!$pp) {
                warn "Failed to create LSF job for $ace_path";
                sleep(10);
            }
            else {
                push @jobs, $pp;
            }
        }
    }
    exit unless @jobs;
    my $scheduler = new PP::JobScheduler(
        job_list => \@jobs,
        day_max => $day_max,
        night_max => $night_max,
        refresh_interval => 120,
    );
    $scheduler->start();
    exit;
}

my $contig_count = 0;

my $importer = new Finishing::Assembly::DBIx::AssemblyImporter(  
    store_reads => $store_reads,
    store_tags  => $store_tags,
    store_ace   => $store_ace,
);

foreach my $ace_path (@acefiles) {

    main->info_msg("importing $ace_path");
    $ace_path = transfer_file($ace_path); # move the ace file to the local /tmp and unzip it if necessary
    my $ace_factory = Finishing::Assembly::Factory->connect('ace', $ace_path);
    my $ace_assembly = $ace_factory->get_assembly;
    $importer->import_ace_assembly(
        ace_assembly=>$ace_assembly,
        dbix_assembly=>$assembly->proxy->source,
    );
    $ace_factory->disconnect;
    unlink $ace_path;
}


sub transfer_file {
    my $path = shift;
    my ($file_name) = $path =~ /([^\/]+)$/;
    my $local_path = "/tmp/ace_files/$file_name";
    unless (-d '/tmp/ace_files') {
        `mkdir /tmp/ace_files`;
    }
    `cp $path $local_path`;
    if ($local_path =~ /\.gz$/) {
        `gunzip $local_path`;
    }
    $local_path =~ s/\.gz//;
    return $local_path;
}

