package Genome::Model::Tools::Ber::BerRunFinish;

use strict;
use warnings;

use Genome;
use Command;

use Carp;
use English;

use BAP::DB::Sequence;
use BAP::DB::SequenceSet;
use BAP::DB::CodingGene;
use Ace;
use Ace::Object;
use Ace::Sequence;

use Data::Dumper;
use IO::File;
use IPC::Run qw/ run timeout /;
use Time::HiRes qw(sleep);
use DateTime;
use Genome::Utility::Email;

use File::Slurp;    # to replace IO::File access...
use File::Copy;
use File::Basename;
use File::stat;
use File::Find::Rule;
use File::stat;

use Bio::SeqIO;

use Cwd;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [
        'locus_tag' => {
            is  => 'String',
            doc => "Locus tag for project, containing DFT/FNL postpended",
        },
        'outdirpath' => {
            is  => 'String',
            doc => "output directory for the ber product naming software",
        },
        'sqlitedatafile' => {
            is  => 'String',
            doc => "Name of sqlite output .dat file",
        },
        'sqliteoutfile' => {
            is  => 'String',
            doc => "Name of sqlite output .out file",
        },
        'acedb_version' => {
            is  => 'String',
            doc => "Current acedb version",
        },
        'amgap_path' => {
            is  => 'String',
            doc => "Current path to AMGAP data",
        },
        'pipe_version' => {
            is  => 'String',
            doc => "Current pipeline version running",
        },
        'project_type' => {
            is  => 'String',
            doc => "Current project type",
        },
        'org_dirname' => {
            is  => 'String',
            doc => "Current organism directory name",
        },
        'assembly_name' => {
            is  => 'String',
            doc => "Current assembly name",
        },
        'assembly_version' => {
            is  => 'String',
            doc => "Current assembly version",
        },
        'sequence_set_id' => {
            is  => 'String',
            doc => "Current sequence set id",
        },
        'locustagdir' => {
            is  => 'String',
            doc => "Locus tag directory",
        },
        'berbasedir' => {
            is  => 'String',
            doc => "BER base directory",
        },
    ]
);

sub help_brief
{
    "Tool for making the final BER ace file, write new parse script and parse into acedb via tace, gather stats from phase5 ace file and acedb, writes the rt file and mails when finished ";
}

sub help_synopsis
{
    return <<"EOS"
      Tool for making the final BER ace file, write new parse script and parse into acedb via tace, gather stats from phase5 ace file and acedb, writes the rt file and mails when finished.
EOS
}

sub help_detail
{
    return <<"EOS"
Tool for making the final BER ace file, write new parse script and parse into acedb via tace, gather stats from phase5 ace file and acedb, writes the rt file and mails when finished.
EOS
}

sub execute
{
    my $self          = shift;
    my $locus_tag     = $self->locus_tag;
    my $outdirpath    = $self->outdirpath;
    my $sqlitedata    = $self->sqlitedatafile;
    my $sqliteout     = $self->sqliteoutfile;
    my $acedb_ver     = $self->acedb_version;
    my $amgap_path    = $self->amgap_path;
    my $pipe_version  = $self->pipe_version;
    my $project_type  = $self->project_type;
    my $org_dirname   = $self->org_dirname;
    my $assembly_name = $self->assembly_name;
    my $assembly_version = $self->assembly_version;
    my $ssid          = $self->sequence_set_id;
    my $locustagdir   = $self->locustagdir;
    my $berbasedir   = $self->berbasedir;

	my %mgap_genome = ();
	$mgap_genome{ 'db_version_id' } = 1;
	$mgap_genome{ 'locus_name' } = $locus_tag;
	$mgap_genome{ 'assembly_name' } = $assembly_name;
	$mgap_genome{ 'org_dirname' } = $org_dirname;

	$mgap_genome{ 'project_type' } = 1 if ($project_type eq 'HGMI');
	$mgap_genome{ 'project_type' } = 2 if ($project_type eq 'CORE');

	$mgap_genome{ 'pipe_version' } = 1 if ($pipe_version eq 'Version_1.0');
	$mgap_genome{ 'pipe_version' } = 2 if ($pipe_version eq 'Version_2.0');

	$mgap_genome{ 'ssid' } = $ssid;
	$mgap_genome { 'acedb_ver' } = $acedb_ver;

	if ($locus_tag =~ m /(DFT|DFT\d|FNL|FNL\d|MSI|MSI\d)$/) {
		$mgap_genome{ 'run_type' } = 1;
	} 

	if ($locus_tag =~ m /(TST|TST\d)$/) {
		$mgap_genome{ 'run_type' } = 2;
	} 

	if (Genome::Sys->username)
	{
		$mgap_genome{ 'user' } = Genome::Sys->username;
	}
	elsif (exists($ENV{LOGIN}) )
	{
		$mgap_genome{ 'user' } = $ENV{LOGIN};
	}
	elsif (exists($ENV{USERNAME}) )
	{
		$mgap_genome{ 'user' } = $ENV{USERNAME};
	}

	#populate_mgap_genome_db (\%mgap_genome);
	#exit;

	my $anno_submission = $amgap_path . "/"
		. $org_dirname . "/"
		. $assembly_name . "/"
		. $assembly_version . "/"
		. "Genbank_submission/"
		. $pipe_version . "/Annotated_submission/";

    my $program = "/gsc/scripts/bin/tace";
    my $cwd     = getcwd();
    my $outdir
		#= qq{/gscmnt/ams1102/info/annotation/BER/autoannotate_v2.5/out};
		#= qq{/gscmnt/sata835/info/annotation/BER/autoannotate_v2.5/out/};
		= $berbasedir . "/out/";
    unless ( $cwd eq $outdir )
    {
        chdir($outdir)
            or die
            "Failed to change to '$outdir'...  from BerRunFinish.pm: $OS_ERROR\n\n";
    }

    my $sqlitedatafile = qq{$outdirpath/$sqlitedata};
    # recent versions of BER try to be helpful and append '.dat' to dat files
    # then we get .dat.dat on the ends of the files.
    if(-f $sqlitedatafile.".dat") {
        # rename file.
        rename($sqlitedatafile.".dat",$sqlitedatafile);
    }
    
    ## get the latest filename
    (my $sqlitedatafilename, $outdirpath) = fileparse($sqlitedatafile);

    ## Before copying sqlite file we will need to delete sqlite*.dat && *.dat.dat file in the Annotated_submission directory if already exists
	my @sqliteDataFiles = find (
							file =>
							name => [ qw/ *.dat *.dat.dat / ],
							in	 => $anno_submission
						  );

    ## Copy sqlite.dat file Annotated_submission directory
	unlink $_ or croak qq{ \n\n Error removing file $_: $! \n\n }
		for (@sqliteDataFiles);

    copy($sqlitedatafile, $anno_submission.$sqlitedatafilename) || croak qq{\n\n Copying of $sqlitedatafile to $anno_submission.$sqlitedatafilename failed ...  from BerRunFinish.pm: $OS_ERROR\n\n };

    my $sqliteoutfile  = qq{$outdirpath/$sqliteout};
    unless ( ( -e $sqlitedatafile ) and ( !-z $sqlitedatafile ) )
    {
        croak
            qq{\n\n NO file,$sqlitedatafile, found for or empty ... from BerRunFinish.pm: $OS_ERROR\n\n };
    }

    my $acedb_short_ver = $self->version_lookup( $self->acedb_version );
    my $acedb_data_path = $self->{amgap_path}
        . "/Acedb/"
        . $acedb_short_ver
        . "/ace_files/"
        . $self->locus_tag . "/"
        . $self->pipe_version;
    unless ( -d $acedb_data_path )
    {
        croak
            qq{\n\n NO acedb_dir_path, $acedb_data_path, found ... from BerRunFinish.pm: $OS_ERROR\n\n };
    }
    ################################
    # parse the sqlite data file
    ################################

    my $bpnace_fh = IO::File->new();
    my $bpname    = qq{_BER_product_naming.ace};
    my $bpn_file  = qq{$acedb_data_path/$locus_tag$bpname};
    $bpnace_fh->open("> $bpn_file")
        or die
        "Can't open '$bpn_file', bpn_file for writing ...from BerRunFinish.pm: $OS_ERROR\n\n";

    my $data_fh = IO::File->new();
    $data_fh->open("< $sqlitedatafile")
        or die
        "Can't open '$sqlitedatafile',sqlite data file ...from BerRunFinish.pm: $OS_ERROR\n\n";

    # FIXME
    # this could be cleaned up with read_file/write_file and Text::CSV_XS...
    while (<$data_fh>)
    {
        my ($featname,   $proteinname, $genesymbol, $go_terms,
            $ec_numbers, $t_equivalog, $speciesname
            )
            = split( /\t/, $ARG );
        print $bpnace_fh "Sequence ", $featname, "\n";
        print $bpnace_fh "BER_product " , "\"$proteinname\"", "\n\n";
    }
    $bpnace_fh->close();
    $data_fh->close();
    ########################################################
    # check acedb readlock before proceeding
    ########################################################

    $cwd = getcwd();

    # is $self->{amgap_path} legal here???
    my $acedb_readlock_path = $self->{amgap_path}
        . "/Acedb/"
        . $acedb_short_ver
        . "/database/readlocks";
    unless ( $cwd eq $acedb_readlock_path )
    {
        chdir($acedb_readlock_path)
            or die
            "Failed to change to '$acedb_readlock_path'...  from BerRunFinish.pm: $OS_ERROR\n\n";
    }

    #FIXME
    # oh, shit, not again.
    # File::Find on the acedb_readlock_path, and push into @readlock
    # this should also be in it's own subroutine for easier testing
    while (1)
    {
        opendir( DIR, $acedb_readlock_path )
            or die
            "Can't open '$acedb_readlock_path'... from BerRunFinish.pm: $OS_ERROR\n\n";
        my @readlock = ();
        while ( defined( my $file = readdir DIR ) )
        {
            next if $file =~ m/^\.\.?$/;
            push( @readlock, $file );
        }
        closedir(DIR);
        my ( $readlock_fh, $session, $machine, $gsc, $wustl, $edu, $process );

	    my $current_time = time;
		my $index = 1;
        if (@readlock)
        {
			foreach my $lockfile (@readlock)
			{
				my $stat = stat($lockfile);
				#my $diff = ($current_time - $stat->mtime) / 3600;
				my $diff = ($current_time - $stat->mtime) / 60;

				#if ($diff > 2) { ## 2 hours
				if ($diff > 15) { ## 15 minutes
					$self->status_message("Removing lockfile: $lockfile");
					unlink $lockfile or die "Can't delete lock file: $lockfile...- $OS_ERROR\n\n";
				} else {
					( $session, $machine, $gsc, $wustl, $edu, $process )
						= split /\./, $lockfile;
					$machine = join( '.', $machine, $gsc, $wustl, $edu );
					$readlock_fh = IO::File->new();
					$readlock_fh->open("< $lockfile")
						or die
						"Can't open '$readlock_fh', readlock_fh for reading ...from BerRunFinish.pm: $OS_ERROR\n\n";
				}
			}

			if (defined ($readlock_fh)) {

					my %readlock = ();

					while (<$readlock_fh>)
					{
						chomp $ARG;
						if (   $ARG !~ /^[\#\s\t\n]/
								&& $ARG =~ /^([^\:]+)\:([^\n\#]*)/ )
						{
							my $key   = $1;
							my $value = $2;

							#remove preceding and trailing whitespace
							$value =~ s/^[\s\t]+|[\s\t]+$//g;

							#set value
							$readlock{$key} = defined($value) ? $value : '';
						}
					}
					print
						qq{\n\nACeDB Session number: $session is currently readlocked by: $readlock{User} on $readlock{Created} using machine:$machine (process ID is: $process)\n};
					sleep(300);
					next;
			}
        }
        else
        {
            print qq{\n No readlock detected, continuing to next step....\n};
            last;
        }
    }
    ########################################################
    # check acedb via tace for previous data and remove it
    ########################################################
    # this should be in a separate sub too...
    $cwd = getcwd();
    my $acedb_maindir_path
        = $self->{amgap_path} . "/Acedb/" . $acedb_short_ver;
    unless ( $cwd eq $acedb_maindir_path )
    {
        chdir($acedb_maindir_path)
            or die
            "Failed to change to '$acedb_maindir_path'...  from BerRunFinish.pm: $OS_ERROR\n\n";
    }

    #connecting to acedb database
    my $db = Ace->connect(
        -path    => "$acedb_maindir_path",
        -program => "$program"
        )
        or die
        "ERROR: cannot connect to acedb ... from BerRunFinish.pm: $OS_ERROR\n \n";

    #mining data from acedb
    my @ace_objects = $db->fetch(
        -name  => "$locus_tag*",
        -class => 'Sequence',
    );

    #my @ace_objects = $db->fetch(
    #	                           -name    => "$locus_tag*",
    #				   -class   => 'Sequence',
    #				   -filltag => 'Visible',
    #			          );

    foreach my $ace_stuff (@ace_objects)
    {
        my $aceseq_obj = $db->fetch( 'Sequence' => $ace_stuff );
        my $BER_product = $aceseq_obj->BER_product();
        if ( defined($BER_product) )
        {
            my $result_code = $aceseq_obj->kill();
        }
        else
        {
            next;
        }
    }
    print qq{\n Done checking/removing previous data from ACeDB\n\n};

    ########################################################
    # write new parse script and parse into acedb via tace
    ########################################################
    # separate subroutine again...
    $cwd                = getcwd();
    $acedb_maindir_path = $self->{amgap_path} . "/Acedb/" . $acedb_short_ver;
    unless ( $cwd eq $acedb_maindir_path )
    {
        chdir($acedb_maindir_path)
            or die
            "Failed to change to '$acedb_maindir_path'...  from BerRunFinish.pm: $OS_ERROR\n\n";
    }

    my $acedb_scripts_path = $self->{amgap_path} . "/Acedb/Scripts";
    my $parse_script_name
        = "parsefiles_pap_ber_" . $locus_tag . "_" . $pipe_version . ".sh";
    my $parse_script_full = qq{$acedb_scripts_path/$parse_script_name};
    my $parse_script_fh   = IO::File->new();
    $parse_script_fh->open("> $parse_script_full")
        or die
        "Can't open '$parse_script_full', parse_script_full for writing ...from BerRunFinish.pm: $OS_ERROR\n\n";

    opendir( ACEDATA, $acedb_data_path )
        or die
        "Can't open $acedb_data_path, acedb_data_path from BerRunFinish.pm: $OS_ERROR\n";

    my @acefiles = ();
    @acefiles = readdir(ACEDATA);
    closedir(ACEDATA);

    my $phase5file = $locus_tag . "_phase_5_ssid_";

    my $parse = "parse";
    print $parse_script_fh "#!/bin/sh -x\n\n";
    print $parse_script_fh
        "#if you call script from bash, tace will follow links!\n\n";
    print $parse_script_fh "TACE=/gsc/scripts/bin/tace\n";
    print $parse_script_fh "ACEDB=`pwd`\n\n";
    print $parse_script_fh "export ACEDB\n\n";
    print $parse_script_fh "echo \$acedb\n\n";
    print $parse_script_fh "\$TACE << EOF\n\n";

    my $shortph5file;
    foreach my $acefile (@acefiles)
    {
        next if $acefile =~ /^\.\.?$/;
        next if $acefile =~ /\.gff$|\.fasta$|\.genomic_canonical\.ace$/;
        next if $acefile =~ /\.txt$/;
        next if $acefile =~ m/dead_genes_list/;

        if ( $acefile =~ /$phase5file/ )
        {
            $shortph5file = $phase5file = $acefile;
        }

        next if $acefile =~ /_phase_[0-5]_ssid_/;
        print $parse_script_fh "$parse  $acedb_data_path/$acefile\n";
    }
    print $parse_script_fh "\nsave\n";
    print $parse_script_fh "quit\n\n";
    print $parse_script_fh "EOF\n\n";
#    print $parse_script_fh
#        "echo \"Parsing of HGMI_$locus_tag $pipe_version files, complete.\" | mailx -s \"HGMI_$locus_tag $pipe_version\" ssurulir\n";

    $parse_script_fh->close();

    my $mode = 0775;
    chmod $mode, $parse_script_full;
    my $aceparce_stdout
        = $acedb_data_path . "/STDOUT_" . $locus_tag . "_ace_parse.txt";

    my @aceparcecmd = ( $parse_script_full, );

    my $aceparse_stderr;
    IPC::Run::run( \@aceparcecmd, '>', \$aceparce_stdout, '2>',
        \$aceparse_stderr, )
        or die "problem: $aceparse_stderr";

    ########################################################
    # gather stats from phase5 ace file and acedb
    ########################################################
    $phase5file = qq{$acedb_data_path/$phase5file};

    unless ( ( -e $phase5file ) and ( !-z $phase5file ) )
    {
        croak
            qq{\n\n NO file,$phase5file,(phase5file) found  or else empty ... from BerRunFinish.pm: $OS_ERROR\n\n };
    }

    my $phase5ace_fh = IO::File->new();
    $phase5ace_fh->open("< $phase5file")
        or die
        "Can't open '$phase5file', phase5file for reading ...from BerRunFinish.pm: $OS_ERROR\n\n";

    my @phase5acecount = ();
    while (<$phase5ace_fh>)
    {
        chomp $ARG;
        if ( $ARG =~ /Subsequence/ )
        {
            push( @phase5acecount, $ARG );
        }
    }
    $phase5ace_fh->close();

    my $acefilecount = scalar(@phase5acecount);
    $db = Ace->connect(
        -path    => "$acedb_maindir_path",
        -program => "$program",
        )
        or die
        "ERROR: cannot connect to acedb ...from BerRunFinish.pm: $OS_ERROR\n \n";

    my @trna_all_d
        = $db->fetch( -query => "Find Sequence $locus_tag\_C*.t* & ! Dead" );
    my @rfam_all_d = $db->fetch(
        -query => "Find Sequence $locus_tag\_C*.rfam* & ! Dead" );
    my @rnammer_all_d = $db->fetch(
        -query => "Find Sequence $locus_tag\_C*.rnammer* & ! Dead" );
    my @orfs_d = $db->fetch(
        -query => "Find Sequence $locus_tag\_C*.p5_hybrid* & ! Dead" );

    my $Totals_not_dead      = 0;
    my $Totals_not_dead_rna  = 0;
    my $Totals_not_dead_orfs = 0;

    $Totals_not_dead
        = scalar(@rfam_all_d) + scalar(@rnammer_all_d) + scalar(@trna_all_d) +
        scalar(@orfs_d);
    $Totals_not_dead_rna
        = scalar(@rfam_all_d) + scalar(@rnammer_all_d) + scalar(@trna_all_d);
    $Totals_not_dead_orfs = scalar(@orfs_d);

    my @trna_all = $db->fetch( -query => "Find Sequence $locus_tag\_C*.t*" );
    my @rfam_all
        = $db->fetch( -query => "Find Sequence $locus_tag\_C*.rfam*" );
    my @rnammer_all
        = $db->fetch( -query => "Find Sequence $locus_tag\_C*.rnammer*" );
    my @orfs
        = $db->fetch( -query => "Find Sequence $locus_tag\_C*.p5_hybrid*" );

    my $Totals_with_dead      = 0;
    my $Totals_with_dead_rna  = 0;
    my $Totals_with_dead_orfs = 0;

    $Totals_with_dead
        = scalar(@rfam_all) + scalar(@rnammer_all) + scalar(@trna_all) +
        scalar(@orfs);
    $Totals_with_dead_rna
        = scalar(@rfam_all) + scalar(@rnammer_all) + scalar(@trna_all);
    $Totals_with_dead_orfs = scalar(@orfs);

	my @p5_blastx = $db->fetch( -query => "Find Sequence $locus_tag\_C*.blastx.p5*" );
	my $p5_blastx_genes = scalar(@p5_blastx);
	$mgap_genome{'p5_blastx_count'} = $p5_blastx_genes;

	my $rnammer_all = scalar(@rnammer_all);
	$mgap_genome{'rnammer_count'} = $rnammer_all;

	my @p5_ace_objects = $db->fetch(
			Sequence => "$locus_tag\_C*p5*"
			);

	my $p5_gene_count = scalar(@p5_ace_objects);
	$mgap_genome{'p5_gene_count'} = $p5_gene_count;

	$_ = 0 for my ($keggscan_count, $iprscan_count, $psortb_count, $ber_naming_count, $dead_gene_count, $manual_review );
	foreach my $gene (@p5_ace_objects) {
		$dead_gene_count++ if $gene->Dead();
		$keggscan_count++ if $gene->KEGG();
		$ber_naming_count++ if $gene->BER_product();
		$iprscan_count++ if $gene->Interpro();
		$psortb_count++ if $gene->PSORT_B();
		$manual_review++ if $gene->ManualReview();
	}

	$mgap_genome{'dead_gene_count'} = $dead_gene_count;
	$mgap_genome{'keggscan_count'} = $keggscan_count;
	$mgap_genome{'iprscan_count'} = $iprscan_count;
	$mgap_genome{'psortb_count'} = $psortb_count;
	$mgap_genome{'ber_naming_count'} = $ber_naming_count;
	$mgap_genome{'mr_count'} = $manual_review;

	print "\n\n" . $locus_tag . "\n\n";
	print $acefilecount
		. "\tSubsequence counts from acefile $shortph5file\n\n";
    print $Totals_not_dead
        . "\tp5_hybrid counts from ACEDB  orfs plus RNA's that are NOT dead genes\n";
    print $Totals_not_dead_rna
        . "\tp5_hybrid counts from ACEDB for ALL RNA's that are NOT dead genes\n";
    print $Totals_not_dead_orfs
        . "\tp5_hybrid counts from ACEDB orfs minus RNA's that are NOT dead genes\n\n";
    print $Totals_with_dead
        . "\tp5_hybrid counts from ACEDB orfs plus RNA's with dead genes (should match acefile $shortph5file)\n";
    print $Totals_with_dead_rna
        . "\tp5_hybrid counts from ACEDB for ALL RNA's with dead genes\n";
    print $Totals_with_dead_orfs
        . "\tp5_hybrid counts from ACEDB orfs with dead genes\n\n";

    if ( $acefilecount == $Totals_with_dead )
    {

        print
            "p5_hybrid ace file counts match p5_hybrid counts in ACEDB... Good :) \n\n";

    }
    else
    {
        print $acefilecount, " ", $Totals_with_dead, "\n";
        print
            "HOUSTON, WE HAVE A PROBLEM, p5_hybrid ace file counts DO NOT MATCH p5_hybrid counts in ACEDB (Totals_with_dead)... BAD :(\n\n";

    }

    print qq{\n\nOther information from AceDB:\n};
    print qq{-----------------------------\n};

	print qq{blastx p5_hybrid gene count:\t$p5_blastx_genes\n};
	print qq{genes with rnammer hits:\t $rnammer_all\n\n};

#
    ## We will dump gff for this genome
    unless ( $cwd eq $acedb_data_path )
    {
        chdir($acedb_data_path)
            or die
            "Failed to change to '$acedb_data_path'...  from BerRunFinish.pm: $OS_ERROR\n\n";
    }

    my $gff_dump_file = $locus_tag ."_phase_5_ssid_". $ssid. ".gff";
	my $bap_dump_gff_cmd = "bap_dump_gene_predictions_gff --sequence-set-id $ssid > $gff_dump_file";
#print "Dumping gff dump (cmd): ". $bap_dump_gff_cmd."\n";
	system("$bap_dump_gff_cmd") == 0
			or die "system $bap_dump_gff_cmd failed: $?";

    ## List dead genes
    my $dead_genes_file = $locus_tag."_dead_genes_list";
    my $dead_genes_cmd = " bap_list_dead_genes --sequence-set-id $ssid > $dead_genes_file";
#    print "Running bap_list_dead_genes (cmd): ". $dead_genes_cmd. "\n";
    system("$dead_genes_cmd") == 0
			or die "system $bap_dump_gff_cmd failed: $?";

    ## Count the dead genes - depends on number of lines
     my $dead_genes = `wc -l < $dead_genes_file` ||
            die "wc failed: $?\n";
    chomp($dead_genes);

    ## BER data location
    my $ber_dir = $locustagdir. "/ber";

	## We will find the *btab files with 0 size
	my $zero_size_btab_file_count = 0;
	my @btab_files = <$ber_dir/*btab>;
	foreach my $btab_file (@btab_files) {
		$zero_size_btab_file_count++ if (stat($btab_file)->size == 0);
	}
	$mgap_genome { 'btab_0_size' } = $zero_size_btab_file_count;

    my $location = $amgap_path . "/"
        . $org_dirname . "/"
        . $assembly_name . "/"
        . $assembly_version;

    my $fasta_file = $location. "/"
					. "Sequence/Unmasked/"
					. $locus_tag . ".v1.contigs.newname.fasta";

    ## We will find the sequence length
	my $total_length = 0;
	if (-e $fasta_file) {
		my $seqio = new Bio::SeqIO('-file' => $fasta_file,
				'-format' => 'fasta', );
		while(my $seq = $seqio->next_seq() ) {
			$total_length += $seq->length;
		}
	}
	$mgap_genome{'sequence_length'} = $total_length;

    ########################################################
    # Writing the rt file
    ########################################################

    my $rtfilename = $project_type
        . "_rt_let_"
        . $locus_tag . "_"
        . $pipe_version . ".txt";
    my $rtfileloc
        = $amgap_path . "/Acedb/Scripts/" . $project_type . "_files";
    my $rtfullname = qq{$rtfileloc/$rtfilename};
    my $rtfile_fh  = IO::File->new();
    $rtfile_fh->open("> $rtfullname")
        or die
        "Can't open '$rtfullname', rtfullname for writing ...from BerRunFinish.pm: $OS_ERROR\n\n";

    print $rtfile_fh
        qq{\n$assembly_name, $locus_tag, a $project_type project has finished processing in AMGAP, BER product naming and now ready to be processed for submissions\n\n};

    print $rtfile_fh qq{\nA copy of BER naming file $sqlitedatafilename has been placed in $anno_submission$sqlitedatafilename\n\n};

    print $rtfile_fh qq{\nBlast data for this genome can be found at $ber_dir\n\n};

    my $sequence_set     = BAP::DB::SequenceSet->retrieve($ssid);
    my $software_version = $sequence_set->software_version();
    my $data_version     = $sequence_set->data_version();

    print $rtfile_fh qq{BAP/MGAP Version: $software_version }, "\n";
    print $rtfile_fh qq{Data Version: $data_version},          "\n\n";
    print $rtfile_fh qq{Location:\n\n};

    print $rtfile_fh qq{$location\n\n};
    print $rtfile_fh
        qq{Gene prediction by the following programs has been run via bap_predict_genes:\n\n};
    print $rtfile_fh qq{Glimmer3\n};
    print $rtfile_fh qq{GeneMark\n};
    print $rtfile_fh qq{trnascan\n};
    print $rtfile_fh qq{RNAmmer\n};
    print $rtfile_fh qq{Rfam v8.1, with Rfam_product\n\n};
    print $rtfile_fh
        qq{bap_merge_genes has been run and includes blastx through phase_5\n\n};
    print $rtfile_fh qq{Here are the gene counts from Oracle AMGAP:\n\n};

    my @sequences        = $sequence_set->sequences();
    my $blastx_counter   = 0;
    my $glimmer2_counter = 0;
    my $glimmer3_counter = 0;
    my $genemark_counter = 0;

    foreach my $i ( 0 .. $#sequences )
    {
        my $sequence     = $sequences[$i];
        my @coding_genes = $sequence->coding_genes();
        foreach my $ii ( 0 .. $#coding_genes )
        {
            my $coding_gene = $coding_genes[$ii];
            if ( $coding_gene->source() =~ 'blastx' )
            {
                $blastx_counter++;
            }
            elsif ( $coding_gene->source() =~ 'glimmer3' )
            {
                $glimmer3_counter++;
            }
            else
            {
                $genemark_counter++;
            }
        }
    }
	
	$mgap_genome {'genemark_gene_count'} = $genemark_counter;
	$mgap_genome {'glimmer_gene_count'} = $glimmer3_counter;
	$mgap_genome {'blastx_gene_count'} = $blastx_counter;

    print $rtfile_fh qq{blastx count   =\t $blastx_counter},   "\n";
    print $rtfile_fh qq{GeneMark count =\t $genemark_counter}, "\n";
    print $rtfile_fh qq{Glimmer3 count =\t $glimmer3_counter}, "\n\n";
    print $rtfile_fh
        qq{Protein analysis by the following programs has been run via PAP workflow:\n\n};
   
    #print $rtfile_fh qq{Interpro v4.8 (database v29.0)\n};
    #print $rtfile_fh qq{Keggscan v56\n};
    #print $rtfile_fh qq{psortB v3.0.3\n};
    #print $rtfile_fh qq{BER v2.5\n}

    print $rtfile_fh qq{Interpro \n};
    print $rtfile_fh qq{Keggscan \n};
    print $rtfile_fh qq{psortB \n};
    print $rtfile_fh qq{BER\n};
    print $rtfile_fh qq{Blastp\n\n};
    print $rtfile_fh
        qq{Location of AMGAP ace files can be located, here:\n\n};
    print $rtfile_fh qq{$acedb_data_path\n\n};

    foreach my $acefile (@acefiles)
    {
        next if $acefile =~ /^\.\.?$/;
        next if $acefile =~ /\.gff$/;
        next if $acefile =~ /\.txt$/;
        print $rtfile_fh "$acefile\n";
    }

    print $rtfile_fh
        qq{\n$locus_tag, QC ace file gene counts verses ACEDB gene counts},
        "\n\n";
    print $rtfile_fh qq{$acefilecount\tgenes from acefile $shortph5file},
        "\n\n";
    print $rtfile_fh
        qq{$Totals_not_dead\tp5_hybrid genes from ACEDB orfs plus RNA\'s that are NOT dead genes},
        "\n";
    print $rtfile_fh
        qq{$Totals_not_dead_rna\tp5_hybrid genes from ACEDB for ALL RNA\'s that are NOT dead genes},
        "\n";
    print $rtfile_fh
        qq{$Totals_not_dead_orfs\tp5_hybrid genes from ACEDB orfs minus RNA\'s that are NOT dead genes minus RNA\'s},
        "\n\n";
    print $rtfile_fh
        qq{$Totals_with_dead\tp5_hybrid genes from ACEDB orfs plus RNA\'s with dead genes (should match acefile $shortph5file ) },
        "\n";
    print $rtfile_fh
        qq{$Totals_with_dead_rna\tp5_hybrid genes from ACEDB for ALL RNA\'s with dead genes},
        "\n";
    print $rtfile_fh
        qq{$Totals_with_dead_orfs\tp5_hybrid genes from ACEDB orfs with dead genes minus RNA\'s},
        "\n\n";

    if ( $acefilecount == $Totals_with_dead )
    {
        print $rtfile_fh
            qq{p5_hybrid ace file counts match p5_hybrid counts in ACEDB... Good :) },
            "\n\n";
    }
    else
    {
        print $rtfile_fh
            qq{HOUSTON, WE HAVE A PROBLEM, p5_hybrid ace file counts DO NOT MATCH p5_hybrid counts in ACEDB (Totals_with_dead)... BAD :\(  },
            "\n\n";
    }

    print $rtfile_fh qq{Other information from AceDB:\n};
    print $rtfile_fh qq{-----------------------------\n};

	print $rtfile_fh qq{blastx p5_hybrid gene count:\t$p5_blastx_genes\n};

	print $rtfile_fh qq{genes with rnammer hits:\t $rnammer_all\n\n};

    print $rtfile_fh qq{GFF dump for thus genome can be downloaded from: $acedb_data_path/$gff_dump_file\n\n};
    
    ## Dead gene stuff
    if ($dead_genes) {
    	print $rtfile_fh qq{Further I found $dead_genes gene tagged as 'Dead'\n};
        open (DG, $dead_genes_file) || die ("Error opening $dead_genes_file :$?\n");
		while (<DG>) {
			chomp;
            print $rtfile_fh qq{$_\n};
		}
	}
    
    print $rtfile_fh qq{\nLocation of this file: $rtfullname\n\n};
    print $rtfile_fh qq{I am transferring ownership to Veena.\n\n};
    print $rtfile_fh qq{Thanks,\n\n};
    print $rtfile_fh qq{ \n};

    send_mail( $ssid, $assembly_name, $rtfileloc, $rtfilename, $rtfullname );
	populate_mgap_genome_db (\%mgap_genome);

    return 1;
}

################################################
################################################

sub version_lookup
{
    my $self           = shift;
    my $v              = shift;
    my $lookup         = undef;
    my %version_lookup = (
        'DEV'  => 'Development',
        'DEV1'  => 'Development_1',
        'DEV2'  => 'Development_2',
        'V1' => 'Version_1.0',
        'V2' => 'Version_2.0',
        'V3' => 'Version_3.0',
        'V4' => 'Version_4.0',
        'V5' => 'Version_5.0',
        'V6' => 'Version_6.0',
        'V7' => 'Version_7.0',
        'V8' => 'Version_8.0',
    );

    if ( exists( $version_lookup{$v} ) )
    {
        $lookup = $version_lookup{$v};
    }
    else
    {
        $v =~ s/V(\d+)/Version_$1.0/;
        $lookup = $v;
    }

    return $lookup;
}

sub send_mail
{

    my ( $ssid, $assembly_name, $rtfileloc, $rtfilename, $rtfullname ) = @ARG;
    my $from = Genome::Utility::Email::construct_address(Genome::Sys->username);

    my $subject
        = "Amgap BER Product Naming script mail for AMGAP SSID: $ssid ($assembly_name)";

    my $body = <<BODY;
The Amgap BER Product Naming script has finished running for MGAP SSID: $ssid ($assembly_name).
The information for the rt ticket has been attached:

File: $rtfilename

Path: $rtfileloc
BODY

    Genome::Utility::Email::send(
        from    => $from,
        to      => [ $from, 'kpepin@genome.wustl.edu'],  #FIXME remove hard-coded addresses
        subject => $subject,
        body    => $body,
        attachments => {
            path        => $rtfullname,
        },
    );
    return 1;
}

sub populate_mgap_genome_db
{
	my (%mgap_genome) = %{$_[0]};
	my $db_file = '/gscmnt/temp212/info/annotation/BAP_db/mgap_genome.sqlite';
	my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file",'','',
			{RaiseError => 1, AutoCommit => 1});

	return 0 unless(defined($dbh));

	## Insert data to genome_info
	my $genome_info_sql = <<SQL;
	INSERT INTO genome_info (locus_name, assembly_name, sequence_length, org_dirname)
	VALUES (?, ?, ?, ?);
SQL

	my $sth = $dbh->prepare($genome_info_sql);
	$sth->bind_param(1, $mgap_genome{'locus_name'});
	$sth->bind_param(2, $mgap_genome{'assembly_name'});
	$sth->bind_param(3, $mgap_genome{'sequence_length'});
	$sth->bind_param(4, $mgap_genome{'org_dirname'});

	$sth->execute() or confess "Couldn't execute statement: " . $sth->errstr ;
	$sth->finish;
	my $genome_id = $dbh->func('last_insert_rowid');

	## Insert data to bacterial_run_details
	my $genome_bacterial_run_details_sql = <<SQL;
	INSERT INTO bacterial_run_details (ssid, user, genome_id, db_version_id, acedb_version, project_type_id, run_type_id, pipeline_version_id)
	VALUES (?,?,?,?,?,?,?,?);
SQL

	my $sth2 = $dbh->prepare($genome_bacterial_run_details_sql);
	$sth2->bind_param(1, $mgap_genome{'ssid'});
	$sth2->bind_param(2, $mgap_genome{'user'});
	$sth2->bind_param(3, $genome_id);
	$sth2->bind_param(4, $mgap_genome{'db_version_id'});
	$sth2->bind_param(5, $mgap_genome{'acedb_ver'});
	$sth2->bind_param(6, $mgap_genome{'project_type'});
	$sth2->bind_param(7, $mgap_genome{'run_type'});
	$sth2->bind_param(8, $mgap_genome{'pipe_version'});

	$sth2->execute() or confess "Couldn't execute statement: " . $sth2->errstr ;
	$sth2->finish;

	## Insert data to hit_counts table
	my $hit_counts_sql = <<SQL;
	INSERT INTO hit_counts (p5_gene_count, p5_blastx_count, rnammer_count, 
							genemark_gene_count, glimmer_gene_count,keggscan_count, 
							iprscan_count, psortb_count, ber_naming_count, genome_id, 
							dead_gene_count, blastx_gene_count, btab_0_size, mr_count)
	VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);
SQL

	my $sth3 = $dbh->prepare($hit_counts_sql);
	$sth3->bind_param(1, $mgap_genome{'p5_gene_count'});
	$sth3->bind_param(2, $mgap_genome{'p5_blastx_count'});
	$sth3->bind_param(3, $mgap_genome{'rnammer_count'});
	$sth3->bind_param(4, $mgap_genome{'genemark_gene_count'});
	$sth3->bind_param(5, $mgap_genome{'glimmer_gene_count'});
	$sth3->bind_param(6, $mgap_genome{'keggscan_count'});
	$sth3->bind_param(7, $mgap_genome{'iprscan_count'});
	$sth3->bind_param(8, $mgap_genome{'psortb_count'});
	$sth3->bind_param(9, $mgap_genome{'ber_naming_count'});
	$sth3->bind_param(10, $genome_id);
	$sth3->bind_param(11, $mgap_genome{'dead_gene_count'});
	$sth3->bind_param(12, $mgap_genome{'blastx_gene_count'});
	$sth3->bind_param(13, $mgap_genome{'btab_0_size'});
	$sth3->bind_param(14, $mgap_genome{'mr_count'});

	$sth3->execute() or confess "Couldn't execute statement: " . $sth3->errstr ;
	$sth3->finish;

	$dbh->disconnect if defined($dbh);

	return 1;
}

1;
