package Genome::Model::GenePrediction::Command::Bacterial::Finish;

use strict;
use warnings;

use Genome;

use BAP::DB::SequenceSet;
use BAP::DB::CodingGene;
use English;
use IO::File;
use Cwd;
require "pwd.pl";
use Getopt::Long;
use Pod::Usage;
use Ace;
use Ace::Sequence;
use Carp;
use File::Path;
use Genome::Utility::Email;
use File::Path 'make_path';

class Genome::Model::GenePrediction::Command::Bacterial::Finish {
    is  => 'Command',
    doc => "BAP finish project",
    has => [
        locus_id => {
            is  => 'String',
            doc => "locus id string",
        },
        org_dirname => {
            is  => 'String',
            doc => "organism directory name",
        },
        project_type => {
            is  => 'String',
            doc => "project type; usually HGMI",
        },
        acedb_version => {
            is  => 'String',
            doc => "AceDB version; specify like V1 or Version_1.0",
        },
        assembly_name => {
            is  => 'String',
            doc => "assembly name",

        },
        assembly_version => {
            is  => 'String',
            doc => "assembly version; specify like V1 or Version_1.0",
        },
        pipe_version => {
            is  => 'String',
            doc => "pipeline version; specify like V1 or Version_1.0",
        },
        path_base => {
            is  => 'String',
            doc => "path base for output file locations",
        },
        sequence_set_id => {
            is  => 'Integer',
            doc => "sequence set id for genome assembly",
        },
    ],
    has_optional => [
        dev => {
            is      => 'Boolean',
            doc     => "use development MGAP database",
            default => 0,
        },
        no_acedb => {
            is      => 'Boolean',
            doc     => "skip acedb",
		    #default => 0,
        },
        no_mail => {
            is => 'Boolean',
            doc => "disable sending results email",
            default => 0,
        },
        _tace_location => {
            is      => 'String',
            doc     => "location of tace/acedb terminal interface",
            default => "/gsc/scripts/bin/tace",
        },

    ],
};


sub sub_command_sort_position {40}

sub help_synopsis
{
return <<"EOS"
genome model gene-prediction finish [options]
EOS
}

sub help_brief
{
    "BAP project finish steps";

}

sub help_detail
{
return <<"EOS"
runs project finishing steps for bacterial prediction pipeline
EOS
}

sub execute
{
    my $self = shift;
    my $ssid = $self->sequence_set_id;

    if ( $self->dev ) { $BAP::DB::DBI::db_env = 'dev'; }

    my $program = $self->_tace_location;
    my $cwd     = getcwd();
    my @bugline = split( '/', $cwd );
    my ( $bug, $version, $hgmi_path, $acedb_scripts_path, $acedb_path,
        $acedb_acefile_path, $acedb_scripts_HGMI_files, $acedb_installation_path, $path_base );

#This section chops the directory struture up for each of the project types: HGMI, HMPP, ENTER and BACTER
    $bug       = $self->assembly_name;
    $version   = $self->pipe_version;
    $hgmi_path = $self->path_base;
    $path_base = $self->path_base;

    my $hgmi_dumps_path
        = $hgmi_path . "/"
        . $self->org_dirname . "/"
        . $self->assembly_name . "/"
        . $self->assembly_version . "/BAP/"
        . $self->pipe_version
        . "/Dumps";
    unless ( -d $hgmi_dumps_path )
    {
        confess
            "The directory '$hgmi_dumps_path', hgmi_dumps_path does not exist, in bap_finish_project: $OS_ERROR";
    }
    unless ( $cwd eq $hgmi_dumps_path )
    {
        chdir($hgmi_dumps_path)
            or confess
            "Failded to change to '$hgmi_dumps_path'...from bap_finish_project: $OS_ERROR";
    }

    # FIXME This is so sloppy there are no words to adequately describe it. Holy shit would be a good start.
    # FIXME For now, there is an ace installation directory (for HGMI, at /gscmnt/278/analysis/HGMI) that has
    # these wspec and databases directories that are required for ace upload to work. I'd like to look into
    # getting rid of those in the near future...
    if ($self->project_type =~ /HGMI/) {
        $acedb_installation_path = "/gscmnt/278/analysis/HGMI/Acedb/";
        my $acedb_version_long = $self->version_lookup($self->acedb_version);
        $acedb_installation_path .= $acedb_version_long;

        $acedb_path         = $hgmi_path . "/Acedb/" . $acedb_version_long;
        $acedb_scripts_path = $hgmi_path . "/Acedb/Scripts";
        $acedb_scripts_HGMI_files = $hgmi_path . "/Acedb/Scripts/HGMI_files";

        unless (-d $acedb_path) {
            confess "The directory '$acedb_path' does not exist!";
        }

        $acedb_acefile_path = $acedb_path . "/ace_files/" . $self->locus_id . "/" . $version;

        unless ( -d $acedb_acefile_path) {
            confess "The directory '$acedb_acefile_path', acedb_acefile_path does not exist, in bap_finish_project: $OS_ERROR";
        }
    }
    elsif ( $self->project_type =~ /HMPP/ )
    {

        unless ( $hgmi_path =~ /HMPP$/ )
        {
            $hgmi_path
                = qq{/gscmnt/254/analysis/bacterial_analysis/HMPP};    #HMPP
        }
        #######################
        #acedb_scripts_path
        #######################
#$acedb_scripts_path = qq{$hgmi_path/$bugline[6]/$bugline[7]/$bugline[8]/Acedb/Scripts}; #HMPP
        $acedb_scripts_path
            = $hgmi_path .. "/"
            . $self->org_dirname . "/"
            . $self->assembly_name . "/"
            . $self->assembly_version
            . "/Acedb/Scripts";

        if ( $self->org_dirname eq "E_coli" )
        {
            $acedb_scripts_path
                = qq{$hgmi_path/E_coli/Escherichia_coli_K12_MG1655/Version_1.0/Acedb/Scripts}
                ;    #HMPP_ECOLI
            $acedb_path
                = qq{$hgmi_path/E_coli/Escherichia_coli_K12_MG1655/Version_1.0/Acedb/Version_1.0}
                ;    #HMPP_ECOLI
        }
        elsif ( $self->org_dirname eq "R_sphaeroides" )
        {
            $acedb_scripts_path
                = qq{$hgmi_path/R_sphaeroides/Rhodobacter_sphaeroides_FNL/Version_1.0/Acedb/Scripts}
                ;    #HMPP_R_sphaeroides
            $acedb_path
                = qq{$hgmi_path/R_sphaeroides/Rhodobacter_sphaeroides_FNL/Version_1.0/Acedb/Version_1.0}
                ;    #HMPP_R_sphaeroides
        }
        elsif ( $self->org_dirname eq "S_aureus" )
        {
            $acedb_scripts_path
                = qq{$hgmi_path/S_aureus/Staphylococcus_aureus_usa300/Version_1.0/Acedb/Scripts}
                ;    #HMPP_S_aureus
            $acedb_path
                = qq{$hgmi_path/S_aureus/Staphylococcus_aureus_usa300/Version_1.0/Acedb/Version_1.0}
                ;    #HMPP_S_aureus
        }
        elsif ( $self->org_dirname eq "T_pallidum" )
        {
            $acedb_scripts_path
                = qq{$hgmi_path/T_pallidum/Tpallidum_SamoanD.10x.454AllContigs_DFT/Version_1.0/Acedb/Scripts};
            $acedb_path
                = qq{$hgmi_path/T_pallidum/Tpallidum_SamoanD.10x.454AllContigs_DFT/Version_1.0/Acedb/Version_1.0};
        }
        else
        {

#$acedb_scripts_path = qq{$hgmi_path/$bugline[6]/$bugline[7]/$bugline[8]/Acedb/Scripts}; #HMPP
            $acedb_scripts_path
                = $hgmi_path . "/"
                . $self->org_dirname . "/"
                . $self->assembly_name . "/"
                . $self->assembly_version
                . "/Acedb/Scripts";

#$acedb_path         = qq{$hgmi_path/$bugline[6]/$bugline[7]/$bugline[8]/Acedb/Version_1.0};     #HMPP
            $acedb_path
                = $hgmi_path . "/"
                . $self->org_dirname . "/"
                . $self->assembly_name . "/"
                . $self->assembly_version
                . "/Acedb/"
                . $self->pipe_version;
        }

        $acedb_acefile_path
            = $acedb_path."/ace_files/".$self->locus_id."/$version";    #HMPP
        $acedb_scripts_HGMI_files = qq{$acedb_scripts_path/HMPP_files};  #HMPP

    }
    elsif ( $self->project_type =~ /ENTER/ )
    {
        $hgmi_path
            = qq{/gscmnt/254/analysis/bacterial_analysis/Enterobacter}; #ENTER
        $acedb_scripts_path = qq{$hgmi_path/Acedb/Scripts};             #ENTER
        $acedb_path         = qq{$hgmi_path/Acedb/Version_1.0};         #ENTER
        $acedb_acefile_path
            = qq{$acedb_path/ace_files/}.$self->locus_id.qq{/$version};             #ENTER
        $acedb_scripts_HGMI_files
            = qq{$acedb_scripts_path/ENTER_files};                      #ENTER
    }
    elsif ( $self->project_type =~ /BACTER/ )
    {
        $version = $self->pipe_version;

#$hgmi_path = qq{/gscmnt/277/analysis/bacterial_analysis/Bacteroides};                   #BACTER
        $hgmi_path = qq{/gscmnt/temp212/info/annotation/HGMI_tmp/Bacteroides/}
            ;    #BACTER
        $acedb_scripts_path = qq{$hgmi_path/Acedb/Scripts};        #BACTER
        $acedb_path         = qq{$hgmi_path/Acedb/Version_1.0};    #BACTER
          #$acedb_acefile_path = qq{$acedb_path/ace_files/$locus_id/$version};                     #BACTER
        $acedb_acefile_path
            = "/gscmnt/temp212/info/annotation/HGMI_tmp/Bacteroides/Acedb/Version_1.0/ace_files/" . $self->locus_id . "/$version"
            ;    #BACTER
        $acedb_scripts_HGMI_files
            = qq{$acedb_scripts_path/BACTER_files};    #BACTER
    }
    elsif ( $self->project_type =~ /CORE/ )
    {

        #FIXME: maybe we should check these directories, and make them
        # if they don't exist.
##		$acedb_installation_path = "/gscmnt/254/analysis/bacterial_analysis/CORE/Acedb/"; VEE modified
                $acedb_installation_path = qq{$path_base/Acedb/};
		 
        my $acedb_version_long = $self->version_lookup($self->acedb_version);
		$acedb_installation_path .= $acedb_version_long;

##        $hgmi_path = qq{/gscmnt/254/analysis/bacterial_analysis/CORE};   #CORE vee modified
##        $hgmi_path = qq{$path_base/CORE};   ##vee modified
          $hgmi_path = qq{$path_base};   
        $acedb_scripts_path = qq{$path_base/Acedb/Scripts};              #CORE
             #$acedb_path = qq{$hgmi_path/Acedb/Version_1.0}; #CORE
        
        #$acedb_path = qq{$hgmi_path/Acedb/$acedb_version_long};    #CORE vee modified
        $acedb_path = qq{$path_base/Acedb/$acedb_version_long};    #CORE
        $acedb_acefile_path  = qq{$acedb_path/ace_files/}.$self->locus_id."/".$version;        #CORE
  
        $acedb_scripts_HGMI_files = qq{$acedb_scripts_path/CORE_files};  #CORE
        mkpath( [$acedb_acefile_path], 0, 0775 );
        carp "set path for CORE: $acedb_acefile_path";
    }
    else
    {
        confess
            "project-type: ".$self->project_type.", not set correctly!  Please see documentation.\n\n";
    }

    $cwd = getcwd();
    unless ( -d $acedb_acefile_path )
    {
        carp "acefile path $acedb_acefile_path doesn't exist. Creating...";
        mkpath( [$acedb_acefile_path], 0, 0775 );
    }

    unless ( -d $hgmi_dumps_path )
    {
        confess
            "The directory '$hgmi_dumps_path', hgmi_dumps_path does not exist, in bap_finish_project: $OS_ERROR";
    }
    unless ( $cwd eq $hgmi_dumps_path )
    {
        chdir($hgmi_dumps_path)
            or confess
            "Failded to change to '$hgmi_dumps_path'...from bap_finish_project: $OS_ERROR";
    }
    my @phase_number = ( 0 .. 5 );

    foreach my $phase (@phase_number)
    {

        my $dump_output
            = join( "_", $self->locus_id, 'phase', $phase, 'ssid', $ssid );
        $dump_output = join( ".", $dump_output, 'ace' );

        # Blow up the files we're going to write $outfile to, if they exist

        foreach my $old_dump_output ($dump_output)
        {

            if ( defined($old_dump_output) && ( -e $old_dump_output ) )
            {
                unlink($old_dump_output)
                    || warn "Failed to unlink '$old_dump_output'!";
            }
        }
        my ( $dump_cmd, $dump_dumper );

        $self->status_message("Running ace dump from " . ($self->dev ? "dev" : "prod"));
        $self->status_message("sequence set id, $ssid : phase, $phase : ace file, $dump_output");
        my $rv = Genome::Model::Tools::Bacterial::AceDumpGenes->execute(
            sequence_set_id => $ssid,
            phase => $phase,
            ace_file => $dump_output,
            dev => $self->dev,
        );
        unless ($rv) {
            $self->error_message("Ace dumping failed with sequence set id $ssid");
            confess;
        }

        my $newACEfilelink = qq{$acedb_acefile_path/$dump_output};
        unless ( -s $newACEfilelink )
        {
            carp "current directory: $cwd";

            symlink "$cwd/$dump_output", "$acedb_acefile_path/$dump_output"
                or warn
                "File exists or Can not make symlink for $dump_output : $OS_ERROR\n";
        }
    }
######################################
    # HGMI Ace_Parse_Script and RT_Writer
######################################

    if ( $cwd ne $acedb_scripts_path )
    {
        $self->status_message("Changing directory to $acedb_scripts_path");
        make_path($acedb_scripts_path) unless -d $acedb_scripts_path;
        chdir($acedb_scripts_path);
    }

    opendir( ACEDIR, "$acedb_acefile_path" )
        or confess "Can't open acefile path: $!\n";

    my @acefiles = ();

    @acefiles = readdir(ACEDIR);

    closedir(ACEDIR);

    my $parsefile_name
        = "parsefiles_wens_" . $self->locus_id . "_" . $version . ".sh";

    open( PARSEFILE, "> $parsefile_name" )
        or confess "Can not open new parse file: $!\n";

    my $parse = "parse";

    print PARSEFILE "#!/bin/sh -x\n\n";
    print PARSEFILE
        "#if you call script from bash, tace will follow links!\n\n";
    print PARSEFILE "TACE=/gsc/scripts/bin/tace\n";
    #print PARSEFILE "ACEDB=" . $acedb_scripts_path . "\n\n";
    print PARSEFILE "ACEDB=" . $acedb_path . "\n\n";
    print PARSEFILE "export ACEDB\n\n";
    print PARSEFILE "echo \$acedb\n\n";
    print PARSEFILE "\$TACE << EOF\n\n";

    foreach my $acefile (@acefiles)
    {
        next if $acefile =~ /^\.\.?$/;
        print PARSEFILE "$parse $acedb_acefile_path/$acefile\n";
    }
    print PARSEFILE "\nsave\n";
    print PARSEFILE "quit\n\n";
    print PARSEFILE "EOF\n\n";

#print PARSEFILE "echo \"Parsing of HGMI_$locus_id /$version files, complete.\" | mailx -s \"HGMI_$locus_id /$version\" wnash , pthiru\n";
    #print PARSEFILE
    #    "echo \"Parsing of HGMI_".$self->locus_id." /$version files, complete.\" | mailx -s \"HGMI_".$self->locus_id." /$version\" ssurulir\n";

    close(PARSEFILE);

    my $mode = 0775;
    chmod $mode, $parsefile_name;

    # read ace files into ACEDB
    my $cwd3 = getcwd();

    if ( $cwd3 ne $acedb_path )
    {
        $self->status_message("Changing directory to $acedb_path");
        chdir($acedb_path);
    }

    my $parse_ACEDB = qq{$acedb_scripts_path/$parsefile_name};

    my @trna_all_d    = ();
    my @rfam_all_d    = ();
    my @rnammer_all_d = ();
    my @orfs_d        = ();
    my @trna_all      = ();
    my @rfam_all      = ();
    my @rnammer_all   = ();
    my @orfs          = ();

    my ( $ace_file, $Totals_not_dead, $Totals_not_dead_rna,
        $Totals_not_dead_orfs );
    my ( $ace_file2, $Totals_with_dead, $Totals_with_dead_rna,
        $Totals_with_dead_orfs );
    my $acefilecount = 0;
    my $phase5       = qq{_phase_5_ssid};
    my $file         = $self->locus_id . $phase5;

    unless ($self->no_acedb)
    {
        my $locus_id = $self->locus_id;

        # FIXME: don't use backtics
        #my $pACEDBcmd = `$parse_ACEDB`;
        IPC::Run::run([$parse_ACEDB]);

        # get p5_hybrid and acedb counts
        opendir( ACEDIR2, $acedb_acefile_path )
            or confess "Can't open $acedb_acefile_path path: $OS_ERROR\n";

        my @acefileseek = ();

        @acefileseek = readdir(ACEDIR2);

        closedir(ACEDIR2);

        foreach $ace_file (@acefileseek)
        {

            if ( $ace_file =~ /^\.\.?$/ )
            {
                next;
            }

            if ( $ace_file =~ /$file/ )
            {

                $ace_file2 = $ace_file;

                my $acefile_search = "$acedb_acefile_path/$ace_file";

                if ( -e $acefile_search )
                {

                    open( ACEFILE, $acefile_search )
                        or confess "Can't open $acefile_search: $OS_ERROR\n";
                    my @acecount = ();
                    while ( my $aceline = <ACEFILE> )
                    {
                        chomp $aceline;

                        if ( $aceline =~ /^Subsequence/ )
                        {

                            push( @acecount, $aceline );
                        }
                    }

                    $acefilecount = scalar(@acecount);

                    #connecting to acedb database

                    my $db = Ace->connect(
                        -path    => "$acedb_installation_path",
                        -program => "$program"
                    ) or confess "ERROR: cannot connect to acedb\n";

                    #mining data from acedb
                    # FIXME: metarna goes here.
                    @trna_all_d = $db->fetch(
                        -query => "Find Sequence $locus_id\_C*.t* & ! Dead" );
                    @rfam_all_d = $db->fetch( -query =>
                            "Find Sequence $locus_id\_C*.rfam* & ! Dead" );
                    @rnammer_all_d = $db->fetch( -query =>
                            "Find Sequence $locus_id\_C*.rnammer* & ! Dead" );
                    @orfs_d
                        = $db->fetch( -query =>
                            "Find Sequence $locus_id\_C*.p5_hybrid* & ! Dead"
                        );

                    $Totals_not_dead      = 0;
                    $Totals_not_dead_rna  = 0;
                    $Totals_not_dead_orfs = 0;

                    $Totals_not_dead
                        = scalar(@rfam_all_d) 
                        + scalar(@rnammer_all_d)
                        + scalar(@trna_all_d)
                        + scalar(@orfs_d);
                    $Totals_not_dead_rna
                        = scalar(@rfam_all_d) 
                        + scalar(@rnammer_all_d)
                        + scalar(@trna_all_d);
                    $Totals_not_dead_orfs = scalar(@orfs_d);

                    @trna_all = $db->fetch(
                        -query => "Find Sequence $locus_id\_C*.t*" );
                    @rfam_all = $db->fetch(
                        -query => "Find Sequence $locus_id\_C*.rfam*" );
                    @rnammer_all = $db->fetch(
                        -query => "Find Sequence $locus_id\_C*.rnammer*" );
                    @orfs = $db->fetch(
                        -query => "Find Sequence $locus_id\_C*.p5_hybrid*" );

                    $Totals_with_dead      = 0;
                    $Totals_with_dead_rna  = 0;
                    $Totals_with_dead_orfs = 0;

                    $Totals_with_dead
                        = scalar(@rfam_all) 
                        + scalar(@rnammer_all)
                        + scalar(@trna_all)
                        + scalar(@orfs);
                    $Totals_with_dead_rna
                        = scalar(@rfam_all) 
                        + scalar(@rnammer_all)
                        + scalar(@trna_all);
                    $Totals_with_dead_orfs = scalar(@orfs);

                    print "\n\n" . $self->locus_id . "\n\n";
                    print $acefilecount
                        . "\tSubsequence counts from acefile $ace_file\n\n";
                    print $Totals_not_dead
                        . "\tp5_hybrid counts from ACEDB  orfs plus RNA's that are NOT dead genes\n";
                    print $Totals_not_dead_rna
                        . "\tp5_hybrid counts from ACEDB for ALL RNA's that are NOT dead genes\n";
                    print $Totals_not_dead_orfs
                        . "\tp5_hybrid counts from ACEDB orfs minus RNA's that are NOT dead genes\n\n";
                    print $Totals_with_dead
                        . "\tp5_hybrid counts from ACEDB orfs plus RNA's with dead genes (should match acefile $ace_file)\n";
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

                        print
                            "HOUSTON, WE HAVE A PROBLEM, p5_hybrid ace file counts DO NOT MATCH p5_hybrid counts in ACEDB... BAD :(\n\n";

                    }
                }
            }
        }
    }

    # Writing the rt file

    my $cwd2 = getcwd();

    if ( $cwd2 ne $acedb_scripts_HGMI_files )
    {
        $self->status_message("Changing directory to $acedb_scripts_HGMI_files");
        chdir($acedb_scripts_HGMI_files);
    }
    my $rtfile_name
        = $self->project_type . "_rt_let_" . $self->locus_id . "_" . $version . ".txt";

    open( RTFILE, "> $rtfile_name" )
        or confess "Can not open new RT file: $OS_ERROR\n";
    print RTFILE
        $bug,", ",$self->locus_id," , ",$self->project_type," project has finished in MGAP.\n\n";

    my $sequence_set     = BAP::DB::SequenceSet->retrieve($ssid);
    my $software_version = $sequence_set->software_version();
    my $data_version     = $sequence_set->data_version();

    my $svn_version = '$Revision$';
    $svn_version =~ s/\D+//g;

    print RTFILE qq{BAP/MGAP Version: $software_version ($svn_version)}, "\n";
    print RTFILE qq{Data Version: $data_version}, "\n\n";
    print RTFILE "Location:\n\n";

    my $location;

    if ( $self->project_type =~ /BACTER/ )
    {
        $location = join( "/", @bugline[ 0 .. 9 ] );
    }
    else
    {
        $location = join( "/", @bugline[ 0 .. 7 ] );
        carp "location is $location";
    }
    print RTFILE "$location\n\n";
    print RTFILE
        "The following analysis have been run via bap_predict_genes:\n\n";
    print RTFILE "Glimmer2\n";
    print RTFILE "Glimmer3\n";
    print RTFILE "GeneMark\n";
    print RTFILE "trnascan\n";
    print RTFILE "RNAmmer\n";
    print RTFILE "Meta RNA\n";
    print RTFILE "Rfam v8.1, with Rfam_product\n\n";
    print RTFILE
        "bap_merge_genes has been run and includes blastx through phase_5\n\n";
    print RTFILE "Here are the gene counts from MGAP:\n\n";

    my @sequences = $sequence_set->sequences();

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
            elsif ( $coding_gene->source() =~ 'glimmer2' )
            {
                $glimmer2_counter++;
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

    print RTFILE qq{blastx count   =\t $blastx_counter},   "\n";
    print RTFILE qq{GeneMark count =\t $genemark_counter}, "\n";
    print RTFILE qq{Glimmer2 count =\t $glimmer2_counter}, "\n";
    print RTFILE qq{Glimmer3 count =\t $glimmer3_counter}, "\n\n";
    print RTFILE "Location of MGAP ace files can be located, here:\n\n";
    print RTFILE "$acedb_acefile_path\n\n";
    print RTFILE
        "Here is a list of ace files in the directory.  All of the files have been\n";
    print RTFILE "parsed into ACEdb:\n\n";

    foreach my $acefile (@acefiles)
    {
        next if $acefile =~ /^\.\.?$/;
        print RTFILE "$acefile\n";
    }
    #if ($skip_acedb_flag)
    if ($self->no_acedb)
    {
        $ace_file2             = $file;
        $acefilecount          = 0;
        $Totals_not_dead       = 0;
        $Totals_not_dead_rna   = 0;
        $Totals_not_dead_orfs  = 0;
        $Totals_with_dead      = 0;
        $Totals_with_dead_rna  = 0;
        $Totals_with_dead_orfs = 0;
    }

    print RTFILE
        "\n",$self->locus_id,", QC ace file gene counts verses ACEDB gene counts",
        "\n\n";
    print RTFILE qq{$acefilecount \tgenes from acefile $ace_file2\n\n};
    print RTFILE
        qq{$Totals_not_dead\ tp5_hybrid genes from ACEDB orfs plus RNA\'s that are NOT dead genes},
        "\n";
    print RTFILE
        qq{$Totals_not_dead_rna \tp5_hybrid genes from ACEDB for ALL RNA\'s that are NOT dead genes},
        "\n";
    print RTFILE
        qq{$Totals_not_dead_orfs \tp5_hybrid genes from ACEDB orfs minus RNA\'s that are NOT dead genes minus RNA\'s},
        "\n\n";
    print RTFILE
        qq{$Totals_with_dead \tp5_hybrid genes from ACEDB orfs plus RNA\'s with dead genes (should match acefile $ace_file2 ) },
        "\n";
    print RTFILE
        qq{$Totals_with_dead_rna \tp5_hybrid genes from ACEDB for ALL RNA\'s with dead genes\n};
    print RTFILE
        qq{$Totals_with_dead_orfs \tp5_hybrid genes from ACEDB orfs with dead genes minus RNA\'s},
        "\n\n";

    if ( $acefilecount == $Totals_with_dead )
    {
        print RTFILE
            qq{p5_hybrid ace file counts match p5_hybrid counts in ACEDB... Good :) },
            "\n\n";
    }
    else
    {
        print RTFILE
            qq{HOUSTON, WE HAVE A PROBLEM, p5_hybrid ace file counts DO NOT MATCH p5_hybrid counts in ACEDB... BAD :\(  },
            "\n\n";
    }
    print RTFILE "\nI am transferring ownership to Veena.\n\n";
    print RTFILE "Thanks,\n\n";
    print RTFILE "Bill\n";

    close(RTFILE);

    my $user              = Genome::Sys->username;
    my $sequence_set_name = $sequence_set->sequence_set_name();

    unless ($self->no_mail) {
        $self->send_mail( $ssid, $sequence_set_name, $user, );
    }

    return 1;
}

sub version_lookup
{
    my $self           = shift;
    my $v              = shift;
    my $lookup         = undef;
    my %version_lookup = (
		'DEV'  => 'Development',
		'DEV1'  => 'Development_1',
		'DEV2'  => 'Development_2',
        'V1'  => 'Version_1.0',
        'V2'  => 'Version_2.0',
        'V3'  => 'Version_3.0',
        'V4'  => 'Version_4.0',
        'V5'  => 'Version_5.0',
        'V6'  => 'Version_6.0',
        'V7'  => 'Version_7.0',
        'V8'  => 'Version_8.0',
        'V9'  => 'Version_9.0',
        'V10' => 'Version_10.0',
    );

    if ( exists( $version_lookup{$v} ) )
    {
        $lookup = $version_lookup{$v};
    }
    else
    {
        $lookup = $v;
        $lookup =~ s/V(\d+)/Version_$1.0/;
    }

    return $lookup;
}

sub send_mail
{
    my $self = shift;
    my ( $ss_id, $ss_name, $user ) = @_;

    my $from = Genome::Utility::Email::construct_address($user);

    my $subject = "BAP Finish script mail for MGAP SSID: $ss_id ($ss_name)";

    my $body = <<BODY;
The bap_finish_scripts script has finished running MGAP SSID: $ss_id ($ss_name).
BODY

    Genome::Utility::Email::send(
        from => $from,
        to => [ $from, 'kpepin@watson.wustl.edu' ], # FIXME remove hard-coded addresses
        subject => $subject,
        body => $body,
    );
}

1;
