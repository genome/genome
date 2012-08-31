package Genome::Model::Tools::Hgmi::DirBuilder;

use strict;
use warnings;

use Genome;
use Command;
use Carp;
use English;
use File::Slurp;
use File::Path qw(make_path);
use DateTime;
use List::MoreUtils qw/ uniq /;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [
		   'path'                  => {
					       is  => 'String',
					       doc => "direcotry path"
					      },
		   'org_dirname'           => {
					       is  => 'String',
					       doc => "organism abbreviated name"
					      },
		   'assembly_version_name' => {
					       is  => 'String',
					       doc => "complete assembly name and version"
					      },
		   'assembly_version'      => {
					       is  => 'String',
					       doc => "analysis version"
					      },
		   'pipe_version'         => {
					      is  => 'String',
					      doc => "pipeline version"
					     },
		   'cell_type'            => {
					      is  => 'String',
					      doc => "[BACTERIA|EUKARYOTES|VIRAL]"
					     },
		  ],
			 
			);

sub help_brief
{
    "tool for creating the directory structure for HGMI/Annotation projects";
}

sub help_synopsis
{
    my $self = shift;
    return <<"EOS"
Creates the standard annotation/analysis directory structure.
EOS

}

sub help_detail
{
    my $self = shift;
    return <<"EOS"
This tool creates the standard directory structure for a project.
EOS

}

sub execute
{
    my $self = shift;

    my $date;
    my $dt  = DateTime->now;
    my $ymd = $dt->ymd('');
    $date = $ymd;
    my @moredirs = ();

    if ( ($self->cell_type =~ /BACTERIA/) || 
         ($self->cell_type =~ /ARCHEA/)  ||
         ($self->cell_type =~ /VIRAL/) )
    {
        @moredirs = (
            'Ensembl_pipeline', 'Gene_predictions',
            'Gene_merging',     'Genbank_submission',
            'Sequence',         'Rfam',
            'Repeats',          'Kegg',
            'Cog',              'Interpro',
            'psortB',           'Blastp',
            'BAP'
        );    #Bacteria
    } 
    else {
	die "Can't run dir-builder... check cell_type in config file ... from DirBuilder.pm: $OS_ERROR\n\n";
    }

    # it looks like there's a bunch of extra
    # filepath variables that are built here
    my $dirpatha = undef;
    my $dirpath  = undef;
    $dirpatha = $self->path . "/"
        . $self->org_dirname . "/"
        . $self->assembly_version_name;
    $dirpath = $dirpatha . "/" . $self->assembly_version;
    my $newdir       = $self->path . "/" . $self->org_dirname;
    my $assembly_dir = $self->assembly_version_name;
    my $version      = $self->pipe_version;

    foreach my $file (@moredirs)
    {
        my @cmd = ();

        $self->add_directory(\@cmd,$newdir);
        $self->add_directory(\@cmd,$dirpatha);
        $self->add_directory(\@cmd,$dirpath);
        my $dirpathfile = $dirpath . "/" . $file;
        my $dirpfversion = $dirpath . "/" . $file . "/" . $version;
        my $dirpfScripts = $dirpath . "/" . $file . "/Scripts";
        $self->add_directory(\@cmd,$dirpathfile);
        $self->add_directory(\@cmd,$dirpfversion);
        $self->add_directory(\@cmd,$dirpfScripts);

	my (
	    $dirpfvGlimmer2,$dirpfvGlimmer3,
	    $dirpfvHybrid,$dirpfvGenemark,
	    $dirpfvDumps,$dirpfvSequence,
	    );

        if ( $file =~ 'Ensembl_pipeline' )
        {

            my $dirpfvSequence
                = $dirpath . "/" . $file . "/" . $version . "/Sequence";
            my $dirpfvDumps
                = $dirpath . "/" . $file . "/" . $version . "/Dumps";

            $self->add_directory(\@cmd,$dirpfvSequence);
            $self->add_directory(\@cmd,$dirpfvDumps);

        }
        if ( $file =~ 'Genbank_submission' )
        {
            my $dirpfvDraft
                = qq{$dirpath/$file/$version/Draft_submission_files};
            my $dirpfvAnnotated
                = qq{$dirpath/$file/$version/Annotated_submission};

            $self->add_directory(\@cmd,$dirpfvDraft);
            $self->add_directory(\@cmd,$dirpfvAnnotated);

        }
        if ( $file =~ 'Gene_predictions' )
        {
            my $dirpfvEannot = qq{$dirpath/$file/$version/Eannot};
            my $dirpfvGenewise = qq{$dirpath/$file/$version/Genewise};
            my $dirpfvGenemark = qq{$dirpath/$file/$version/Genemark};
            my $dirpfvGlimmer2 = qq{$dirpath/$file/$version/Glimmer2};
            my $dirpfvGlimmer3 = qq{$dirpath/$file/$version/Glimmer3};
            my $dirpfvSnap = qq{$dirpath/$file/$version/Snap};
            my $dirpfvGenscan = qq{$dirpath/$file/$version/Genscan};
            my $dirpfvFgenesh = qq{$dirpath/$file/$version/Fgenesh};
            my $dirpfvGenefinder = qq{$dirpath/$file/$version/Genefinder};

            $self->add_directory(\@cmd,$dirpfvEannot );
            $self->add_directory(\@cmd,$dirpfvGenewise );
            $self->add_directory(\@cmd,$dirpfvGenemark );
            $self->add_directory(\@cmd,$dirpfvGlimmer2 );
            $self->add_directory(\@cmd,$dirpfvGlimmer3 );
            $self->add_directory(\@cmd,$dirpfvSnap );
            $self->add_directory(\@cmd,$dirpfvGenscan );
            $self->add_directory(\@cmd,$dirpfvFgenesh );
            $self->add_directory(\@cmd,$dirpfvGenefinder );

        }
        if ( $file =~ 'Gene_merging' )
        {
            my $dirpfvHybrid = qq{$dirpath/$file/$version/Hybrid};
            my $dirpfvHIntergenic
                = qq{$dirpath/$file/$version/Hybrid/intergenic};
            $self->add_directory(\@cmd,$dirpfvHybrid );
            $self->add_directory(\@cmd,$dirpfvHIntergenic );
        }
        if ( $file =~ 'Repeats' )
        {
            my $dirpfvRECON = qq{$dirpath/$file/$version/RECON};
            my $dirpfvRQC = qq{$dirpath/$file/$version/RECON/QC};
            my $dirpfvRRepeatmasker
                = qq{$dirpath/$file/$version/RECON/RepeatMasker};
            $self->add_directory(\@cmd,$dirpfvRECON );
            $self->add_directory(\@cmd,$dirpfvRQC );
            $self->add_directory(\@cmd,$dirpfvRRepeatmasker );
        }
        if ( $file =~ 'Sequence' )
        {
            my $dirpfMasked = qq{$dirpath/$file/Masked};
            my $dirpfUnmasked = qq{$dirpath/$file/Unmasked};
            $self->add_directory(\@cmd,$dirpfMasked );
            $self->add_directory(\@cmd,$dirpfUnmasked );

        }
        if ( $file =~ 'Acedb' )
        {
            my $dirpfvwgf = qq{$dirpath/$file/$version/wgf};
            my $dirpfvwspec = qq{$dirpath/$file/$version/wspec};
            my $dirpfvDatabase = qq{$dirpath/$file/$version/database};
            my $dirpfvGff = qq{$dirpath/$file/$version/Gff_files};
            $self->add_directory(\@cmd,$dirpfvwgf );
            $self->add_directory(\@cmd,$dirpfvwspec );
            $self->add_directory(\@cmd,$dirpfvDatabase );
            $self->add_directory(\@cmd,$dirpfvGff );


        }
        if ( $file =~ 'Rfam' )
        {
            #unless ( -e $dirpfScripts ) { push @cmd, qq{$dirpfScripts}; }
            $self->add_directory(\@cmd,$dirpfScripts );

        }
        if ( $file =~ 'Kegg' )
        {
            #unless ( -e $dirpfScripts ) { push @cmd, qq{$dirpfScripts}; }
            $self->add_directory(\@cmd,$dirpfScripts );

            $dirpfvGenemark = qq{$dirpath/$file/$version/Genemark};
            #unless ( -e $dirpfvGenemark ) { push @cmd, qq{$dirpfvGenemark}; }
             $dirpfvGlimmer2 = qq{$dirpath/$file/$version/Glimmer2};
            #unless ( -e $dirpfvGlimmer2 ) { push @cmd, qq{$dirpfvGlimmer2}; }
             $dirpfvGlimmer3 = qq{$dirpath/$file/$version/Glimmer3};
            #unless ( -e $dirpfvGlimmer3 ) { push @cmd, qq{$dirpfvGlimmer3}; }
            $dirpfvHybrid = qq{$dirpath/$file/$version/Hybrid};
            #unless ( -e $dirpfvHybrid ) { push @cmd, qq{$dirpfvHybrid}; }
            $self->add_directory(\@cmd,$dirpfvGenemark );
            $self->add_directory(\@cmd,$dirpfvGlimmer2 );
            $self->add_directory(\@cmd,$dirpfvGlimmer3 );
            $self->add_directory(\@cmd,$dirpfvHybrid );
        }
        if ( $file =~ 'Cog' )
        {
            #unless ( -e $dirpfScripts ) { push @cmd, qq{$dirpfScripts}; }
            $self->add_directory(\@cmd,$dirpfScripts );

             $dirpfvGenemark = qq{$dirpath/$file/$version/Genemark};
            #unless ( -e $dirpfvGenemark ) { push @cmd, qq{$dirpfvGenemark}; }
             $dirpfvGlimmer2 = qq{$dirpath/$file/$version/Glimmer2};
            #unless ( -e $dirpfvGlimmer2 ) { push @cmd, qq{$dirpfvGlimmer2}; }
             $dirpfvGlimmer3 = qq{$dirpath/$file/$version/Glimmer3};
            #unless ( -e $dirpfvGlimmer3 ) { push @cmd, qq{$dirpfvGlimmer3}; }
             $dirpfvHybrid = qq{$dirpath/$file/$version/Hybrid};
            #unless ( -e $dirpfvHybrid ) { push @cmd, qq{$dirpfvHybrid}; }
            $self->add_directory(\@cmd,$dirpfvGenemark );
            $self->add_directory(\@cmd,$dirpfvGlimmer2 );
            $self->add_directory(\@cmd,$dirpfvGlimmer3 );
            $self->add_directory(\@cmd,$dirpfvHybrid );
        }
        if ( $file =~ 'Interpro' )
        {
            unless ( -e $dirpfScripts ) { push @cmd, qq{$dirpfScripts}; }

             $dirpfvGenemark = qq{$dirpath/$file/$version/Genemark};
#            unless ( -e $dirpfvGenemark ) { push @cmd, qq{$dirpfvGenemark}; }
             $dirpfvGlimmer2 = qq{$dirpath/$file/$version/Glimmer2};
#            unless ( -e $dirpfvGlimmer2 ) { push @cmd, qq{$dirpfvGlimmer2}; }
             $dirpfvGlimmer3 = qq{$dirpath/$file/$version/Glimmer3};
#            unless ( -e $dirpfvGlimmer3 ) { push @cmd, qq{$dirpfvGlimmer3}; }
             $dirpfvHybrid = qq{$dirpath/$file/$version/Hybrid};
#            unless ( -e $dirpfvHybrid ) { push @cmd, qq{$dirpfvHybrid}; }
            $self->add_directory(\@cmd,$dirpfvGenemark );
            $self->add_directory(\@cmd,$dirpfvGlimmer2 );
            $self->add_directory(\@cmd,$dirpfvGlimmer3 );
            $self->add_directory(\@cmd,$dirpfvHybrid );

        }
        if ( $file =~ 'psortB' )
        {
            unless ( -e $dirpfScripts ) { push @cmd, qq{$dirpfScripts}; }

             $dirpfvGenemark = qq{$dirpath/$file/$version/Genemark};
            #unless ( -e $dirpfvGenemark ) { push @cmd, qq{$dirpfvGenemark}; }
             $dirpfvGlimmer2 = qq{$dirpath/$file/$version/Glimmer2};
            #unless ( -e $dirpfvGlimmer2 ) { push @cmd, qq{$dirpfvGlimmer2}; }
             $dirpfvGlimmer3 = qq{$dirpath/$file/$version/Glimmer3};
            #unless ( -e $dirpfvGlimmer3 ) { push @cmd, qq{$dirpfvGlimmer3}; }
             $dirpfvHybrid = qq{$dirpath/$file/$version/Hybrid};
            #unless ( -e $dirpfvHybrid ) { push @cmd, qq{$dirpfvHybrid}; }
            $self->add_directory(\@cmd,$dirpfvGenemark );
            $self->add_directory(\@cmd,$dirpfvGlimmer2 );
            $self->add_directory(\@cmd,$dirpfvGlimmer3 );
            $self->add_directory(\@cmd,$dirpfvHybrid );
        }
        if ( $file =~ 'Blastp' )
        {
            unless ( -e $dirpfScripts ) { push @cmd, qq{$dirpfScripts}; }

             $dirpfvHybrid = qq{$dirpath/$file/$version/Hybrid};
#            unless ( -e $dirpfvHybrid ) { push @cmd, qq{$dirpfvHybrid}; }
            $self->add_directory(\@cmd,$dirpfvHybrid );
        }
        if ( $file =~ 'BAP' )
        {
            $dirpfvSequence = qq{$dirpath/$file/$version/Sequence};
#            unless ( -e $dirpfvSequence ) { push @cmd, qq{$dirpfvSequence}; }
             $dirpfvDumps = qq{$dirpath/$file/$version/Dumps};
#            unless ( -e $dirpfvDumps ) { push @cmd, qq{$dirpfvDumps}; }
            $self->add_directory(\@cmd,$dirpfvSequence );
            $self->add_directory(\@cmd,$dirpfvDumps );
        }
        if ( $file =~ 'GAP' )
        {
             $dirpfvSequence = qq{$dirpath/$file/$version/Sequence};
#            unless ( -e $dirpfvSequence ) { push @cmd, qq{$dirpfvSequence}; }
             $dirpfvDumps = qq{$dirpath/$file/$version/Dumps};
#            unless ( -e $dirpfvDumps ) { push @cmd, qq{$dirpfvDumps}; }
            $self->add_directory(\@cmd,$dirpfvDumps );
            $self->add_directory(\@cmd,$dirpfvSequence );
        }

        @cmd = uniq @cmd;    # make sure there are now duplicates.

        foreach my $cmd (@cmd)
        {
            #mkdir $cmd or croak "can't mkdir $cmd : $!\n";
            make_path( $cmd ,{verbose => 1, mode => 0775} ) or croak "can't mkdir $cmd : $!\n";
        }
        my $readme = $newdir . "/"
            . $assembly_dir
            . "/README_"
            . $self->assembly_version_name;
        my $message = " This file was created on $date for "
            . $self->org_dirname
            . "using assembly: "
            . $self->assembly_version_name . "\n";
        write_file( $readme, ($message) );

    }

    return 1;
}

sub add_directory
{
    my $self = shift;
    my $list = shift;
    my $dir = shift;
    unless( -e $dir)
    {
        push(@$list,$dir);
    }
    
    return 1;
}



1;

# $Id$
