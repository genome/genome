package Genome::Model::Tools::Ber::BerRunBlastphmmpfam;

use strict;
use warnings;

use Genome;
use Command;

use Carp;
use English;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::BioDB;
use Bio::DB::Query::BioQuery;

use Data::Dumper;
use IO::File;
use IPC::Run qw/ run timeout /;
use Time::HiRes qw(sleep);

use Cwd;

UR::Object::Type->define(
			 class_name => __PACKAGE__,
			 is         => 'Command',
			 has        => [
					'locus_tag'       => {
							      is  => 'String',
							      doc => "Locus tag for project, containing DFT/FNL postpended",
							     },
					'proteinnamefile' => {
							      is  => 'String',
							      doc => "Protein name file for all genes dumped from BioSQL",
							     },
					'fastadirpath'    => {
							      is  => 'String',
							      doc => "fasta directory path",
							     },
					'berdirpath'      => {
							       is  => 'String',
							       doc => "blastp directory path",
							     } ,
					'hmmdirpath'      => {
					                       is  => 'String',
							       doc => "hmmpfam directory path",
                                                             },
					'bsubfiledirpath' => {
							      is  => 'String',
							      doc => "bsub error files directory path",
							     } ,
					'blpqueryfile'    => {
                                                              is => 'String',
							      doc => "blastp query file",
                                                             } ,
					'hmmdatabase'     => {
					                      is => 'String',
							      doc => "hmmpfam query file",
					                      } ,
				       ]
			);

sub help_brief
  {
    "Tool for Running Blastp and hmmpfam for BER product naming pipeline";
  }

sub help_synopsis
  {
    return <<"EOS"
      Tool for Running Blastp and hmmpfam for BER product naming pipeline.
EOS
  }

sub help_detail
  {
    return <<"EOS"
Tool for Running Blastp and hmmpfam for BER product naming pipeline.
EOS
  }


sub execute
  {
    my $self = shift;

    my @berrunblastphmmpfam = $self->gather_details();

    my $i = 0;
    foreach my $runs (@berrunblastphmmpfam) {
	$i++;
	#print qq{\n\n@$runs\n\n};

	#print Dumper(@$runs);

	#print qq{\n\n$i\n\n};
	IPC::Run::run(@$runs) || croak "\n\nCan't run BerRunBlastphmmpfam.pm from BerRunBlastphmmpfam.pm: $OS_ERROR\n\n";

    }

    print qq{\n\nTotal Jobs for Blastp/Hmmpfam:\t$i\n\n};

    while ( 1 ) {

	my $berdirpath = $self->berdirpath;
	my $hmmdirpath = $self->hmmdirpath;

	opendir (BERDIR, "$berdirpath")or die "Can't open '$berdirpath' ...  from BerRunBlastphmmpfam.pm : $OS_ERROR\n";
	my @berfile = readdir (BERDIR);
	closedir (BERDIR);

	my @storageber = ( );
	foreach my $bername (@berfile) {
	    next if $bername =~ m/^\.\.?$/;
	    my $nr           = qq{.blastp};

	    if ( $bername =~ /$nr$/ ) {
		my $bernamecheck = qq{$berdirpath/$bername};
		if ( (-e $bernamecheck) && (! -z $bernamecheck)) {
		    my $size = -s $bernamecheck;
		    if ( $size > 550 ) {
			push (@storageber, $bername);
		    }
		}
	    }
	}

	my $bercount = $#storageber + 1;

	opendir (HMMDIR, "$hmmdirpath")or die "Can't open '$hmmdirpath' ...  from BerRunBlastphmmpfam.pm : $OS_ERROR\n";
	my @hmmfile = readdir (HMMDIR);
	closedir (HMMDIR);

	my @storagehmm = ( );
	foreach my $hmmname (@hmmfile) {
	    next if $hmmname =~ m/^\.\.?$/;
	    my $hmmpfam      = qq{.hmmpfam};
	    if  ($hmmname =~ /$hmmpfam$/) {
		my $hmmnamecheck = qq{$hmmdirpath/$hmmname};
		if (( -e $hmmnamecheck) && ( ! -z $hmmnamecheck )){
		    my $size = -s $hmmnamecheck;
		    if ( $size > 550 ){
			push (@storagehmm, $hmmname);
		    }

            # need to check if there is an error here, then try to resubmit or something???
            my $chkcmd = "grep -H \"Exited with exit code\" $hmmnamecheck";
            my $rc = system($chkcmd);
            if($rc == 0) {
                warn "errors in $hmmnamecheck - this portion may need to be rerun";
            }
		}
	    }
	}

	my $hmmcount = $#storagehmm + 1;

	my $totalcount = $bercount + $hmmcount;
	unless ( $i == $totalcount ) {
	    sleep 900;
	    next;
	}
	else {
	    sleep 120;
	    last;
	}
    }

    return 1;

  }

sub gather_details
  {
    my $self = shift;
    my $locus_tag       = $self->locus_tag;
    my $proteinnamefile = $self->proteinnamefile;
    my $fastadirpath    = $self->fastadirpath;
    my $berdirpath      = $self->berdirpath;
    my $hmmdirpath      = $self->hmmdirpath;
    my $bsubfiledirpath = $self->bsubfiledirpath;
    my $blpqueryfile    = $self->blpqueryfile;
    my $hmmdatabase     = $self->hmmdatabase;

    my @allipc;
    my $cwd = getcwd();
    unless ( $cwd eq $fastadirpath ) {
	chdir($fastadirpath)  or croak "Failed to change to '$fastadirpath', from BerRunBlastphmmpfam.pm: $OS_ERROR";
    }

    # FIXME: this sort of hard codes the path.  should test if 
    # proteinnamefile is a full path to valid file.
    $proteinnamefile = qq{$fastadirpath/$proteinnamefile};
    unless ( -e $proteinnamefile ) {
      croak qq{\n\n NO, $proteinnamefile file found, from BerRunBlastphmmpfam.pm : $OS_ERROR \n\n};
    }
    my $proteinnamefile_fh = IO::File->new();
    $proteinnamefile_fh->open("<$proteinnamefile")
      or croak "Can't open '$proteinnamefile' from BerRunBlastphmmpfam.pm : $OS_ERROR";;

    my $blpcount  = 0;
    my $hmpfcount = 0;

    while (my $line = <$proteinnamefile_fh>) {

      chomp($line);
      my $file = qq{$fastadirpath/$line};

      unless ( -e $file ) {
	  croak qq{\n\n NO, $file, found for blastp from BerRunBlastphmmpfam.pm : $OS_ERROR \n\n };
      }

      #blastp send to lsf

	  $line =~ s/\.fasta$//;
      my $blperr        = qq{$bsubfiledirpath/bsub.err.blp.$line};
      my $blpout        = qq{$berdirpath/$line.blastp};
     
      my @blastpcmd = (
		       'blastp',
		       $blpqueryfile,
		       $file,
		      );

      my $Rbp    = qq {rusage[mem=4096]};

      my @bsubcmdbp = (
		       'bsub',
		       '-o',
		       $blpout,
		       '-e',
		       $blperr,
		       '-q',
		       $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
		       '-n',
		       '1',
		       '-R',
		       $Rbp,
		       @blastpcmd,
		      );

      $blpcount++;

      my @bpipc = (
	          \@bsubcmdbp,
	          \undef,
	          '2>&1',
	          );

      push(@allipc, \@bpipc);

      #hmmpfam send to lsf

      my $hmpfmout  = qq{$hmmdirpath/$line.hmmpfam};
      my $hmpfmerr  = qq{$bsubfiledirpath/bsub.err.hmm.$line};
     
      my @hmmpfcmd = (
	             'hmmpfam',
                 '--cpu', '1',
		     $hmmdatabase,
		     $file,
		     );

      my @bsubcmdhf = (
	              'bsub',
		      '-o',
		      $hmpfmout,
		      '-e',
		      $hmpfmerr,
		      '-q',
		      $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
		      '-n',
		      '1',
		      '-R',
		      $Rbp,
		      @hmmpfcmd,
		      );

      my @hmmpfipc = (
	           \@bsubcmdhf,
		   \undef,
		   '2>&1',
                   );

      push(@allipc, \@hmmpfipc );

  }

    return @allipc;

}

