package Genome::Model::Tools::Ber::BerRunAnnoSqlite;

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
use DateTime;

use Cwd;

UR::Object::Type->define(
			 class_name => __PACKAGE__,
			 is         => 'Command',
			 has        => [
					'locus_tag'       => {
							      is  => 'String',
							      doc => "Locus tag for project, containing DFT/FNL postpended",
							     },
					'srcdirpath'      => {
					                      is  => 'String',
							      doc => "src directory for the ber product naming software",
                                                             },
					'outdirpath'      => {
					                      is  => 'String',
							      doc => "output directory for the ber product naming software",
                                                             },
					'sqlitedatafile'  => { is  => 'String',
							       doc => "Name of sqlite output .dat file",
							       is_optional => 1,
							   },
					'sqliteoutfile'  => { is  => 'String',
							       doc => "Name of sqlite output .out file",
							       is_optional => 1,
							   },
				       ]
			);

sub help_brief
  {
    "Tool for Running anno-sqlite.bash for BER product naming pipeline";
  }

sub help_synopsis
  {
    return <<"EOS"
      Tool for Running anno-sqlite.bash for BER product naming pipeline.
EOS
  }

sub help_detail
  {
    return <<"EOS"
Tool for Running anno-sqlite.bash for BER product naming pipeline.
EOS
  }


sub execute
  {
      my $self = shift;

      my @berrunannosqlite = $self->gather_details();

      my $locus_tag             = $self->locus_tag;
      my $outdirpath            = $self->outdirpath;
      my $bsubanno_screenoutput = qq{$outdirpath/bsub.screenoutput.$locus_tag.$<};

      eval {

	  IPC::Run::run(
	                @berrunannosqlite,
			'>',
			$bsubanno_screenoutput,
	       		timeout ( 5 ),
		    ) || croak "\n\nCan't run BerRunAnnoSqlite.pm from BerRunAnnoSqlite.pm: $OS_ERROR\n\n";
      };

      my $screenout_fh = IO::File->new();
      $screenout_fh->open("< $bsubanno_screenoutput")
	  or die "Can't open $bsubanno_screenoutput  from BerRunAnnoSqlite.pm: $OS_ERROR\n";
      
     # my $bjobs_fh = IO::File->new();
     # $bjobs_fh->open("> bjobs.txt")
#	  or die "Can't open bjobs.txt  from BerRunAnnoSqlite.pm: $OS_ERROR\n";

      while (<$screenout_fh>) {

	  my ( $job, $jobid, $word1, $word2, $word3, $word4, $word5 ) = split(' ', $ARG);

	  $jobid =~ s/^<//;
	  $jobid =~ s/>$//;
	  
	SCREEN: while ( 1 ) {
	      my $bjobsfile = qq{bjobs_txt\_$locus_tag\.$PID};

	      my $bjobs_fh = IO::File->new();
	      $bjobs_fh->open("> $bjobsfile")
		  or die "Can't open bjobs.txt  from BerRunAnnoSqlite.pm: $OS_ERROR\n";
	      
	      my @bjobscmd = (
		  'bjobs',
		  $jobid,
	      );
	      
	      IPC::Run::run(
		  \@bjobscmd,
		  '>',
		  $bjobs_fh
	      );
	      
	      $bjobs_fh->close();

	      my $bjobs_open_fh = IO::File->new();
	      $bjobs_open_fh->open("< $bjobsfile")
		  or die "Can't open $bjobs_open_fh (bjobs file)  from BerRunAnnoSqlite.pm: $OS_ERROR\n";

	      BJOB: while (<$bjobs_open_fh>) {
	      
		  unless ($ARG =~ /^JOBID/) {

		      #print qq{\n$ARG\n};

		      my ($jobid2, $user, $stat, $queue, $host, $host2, $jobname, $submit) = split(' ', $ARG);

		      #print qq{\n$stat\n};

		      if ( ($stat eq "RUN" ) || ($stat eq "PEND")) {
			  sleep 600;
			  $bjobs_open_fh->close();
			  next SCREEN;
		      }
		      elsif ( $stat eq "DONE" ) {
			  $bjobs_open_fh->close();
			  unlink("$bjobsfile");
			  last SCREEN;
		      }
		  }
	      }
	  }

      }

      return 1;

  }

sub gather_details
  {
      my $self                  = shift;

      my $locus_tag             = $self->locus_tag;
      my $srcdirpath            = $self->srcdirpath;
      my $outdirpath            = $self->outdirpath;
      my $bsubannoout           = qq{$outdirpath/bsub.out.$locus_tag};
      my $bsubannoerr           = qq{$outdirpath/bsub.err.$locus_tag};
      my $bsubanno_screenoutput = qq{$outdirpath/bsub.screenoutput.$locus_tag};

      my $cwd = getcwd();
      unless ($cwd eq $srcdirpath) {
	  chdir($srcdirpath) or die "Failed to change to '$srcdirpath'...  from BerRunBtabhmmtab.pm: $OS_ERROR\n\n";
      }

      my $now      = DateTime->now(time_zone => 'America/Chicago');
      my $datecode = $now->ymd('_');

      my @annosqlitecmd = (
	                   './anno-sqlite.bash',
			   $locus_tag,
			   $datecode,
		          );

      my $sqlitedata = qq{sqlite-$locus_tag-$datecode.dat};
      $self->sqlitedatafile($sqlitedata);

      my $sqliteout  = qq{sqlite-$locus_tag-$datecode.out};
      $self->sqliteoutfile($sqliteout);

      my $Rbp    = qq{'select[mem>8000] rusage[mem=8000]'};

      my @bsubcmdb = (
	              'bsub',
		      '-o',
		      $bsubannoout,
		      '-e',
		      $bsubannoerr,
		      '-q',
		      'long',
              '-M',
              '8000000',
		      '-n',
		      '1',
		      '-R',
		      $Rbp,
		      @annosqlitecmd,
		     );
      
      my @annoipc = (
	  \@bsubcmdb,
	  \undef,
	  '2>&1',
      );
      
      
      return @annoipc;
      
  }
