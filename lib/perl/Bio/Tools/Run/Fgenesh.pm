# $Id$
#
# BioPerl module for Bio::Tools::Run::Fgenesh
#
# Cared for by Bioperl
#
# Copyright Bioperl, Mark Johnson <mjohnson-at-watson-dot-wustl-dot-edu>
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::Fgenesh - Wrapper for local execution of
FGENESH/FGENESH+

=head1 SYNOPSIS

  # FGENESH
  my $factory =
      Bio::Tools::Run::Fgenesh->new('-program' => 'fgenesh',
                                    '-param'   => 'C_elegans');

  # FGENESH+
  my $factory =
      Bio::Tools::Run::Fgenesh->new('-program' => 'fgenesh+',
                                    '-param'   => 'C_elegans',
                                    '-O'       => 'prot_map.CFG');

  # Pass the factory Bio::Seq objects
  # returns a Bio::Tools::Fgenesh object
  my $fgenesh = $factory->run($seq);

=head1 DESCRIPTION

Wrapper module for the FGENESH and FGENESH+ gene predictors.

General information about FGENESH is available at
L<http://www.softberry.com/berry.phtml?topic=fgenesh&group=help&subgroup=gfind>.

General information about FGENESH+ is available at
L<http://http://www.softberry.com/berry.phtml?topic=fgenes_plus&group=help&subgroup=gfs>.

Contact information for licensing inquiries is available at:
L<http://www.softberry.com/berry.phtml?topic=contact>.

Note that FGENESH/FGENESH+ will throw an error if it encounters a
multi-fasta file (if you run() more than one sequence at a time,
only the first will be processed).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Mark Johnson 

 Email: mjohnson-at-watson-dot-wustl-dot-edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::Run::Fgenesh;

use strict;
use warnings;

use Bio::SeqIO;
use Bio::Root::Root;
use Bio::Tools::Run::WrapperBase;
use Bio::Tools::Fgenesh;
use English;
use IPC::Run;     # Should be okay on WIN32 (See IPC::Run Docs)

use base qw(Bio::Root::Root Bio::Tools::Run::WrapperBase);

our @params           = (qw(program param));
our @fgenesh_params   = (qw(i m p));
our @fgenesh_switches = (qw(a n r));

=head2 program_name

 Title   : program_name
 Usage   : $factory>program_name()
 Function: gets/sets the program name
 Returns:  string
 Args    : string

=cut

sub program_name {
    
    my ($self, $val) = @_;
    
    $self->program($val) if $val;
    
    return $self->program();

}

=head2 program_dir

 Title   : program_dir
 Usage   : $factory->program_dir()
 Function: gets/sets the program dir
 Returns:  string
 Args    : string

=cut

sub program_dir {
    
    my ($self, $val) = @_;
    
    $self->{'_program_dir'} = $val if $val;
    
    return $self->{'_program_dir'};
    
}

=head2 new

 Title   : new
 Usage   : $fgenesh->new(@params)
 Function: creates a new Fgenesh factory
 Returns:  Bio::Tools::Run::Fgenesh
 Args    : 

=cut

sub new {
    
       my ($class,@args) = @_;
       my $self = $class->SUPER::new(@args);
       
       $self->io->_initialize_io();
       
       $self->_set_from_args(
                             \@args,
                             -methods => [
                                          @params,
                                          @fgenesh_params,
                                          @fgenesh_switches,
                                         ],
                             -create  => 1,
                         );

       my $program = $self->program();
       
       unless (defined($program)) {
           $self->throw('Must specify program');
       }
       
       unless (defined($self->param())) {
           $self->throw('Must specify param (parameters file)');
       }

       if ($program eq 'fgenesh+') {
           unless (defined($self->O())) {
               $self->throw('Must specify -O (prot_map) configuration file for fgenesh+');
           }
       }
       
       return $self;
       
}

=head2 run

 Title   :   run
 Usage   :   $obj->run($seq_file)
 Function:   Runs Fgenesh 
 Returns :   A Bio::Tools::Fgenesh object
 Args    :   An array of Bio::PrimarySeqI objects

=cut

sub run {
    
    my ($self, @seq) = @_;

    unless (@seq) {
        $self->throw("Must supply at least one Bio::PrimarySeqI");
    }
    
    foreach my $seq (@seq) {
        
        unless ($seq->isa('Bio::PrimarySeqI')) {
            $self->throw("Object does not implement Bio::PrimarySeqI");
        }
        
    }

    my $program_name = $self->program_name();
    my $file_name    = $self->_write_seq_file(@seq);

    # fgenesh (and fgenesh+) don't do multi-FASTA.
    if (@seq > 1) {
        $self->warn("Program $program_name processes one sequence at a time");
    }
    
    return $self->_run($file_name);
    
}

=head2 _run

 Title   :   _run
 Usage   :   $obj->_run()
 Function:   Internal(not to be used directly)
 Returns :   An instance of Bio::Tools::Fgenesh
 Args    :   file name, sequence identifier (optional)

=cut

sub _run {
    
    my ($self, $seq_file_name) = @_;

    my ($temp_fh, $temp_file_name) =
        $self->io->tempfile(-dir=>$self->tempdir());
    close($temp_fh);

    my $param_file = $self->param();
    
    # IPC::Run wants an array where the first element is the executable
    my @cmd = (
               $self->executable(),
               $param_file,
               $seq_file_name,
               split(/\s+/, $self->_setparams()),
           );

    my $cmd = join(' ', @cmd);
    $self->debug("Fgenesh Command = $cmd");

    # Run the program via IPC::Run so:
    # 1) The console doesn't get cluttered up with the program's STDERR/STDOUT
    # 2) We don't have to embed STDERR/STDOUT redirection in $cmd
    # 3) We don't have to deal with signal handling (IPC::Run should take care
    #    of everything automagically. 
    my $program_stderr;
    
    eval {
        IPC::Run::run(
                      \@cmd,
                      \undef,
                      '>',
                      $temp_file_name,
                      '2>',
                      \$program_stderr,
                  ) || die $CHILD_ERROR;
        
    };
    
    if ($EVAL_ERROR) {
        $self->throw("Fgenesh call crashed: $EVAL_ERROR"); 
    }

    $self->debug(join("\n", 'Fgenesh STDERR:', $program_stderr)) if $program_stderr;
                 
    return Bio::Tools::Fgenesh->new(-file => $temp_file_name);
    
}

sub _setparams {
    
    my ($self) = @_;
    
    my $param_string = $self->SUPER::_setparams(
                                                -params   => [@fgenesh_params],
                                                -switches => [@fgenesh_switches],
                                                -dash     => 1,
                                            );

    # Kill leading and trailing whitespace
    $param_string =~ s/^\s+//g;
    $param_string =~ s/\s+$//g;

    return $param_string;
    
}

=head2 _write_seq_file

 Title   :   _write_seq_file
 Usage   :   obj->_write_seq_file($seq) or obj->_write_seq_file(@seq)
 Function:   Internal(not to be used directly)
 Returns :   Name of a temp file containing program output
 Args    :   One or more Bio::PrimarySeqI objects

=cut

sub _write_seq_file {

    my ($self, @seq) = @_;
    
    my ($fh, $file_name) = $self->io->tempfile(-dir=>$self->tempdir());
    my $out              = Bio::SeqIO->new(-fh => $fh , '-format' => 'Fasta');

    foreach my $seq (@seq){
	$out->write_seq($seq);
    }

    close($fh);
    $out->close();
    
    return $file_name;

}

1;
