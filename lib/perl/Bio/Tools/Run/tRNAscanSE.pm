# BioPerl module for Bio::Tools::Run::tRNAscanSE
#
# Cared for by Bioperl
#
# Copyright Bioperl, Mark Johnson <mjohnson-at-watson-dot-wustl-dot-edu>
#
# Special thanks to Chris Fields, Sendu Bala
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::tRNAscanSE - Wrapper for local execution of tRNAscan-SE

=head1 SYNOPSIS

  my $factory = Bio::Tools::Run::tRNAscanSE->new('-program' => 'tRNAscan-SE');
 
  # Pass the factory Bio::Seq objects
  # returns a Bio::Tools::tRNAscanSE object
  my $factory = $factory->run($seq);
  or
  my $factory = $factory->run(@seq);

=head1 DESCRIPTION

Wrapper module for tRNAscan-SE.

tRNAscan-SE is open source and available at
L<http://selab.wustl.edu/cgi-bin/selab.pl?mode=software#trnascan/>.

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

package Bio::Tools::Run::tRNAscanSE;

use strict;
use warnings;

use Bio::SeqIO;
use Bio::Root::Root;
use Bio::Tools::Run::WrapperBase;
use Bio::Tools::tRNAscanSE;
use English;
use IPC::Run;     # Should be okay on WIN32 (See IPC::Run Docs)

use base qw(Bio::Root::Root Bio::Tools::Run::WrapperBase);

our @params              = (qw(program));
our @tRNAscanSE_switches = (qw(A B C G O P));

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
 Usage   : $tRNAscanSE->new(@params)
 Function: creates a new tRNAscanSE factory
 Returns:  Bio::Tools::Run::tRNAscanSE
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
                                          @tRNAscanSE_switches,
                                         ],
                             -create =>  1,
                         );

       unless (defined($self->program())) {
           $self->throw('Must specify program');
       }

       return $self;
       
}

=head2 run

 Title   :   run
 Usage   :   $obj->run($seq_file)
 Function:   Runs tRNAscan-SE
 Returns :   A Bio::Tools::tRNAscanSE object
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

    return $self->_run($file_name);
    
}

=head2 _run

 Title   :   _run
 Usage   :   $obj->_run()
 Function:   Internal(not to be used directly)
 Returns :   An instance of Bio::Tools::tRNAscanSE 
 Args    :   file name

=cut

sub _run {
    
    my ($self, $seq_file_name) = @_;
    
    my @cmd = (
               $self->executable(),
               split(/\s+/, $self->_setparams()),
               $seq_file_name,
           );

    my $cmd = join(' ', @cmd);
    $self->debug("tRNAscan-SE Command = $cmd");
    
    my $program_name = $self->program_name();
    my ($program_stderr);

    my ($output_fh, $output_file_name) = $self->io->tempfile(-dir=> $self->tempdir());
    
    
    my @ipc_args = (\@cmd, \undef, '>', $output_file_name, '2>', \$program_stderr);

    # Run the program via IPC::Run so:
    # 1) The console doesn't get cluttered up with the program's STDERR/STDOUT
    # 2) We don't have to embed STDERR/STDOUT redirection in $cmd
    # 3) We don't have to deal with signal handling (IPC::Run should take care
    #    of everything automagically. 

    eval {
        IPC::Run::run(@ipc_args) || die $CHILD_ERROR;
    };

    if ($EVAL_ERROR) {
        $self->throw("tRNAscan-SE call crashed: $program_stderr"); 
    }

    $self->debug(join("\n", 'tRNAscanSE STDERR:', $program_stderr)) if $program_stderr;
        
    return Bio::Tools::tRNAscanSE->new(-file => $output_file_name);
    
}

sub _setparams {

    my ($self) = @_;

    my $param_string = $self->SUPER::_setparams(
                                                -params   => [ ],
                                                -switches => [
                                                              @tRNAscanSE_switches,
                                                             ],
                                                -dash     => 1

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
