# BioPerl module for Bio::Tools::Run::Glimmer
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

Bio::Tools::Run::Glimmer - Wrapper for local execution of Glimmer,
GlimmerM and GlimmerHMM.

=head1 SYNOPSIS

  # glimmer2
  my $factory =
      Bio::Tools::Run::Glimmer->new('-program' => 'glimmer3',
                                    '-model'   => 'model.icm');
  # glimmer3
  my $factory =
      Bio::Tools::Run::Glimmer->new('-program' => 'glimmer2',
                                    '-model'   => 'model.icm');
  # glimmerm
  my $factory =
      Bio::Tools::Run::Glimmer->new('-program' => 'glimmerm');

  # glimmerHMM
  my $factory =
      Bio::Tools::Run::Glimmer->new('-program' => 'glimmerHMM');

  # Pass the factory Bio::Seq objects
  # returns a Bio::Tools::Glimmer object
  my $glimmer = $factory->run($seq);
  or
  my $glimmer = $factor->run(@seq);

=head1 DESCRIPTION

Wrapper module for the Glimmer family of programs.  Should work with
all currently available flavors: Glimmer, GlimmerM and GlimmerHMM.
However, only Glimmer 2.X and 3.X have been tested.

Glimmer is open source and available at
L<http://www.cbcb.umd.edu/software/glimmer/>.

GlimmerM is open source and available at
L<http://www.tigr.org/software/glimmerm/>.

GlimmerHMM is open source and available at
L<http://www.cbcb.umd.edu/software/GlimmerHMM/>.

Note that Glimmer 2.X will only process the first
sequence in a fasta file (if you run() more than one
sequence at a time, only the first will be processed).

Note that Glimmer 3.X produces two output files.  This 
module only passes the .predict file to Bio::Tools::Glimmer.

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

package Bio::Tools::Run::Glimmer;

use strict;
use warnings;

use Bio::SeqIO;
use Bio::Root::Root;
use Bio::Tools::Run::WrapperBase;
use Bio::Tools::Glimmer;
use English;
use IPC::Run;     # Should be okay on WIN32 (See IPC::Run Docs)

use base qw(Bio::Root::Root Bio::Tools::Run::WrapperBase);

our @params              = (qw(program model));

our @glimmer2_params     = (qw(C L g i o p q s t w));
our @glimmer2_switches   = (qw(M X f l r));

our @glimmer3_params     = (qw(A C E L M P Z b g i t z));
our @glimmer3_switches   = (qw(X f l o q r));

our @glimmerM_params     = (qw(d g t));
our @glimmerM_switches   = (qw(5 3 f r s));

our @glimmerHMM_params   = (qw(d n p));
our @glimmerHMM_switches = (qw(f h v));

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

=head2 model

 Title   : model
 Usage   : $factory>model()
 Function: gets/sets the name of the model (icm) file
 Returns:  string
 Args    : string

=cut

sub model {
    
    my ($self, $val) = @_;
    
    $self->{'_model'} = $val if $val;
    
    return $self->{'_model'};

}

=head2 new

 Title   : new
 Usage   : $glimmer->new(@params)
 Function: creates a new Glimmer factory
 Returns:  Bio::Tools::Run::Glimmer
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
                                          @glimmer2_params,
                                          @glimmer2_switches,
                                          @glimmer3_params,
                                          @glimmer3_switches,
                                          @glimmerM_params,
                                          @glimmerM_switches,
                                          @glimmerHMM_params,
                                          @glimmerHMM_switches
                                         ],
                             -create =>  1,
                         );

       unless (defined($self->program())) {
           $self->throw('Must specify program');
       }

       unless (defined($self->model())) {
           $self->throw('Must specify model');
       }
       
       return $self;
       
}

=head2 run

 Title   :   run
 Usage   :   $obj->run($seq_file)
 Function:   Runs Glimmer/GlimmerM/GlimmerHMM
 Returns :   A Bio::Tools::Glimmer object
 Args    :   An array of Bio::PrimarySeqI objects

=cut

sub run{
    
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

    my @run_args = ( $file_name );

    # Glimmer 2.X ignores sequences after the first in a fasta file
    # Glimmer 3.X will process multiple sequences at once
    if ($program_name eq 'glimmer2') {
        if (@seq > 1) {
            $self->warn("Program $program_name processes one sequence at a time");
        }
        push @run_args, $seq[0]->display_id();
        push @run_args, $seq[0]->length();
    }
    
    return $self->_run(@run_args);
    
}

=head2 _run

 Title   :   _run
 Usage   :   $obj->_run()
 Function:   Internal(not to be used directly)
 Returns :   An instance of Bio::Tools::Glimmer 
 Args    :   file name, sequence identifier (optional)

=cut

sub _run {
    
    my ($self, $seq_file_name, $seq_id, $seq_length) = @_;
    
    my @cmd = (
               $self->executable(),
               $seq_file_name,
               $self->model(),
               split(/\s+/, $self->_setparams()),
           );

    my $cmd = join(' ', @cmd);
    $self->debug("Glimmer Command = $cmd");
    
    my $program_name = $self->program_name();
    my ($output_fh, $output_file_name, $detail_file_name);
    my ($program_stdout, $program_stderr);
        
    my @ipc_args = (\@cmd, \undef);
    
    # No STDOUT option for glimmer3, it takes a
    # 'tag' argument, and outputs tag.predict and
    # tag.detail.  It seems that tag can be a path,
    # which is handy.  
    if ($program_name eq 'glimmer3') {
        
        my $temp_dir     = $self->tempdir();
        my $glimmer3_tag = "$temp_dir/glimmer3";

        push @cmd, $glimmer3_tag;     
        $output_file_name = "$glimmer3_tag.predict";
        $detail_file_name = "$glimmer3_tag.detail";
        push @ipc_args, \$program_stdout, \$program_stderr;
        
    }
    else {
        
        ($output_fh, $output_file_name) = $self->io->tempfile(-dir=>$self->tempdir());
        close($output_fh);
        push @ipc_args, '>', $output_file_name;
        push @ipc_args, '2>', \$program_stderr;
        
    }

    # Run the program via IPC::Run so:
    # 1) The console doesn't get cluttered up with the program's STDERR/STDOUT
    # 2) We don't have to embed STDERR/STDOUT redirection in $cmd
    # 3) We don't have to deal with signal handling (IPC::Run should take care
    #    of everything automagically. 

    eval {
        IPC::Run::run(@ipc_args) || die $CHILD_ERROR;;
    };

    if ($EVAL_ERROR) {
        $self->throw("Glimmer call crashed: $EVAL_ERROR"); 
    }

    $self->debug(join("\n", 'Glimmer STDOUT:', $program_stdout)) if $program_stdout;
    $self->debug(join("\n", 'Glimmer STDERR:', $program_stderr)) if $program_stderr;
    
    my %parser_args = (-file => $output_file_name);

    # Pass along $seq_id and $seq_length if they were provided
    # (only should be for glimmer2).
    if (defined($seq_id))     { $parser_args{-seqname  } = $seq_id;     }
    if (defined($seq_length)) { $parser_args{-seqlength} = $seq_length; }

    # Pass along the name of extra output file, with handy information about
    # sequence lengths (only produced by glimmer3)
    if (defined($detail_file_name)) { $parser_args{-detail} = $detail_file_name; }
    
    return Bio::Tools::Glimmer->new(%parser_args);
    
}

sub _setparams {

    my ($self) = @_;

    my $param_string = $self->SUPER::_setparams(
                                                -params   => [
                                                              @glimmer2_params,
                                                              @glimmer3_params,
                                                              @glimmerM_params,
                                                              @glimmerHMM_params,
                                                             ],
                                                -switches => [
                                                              @glimmer2_switches,
                                                              @glimmer2_switches,
                                                              @glimmerM_switches,
                                                              @glimmerHMM_switches,
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
