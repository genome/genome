#$Id$

package Genome::Model::GenePrediction::Command::Pap::BtabHmmtab;

use strict;
use warnings;


use English;
use File::Basename;
use File::chdir;
use File::Temp;
use File::Slurp;
use IO::File;
use IO::Dir;
use IPC::Run;
use Carp;

class Genome::Model::GenePrediction::Command::Pap::BtabHmmtab {
    is  => ['Command::V1'],
    has => [
        locus_tag => {
            is  => 'SCALAR',
            doc => 'locus tag name',
            is_input => 1,
        },
        hmmdirpath => {
            is  => 'SCALAR',
            doc => 'hmm result files path',
            is_input => 1,
        },
        berdirpath => {
            is  => 'SCALAR',
            doc => 'ber/blastp result files path',
            is_input => 1,
        },
        srcdirpath => {
            is  => 'SCALAR',
            doc => 'ber annotation base path (scripts/libs located here)',
            is_input => 1,
        },
        bsubfiledirpath => {
            is  => 'SCALAR',
            doc => 'bsub error files path',
            is_input => 1,
        },
        fastadir => {
            is  => 'SCALAR',
            doc => 'directory containing fasta files',
            is_input => 1,
        },
        sequence_names => {
            is  => 'ARRAY',
            doc => 'a list of sequence names for running hmmpfam on',
            is_input => 1,
        },
        success => {
            is          => 'SCALAR',
            doc         => 'success flag',
            is_optional => 1,
            is_output => 1,
        },
        lsf_queue => {
            is_param => 1,
            default_value => $ENV{GENOME_LSF_QUEUE_SHORT},
        },
        lsf_resource => {
            is_param => 1,
            default_value => '-R "rusage[tmp=100]" -n 1'
        },
    ],
};


sub sub_command_sort_position {10}

sub help_brief
{
    "Run run the BER hmmpfam step";
}

sub help_synopsis
{
    return <<"EOS"
    my \$w = Workflow::Model->create_from_xml("some xml...");
    my \$result = \$w->execute( 'locus tag'       => "locustag",
                              'fastadir'        => \$self->fastadirpath,
                              'berdirpath'      => \$self->berdirpath,
                              'hmmdirpath'      => \$self->hmmdirpath,
                              'srcdirpath'      => \$self->srcdirpath,
                              'bsubfiledirpath' => \$self->bsubfiledirpath,
                              'seq names'       => \@sequence_names, );

EOS
}

sub help_detail
{
    return <<"EOS"
Need documenation here.
EOS
}

sub execute
{

    my $self = shift;

    # convert ber to btab stuff
    $self->ber2btab();
    # convert hmm to htab stuff.
    $self->hmm2htab();
    $self->success(1);
    return 1;

}

sub get_fasta_store
{
    my $self = shift;
    my $locus_tag = $self->locus_tag;
    my $d = IO::Dir->new($self->fastadir);
    my $error_state = $OS_ERROR;
    my @fasta_store = ();
    my %tmpstore;
    if(defined $d)
    {
        #@fasta_store = $d->read;
        %tmpstore = map { $_ => 1 } $d->read;
        $d->close;
    }
    else
    {
        croak "can't open ".$self->fastadir. " : ". $error_state;
    }
    delete $tmpstore{'.'}; 
    delete $tmpstore{'..'}; 
    foreach my $key (sort keys %tmpstore)
    {
        if($key =~ /fof$/)
        {
            delete $tmpstore{$key};
            next;
        }
        if($key =~ /^$locus_tag/)
        {
            push(@fasta_store,$key);
        }
    }
    return @fasta_store;
}

sub ber2btab
{
    my $self = shift;
    {
        local $CWD = $self->srcdirpath;
        my @fasta_store = $self->get_fasta_store();
        # get all the files from $self->berdirpath that end in .nr
        my @ber_store = $self->get_ber_store();

        my $fasta_count = scalar @fasta_store;
        my $ber_count = scalar @ber_store;
        unless($fasta_count == $ber_count)
        {
            croak "mismatch between the fasta count and the blastp result count ( $fasta_count , $ber_count )";
        }

        my $btabcount = 0;
        foreach my $btabfile (@ber_store)
        {

            my $datestamp = join('.',time,$PID);
            my $perllib = $self->srcdirpath."/lib";
            my $btabout = $self->berdirpath."/".$btabfile.".btab"; 
            my $btabin  = $self->berdirpath."/".$btabfile;
            my $btab    = $self->srcdirpath."/wu-blast2btab.pl";
            my $btablog = $self->srcdirpath."/".$self->locus_tag.".btab.".$datestamp.".log";

            my @ber2btabcmd = (
                               'perl',
                               '-I', $perllib,
                               $btab,
                               '--input', $btabin, '--output', $btabout,
                               '--log',   $btablog,
                               );

            #print join(' ',@ber2btabcmd),"\n";
            IPC::Run::run(\@ber2btabcmd,
                          \undef,
                          '2>&1',) or 
                    croak "can't run btab conversion ". 
                          join(' ',@ber2btabcmd). " : ". $CHILD_ERROR;
            $btabcount++;
        }

        # check $btabcount == $ber_count ?
       
    }
    return 1;
}

sub hmm2htab
{
    my $self = shift;

    {
        local $CWD = $self->srcdirpath;
        my @fasta_store = $self->get_fasta_store();
        # get all the files from $self->berdirpath that end in .nr
        my @hmm_store = $self->get_hmm_store();

        my $fasta_count = scalar @fasta_store;
        my $hmm_count = scalar @hmm_store;
        unless($fasta_count == $hmm_count)
        {
            croak "mismatch between the fasta count and the hmmpfam result count ( $fasta_count , $hmm_count )";
        }
        my $htabcount = 0;
        my $htabinfo = $self->srcdirpath."/hmm_info.txt";
        # check for htabinfo???

        foreach my $htabfile (@hmm_store)
        {

            my $datestamp = join('.',time,$PID);
            my $htabout = $self->hmmdirpath."/".$htabfile.".htab"; 
            my $htabin  = $self->hmmdirpath."/".$htabfile;
            my $htab    = $self->srcdirpath."/hmmToHtab.pl";
            my $htablog = $self->srcdirpath."/".$self->locus_tag.".htab.".$datestamp.".log";

            my @hmm2htab = (
                               $htab,
                               '-H', $htabinfo,
                               '-f',$htabin, );

            #print join(' ',@ber2btabcmd),"\n";
            my ($htab_out,$htab_err);
            IPC::Run::run(\@hmm2htab,
                          \undef,
                          #'2>&1',
                          '2>',
                          \$htab_err,
                          '>',
                          \$htab_out,
                          ) or 
                    croak "can't run htab conversion ". 
                          join(' ',@hmm2htab). " : $htab_err : ". $CHILD_ERROR;
            print "#",$htab_out;
            $htabcount++;
            # write out the output to logs...
        }

        # check $htabcount == $hmm_count ?

    }
    return 1;
}

sub get_ber_store
{
    my $self = shift;
    my $locus_tag = $self->locus_tag;
    my $d = IO::Dir->new($self->berdirpath);
    my $error_state = $OS_ERROR;
    my @ber_store = ();
    my %tmpstore;
    if(defined $d)
    {
        #@fasta_store = $d->read;
        %tmpstore = map { $_ => 1 } $d->read;
        $d->close;
    }
    else
    {
        croak "can't open ".$self->berdirpath. " : ". $error_state;
    }
    delete $tmpstore{'.'}; 
    delete $tmpstore{'..'}; 
    foreach my $key (sort keys %tmpstore)
    {
        if($key =~ /\.nr$/)
        {
            my $file = $self->berdirpath ."/".$key;
            if((-e $file) && (! -z $file))
            {
                push(@ber_store,$key);
            }
        }
    }
    return @ber_store;
}

sub get_hmm_store
{
    my $self = shift;
    my $locus_tag = $self->locus_tag;
    my $d = IO::Dir->new($self->hmmdirpath);
    my $error_state = $OS_ERROR;
    my @hmm_store = ();
    my %tmpstore;
    if(defined $d)
    {
        #@fasta_store = $d->read;
        %tmpstore = map { $_ => 1 } $d->read;
        $d->close;
    }
    else
    {
        croak "can't open ".$self->hmmdirpath. " : ". $error_state;
    }
    delete $tmpstore{'.'}; 
    delete $tmpstore{'..'}; 
    foreach my $key (sort keys %tmpstore)
    {
        if($key =~ /\.hmmpfam$/)
        {
            my $file = $self->hmmdirpath ."/".$key;
            if((-e $file) && (! -z $file))
            {
                push(@hmm_store,$key);
            }
        }
    }
    return @hmm_store;
}



1;
