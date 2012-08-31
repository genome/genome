package GAP::Job::RNAmmer;

use strict;
use warnings;

use GAP::Job;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Carp 'confess';
use English;
use File::Temp qw/tempdir/;
use IO::File;
use IPC::Run;

use base qw(GAP::Job);


sub new {
    my ($class, $seq, $domain, $job_id, $version) = @_;

    my $self = {};
    bless $self, $class;
    
    unless (defined($seq)) {
        confess 'missing seq object!';
    }
    
    unless ($seq->isa('Bio::PrimarySeqI')) {
        confess 'seq object is not a Bio::PrimaySeqI!';
    }
    $self->{_seq} = $seq;

    unless (defined($domain)) {
        confess 'missing domain';
    }
    unless (($domain eq 'archaea') || ($domain eq 'bacteria') || ($domain eq 'eukaryota')) {
        confess "invalid domain '$domain'";
    }
    $self->{_domain} = $domain;
    
    unless (defined($job_id)) {
        confess 'missing job id';
    }
    $self->job_id($job_id);
    
    unless (defined $version) {
        $version = '1.2.1';
    }
    $self->{_version} = $version;

    return $self;
    
}

sub execute {
    my ($self) = @_;

    $self->SUPER::execute(@_);
   
    my $seq    = $self->{_seq};
    my $domain = $self->{_domain};
    my $seq_fh = $self->_write_seqfile($seq);
    
    my ($rnammer_stdout, $rnammer_stderr);

    my $temp_fh       = File::Temp->new();
    my $temp_filename = $temp_fh->filename();
    close($temp_fh);
    
    my $tempdir = tempdir(CLEANUP => 1);

    my $rnammer_path = $self->get_rnammer_path;
    confess "Could not determine rnammer path for version " . $self->{_version} unless defined $rnammer_path;

    ## For RNAmmer 1.2 at least, the default
    ## behaviour when -m is not specified
    ## is the same as -m tsu,lsu,ssu.
    my @cmd = (
        $rnammer_path,
        '-S',
        substr($domain, 0, 3),
        '-T',
        $tempdir,
        '-gff',
        $temp_filename,
    );

    eval {
        IPC::Run::run(
            \@cmd,
            IO::File->new($seq_fh->filename()),
            '>',
            \$rnammer_stdout,
            '2>',
            \$rnammer_stderr, 
      ) || die $CHILD_ERROR;
    };
    if ($EVAL_ERROR) {
        die "Failed to exec rnammer: $EVAL_ERROR";
    }

    my $gff = Bio::Tools::GFF->new(-file => $temp_filename, -gff_version => 1);
    while (my $feature = $gff->next_feature()) {
        $seq->add_SeqFeature($feature);
    }
        
}

sub get_rnammer_path {
    my $self = shift;
    my $version = $self->{_version};
    confess "No rnammer version provided!" unless defined $version;

    my $base_path = '/gsc/pkg/bio/rnammer/';
    my $rnammer_path = $base_path . "rnammer-$version/rnammer";

    unless (-e $rnammer_path) {
        confess "No rnammer executable found at expected location $rnammer_path";
    }

    return $rnammer_path;
}

sub genes {

    my ($self) = @_;


    return $self->{_genes};

}

sub _write_seqfile {

    my ($self, @seq) = @_;


    my $seq_fh = File::Temp->new();

    my $seqstream = Bio::SeqIO->new(
                                    -fh => $seq_fh,
                                    -format => 'Fasta',
                                );

    foreach my $seq (@seq) {
        $seqstream->write_seq($seq);
    }

    close($seq_fh);
    $seqstream->close();

    return $seq_fh;
    
}

1;
