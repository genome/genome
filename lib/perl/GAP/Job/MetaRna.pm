package GAP::Job::MetaRna;

use strict;
use warnings;

use GAP::Job;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Carp;
use English;
use File::Temp qw/tempdir/;
use IO::File;
use IPC::Run;

use base qw(GAP::Job);


sub new {

    my ($class, $seq, $domain, $job_id) = @_;

    
    my $self = { };
    bless $self, $class;
    
    unless (defined($seq)) {
        croak 'missing seq object!';
    }
    
    unless ($seq->isa('Bio::PrimarySeqI')) {
        croak 'seq object is not a Bio::PrimaySeqI!';
    }
    
    $self->{_seq} = $seq;

    unless (defined($domain)) {
        croak 'missing domain';
    }

    unless (
            ($domain eq 'archaea'  )
            ||
            ($domain eq 'bacteria' )
            ||
            ($domain eq 'eukaryota')
        ) {
        croak "invalid domain '$domain'";
    }

    $self->{_domain} = $domain;
    
    unless (defined($job_id)) {
        croak 'missing job id';
    }

    $self->job_id($job_id);
    
    return $self;
    
}

sub execute {
    
    my ($self) = @_;

    $self->SUPER::execute(@_);
   
   
    my $seq    = $self->{_seq};
    my $domain = $self->{_domain};
    my $seq_fh = $self->_write_seqfile($seq);
    
    my ($metarna_stdout, $metarna_stderr);

    my $temp_fh       = File::Temp->new();
    my $temp_filename = $temp_fh->filename();
    
    close($temp_fh);
    
    my $tempdir = tempdir(CLEANUP => 1);

    # may need to set the HMM dir if we are running outside of the dir...
    my @cmd = (
            '/gsc/pkg/bio/meta_rna/rRNA_hmm_fs/rna_hmm.py',
            '-i',
            $seq_fh->filename(),
            '-o',
            $temp_filename,
            '-L',
            '/gsc/pkg/bio/meta_rna/rRNA_hmm_fs/HMMs/',
            '-e',
            '10e-5',
        );
    eval {
        
        IPC::Run::run(
                      \@cmd,
                      IO::File->new($seq_fh->filename()),
                      '>',
                      \$metarna_stdout,
                      '2>',
                      \$metarna_stderr, 
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
