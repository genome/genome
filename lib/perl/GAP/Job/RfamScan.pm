package GAP::Job::RfamScan;

use strict;
use warnings;

use GAP::Job;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Carp;
use English;
use File::Temp;
use IO::File;
use IPC::Run;
use Bio::AlignIO;

use base qw(GAP::Job);


sub new {

    my ($class, $seq, $job_id) = @_;

    
    my $self = { };
    bless $self, $class;
        

    unless (defined($job_id)) {
        croak 'missing job id';
    }

    $self->job_id($job_id);
    
    unless (defined($seq)) {
        croak 'missing seq object!';
    }

    unless ($seq->isa('Bio::PrimarySeqI')) {
        croak 'seq object is not a Bio::PrimaySeqI!';
    }
    
    $self->{_seq} = $seq;
    
    return $self;
    
}

sub execute {
    
    my ($self) = @_;

    $self->SUPER::execute(@_);
   
   
    my $seq = $self->{_seq};

    my $seq_fh = $self->_write_seqfile($seq);
    
    my ($rfam_stdout, $rfam_stderr);

    my $temp_fh       = File::Temp->new();
    my $temp_filename = $temp_fh->filename();

    close($temp_fh);
    
    my @cmd = (
               'rfam_scan',
               '-o',
               $temp_filename,
               $seq_fh->filename(),
           );
    eval {
        
        IPC::Run::run(
                      \@cmd,
                      undef,
                      '>',
                      \$rfam_stdout,
                      '2>',
                      \$rfam_stderr, 
                  ) || die $CHILD_ERROR;
        
    };

    if ($EVAL_ERROR) {
        die "Failed to exec rfam_scan: $EVAL_ERROR";
    }

    my $rfam_seed_file = '/gsc/pkg/bio/rfam/installed/Rfam.seed';

    unless (-e $rfam_seed_file ){
    
	die " Missing rfam seed file $rfam_seed_file! ";
    
    }

    my %rfam_basket = ( );
    
    my $io = Bio::AlignIO->new(-file => $rfam_seed_file, -format => 'stockholm' );
    
    while ( my $aln = $io->next_aln() ) {
	
	my $annotation = $aln->annotation();
	
	my ( $entry_type ) = $annotation->get_Annotations('entry_type');
	
	$rfam_basket{$aln->id()} = $aln->description(); 
	
    }
    
    my $rfam_fh = IO::File->new();
    $rfam_fh->open($temp_filename)
	or die "Can't open '$temp_filename': $OS_ERROR";
    
    while (my $line = <$rfam_fh>) {
        
	chomp $line;
        
	my (
	    $seq_id,
	    $seq_start,
	    $seq_end,
	    $rfam_acc,
	    $model_start,
	    $model_end,
	    $bit_score,
	    $rfam_id,
	    ) = split /\s+/, $line;
	
	my $seq_strand = $seq_start > $seq_end ? -1 : 1;

#some rfam id's (RF00079, BACEGGDFTB_Contig301.1 ) are coming back with a semicolon on the end messing up pattern matching.

	if ( $rfam_id =~ /\;/ ) {

	    $rfam_id =~ s/\;//;
	
	}
	
	my $rfam_prod = $rfam_basket{$rfam_id};
	
	if ($rfam_prod =~ /Bacterial RNase P class A/i){ #genbank wanted us to replace (bacterial RNase P class A) with (RNA component of RNase P);06/15/07 per Veena Bhonagiri)
	    
	    $rfam_prod = "RNA component of RNase P";
	    
	}
	
	my $feature = Bio::SeqFeature::Generic->new(
						    -seq_id => $seq_id,
						    -start  => $seq_start,
						    -end    => $seq_end,
						    -strand => $seq_strand,
						    -source => 'Infernal',
						    -score  => $bit_score,
						    -tag => {
							acc       => $rfam_acc,
							id        => $rfam_id,
							rfam_prod => $rfam_prod,
						    },
						    );
	
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
