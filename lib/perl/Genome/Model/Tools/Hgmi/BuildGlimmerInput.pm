package Genome::Model::Tools::Hgmi::BuildGlimmerInput;

use strict;
use warnings;

use Genome;
use lib "/gsc/scripts/opt/bacterial-bioperl";

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Glimmer;
use Bio::Tools::Run::Glimmer;

use English;
use IO::File;
use IPC::Run;

class Genome::Model::Tools::Hgmi::BuildGlimmerInput {
    is => 'Command',
    has => [
        fasta_files => { is => 'ARRAY',
                         doc => 'array of fasta file names',
                         is_input => 1, },
        model_file => { is => 'SCALAR',
                        is_optional => 1,
                        doc => 'absolute path to the model file for this fasta',
                        is_output => 1, },
        pwm_file => { is => 'SCALAR' ,
                      is_optional => 1,
                      doc => 'absolute path to the pwm file for this fasta',
                      is_output => 1, }
    ],
};

sub help_brief {
    "Write a set of fasta files for an assembly";
}

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}

sub execute {

    my $self = shift;

    ##FIXME: This should not be hardcoded, the necessary binaries should be deployed to /gsc.
    my $glimmer_bin_dir = $ENV{GENOME_SW} . '/glimmer/glimmer3.02/bin';

    ##FIXME: This should not be hardcoded either, should be in a superclass or it should be a
    ##       parameter.
    my $network_temp_dir = '/gscmnt/temp212/info/annotation/MGAP_tmp';

    my %network_temp_params = (
                               'DIR'     => $network_temp_dir,
                               'SUFFIX'  => 'g3.icm.tmp',
                               'TEMPLATE'=> 'MGAP_XXXXXXXX',
                               'UNLINK'  => 0,
                           );

    ## Semi-ugly hack to handle multi-fasta input outside
    ## of a Workflow
    my @input_fasta_files = @{$self->fasta_files()};
    my @fasta_handles     = ( );
    my @fasta_files       = ( );

    foreach my $fasta_file (@input_fasta_files) {

        my $seqin = Bio::SeqIO->new(-file => $fasta_file, -format => 'Fasta');

        while (my $seq = $seqin->next_seq()) {
        
            my $tmp_fh = File::Temp->new();
            my $tmp_fn = $tmp_fh->filename();

            $tmp_fh->close();
       
            my $seqout = Bio::SeqIO->new(-file => ">$tmp_fn", -format => 'Fasta');
            
            $seqout->write_seq($seq);

            push @fasta_handles, $tmp_fh;
            push @fasta_files,   $tmp_fn;
            
        }
    
    }
    
    my $train_fh    = File::Temp->new();
    my $train_seqio = Bio::SeqIO->new(-fh => $train_fh, -format => 'Fasta');
    
    foreach my $fasta_file (@fasta_files) {

        my $contig_seq = Bio::SeqIO->new(-file => $fasta_file, -format => 'Fasta')->next_seq();
        
        my $lo_stdout = File::Temp->new();
        
        my $lo_stderr;
        
        my @long_orfs_command = (
                                 "$glimmer_bin_dir/long-orfs",
                                 '-n',
                                 '-t',
                                 '1.15',
                                 "$fasta_file",
                                 '-',
                             );
        
        IPC::Run::run(
                      \@long_orfs_command, 
                      \undef, 
                      '>', 
                      $lo_stdout->filename(), 
                      '2>', 
                      \$lo_stderr, 
                     ); 

        my $lo_fh = IO::File->new();
        $lo_fh->open($lo_stdout->filename()) or die "Can't open '$lo_stdout': $OS_ERROR";

        while (my $orf = <$lo_fh>) {
            
            chomp $orf;
            
            my ($orf_no, $start, $end, $frame, $score) = split /\s+/, $orf;
            my $strand = substr($frame, 0, 1) eq '+' ? '1' : '-1';
            
            if ($start > $end) { ($start, $end) = ($end, $start); }
            
            my $orf_seq = $contig_seq->trunc($start, $end);
            
            unless ($strand > 0) { $orf_seq = $orf_seq->revcom(); }
            
            $orf_seq->description("from $start to $end on strand $strand");
            
            $train_seqio->write_seq($orf_seq);
            
        }
        
        $lo_fh->close();
        
    }
    
    $train_seqio->close();
    $train_fh->close();
    
    my $icm_fh = File::Temp->new(%network_temp_params);
                                
    $icm_fh->close();
    
    ## Use the training file as input to build-icm to create an icm file

    my ($bi_stdout, $bi_stderr);
    
    my @cmd = (
               "$glimmer_bin_dir/build-icm",
               '-r',
               $icm_fh->filename(),
           );

    IPC::Run::run(\@cmd, '<', $train_fh->filename(),'>', \$bi_stdout, '2>', \$bi_stderr);

    $self->model_file($icm_fh->filename());
    
    ## Run glimmer3 on each contig, extract sequence in front of the 5' end of
    ## each predicted gene, output to the upstream file
    my $upstream_fh    = File::Temp->new();
    my $upstream_seqio = Bio::SeqIO->new(-fh => $upstream_fh, -format => 'Fasta');
    
    foreach my $fasta_file (@fasta_files) {

        my $contig_seq = Bio::SeqIO->new(-file => $fasta_file, -format => 'Fasta')->next_seq();
    
         my $factory = Bio::Tools::Run::Glimmer->new(
                                                     -program => 'glimmer3',
                                                     -model   => $icm_fh->filename(),
                                                     -o       => 50,
                                                     -g       => 110,
                                                     -t       => 30,
                                                 );
        
        my $glimmer = $factory->run($contig_seq);

        GENE: while (my $gene = $glimmer->next_prediction()) {

             my $location = $gene->location();

             ## For now, keep it simple and only deal with
             ## predictions with simple locations.  We
             ## definately don't want fragments, if they
             ## don't have a good start, we'll screw up the
             ## PWM for the Ribosomal Binding Sites by building
             ## it on the wrong sequence.  Wraparound genes
             ## would be okay, but there aren't going to be that
             ## many, at most one per contig, if the contig is
             ## even linear.  For now skip 'em.  Should
             ## refactor this to be able to handle wraparound
             ## genes without making the logic a mess
             unless ($location->isa('Bio::Location::Simple')) {
                 next GENE;
             }

             my $start  = $gene->start();
             my $end    = $gene->end();
             my $strand = $gene->strand();
             
             if ($strand > 0) {
                 ($start, $end) = ($start - 25, $start - 1);
             }
             else {
                 ($start, $end) = ($start + 1, $start + 25);
             }
             
             foreach my $coord ($start, $end) {
                 
                 if ($coord < 1) {
                     next GENE;
                 }
                 
                 if ($coord > $contig_seq->length()) {
                     next GENE;
                 }
                 
             }
             
             my $upstream_seq = $contig_seq->trunc($start, $end);

             if ($strand < 1) {
                 $upstream_seq = $upstream_seq->revcom()
             }
             
             $upstream_seq->description("from $start to $end on strand $strand");
             $upstream_seqio->write_seq($upstream_seq);
             
         }

    }

    $upstream_fh->close();
    $upstream_seqio->close();

    ## Run elph on the upstream file, extract the motif count section into
    ## motif_file
    my $motif_fh = File::Temp->new(%network_temp_params);

    print $motif_fh "6", "\n";

    my $found_motif_counts = 0;
    my $motif_line_count   = 0;
    
    my $elph_stdout = File::Temp->new();
    my $elph_stderr;

    my @elph_command = (
                        'elph',
                        $upstream_fh->filename(),
                        'LEN=6',
                    );

    IPC::Run::run(\@elph_command, \undef, '>', $elph_stdout, '2>', \$elph_stderr);

    my $elph_fh = IO::File->new();
    $elph_fh->open($elph_stdout) or die "Can't open '$elph_stdout': $OS_ERROR";

    while (my $line = <$elph_fh>) {
        
        chomp $line;
        
        if ($line =~ /^Motif counts:/) {
            $found_motif_counts = 1;
            next;
        }
        
        unless ($found_motif_counts) { next; }
        
        if ($line =~ /^[acgt]:/) {
            
            $motif_line_count++;
            $line =~ s/://g;
            
            my @cols = split /\s+/, $line;
            
            print $motif_fh join("\t", @cols), "\n";
            
            if ($motif_line_count == 4) { last; }
            
        }
        
    }

    $motif_fh->close();

    $self->pwm_file($motif_fh->filename());
    return 1;
    
}
 
1;
