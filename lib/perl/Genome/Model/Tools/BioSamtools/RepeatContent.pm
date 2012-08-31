package Genome::Model::Tools::BioSamtools::RepeatContent;

use strict;
use warnings;

use Genome;
use Bio::Tools::RepeatMasker;

my $DEFAULT_REPEATS_FILE = '/gscmnt/sata132/techd/solexa/jwalker/RepeatReducedCapture/RepeatMasker/hg18.fa.out';
class Genome::Model::Tools::BioSamtools::RepeatContent {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        bam_file => { doc => 'A BAM format file of alignment data' },
        output_file => { doc => 'A tab delimited RepeatMasker style table' },
        repeats_file => {
            doc => 'A RepeatMasker table of repeats for the genome.',
            is => 'Text',
            default_value => $DEFAULT_REPEATS_FILE,
            is_optional => 1,
        }
    ],
};

sub execute {
    my $self = shift;
    my $output_fh = IO::File->new($self->output_file,'w');
    unless ($output_fh) {
        die('Failed to open output file '. $self->output_file);
    }
    my $bam  = Bio::DB::Bam->open( $self->bam_file );
    unless ($bam) {
        die('Failed to load bam file '. $self->bam_file);
    }
    my $bai_file = $self->bam_file .'.bai';
    my @symlinks;
    if (-e $bai_file) {
        while (-l $bai_file) {
            push @symlinks, $bai_file;
            $bai_file = readlink($bai_file);
        }
        my $bam_mtime = (stat($self->bam_file))[9];
        my $bai_mtime = (stat($bai_file))[9];
        if ($bam_mtime > $bai_mtime) {
            unless (unlink $bai_file) {
                die('Failed to remove old bai file'. $bai_file);
            }
        }
    }
    unless (-e $bai_file) {
        Bio::DB::Bam->index_build($self->bam_file);
        if (@symlinks) {
            @symlinks = reverse(@symlinks);
            my $to_file = $bai_file;
            for my $from_file (@symlinks) {
                unless (symlink($to_file,$from_file)) {
                    die('Failed to create symlink '. $from_file .' => '. $to_file);
                }
                $to_file = $from_file;
            }
        }
    }

    unless (-e $bai_file) {
        die('Failed to create BAM index file!');
    }

    my $index  = Bio::DB::Bam->index_open( $self->bam_file );
    unless ($index) {
        die('Failed to find index for bam '. $self->bam_file);
    }

    my $header = $bam->header;
    my $n_targets = $header->n_targets;
    my $target_name = $header->target_name;
    my $target_len = $header->target_len;
    
    my %target_name_index;
    my $i = 0;
    for my $target_name (@{ $target_name }) {
        $target_name_index{$target_name} = $i++;
    }
    
    unless ($n_targets == $i) {
        die 'Expected '. $n_targets .' targets but counted '. $i .' indices';
    }
    print 'Start loop of all reads: '. `date` ."\n";
    my $aligned_reads = 0;
    my $total_count = 0;
    my $total_bp = 0;
    while (my $align = $bam->read1()) {
        $total_count++;
        $total_bp += $align->l_qseq;
        unless ($align->unmapped) {
            $aligned_reads++;
        }
    }
    my $parser = Bio::Tools::RepeatMasker->new(-file=>$self->repeats_file);
    my %repeats;
    my $repeat_aligned_reads = 0;
    print 'Start loop of all repeat regions: '. `date` ."\n";
    while (my $result = $parser->next_result) {
        my $tag = $result->primary_tag;
        my ($family,$class) = split("/",$tag);
        if ($family =~ /RNA/) {
            $family = 'Small RNA';
        }
        my $target = $result->seq_id;
        my $start = $result->start;
        my $end = $result->end;
        $target =~ s/chr//;
        # Here we get the $tid from the $gene_name or $seq_id
        my $tid = $target_name_index{$target};
        unless (defined $tid) { warn('Failed to get tid for target '. $target); next; }
        
        my $coverage = $index->coverage( $bam, $tid, $start - 1, $end );
        for my $depth (@{$coverage}) {
            $repeats{$family}{base_pair} += $depth;
            if ($class) {
                $repeats{$family}{$class}{base_pair} += $depth;
            }
        }
        my $callback = sub {
            my $align = shift;
            $repeat_aligned_reads++;
            $repeats{$family}{elements}++;
            if ($class) {
                $repeats{$family}{$class}{elements}++;
            }
        };
        $index->fetch( $bam, $tid, $start - 1, $end, $callback );;
    }
    print 'End loop of all repeat regions: '. `date` ."\n";
    my $masked_bp = 0;
    my $string = '';
    for my $family (sort keys %repeats) {
        my $family_bp = delete($repeats{$family}{base_pair}) || 0;
        $masked_bp += $family_bp;
        my $family_elements = delete($repeats{$family}{elements}) || 0;
        my $family_pc = sprintf("%.02f",(($family_bp / $total_bp ) * 100)) .'%';
        $string .= $family .":\t". $family_elements ."\t". $family_bp ."\t". $family_pc."\n";
        for my $class (sort keys %{$repeats{$family}}) {
            my $class_elements = $repeats{$family}{$class}{elements} || 0;
            my $class_bp = $repeats{$family}{$class}{base_pair} || 0;
            my $class_pc = sprintf("%.02f",(($class_bp / $total_bp ) * 100)) .'%';
            if ($class) {
                $string .= "\t". $class .":\t". $class_elements ."\t". $class_bp ."\t". $class_pc."\n";
            }
        }
        $string .= "\n";
    }
    
    print $output_fh "sequences:\t". $total_count ."\n";
    print $output_fh "total length:\t". $total_bp ."\n";
    print $output_fh "aligned:\t". $aligned_reads ."\n";
    print $output_fh "repeat aligned:\t". $repeat_aligned_reads ."\n";
    print $output_fh "masked:\t". $masked_bp ." bp ( ". sprintf("%02f",(($masked_bp / $total_bp ) * 100)) ." %) \n";
    print $output_fh $string ."\n";
    $output_fh->close;
    return 1;
}

1;
