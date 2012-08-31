package Genome::Model::Tools::BioSamtools::Breakdown;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::Breakdown {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        bam_file => { doc => 'BAM format file(s) of alignment data' },
        output_file => { doc => 'The tsv output file' },
    ],
    has_optional => [
        _tsv_fh => { },
        _annotation => { },
        _total => {
            default_value => 0,
        },
    ],
};

my @TAGS = qw/TRANSCRIPT CHEMISTRY MITOCHONDRIAL RIBOSOME-PROTEIN RIBOSOME/;

sub create {
    my $class = shift;
    my %params = @_;
    my $bam_file = delete($params{bam_file});
    my $self = $class->SUPER::create(%params);
    unless ($self) { return; }
    $self->bam_file($bam_file);
    return $self;
}

sub execute {
    my $self = shift;
    my $tsv_fh = IO::File->new($self->output_file,'w');
    unless ($tsv_fh) {
        die('Failed to create output tsv file handle for file '. $self->output_file);
    }
    $self->_tsv_fh($tsv_fh);
    my $bam_files = $self->bam_file;

    my @bam_files;
    if (ref($bam_files) eq 'ARRAY') {
        @bam_files = @{$bam_files};
    } else {
        push @bam_files, $bam_files;
    }
    my $annotation = {};
    my $total = 0;
    for my $bam_file (@bam_files) {
        my $bam = Bio::DB::Bam->open($bam_file);
        my $bai_file = $bam_file .'.bai';
        if (-e $bai_file) {
            my $bam_mtime = (stat($bam_file))[9];
            my $bai_mtime = (stat($bai_file))[9];
            if ($bam_mtime > $bai_mtime) {
                unless (unlink $bai_file) {
                    die('Failed to remove old bai file'. $bai_file); 
                }
            }
        }
        unless (-e $bai_file) {
            my $index_cmd = "samtools index $bam_file";
            my $index_rv = system($index_cmd);
            unless ($index_rv == 0) {
                die('non-zero return value '. $index_rv .' from command: '. $index_cmd .":  $!");
            }
        }

        my $index = Bio::DB::Bam->index_open($bam_file);
        unless ($index) {
            die('Failed to find index for bam '. $bam_file);
        }
        my $header = $bam->header;
        my $n_targets = $header->n_targets;
        my $target_name = $header->target_name;
        my $target_len = $header->target_len;
        my $callback = sub {
            my $alignment = shift;
            my $found;
            my $tid = $alignment->tid;
            my $target = $target_name->[$tid];
            my $length = $target_len->[$tid];
            for my $tag (@TAGS) {
                if ($target =~ /^$tag/) {
                    $self->_update_totals( $alignment->query->dna, $tag );
                    $found = 1;
                }
            }
            unless($found) {
                $self->_update_totals( $alignment->query->dna, 'UNKNOWN' )
            }
            $total++;
        };
        for (my $tid = 0; $tid < $n_targets; $tid++) {
            my $length = $target_len->[$tid];
            $index->fetch($bam,$tid,0,$length-1,$callback);
        }
    }
    
    $self->_total($total);
    for my $tag (@TAGS) {
        $self->_print_category_report( $tag );
    }
    $self->_print_category_report( 'UNKNOWN' );
    $self->_tsv_fh->close;
    return 1;
}


sub _update_totals {
    my ($self, $seq, $type) = @_;
    my $annotation = $self->_annotation;
    $annotation->{ $type }->{total}++;
    my ($polyAT_count, $polyA_count, $polyT_count) = &_assess_homopolymer_runs( $seq );
    if ($polyAT_count > 0) { $annotation->{ $type }->{polyAT}++ }
    if ($polyA_count  > 0) { $annotation->{ $type }->{polyA}++  }
    if ($polyT_count  > 0) { $annotation->{ $type }->{polyT}++  }
    $self->_annotation($annotation);
    return 1;
}

sub _assess_homopolymer_runs {
    # Loose definition is a run of 10 or more A's or T's for "long" homopolymer
    # runs.
    my $seq = shift;
    map {$_ = 0 } my ($polyAT_count, $polyA_count, $polyT_count);
    if ($seq =~ /TTTTTTTTTT/i || $seq =~ /AAAAAAAAAA/i) {
        if ($seq =~ /TTTTTTTTTT/i && $seq =~ /AAAAAAAAAA/i) {
            $polyAT_count++;
        }
        elsif ($seq =~ /AAAAAAAAAA/i) {
            $polyA_count++;
        }
        elsif ($seq =~ /TTTTTTTTTT/i) {
            $polyT_count++;
        }
    }
    return ($polyAT_count, $polyA_count, $polyT_count);
}

sub _print_category_report {
    my $self = shift;
    my $category = shift;
    my $annotation = $self->_annotation;
    my ($percent_category_reads, $percent_polyAT, $percent_polyA, $percent_polyT);

    # Means that the category is un-populated, so set it to 0 strictly
    # for calculation purposes--i.e., won't throw an error.
    if (!$annotation->{ $category }->{total}) { $annotation->{ $category }->{total} = 0 }

    # Calculations:
    $percent_category_reads = (($annotation->{ $category }->{total}  / $self->_total) * 100);

    if ($annotation->{$category }->{polyAT}) {
	$percent_polyAT = ($annotation->{$category }->{polyAT} == 0) ? 0 : (($annotation->{ $category }->{polyAT} / $annotation->{ $category }->{total}) * 100);
    }
    else {
	$percent_polyAT = 0;
	$annotation->{$category }->{polyAT} = 0;
    }

    if ($annotation->{$category }->{polyA}) {
	$percent_polyA  = ($annotation->{$category }->{polyA}  == 0) ? 0 : (($annotation->{ $category }->{polyA}  / $annotation->{ $category }->{total}) * 100);
    }
    else {
	$percent_polyA = 0;
	$annotation->{$category }->{polyA} = 0;
    }

    if ($annotation->{$category }->{polyT}) {
	$percent_polyT = ($annotation->{$category }->{polyT}  == 0) ? 0 : (($annotation->{ $category }->{polyT}  / $annotation->{ $category }->{total}) * 100);
    }
    else {
	$percent_polyT = 0;
	$annotation->{$category }->{polyT} = 0;
    }

    # TAB-DELIMITED FIELDS:
    # [0]  category
    # [1]  # reads in category
    # [2]  # total reads
    # [3]  % category reads in total reads
    # [4]  # polyAT reads in category
    # [5]  % polyAT reads in category
    # [6]  # polyA  reads in category
    # [7]  % polyA  reads in category
    # [8]  # polyT  reads in category
    # [9]  % polyT  reads in category
    my $tsv_fh = $self->_tsv_fh;
    unless ($annotation->{ $category }->{total} == 0) {
	print $tsv_fh join (
            "\t",
            $category,
            $annotation->{ $category }->{total},
            $self->_total,
            sprintf( "%.2f", $percent_category_reads ),
            $annotation->{ $category }->{polyAT},
            sprintf( "%.2f", $percent_polyAT ),
            $annotation->{ $category }->{polyA},
            sprintf( "%.2f", $percent_polyA ),
            $annotation->{ $category }->{polyT},
            sprintf( "%.2f", $percent_polyT ),
        ) . "\n";
    }
    $self->_annotation($annotation);
    return 1;
}

__END__
1;
