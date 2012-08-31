#:eclark This report has a lot of business logic that could be refactorered for readability, the biggest immediately concern is the hardcoded path to human build36

package Genome::Model::ReferenceAlignment::Report::HqPos;

use strict;
use warnings;

use Genome;
use Genome::Info::IUB;
use Sort::Windowed;
use Sort::Merge;

class Genome::Model::ReferenceAlignment::Report::HqPos {
    is => 'Genome::Model::Command::BaseDeprecated',
    has => [
        min_alignment_quality   =>  { is => 'Integer', default_value => 15, is_optional => 1,
                                        doc => "minimum alignmnet quality", },
        min_base_quality        =>  { is => 'Integer', default_value => 20, is_optional => 1,
                                        doc => "minimum base quality at discrep. point" },
        min_read_count          =>  { is => 'Integer', default_value => 2,  is_optional => 1,
                                        doc => "minimum number of supporting variant reads" },
        chromosomes             =>  { is => 'List', is_optional => 1, default_value => join(',',1..22,'X','Y'), is_optional => 1,
                                        doc => "the list of reference sequences to check" },
        duplicate_weight        =>  { is => 'Decimal', default_value => .02,
                                        doc => 'the liklihood that a same start-site read is a distinct fragment NOT ENABLED' },
    ],
    doc => "scan alignments for positions at which there is hiqh quality variant evidence"
};

sub help_detail {
    return <<EOS;
genome-model report hq-pos tumor%v0b
genome-model report hq-pos tumor%v0b --chromosomes 22,X --min-alignment-quality=15 --min-base-quality=20
EOS
}

sub execute {
    $DB::single = $DB::stopper;
    my $self = shift;
    my $model = $self->model;
    my $chromosomes = $self->chromosomes;
    my @chromosomes = split(/,/,$chromosomes);
    my $n;
    for my $chrom (@chromosomes) {
        $self->status_message("processing reference sequence $chrom");
        my $i = $self->_refseq_position_summary_iterator($chrom);
        while (my $data = $i->()) {
            no warnings;
            $n++;
            print join("\t",@$data),"\n";
        }
    }
    $self->status_message("found $n positions with high-quality discrepancies");
    return 1;
}

sub _refseq_position_summary_iterator {
    my ($self,$chrom) = @_;
    my $ii = $self->_refseq_library_merge_position_iterator($chrom);
    my $min_read_count = $self->min_read_count;
    my @libraries = sort $self->model->libraries;
    my $next;
    my $i = sub {
        my $result;
        until ($result) {
            $next ||= $ii->();
            unless ($next) {
                return;
            }

            my @all_at_this_position = ($next);
            my ($chrom1,$p1,$r1,$a1,$c1,$g1,$t1,$detail,$lib) = @$next;
            while(1) {
                $next = $ii->();
                last unless ($next);
                my ($chrom2,$p2,$r2,$a2,$c2,$g2,$t2,$detail,$lib) = @$next;
                if (($chrom2 ne $chrom1) or ($p2 != $p1)) {
                    #stop collecting lines
                    #$next will be around for the next iteration
                    last;
                }
                else {
                    #the position and chromosome matches
                    push @all_at_this_position, $next;
                }
            }

            #assert: @all_at_this_position now contains every line dealing with a particular snp
            my %total_read_count;
            my %lib;
            for my $hqd (@all_at_this_position) {
                my ($chrom,$p,$r,$a,$c,$g,$t,$detail,$lib) = @$hqd;
                $total_read_count{A}+=$a;
                $total_read_count{C}+=$c;
                $total_read_count{G}+=$g;
                $total_read_count{T}+=$t;
                $lib{$lib} = [$a,$c,$g,$t];
            }

            # NOTE!!!: this call is really just filler since the DTR does the real call
            # If this call were used, we'd need a better algorithm (like the consensus caller)
            my $call;
            my $max_count = 0;
            for my $base (qw/A C G T/) {
                my $count = $total_read_count{$base} || 0;
                if ($count >= $max_count) {
                    $call = $base;
                    $max_count = $count;
                }
            }
            my @bcall_fields = (3,4,5,6);
            if ($max_count >= $min_read_count) {
                $result = [
                    $chrom1,$p1,$r1,Genome::Info::IUB->iub_for_alleles($r1,$call),
                    $total_read_count{A},
                    $total_read_count{C},
                    $total_read_count{G},
                    $total_read_count{T},

                ];
                for my $lib (@libraries) {
                    if ($lib{$lib}) {
                        push @$result, @{$lib{$lib}};
                    }
                    else {
                        push @$result, 0,0,0,0;
                    }
                }
            }
        }
        return $result;
    };
    return $i;
}

sub _refseq_library_merge_position_iterator {
    my $self = shift;
    my $refseq = shift;

    my @libraries = $self->model->libraries;
    my @inputs =
        map {
            $self->_refseq_library_position_iterator(
                $refseq,
                $_
            );
        } @libraries;

    my @next = map { $_->() } @inputs;
    my $i = sub {
        my $first_n = -1;
        my $first_p = -1;
        for (my $n=0; $n<@next; $n++) {
            next unless $next[$n];
            if ($next[$n][1] < $first_p or $first_n == -1) {
                $first_n = $n;
                $first_p = $next[$n][1];
                #print "first item $first_p in lib $first_n\n";
            }
        }
        return if $first_n == -1;
        my $this = $next[$first_n];
        $next[$first_n] = $inputs[$first_n]->();
        push @$this,$libraries[$first_n];
        return $this;
    };
    return $i;
}

sub _refseq_library_position_iterator {
    my ($self,$chrom,$lib) = @_;

    my ($ma_event) = $self->model->events("event_type like" => '%merge-ali%', ref_seq_id => $chrom);
    my $map = $ma_event->resolve_accumulated_alignments_filename(ref_seq_id => $chrom, library_name => $lib);

    $self->status_message("map for $chrom $lib is $map\n");

    my $mapview_iterator = $self->_mapview_iterator($map);
    unless ($mapview_iterator) {
        die "Failed to get a mapview iterator for $map!";
    }

    my $hq_discrepancy_iterator = $self->_extract_hq_discrepancies_iterator($mapview_iterator);

    my $sorted_hq_discrepancy_iterator = Sort::Windowed::new_iterator(
        sub {
            my $d = $hq_discrepancy_iterator->();
            return unless defined $d;
            return ($d->[1],$d);
        },
        100
    );

    my $library_position_iterator = $self->_library_position_summary_iterator_for_hq_discrepancy_iterator(
        $sorted_hq_discrepancy_iterator
    );

    return $library_position_iterator;
}

sub _library_position_summary_iterator_for_hq_discrepancy_iterator {
    my ($self,$sorted_hq_discrepancy_iterator) = @_;
    my $min_read_count = $self->min_read_count;
    my $next;
    my $i = sub {
        my $result;
        until ($result) {
            ##are we using LAST to store any previous thing, or just filling next and not doing anything with it?
            $next ||= $sorted_hq_discrepancy_iterator->();
            unless ($next) {
                return;
            }

            my @all_at_this_position = ($next);
            my ($c1,$p1,$r1,$v1,$aq1,$bp1,$rn1,$cc1,$ss1,$o1) = @$next;
            while(1) {
                $next = $sorted_hq_discrepancy_iterator->();
                last unless ($next);
                my ($c2,$p2,$r2,$v2,$aq2,$bp2,$rn2,$cc2,$ss2,$o2) = @$next;
                if (($c2 ne $c1) or ($p2 != $p1)) {
                    #this line doesn't concern the snp we are about to generate information about.
                    #stop collecting lines
                    last;
                }
                else {
                    #the position and chromosome matches
                    push @all_at_this_position, $next;
                }
            }

            #assert: @all_at_this_position now contains every line dealing with a particular snp
            my %genotype_count;
            for my $hqd (@all_at_this_position) {
                my ($c1,$p1,$r1,$v1,$aq1,$bp1,$rn1,$cc1,$ss1,$o1) = @$hqd;
                $genotype_count{$v1}{$ss1.$o1}++;
            }

            # todo: incorporate weighting of duplicates here
            my $a = scalar(keys(%{$genotype_count{A}}));
            my $c = scalar(keys(%{$genotype_count{C}}));
            my $g = scalar(keys(%{$genotype_count{G}}));
            my $t = scalar(keys(%{$genotype_count{T}}));

            $result = [
                $c1,$p1,$r1,
                $a,$c,$g,$t,
                \@all_at_this_position
            ];
        }
        return $result;
    };
    return $i;
}

sub _mapview_iterator {
    # return a mapview iterator for a map file
    my ($self,$map) = @_;
    my $fh = IO::File->new("maq mapview $map |");
    unless ($fh) {
        die "failed to execute maq mapview: $!";
    }
    my $nextline;
    my @cols;
    return
        sub {
            $nextline = $fh->getline;
            if (not defined $nextline) {
                $fh->close;
                return;
            }
            @cols = split(/\t/,$nextline);
            #print "mvi: @cols\n";
            return \@cols;
        };
}

sub _extract_hq_discrepancies_iterator {
    # return a hq discrepancy iterator for a mapview iterator
    my ($self,$mapview_iterator) = @_;
    my $minaq = $self->min_alignment_quality;
    my $minbq = $self->min_base_quality;
    my %rh;
    my @out;
    my $show;
    my $i = sub {
        if (@out) {
            ##print "hqd cache: @{$out[0]}\n";
            return shift @out;
        }

        my ($mapview_data,$c,$p,$aq);
        while($mapview_data = $mapview_iterator->()) {
            $c = $mapview_data->[1];
            $p = $mapview_data->[2];
            $aq = $mapview_data->[6];
            if($aq >= $minaq) {
                my $l = $mapview_data->[13];
                my $s = uc($mapview_data->[14]);
                my @bq = map { ord($_) } split(//,$mapview_data->[15]);
                my $rh = $rh{$c} ||= IO::File->new(Genome::Config::reference_sequence_directory() . '/NCBI-human-build36-flat/$c.bases');
                seek($rh,$p-1,0);
                my $rs;
                sysread($rh,$rs,$l*2,0);
                my $m = 0;
                for(my $n=0;$n<$l;$n++) {
                    if ($bq[$n]>=$minbq and (substr($s,$n,1) ne substr($rs,$n,1))) {
                        if ($show) {
                            print STDERR $rs,"\n" unless $m;
                            print STDERR ($s,"\n") unless $m;
                            $m = 1;
                            print STDERR ((" " x $n),"*","\n");
                        }
                        push @out, [$c,$p+$n,substr($rs,$n,1),substr($s,$n,1),$aq,$bq[$n],$n,@$mapview_data];
                    }
                }
            }
            last if @out;
        }
        return unless @out;
        ##print "hqd new: @{$out[0]}\n";
        return shift @out;
    };
    return $i;
}

1;
