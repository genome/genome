package BAP::GeneMerge;

use strict;
use warnings;

use Carp;
use Graph;


sub find_overlaps {
    
    my ($callback, $feature_ref) = @_;
    
    
    unless (ref($callback) eq 'CODE') {
        croak 'callback arg is not a CODE ref';
    }

    foreach my $feature (@{$feature_ref}) {
        
        unless ($feature->isa('Bio::SeqFeatureI')) {
            croak 'feature does not implement Bio::SeqFeatureI';
        }
        
    }
    
    unless (@{$feature_ref} > 1) {
        return;
    }

    my @sorted_features = sort {
                                $a->start() <=> $b->start()
                                    ||
                                $a->length() <=> $b->length()
                               } @{$feature_ref};
    
  OUTER: foreach my $o (0..$#sorted_features) {
        
        my $f = $sorted_features[$o];
        
        if ($o < $#sorted_features) {
            
          INNER1: foreach my $i ($o+1..$#sorted_features) {
                
                my $g = $sorted_features[$i];
                
                if ($g->start() > $f->end()) {
                    last INNER1;
                }
                else {
                    $callback->($f, $g);
                }
                
            }
            
        }
        
        if ($o > 0) {
            
          INNER2: foreach my $i (reverse(0..$o-1)) {
                
                my $g = $sorted_features[$i];
                
                if ($g->end() < $f->start()) {
                    last INNER2;
                }
                else {
                    $callback->($g, $f);
                }
                
            }
            
        }
        
    }
    
}

sub graph_overlaps {
    
    my (
        $overlap_threshold_bases,
        $overlap_threshold_percent,
        $feature_ref
    ) = @_;
    
    
    my $graph = Graph::Undirected->new();
    
    my $overlap = sub {
        
        my ($f, $g) = @_;

        
        my $overlap_count;
        
        if ($g->end() < $f->end()) {
            $overlap_count = $g->length();
        }
        else {
            $overlap_count = abs($f->end() - $g->start()) + 1;
        }
        
        my $overlap_percent_f = sprintf("%.1f", ($overlap_count / $f->length) * 100);
        my $overlap_percent_g = sprintf("%.1f", ($overlap_count / $g->length) * 100);
        
        if (
            ($overlap_count     > $overlap_threshold_bases) ||
            ($overlap_percent_f > $overlap_threshold_percent) ||
            ($overlap_percent_g > $overlap_threshold_percent)
        ) {

            $graph->add_edge($f->display_name(), $g->display_name());
        
        }
        
    };
    
    find_overlaps($overlap, $feature_ref);

    return $graph;
    
}

sub better_gene {

    my ($gene_a, $gene_b) = @_;

    
    ## pfam (iprscan) evidence > blastp evidence
    ##
    ## Truth table for best evidence for two genes
    ##
    ## P = pfam (iprscan) evidence
    ## B = blastp evidence
    ## - = no evidence
    ##
    ## Gene A | Gene B |  Result
    ## ---------------------------
    ##    P        P    keep both
    ##    P        B    toss B
    ##    B        P    toss A
    ##    B        B    keep longest
    ##    -        -    keep longest
    ##    B        -    toss B
    ##    -        B    toss A

    if ($gene_a->has_tag('pfam_evidence')) {
        return 0;
    }
    elsif ($gene_b->has_tag('pfam_evidence')) {
        return 1;
    }
    elsif ($gene_b->has_tag('blastp_evidence')) {

        if ($gene_a->has_tag('blastp_evidence')) {

            if ($gene_b->length() > $gene_a->length()) {
                return 1;
            }
            else {
                return 0;
            }
            
        }
        else {
            return 1;
        }
    }
    elsif ($gene_a->has_tag('blastp_evidence')) {
        return 0;
    }
    else {
        
        if ($gene_b->length() > $gene_a->length()) {
            return 1;
        }
        else {
            return 0;
        }
        
    }
    
}

sub fancy_tag_overlapping {
    
    my ($graph, $feature_ref) = @_;
    
    
    my %feature_index = map { $_->display_name() => $_ } @{$feature_ref};
    
    foreach my $vertex ($graph->vertices()) {
        
        my $vertex_feature = $feature_index{$vertex};
        
        my $better_neighbor;
        
        if ($better_neighbor = better_neighbor($graph, $vertex, \%feature_index)) {

            $vertex_feature->add_tag_value('delete_overlap' => 1);
            $vertex_feature->add_tag_value('delete_for' => $better_neighbor);
            $graph->delete_vertex($vertex);
            
        }

    }
    
}

sub better_neighbor {
    
    my ($graph, $vertex, $feature_ref) = @_;
    
    
    my $vertex_feature  = $feature_ref->{$vertex};
    
    foreach my $neighbor ($graph->neighbors($vertex)) {

        my $neighbor_feature = $feature_ref->{$neighbor};
        
        if (better_gene($vertex_feature, $neighbor_feature)) {
            unless (better_neighbor($graph, $neighbor, $feature_ref)) {
                return $neighbor;
            }
            
        }
        
    }
    
    return undef;
    
}

sub tag_redundant_rfam {

    my ($feature_ref) = @_;


    my $overlap = sub {

        my ($rfam)    = grep { $_->source_tag() eq 'rfam'    } @_;
        my ($rnammer) = grep { $_->source_tag() eq 'rnammer' } @_;
        
        if (
            defined($rfam) && defined($rnammer)
        ) {
            unless ($rfam->has_tag('redundant')) {
                $rfam->add_tag_value('redundant', 1);
            }
        }
        
    };
    
    find_overlaps($overlap, $feature_ref);
    
}

sub tag_rna_overlap {
    
    my ($feature_ref) = @_;
    
    
    my $overlap = sub {
        
        my ($f, $g) = @_;
        
        my ($f_gene_type) = $f->each_tag_value('type');
        my ($g_gene_type) = $g->each_tag_value('type');
        
        my $overlap_count;
        
        if ($g->end() < $f->end()) {
            $overlap_count = $g->length();
        }
        else {
            $overlap_count = abs($f->end() - $g->start()) + 1;
        }
        
        my $overlap_percent_f = sprintf("%.1f", ($overlap_count / $f->length) * 100);
        my $overlap_percent_g = sprintf("%.1f", ($overlap_count / $g->length) * 100);
        
        foreach my $h ($f, $g) {
            
            my ($h_gene_type) = $h->each_tag_value('type');
            
            if ($h_gene_type eq 'rRNA') {
                
                if ($h->source_tag() eq 'rfam') {
                    
                    my $threshold = 10;
                    
                    if ($h->has_tag('overlap_50')) {
                        $threshold = 50;
                    }
                    
                    if (
                        ($overlap_percent_f > $threshold)
                        ||
                        ($overlap_percent_g > $threshold)
                    ) {
                        
                        if ($f_gene_type eq 'coding') {
                            unless ($f->has_tag('delete_rrna_overlap')) {
                                $f->add_tag_value('delete_rrna_overlap', 1);
                            }
                        }
                        
                        if ($g_gene_type eq 'coding') {
                            unless ($g->has_tag('delete_rrna_overlap')) {
                                $g->add_tag_value('delete_rrna_overlap', 1);
                            }
                        }
                        
                    }
                    
                }
                elsif ($h->source_tag() eq 'rnammer') {
                    
                    if ($f_gene_type eq 'coding') {
                        unless ($f->has_tag('delete_rrna_overlap')) {
                            $f->add_tag_value('delete_rrna_overlap', 1);
                        }
                    }
                    
                    if ($g_gene_type eq 'coding') {
                        unless ($g->has_tag('delete_rrna_overlap')) {
                            $g->add_tag_value('delete_rrna_overlap', 1);
                        }
                    }
                    
                }
                
            }
            elsif ($h_gene_type eq 'tRNA') {
                
                if ($f->strand() eq $g->strand()) {
                    
                    if ($f_gene_type eq 'coding') {
                        unless ($f->has_tag('delete_trna_overlap')) {
                            $f->add_tag_value('delete_trna_overlap', 1);
                        }
                    }
                    
                    if ($g_gene_type eq 'coding') {
                        unless ($g->has_tag('delete_trna_overlap')) {
                                $g->add_tag_value('delete_trna_overlap', 1);
                            }
                    }
                    
                }
                else {
                    
                    if ($f_gene_type eq 'coding') {
                        unless ($f->has_tag('check_trna_overlap')) {
                            $f->add_tag_value('check_trna_overlap', 1);
                        }
                    }
                    
                    if ($g_gene_type eq 'coding') {
                        unless ($g->has_tag('check_trna_overlap')) {
                            $g->add_tag_value('check_trna_overlap', 1);
                        }
                    }
                    
                }
                
            }
            
        }
        
    };
    
    find_overlaps(
                  $overlap,
                  $feature_ref,
              );
    
}

1;
