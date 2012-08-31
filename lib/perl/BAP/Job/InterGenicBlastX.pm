package BAP::Job::InterGenicBlastX;

use strict;
use warnings;

use Bio::Tools::Run::StandAloneBlast;
use Carp;
use Sys::Hostname;
use base qw(GAP::Job);

sub new {
    my ($class, $job_id, $seq, $mask_ref, $db, $core_num, $expansion, $mask_char, $bit_score) = @_;
    
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

    unless (defined($mask_ref)) {
        croak 'missing mask ref!';
    }
    unless (ref($mask_ref) eq 'ARRAY') {
        croak 'mask ref is not an array ref!'
    }
    $self->{_masks} = $mask_ref;
    
    unless (defined($db)) {
        croak 'missing db!';
    }
    $self->{_db} = $db;

    unless (defined($core_num)) {
	    croak 'missing number of cores to run blast in Job!';
    }
    $self->{_core_num} = $core_num;
    
    unless (defined($expansion)) {
        $expansion = 300;
    }
    $self->{_expansion} = $expansion;
    
    unless (defined($mask_char)) {
        $mask_char = 'N';
    }

    $self->{_mask_char} = $mask_char;

    unless (defined($bit_score)) {
        $bit_score = 60;
    }

    $self->{_bit_score} = $bit_score;
    
    $self->{_genes} = [ ];
    
    return $self;
    
}

sub execute {
    
    my ($self) = @_;
 
    $self->SUPER::execute(@_);
    
    $self->_mask_seq();

    my $core_num = $self->{_core_num};


    my $factory = Bio::Tools::Run::StandAloneBlast->new(
                                                        -database => $self->{_db},
                                                        -expect   => '1e-6',
                                                        -p        => 'blastx',
                                                        -a        => $core_num,
                                                        -b        => 1500,
                                                        -v        => 1500,
                                                    );

    my $report = $factory->blastall($self->{_masked_seq});

    my $result     = $report->next_result();

    unless (defined($result)) { return; }

    my $query_name = $result->query_name();
    
    my %seen_coords = ( );
    
    while (my $hit = $result->next_hit()) {
        
        while (my $hsp = $hit->next_hsp()) {
            
            if ($hsp->bits() >= $self->{_bit_score}) {
                
                my $coverage = (($hsp->length('hsp') / $hsp->length('query')) * 100);
                
                if (
                    ($coverage >= 30) &&
                    ($hsp->percent_identity() >= 30)
                ) {
                    
                    my $start  = $hsp->start('query');
                    my $end    = $hsp->end('query');
                    
                    if ($start > $end) { ($start, $end) = ($end, $start) }
                    
                    ## No point in processing more than one alignment to the same
                    ## query coordinates, and it's cheapest to discard them here
                    unless (exists($seen_coords{$start}{$end})) {
                        
                        my $feature = Bio::SeqFeature::Generic->new(
                                                                    -seq_id => $self->{_seq}->display_id(), 
                                                                    -start  => $start,
                                                                    -end    => $end,
                                                                    -strand => $hsp->strand(),
                                                                    -source => 'blastx',
                                                                    -score  => $hsp->evalue(),
                                                                );
                        
                        $self->{_seq}->add_SeqFeature($feature);
                        
                        $seen_coords{$start}{$end} = 1;
                        
                    }
                    
                }   
                
            }
            
        }

    }
    
    
    ## Pitch any that do not have good translations
    $self->_check_translation();

    ## Get rid of any features contained by other features before we extend them
    $self->_remove_redundant();

    ## Extend the features, looking for start and stop codons
    $self->_extend_features();

    ## Get rid of any features contained by other features (again)
    $self->_remove_redundant();

    @{$self->{_genes}} = $self->{_seq}->all_SeqFeatures();
    $self->{_seq}->flush_SeqFeatures();
    
}

sub genes {

    my ($self) = @_;


    return $self->{_genes};

}

sub _mask_seq {

    my ($self) = @_;


    my @masks = (@{$self->{_masks}});
    $self->{_masks} = [ ];
    
    unless (@masks > 0) {
        $self->{_masked_seq} = $self->{_seq};
        return;
    }
    
    foreach my $mask (@masks) {
        
        my $start  = $mask->start();
        my $end    = $mask->end();
        
        if ($start > $end) { ($start, $end) = ($end, $start) }
        
        $mask->start($start);
        $mask->end($end);
        
    }

    @masks = sort {
        $a->start() <=> $b->start()
            ||
        $a->length() <=> $b->length()
    } @masks;

    my @composite_masks = ( );
    
    my $composite = $masks[0];     ## Prime the pump
    
  MASK: foreach my $mask (@masks) {
        
        ## If we hit a gap, stop extending, store the composite mask and
        ## move on to the next mask.
        if ($mask->start() > ($composite->end() + 1)) {
            
            push @composite_masks, $composite;
            $composite = $mask;
            next MASK;
        }
        
        ## No gap, extend the composite to include the current mask.
        if ($mask->end() > $composite->end()) {
            $composite->end($mask->end());
        }
        
    }
    
    ## After the loop, there will always be a leftover mask that had
    ## no masks after it to extend it.
    push @composite_masks, $composite;
    
    
    @composite_masks = sort {
        $a->start()  <=> $b->start()
            ||
        $a->length() <=> $b->length()
    } @composite_masks;

    my @intergenic = ( );
    
    ## If there is a gap between the start of the sequence and the first mask,
    ## create an intergenic region to cover it.
    my $first_mask = $composite_masks[0];
    
    if ($first_mask->start() > 1) {
        
        my $start  = 1;
        my $end    = $first_mask->start() - 1;

        my $feature = Bio::SeqFeature::Generic->new(
                                                    -start => $start,
                                                    -end   => $end,
                                                );

        push @intergenic, $feature;
        
    }

    ## If there is a gap between the last mask and the end of the sequence,
    ## create an intergenic region to cover it.
    my $last_mask = $composite_masks[$#composite_masks];

    if ($last_mask->end() < $self->{_seq}->length()) {
        
        my $start  = $last_mask->end() + 1;
        my $end    = $self->{_seq}->length();;
     
        my $feature = Bio::SeqFeature::Generic->new(
                                                    -start => $start,
                                                    -end   => $end,
                                                );

        push @intergenic, $feature;

    }
    
    ## All the composite masks are disjoint, by virtue of the way they
    ## were created, so create intergenic regions in the gaps between
    ## them.
    if (@composite_masks >= 2) {

        my $first = shift @composite_masks;

        while (my $second = shift @composite_masks) {

            my $start  = $first->end() + 1;
            my $end    = $second->start() - 1;

            my $feature = Bio::SeqFeature::Generic->new(
                                                        -start => $start,
                                                        -end   => $end,
                                                    );
            
            push @intergenic, $feature;
            
            $first = $second;
            
        }
        
    }

    ## Expand each intergenic region by $expansion base pairs on each
    ## side (for a total expansion of 2 X $expansion base pairs), but
    ## don't expand off the ends of the sequence.
    foreach my $intergenic (@intergenic) {
        
        if (($intergenic->start() - $self->{_expansion}) < 1) {
            $intergenic->start(1);
        }
        else {
            $intergenic->start($intergenic->start() - $self->{_expansion});
        }

        if ($intergenic->end() + $self->{_expansion} > $self->{_seq}->length()) {
            $intergenic->end($self->{_seq}->length());
        }
        else {
            $intergenic->end($intergenic->end() + $self->{_expansion});
        }
        
    }
    
    @intergenic = sort {
        $a->start()  <=> $b->start()
            ||
        $a->length() <=> $b->length()
    } @intergenic;

    my $seq_string        = $self->{_seq}->seq();
    my $masked_seq_string = $self->{_mask_char} x length($seq_string);
    
    ## For each intergenic region, copy the sequence from the sequence
    ## to the 'masked' sequence, thereby 'unmasking' it.
    foreach my $intergenic (@intergenic) {

        my $start  = $intergenic->start();
        my $length = $intergenic->length();

        ## Offset to substr are 0-based, our coordinates are 1-based
        substr($masked_seq_string, $start - 1, $length, substr($seq_string, $start - 1, $length));
        
    }
    
    $self->{_masked_seq} = Bio::Seq->new(
                                         -seq => $masked_seq_string,
                                         -id  => $self->{_seq}->display_id(),
                                     );
    
}

sub _remove_redundant {

    my ($self) = @_;
    

    my @features = $self->{_seq}->all_SeqFeatures();
    
    unless (@features > 1) { return; }

    @features = sort {
        $a->start() <=> $b->start()
            ||
        $a->length() <=> $b->length()
    } @features;

    my @filtered_features = ( );
    
  OUTER: foreach my $o (0..$#features) {
        
        my $f = $features[$o];
        
        if ($o < $#features) {
            
          INNER1: foreach my $i ($o+1..$#features) {
                
                my $g = $features[$i];
                
                if ($g->start() > $f->end()) {
                    
                    last INNER1;
                    
                }
                
                if (
                    ($f->start() >= $g->start())
                    &&
                    ($f->end() <= $g->end())
                    &&
                    ($f->length() <= $g->length())
                    &&
                    (!$g->has_tag('deleted'))
                ) {
                    
                    $f->add_tag_value('deleted', 1);
                    next OUTER;
                    
                }
                
            }
            
        }
        
        if ($o > 0) {
            
          INNER2: foreach my $i (reverse(0..$o-1)) {
                
                my $g = $features[$i];
                
                if ($g->end() < $f->start()) {
                    
                    last INNER2;
                    
                }
                
                if (
                    ($f->start() >= $g->start())
                    &&
                    ($f->end() <= $g->end())
                    &&
                    ($f->length() <= $g->length())
                    &&
                    (!$g->has_tag('deleted'))
                ) {
                    
                    $f->add_tag_value('deleted', 1);
                    next OUTER;
                    
                }
                
            }
            
        }
        
        push @filtered_features, $f;
        
    }

    $self->{_seq}->flush_SeqFeatures();
    
    if (@filtered_features) {
        $self->{_seq}->add_SeqFeature(@filtered_features);
    }
        
}

sub _extend_features {

    my ($self) = @_;


    my @features          = $self->{_seq}->all_SeqFeatures();
    my @filtered_features = ( );
    
  FEAT: foreach my $feature (@features) {

        ## There should be other checks to filter out candidate
        ## genes shorter than a certain threshold.  However, a
        ## last, stopgap sanity check is still in order.  Make
        ## sure we have at least one codon.  
        unless ($feature->length() >= 3) {
            carp 'skipping short (<3bp) gene candidate';
            next FEAT;
        }

        ## Should not get passed a candidate gene that isn't a
        ## multiple of 3, but check and complain anyway
        unless (($feature->length() % 3) == 0) {
            carp 'possibly extending gene candidate with length not a multiple of 3';
        }

        my $first_codon = substr($feature->seq->seq(), 0, 3);

        while (!$self->_is_start($first_codon)) {

            ## Make sure we don't hit a stop codon
            if ($self->_is_stop($first_codon)) { next FEAT; }

            ## No Ns or Xs allowed
            if ($first_codon =~ /N/io) { next FEAT; }
            if ($first_codon =~ /X/io) { next FEAT; }

            ## Make sure we don't run off the start of the sequence
            if ($feature->strand > 0) {
                if (($feature->start() - 3) < 1) { next FEAT; }
                $feature->start($feature->start() - 3);
            }
            else {
                if ($feature->end() + 3 > $self->{_seq}->length()) { next FEAT; }
                $feature->end($feature->end() + 3);
            }

            $first_codon = substr($feature->seq->seq(), 0, 3)
            
        }

        my $last_codon = substr($feature->seq->seq(), -3, 3); 

        while (!$self->_is_stop($last_codon)) {

            ## No Ns or Xs allowed
            if ($last_codon =~ /N/io) { next FEAT; }
            if ($last_codon =~ /X/io) { next FEAT; }

            ## Make sure we don't run off the end of the sequence
            if ($feature->strand > 0) {
                if (($feature->end() + 3) > $self->{_seq}->length()) { next FEAT; }
                $feature->end($feature->end() + 3);
            }
            else {
                if (($feature->start() - 3) < 1) { next FEAT; }
                $feature->start($feature->start() - 3);
            }

            $last_codon = substr($feature->seq->seq(), -3, 3);

        }

        ## However far we had to go, we found start and stop codons, so
        ## keep the feature
        push @filtered_features, $feature;

    }

    $self->{_seq}->flush_SeqFeatures();
    
    if (@filtered_features) {
        $self->{_seq}->add_SeqFeature(@filtered_features);
    }
    
}

sub _check_translation {
    
    my ($self) = @_;
    
    
    my @features          = $self->{_seq}->all_SeqFeatures(); 
    my @filtered_features = ( );
    
  FEAT: foreach my $feature (@features) {
        
        ## Should not be possible, but get rid of anything
        ## not an even multiple of 3
        unless (($feature->length() % 3) == 0) {
            carp 'discarding gene candidate with length not a multiple of 3';
            next FEAT;
        }

        my $feature_seq    = $feature->seq();
        my $feature_string = $feature_seq->seq();
        
        ## Get rid of anything with Ns or Xs
        if (
            ($feature_string =~ /N/io)
            ||
            ($feature_string =~ /X/io)
        ){
            next FEAT;
        }

        my $protein_seq    = $feature_seq->translate();
        my $protein_string = $protein_seq->seq();

        ## Discard internal stops
        if ($protein_string =~ /\*.+/o) {
            next FEAT;
        }

        ## Feature must be OK, no Ns, Xs or internal stops
        push @filtered_features, $feature;

    }

    $self->{_seq}->flush_SeqFeatures();
    
    if (@filtered_features) {
        $self->{_seq}->add_SeqFeature(@filtered_features);
    }
    
}

sub _is_start {

    my ($self, $codon) = @_;


    my $lc_codon = lc($codon);

    if ($lc_codon eq 'atg') { return 1; }
    if ($lc_codon eq 'ctg') { return 1; }
    if ($lc_codon eq 'gtg') { return 1; }
    if ($lc_codon eq 'ttg') { return 1; }

    return 0;

}

sub _is_stop {

    my ($self, $codon) = @_;
    
    
    my $lc_codon = lc($codon);
    
    if ($lc_codon eq 'taa') { return 1; }
    if ($lc_codon eq 'tag') { return 1; }
    if ($lc_codon eq 'tga') { return 1; }
    
    return 0;
    
}

1;
