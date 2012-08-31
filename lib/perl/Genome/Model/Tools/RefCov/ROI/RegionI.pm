package Genome::Model::Tools::RefCov::ROI::RegionI;

use strict;
use warnings;

use Carp;
use integer;
use vars qw(%STRAND_OPTIONS);

BEGIN {
# STRAND_OPTIONS contains the legal values for the strand-testing options
    %STRAND_OPTIONS = map { $_, '_' . $_ }
    (
     'strong', # regions must have the same strand
     'weak',   # regions must have the same strand or no strand
     'ignore', # ignore strand information
     );
}

class Genome::Model::Tools::RefCov::ROI::RegionI {
    has => [
        start => {},
        end => {},
        strand => {
            is_optional => 1,
            default_value => 0,
        },
    ]
};

sub length {
    my $self = shift;
    return $self->end - $self->start + 1;
}


# utility methods
#

# returns true if strands are equal and non-zero
sub _strong {
    my ($r1, $r2) = @_;
    my ($s1, $s2) = ($r1->strand(), $r2->strand());

    return 1 if $s1 != 0 && $s1 == $s2;
}

# returns true if strands are equal or either is zero
sub _weak {
    my ($r1, $r2) = @_;
    my ($s1, $s2) = ($r1->strand(), $r2->strand());
    return 1 if $s1 == 0 || $s2 == 0 || $s1 == $s2;
}

# returns true for any strandedness
sub _ignore {
    return 1;
}

# works out what test to use for the strictness and returns true/false
# e.g. $r1->_testStrand($r2, 'strong')
sub _testStrand() {
    my ($r1, $r2, $comp) = @_;
    return 1 unless $comp;
    my $func = $STRAND_OPTIONS{$comp};
    return $r1->$func($r2);
}


=head1 Boolean Methods

These methods return true or false. They throw an error if start and
end are not defined.

  $region->overlaps($otherRegion) && print "Regions overlap\n";

=head2 overlaps

  Title   : overlaps
  Usage   : if($r1->overlaps($r2)) { do stuff }
  Function: tests if $r2 overlaps $r1
  Args    : arg #1 = a region to compare this one to (mandatory)
            arg #2 = optional strand-testing arg ('strong', 'weak', 'ignore')
  Returns : true if the regions overlap, false otherwise

=cut

sub overlaps {
    my ($self, $other, $so) = @_;

    $self->throw("start is undefined") unless defined $self->start;
    $self->throw("end is undefined") unless defined $self->end;
    $self->throw("not a Genome::Model::Tools::RefCov::ROI::RegionI object") unless defined $other &&
        $other->isa('Genome::Model::Tools::RefCov::ROI::RegionI');
    $other->throw("start is undefined") unless defined $other->start;
    $other->throw("end is undefined") unless defined $other->end;
    return
        ($self->_testStrand($other, $so)
             and not (
                 ($self->start() > $other->end() or
                      $self->end() < $other->start()   )
             ));
}

=head2 contains

  Title   : contains
  Usage   : if($r1->contains($r2) { do stuff }
  Function: tests whether $r1 totally contains $r2
  Args    : arg #1 = a region to compare this one to (mandatory)
	             alternatively, integer scalar to test
            arg #2 = optional strand-testing arg ('strong', 'weak', 'ignore')
  Returns : true if the argument is totally contained within this region

=cut

sub contains {
    my ($self, $other, $so) = @_;
    $self->throw("start is undefined") unless defined $self->start;
    $self->throw("end is undefined") unless defined $self->end;
    if(defined $other && ref $other) { # a region object?
        $other->throw("Not a Genome::Model::Tools::RefCov::ROI::RegionI object: $other") unless  $other->isa('Genome::Model::Tools::RefCov::ROI::RegionI');
        $other->throw("start is undefined") unless defined $other->start;
        $other->throw("end is undefined") unless defined $other->end;
        return ($self->_testStrand($other, $so)      and
                    $other->start() >= $self->start() and
                        $other->end() <= $self->end());
    } else { # a scalar?
        $self->throw("'$other' is not an integer.\n") unless $other =~ /^[-+]?\d+$/;
        return ($other >= $self->start() and $other <= $self->end());
    }
}

=head2 equals

  Title   : equals
  Usage   : if($r1->equals($r2))
  Function: test whether $r1 has the same start, end, length as $r2
  Args    : arg #1 = a region to compare this one to (mandatory)
            arg #2 = optional strand-testing arg ('strong', 'weak', 'ignore')
  Returns : true if they are describing the same region

=cut

sub equals {
    my ($self, $other, $so) = @_;

    $self->throw("start is undefined") unless defined $self->start;
    $self->throw("end is undefined") unless defined $self->end;
    $other->throw("Not a Genome::Model::Tools::RefCov::ROI::RegionI object") unless  $other->isa('Genome::Model::Tools::RefCov::ROI::RegionI');
    $other->throw("start is undefined") unless defined $other->start;
    $other->throw("end is undefined") unless defined $other->end;

    return ($self->_testStrand($other, $so)   and
                $self->start() == $other->start() and
                    $self->end()   == $other->end()       );
}

=head1 Geometrical methods

These methods do things to the geometry of regions, and return
Genome::Model::Tools::RefCov::ROI::RegionI compliant objects or triplets (start, stop, strand) from
which new regions could be built.

=head2 intersection

 Title   : intersection
 Usage   : ($start, $stop, $strand) = $r1->intersection($r2); OR
           ($start, $stop, $strand) = Genome::Model::Tools::RefCov::Region->intersection(\@regions); OR
           my $containing_region = $r1->intersection($r2); OR
           my $containing_region = Genome::Model::Tools::RefCov::Region->intersection(\@regions);
 Function: gives the region that is contained by all regions
 Returns : undef if they do not overlap, or
           the region that they do overlap (in the form of an object
            like the calling one, OR a three element array)
 Args    : arg #1 = [REQUIRED] a region to compare this one to,
                    or an array ref of regions
           arg #2 = optional strand-testing arg ('strong', 'weak', 'ignore')

=cut

sub intersection {
    my ($self, $given, $so) = @_;
    $self->throw("missing arg: you need to pass in another feature") unless $given;

    my @regions;
    if ($self eq "Genome::Model::Tools::RefCov::ROI::RegionI") {
        $self = "Genome::Model::Tools::RefCov::Region";
        $self->warn("calling static methods of an interface is deprecated; use $self instead");
    }
    if (ref $self) {
        push(@regions, $self);
    }
    ref($given) eq 'ARRAY' ? push(@regions, @{$given}) : push(@regions, $given);
    $self->throw("Need at least 2 regions") unless @regions >= 2;

    my $intersect;
    while (@regions > 0) {
        unless ($intersect) {
            $intersect = shift(@regions);
            $self->throw("Not an object: $intersect") unless ref($intersect);
            $self->throw("Not a Genome::Model::Tools::RefCov::ROI::RegionI object: $intersect") unless $intersect->isa('Genome::Model::Tools::RefCov::ROI::RegionI');
            $self->throw("start is undefined") unless defined $intersect->start;
            $self->throw("end is undefined") unless defined $intersect->end;
        }

        my $compare = shift(@regions);
        $self->throw("Not an object: $compare") unless ref($compare);
        $self->throw("Not a Genome::Model::Tools::RefCov::ROI::RegionI object: $compare") unless $compare->isa('Genome::Model::Tools::RefCov::ROI::RegionI');
        $self->throw("start is undefined") unless defined $compare->start;
        $self->throw("end is undefined") unless defined $compare->end;
        return unless $compare->_testStrand($intersect, $so);

        my @starts = sort {$a <=> $b} ($intersect->start(), $compare->start());
        my @ends   = sort {$a <=> $b} ($intersect->end(), $compare->end());

        my $start = pop @starts; # larger of the 2 starts
        my $end = shift @ends;   # smaller of the 2 ends

        my $intersect_strand;    # strand for the intersection
        if (defined($intersect->strand) && defined($compare->strand) && $intersect->strand == $compare->strand) {
            $intersect_strand = $compare->strand;
        }
        else {
            $intersect_strand = 0;
        }

        if ($start > $end) {
            return;
        }
        else {
            $intersect = $self->class->create(
                start  => $start,
                end    => $end,
                strand => $intersect_strand
            );
        }
    }

    if (wantarray()) {
        return ($intersect->start, $intersect->end, $intersect->strand);
    }
    else {
        return $intersect;
    }
}

=head2 union

   Title   : union
    Usage   : ($start, $stop, $strand) = $r1->union($r2);
            : ($start, $stop, $strand) = Genome::Model::Tools::RefCov::Region->union(@regions);
              my $newregion = Genome::Model::Tools::RefCov::Region->union(@regions);
    Function: finds the minimal Region that contains all of the Regions
    Args    : a Region or list of Region objects
    Returns : the region containing all of the region
              (in the form of an object like the calling one, OR
              a three element array)

=cut

sub union {
	my $self = shift;
	my @regions = @_;
	if ($self eq "Genome::Model::Tools::RefCov::ROI::RegionI") {
		$self = "Genome::Model::Tools::RefCov::Region";
		$self->warn("calling static methods of an interface is deprecated; use $self instead");
	}
	if(ref $self) {
		unshift @regions, $self;
	}

	my @start = sort {$a<=>$b}
	  map( { $_->start() } @regions);
	my @end   = sort {$a<=>$b}
	  map( { $_->end()   } @regions);

	my $start = shift @start;
	while( !defined $start ) {
		$start = shift @start;
	}

	my $end = pop @end;

	my $union_strand;  # Strand for the union region object.

	foreach(@regions) {
		if(! defined $union_strand) {
			$union_strand = $_->strand;
			next;
		} else {
			if(not defined $_->strand or $union_strand ne $_->strand) {
				$union_strand = 0;
				last;
			}
		}
	}
	return unless $start or $end;
	if( wantarray() ) {
		return ( $start,$end,$union_strand);
	} else {
		return $self->class->create(
                    start => $start,
                    end => $end,
                    strand => $union_strand,
                );
	}
}

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : ($a_unique,$common,$b_unique) = $a->overlap_extent($b)
 Function: Provides actual amount of overlap between two different
           regions
 Example :
 Returns : array of values containing the length unique to the calling
           region, the length common to both, and the length unique to
           the argument region
 Args    : a region

=cut

sub overlap_extent{
	my ($a,$b) = @_;

	$a->throw("start is undefined") unless defined $a->start;
	$a->throw("end is undefined") unless defined $a->end;
	$b->throw("Not a Genome::Model::Tools::RefCov::ROI::RegionI object") unless  $b->isa('Genome::Model::Tools::RefCov::ROI::RegionI');
	$b->throw("start is undefined") unless defined $b->start;
	$b->throw("end is undefined") unless defined $b->end;

	if( ! $a->overlaps($b) ) {
		return ($a->length,0,$b->length);
	}

	my ($au,$bu) = (0, 0);
	if( $a->start < $b->start ) {
		$au = $b->start - $a->start;
	} else {
		$bu = $a->start - $b->start;
	}

	if( $a->end > $b->end ) {
		$au += $a->end - $b->end;
	} else {
		$bu += $b->end - $a->end;
	}

	my $intersect = $a->intersection($b);
	my $ie = $intersect->end;
	my $is = $intersect->start;

	return ($au,$ie-$is+1,$bu);
}

=head2 disconnected_regions

    Title   : disconnected_regions
    Usage   : my @disc_regions = Genome::Model::Tools::RefCov::Region->disconnected_regions(@regions);
    Function: finds the minimal set of regions such that each input region
              is fully contained by at least one output region, and none of
              the output regions overlap
    Args    : a list of regions
    Returns : a list of objects of the same type as the input
              (conforms to RegionI)

=cut

sub disconnected_regions {
    my $self = shift;
    if ($self eq "Genome::Model::Tools::RefCov::ROI::RegionI") {
	$self = "Genome::Model::Tools::RefCov::Region";
	$self->warn("calling static methods of an interface is deprecated; use $self instead");
    }
    my @inregions = @_;
    if(ref $self) {
	unshift @inregions, $self;
    }

    my @outregions = (); # disconnected regions

    # iterate through all input regions $inregion,
    # adding each input region to the set of output regions @outregions,
    # provided $inregion does not overlap ANY region in @outregions
    # - if it does overlap an outregion, then merge it
    foreach my $inregion (@inregions) {
	my $intersects = 0;
	my @outregions_new = ();
	my @intersecting_regions = ();

        # iterate through all @outregions, testing if it intersects
        # current $inregion; if it does, merge and add to list
        # of @intersecting_regions, otherwise add $outregion to
        # the new list of outregions that do NOT intersect
	for (my $i=0; $i<@outregions; $i++) {
	    my $outregion = $outregions[$i];
	    my $intersection = $inregion->intersection($outregion);
	    if ($intersection) {
		$intersects = 1;
		my $union = $inregion->union($outregion);
		push(@intersecting_regions, $union);
	    }
	    else {
		push(@outregions_new, $outregion);
	    }
	}
	@outregions = @outregions_new;
        # @outregions now contains a list of non-overlapping regions
        # that do not intersect the current $inregion

	if (@intersecting_regions) {
	    if (@intersecting_regions > 1) {
		# this sf intersected > 1 region, which means that
		# all the regions it intersects should be joined
		# together in a new region
                my $merged_region =
                  $self->union(@intersecting_regions);
		push(@outregions, $merged_region);

	    }
	    else {
		# exactly 1 intersecting region
		push(@outregions, @intersecting_regions);
	    }
	}
	else {
	    # no intersections found - new region
	    push(@outregions,
		 $self->class->create(
                     start => $inregion->start,
                     end => $inregion->end,
                     strand => $inregion->strand,
                 ));
	}
    }
    return @outregions;
}

=head2 offsetStranded

    Title    : offsetStranded
    Usage    : $rnge->ofsetStranded($fiveprime_offset, $threeprime_offset)
    Function : destructively modifies RegionI implementing object to
               offset its start and stop coordinates by values $fiveprime_offset and
               $threeprime_offset (positive values being in the strand direction).
    Args     : two integer offsets: $fiveprime_offset and $threeprime_offset
    Returns  : $self, offset accordingly.

=cut

sub offsetStranded {
  my ($self, $offset_fiveprime, $offset_threeprime) = @_;
  my ($offset_start, $offset_end) = $self->strand() eq -1 ? (- $offset_threeprime, - $offset_fiveprime) : ($offset_fiveprime, $offset_threeprime);
  $self->start($self->start + $offset_start);
  $self->end($self->end + $offset_end);
  return $self;
};

=head2 subtract

  Title   : subtract
  Usage   : my @subtracted = $r1->subtract($r2)
  Function: Subtract region r2 from region r1
  Args    : arg #1 = a region to subtract from this one (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : undef if they do not overlap or r2 contains this RegionI,
            or an arrayref of Region objects (this is an array since some
            instances where the subtract region is enclosed within this region
            will result in the creation of two new disjoint regions)

=cut

sub subtract() {
   my ($self, $region, $so) = @_;
    $self->throw("missing arg: you need to pass in another feature")
      unless $region;
    return unless $self->_testStrand($region, $so);

    if ($self eq "Genome::Model::Tools::RefCov::ROI::RegionI") {
	$self = "Genome::Model::Tools::RefCov::Region";
	$self->warn("calling static methods of an interface is
deprecated; use $self instead");
    }
    $region->throw("Input a Genome::Model::Tools::RefCov::ROI::RegionI object") unless
$region->isa('Genome::Model::Tools::RefCov::ROI::RegionI');

    if (!$self->overlaps($region)) {
        return undef;
    }

    ##Subtracts everything
    if ($region->contains($self)) {
        return undef;   
    }

    my ($start, $end, $strand) = $self->intersection($region, $so);
    ##Subtract intersection from $self region

    my @outregions = ();
    if ($self->start < $start) {
        push(@outregions,
             $self->class->create(
                 start => $self->start,
                 end => $start - 1,
                 strand => $self->strand,
             ))
    }
    if ($self->end > $end) {
        push(@outregions, 
             $self->class->create(
                 start => $end + 1,
                 end => $self->end,
                 strand => $self->strand,
             ));
    }
    return \@outregions;
}

# THIS MODULE IS A TRIMMED DOWN VERSION OF THE BIOPERL BIO::RANGE CLASS


=head2 unions

 Title   : unions
 Usage   : @unions = Genome::Region->unions(@regions);
 Function: generate a list of non-intersecting Genome::Region objects
           from a list of Genome::Region objects which may intersect
 Returns : a list of Genome::Region objects
 Args    : a list of Genome::Region objects


=cut

sub unions {
  my ($class,@i) = @_;
  die('You probably wanted the method disconnected_regions in Genome::Model::Tools::RefCov::ROI::RegionI');
  if (ref($class)) {
      $class = ref($class);
  }
  @i = sort{ $a->start <=> $b->start } @i;

  for (my $j = 0; $j < scalar(@i); $j++) {
      for (my $k = 0; $k < scalar(@i); $k++) {
          #it may have been replaced by a union under the key of
          #the overlapping range, we are altering the hash in-place
          next unless $i[$j];
          next if $i[$k]->end   < $i[$j]->start;
          last if $i[$k]->start > $i[$j]->end;
          if($i[$j]->overlaps($i[$k])){
              if ($i[$j]->equals($i[$k])) {next;}
              my ($start,$end,$strand) = $i[$j]->union($i[$k]);
              splice(@i,$k,1);
              $i[$j] = $class->create( start => $start , end => $end , strand => $strand );
          }
      }
  }

  return @i;
}


=head2 toString

  Title   : toString
  Function: stringifies this region
  Example : print $range->toString(), "\n";
  Returns : a string representation of this range

=cut

sub toString {
  my $self = shift;
  return  "(${\$self->start}, ${\$self->end}) strand=${\$self->strand}";
}



sub throw {
    my $self = shift;
    my $error_message = shift;
    confess($error_message);
}

1;
