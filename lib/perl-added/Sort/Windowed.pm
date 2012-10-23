package Sort::Windowed;
use strict;
use warnings;

sub new_iterator {
    my ($in, $bsize) = @_;
    my @buffer = map { [] } (1..$bsize);
    my ($p,$r);
    my ($minp, $last_p);
    my ($minbp, $lastbp);
    my $initialized;

    my $in_complete;
    my @out;    
    my $i = sub {
        unless ($initialized) {
            my ($p,$r) = $in->();
            return if not defined $r;
            die "undefined sort key returned from input iterator for item: $r.  expected (key,item) from the iterator" if not defined $p;
            $minp = $p;
            $minbp = $minp % $bsize;
            push @{$buffer[$minbp]}, $r;
            $last_p = 0;
            $initialized = 1;
        }
        until (@out or $in_complete) {
            ($p,$r) = $in->();
            if (defined $r) {
                #print "got $p: $r\n";
                if ($p < $last_p) {
                    if ($last_p - $p > $bsize) {
                        die "index $p found after index $last_p: difference is greater than window size $bsize!";
                    }
                }
                my $bp = ($p % $bsize);

                my $flushes=0;
                while ($p - $minp >= $bsize) {
                    if (@{$buffer[$minbp]}) {
                        ##print "$p requires $minp be flushed rom $minbp\n";
                        for my $e (sort @{$buffer[$minbp]}) {
                            #print "printing $e\n";
                            push @out, $e;
                        }
                        @{$buffer[$minbp]} = ();
                    }
                    $minp++;
                    $minbp++;
                    $minbp = 0 if $minbp == $bsize;    
                    $flushes++;
                    if ($flushes >= $bsize) {
                        # no point in going around the buffer more than once: it's empty
                        $minp = $p;
                        $minbp = $minp % $bsize;
                        #print "set minp $minp minbp $minbp\n";
                        last;
                    }
                }
                if ($p < $minp) {
                    $minp = $p;
                    $minbp = $minp % $bsize;
                }
                #print "storing $p at $bp: min is $minbp\n";
                push @{$buffer[$bp]}, $r;
                $last_p = $p;
            }
            else { 
                #print "at final\n";
                #print "final\n";
                for my $bp (($minbp..$bsize-1), (0..$minbp-1)) {
                    for my $e (sort @{$buffer[$bp]}) {
                        #print "printing $e\n";
                        push @out, $e;
                    }
                }
                $in_complete = 1;
            } 
        }
        #print (ref($out[0]) ? "s: @{$out[0]}\n" : "s: $out[0]\n");
        return shift @out;
    };
    return $i;
}

sub sort {
    my ($in,$out,$bsize) = @_;

    my @buffer = map { [] } (1..$bsize);

    my ($p,$r) = $in->();
    return if not defined $r;
    die "undefined sort key returned from input iterator for item: $r.  expected (key,item) from the iterator" if not defined $p;
    my $minp = $p;
    my $minbp = $minp % $bsize;
    push @{$buffer[$minbp]}, $r;
    
    my $last_p = 0;
    for (($p,$r)=$in->(); defined($r); ($p,$r) = $in->()) {
        #print "got $p: $r\n";
        if ($p < $last_p) {
            if ($last_p - $p > $bsize) {
                die "index $p found after index $last_p: difference is greater than window size $bsize!";
            }
        }
        my $bp = ($p % $bsize);

        my $flushes=0;
        while ($p - $minp >= $bsize) {
            if (@{$buffer[$minbp]}) {
                #print "$p requires $minp be flushed rom $minbp\n";
                for my $e (sort @{$buffer[$minbp]}) {
                    #print "printing $e\n";
                    $out->($e)
                }
                @{$buffer[$minbp]} = ();
            }
            $minp++;
            $minbp++;
            $minbp = 0 if $minbp == $bsize;    
            $flushes++;
            if ($flushes >= $bsize) {
                # no point in going around the buffer more than once: it's empty
                $minp = $p;
                $minbp = $minp % $bsize;
                #print "set minp $minp minbp $minbp\n";
                last;
            }
        }
        if ($p < $minp) {
            $minp = $p;
            $minbp = $minp % $bsize;
        }
        #print "storing $p at $bp: min is $minbp\n";
        push @{$buffer[$bp]}, $r;
        $last_p = $p;
    }

    #print "final\n";
    for my $bp (($minbp..$bsize-1), (0..$minbp-1)) {
        for my $e (sort @{$buffer[$bp]}) {
            #print "printing $e\n";
            $out->($e);
        }
    }
}
1;

=pod

=head1 NAME

Sort::Windowed - sort items with limited disorder 

=head1 SYNOPSIS

    use Sort::Windowed;
    use IO::File;

    my $fin = IO::File->new("in.tsv");
    my $fout = IO::File->new(">out.tsv");

    my $in = sub { my $r = $fin->getline(); my @f = split(/\t/,$r); return ($f[1],$r); };
    my $out = sub { $fout->print($_[0]) };

    Sort::Windowed::sort($in,$out,50);

    # or

    my $i = Sort::Windowed::new_iterator($in,50);
    while (my $next = $i->()) {
        $fout->print($next);
    }

=head1 DESCRIPTION

This sorter takes an input iterator and an output iterator, and returns sorted data as long as
there is a limited amount of disorder in the input.  The 3rd parameter is a window size for a 
circular buffer used in sorting.  If any item has a position _before_ a prior input it the position
difference must be less than the window size or the sorter will throw an exception.

=head1 COPYRIGHT AND DISCLAIMERS

Copyright (c) 2008 Scott Smith and Washington University Genome Center

This library is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.

=head1 AUTHOR

Scott Smith ssmith@genome.wustl.edu

=cut
