package Genome::Model::Tools::Sword::SlidingWindowPileup;

use strict;
use warnings;

use Bio::DB::Sam;
#use POSIX;

class Genome::Model::Tools::Sword::SlidingWindowPileup {
    is => 'Genome::Model::Tools::Sword::Base',
    doc => 'Add documentation here',
    has => [
        bam_file_fof => {
            is => 'Text',
            doc => 'I need a BAM file',
        },
        scaling => {
            is => 'Text',
            doc => 'A file for scaling...',
        },
        window_size => {
            is => 'Integer',
            is_optional => 1,
            doc => 'The size of the window',
            default_value => 1000,
            valid_values => [100,1000,10000],
        },
        threshold => {

        },
        region => {

        },
        output_file => {
        },
    ],
};

sub execute {
    my $self = shift;
    
    if ($self->threshold < 0) {
        $self->error_message('Threshold must be positive, but you input:'. $self->threshold);
        die($self->error_message);
    }
    ###### Read in sample weights ########
    my $scaling_fh = Genome::Sys->open_file_for_reading($self->scaling);
    my @scaling;
    while ($scaling_fh->getline) {
        chomp;
        push @scaling, $_;
    }
    ###### Calculate number of bins #####
    open OUTFILE, ">" . $self->output_file or die $!;	# open output file

    my $chr = (split /:/, $self->region)[0] ;
    my $temp =  (split /:/, $self->region)[1] ;
    my $startpos = (split /-/, $temp)[0] ;
    my $endpos = (split /-/, $temp)[1] ;
    my $c = $endpos - $startpos ;
    my $bin = ceil( $c / $self->window_size ) ;
    
    print "# Bins = $bin\n";
    for my $b (1..$bin) {
        my $sp = $startpos + ($b-1)*$self->window_size ;
        print OUTFILE "$sp\t";
    }
    print OUTFILE "\n";

#### Extract and print average coverage ## 
open SAMP, "<" . $self->bam_file_fof or die$!;
while (<SAMP>) {
    chomp ;
    my $bfilet = $_ ;
    my $bamt = Bio::DB::Sam->new(-bam => $bfilet);

    $_ = <SAMP>;
    chomp ;
    my $bfilen = $_ ;		
    my $bamn = Bio::DB::Sam->new(-bam => $bfilen);

    my $scale = shift @scaling; 

    my $countt;
    my $countn;
    for my $b (1..$bin) {
        my $sp = $startpos + ($b-1)*$self->window_size ;
        my $ep = $startpos + $b*$self->window_size - 1 ;
            ($ep = $endpos) if ($ep > $endpos) ;
	my $seqid = $chr . ":" . $sp . "-" . $ep ;

        if ( ($chr =~ "11") & ($sp >= 65265000) & ($ep <= 65274000) ) {	
            # skip MALAT1 region
            print OUTFILE "0\t";
            next ;
        } elsif ( ($chr =~ "chr11") & ($sp >= 65265000) & ($ep <= 65274000) ) {	
            # skip MALAT1 region
            print OUTFILE "0\t";
            next ;
        }


        my @alignt = $bamt->get_features_by_location($seqid);
        my @alignn = $bamn->get_features_by_location($seqid);

        $countt = 0;
        for my $a (@alignt) {
            my $NH = $a->aux_get("NH");
            my $pos = $a->pos;
            my $mate_is_unmapped = $a->munmapped;

            next if $mate_is_unmapped;
            next if ($NH > 1) ;
            next unless ( ($pos >= $sp) & ($pos <= $ep) ) ;

            $countt++;
        } 

        $countn = 0;
        for my $a (@alignn) {
            my $NH = $a->aux_get("NH");
            my $pos = $a->pos;
            my $mate_is_unmapped = $a->munmapped;

            next if $mate_is_unmapped;
            next if ($NH > 1) ;
            next unless ( ($pos >= $sp) & ($pos <= $ep) ) ;

            $countn++;
        } 

        my $logr = log2( (($countt + 1) / ($countn + 1)) / $scale ) ;
        my $cov = ($countn + $countt);
	if ($cov < ($self->threshold*$self->window_size)) { 
		$logr = $logr * ($cov/($self->threshold*$self->window_size)) ;
	} 
        print OUTFILE sprintf('%0.3f',$logr) . "\t";

    }

    print OUTFILE "\n";			# start new line for next sample

}

    close OUTFILE;
    close SAMP;
    return 1;
}


######## Subroutines #################

sub log2 {
    my $n = shift;
    return log($n)/log(2);
}

1;
