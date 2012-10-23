package Genome::Model::GenePrediction::Command::Pap::Blast::BlastTopHitLogic;

# bdericks: This module was originally written by Todd Wylie and named 
# BlastTopHitLogic. As of August 27, 2010, the original module can be 
# found at /gscmnt/233/analysis/sequence_analysis/lib. I wasn't able to 
# find it in version control anywhere. I've made some changes to turn 
# this into a command object and put in the PAP namespace

use PAP;
use IO::File;

# FIXME Really need to remove this dependency, but porting this script
# into the PAP namespace may not be straightforward...
use lib "/gsc/scripts/gsc/compbio/lib";
use BPdeluxe;

class Genome::Model::GenePrediction::Command::Pap::Blast::BlastTopHitLogic {
    is => 'Command::V1',
    has => [
        report => {
            is => 'Path',
            doc => 'Path to report file',
        },
        return_type => {
            is => 'Text',
            valid_values => ['top_hit_report', 'full_report'],
            doc => 'Type of report to produce',
        },
    ],
    has_optional => [
        generated_report => {
            is => 'SCALAR',
            doc => 'Hash containing report information',
        },
    ],
};

sub execute {
    my $self = shift;
    my $report = $self->report;
    my $return_type = $self->return_type;

    my $report_in = IO::File->new($report);

    # FIXME BPdeluxe is a deployed module. This dependency needs to be removed!
    my $BP = new BPdeluxe::Multi($report_in);
    my %queries;
    while(my $multi = $BP->nextReport) {
        my $query = $multi->query;
        while(my $sbjct = $multi->nextSbjct) {	
            my $group_hashref = $sbjct->group_list("$query","return_hashref");
            push(@{$queries{$query}}, $group_hashref);	
        }
    }
    close($report_in);

    # ----------------------------------------------------------------------
    # F I L T E R   L O G I C 
    # ----------------------------------------------------------------------
    # The following logic filters the match list. Top hits are determined
    # on a one-to-one basis (1 group hit per query). The rules for hits
    # are:
    # 1) Lowest p_value determines top match.
    # 2) If other HSPs tie the lowest p_value, highest bit score prevails.
    # 3) If p_value, bit scores are identical, then choosing any 1 of the hits
    #    is equally representative.
    my (%full_report, %top_hit_report);
    my $query_count = 0;
    foreach my $query (sort keys %queries) {
        my %grouped;
        foreach my $group (@{$queries{$query}}) {
            $grouped{ $group->{subject_gi} } = $group;
        }
        my $i = 0;
        my $top_p_value;
        my %final_vals;
        foreach my $val (sort {$grouped{$a}->{p_value} <=> $grouped{$b}->{p_value}} keys %grouped) {
            $i++;
            if ($return_type eq "full_report") {
                $query_count++;
                $full_report{$query_count} = [ $query, $val, $grouped{$val} ];
            }
            $top_p_value = $grouped{$val}->{p_value} if ($i == 1);
            $final_vals{$val} = $grouped{$val} if ($grouped{$val}->{p_value} eq $top_p_value);	
        }
        my $ii = 0;
        foreach my $top_p (sort {$final_vals{$b}->{bit_score} <=> $final_vals{$a}->{bit_score}} keys %final_vals) {
            $ii++;
            if (($ii == 1) && ($return_type eq "top_hit_report")) {    
                $top_hit_report{$query} = [ $query, $top_p, $final_vals{$top_p} ]
            }
        }
    }

    # ----------------------------------------------------------------------
    # R E T U R N   H A S H R E F   T O   D A T A 
    # ----------------------------------------------------------------------
    if ($return_type eq "full_report") {
        $self->generated_report(\%full_report);
    }
    elsif ($return_type eq "top_hit_report") {
        $self->generated_report(\%top_hit_report);
    }
    return 1; 
}

