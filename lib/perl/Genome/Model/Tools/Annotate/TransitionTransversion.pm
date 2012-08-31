package Genome::Model::Tools::Annotate::TransitionTransversion;

use strict;
use warnings;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Annotate::TransitionTransversion {
    is => 'Command',                       
    has => [                                # specify the command's single-value properties (parameters) <--- 
        referance_allele      => { is => 'String',    doc => "give the referance allele" },
        variant_allele        => { is => 'String',    doc => "give the variant allele" },
        trans_type        => { is => 'String',    doc => "this is the result", is_optional => 1 },

    ], 
};

        
sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "give the referance and variant alleles in the same orientation"                 
}

sub help_synopsis {                         # replace the text below with real examples <---
    return <<EOS
gmt annotation transition-transversion --referance-allele=A --variant-allele=C
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
    give the referance and variant alleles in the same orientation and get returned the Transition/Transversion annotation
    for example gmt annotation transition-transversion --referance-allele=A --variant-allele=C
    should result in T:A>G:C(transversion)
EOS
}

sub execute {                               # replace with real execution logic.

    my $self = shift;

    my $ref_allele = $self->referance_allele;
    my $variant_allele = $self->variant_allele;

    my $alleles = "$ref_allele$variant_allele";

    my $trans;

    if ($alleles eq "AC") {
	$trans = "T:A" . ">" . "G:C(transversion)";
	#tg_ac
    } elsif ($alleles eq "AT") {
	$trans = "T:A" . ">" . "A:T(transversion)";
	#ta_at
    } elsif ($alleles eq "TA") {
	$trans = "T:A" . ">" . "A:T(transversion)";
	#ta_at
    } elsif ($alleles eq "TG") {
	$trans = "T:A" . ">" . "G:C(transversion)";
	#tg_ac
    } elsif ($alleles eq "CA") {
	$trans = "C:G" . ">" . "A:T(transversion)";
	#ca_gt
    } elsif ($alleles eq "CG") {
	$trans = "C:G" . ">" . "G:C(transversion)";
	#cg_gc
    } elsif ($alleles eq "GC") {
	$trans = "C:G" . ">" . "G:C(transversion)";
	#cg_gc
    } elsif ($alleles eq "GT") {
	$trans = "C:G" . ">" . "A:T(transversion)";
	#ca_gt
    } elsif ($alleles eq "AG") {
	$trans = "T:A" . ">" . "C:G(transition)";
	#tc_ag
    } elsif ($alleles eq "TC") {
	$trans = "T:A" . ">" . "C:G(transition)";
	#tc_ag
    } elsif ($alleles eq "CT") {
	$trans = "C:G" . ">" . "T:A(transition)";
	#ct_ga
    } elsif ($alleles eq "GA") {
	$trans = "C:G" . ">" . "T:A(transition)";
	#ct_ga
    } else {
	$trans = "ambiguos";
    }
    
    print qq($trans\n);
    $self->trans_type($trans);
    return 1;

}

1;
