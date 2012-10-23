package SnpDom;

use strict;
use warnings;
# comments
# FILE: find-muts-pfdom.pl
# VERSION: $Revision: 38920 $ #

use Carp;
use IO::File;
use Text::CSV_XS;
use File::Slurp;
use List::MoreUtils qw/ uniq /;

# moved the use statements for MPSampleData out of here.
# apparently anytime you try to use the help from genome-model
# it's going to connect to the database, due to this fun module.

use SnpDom::Mutation;

sub new {
    my $class = shift;
    my $arghash = shift;
    my %args = @_;
    my $self = { '_mcounts' => undef,
                 '_mutations' => undef,
                 '_ts_flag' => undef, };
    bless($self,$class);
    if(exists($arghash->{'-inc-ts'}) ) {
        $self->{'_ts_flag'} = 1;
    }
    return $self;
}


# subs start here 
sub read_muts {
    my ($self,$f,$cols) = @_;
    my $fh = new IO::File;
    my $c = new Text::CSV_XS({sep_char => "\t"});

    my ($g,$aa_change) = split(/,/x,$cols); # should be like 2,3
    my $retval = $fh->open($f); # goofy.  doesn't work in the if() construct...
    if(!defined($retval) ) {
        carp "can't open $f : $!";
        $fh->close;
        return 0;
    }
    while(<$fh>) {
        chomp;
        $c->parse($_);
        my @fields = $c->fields();
        #next if $fields[$aa_change] eq "NULL";
        # should do:
        # $self->add_mutation($gene,$ts,$aa_change); instead
        if(exists($self->{_mutations}->{$fields[$g]})) {
            my $sdm = $self->{_mutations}->{$fields[$g]};
            $sdm->add_mutation($fields[$aa_change]);
            $self->{_mutations}->{$fields[$g]} = $sdm;
        }
        else {
            my $sdm = new SnpDom::Mutation('gt' => $fields[$g]);
            $sdm->add_mutation($fields[$aa_change]);
            $self->{_mutations}->{$fields[$g]} = $sdm;

        }
    }
    $fh->close;
    return 1;
}

sub get_mut_obj {
    my ($self,$name) = @_;
    if(exists($self->{_mutations}->{$name})) {
        return $self->{_mutations}->{$name};
    }
    return;
}

sub get_all_mutations {
    my ($self) = @_;
    my %hash;
    foreach my $key (sort keys %{$self->{_mutations}}) {
        my $sdm = $self->{_mutations}->{$key};
        $hash{$key} = $sdm->get_all_mutations();
    }
    return \%hash;
}

# this one is fugly.
sub mutation_in_dom {
    # These will connect to the database on use right now...
    eval qq| 
        use MPSampleData::DBI;
        use MPSampleData::Gene;
        use MPSampleData::Transcript;
        use MPSampleData::Protein;
        use MPSampleData::IproGeneTranscriptXref;
        use MPSampleData::IproResults;
    |;
    die $@ if $@;

    my ($self,$dlen,$filter) = @_;
    #MPSampleData::DBI::myinit("dbi:Oracle:dwrac","mguser_prd");
    foreach my $gt (sort keys %{$self->{_mutations}}) {
        my %cache;
        my ($transcript,$gene) = split(/,/x,$gt);
        my ($gobj) = MPSampleData::Gene->search({hugo_gene_name => $gene});
        my ($tobj) = MPSampleData::Transcript->search({transcript_name => $transcript});
        unless(defined($gobj) && defined($tobj)) {
            next;
        }
        my (@xref) = MPSampleData::IproGeneTranscriptXref->search({gene_id => $gobj->gene_id,
                                                          transcript_id => $tobj->transcript_id});
        my %iprosdone;
IPRO:   foreach my $ref (@xref) {
            my $ipro = MPSampleData::IproResults->search({ipro_id => $ref->ipro_id});
            next unless defined($ipro);
RECORD:     while(my $i = $ipro->next) {
            if(!exists($iprosdone{$ref->ipro_id}) ) {
                #$iprosdone{$ref->ipro_id} = 1;
            }
            else {
                next IPRO;
            }
            my $gene_ts = $transcript . ",". $gene;
            my $gname = $gene;
            if(defined($self->{_ts_flag})) {
                $gname = $gene_ts;
            }
#            next if(!defined($ipro->parid));
            next RECORD unless $i->setid =~ /match_part/x;
            if(defined($filter))
            {
                next RECORD unless $i->setid =~ /$filter/x;
            }
            my $s = $i->start_;
            my $e = $i->end;
            $dlen->{$i->name}{$gname} = $e - $s + 1; 
            my $gene_mutations = $self->{_mutations}->{$gt};
            if(exists($cache{$i->name}{$i->start_}{$i->end}) ) {
                next RECORD;
            }
            else {
                $cache{$i->name}{$i->start_}{$i->end} = 1;
            }
AACHG:      foreach my $aachg ( @{$gene_mutations->get_all_mutations()} ) {
                #
                my $orig = undef;
                my $pos = undef;
                my $new = undef;
                if($aachg =~ m/p\.([A-Z*]+)(\d+)([A-Z*]{0,1}|ins|del|insertion|deletion)/x) {
                    $orig = $1;
                    $pos = $2;
                    $new = $3;
                }
                else {
                    next AACHG;
                }
                if(!exists($self->{_mcounts}->{$gname}->{$i->name})) {
                    $self->{_mcounts}->{$gname}->{$i->name} = 0;
                }

                if(($pos >= $s) && ($pos <= $e)){
                    #print $g,"\t", $aachg, "\t", $p->{name},"\tYES\n";
                    $self->{_mcounts}->{$gname}->{$i->name} += 1;
                    $gene_mutations->add_domain($aachg,$i->name);
                }

            }

            } # next ipro
            
        }
    }

    # I'm wondering - does your connection go out of scope here
    # and automagically disconnect? explicitly calling 
    # disconnect on the db handle blows up with UR saying there ain't
    # no data source.

    return 1;
}

# BROKEN! HA!?
sub get_mutation_counts {
    my ($self) = @_;
    return $self->{_mcounts};
}

# BROKEN! HA!?
sub get_mutation_counts_by_gene {
    my ($self,$gene) = @_;
    if(exists($self->{_mcounts}->{$gene})) {
        return $self->{_mcounts}->{$gene};
    }
    else {
        carp "no counts for $gene";
    }
    return;
}

sub get_gene_ids {
    my ($self) = @_;
    my @ids = sort keys %{$self->{_mcounts}};
    return \@ids;
}

sub get_all_domains
{
    my $self = shift;
    eval qq|
        use MPSampleData::DBI;
        use MPSampleData::Gene;
        use MPSampleData::Transcript;
        use MPSampleData::Protein;
        use MPSampleData::IproGeneTranscriptXref;
        use MPSampleData::IproResults;
    |;
    die $@ if $@;
    my $gene = shift;
    my $transcript = shift;

    my ($t) = MPSampleData::Transcript->search(transcript_name => $transcript);
    return unless $t;
    my @iprx = MPSampleData::IproGeneTranscriptXref->search(transcript_id => $t->transcript_id);
    my @domains;
    foreach my $ipr (@iprx)
    {
        # do stuff
        my $ipro = MPSampleData::IproResults->search(ipro_id => $ipr->ipro_id);
        while(my $i = $ipro->next)
        {
            push(@domains, $i->name);
        }
    }
    @domains = uniq @domains;
    return @domains;
}

sub get_lengths_from_db {
    eval qq|
        use MPSampleData::DBI;
        use MPSampleData::Gene;
        use MPSampleData::Transcript;
        use MPSampleData::Protein;
        use MPSampleData::IproGeneTranscriptXref;
        use MPSampleData::IproResults;
    |;
    die $@ if $@;

    my ($self) = @_;
    my @tgs = sort keys %{$self->{_mutations}};
    my %tnames = map { (split(/,/x,$_))[0] => $_  } @tgs;
    my %lengths;
    foreach my $tname (sort keys %tnames) {
        my ($t) = MPSampleData::Transcript->search({transcript_name => $tname});
        if(defined($t))
        {
            my ($p) = MPSampleData::Protein->search({transcript_id => $t->transcript_id});
            if($p) {
            $lengths{$tnames{$tname}} = length($p->amino_acid_seq) or 0; 
            }
            else {
                carp "no pep for $tname";
            }
        }
        else
        {
            $lengths{$tnames{$tname}} = 0;
            print STDERR "no length found for ", $tname, 
                         " ", $tnames{$tname},"\n";
        }
    }
    # I'm wondering - does your connection go out of scope here
    # and automagically disconnect? explicitly calling 
    # disconnect on the db handle blows up with UR saying there ain't
    # no data source.
    return \%lengths;
}

sub add_aachange {
    my $self = shift;
    my %args = @_;
    my $gt = undef;
    if(exists($args{'-gt'})) {
        $gt = $args{'-gt'};
    }
    elsif(exists($args{'-gene'}) && exists($args{'-transcript'}) ) {
        $gt = $args{'-transcript'}. "," . $args{'-gene'};
    }
    else {
        carp "wrong args!";
        return;
    }

    my $sdm = new SnpDom::Mutation('gt' => $gt);
    if(exists($self->{_mutations}->{$gt})) {
        my $sdm = $self->{_mutations}->{$gt};
        $sdm->add_mutation($args{'-aachange'});
        $self->{_mutations}->{$gt} = $sdm;
    }
    else {
        my $sdm = new SnpDom::Mutation('gt' => $gt);
        $sdm->add_mutation($args{'-aachange'});
        $self->{_mutations}->{$gt} = $sdm;
    }
    return 1;
}

sub add_mutation
{
    my ($self,$gene,$transcript,$aa_change) = @_;
    if(exists($self->{_mutations}->{$transcript.",".$gene})) {
        my $sdm = $self->{_mutations}->{$transcript.",".$gene};
        $sdm->add_mutation($aa_change);
        $self->{_mutations}->{$transcript.",".$gene} = $sdm;
    }
    else {
        my $sdm = new SnpDom::Mutation('gt' => $transcript.",".$gene);
        $sdm->add_mutation($aa_change);
        $self->{_mutations}->{$transcript.",".$gene} = $sdm;
    }
    return 1;
}

1;
# $Id: SnpDom.pm 38920 2008-09-22 19:58:03Z josborne $
