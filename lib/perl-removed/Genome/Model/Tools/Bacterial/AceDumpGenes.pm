package Genome::Model::Tools::Bacterial::AceDumpGenes;

use strict;
use warnings;
use Genome;

use BAP::DB::DBI;
use BAP::DB::SequenceSet;
use BAP::DB::CodingGene;
use BAP::DB::Sequence;
use Data::Dumper;
use IO::File;
use IO::Handle;
use Carp;

class Genome::Model::Tools::Bacterial::AceDumpGenes {
    is => 'Command',
    has => [
        sequence_set_id => { 
            is => 'Integer',
            doc => "sequence set id of genes to dump out",
        },
    ],
    has_optional => [
        phase => {
            is => 'Integer',
            doc => "specify which phase of gene merging to dump from",
            default => 5,
                 },
        dev => {
            is => 'Boolean',
            doc => "use development database",
            default => 0,
               },
        ace_file => {
            is => 'Scalar',
            doc => "optional output file;defaults to STDOUT",
        },
    ],
};

sub help_brief {
    return "Used to dump ace file from the orcale database in BAP/MGAP";
}

sub help_detail {
    return <<EOS
This script is intented to be used with the HGMI/bacterial annotation pipeline (bap/mgap).  This script dumps data in ace format from the oracle database for a given sequence-set-id.
EOS
}

sub help_synopsis {
    return "gmt bacterial ace-dump-genes --sequence-set-id 123";
}

sub execute {
    my $self = shift;
    my $phase = $self->phase;
    my $sequence_set_id = $self->sequence_set_id;
    my $dev_flag = $self->dev;

    my $fh;
    if(defined($self->ace_file)) {
        $self->debug_message("writing ace file to". $self->ace_file);
        $fh = IO::File->new(">".$self->ace_file);
    }
    else
    {
        $self->debug_message("ace data will be written to stdout");
        $fh = IO::Handle->new;
        $fh->fdopen(fileno(STDOUT),"w");
    }
    if ($dev_flag) { $BAP::DB::DBI::db_env = 'dev'; }

    my $gene_phase_name = join '', 'p', $phase, '_hybrid';

    $phase = "phase_$phase";

    my $sequence_set = BAP::DB::SequenceSet->retrieve($sequence_set_id);

    unless($sequence_set) {
        $self->debug_message("can't retrieve anything for sequence set $sequence_set_id");
        croak;
    }


    my $software_version = $sequence_set->software_version();
    my $data_version     = $sequence_set->data_version();

    my $svn_version = '$Revision: 60465 $';
    $svn_version    =~ s/\D+//g;

    print $fh qq{//BAP/MGAP Version $software_version ($svn_version)}, "\n";
    print $fh qq{//Data Version $data_version}, "\n";

    if ($software_version eq 'unknown') {
        $software_version = 'v0.0';
    }
    else {
        $software_version = join('v', $software_version);
    }
    $svn_version      = join('svn', $svn_version);

    if ($data_version eq 'unknown') {
        $data_version = '700101';
    }
    else {
        my ($m, $d, $y) = split(/\-/, $data_version);
        $data_version = join('', substr($y, 2, 2), $m, $d);
    }

    my @sequences = $sequence_set->sequences();

    foreach my $i (0..$#sequences) {

        my $sequence      = $sequences[$i];
        my $sequence_name = $sequence->sequence_name();

        my @coding_genes = $sequence->coding_genes($phase => 1);
        #@coding_genes = purge_dead_genes(@coding_genes);

        my @trna_genes   = $sequence->trna_genes();
        my @rna_genes    = $sequence->rna_genes();
        @rna_genes = grep { !$_->redundant() } @rna_genes;
        my @rfam_genes    = grep { $_->source =~ /rfam/i;    } @rna_genes;
        my @rnammer_genes = grep { $_->source =~ /rnammer/i; } @rna_genes;

        unless (
            (@coding_genes > 0) ||
            (@trna_genes > 0)   ||
            (@rfam_genes > 0)   ||
            (@rnammer_genes > 0)
        ) {
            next;
        }

        print $fh qq{Sequence $sequence_name}, "\n";
        print $fh qq{Database $software_version $data_version $svn_version}, "\n\n";
        print $fh qq{Sequence : "$sequence_name"}, "\n";

        foreach my $gene (@coding_genes, @trna_genes, @rna_genes) {

            my $gene_name = $gene->gene_name();

            if ($gene->isa('BAP::DB::CodingGene') && ($phase ne 'phase_0')) {
                $gene_name = rename_gene($gene_name,$gene_phase_name);
            }

            my $start     = $gene->start();
            my $end       = $gene->end();

            # AceDB Ace format specifies minus strand with (start > end)
            # DB should always be (start < end) regardless of strand (but don't trust it)
            if ($gene->strand() < 0) {
                if ($start < $end) {
                    ($start, $end) = ($end, $start);
                }
            }
            else {
                if ($start > $end) {
                    ($start, $end) = ($end, $start);
                }

            }

            print $fh qq{Subsequence "$gene_name" $start $end}, "\n";

        }

        if (@coding_genes) { print $fh "\n"; }

        my %method_fixup = (
            'glimmer2' => 'Glimmer2',
            'glimmer3' => 'Glimmer3',
            'genemark' => 'GeneMark',
            'blastx'   => 'blastx',
        );

        foreach my $i (0..$#coding_genes) {

            my $coding_gene = $coding_genes[$i];

            my $gene_name = $coding_gene->gene_name();

            if ($phase ne 'phase_0') { $gene_name = rename_gene($gene_name,$gene_phase_name); }

            my $start       = $coding_gene->start();
            my $end         = $coding_gene->end();
            my $source      = $sequence_name;
            my $method      = $coding_gene->source();

            if ($phase ne 'phase_0') { $method = $gene_phase_name; }

            my $length      = (abs($start - $end)) + 1;

            my $exon        = "1 $length";
            my $cds         = "1 $length";

            if (exists($method_fixup{$method})) {
                $method = $method_fixup{$method};
            }

            print $fh qq{Sequence : "$gene_name"}, "\n";
            print $fh qq{Source "$source"}, "\n";
            print $fh qq{Method "$method"}, "\n";
            print $fh qq{Source_Exons\t$exon}, "\n";
            print $fh qq{CDS\t$cds}, "\n";

            if ($coding_gene->missing_start()) {
                print $fh qq{Start_not_found}, "\n";
            }
            if ($coding_gene->missing_stop()) {
                print $fh qq{End_not_found}, "\n";
            }

            if (
                $coding_gene->blastp_evidence() &&
                $coding_gene->pfam_evidence()
            ) {
                print $fh qq{Protein_evidence "blastp and pfam"}, "\n";
            }
            elsif ($coding_gene->blastp_evidence()) {
                print $fh qq{Protein_evidence "blastp"}, "\n";
            }
            elsif ($coding_gene->pfam_evidence()) {
                print $fh qq{Protein_evidence "pfam"},"\n";
            }

            if ($phase eq 'phase_5' and ($coding_gene->gene_id == 6093590 or $coding_gene->gene_id == 6093636 or $coding_gene->gene_id == 6139539)) {
            }

            my @coding_gene_tags = $coding_gene->tag_names;
            if(@coding_gene_tags) {
                foreach my $tag (@coding_gene_tags) {
                    print $fh $tag->tag_name," \"", $tag->tag_value,"\"\n";
                }
            }

            unless ($i == $#coding_genes) {
                print $fh "\n";
            }

        }

        if (@trna_genes) { print $fh "\n"; }

        foreach my $i (0..$#trna_genes) {

            my $trna_gene  = $trna_genes[$i];
            my $gene_name  = $trna_gene->gene_name();
            my $source     = $sequence_name;
            my $score      = $trna_gene->score();
            my $codon      = $trna_gene->codon();
            my $aa         = $trna_gene->aa();
            my $aa_code    = substr $aa, 0, 1;

            my $method     = 'tRNAscan';
            my $remark     = "tRNA-$aa Sc=$score";
            my $transcript = qq{tRNA "$codon $aa $aa_code"}; 

            print $fh qq{Sequence : "$gene_name"}, "\n";
            print $fh qq{Source "$source"}, "\n";
            print $fh qq{Method "$method"}, "\n";
            print $fh qq{Remark "$remark"}, "\n";
            print $fh qq{Transcript $transcript}, "\n";

            unless ($i == $#trna_genes) {
                print $fh "\n";
            }
        }

        if (@rna_genes) { print $fh "\n"; }

        foreach my $i (0..$#rna_genes) {

            my $rna_gene   = $rna_genes[$i];
            my $gene_name  = $rna_gene->gene_name();
            my $source     = $sequence_name;
            my $score      = $rna_gene->score();
            my $desc       = $rna_gene->desc();
            my $model      = $rna_gene->acc();
            my $start      = $rna_gene->start();
            my $end        = $rna_gene->end();
            my $rfam_prod  = $rna_gene->rfam_prod();

            my $length     = (abs($start - $end))+1;

            my $exon       = "1 $length";

            my ($method, $remark);

            if ($rna_gene->source() =~ /rfam/i   ){
                $method = 'Rfam';
                $remark = "Predicted by Rfam ($model), score $score";
                $score  = "Rfam $score";
            }
            elsif ($rna_gene->source() =~ /rnammer/) {
                $method = 'RNAmmer';
                $remark = "Predicted by RNAmmer, score $score";
                $score =  "RNAmmer $score";
            }
            else {
                $method = $rna_gene->source();
                $remark = "Predicted by $method, score $score";
                $score  = "$method $score";
            }

            my $locus = "$desc";

            print $fh qq{Sequence : "$gene_name"}, "\n";
            print $fh qq{Source "$source"}, "\n";
            print $fh qq{Method "$method"}, "\n";

            if ($method eq 'RNAmmer') {

                print $fh qq{Source_Exons\t$exon}, "\n";

            }

            print $fh qq{Remark "$remark"}, "\n";
            print $fh qq{Locus $locus}, "\n";

            if ($method eq 'Rfam') {

                print $fh qq{Rfam_product "$rfam_prod"}, "\n";

            }

            print $fh qq{$score}, "\n";

            unless ($i == $#rna_genes) {
                print $fh "\n";
            }
        }

        # Sprinkle a little bit of whitespace in between contigs
        # (unless this is the last contig)
        unless ($i == $#sequences) {
            print $fh "\n";
        }

    }


    return 1;
}


sub rename_gene {

    my ($old_name,$gene_phase_name) = @_;

    unless (defined($old_name)) { croak "encountered undefined gene name"; }


    my @cols = split /\./, $old_name;

    splice @cols, -1, 1, $gene_phase_name, $cols[$#cols];

    return join '.', @cols;

}

sub purge_dead_genes {

    my @coding_genes = @_;
    my @nondead_coding_genes = ( );
    foreach my $gene_id (@coding_genes) {
        my $dead_tag = 0;
        my @tags = $gene_id->tag_names;
        foreach my $t (@tags) {
            if ( ($t->tag_name eq 'Dead') && ($t->tag_value eq 'rrna hit')) {
                $dead_tag = 1;
                last;
            }
        }
        unless ($dead_tag == 1) {
            push @nondead_coding_genes, $gene_id;
        }
    }
    @coding_genes = ( );
    @coding_genes = @nondead_coding_genes;

    return @coding_genes;

}




1;
