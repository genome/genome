package Genome::Model::Tools::Ber::AmgapPrepareBER;

use strict;
use warnings;

use Genome;
use Command;

use Carp;

use BAP::DB::CodingGene;
use BAP::DB::Sequence;
use BAP::DB::SequenceSet;

use Carp;
use Getopt::Long;
use Cwd;
use Pod::Usage;
use English;
use IO::File;

# FIXME: add output directory parameter so we can execute this wherever we want to.
UR::Object::Type->define(
			 class_name => __PACKAGE__,
			 is         => 'Command',
			 has        => [
					'locus_tag'       => {
							      is => 'String',
							      doc => "Locus tag for project, containing DFT/FNL postpended",
							     },
					'sequence_set_id' => {
							      is  => 'Integer',
							      doc => "Sequence set id in MGAP database",
							     },
					'phase'           => {
							      is  => 'Interger',
							      doc => "Phase of protein to dump from Oracle.  Default is phase 5. ",
							      default => 5,
							      is_optional => 1,
							     },

				       ]
			);


sub help_brief
  {
    "Tool for creating config files needed for BER product naming pipeline";
  }

sub help_synopsis
  {
    return <<"EOS"
      Tool for creating config files needed for BER product naming pipeline.  Creates the 4 config files and dumps the gene name from oracle.  Then converts the 4 config files from unix2dos format. 
EOS
  }

sub help_detail
  {
    return <<"EOS"
Tool for creating config files needed for BER product naming pipeline.  Creates the 4 config files and dumps the gene name from oracle.  Then converts the 4 config files from unix2dos format.
EOS
  }


sub execute {
    my $self = shift;
    my $locus_tag      = $self->locus_tag;

    # FIXME: files opened here require you to be in the correct directory.
    # this probably isn't a good idea.

    my $asm_feature_fh = IO::File->new();
    my $asm_feature    = qq{$locus_tag\_asm_feature};
    $asm_feature_fh->open("> $asm_feature") 
        or croak "Can not open new asm_feature file, from AmgapPrepareBER.pm: $OS_ERROR\n";

    print $asm_feature_fh qq{asmbl_id\tend3\tend5\tfeat_name\tfeat_type}, "\n";

    my $ident2_fh = IO::File->new();
    my $ident2    = qq{$locus_tag\_ident2};
    $ident2_fh->open (">$ident2")
        or croak "Can not open new ident2 file, from AmgapPrepareBER.pm: $OS_ERROR\n";

    print $ident2_fh qq{complete\tfeat_name\tlocus}, "\n";

    my $asmbl_data_fh = IO::File->new();
    my $asmbl_data    = qq{$locus_tag\_asmbl_data};
    $asmbl_data_fh->open (">$asmbl_data")
        or croak "Can not open new asmbl_data file, from AmgapPrepareBER.pm: $OS_ERROR\n";

    print $asmbl_data_fh qq{id\tname\ttype}, "\n";

    my $stan_fh = IO::File->new();
    my $stan = qq{$locus_tag\_stan};
    $stan_fh->open(">$stan")
        or croak "Can not open new stan file, from AmgapPrepareBER.pm: $OS_ERROR\n";

    print $stan_fh qq{asmbl_data_id\tasmbl_id\tiscurrent}, "\n";

    my $phase           = $self->phase;
    my $sequence_set_id = $self->sequence_set_id;

    my $gene_phase_name = join '', 'p', $phase, '_hybrid';

    $phase = "phase_$phase";

    my $sequence_set = BAP::DB::SequenceSet->retrieve($sequence_set_id);

    my @sequences = $sequence_set->sequences();

    @sequences = sort { $a <=> $b } @sequences;

    my $asmbl_id_count = 1;
    my $ident2_count   = 1;
    my $iscurrent      = 1;
    my $genecountcheck = 1;

# some fixes to filter out the rrna hits...
    SEQUENCE: foreach my $i (0..$#sequences) {
        my $sequence      = $sequences[$i];
        my $sequence_name = $sequence->sequence_name();
        my @coding_genes  = $sequence->coding_genes($phase => 1);
        @coding_genes     = sort { $a <=> $b } @coding_genes;

        unless (
            (@coding_genes > 0)
        ) {
            next;
        }
        my %method_fixup = (
            'glimmer2' => 'Glimmer2',
            'glimmer3' => 'Glimmer3',
            'genemark' => 'GeneMark',
            'blastx'   => 'blastx',
        );

        $genecountcheck++;

        GENE: foreach my $i (0..$#coding_genes) {
            my $coding_gene = $coding_genes[$i];
            my $gene_name = $coding_gene->gene_name();

            # rrna hits don't get pulled.  other dead genes or anything else tagged,
            # may pass.
            my @tags = $coding_gene->tag_names();
            foreach my $tag (@tags) {
                if(($tag->tag_name eq 'Dead')) {
                    #print $coding_gene->gene_name()."\t".$tag->tag_name."\t".$tag->tag_value."\n";
                    next GENE;
                }
            }

            if ($phase ne 'phase_0') { 
                $gene_name = rename_gene($gene_name,$gene_phase_name);
            }

            my $start       = $coding_gene->start();
            my $end         = $coding_gene->end();
            my $source      = $sequence_name;
            my $method      = $coding_gene->source();
            my $complete    = qq{};

            if ($phase ne 'phase_0') { $method = $gene_phase_name; }

            if (exists($method_fixup{$method})) {
                $method = $method_fixup{$method};
            }

            print $asm_feature_fh qq{$asmbl_id_count\t$end\t$start\t$gene_name\tORF}, "\n";
            printf $ident2_fh qq{$complete\t$gene_name\t$locus_tag%04d\n}, $ident2_count;
            $ident2_count++;
        }
        print $stan_fh qq{$asmbl_id_count\t$asmbl_id_count\t$iscurrent}, "\n";
        print $asmbl_data_fh qq{$asmbl_id_count\tContig\tcontig}, "\n";
        $asmbl_id_count++;
    }

    unless ($asmbl_id_count == $genecountcheck) {

        warn qq{\n\nThe gene counts and the loop counts do not match, from AmgapPrepareBER.pm :$asmbl_id_count\t$genecountcheck:  $OS_ERROR\n\n};

    }

    $asm_feature_fh->close();
    $ident2_fh     ->close();
    $stan_fh       ->close();
    $asmbl_data_fh ->close();

    my @file2convert = ($asm_feature,$ident2,$asmbl_data,$stan);

    foreach my $file (@file2convert) {

        #unless (system(qq{/usr/bin/unix2dos -u -o -p '$file'})==0) {
        unless (system(qq{/usr/bin/recode ..pc '$file'})==0) {

            #die "unix2dos died for '$file'  , from AmgapPrepareBER.pm\n";
            die "recode died while reformating unix to a dos file format for,  '$file'  , from AmgapPrepareBER.pm\n";

        }
    }

    return 1;

}

sub rename_gene {

    my ($old_name, $gene_phase_name) = @_;

    unless (defined($old_name)) { croak "encountered undefined gene name, from AmgapPrepareBER.pm\n"; }

    my @cols = split /\./, $old_name;

    splice @cols, -1, 1, $gene_phase_name, $cols[$#cols];

    return join '.', @cols;

}



1;

# $Id$
