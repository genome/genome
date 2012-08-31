package Genome::Model::Tools::WuBlast::Parse;

use strict;
use warnings;

use Genome;
use IO::File;
use Bio::SearchIO;


class Genome::Model::Tools::WuBlast::Parse {
    is  => 'Command',
    has => [
        blast_outfile => {
            type     => 'String',
            is_input => 1,
            doc      => 'File name of a blast output file including path',
        },
    ],
    has_optional => [
        parse_outfile => {
            type     => 'String',
            is_input => 1,
            doc      => 'File name of the blast-output parsing file including path',
        },
        percent_threshold => {
            type     => 'Number',
            is_input => 1,
            doc      => 'percent_identity value to pass',
        },
        length_threshold  => {
            type     => 'Number',
            is_input => 1,
            doc      => 'match sequence length value to pass',
        },

    ],
};


sub help_brief {
    "Wrapper for parsing blast output" 
}


sub help_detail {  
    return <<EOS
This should work on all WuBlast reports, including blastn, blastp, blastx as 
well as blat with blast-format output. For normal parsing purpose, this will
write output to parse_outfile option, otherwise it will return BioPerl High
Scoring Pairs (HSP) object arrays to other scripts/modules
EOS
}

    
sub execute {
    my $self = shift;
    my $out_file = $self->blast_outfile;

    unless (-s $out_file) {
        $self->error_message('Blast output file '.$out_file.' does not exist');
        return;
    }

    my $p_thresh = $self->percent_threshold || '';
    my $l_thresh = $self->length_threshold  || '';

    my $io = Bio::SearchIO->new(-file => $out_file, -format => 'blast');
    my @hsps = ();

    while (my $result = $io->next_result) {
	    while (my $hit = $result->next_hit) {
	        while (my $hsp = $hit->next_hsp) {
		        $hsp->hit->seqlength($hit->length);
		        push @hsps, $hsp;
            }
        }
    }    

    @hsps = grep{$_->length > $l_thresh} @hsps if $l_thresh;
    @hsps = grep{$_->percent_identity > $p_thresh} @hsps if $p_thresh;

    #unlink $output if $output =~ /^\/tmp\//;
    if ($self->parse_outfile) {
        my $parse_file = $self->parse_outfile;
        my $out_io = IO::File->new(">$parse_file") 
            or ($self->error_message("can't write to $parse_file") and return);
        _tabular_output($out_io, @hsps);
        $out_io->close;
    }
    
    return \@hsps;
}


sub _tabular_output {
    my ($io, @hsps) = @_;

    my $out  = "\n        query_id  q_start   q_end   q_length          hit_id        h_start   h_end  h_length  str  score  percent\n";
       $out .= '='x116;
       $out .= "\n";
    
    while (my $hsp = shift @hsps) {
	    $out .= sprintf("%16s%9d%9d%9d%22s%9d%9d%9d%5s%7d    %.1f", $hsp->query->seq_id, $hsp->query->start, $hsp->query->end,$hsp->query->seqlength, $hsp->hit->seq_id, $hsp->hit->start, $hsp->hit->end, $hsp->hit->seqlength, $hsp->query->strand, $hsp->score, $hsp->percent_identity);
	    $out .= "\n";
    }
    $io->print($out);
    return 1;
}

1;

#$HeadURL$
#$Id$
