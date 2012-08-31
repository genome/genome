package Genome::Model::Tools::Annotate::Adaptor::Sniper;
#adukes: 1/29/10, I started using this to adapt the snps.indelpe.filtered file in refalign pipe to annotation input.  We might want to rename this module to something not Sniper specific like FourColumnAdaptor or something.  This would require some more naive column printing after the first 4 columns, but I think that's okay.

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Annotate::Adaptor::Sniper {
    is => 'Genome::Model::Tools::Annotate',
    has => [
        somatic_file => {
            is  => 'String',
            is_input  => 1,
            doc => 'The somatic file output from sniper to be adapted',
        },
        output_file => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => "Store output in the specified file instead of sending it to STDOUT."
        },
        skip_if_output_present => {
            is => 'Boolean',
            is_input => 1,
            is_optional => 1,
            default => 0,
            doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
        },
    ],
};

sub help_brief {
    "Converts somatic sniper output into a gmt annotate transcript-variants friendly input format",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt annotate adaptor sniper ...    
EOS
}

sub help_detail {                           
    return <<EOS 
Converts somatic sniper output into a gmt annotate transcript-variants friendly input format
EOS
}

# For now assume we have a somatic file already made from the normal and tumor bam files... separate tool for this perhaps
sub execute {
    my $self = shift;

    unless (-s $self->somatic_file) {
        $self->error_message("somatic_file: " . $self->somatic_file . " does not exist or has no size");
        die;
    }
    
    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    my $somatic_fh = IO::File->new($self->somatic_file);

    # establish the output handle for the transcript variants
    my $output_fh;
    if (my $output_file = $self->output_file) {
        $output_fh = $self->_create_file($output_file);
    }
    else {
        $output_fh = 'STDOUT';
    }

    my @return;   
    while (my $line=$somatic_fh->getline) {
        chomp $line;
        if($line =~ m/\*/) {
            #INDELOL!!!
            @return = $self->parse_indel_line($line);
            for my $ret (@return) {
                $output_fh->print(join("\t",@{$ret}) . "\n");
            }
       }
       else { #CONSNPTION FIT!
            my ($chr, $start, $ref_base, $variant_base, $somatic_score, $consensus_quality, $snp_quality, $max_map_q, $depth_tumor, $depth_normal) = split("\t", $line);
            $chr =~ s/chr//;
            $chr =~ s/M[^T]?$/MT/;
            $output_fh->print("$chr\t$start\t$start\t$ref_base\t$variant_base\tSNP\t$somatic_score\t$consensus_quality\t$snp_quality\t$max_map_q\t$depth_tumor\t$depth_normal\n");
        }
    }
    
    $output_fh->close;
    return 1;
}


sub parse_indel_line {
    my $self=shift;
    my ($line) = @_;
    my @return;
    #$self->status_message("(SEARCHING FOR: $line)");

    my %indel1;
    my %indel2;
    my ($chr,
        $start_pos,
        $star, 
        $somatic_score,
    );
    my @rest_of_fields;
    ($chr,
        $start_pos,
        $star, 
        $somatic_score,
        $indel1{'sequence'},
        $indel2{'sequence'}, 
        $indel1{'length'},
        $indel2{'length'},
        @rest_of_fields
    ) = split /\s+/, $line; 
    $chr =~ s/chr//;
    $chr =~ s/M[^T]?$/MT/;
    my @indels;
    push(@indels, \%indel1);
    push(@indels, \%indel2);
    for my $indel(@indels) {
        
        if ($indel->{'sequence'} eq '*') { next; }
        my $hash;
        my $stop_pos;
        my $start;
        if($indel->{'length'} < 0) {
            #it's a deletion!
            $hash->{variation_type}='DEL';
            $start= $start_pos+1;
            $stop_pos = $start_pos + abs($indel->{'length'});
            $hash->{reference}=$indel->{'sequence'};
            $hash->{variant}=0;
        }
        else {
            #it's an insertion
            $hash->{variation_type}='INS';
            $start=$start_pos;
            $stop_pos = $start_pos+1;
            $hash->{reference}=0;
            $hash->{variant}=$indel->{'sequence'};

        }

        $hash->{chromosome}=$chr;
        $hash->{start}=$start;
        $hash->{stop}=$stop_pos;
        #$hash->{num_reads}=$num_reads_across;
        push @return, [@$hash{"chromosome","start","stop","reference","variant","variation_type"}, $somatic_score, @rest_of_fields];
    }
    return @return;
}
1;

