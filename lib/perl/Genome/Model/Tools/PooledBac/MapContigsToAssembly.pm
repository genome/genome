package Genome::Model::Tools::PooledBac::MapContigsToAssembly;

use strict;
use warnings;

use Genome;

use Data::Dumper;
use Genome::Model::Tools::Pcap::Ace;
use Genome::Model::Tools::Pcap::Phd;
use List::Util qw(max min);

class Genome::Model::Tools::PooledBac::MapContigsToAssembly {
    is => 'Command',
    has => 
    [        
        pooled_bac_dir =>
        {
            type => 'String',
            doc => "Pooled BAC Assembly Directory",    
        },
        ace_file_name =>
        {
            type => 'String',
            doc => "Ace file containing pooled bac contigs"
        },
        project_dir =>
        {
            type => 'String',
            doc => "output dir for separate pooled bac projects"        
        },
        contig_map_file =>
        {
            type => 'String',
            is_optional => 1,
            doc => "this file contains a list of contigs and where they map to",
        },
        percent_overlap => 
        {
            type => 'Number',
            is_optional => 1,
            default_value => 50,
            doc => "this is the percent overlap, default is 50%",
        },
        percent_identity =>
        {
            type => 'Number',
            is_optional => 1,
            default_value => 85,
            doc => "this is the percent identity, default is 85%",
        },
    ]
};

sub help_brief {
    "Move Pooled BAC assembly into separate projects"
}   

sub help_synopsis { 
    return;
}
sub help_detail {
    return <<EOS 
    Move Pooled BAC Assembly into separate projects
EOS
}

############################################################
sub execute { 
    my $self = shift;
    print "Finding Matching Contigs...\n";
    my $pooled_bac_dir = $self->pooled_bac_dir;
    my $project_dir = $self->project_dir;
    my $blastfile = $project_dir."/bac_region_db.blast";
    $self->error_message("$blastfile does not exist") and die unless (-e $blastfile);
    my $out = Genome::Model::Tools::WuBlast::Parse->execute(blast_outfile => $blastfile);   
    $self->error_message("Failed to parse $blastfile") and die unless defined $out;

    my $percent_overlap = $self->percent_overlap / 100;
    my $percent_identity = $self->percent_identity / 100;

    my $ace_file = $pooled_bac_dir.'/consed/edit_dir/'.$self->ace_file_name;
    $self->error_message("Ace file $ace_file does not exist") and die unless (-e $ace_file);    

    my $ut = Genome::Model::Tools::PooledBac::Utils->create;
    $self->error_message("Genome::Model::Tools::PooledBac::Utils->create failed.\n") unless defined $ut;

    my $contig_map_file = $self->contig_map_file || "CONTIG_MAP";
    $contig_map_file = $project_dir.'/'.$contig_map_file;
    my $contig_hsps = $ut->get_matching_contigs_list($out->{result});$out=undef;

    my %contig_map;
    my $contig_names = $ut->get_contig_names($ace_file);
    foreach my $contig (@{$contig_names})
    {
        $contig_map{$contig} = {maps_to => 'orphan_project', module => 'MapContigsToAssembly'};        
    }
    
    # Go thru HSPs, and collect regions hit for each contig
    # Bac = HIT
    # Contig = QUERY
    CONTIG: foreach my $contig_hsp ( @$contig_hsps ) {
        my %bacs;
        HSP: for my $hsp ( @$contig_hsp ) {
            # check that the hsp is above our thresholds
            #next unless $hsp->{HSP_LENGTH} >= 1000; # param??
            next unless ( $hsp->{IDENTICAL} / $hsp->{HSP_LENGTH} ) >= $percent_identity;
            
            # set the coverage on contig
            my $bac_name = $hsp->{HIT_NAME};
            unless ( exists $bacs{$bac_name} ){ 
                $bacs{$bac_name} = [];
            }

            my ($query_start, $query_end) = ( $hsp->{QUERY_START} < $hsp->{QUERY_END} ) 
            ? ( $hsp->{QUERY_START}, $hsp->{QUERY_END} )
            : ( $hsp->{QUERY_END}, $hsp->{QUERY_START} );
            print join(' ', $bac_name, $hsp->{QUERY_NAME}, $query_start, $query_end, $hsp->{QUERY_LENGTH})."\n";

            # see if this hsp covers an existing overlap, and combine it
            for my $overlap ( @{$bacs{$bac_name}} ) {
                next if $query_start > $overlap->[1] 
                    or $query_end < $overlap->[0];
                $overlap->[0] = $query_start if $query_start < $overlap->[0];
                $overlap->[1] = $query_end if $query_end > $overlap->[1];
                next HSP; # done, go to next HSP
            }

            # add the overlap
            push @{$bacs{$bac_name}}, [
                $query_start, $query_end,
            ];
        }

        next CONTIG unless %bacs;

        # Calc coverage for contig to each bac
        my %bac_coverage;
        BAC: for my $bac_name ( keys %bacs ) {
            # Get the coverage of overlaps
            my $coverage = 0;
            for my $overlap ( @{$bacs{$bac_name}} ) {
                $coverage += $overlap->[1] - $overlap->[0] + 1;
            }
            next BAC unless ( $coverage / $contig_hsp->[0]->{QUERY_LENGTH} ) >= $percent_overlap;
            $bac_coverage{ $bac_name } = $coverage;
        }

        # If we have a best bac, assign it to this contig
        next CONTIG unless %bac_coverage;
        my ($best_bac) = sort { $bac_coverage{$a} <=> $bac_coverage{$b} } keys %bac_coverage;
        $contig_map{ $contig_hsp->[0]->{QUERY_NAME} }->{maps_to} = $best_bac;
    }

    $ut->write_contig_map(\%contig_map, $contig_map_file);

    return 1;
}

1;

