package Genome::Model::Tools::PooledBac::AddLinkingContigs;

use strict;
use warnings;

use Genome;
use Genome::Model::Tools::Pcap::Ace;
use Genome::Model::Tools::Pcap::Phd;
use Genome::Sys;
use Genome::Model::Tools::PooledBac::Utils;

class Genome::Model::Tools::PooledBac::AddLinkingContigs {
    is => 'Command',
    has => 
    [
        contig_map_file =>
        {
                type => 'String',
                is_optional => 1,
                doc => "this file contains a list of contigs and where they map to",
        },
        project_dir =>
        {
            type => 'String',
            is_optional => 0,
            doc => "output dir for separate pooled bac projects"        
        }
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
    print "Adding Linking Contigs...\n";
    my $project_dir = $self->project_dir;
    my $ut = Genome::Model::Tools::PooledBac::Utils->create;
    my $contig_map_file = $self->contig_map_file || "CONTIG_MAP";
    $contig_map_file = $project_dir.'/'.$contig_map_file;
    $self->error_message("Contig map file, $contig_map_file, does not exist.\n") and die unless (-e $contig_map_file);
    my $contig_map = $ut->open_contig_map($contig_map_file);
    my $match_list;
    my $orphan_list;
    Genome::Sys->create_directory($project_dir.'/reports');
    my $fh = Genome::Sys->open_file_for_writing( $project_dir.'/reports/contigs_that_link_to_matching_contigs' );
    ($match_list, $orphan_list) = $ut->create_match_and_orphan_lists($contig_map);
    foreach my $orphan (keys %{$orphan_list})
    {
        my ($sc_num, $ct_num) = $orphan =~ /Contig(\d+)\.(\d+)/;
        next unless (defined $sc_num && defined $ct_num);        
        my $pre_ctg = "Contig$sc_num.".($ct_num -1);
        my $aft_ctg = "Contig$sc_num.".($ct_num+1);
        if(exists $match_list->{$pre_ctg} && 
           exists $match_list->{$aft_ctg} &&
           $match_list->{$pre_ctg}{maps_to} eq $match_list->{$aft_ctg}{maps_to})        
        {
            $contig_map->{$orphan}->{maps_to} = $match_list->{$pre_ctg}->{maps_to};
            $contig_map->{$orphan}->{module} = 'AddLinkingContigs';
        }
        if(exists $match_list->{$pre_ctg})
        {
            $fh->print("$orphan links to $pre_ctg in $match_list->{$pre_ctg}{maps_to}\n");
        }
        if(exists $match_list->{$aft_ctg})
        {
            $fh->print("$orphan links to $aft_ctg in $match_list->{$aft_ctg}{maps_to}\n");
        }    
    }
    $fh->close;
    $ut->write_contig_map($contig_map,$contig_map_file);
    return 1;
}

1;
