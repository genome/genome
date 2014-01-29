package Genome::Model::Tools::PooledBac::GenerateReports;

use strict;
use warnings;

use Genome;
use Genome::Model::Tools::Pcap::Ace;
use Genome::Model::Tools::Pcap::Phd;
use Genome::Model::Tools::PooledBac::Utils;
use Genome::Sys;
use List::Util qw(max min);

class Genome::Model::Tools::PooledBac::GenerateReports {
    is => 'Command',
    has => 
    [        
        project_dir =>
        {
            type => 'String',
            is_optional => 0,
            doc => "output dir for separate pooled bac projects"        
        },
        contig_map_file =>
        {
            type => 'String',
            is_optional => 1,
            doc => "this file contains a list of contigs and where they map to",
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

sub print_matching_contigs_report
{
    my ($self, $list, $report_name) = @_;
    #print list
    my $fh = IO::File->new('>'.$report_name);
    $self->error_message("Failed to open $report_name for writing.") and die unless defined $fh;
    foreach my $result (@$list)
    {
        foreach my $res (@{$result})
        {
            $fh->print ($res->{QUERY_NAME}, ' ',$res->{HIT_NAME},' ',$res->{HSP_LENGTH},' ',$res->{IDENTICAL}/$res->{HSP_LENGTH},"\n");        
        }    
    }
}

sub print_multiple_hits_report
{
    my ($self, $list, $report_name) = @_;
    #print list
    my $fh = IO::File->new('>'.$report_name);
    $self->error_message("Failed to open $report_name for writing.")  and die unless defined $fh;
    foreach my $result (@$list)
    {
        next unless (@{$result} > 1);
        
        for(my $i=0;$i<@{$result};$i++)
        {
            my $res = $result->[$i];
            if($i==0)
            {
                $fh->print ($res->{QUERY_NAME});        
            }
            else
            {
                my $name = $res->{QUERY_NAME};
                $name =~ s/./ /g;
                $fh->print ($name);
            }
            $fh->print ("\t",$res->{HIT_NAME},' ',$res->{HSP_LENGTH},' ',$res->{IDENTICAL}/$res->{HSP_LENGTH},"\n");
        }    
    }
}

sub print_close_match_report
{
    my ($self, $list, $report_name) = @_;
    #print list
    my $l_pcutoff = 0.05;#length percent difference cutoff, if there is a 5% or less diff in length
    my $m_pcutoff = 0.05;#matching percent difference cutoff, if there is a 5% or less diff in identitiy
    my $fh = IO::File->new('>'.$report_name);
    $self->error_message("Failed to open $report_name for writing.")  and die unless defined $fh;
    foreach my $result (@$list)
    {
        next unless (@{$result} > 1);
        my $res0 = $result->[0];
        my $res1 = $result->[1];
        my $max_length = max($res0->{HSP_LENGTH},$res1->{HSP_LENGTH});
        my $l_pdiff = abs (($res0->{HSP_LENGTH}-$res1->{HSP_LENGTH})/$max_length);
        my $pid0 = $res0->{IDENTICAL}/$res0->{HSP_LENGTH};
        my $pid1 = $res1->{IDENTICAL}/$res1->{HSP_LENGTH};
        my $max_id = max($pid0,$pid1);
        my $m_pdiff = abs(($pid0-$pid1)/$max_id);
        next unless (($l_pdiff < $l_pcutoff) && ($m_pdiff < $m_pcutoff));
        for(my $i=0;$i<@{$result};$i++)
        {
            my $res = $result->[$i];
            if($i==0)
            {
                $fh->print ($res->{QUERY_NAME});        
            }
            else
            {
                my $name = $res->{QUERY_NAME};
                $name =~ s/./ /g;
                $fh->print ($name);
            }
            $fh->print ("\t",$res->{HIT_NAME},' ',$res->{HSP_LENGTH},' ',$res->{IDENTICAL}/$res->{HSP_LENGTH},"\n");
        }    
    }
}

sub print_orphan_contigs_report
{
    my ($self, $orphan_contigs, $report_name) = @_;
    my $fh = IO::File->new('>'.$report_name);
    $self->error_message("Failed to open $report_name for writing.")  and die unless defined $fh;
    foreach my $name (@{$orphan_contigs})
    {
        $fh->print ($name,"\n");
    } 
}

sub print_used_contigs_report
{
    my ($self, $contig_map, $report_name) = @_;
    my $fh = IO::File->new('>'.$report_name);
    $self->error_message("Failed to open $report_name for writing.") and die unless defined $fh;
    foreach my $contig_name (keys %{$contig_map})
    {
        my $hash_ref = $contig_map->{$contig_name};
        print $fh "$contig_name $hash_ref->{maps_to}\n";    
    }  
}

sub combined_used_and_orphan_lists {
    my $self = shift;

    my $report_dir = $self->project_dir.'/reports';
    my $out_file = $report_dir.'/complete_contig_list_with_orphan_contigs';

    unlink $out_file;

    my $fh_out = Genome::Sys->open_file_for_writing( $out_file );
    my $fh_used_in = Genome::Sys->open_file_for_reading( $report_dir.'/complete_contig_list' );
    my $fh_orph_in = Genome::Sys->open_file_for_reading( $report_dir.'/orphan_contigs' );

    while ( my $line = $fh_used_in->getline ) {
	$fh_out->print( $line );
    }
    while ( my $line = $fh_orph_in->getline ) {
	chomp $line;
	$fh_out->print( $line." orphan \n" );
    }

    $fh_out->close;
    $fh_used_in->close;
    $fh_orph_in->close;

    return 1;
}
############################################################
sub execute { 
    my $self = shift;
    print "Generating Reports...\n";
    my $project_dir = $self->project_dir;
    my $blastfile = $project_dir."/bac_region_db.blast";
    my $reports_dir = $project_dir."/reports/";
    $self->error_message("Failed to create directory $reports_dir")  and die unless Genome::Sys->create_directory($reports_dir);
    #`mkdir -p $reports_dir`;
    my $out = Genome::Model::Tools::WuBlast::Parse->execute(blast_outfile => $blastfile, parse_outfile => $reports_dir."blast_report");
    $self->error_message("Failed to parse $blastfile")  and die unless defined $out; 

    my $ut = Genome::Model::Tools::PooledBac::Utils->create;
    $self->error_message("Genome::Model::Tools::PooledBac::Utils->create failed.\n") unless defined $ut;

    my $list = $ut->get_matching_contigs_list($out->{result});$out=undef;
    $self->print_matching_contigs_report($list, $reports_dir."matching_contigs");
    $self->print_close_match_report($list,$reports_dir."ambiguous_matching_contigs");
    $self->print_multiple_hits_report($list,$reports_dir."contigs_with_multiple_hits");
    
    

    my $contig_map_file = $self->contig_map_file || "CONTIG_MAP";
    $contig_map_file = $project_dir.'/'.$contig_map_file;    
    my $contig_map = $ut->open_contig_map($contig_map_file);
    
    my ($match_list, $orphan_list) = $ut->create_match_and_orphan_lists($contig_map);
    my $orphan_contig_names = [keys %{$orphan_list}];
    $self->print_orphan_contigs_report($orphan_contig_names, $reports_dir."orphan_contigs");
    $self->print_used_contigs_report($match_list, $reports_dir."complete_contig_list");

    unless ( $self->combined_used_and_orphan_lists() ) {
	$self->debug_message("Failed to create complete_contig_list_with_orphan_contigs list");
    }
    return 1;
}



1;
