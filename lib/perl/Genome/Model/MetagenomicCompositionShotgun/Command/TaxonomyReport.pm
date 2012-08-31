package Genome::Model::MetagenomicCompositionShotgun::Command::TaxonomyReport;

use strict;
use warnings;

use Genome;

use File::Path;
use File::Find;

class Genome::Model::MetagenomicCompositionShotgun::Command::TaxonomyReport{
    is => 'Genome::Command::Base',
    doc => 'Generate metagenomic reports for a MetagenomicCompositionShotgun build.',
    has => [
        build => {
            is => 'Genome::Model::Build::MetagenomicCompositionShotgun',
            doc => 'MCS build',
        },
    ],
};

sub execute {
    my $self = shift;

    my $build = $self->build;
    $self->status_message('Metagenomic report for '.$build->__display_name__);

    if ( not -d $self->build->reports_directory ) {
        $self->error_message("Build report directory does not exist");
        return;
    }
    $self->status_message("Report dir: " . $self->build->reports_directory);

    my $refcov_output = $build->refcov_output;
    if ( not -s $refcov_output ) {
        $self->error_message("Refcov output ($refcov_output) does not exist");
        return;
    }

    my $taxonomy_file = $self->build->metagenomic_reference_taxonomy_file;
    if ( not -s $taxonomy_file ) {
        $self->error_message("Taxonomy file ($taxonomy_file) does not exist");
        return;
    }

    my $taxonomy_report = $self->_generate_taxonomy_report;
    return if not $taxonomy_report;

    my $refcov_summary = $self->_generate_refcov_summary;
    return if not $refcov_summary;

    return 1;
}

sub _generate_taxonomy_report {
    my $self = shift;

    $self->status_message("Generate taxonomy report");

    $self->status_message('Load reference hits');
    my $ref_counts = $self->_load_reference_hits;
    return if not $ref_counts;
    $self->status_message('Load reference hits...OK');

    $self->status_message('Load taxonomy');
    my $taxonomy = $self->_load_basic_taxonomy;
    return if not $taxonomy;
    $self->status_message('Load taxonomy...OK');

    my $viral_taxonomy = {};
    if ( -s $self->build->metagenomic_reference_viral_taxonomy_file ) {
        $self->status_message('Load viral taxonomy');
        $viral_taxonomy = $self->_load_viral_taxonomy;
        return if not $viral_taxonomy;
        $self->status_message('Load viral taxonomy...OK');
    }

    my %species_counts_hash;
    my %phyla_counts_hash;
    my %genus_counts_hash;
    my %viral_family_counts_hash;
    my %viral_subfamily_counts_hash;
    my %viral_genus_counts_hash;
    my %viral_species_counts_hash;

    $self->status_message('Open read count file: '.$self->build->read_count_file);
    unlink $self->build->read_count_file if -e $self->build->read_count_file;
    my $read_cnt_o = eval{ Genome::Sys->open_file_for_writing($self->build->read_count_file); };
    if ( not $read_cnt_o ) {
        $self->error_message('Failed to open read count file: '.$@);
        return;
    }

    print $read_cnt_o "Reference Name\t#Reads with hits\tSpecies\tPhyla\tHMP genome\n";
    $self->status_message('Going through references');
    for my $ref_id (sort keys %$ref_counts){
        if (($ref_id =~ /^BACT/) or ($ref_id =~ /^ARCH/) or ($ref_id =~ /^EUKY/)){
            my $species= $taxonomy->{$ref_id}->{species} || '';
            $species_counts_hash{$species}+=$ref_counts->{$ref_id};
            my $phyla=$taxonomy->{$ref_id}->{phyla} || '';
            $phyla_counts_hash{$phyla}+=$ref_counts->{$ref_id};
            my $genus=$taxonomy->{$ref_id}->{genus} || '';
            $genus_counts_hash{$genus}+=$ref_counts->{$ref_id};
            my $hmp_flag=$taxonomy->{$ref_id}->{hmp}|| '';	
            my $order=$taxonomy->{$ref_id}->{order} || '';
            print $read_cnt_o "$ref_id\t$ref_counts->{$ref_id}\t$species\t$phyla\t$genus\t$order\t$hmp_flag\n";
        }elsif ($ref_id =~ /^VIRL/){ #produce reports for viral taxonomy if available
            my ($gi) = $ref_id =~/^VIRL_(\d+)$/;
            if ($viral_taxonomy->{$gi}){
                my $species = $viral_taxonomy->{$gi}->{species} || '';
                $viral_species_counts_hash{$species}+=$ref_counts->{$ref_id};
                my $genus = $viral_taxonomy->{$gi}->{genus} || '';
                $viral_genus_counts_hash{$genus}+=$ref_counts->{$ref_id};
                my $family = $viral_taxonomy->{$gi}->{family} || '';
                $viral_family_counts_hash{$family}+=$ref_counts->{$ref_id};
                my $subfamily = $viral_taxonomy->{$gi}->{subfamily} || '';
                $viral_subfamily_counts_hash{$subfamily}+=$ref_counts->{$ref_id};
                print $read_cnt_o "$ref_id\t$ref_counts->{$ref_id}\t$species\t\t\t\t\n";
            }else{
                print $read_cnt_o "$ref_id\t$ref_counts->{$ref_id}\t\t\t\t\t\n";
            }
        }else{
            print $read_cnt_o "$ref_id\t$ref_counts->{$ref_id}\t\t\t\t\t\n";
        }
    }
    $read_cnt_o->close;
    $self->status_message('Going through references...OK');

    $self->status_message("Write taxonomy counts");
    my $species_output_file = $self->build->reports_directory . '/species_count';
    $self->_write_count_and_close($species_output_file, "Species", \%species_counts_hash);
    my $phyla_output_file = $self->build->reports_directory . '/phyla_count';
    $self->_write_count_and_close($phyla_output_file, "Phyla", \%phyla_counts_hash);
    my $genus_output_file = $self->build->reports_directory . '/genus_count';
    $self->_write_count_and_close($genus_output_file, "Genus", \%genus_counts_hash);
    my $viral_species_output_file = $self->build->reports_directory . '/viral_species_count';
    $self->_write_count_and_close($viral_species_output_file, "Viral Species", \%viral_species_counts_hash);
    my $viral_genus_output_file = $self->build->reports_directory . '/viral_genus_count';
    $self->_write_count_and_close($viral_genus_output_file, "Viral Genus", \%viral_genus_counts_hash);
    my $viral_family_output_file = $self->build->reports_directory . '/viral_family_count';
    $self->_write_count_and_close($viral_family_output_file, "Viral Family", \%viral_family_counts_hash);
    my $viral_subfamily_output_file = $self->build->reports_directory . '/viral_subfamily_count';
    $self->_write_count_and_close($viral_subfamily_output_file, "Viral Subfamily", \%viral_subfamily_counts_hash);
    $self->status_message("Write taxonomy counts...OK");

    $self->status_message("Generate taxonomy report...OK");

    return 1;
}

sub _load_taxonomy {
    my ($self, $filename, $header_ignore_str, $taxon_map_ref) = @_;

    my $fh = Genome::Sys->open_file_for_reading($filename);
    my $taxonomy = {};
    my $header = <$fh>;
    unless ($header =~ /^$header_ignore_str/) {
        die "unexpected header $header!  expected =~ $header_ignore_str";
    }
    while (<$fh>) {
        chomp;
        if (/^$header_ignore_str/) {
            die "duplicated header?!?!: $_\n";
        }
        my @fields = split(/\t/, $_);
        for (@fields) {
            s/^\s+//; 
            s/\s+$//;
        }
        # todo: this is a one-line hash slice -ss
        my $ref_id = $fields[0]; 
        for my $taxon (keys %$taxon_map_ref) {
            $taxonomy->{$ref_id}{$taxon} = $fields[$taxon_map_ref->{$taxon}];
        }
    }
    return $taxonomy;
}

sub _load_basic_taxonomy {
    my $self = shift;

    my $taxonomy = $self->_load_taxonomy($self->build->metagenomic_reference_taxonomy_file, 'Species', {
            species => '1',
            phyla   => '2',
            genus   => '3',
            order   => '4',
            hmp     => '5',
        }
    );

    unless(%$taxonomy) {
        $self->error_message("No taxonomy data loaded from " . $self->build->metagenomic_reference_taxonomy_file . "!");
        return;
    }

    return $taxonomy;
}

sub _load_viral_taxonomy {
    my $self = shift;

    my $viral_taxonomy = $self->_load_taxonomy(
        $self->build->metagenomic_reference_viral_taxonomy_file, 
        'gi', 
        {
            species    => '1',
            genus      => '2',
            subfamily  => '3',
            family     => '4',
            infraorder => '5',
            suborder   => '6',
            superorder => '7',
        },
    );

    unless(%$viral_taxonomy) {
        $self->error_message("No viral taxonomy data loaded from " . $self->build->metagenomic_reference_viral_taxonomy_file . "!");
        return;
    }

    return $viral_taxonomy;
}

sub _load_reference_hits {
    my $self = shift;

    # Count Reference Hits
    my %ref_counts_hash;
    my $ignore_unmapped;
    my $sorted_bam = $self->build->_final_metagenomic_bam;
    my $fh = IO::File->new("samtools view $sorted_bam |");
    while (<$fh>){
        my @fields = split(/\t/, $_);
        my $bitflag = $fields[1];
        if ($bitflag & 0x0004){
            $ignore_unmapped++;
            next;
        }

        my ($ref_name, $null, $gi) = split(/\|/, $fields[2]);
        if ($ref_name eq "VIRL"){
            $ref_name .= "_$gi";
        }
        $ref_counts_hash{$ref_name}++;
    }

    $self->status_message("Skipped $ignore_unmapped reads without a metagenomic mapping");

    return \%ref_counts_hash;
}

sub _generate_refcov_summary {
    my $self = shift;

    $self->status_message('Generate ref cov summary');

    my $read_counts_fh = eval{ Genome::Sys->open_file_for_reading($self->build->read_count_file); };
    if ( not $read_counts_fh ) {
        $self->error_message('Failed to open read count file ('.$self->build->read_count_file.'): '.$@);
        return;
    }
    my $ref_data;
    while (my $line = $read_counts_fh->getline) {
        chomp $line;
        next if ($line =~ /^Reference/);
        my @array=split(/\t/,$line);
        my $ref = $array[0];
        $ref = 'VIRL' if $ref =~/VIRL/;
        if ($ref eq 'VIRL'){
            $ref_data->{$ref}->{reads}+=$array[1];
        }else{
            $ref_data->{$ref}->{reads}=$array[1];
            $ref_data->{$ref}->{species}=$array[2];
            $ref_data->{$ref}->{phyla}=$array[3];
            $ref_data->{$ref}->{genus}=$array[4];
            $ref_data->{$ref}->{order}=$array[5];
            $ref_data->{$ref}->{hmp}=$array[6];
        }
    }
    $read_counts_fh->close;

    my $taxonomy_fh = eval{ Genome::Sys->open_file_for_reading($self->build->metagenomic_reference_taxonomy_file); };
    if ( not $taxonomy_fh ) {
        $self->error_message('Failed to open taxonomy file ('.$self->build->metagenomic_reference_taxonomy_file.'): '.$@);
        return;
    }
    my %header_hash;
    while (my $line = $taxonomy_fh->getline) {
        chomp $line;
        my ($ref, $species) = split(/\t/,$line);
        my ($gi) = split(/\|/, $ref);
        ($gi) = $gi =~ /([^>]+)/;
        $header_hash{$gi}=$species;
    }
    $taxonomy_fh->close;

    if ( -s $self->build->metagenomic_reference_viral_headers_file ) {
        my $viral_taxonomy_fh = Genome::Sys->open_file_for_reading($self->build->metagenomic_reference_viral_headers_file);
        while (my $line = $viral_taxonomy_fh->getline) {
            chomp $line;
            my ($gi, @species) = split(/\s+/,$line);
            my $species = "@species";
            $gi = "VIRL_$gi";
            $header_hash{$gi}=$species;
        }
        $viral_taxonomy_fh->close;
    }

    my $refcov_fh = eval{ Genome::Sys->open_file_for_reading($self->build->refcov_output); };
    if ( not $refcov_fh ) {
        $self->error_message('Failed to open refcov output ('.$self->build->refcov_output.'): '.$@);
        return;
    }
    my $data;
    while(my $line = $refcov_fh->getline) {
        chomp $line;
        my (@array)=split(/\t/,$line);
        my ($ref)  =split(/\|/, $array[0]);

        my $species = $header_hash{$ref};

        #Assuming that average coverage is calculated over the whole reference instead of just the covered reference. 
        my $cov=$array[2]*$array[5];#2 is total ref bases 5 is avg coverage

        #Refcov fields
        $data->{$ref}->{cov}+=$cov;
        $data->{$ref}->{tot_bp}+=$array[2];	    	
        $data->{$ref}->{cov_bp}+=$array[3];
        $data->{$ref}->{missing_bp}+=$array[4];
    }
    $refcov_fh->close;

    my $summary_report = $self->build->reports_directory."/metagenomic_refcov_summary.txt";
    my $summary_report_fh = eval{ Genome::Sys->open_file_for_writing($summary_report); };
    if ( not $summary_report_fh ) {
        $self->error_message("Failed to open summary report ($summary_report): $@");
        return;
    }
    print $summary_report_fh "Reference Name\tPhyla\tGenus\tOrder\tHMP flag\tDepth\tBreadth\tTotal reference bases\tBases not covered\t#Reads\n";
#foreach my $s (keys%{$data}){
    for my $s (sort {$a cmp $b} keys%{$data}){
        my $desc=$header_hash{$s};
        $desc ||= $s;
        next if $desc =~/^gi$/;
        my ($phy, $hmp, $reads, $gen, $ord);
        if ( $ref_data->{$s}->{reads}){
            $phy=$ref_data->{$s}->{phyla};
            $gen=$ref_data->{$s}->{genus};
            $ord=$ref_data->{$s}->{order};
            $hmp=$ref_data->{$s}->{hmp};
            $reads=$ref_data->{$s}->{reads};
        }
        $phy ||= '-';
        $hmp ||= 'N';
        $reads ||= 0;
        $ord ||= '-';
        $gen ||= '-';

        my $new_avg_cov=$data->{$s}->{cov}/$data->{$s}->{tot_bp};
        my $new_avg_breadth=$data->{$s}->{cov_bp}*100/$data->{$s}->{tot_bp};
        my $total_bp = $data->{$s}->{tot_bp};
        my $missing_bp = $data->{$s}->{missing_bp};
        print $summary_report_fh "$desc\t$phy\t$gen\t$ord\t$hmp\t$new_avg_cov\t$new_avg_breadth\t$total_bp\t$missing_bp\t$reads\n";
    }

    $self->status_message("metagenomic report successfully completed");

    system("touch ".$self->build->reports_directory."/FINISHED");

    return 1;
}

sub _write_count_and_close {
    my($self, $filename, $title, $counts_ref) = @_;
    unlink $filename if -e $filename;
    my $file_o=Genome::Sys->open_file_for_writing($filename);
    print $file_o "$title Name\t#Reads with hits\n";
    for my $name (keys %$counts_ref){
        next if (($name eq "") or ($name =~ /^\s+$/));
        print $file_o "$name\t$counts_ref->{$name}\n";
    }
    $file_o->close;
}

1;


