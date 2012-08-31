package Genome::Model::Tools::Analysis::MutationSpectrum;

#####################################################################################################################################
# MutationSpectrum - Given an annotation file, gives an output of transition/transversion, cpg islands, and cpgs within cpg islands.
#					
#	AUTHOR:		Will Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	3/05/2010 by W.S.
#	MODIFIED:	3/05/2010 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use warnings;
use strict;

use Genome;
use Workflow;
use Carp;
use FileHandle;
use Data::Dumper;
use List::Util qw( max );
use IO::File;
use Genome::Info::IUB;
use DBI;
use Cwd qw( abs_path );
class Genome::Model::Tools::Analysis::MutationSpectrum {
    is => ['Command'],
    has => [
    fasta_file => { 
        is  => 'String',
        is_input=>1, 
        doc => 'The input fasta file.',
    },
    mutation_file => { 
        is  => 'String',
        is_input=>1, 
        doc => 'The input annotated mutation file.',
    },
    chr_col => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        default_value => '1',
        doc => 'The number (1-based) of the column that points to Chromosome',
    },
    start_col => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        default_value => '2',
        doc => 'The number (1-based) of the column that points to Start Position',
    },
    stop_col => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        default_value => '3',
        doc => 'The number (1-based) of the column that points to Stop Position',
    },
    ref_col => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        default_value => '4',
        doc => 'The number (1-based) of the column that points to Reference Base',
    },
    var_col => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        default_value => '5',
        doc => 'The number (1-based) of the column that points to Tumor Base #1',
    },
    type_col => { 
        is  => 'String',
        is_input=>1, 
        is_optional => 1,
        default_value => '14',
        doc => 'The number (1-based) of the column that points to syn/silent type column',
    },
    output_trans_file => {
        is  => 'String',
        is_input => '1',
        is_output => '1',
        doc => 'The Transition/Transversion output file.',
    },
    output_CpGFinder_file => {
        is  => 'String',
        is_input => '1',
        is_output => '1',
        doc => 'The CpG Island Finder output file.',
    },
    output_CpGSites_file => {
        is  => 'String',
        is_input => '1',
        is_output => '1',
        doc => 'The CpG Within CpG output file.',
    },
    plot_spectrum_file_name => {
        is => 'String',
        is_optional => 1,
        is_input => 1,
        is_output => 1,
        doc => 'Filename of the pdf to output of the spectrum bar graph',
    },
    plot_spectrum_genome_name => {
        is => 'String',
        is_optional => 1,
        is_input => 1,
        is_output => 1,
        default => "",
        doc => 'name to put into the title of the spectrum plot',
    },   
    absolute_axis => {
        is => 'Boolean',
        is_optional => 1,
        is_input => 1,
        is_output => 1,
        default => 1,
        doc => 'Whether or not to force the plot y-axis to be 0-100',
    },   
    #below here are variables with which to store results
    #hopefully these can then be judiciously used to write cross-comparison scripts

    _total_lines => {
        is => 'Integer',
        is_optional => 1,
        doc => "Total number of lines in the file",
    },   
    _num_indels => {
        is => 'Integer',
        is_optional => 1,
        doc => "Total number of indels in the file",
    },   
    _num_mnps => {
        is => 'Integer',
        is_optional => 1,
        doc => "Total number of mnp (likely dnps) in the file",
    },   
    _num_snvs => {
        is => 'Integer',
        is_optional => 1,
        doc => "Total number of snvs in the file",
    },   
    _transitions_transversion_count => {
        is => 'HashRef',
        is_optional => 1,
        doc => "Hash containing the counts for the various base changes as hash{base1}{base2}",
    },   
    _num_transitions => {
        is => 'Integer',
        is_optional => 1,
        doc => "Number of transitions",
    },
    _num_transversions => {
        is => 'Integer',
        is_optional => 1,
        doc => "Number of transversions",
    },

    _synonymous_transitions_transversion_count => {
        is => 'HashRef',
        is_optional => 1,
        doc => "Hash containing the counts for the various base changes for synonymous changes as hash{base1}{base2}",
    },   
    #not dealing with CpG/NpG stuff yet

    # Make workflow choose 64 bit blades
    lsf_resource => {
        is_param => 1,
        default_value => 'rusage[mem=4000] select[type==LINUX64] span[hosts=1]',
    },
    lsf_queue => {
        is_param => 1,
        default_value => 'long'
    }, 
    ],
};

sub help_brief {
    "Given an annotation file, gives an output of transition/transversion, cpg islands, and cpgs within cpg islands.";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis mutation-spectrum --fasta-file --mutation-file --chr-col --start-col --stop-col --ref-col --var-col --type-col --output-trans-file --output-CpGFinder-file --output-CpGSites-file

EXAMPLE:
gmt analysis mutation-spectrum --fasta-file /path/to/your/reference_sequences/NCBI-human-build36/all_sequences.fa --mutation-file /gscuser/charris/testing/TSP_Lung_Adenocarcinoma_Mutation_Table2-29aug2008i.tsv --chr-col 3 --start-col 4 --stop-col 5 --ref-col 12 --var-col 14 --type-col 6 --output-trans-file exampleout.transitiontransversion.csv --output-CpGFinder-file exampleout.CpGFinder.csv --output-CpGSites-file exampleout.CpGIsland.csv
EOS
}

sub help_detail {                           
    return <<EOS 
    Given an annotation file, gives an output of transition/transversion, cpg islands, and cpgs within cpg islands.
EOS
}

sub execute {
    my $self = shift;

    my %base_complement;
    $base_complement{'A'} = 'T';
    $base_complement{'T'} = 'A';
    $base_complement{'G'} = 'C';
    $base_complement{'C'} = 'G';

    my $indels=0;
    my $DNP=0;
    my $lines = 0;
    #first command fasta, next is variant file
    my $fasta = $self->fasta_file;
    unless(-s $self->fasta_file) {
        $self->error_message(<'Could not open fasta input file '$fasta' for reading'>);
        return;
    }
    my $fh = IO::File->new($self->mutation_file);
    unless($fh) {
        my $errorfile = $self->mutation_file;
        $self->error_message(<'Could not open mutation input file '$errorfile' for reading'>);
        return;
    }

    #next columns are pointer to new columns
    my $chr_col = ($self->chr_col - 1);
    my $start_col = ($self->start_col - 1);
    my $stop_col = ($self->stop_col - 1);
    my $ref_col = ($self->ref_col - 1);
    my $var_col = ($self->var_col - 1);
    #next column points to syn/silent type column
    my $type_col = ($self->type_col - 1);

    my $delimiter="\t";
    my %nonsyn_type = ( A => { C => 0, G => 0, T => 0 } , C => { A => 0, G => 0, T => 0 } );
    my %synonomous_type = ( A => { C => 0, G => 0, T => 0 } , C => { A => 0, G => 0, T => 0 } );
    my $db = "ucsc";
    my $user = "mgg_admin";
    my $password = "c\@nc3r"; 
    my $dataBase = "DBI:mysql:$db:mysql2";
    my $dbh = DBI->connect($dataBase, $user, $password) || die "ERROR: Could not connect to database: $! \n";
    my %cpg_count;
    my $cpgs_in_islands;
    my $cpgs_not_in_islands;
    my %npg_count;
    my %cpg_island_count;
    my %cpg_finder_count;
    my %npg_island_count;
    my %npg_finder_count;

    $cpgs_in_islands=2044537;
    $cpgs_not_in_islands=25594255;
    print "cpgs in islands: $cpgs_in_islands\n";
    my $query = "SELECT name from cpgIslandExt where chrom=? and chromStart <= ? and chromEnd >= ?";
    my $cpg_statement=$dbh->prepare($query);

    my $basename |= 'exampleout';
    my $output_trans = $self->output_trans_file || $basename."transitiontransversion.csv";
    my $output_island = $self->output_CpGSites_file || $basename."CpGIsland.csv";
    my $output_finder = $self->output_CpGFinder_file || $basename."CpGFinder.csv";
    # Open Output Transition/Transversion
    unless (open(TRANS,">$output_trans")) {
        die "Could not open output file '$output_trans' for writing";
    }
    # Open Output CpG Island
    unless (open(CPG_ISLAND,">$output_island")) {
        die "Could not open output file '$output_island' for writing";
    }
    # Open Output CpG Finder
    unless (open(CPG_FINDER,">$output_finder")) {
        die "Could not open output file '$output_finder' for writing";
    }
    my $i = 0;
    while (my $line=$fh->getline) {
        $lines++;
        chomp $line;
        my @fields = split $delimiter, $line;
        my $chr = $fields[$chr_col];
        unless ($chr =~ m/^\d+/ || $chr =~ m/^X+/ || $chr =~ m/^Y+/ || $chr =~ m/^MT+/ ){
            print "Header: $line\n";
            next;
        }
        my $ref_base = $fields[$ref_col];
        my $var_base = $fields[$var_col];
        my $ref = $ref_base;
        my $var = $var_base;
        my $start = $fields[$start_col] - 1;
        my $stop = $fields[$stop_col] + 1;
        my $type = $fields[$type_col];

        unless($ref_base =~ m/^[ACGT]$/) {
            ($ref_base, $var_base) = ($ref_base =~ m/([ACGT])>([ACGT])/);
        }
        if($type =~ m/del/i || $type =~ m/ins/i) {
            $indels++;
            next;
        }
        if($ref =~ m/\D{2,}/i || $var =~ m/\D{2,}/i) {
            print "Skipped (likely DNP): $line\n";
            $DNP++;
            next;
        }
        if (defined $ref_base && defined $var_base) {
            if($ref_base eq $var_base) {
                $var_base = $fields[$var_col+1];
            }
            if($ref_base eq '-' || $var_base eq '-') {
                $indels++;
                next;
            }
        }
        else {
            print "Unidentified variant type in line $lines\n";
        }

        my @variant_alleles = Genome::Info::IUB->variant_alleles_for_iub($ref_base, $var_base);
        unless (@variant_alleles) {
            print "Error reading variant alleles for: $line\n";
        }

        for my $var_base ( @variant_alleles) {
            if($ref_base eq 'G' || $ref_base eq 'T') {
                $ref_base = $base_complement{$ref_base};
                $var_base = $base_complement{$var_base};
            }
            if($type) {
                if($fields[$type_col] =~ m/(Synonymous|Silent)/i) {
                    $synonomous_type{$ref_base}{$var_base}+=1;          
                }
            }
            $nonsyn_type{$ref_base}{$var_base}+=1;   
        }

        my @var_alleles = Genome::Info::IUB::variant_alleles_for_iub($ref, $var);
        #   changed commands for better functionality
        #   my $cmd = "xdget -n $fasta.db $chr -a $start -b $stop";
        my $cmd = "samtools faidx $fasta $chr:$start-$stop";
        my @return = `$cmd`;
        my @var_alleles2 = @var_alleles;
        my @return2 = @return;
        next unless ($ref eq 'G' || $ref eq 'C');
        for my $variant(@var_alleles) {
            for my $line (@return) {
                $cpg_statement->execute("chr$chr",$fields[$start_col], $fields[$stop_col]);
                if($cpg_statement->fetchrow_array) {
                    #cpg island cpg
                    $cpg_island_count{$ref}{$variant}+=1 while $line =~ /cg/icg;
                }
                else {
                    #non cpg island cpg
                    $npg_island_count{$ref}{$variant}+=1 while $line =~ /cg/icg;
                }
            }
        }
        for my $variant(@var_alleles2) {
            for my $line (@return2) {
                #cpg island
                $cpg_finder_count{$ref}{$variant}+=1 while $line =~ /cg/icg;
                #non cpg island
                $npg_finder_count{$ref}{$variant}+=1 while ($line =~ /(AG|TG|GG|CT|CA|CC)/icg); 
            }
        }
    }
    $fh->close;

    #Headers
    print TRANS "base->base_change\tAll\tSynonomous\n";
    print CPG_ISLAND "base->base_change\tcpg_island_changes\tnpg_island_changes\t(cpg_island_changes+npg_island_changes)\tcpg_island_changes/(cpg_island_changes+npg_island_changes)\n";
    print CPG_FINDER "base->base_change\tcpg_finder_changes\tnpg_finder_changes\t(cpg_finder_changes+npg_finder_changes)\tcpg_finder_changes/(cpg_finder_changes+npg_finder_changes)\n";

    my $transitions = 0;
    my $transversions = 0;
    for my $base (sort keys %nonsyn_type) {
        for my $base_change (sort keys %{$nonsyn_type{$base}} ) { 
            unless (exists ($nonsyn_type{$base}{$base_change})) {
                $nonsyn_type{$base}{$base_change} = 0;
            }
            unless (exists ($synonomous_type{$base}{$base_change})) {
                $synonomous_type{$base}{$base_change} = 0;
            }
            #tally as transition or transversion
            #only two allowed transitions so test explicitly
            #Nonsynonymous appears to actually be ALL changes
            if(($base eq 'A' && $base_change eq 'G') || ($base eq 'C' && $base_change eq 'T')) {
                $transitions += $nonsyn_type{$base}{$base_change};
            }
            else {
                $transversions += $nonsyn_type{$base}{$base_change};
            }
            print TRANS "$base->$base_change\t" . $nonsyn_type{$base}{$base_change} . "\t" . $synonomous_type{$base}{$base_change} .  "\n";
        }
    } 
    print TRANS "Transitions\t$transitions\nTransversions\t$transversions\n";
    my $total_snvs = $transversions + $transitions;
    print TRANS "SNVs\t$total_snvs\nIndels\t$indels\nDNP$DNP\nTotal\t$lines";
    close TRANS;

    #set class values
    $self->_num_transversions($transversions);
    $self->_num_transitions($transitions);
    $self->_synonymous_transitions_transversion_count(\%synonomous_type);
    $self->_transitions_transversion_count(\%nonsyn_type);
    $self->_num_snvs($total_snvs);
    $self->_num_mnps($DNP);
    $self->_num_indels($indels);
    $self->_total_lines($lines);

    for my $base (sort keys %cpg_island_count) {
        for my $base_change (sort keys %{$cpg_island_count{$base}} ) { 
            my $cpg_island_changes = $cpg_island_count{$base}{$base_change} ? $cpg_island_count{$base}{$base_change} : 0;
            my $npg_island_changes = $npg_island_count{$base}{$base_change} ? $npg_island_count{$base}{$base_change} : 0;
            $cpg_island_changes += $cpg_island_count{$base_complement{$base}}{$base_complement{$base_change}} if exists($cpg_island_count{$base_complement{$base}}{$base_complement{$base_change}});
            $npg_island_changes += $npg_island_count{$base_complement{$base}}{$base_complement{$base_change}} if exists($npg_island_count{$base_complement{$base}}{$base_complement{$base_change}});
            print CPG_ISLAND "$base->$base_change\t$cpg_island_changes\t$npg_island_changes\t" . ($cpg_island_changes+$npg_island_changes) . "\t" .  $cpg_island_changes/($cpg_island_changes + $npg_island_changes) . "\n";
        }
    } 

    for my $base (sort keys %cpg_finder_count) {
        for my $base_change (sort keys %{$cpg_finder_count{$base}} ) { 
            my $cpg_finder_changes = $cpg_finder_count{$base}{$base_change};
            my $npg_finder_changes = $npg_finder_count{$base}{$base_change};
            $cpg_finder_changes += $cpg_finder_count{$base_complement{$base}}{$base_complement{$base_change}} if exists($cpg_finder_count{$base_complement{$base}}{$base_complement{$base_change}});
            $npg_finder_changes += $npg_finder_count{$base_complement{$base}}{$base_complement{$base_change}} if exists($npg_finder_count{$base_complement{$base}}{$base_complement{$base_change}});
            print CPG_FINDER "$base->$base_change\t$cpg_finder_changes\t$npg_finder_changes\t" . ($cpg_finder_changes+$npg_finder_changes) . "\t" .  $cpg_finder_changes/($cpg_finder_changes + $npg_finder_changes) . "\n";
        }
    } 
    my $plot_spectrum_file_name = $self->plot_spectrum_file_name;
    if($plot_spectrum_file_name) {
        #assume that the spectrum plot is requested
        unless($plot_spectrum_file_name =~ /\.pdf$/) {
            $plot_spectrum_file_name .= ".pdf";
        }
        #statistics::R is horrible so try to prevent people from causing themselves problems
        my $abs_filename = abs_path($plot_spectrum_file_name);
        my $abs_output_trans = abs_path($output_trans);
        my $genome = $self->plot_spectrum_genome_name;
        my $absolute_axis_boolean;
        if($self->absolute_axis) {
            $absolute_axis_boolean = "T";
        }
        else {
            $absolute_axis_boolean = "F";
        }

        my $plot_cmd = qq{ plot_spectrum("$abs_output_trans",output_file="$abs_filename",genome="$genome",absolute_axis=$absolute_axis_boolean) };
        my $call = Genome::Model::Tools::R::CallR->create(command=>$plot_cmd, library=> "MutationSpectrum.R");
        $call->execute;
    }


    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

