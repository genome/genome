package Genome::Model::Tools::Predictor::Interproscan;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Predictor::Interproscan {
    is => 'Genome::Model::Tools::Predictor::Base',
    doc => 'execute interproscan',
};

sub requires_chunking {
    return 0;
}

sub run_predictor {
    my $self = shift;

    my $unsorted_file = $self->raw_output_path . '.unsorted';
    my $sorted_file = $self->raw_output_path;

    my $tool_path = $self->tool_path_for_version($self->version);
    my $cmd = join(' ', $tool_path, $self->parameters, 
        '-i '. $self->input_fasta_file, '-o ' . $unsorted_file .
        '> ' . $self->debug_output_path, '2>&1');

    my $rv = Genome::Sys->shellcmd(
        cmd => $cmd,
    );
    unless ($rv == 1) {
        die "Could not execute interpro command $cmd";
    }

    Genome::Sys->shellcmd(
        cmd => "sort -k 1,1 -k 4,4 -k 12,12 < $unsorted_file > $sorted_file",
    );

    return 1;
}

sub parse_output {
    my $self = shift;
    
    my $output_fh = IO::File->new($self->raw_output_path, 'r');
    unless ($output_fh) {
        die "Could not get file handle for " . $self->raw_output_path;
    }

    my @features;
    while (my $line = $output_fh->getline) {
        chomp $line;
        
        # Raw format described here:
        # ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/README.html#3
        my (
            $protein_name,
            $checksum,
            $length,
            $analysis_method,
            $database_entry,
            $database_member,
            $start,
            $end,
            $evalue,
            $status,
            $run_date,
            $ipr_number,
            $ipr_description,
            $go_description,
        ) = split /\t/, $line;

        # Selection criteria taken from Prat, refactored to be easier to read (according to whoever 
        # wrote this interpro wrapper originally...)
        if (($ipr_number ne 'NULL') and (($evalue eq 'NA') or ($evalue <= 0.01))) {
            if ($evalue eq 'NA') {
                $evalue = 1e10;
            }

            my $feature = Bio::SeqFeature::Generic->new(
                -display_name => $protein_name,
                -primary      => $analysis_method,
                -source_tag   => 'InterPro',
                -start        => $start,
                -end          => $end,
                -score        => $evalue,
            );
            $feature->add_tag_value('interpro_analysis',    $analysis_method);
            $feature->add_tag_value('interpro_evalue',      $evalue);
            $feature->add_tag_value('interpro_description', $ipr_description);
            
            my $dblink = Bio::Annotation::DBLink->new(
                -database   => 'InterPro',
                -primary_id => $ipr_number,
            );
            $feature->annotation->add_Annotation('dblink', $dblink);

            push @features, $feature;
        }
    }

    $self->unfiltered_bio_seq_features(\@features);
    return 1;
}

sub filter_results {
    my $self = shift;

    # Only interested in keeping the features with the lowest evalue per gene
    my %least_evalue;
    for my $feature ($self->unfiltered_bio_seq_features) {
        my $protein_name = $feature->display_name;
        my $analysis_method = $feature->primary_tag;
        my $evalue = $feature->score;

        my ($dblink) = $feature->annotation->get_Annotations('dblink');
        my $ipr_number = $dblink->primary_id;

        # Yet more logic taken from Prat. The original author of this wrapper notes that the <= may
        # be a bug, but his goal was replication of results, so he didn't take the time to determine for sure.
        if (exists($least_evalue{$protein_name}{$analysis_method}{$ipr_number})) {
            if ($evalue <= $least_evalue{$protein_name}{$analysis_method}{$ipr_number}->score()) {
                $least_evalue{$protein_name}{$analysis_method}{$ipr_number} = $feature;
            }
        }
        else {
            $least_evalue{$protein_name}{$analysis_method}{$ipr_number} = $feature
        }
    }

    # Only interested in the result for each gene with the lowest evalue
    my @features;
    foreach my $gene (keys %least_evalue) {
        foreach my $analysis (keys %{$least_evalue{$gene}}) {
            foreach my $ipr (keys %{$least_evalue{$gene}{$analysis}}) {
                push @features, $least_evalue{$gene}{$analysis}{$ipr};
            }
        }
    }
    
    $self->bio_seq_features(\@features);
    return 1;
}

# Originally stolen from /gsc/scripts/gsc/annotation/interpro2ace, but it was discovered
# that script didn't work correctly. So now this logic is drawn from biosql2ace.
sub create_ace_file {
    my $self = shift;

    my $ace_file_fh = Genome::Sys->open_file_for_writing($self->ace_file_path);
    my %ipr;

    for my $feature ($self->bio_seq_features) {
        my $display_name = $feature->display_name();

        my ($dblink) = grep { $_->database() eq 'InterPro' } $feature->annotation->get_Annotations();
        next unless $dblink;

        my $ipr_number = $dblink->primary_id();
        $ipr{$display_name}{$ipr_number} = 1;
        my ($analysis) = $feature->each_tag_value('interpro_analysis');
        my ($evalue)   = $feature->each_tag_value('interpro_evalue');
        my ($desc)     = $feature->each_tag_value('interpro_description');

        ## Bug for bug replication...
        if ($evalue == 1e10) { $evalue = ''; }

        if($evalue) { $evalue = uc sprintf('%.2e', $evalue); }

        $ace_file_fh->print("Sequence $display_name\n");
        $ace_file_fh->print("Interpro   \"$analysis : $ipr_number $desc : pval $evalue\"\n\n");
    }
    $ace_file_fh->close;

    my $ipr_fh = Genome::Sys->open_file_for_writing($self->ace_file_path . '.ipr');
    for my $gene (sort keys %ipr) {
        my $ipr_string = join " ", (sort { $a cmp $b } keys %{$ipr{$gene}});
        $ipr_fh->print("Sequence $gene\n");
        $ipr_fh->print("IPR_ID \"$ipr_string\"\n\n");
    }
    $ipr_fh->close;

    return 1;
}

sub tool_path_for_version {
    my ($self, $version) = @_;
    unless ($version) {
        die "Not given interproscan version!";
    }
    my $path = "/gsc/scripts/pkg/bio/iprscan/iprscan-$version/bin/iprscan";
    unless (-e $path) {
        die "Nothing found at expected interproscan path $path!";
    }
    unless (-x $path) {
        die "Cannot execute $path!";
    }
    return $path;
}

sub raw_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'interpro.output.raw');
}

sub debug_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'interpro.output.error');
}

sub dump_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'interpro.results.dump');
}

sub ace_file_path {
    my $self = shift;
    return join('/', $self->output_directory, 'interpro.ace');
}

1;

