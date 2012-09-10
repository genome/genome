package Genome::Model::Tools::Cufflinks::Cummerbund;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Cufflinks::Cummerbund {
    is => 'Genome::Model::Tools::Cufflinks',
    has_input => [
        cuffdiff_directory => { },
        output_directory => {
            is_optional => 1,
        },
        genes_of_interest => {
            is_optional => 1,
        },
        image_format => {
            is_optional => 1,
            default_value => 'png',
            valid_values => ['png','jpeg'],
        },
        image_width => {
            is_optional => 1,
            default_value => 800,
        },
        image_height => {
            is_optional => 1,
            default_value => 600,
        },
        image_background => {
            is_optional => 1,
            default_value => 'transparent',
        },
    ],
};

sub execute {
    my $self = shift;

    unless ($self->output_directory) {
        $self->output_directory($self->cuffdiff_directory);
    }

    my $rscript_file = $self->output_directory.'/cummerbund.R';

    my $tmp_fh = Genome::Sys->open_file_for_writing($rscript_file);
    unless ($tmp_fh) {
        die('Failed to open Rscript file for writing: '. $rscript_file);
    }

    my $fpkm_tracking_file = $self->cuffdiff_directory .'/genes.fpkm_tracking';
    my $fpkm_tracking_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $fpkm_tracking_file,
        separator => "\t",
    );

    my @fpkm_tracking_headers = @{$fpkm_tracking_reader->headers};
    my @fpkm_headers = grep {$_ =~ /_FPKM$/} @fpkm_tracking_headers;
    my %labels;
    for my $fpkm_header (@fpkm_headers) {
        if ($fpkm_header =~ /^(.*)_FPKM/) {
            my $label = $1;
            $label =~ s/-/_/g;
            if ($labels{$label}) {
                die('Multiple labels found for '. $label);
            }
            $labels{$label} = 1;
        } else {
            die('Failed to match FPKM header : '. $fpkm_header);
        }
    }
    my @labels = sort keys %labels;
    my @comparisons;
    for (my $i = 0; $i < scalar(@labels); $i++) {
        for (my $j = ($i+1); $j < scalar(@labels); $j++) {
            push @comparisons, [$labels[$i], $labels[$j]];
        }
    }
    #print Data::Dumper::Dumper(@comparisons);

    # Load Library and Data
    print $tmp_fh 'library("cummeRbund")' ."\n";
    print $tmp_fh 'cuff <- readCufflinks(dir="'. $self->cuffdiff_directory .'")' ."\n";

    # Create Gene Set
    if ($self->genes_of_interest) {
        my @goi = split(',',$self->genes_of_interest);
        my $goi_string = '"'. join('","',@goi) .'"';
        print $tmp_fh 'myGeneIds <- c('. $goi_string .')' ."\n";
        print $tmp_fh 'myGenes <- getGenes(cuff, myGeneIds)' ."\n";
        print $tmp_fh 'myGenes' ."\n";
    }

    # Over-dispersion
    print $tmp_fh $self->resolve_image_device_string('overdispersion') ."\n";
    print $tmp_fh 'disp<-dispersionPlot(genes(cuff))' ."\n";
    print $tmp_fh 'disp' ."\n";
    print $tmp_fh 'dev.off()' ."\n";
    
    # Gene FPKM Density Histogram
    print $tmp_fh $self->resolve_image_device_string('gene_fpkm_density') ."\n";
    print $tmp_fh 'dens <- csDensity(genes(cuff))'. "\n";
    print $tmp_fh 'dens' ."\n";
    print $tmp_fh 'dev.off()' ."\n";

    # Gene FPKM Density Histogram as Replicates
    print $tmp_fh $self->resolve_image_device_string('gene_fpkm_density_as_replicates') ."\n";
    print $tmp_fh 'densRep<-csDensity(genes(cuff),replicates=T)'. "\n";
    print $tmp_fh 'densRep' ."\n";
    print $tmp_fh 'dev.off()' ."\n";
    
    # Gene FPKM Boxplot
    print $tmp_fh $self->resolve_image_device_string('gene_fpkm_boxplot') ."\n";
    print $tmp_fh 'b <- csBoxplot(genes(cuff))' . "\n";
    print $tmp_fh 'b' ."\n";
    print $tmp_fh 'dev.off()' ."\n";

    # Gene FPKM Boxplot as Replicates
    print $tmp_fh $self->resolve_image_device_string('gene_fpkm_boxplot_as_replicates') ."\n";
    print $tmp_fh 'brep<-csBoxplot(genes(cuff),replicates=T)' . "\n";
    print $tmp_fh 'brep' ."\n";
    print $tmp_fh 'dev.off()' ."\n";

    # Gene FPKM Dendrogam
    print $tmp_fh $self->resolve_image_device_string('gene_fpkm_dendrogram') ."\n";
    print $tmp_fh 'den <- csDendro(genes(cuff))' ."\n";
    print $tmp_fh 'den' ."\n";
    print $tmp_fh 'dev.off()' ."\n";

    # Gene FPKM Dendrogam as Replicates
    print $tmp_fh $self->resolve_image_device_string('gene_fpkm_dendrogram_as_replicates') ."\n";
    print $tmp_fh 'dend.rep<-csDendro(genes(cuff),replicates=T)' ."\n";
    print $tmp_fh 'dend.rep' ."\n";
    print $tmp_fh 'dev.off()' ."\n";

    if ($self->genes_of_interest) {
        # GOI Gene FPKM Heatmap
        print $tmp_fh $self->resolve_image_device_string('gene_fpkm_goi_heatmap') ."\n";
        print $tmp_fh 'h.goi <- csHeatmap(myGenes, cluster = "both")' ."\n";
        print $tmp_fh 'h.goi'. "\n";
        print $tmp_fh 'dev.off()' ."\n";

        # GOI Gene FPKM Heatmap as Replicates
        print $tmp_fh $self->resolve_image_device_string('gene_fpkm_goi_heatmap_as_replicates') ."\n";
        print $tmp_fh 'h.rep<-csHeatmap(myGenes,cluster=\'both\',replicates=T)'."\n";
        print $tmp_fh 'h.rep'. "\n";
        print $tmp_fh 'dev.off()'. "\n";

        # GOI Gene FPKM Barplot
        print $tmp_fh $self->resolve_image_device_string('gene_fpkm_goi_barplot') ."\n";
        print $tmp_fh 'b.goi <- expressionBarplot(myGenes)' ."\n";
        print $tmp_fh 'b.goi' ."\n";
        print $tmp_fh 'dev.off()' ."\n";

        # GOI Gene FPKM Dendrogram
        print $tmp_fh $self->resolve_image_device_string('gene_fpkm_goi_dendrogram') ."\n";
        print $tmp_fh 'den.goi<-csDendro(myGenes)' ."\n";
        print $tmp_fh 'den.goi' ."\n";
        print $tmp_fh 'dev.off()'. "\n";

        # Generate per-gene plots
        my @goi = split(',',$self->genes_of_interest);
        for my $goi (@goi) {
            print $tmp_fh 'myGeneId <- "'. $goi .'"' ."\n";
            print $tmp_fh 'myGene <- getGene(cuff,myGeneId)'. "\n";
            
            # GOI FPKM Expression Line Plot
            print $tmp_fh $self->resolve_image_device_string($goi .'_fpkm_expression_lineplot') ."\n";
            print $tmp_fh 'gl <- expressionPlot(myGene)' ."\n";
            print $tmp_fh 'gl'. "\n";
            print $tmp_fh 'dev.off()'. "\n";

            # GOI FPKM Expression Line Plot as Replicates
            print $tmp_fh $self->resolve_image_device_string($goi .'_fpkm_expression_lineplot_as_replicates') ."\n";
            print $tmp_fh 'gl.rep<-expressionPlot(myGene,replicates=TRUE)' ."\n";
            print $tmp_fh 'gl.rep'. "\n";
            print $tmp_fh 'dev.off()'. "\n";

            # GOI FPKM Expression Barplot
            print $tmp_fh $self->resolve_image_device_string($goi .'_fpkm_expression_barplot') ."\n";
            print $tmp_fh 'gb<-expressionBarplot(myGene)' ."\n";
            print $tmp_fh 'gb'. "\n";
            print $tmp_fh 'dev.off()'. "\n";

            # GOI FPKM Expression Barplot as replicates
            print $tmp_fh $self->resolve_image_device_string($goi .'_fpkm_expression_barplot_as_replicates') ."\n";
            print $tmp_fh 'gb.rep<-expressionBarplot(myGene,replicates=T)' ."\n";
            print $tmp_fh 'gb.rep'. "\n";
            print $tmp_fh 'dev.off()'. "\n";
        }
    }
    
    # Individual Comparisons between experimental conditions
    for my $comparison (@comparisons) {
        my $condition_1 = $comparison->[0];
        my $condition_2 = $comparison->[1];

        # Gene FPKM Scatter Plots
        my $scatter_basename = 'gene_fpkm_scatter_'. $condition_1 .'-vs-'. $condition_2;
        print $tmp_fh $self->resolve_image_device_string($scatter_basename) ."\n";
        print $tmp_fh 's <- csScatter(genes(cuff), "'. $condition_1 .'", "'. $condition_2 .'", smooth = T)' ."\n";
        print $tmp_fh 's' ."\n";
        print $tmp_fh 'dev.off()' ."\n";

        # MAPlot
        print $tmp_fh $self->resolve_image_device_string('maplot_'. $condition_1 .'-vs-'. $condition_2) ."\n";
        print $tmp_fh 'm<-MAplot(genes(cuff),"'. $condition_1 .'","'. $condition_2  .'")' ."\n";
        print $tmp_fh 'm' ."\n";
        print $tmp_fh 'dev.off()'. "\n";

        # MAPlot Count based
        print $tmp_fh $self->resolve_image_device_string('maplot_count_'. $condition_1 .'-vs-'. $condition_2) ."\n";
        print $tmp_fh 'mCount<-MAplot(genes(cuff),"'. $condition_1 .'","'. $condition_2 .'",useCount=T)'. "\n";
        print $tmp_fh 'mCount' ."\n";
        print $tmp_fh 'dev.off()'. "\n";
        
        # Gene FPKM Volcano Plots
        my $volcano_basename = 'gene_fpkm_volcano_'. $condition_1 .'-vs-'. $condition_2;
        print $tmp_fh $self->resolve_image_device_string($volcano_basename) ."\n";
        print $tmp_fh 'v <- csVolcano(genes(cuff), "'. $condition_1 .'", "'. $condition_2 .'")' ."\n";
        print $tmp_fh 'v' ."\n";
        print $tmp_fh 'dev.off()' ."\n";
        if ($self->genes_of_interest) {
            # GOI Gene FPKM Scatter Plots
            my $goi_scatter_basename = 'gene_fpkm_goi_scatter_'. $condition_1 .'-vs-'. $condition_2;
            print $tmp_fh $self->resolve_image_device_string($goi_scatter_basename) ."\n";
            print $tmp_fh 's <- csScatter(myGenes, "'. $condition_1 .'", "'. $condition_2 .'", smooth = T)' ."\n";
            print $tmp_fh 's' ."\n";
            print $tmp_fh 'dev.off()' ."\n";

            # GOI Gene FPKM Volcano Plots
            my $goi_volcano_basename = 'gene_fpkm_goi_volcano_'. $condition_1 .'-vs-'. $condition_2;
            print $tmp_fh $self->resolve_image_device_string($goi_volcano_basename) ."\n";
            print $tmp_fh 'v <- csVolcano(myGenes, "'. $condition_1 .'", "'. $condition_2 .'")' ."\n";
            print $tmp_fh 'v' ."\n";
            print $tmp_fh 'dev.off()' ."\n";
        }
    }

    # TODO: For each GOI add a set of gene-level plots
    # TODO: we really need a way to validate the genes_of_interest list with actual gene IDs in cuffdiff
    
    $tmp_fh->close;
    my $cmd = 'Rscript '. $rscript_file;
    eval {
        Genome::Sys->shellcmd(
            cmd => $cmd,
        );
    };
    if ($@) {
        warn('The Rscript on '. $rscript_file .' did not run correctly: '. $@);
        return 0;
    }
    return 1;
}

sub resolve_image_device_string {
    my $self = shift;
    my $basename = shift;
    my $file = $self->resolve_image_file($basename);
    return $self->open_image_device_string($file);
}

sub resolve_image_file {
    my $self = shift;
    my $basename = shift;

    return $self->output_directory .'/'. $basename .'.'. $self->image_format;
}

sub open_image_device_string {
    my $self = shift;
    my $output_file = shift;
    
    my $string = $self->image_format .'("'. $output_file .'", bg="'. $self->image_background .'", width='. $self->image_width .', height='. $self->image_height .')';
    return $string;
}


1;
