package Genome::Model::Tools::Cufflinks::Cummerbund;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Cufflinks::Cummerbund {
    is => 'Genome::Model::Tools::Cufflinks',
    has_input => [
        cuffdiff_directory => { },
        false_discovery_rate => {
            default_value => '0.05',
            is_optional => 1,
        },
        output_directory => {
            is_optional => 1,
        },
        genes_of_interest => {
            is_optional => 1,
        },
        pdf_report_file => {
            is => 'Text',
            doc => 'Redirect graphic output to the defined path.',
            is_optional => 1,
        },
        image_format => {
            is_optional => 1,
            doc => 'The output format for each image',
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
    print $tmp_fh 'library("IRanges")' . "\n";
    print $tmp_fh 'library("ggplot2")' . "\n";
    print $tmp_fh 'library("cummeRbund")' ."\n";
    if ($self->pdf_report_file) {
        print $tmp_fh 'pdf("'. $self->pdf_report_file .'")' ."\n";
    }
    print $tmp_fh 'cuff <- readCufflinks(dir="'. $self->cuffdiff_directory .'")' ."\n";

    # Create Gene Set
    # TODO: we really need a way to validate the genes_of_interest list with actual gene IDs in cuffdiff
    if ($self->genes_of_interest) {
        my @goi = split(',',$self->genes_of_interest);
        my $goi_string = '"'. join('","',@goi) .'"';
        print $tmp_fh 'myGeneIds <- c('. $goi_string .')' ."\n";
        print $tmp_fh 'myGenes <- getGenes(cuff, myGeneIds)' ."\n";
        print $tmp_fh 'myGenes' ."\n";
    }

    # Over-dispersion
    $self->generate_r_plot_code($tmp_fh,'overdispersion','disp','dispersionPlot(genes(cuff))');
    
    # Gene FPKM Density Histogram
    $self->generate_r_plot_code($tmp_fh,'gene_fpkm_density','dens','csDensity(genes(cuff))');

    # Gene FPKM Density Histogram as Replicates
    $self->generate_r_plot_code($tmp_fh,'gene_fpkm_density_as_replicates','densRep','csDensity(genes(cuff),replicates=T)');
    
    # Gene FPKM Boxplot
    $self->generate_r_plot_code($tmp_fh,'gene_fpkm_boxplot','b','csBoxplot(genes(cuff))');

    # Gene FPKM Boxplot as Replicates
    $self->generate_r_plot_code($tmp_fh,'gene_fpkm_boxplot_as_replicates','brep','csBoxplot(genes(cuff),replicates=T)');

    # Gene FPKM Dendrogam
    $self->generate_r_plot_code($tmp_fh,'gene_fpkm_dendrogram','den','csDendro(genes(cuff))');

    # Gene FPKM Dendrogam as Replicates
    $self->generate_r_plot_code($tmp_fh,'gene_fpkm_dendrogram_as_replicates','den.rep','csDendro(genes(cuff),replicates=T)');

    # Global DE Genes
    print $tmp_fh 'mySigGeneIds <- getSig(cuff,alpha='. $self->false_discovery_rate .',level=\'genes\')' ."\n";
    print $tmp_fh 'mySigGenes <- getGenes(cuff, mySigGeneIds)' ."\n";
    print $tmp_fh 'mySigGenes' ."\n";

    $self->generate_r_plot_code($tmp_fh,'sig_gene_fpkm_heatmap','h.sig','csHeatmap(mySigGenes,cluster=\'both\',labRow=F)');
    $self->generate_r_plot_code($tmp_fh,'sig_gene_fpkm_heatmap_as_replicates','h.sig.rep','csHeatmap(mySigGenes,cluster=\'both\',replicates=T,labRow=F)');
    
    if ($self->genes_of_interest) {
        # GOI Gene FPKM Heatmap
        $self->generate_r_plot_code($tmp_fh,'gene_fpkm_goi_heatmap','h.goi','csHeatmap(myGenes, cluster = "both")');

        # GOI Gene FPKM Heatmap as Replicates
        $self->generate_r_plot_code($tmp_fh,'gene_fpkm_goi_heatmap_as_replicates','h.goi.rep','csHeatmap(myGenes,cluster=\'both\',replicates=T)');
        
        # GOI Gene FPKM Barplot
        $self->generate_r_plot_code($tmp_fh,'gene_fpkm_goi_barplot','b.goi','expressionBarplot(myGenes)');
        
        # GOI Gene FPKM Dendrogram
        $self->generate_r_plot_code($tmp_fh,'gene_fpkm_goi_dendrogram','den.goi','csDendro(myGenes)');

        # GOI Gene FPKM Dendrogram as Replicates
        $self->generate_r_plot_code($tmp_fh,'gene_fpkm_goi_dendrogram_as_replicates','den.goi.rep','csDendro(myGenes,replicates=T)');

        # Generate per-gene plots
        my @goi = split(',',$self->genes_of_interest);
        for my $goi (@goi) {
            print $tmp_fh 'myGeneId <- "'. $goi .'"' ."\n";
            print $tmp_fh 'myGene <- getGene(cuff,myGeneId)'. "\n";
            
            # GOI FPKM Expression Line Plot
            $self->generate_r_plot_code($tmp_fh,$goi .'_fpkm_expression_lineplot','gl','expressionPlot(myGene)');

            # GOI FPKM Expression Line Plot as Replicates
            $self->generate_r_plot_code($tmp_fh,$goi .'_fpkm_expression_lineplot_as_replicates','gl.rep','expressionPlot(myGene,replicates=TRUE)');

            # GOI FPKM Expression Barplot
            $self->generate_r_plot_code($tmp_fh,$goi .'_fpkm_expression_barplot','gb','expressionBarplot(myGene)');

            # GOI FPKM Expression Barplot as Replicates
            $self->generate_r_plot_code($tmp_fh,$goi .'_fpkm_expression_barplot_as_replicates','gb.rep','expressionBarplot(myGene,replicates=T)');
        }
    }
    
    # Individual Comparisons between experimental conditions
    for my $comparison (@comparisons) {
        my $condition_1 = $comparison->[0];
        my $condition_2 = $comparison->[1];

        # Gene FPKM Scatter Plots
        my $scatter_basename = 'gene_fpkm_scatter_'. $condition_1 .'-vs-'. $condition_2;
        $self->generate_r_plot_code($tmp_fh,$scatter_basename,'s','csScatter(genes(cuff), "'. $condition_1 .'", "'. $condition_2 .'", smooth = T)');

        # MAPlot
        $self->generate_r_plot_code($tmp_fh,'MAplot_'. $condition_1 .'-vs-'. $condition_2,'m','MAplot(genes(cuff),"'. $condition_1 .'","'. $condition_2  .'")');

        # MAPlot Count based
        $self->generate_r_plot_code($tmp_fh,'MAplot_count_'. $condition_1 .'-vs-'. $condition_2,'mCount','MAplot(genes(cuff),"'. $condition_1 .'","'. $condition_2 .'",useCount=T)');
        
        # Gene FPKM Volcano Plots
        my $volcano_basename = 'gene_fpkm_volcano_'. $condition_1 .'-vs-'. $condition_2;
        $self->generate_r_plot_code($tmp_fh,$volcano_basename,'v','csVolcano(genes(cuff), "'. $condition_1 .'", "'. $condition_2 .'")');
        
        if ($self->genes_of_interest) {
            # GOI Gene FPKM Scatter Plots
            my $goi_scatter_basename = 'gene_fpkm_goi_scatter_'. $condition_1 .'-vs-'. $condition_2;
            $self->generate_r_plot_code($tmp_fh,$goi_scatter_basename,'s','csScatter(myGenes, "'. $condition_1 .'", "'. $condition_2 .'", smooth = T)');

            # GOI Gene FPKM Volcano Plots
            my $goi_volcano_basename = 'gene_fpkm_goi_volcano_'. $condition_1 .'-vs-'. $condition_2;
            $self->generate_r_plot_code($tmp_fh,$goi_volcano_basename,'v','csVolcano(myGenes, "'. $condition_1 .'", "'. $condition_2 .'")');
        }
    }
    if ($self->pdf_report_file) {
        print $tmp_fh 'dev.off()' ."\n";
    }
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


sub generate_r_plot_code {
    my $self = shift;
    my $fh = shift;
    my $file_basename = shift;
    my $r_variable = shift;
    my $r_routine = shift;

    unless ($self->pdf_report_file) {
        # Write a command that will open an output image device
        print $fh $self->resolve_image_device_string($file_basename) ."\n";
    }
    # Print the command to generate the image
    print $fh $r_variable .' <- '. $r_routine."\n";
    # Print the command to render the image
    print $fh $r_variable ."\n";
    unless ($self->pdf_report_file) {
        # Close the image output device
        print $fh 'dev.off()' ."\n";
    }
    return 1;
}

1;
