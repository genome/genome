package Genome::Model::Tools::Rna::GraphCufflinks;

use warnings;
use strict;
use Genome;
use Cwd;
use Statistics::R;
use Genome::ModelGroup;
use File::Basename;
require Genome::Sys;

class Genome::Model::Tools::Rna::GraphCufflinks {
    is => 'Command',
    has => [
    model_group => {
        is => 'Text',
        is_optional => 0,
        is_input => 1,
        doc => 'Model group to graph',
    },
    chrom => {
        is_optional=>0,
        doc=>"chromosome to graph"
    },
    start => {
        is_optional=>0,
        doc=>"start of graphing region",
    },
    stop => {
        is_optional=>0,
        doc=>"stop of graphing region",
    },
    output_dir => {
        is_optional=>0,
    },
    rows_columns => {
        is_optional=>1,
        is=>'Text',
        default_value=>'2,3',
        doc=>'comma separated values eg:"2,3", but not implemented yet',
    },
    ]
};

sub help_brief {
    "Plot Tophat and Cufflinks output";
}
sub help_detail {
    "Plot Tophat and Cufflinks output";
}

#this shoudl eventually be derivable from an annotation build
my $gene_annotation_file='/gscmnt/gc2105/info/medseq/jbrea/rnaseq/annotation/NCBI-human.combined-annotation_54_36p_v2_102210.bed';
my $transcript_annotation_file = '/gscmnt/gc2105/info/medseq/jbrea/rnaseq/annotation/Homo_sapiens.GRCh36.61.gtf.sorted.mod.sorted';
my $repeat_mask_dir = '/gscmnt/sata194/info/sralign/UCSC/data/';

sub execute {
    my $self = shift;

    my $gene_annotation_processed_file = $self->gene_annotation($gene_annotation_file);
    my $transcript_annotation_processed_file = $self->transcript_annotation($transcript_annotation_file);
    my $repeat_mask_file = '/gscmnt/sata194/info/sralign/UCSC/data/chr'.$self->chrom.'_rmsk.txt';
    die $self->error_message("Couldn't find repeat mask file $repeat_mask_file for chrom ".$self->chrom) unless -e $repeat_mask_file;
    my $repeat_mask_processed_file = $self->repeat_mask($repeat_mask_file);

    my $rows;
    my $columns;
    if($self->rows_columns =~ /^(\d+),(\d+)$/)
    {
        $rows = 2;
        $columns = 3;
    }
    else
    {
        die $self->error_message("wrong input for rows and columns");
    }

    my($align_fh, $alignments) = Genome::Sys->create_temp_file("alignments");
    my($junction_fh, $junctions) = Genome::Sys->create_temp_file("junctions");
    my($transcripts_fh, $transcripts)= Genome::Sys->create_temp_file("transcripts");
    my ($subject_name_fh, $subject_name_file) = Genome::Sys->create_temp_file("subjects_names");
    my $temp_dir = Genome::Sys->create_temp_directory();
    my $output_dir= $self->output_dir;

    #$DB::single = 1; #breakpoint
    my $model_group = Genome::ModelGroup->get($self->model_group);
    for my $model ($model_group->models) {
        my $latest_build = $model->last_succeeded_build;
        if($latest_build)
        {
            my $subject_name = $model->subject_name;
            my $alignments_file = $latest_build->accumulated_alignments_directory . "/all_reads_merged.bam";
            my $junction_file = $latest_build->accumulated_alignments_directory . "/junctions.bed";
            my $junction_processed_file = $self->junctions($junction_file);
            my $transcript_file = $latest_build->accumulated_expression_directory . "/transcripts.gtf";
            my $transcript_processed_file = $self->cufflinks($transcript_file);
            $subject_name_fh->print($subject_name ."\n");
            $align_fh->print($alignments_file ."\n");
            $junction_fh->print($junction_processed_file ."\n");
            $transcripts_fh->print($transcript_processed_file ."\n");
        }else{
            $self->debug_message("No succeeded build for model ".$model->name.", skipping.");
        }

    } 
    $subject_name_fh->close;
    $align_fh->close;
    $junction_fh->close;
    $transcripts_fh->close;

    my$chr = $self->chrom;
    my$start = $self->start;
    my$stop = $self->stop;
    my($config_fh, $config_file) = Genome::Sys->create_temp_file("configuration");
    print "$config_file\n";
    $config_fh->print("tophat_alignment_path_files: $alignments\n");
    $config_fh->print("tophat_junctions_path_file: $junctions\n");
    $config_fh->print("cufflinks_transcripts_path_file: $transcripts\n");
    $config_fh->print("locus: $chr $start $stop\n");
    $config_fh->print("num_plot_in_each_fig_row: $rows\n");
    $config_fh->print("num_plot_in_each_fig_col: $columns\n");
    $config_fh->print("plot__bamreadcounts: yes\n");
    $config_fh->print("plot_gene_annotation: yes\n");
    $config_fh->print("plot_transcript_annotation: yes\n");
    $config_fh->print("plot_repeat_mask: yes\n");
    $config_fh->print("plot_cufflinks_junctions: yes\n");
    $config_fh->print("plot_cufflinks_transcript_assembly: yes\n");
    $config_fh->print("temp_directory: $temp_dir\n");
    $config_fh->print("output_directory: $output_dir"."/","\n");
    $config_fh->print("gene_annotation_file: $gene_annotation_processed_file\n");
    $config_fh->print("transcript_annotation_file: $transcript_annotation_processed_file\n");
    $config_fh->print("mask_annotation_file: $repeat_mask_processed_file\n");
    $config_fh->print("sample_name:  $subject_name_file\n");
    $config_fh->close();

    print "$output_dir"."/"."\n";
    my($script_fh, $script_file) = Genome::Sys->create_temp_file("matlab_script.m");

    my $dirname = dirname(__FILE__);
    $script_fh->print("P=path;\n");
    $script_fh->print("path(P, '$dirname')\n");
    $script_fh->print("rnaseq_plot_bams('$config_file');\n");
    $script_fh->close();
    #$DB::single=1;
    print `matlab -nodisplay < $script_file`;
}

sub gene_annotation{
    my $self = shift;
    my ($annot) = @_;
    $self->debug_message("processing $annot");
    my $annoth = IO::File->new($annot);
    my $chr = $self->chrom;
    my $start = $self->start;
    my $stop = $self->stop;

    my ($outfileh, $outfile) = Genome::Sys->create_temp_file('gene_annotation');
    while(my $annotation = $annoth->getline)
    {
        chomp($annotation);
        my ($annot_chr,$annot_start,$annot_stop,$region) = split /\s+/,$annotation;

        if(($annot_chr eq $chr) && ($annot_start >= $start) && ($annot_stop <= $stop))
        {
            my ($gene,$trans,$type) = split /:/, $region; if(($type eq "utr_exon") || ($type eq "cds_exon")) {
                my $exon_type;
                if($type eq "utr_exon")
                {
                    $exon_type = 1;
                }
                else
                {
                    $exon_type = 0;
                }
                if($annot_chr eq "X")
                {
                    $annot_chr = 22;
                }
                elsif($annot_chr eq "Y")
                {
                    $annot_chr = 23;
                }
                print $outfileh $annot_chr,"\t",$annot_start,"\t",$annot_stop,"\t",$exon_type,"\n";
            }

        }
    }    
    $annoth->close;
    $outfileh->close;
    return $outfile
}

sub transcript_annotation{
    my $self = shift;
    my ($annot) = @_;
    $self->debug_message("processing $annot");
    my $annoth = IO::File->new($annot);

    my $locus_chrom = $self->chrom;
    my $locus_start = $self->start;
    my $locus_stop = $self->stop;

    my ($outfileh, $outfile) = Genome::Sys->create_temp_file('transcript_annotation');


    $locus_chrom = 'chr'.$locus_chrom;

    while(my $transcript = $annoth->getline)
    {  
        my ($chr,$source,$exon_type,$start,$stop,$score,$strand,$exon_type_number,
            $gene_id_label,$gene_id,$trans_id_label,$trans_id)= split /\s+/,$transcript;
        if($locus_chrom eq $chr)
        {
            if($start > $locus_start && $stop < $locus_stop)
            {  
                print $outfileh $chr,"\t",$source,"\t",$exon_type,"\t",$start,"\t",$stop,"\t",$strand,"\t",
                $exon_type_number,"\t",$gene_id_label,"\t",$gene_id,"\t",$trans_id_label,"\t",$trans_id,"\n";
            }
        }
    }
    return $outfile;
}

sub repeat_mask{
    my $self = shift;
    my ($mask) = @_;
    $self->debug_message("processing $mask");
    my $maskh = IO::File->new($mask);

    my $start = $self->start;
    my $stop = $self->stop;

    my ($outfileh, $outfile) = Genome::Sys->create_temp_file('repeat_mask');

    while(my $annot = $maskh->getline)
    {
        chomp($annot);
        my ($var1,$var2,$var3,$var4,$var5,$annotline_chr,$annotline_start,$annotline_stop,$var6,$var7,$mask_type) = split /\s+/,$annot;

        if(($annotline_start >= $start) && ($annotline_stop <= $stop))# && ($mask_type == "LINE" || $mask_type == "SINE"))
        {
            if($annotline_chr eq "X")
            {
                $annotline_chr = 22;
            }
            elsif($annotline_chr eq "Y")
            {
                $annotline_chr = 23;
            }
            $annotline_chr =~ s/^[A-Za-z0-9]{3}//;
            print $outfileh $annotline_start,"\t",$annotline_stop,"\n";
        }
    }
    $maskh->close;
    $outfileh->close;
    return $outfile;
}

sub junctions{
    my $self = shift;
    my ($junct) = @_;
    $self->debug_message("processing $junct");
    my $juncth = IO::File->new($junct);

    my $chr = $self->chrom;
    my $start = $self->start;
    my $stop = $self->stop;

    my ($outfileh, $outfile) = Genome::Sys->create_temp_file();
    $juncth->getline;
    while (my $junction = $juncth->getline)
    {
        chomp($junction);
        my ($junct_chrom,$start_block_feature,$end_block_feature,$junc_name,
            $cov,$strand,$start_block_pos,$end_block_pos,$rgb,$num_blocks,
            $block_sizes,$block_starts) = split /\s+/,$junction;
        if($junct_chrom eq $chr)
        {
            my ($block_size1,$block_size2) = split /\,/, $block_sizes;  
            my $junct_start = $start_block_feature + $block_size1; 
            my $junct_end = $end_block_feature - $block_size2;
            if($junct_start >= $start && $junct_end <= $stop) 
            {
                if($junct_chrom eq "X")
                {
                    $junct_chrom = 22;
                }
                elsif($junct_chrom eq "Y")
                {
                    $junct_chrom = 23;
                }
                print $outfileh $junct_chrom,"\t",$junct_start, "\t", $junct_end,"\t",$cov,"\n";
            }
        }
    }
    $juncth->close;
    $outfileh->close;
    return $outfile;
}

sub cufflinks{
    my $self = shift;
    my ($cuff) = @_;
    $self->debug_message("processing $cuff");
    my $cuffh = IO::File->new($cuff);
    
    my $chr = $self->chrom;
    my $start = $self->start;
    my $stop = $self->stop;

    my ($outfileh, $outfile) = Genome::Sys->create_temp_file();
    my $exon_sum_reads;
    my $exon_length;
    my $exon_av_cov;
    my ($temp_feature);
    my $gene_count = 0;
    my $gene_track = "";
    while(my $transcript = $cuffh->getline)
    {  
        (undef,undef,$temp_feature)= split /\s+/,$transcript;
        if($temp_feature eq "exon")
        {
            my ($cuff_chrom,$source_name,$feature,$cuff_start,$cuff_stop,$score, $cuff_strand,$frame,$gene_id_label,
                $gene_id,$trans_id_label,$trans_id,$exon_number_label,$exon,$fpkm_label,$fpkm,$frac_label,
                $frac,$conf_lo_label,$conf_lo,$conf_hi_label,$conf_hi,$cov_label,$cuff_cov)= split /\s+/,$transcript;
            my ($gene_id_pref1,$gene_id_name,$gene_id_suff1) = split /"/,$gene_id;
            my ($gene_id_pref2,$gene_id_name_num) = split /\./,$gene_id_name;
            my ($trans_id_pref1,$trans_id_name,$trans_id_suff1) = split /"/,$trans_id;
            my ($trans_id_pref2,$gene_trans_id_name_num,$trans_id_name_num) = split /\./,$trans_id_name;
            my ($exon_pref,$exon_num,$exon_suff) = split /"/,$exon;
            my ($fpkm_pref,$fpkm_num,$fpkm_suff) = split /"/,$fpkm;
            my ($frac_pref,$frac_num,$frac_suff) = split /"/,$frac;
            my ($conf_lo_pref,$conf_lo_num,$conf_lo_suff) = split /"/,$conf_lo;
            my ($conf_hi_pref,$conf_hi_num,$conf_hi_suff) = split /"/,$conf_hi;
            my ($cov_pref,$cov_num,$cov_suff) = split /"/,$cuff_cov;

            if(($cuff_chrom eq $chr)&&($cuff_start>=$start)&&($cuff_stop<=$stop))
            {
                if($gene_id ne $gene_track)
                {
                    $gene_count++;
                    $gene_track = $gene_id;
                }
                if($cuff_chrom eq "X")
                {
                    $cuff_chrom = 22;
                }
                elsif($cuff_chrom eq "Y")
                {
                    $cuff_chrom = 23;
                }
                #print $cuff_chrom,"\t",$cuff_start,"\t",$cuff_stop,"\t",$fpkm_num,"\n";
                print $outfileh $gene_id_name_num,"\t",$gene_count,"\t",$trans_id_name_num,"\t",$exon_num,"\t";
                print $outfileh $cuff_chrom,"\t",$cuff_start,"\t",$cuff_stop,"\t",$fpkm_num,"\t",$conf_lo_num,"\t";
                print $outfileh $conf_hi_num,"\t",$cov_num,"\n";
            }
        }  
    }
    $cuffh->close;
    $outfileh->close;
    return $outfile;
}



