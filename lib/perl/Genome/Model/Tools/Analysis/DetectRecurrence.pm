package Genome::Model::Tools::Analysis::DetectRecurrence;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use FileHandle;

class Genome::Model::Tools::Analysis::DetectRecurrence {
    is => 'Command',
    has => [
    major_group => { 
        type => 'String',
        is_optional => 1,
        doc => "name of the major model groups to process, separated by comma",
    },
    major_files => {
        type => 'String',
        is_optional => 1,
        doc => "files' name (in annotation format), including path, could contain '*'",
    },
    major_files_regex => {
        type => 'String',
        is_optional=> 1,
        doc => "project name of major files, eg. \"AML\\d+\"",
    },
    second_group => {
        type => 'String',
        is_optional => 1,
        doc => "name of the second group of models, separated by comma",
    },
    second_files => {
        type => 'String',
        is_optional => 1,
        doc => "files' name (in annotation format), including path, could contain '*'",
    },    
    second_files_regex => {
        type => 'String',
        is_optional => 1,
        doc => "project name of second files, eg. \"AML\\d+\"",
    },
    include_silent => {
        type => 'Boolean',
        is_optional => 1,
        doc => "whether or not to include silent mutations in the recurrence",
        default => 0,
    },
    summary_output_file => {
        type => 'String',
        is_optional => 1,
        doc => "summary output file name, default is summary.recurrence.csv in your working directory",
    },
    detailed_output_file => {
        type => 'String',
        is_optional => 1,
        doc => "detailed output file (containing annotation information) name, default is detailed.recurrence.csv in your working directory",
    },
    merge_major_groups => {
        type => 'Boolean', 
        is_optional => 1,
        doc => "whether or not to merge the major_groups into one group when reporting",
        default => 0,
    }, 
    
    ]
};

#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-group /gscmnt/sata423/info/medseq/analysis/AML*/curated*/varScan.output.snp.targeted.tier1.snvs.pass_strand.post_filter.curated.anno --second-group /gscmnt/sata423/info/medseq/analysis/AML*/curated_validation_results/varScan.output.snp.targeted.pass_filter.pvalue_filtered.post_reviewed.curated.anno.sorted.tier1 --third-group TCGA-LAML-Exomes-Somatic

#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-group TCGA-LAML-Exomes-Somatic --second-group TCGA-AML-Exome-Imported-Somatic,TCGA-AML-WashU-SomaticCapture
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-group TCGA-AML-WashU-SomaticCapture,TCGA-LAML-Exomes-Somatic --second-group TCGA-AML-Exome-Imported-Somatic,TCGA-AML-WashU-SomaticCapture
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-files /gscmnt/sata423/info/medseq/analysis/AML*/curated_validation_results/varScan.output.snp.targeted.pass_filter.pvalue_filtered.post_reviewed.curated.anno.sorted.tier1 --major-files-project TCGA-AML --second-files /gscmnt/sata423/info/medseq/analysis/AML*/curated*/varScan.output.snp.targeted.tier1.snvs.pass_strand.post_filter.curated.anno --second-files-project AML-M1M3

#perl -d:ptkdb -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated_validation_results/varScan.output.snp.targeted.pass_filter.pvalue_filtered.post_reviewed.curated.anno.sorted.tier1" --second-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated*/varScan.output.snp.targeted.tier1.snvs.pass_strand.post_filter.curated.anno"

########################
#test cases
########################

#major & second files
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated_validation_results/varScan.output.snp.targeted.pass_filter.pvalue_filtered.post_reviewed.curated.anno.sorted.tier1" --second-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated*/varScan.output.snp.targeted.tier1.snvs.pass_strand.post_filter.curated.anno"
#passed

#major file only
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated_validation_results/varScan.output.snp.targeted.pass_filter.pvalue_filtered.post_reviewed.curated.anno.sorted.tier1" 
#passed

#one major groups only
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-group TCGA-AML-WashU-SomaticCapture
#passed

#major groups
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-group TCGA-AML-WashU-SomaticCapture,TCGA-AML-Exome-Imported-Somatic

#major groups & one second group
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-group TCGA-AML-WashU-SomaticCapture,TCGA-AML-Exome-Imported-Somatic --second-group TCGA-LAML-Exomes-Somatic-all
#passed

#major groups & second groups
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-group TCGA-AML-WashU-SomaticCapture,TCGA-AML-Exome-Imported-Somatic --second-group TCGA-LAML-Exomes-Somatic-all,Ley-AML-50-somatic

#major groups, second groups, major files and second files
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-group TCGA-AML-WashU-SomaticCapture,Ley-AML-50-somatic --second-group TCGA-LAML-Exomes-Somatic-all,TCGA-AML-Exome-Imported-Somatic --major-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated_validation_results/varScan.output.snp.targeted.pass_filter.pvalue_filtered.post_reviewed.curated.anno.sorted.tier1" --second-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated*/varScan.output.snp.targeted.tier1.snvs.pass_strand.post_filter.curated.anno"
#passed

#major file only
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated_validation_results/varScan.output.snp.targeted.pass_filter.pvalue_filtered.post_reviewed.curated.anno.sorted.tier1" --major-files-regex "ank\d+"

#major/second files/groups
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-group TCGA-AML-WashU-SomaticCapture,Ley-AML-50-somatic --second-group TCGA-LAML-Exomes-Somatic-all,TCGA-AML-Exome-Imported-Somatic --major-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated_validation_results/varScan.output.snp.targeted.pass_filter.pvalue_filtered.post_reviewed.curated.anno.sorted.tier1" --major-files-regex "AML\d+" --second-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated*/varScan.output.snp.targeted.tier1.snvs.pass_strand.post_filter.curated.anno" --second-files-regex "AML\d+" 

#testing the merger-major-groups function
#perl -I /gscuser/llin/git/Genome/lib/perl/ `which gmt` analysis detect-recurrence --major-group TCGA-AML-WashU-SomaticCapture,Ley-AML-50-somatic --second-group TCGA-LAML-Exomes-Somatic-all,TCGA-AML-Exome-Imported-Somatic --major-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated_validation_results/varScan.output.snp.targeted.pass_filter.pvalue_filtered.post_reviewed.curated.anno.sorted.tier1" --major-files-regex "AML\d+" --second-files "/gscmnt/sata423/info/medseq/analysis/AML*/curated*/varScan.output.snp.targeted.tier1.snvs.pass_strand.post_filter.curated.anno" --second-files-regex "AML\d+"  --merge-major-groups
#################
#model groups

#my $model_group="TCGA-AML-WashU-SomaticCapture"; #somatic capture
#my $model_group="Ley-AML-50-somatic"; #somatic
#my $model_group="TCGA-LAML-Exomes-Somatic"; #somatic capture this one includes TCGA-AML-Exome-Imported-Somatic and TCGA-AML-WashU-SomaticCapture (small group) 
#my $model_group="TCGA-AML-Exome-Imported-Somatic"; #somatic capture
#TCGA-LAML-Exomes-Somatic-all (4288) the complete group


sub execute 
{
        my $self=shift;
    
        my @arr_majormodel_hash=();
        my @arr_secmodel_hash=();
        my %majorfile_detect_hash;
        my %secfile_detect_hash;
        
        my %merged_major_group_hash;
        my @arr_majormodel_hash_merged;
        if ($self->major_group)
        {
            my @major_groups=split/\,/, $self->major_group;
            foreach my $major_group (@major_groups)
            {    
                my %majorgroup_var_hash=%{$self->build_var_hash($major_group)};
                my %majorgroup_detect_hash=%{$self->build_detect_hash(\%majorgroup_var_hash)};
                push @arr_majormodel_hash, \%majorgroup_detect_hash;        
                %merged_major_group_hash=(%merged_major_group_hash, %majorgroup_detect_hash);
            }
            push @arr_majormodel_hash_merged, \%merged_major_group_hash;
        }
    
        if($self->second_group) 
        {
                my @sec_groups=split/\,/, $self->second_group;
                foreach my $sec_group (@sec_groups)
                {       
                        my %secondgroup_var_hash =%{$self->build_var_hash($sec_group)};
                        my %secondgroup_detect_hash=%{$self->build_detect_hash(\%secondgroup_var_hash)};
                        push @arr_secmodel_hash, \%secondgroup_detect_hash;
                }
        }
 
         
        if($self->major_files)
        {
                unless ($self->major_files_regex)
                {
                        print "Please input a regex for major files to detect sample name, eg. \"AML\\d+\"\n";
                        return 1;
                }
                
                my %majorfile_var_hash=%{$self->build_var_hash_from_files($self->major_files, $self->major_files_regex)};
                %majorfile_detect_hash=%{$self->build_detect_hash(\%majorfile_var_hash)};
        }
        
        if ($self->second_files)
        {
                unless ($self->second_files_regex)
                {
                        print "Please input a regex for second files to detect sample name, eg. \"AML\\d+\"\n";
                        return 1;
                }    
                my %secfile_var_hash=%{$self->build_var_hash_from_files($self->second_files, $self->second_files_regex)};
                %secfile_detect_hash=%{$self->build_detect_hash(\%secfile_var_hash)};                  
        }
           
        

        if ($self->major_group || $self->major_files)
        {
                if ($self->merge_major_groups)
                {
                        $self->report_recur(\@arr_majormodel_hash_merged, \@arr_secmodel_hash, \%majorfile_detect_hash, \%secfile_detect_hash);
                }
                else
                {
               
                        $self->report_recur(\@arr_majormodel_hash, \@arr_secmodel_hash, \%majorfile_detect_hash, \%secfile_detect_hash);
                }
        }
        else
        {
                print "Please input at least a major_group or major_files\n";
        }
        return 1;
}

1;

sub help_brief {
    "detect recurrence among model groups/files"
}

sub help_detail {
    <<'HELP';
Hopefully this script will detect recurrence within first model group, and detect overlap between the first and second model group/files
HELP
}

sub build_var_hash_from_files
{
    my ($self, $file_name, $regex)=@_;
    $regex=~s/\//\/\//;  #\ => \\
    my $compiled_regex = qr{($regex)};
    my @files=glob($file_name);
    my %variants;
    foreach my $file (@files)
    {
        my $fh=new FileHandle($file);
        #my ($sample_name)=$file=~/\/(AML\d+)\//;        #need help here
        my $sample_name;
        eval { ($sample_name)=$file=~$compiled_regex;
        print "An error occurred: $@" if $@; 
        };
        while(<$fh>)
        {
            chomp;
            my ($chr,$start,$stop,$ref,$var,$mut_type,$gene,@others)=split/\t/, $_;
            my ($mutation_class) = $others[6];

            next if ($gene eq "-");
            next if (!$self->include_silent && $mutation_class eq 'silent');
            if ($sample_name)
            {
                $variants{$gene}{project_name} .= "$sample_name;";
            
                $variants{$gene}{lines} .= "$_;";
            }
            else
            {
                print "cannot detect sample name, please check your input regex, existing...\n";
                last;
            }
        }
    }
    return \%variants;
}

sub report_recur
{
   
    my ($self,$majorgroup_refarr, $secondgroup_refarr, $majorfile_refhash, $secondfile_refhash)=@_;
    my ($major_group, @major_groups, $second_group, @second_groups, %genelist);
        
    my %majorfile_hash=%{$majorfile_refhash};
    my %secondfile_hash=%{$secondfile_refhash};
    
    ######### 
    #open files for output
    if ($self->summary_output_file)
    {
        my $sum_file=$self->summary_output_file;
        open SUM, ">$sum_file";
    }
    else
    {
        open SUM, ">summary.recurrence.csv";
    }
   
    if ($self->detailed_output_file)
    {
        my $detail_file=$self->detailed_output_file;
        open DETAIL, ">$detail_file";
    }
    else
    {
        open DETAIL, ">detailed.recurrence.csv";
    }
    #################
    
    print SUM "Gene\tTotal\t";
    
    if ($self->major_group) 
    {
        $major_group = $self->major_group;
       
        @major_groups=split/,/, $major_group;
        my $major_group_str=$major_group;
        $major_group_str=~s/\,/\t/;
        
        if ($self->merge_major_groups)
        {
                $major_group_str=$major_group;
                @major_groups=($major_group);
        }
        
        $self->build_genelist_hash_from_groups(\@major_groups, $majorgroup_refarr, \%genelist);
        print SUM $major_group_str; #print the header for summary file
        
    }
        
    if ($self->second_group) 
    {
        $second_group=$self->second_group;
        @second_groups=split/,/, $second_group;
        my $second_group_str=$second_group;
        $second_group_str=~s/\,/\t/;
        print SUM "\t";   #print the header for summary file
        print SUM $second_group_str;
        $self->build_genelist_hash_from_groups(\@second_groups, $secondgroup_refarr, \%genelist);
    }
    
    if ($self->major_files)  
    {
        my $major_files=$self->major_files;
        print SUM "\t$major_files"; #print the header for summary file
        
        $self->build_genelist_hash_from_files($majorfile_refhash, $self->major_files, \%genelist);
        push @major_groups, $self->major_files;  #now @major_groups includes both groups and files
    
    }
    
    if ($self->second_files)  
    {
        my $second_files=$self->second_files;
        print SUM "\t$second_files";  #print the header for summary file
        
        $self->build_genelist_hash_from_files($secondfile_refhash, $self->second_files, \%genelist);
        push @second_groups, $self->second_files; #now @second_groups includes both groups and files
    }
    print SUM "\n";
    
    #need to make a gene list that contains genes found in all major/second groups
    #genelist{gene}{group}=num_samples

    my @major_arr=@{$majorgroup_refarr};
    my @second_arr=@{$secondgroup_refarr};

    my %print_hash; #for printing purpose totalnumber_gene=>\@sample_numbs
    foreach my $gene (sort keys %genelist)
    {
        my @groups= keys %{$genelist{$gene}};
        my $group_number=@groups;
        my $flag=0;
        
        if ($group_number<2) #need to be in major group and is recurrent in order to be reported
        {
                foreach my $group (@groups)
                {
                        #if (grep /^$group$/, @major_groups)  #not working
                        foreach my $major_g (@major_groups)
                        {
                                if ($major_g eq $group) #if in major group
                                {
                                        my $sample_num=$genelist{$gene}{$group};
                                        if ($sample_num>1) 
                                        {
                                                #report
                                                $flag=1;
                                        
                                                last;
                                        }
                                }
                        }
                     
                }
        }
        else #more than one groups have this gene, can report
        {
                $flag=1;
        }
        
        next unless ($flag==1); #report 
        
        #print SUM "$gene\t";
        my @whole_groups=(@major_groups,@second_groups);
        my $total_number=0;
        my @sample_numbs;
        foreach my $group (@whole_groups)
        {
                my $sample_num=0;
                
                if (defined $genelist{$gene}{$group})
                {
                        $sample_num=$genelist{$gene}{$group};
                }
                push @sample_numbs, $sample_num;
               $total_number+=$sample_num;
        }
        
        $print_hash{$total_number."_".$gene}=\@sample_numbs;
        #print SUM "$total_number\t";
        #print SUM join "\t", @sample_numbs;
        #print SUM "\n";
        
        #print detailed annotation line, need to include everyhash here
        foreach my $major_refhash (@major_arr)
        {
                $self->print_details($major_refhash, $gene);           
        }
        
        foreach my $second_refhash (@second_arr)
        {
                $self->print_details($second_refhash, $gene);
        }
        
        $self->print_details($majorfile_refhash, $gene);
        $self->print_details($secondfile_refhash, $gene);

    }
    
    ##################
    #print sorted summary file
    foreach my $k (reverse sort numer keys %print_hash)
    {
        my ($total_number, $gene)=split/_/, $k;
        my @sample_numbs=@{$print_hash{$k}};
        print SUM "$gene\t$total_number\t";
        print SUM join "\t", @sample_numbs;
        print SUM "\n";
    }

    
}

sub build_genelist_hash_from_groups
{
        my ($self, $arr_ref, $group_refhash_arr, $refhash)=@_;
        my @group_name_arr=@{$arr_ref};
        my @group_arr=@{$group_refhash_arr};
        #my %genelist=%{$refhash};
        for (my $i=0; $i<=$#group_name_arr;$i++)
        {
                my $refhash_i=$group_arr[$i];
                my %hash_i=%{$refhash_i};
                my @genes=sort keys %{$refhash_i};
                foreach my $gene (@genes)
                {
                        my @samples=keys %{$hash_i{$gene}};
                        my $sample_num=@samples;
                        my $group_name=$group_name_arr[$i];
                        $refhash->{$gene}->{$group_name}=$sample_num;
                }
        }        
        
}


sub print_details
{
        my ($self, $hashref, $gene)=@_;
        my %Hash=%{$hashref};
        foreach my $sample (keys %{$Hash{$gene}})
        {
                my @lines=@{$Hash{$gene}{$sample}};
                foreach my $line (@lines)
                {
                        print DETAIL "$gene\t$sample\t$line\n";
                }
        }
}

sub numer 
{

        my ($count_a) = $a =~ /^(\d+)/;
        my ($count_b) = $b =~ /^(\d+)/;

        return $count_a <=> $count_b;
}

sub build_genelist_hash_from_files
{
        my ($self,$refhash,$project,$genelist_refhash)=@_;
        my %file_hash=%{$refhash};
        
        foreach my $gene (sort keys %file_hash)
        {
                my @samples=sort keys %{$file_hash{$gene}};
                my $sample_num=@samples;
                $genelist_refhash->{$gene}->{$project}=$sample_num;
        }
}
    

sub build_var_hash
{
    my ($self,$model_group)=@_;
    my %pos_hash;
    my %variants;
    my ($indel_anno, $snp_anno, $t1_hc_snp, $t1_snp, $t1_indel);
    my @builds = Genome::ModelGroup->get(name => $model_group)->builds; 
    for my $build (@builds) 
    {
        my $model=$build->model;
        my $type=$model->type_name;
        my $full_tcga_id_name=$model->tumor_model->subject->extraction_label;
        my ($project_name) = $full_tcga_id_name =~ /(TCGA-\S{2}-\d{4})/;
        unless ($project_name)
        {
                $project_name=$model->subject->common_name;
        }
        
        my $data_directory = $build->data_directory; 
        #print "$data_directory\t$project_name\n";
        $indel_anno="$data_directory/ani_annotated_indel.csv" if ($type eq "somatic");
        $snp_anno="$data_directory/anv_annotated_snp.csv" if ($type eq "somatic");
        $t1_hc_snp="$data_directory/hc1_tier1_snp_high_confidence.csv" if ($type eq "somatic");
        $t1_snp="$data_directory/t1v_tier1_snp.csv" if ($type eq "somatic");
        $t1_indel="$data_directory/t1i_tier1_indel.csv" if ($type eq "somatic");

        $snp_anno="$data_directory/annotation.somatic.snp.transcript" if ($type eq "somatic capture");
        $t1_hc_snp="$data_directory/merged.somatic.snp.filter.novel.tier1.hc" if ($type eq "somatic capture");

        my $t1_hc_snp_fh = IO::File->new($t1_hc_snp,"r") or die "Can't open $t1_hc_snp\n";
        while(<$t1_hc_snp_fh>)
        {
            chomp;
            my ($chr,$start,$stop, @others) = split /\t/;
            my $pos = $chr."_".$start."_".$stop;
            $pos_hash{$pos}=1;
        }
        $t1_hc_snp_fh->close;

        my $snp_anno_fh=IO::File->new($snp_anno,"r") or die "Can't open $snp_anno\n";

        while(<$snp_anno_fh>)
        {
            chomp;
            my ($chr,$start,$stop,$ref,$var,$mut_type,$gene,@others)=split/\t/, $_;
            next if ($gene eq "-");

            my $pos=$chr."_".$start."_".$stop;
            next unless (defined ($pos_hash{$pos}));
            $variants{$gene}{project_name} .= "$project_name;";
            $variants{$gene}{lines} .= "$_;";
        }
        $snp_anno_fh->close;
    }   
    return \%variants;
} 

sub build_detect_hash
{
    my ($self,$var_refhash)=@_;
    my %variants=%{$var_refhash};
    my %detect_gene;  #detect_gene{gene}{sample}=line
    foreach my $gene (sort keys %variants) 
    {
        my @samples_in = split /;/, $variants{$gene}{project_name};
        my @lines_in = split /;/, $variants{$gene}{lines};
        #my %samples = map {$_ => 1} @samples_in;

        #print "$gene:\n";

        foreach my $sample (@samples_in) 
        {
            my $line = shift @lines_in;

            if (defined $detect_gene{$gene}{$sample})
            {
                unless (grep /^$line$/, @{$detect_gene{$gene}{$sample}})
                {
                    push @{$detect_gene{$gene}{$sample}}, $line;
                }
            }
            else
            {
                my @arr=($line);
                $detect_gene{$gene}{$sample}=\@arr;
            }

            #print "\t$line\t$sample\n";
        }
    }
    return \%detect_gene;
}

