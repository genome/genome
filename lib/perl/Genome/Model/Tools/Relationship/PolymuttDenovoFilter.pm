package Genome::Model::Tools::Relationship::PolymuttDenovoFilter;

use strict;
use warnings;
use Data::Dumper;
use Genome;
use Genome::Info::IUB;
use POSIX;
our $VERSION = '0.01';
use Cwd;
use File::Basename;
use File::Path;

class Genome::Model::Tools::Relationship::PolymuttDenovoFilter {
    is => 'Command',
    has_optional_input => [
    model_group_id => {
        is_optional=>0,
        doc=>'id of model group for family',
    },
    denovo_vcf=> {
        is_optional=>0,
        doc=>'denovo vcf output for same family',
    },
    output_file=> {
        is_optional=>0,
        doc=>'binomial test outputs for each individual in the family',
    },
    min_read_qual=> {
        is_optional=>1,
        doc=>'the lowest quality reads to use in the calculation. default q20',
        default=>"20",
    },
    min_unaffected_pvalue=> {
        is_optional=>1,
        doc=>"the minimum binomial test result from unaffected members to pass through the filter",
        default=>"1.0e-4",
    },
    ],


};

sub help_brief {
    "blahblahblah"
}

sub help_detail {
}

sub execute {
    $DB::single=1;
    my $self=shift;
    my $mg_id = $self->model_group_id;
    my $vcf = $self->denovo_vcf;

    my $sites_file = Genome::Sys->create_temp_file_path();
    my $cat_cmd = "cat";
    my $vcf_fh = open_vcf_file($vcf);
    if(Genome::Sys->file_is_gzipped($vcf)) {
        $cat_cmd = "zcat";
    }

    my $sites_cmd = qq/ $cat_cmd $vcf | cut -f1,2 | grep -v "^#" | awk '{OFS="\t"}{print \$1, \$2, \$2}' > $sites_file/;
    my $header_cmd = qq/ $cat_cmd $vcf | grep "^#CHR"/;
    unless(-s $sites_file) {
        print STDERR "running $sites_cmd\n";
        `$sites_cmd`;
    }
    my $mg = Genome::ModelGroup->get($mg_id);
    my $header_line= `$header_cmd`;
    $DB::single=1;
    my @models = $self->sort_models_by_header($header_line, $mg->models);

###prepare readcounts
    my $ref_fasta=$models[0]->reference_sequence_build->full_consensus_path("fa");
    my $readcount_file_ref = $self->prepare_readcount_files($ref_fasta, $sites_file, \@models, $header_line);
    my $r_input_path = $self->prepare_r_input($vcf_fh, $readcount_file_ref);
    my ($r_fh, $r_script_path) = Genome::Sys->create_temp_file();
    my ($r_output) = Genome::Sys->create_temp_file_path();
    $r_fh->print($self->r_code($r_input_path, $r_output));
    $r_fh->close;
    `R --vanilla < $r_script_path`;
    $self->output_passing_vcf($self->denovo_vcf, $r_output, $self->min_unaffected_pvalue, $self->output_file);
    return 1;
}

sub sort_models_by_header {
    my ($self, $header_line, @models)= @_;
    chomp($header_line);
    my @fields = split "\t", $header_line;
    splice(@fields, 0, 9);
    my @return_models;
    for my $subject_id (@fields) {
        for my $model (@models) {
            if ($model->subject->name eq $subject_id) {
                push @return_models, $model;
            }
        }
    }
    if(scalar(@return_models) != scalar(@fields)) {
        die "can't match this model group the the denovo output";
    }
    return @return_models;
}



sub output_passing_vcf {
    my($self, $denovo_vcf, $r_output, $pvalue_cutoff, $output_filename) = @_;
    my $pvalues = Genome::Sys->open_file_for_reading($r_output);
    my $output_file = IO::File->new($output_filename, ">");
    my $cat_cmd = "cat";
    if(Genome::Sys->file_is_gzipped($denovo_vcf)) {
        $cat_cmd = "zcat";
    }


    my @lines = `$cat_cmd $denovo_vcf | grep "^#" `;
    $output_file->print(@lines);
    my (@unaffected, $chr, $pos, $affected);
    while(my $line = $pvalues->getline) {
        ($chr,$pos, $unaffected[0], $unaffected[1], $affected) = split "\t", $line; ###generic method for any kind of family later
        my $output=1;
        for my $unaff_pvalue (@unaffected) {
            if($unaff_pvalue < $pvalue_cutoff) {
                $output=0;
            }
        }
        if($output) {
            my $line = `$cat_cmd $denovo_vcf | grep "^$chr\t$pos" `;
            $output_file->print($line);
        }
    }
}



sub prepare_readcount_files {
    my $self = shift;
    my $ref_fasta = shift;
    my $sites_file = shift;
    my $model_ref = shift;
    my @readcount_files;
    my $qual = $self->min_read_qual;
    for my $model (@$model_ref) {
        my $readcount_out = Genome::Sys->create_temp_file_path($model->subject->name . ".readcount.output");
        my $bam = $model->last_succeeded_build->whole_rmdup_bam_file;
        push @readcount_files, $readcount_out;
        unless(-s $readcount_out) {
            print STDERR "running bam-readcount";
            my $rv = Genome::Model::Tools::Sam::Readcount->execute(
                bam_file => $bam,
                minimum_mapping_quality => $qual,
                output_file => $readcount_out,
                reference_fasta => $ref_fasta,
                region_list => $sites_file,
            );
            unless ($rv) {
                $self->error_message("Failed to run readcount");
                return;
            }
        }
    }
    return \@readcount_files;
}




1;
sub prepare_r_input {
    my $self = shift;
    my $vcf_fh=shift;
    my $readcount_files = shift;
    my ($r_input_fh,$r_input_path) = Genome::Sys->create_temp_file();
    while (my $line = $vcf_fh->getline) {
        next if ($line =~ m/^#/);
        chomp($line);
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split "\t", $line;
        my @possible_alleles = ($ref);
        my @alts = split ",", $alt;
        push @possible_alleles, @alts;
        my $novel=0;
        for (my $i = 0; $i < @samples; $i++) {
            my ($gt, $rest) = split ":", $samples[$i];
            my ($all1, $all2) = split "[|/]", $gt;
            unless(grep {/$all1/} @possible_alleles) {
                $novel = $all1;
            }
            unless(grep {/$all2/} @possible_alleles) {
                $novel = $all2;
            }

        }
        if($novel) {
            $r_input_fh->print("$chr\t$pos");
            for (my $i = 0; $i < @samples; $i++) {
                my $readcount_file = $readcount_files->[$i];
                chomp(my $readcount_line = `grep "^$chr\t$pos" $readcount_file`);
                my ($chr, $pos, $rc_ref, $depth, @fields) = split "\t", $readcount_line;
                my $prob;
                my ($gt, $rest) = split ":", $samples[$i];
                my @ref_bases = split "[|/]", $gt;
                if(grep /$novel/, @ref_bases) {
                    $prob = .5;
                }
                else {
                    $prob = .01;
                }
                my $ref_depth=0;
                my $var_depth=0;
                if($readcount_line) {
                    for my $fields (@fields) { 
                        my ($base, $depth, @rest)  = split ":", $fields;
                        next if($base =~m/\+/);
                        if($base eq $novel) {
                            $var_depth+=$depth;
                        }
                        else {
                            $ref_depth+=$depth;
                        }
                    }
                }
                $r_input_fh->print("\t$var_depth\t$ref_depth\t$prob");
            }
            $r_input_fh->print("\n");
        }
    }
    $r_input_fh->close;
    return $r_input_path;
}


sub r_code {
    my $self = shift;
    my $input_name = shift;
    my $output_name = shift;
    return <<EOS
options(stringsAsFactors=FALSE);
mytable=read.table("$input_name");
num_indices=(ncol(mytable)-2)/3
pvalue_matrix=matrix(nrow=nrow(mytable), ncol=num_indices);
indices = seq(from=3,to=num_indices*3, by=3)
for (i in 1:nrow(mytable)) {  
    for (j in indices) { 
        pvalue_matrix[i,(j/3)]=binom.test(as.vector(as.matrix(mytable[i,j:(j+1)])), 0, mytable[i,(j+2)], "t", .95)\$p.value; 
    }
}
chr_pos_matrix=cbind(mytable\$V1,mytable\$V2, pvalue_matrix);
write(t(chr_pos_matrix), file="$output_name", ncolumns=(num_indices+2), sep="\t");
EOS
}

