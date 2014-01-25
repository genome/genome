package Genome::Model::RnaSeq::Command::RawcountMatrix;

# Written by: Zachary Skidmore, Avinash Ramu

######################################################################
####################### set up modules ###############################
######################################################################

use strict;
use warnings;
use Genome;

#######################################################################
####################### set up input parameters #######################
#######################################################################

class Genome::Model::RnaSeq::Command::RawcountMatrix
{
    is => 'Genome::Command::Base',
    has_input =>
    [
        models =>
            {
            is => 'Genome::Model::RnaSeq',
            is_many => 1,
            shell_args_position => 1,
            doc => 'RNAseq models to generate a raw count matrix.',
            },
        output =>
        {
            doc => 'The output tsv file of raw count values obtained from HTseq tool per model.',
        },
        edgeR =>
        {
            doc => 'An optional Boolean, if true (i.e. 1) performs a simple differential expression analysis',
            is_optional => 1,
        },
    ],
};

######################################################################
###################### set up various help sub-routines ##############
######################################################################

sub help_synopsis
{
    return <<"EOS"
    genome model rna-seq rawcount-Matrix --output=FILE.tsv --edgeR=(1 or 0)(optional)  MODEL_GROUP_ID
EOS
}

sub help_brief
{
    return "accumulate Raw Count values obtained from (presumably) the HTseq tool into a matrix, optionally performs DE analysis with edgeR.";
}

sub help_detail
{
    return <<EOS
    Accumulate raw count values for genes and isoforms across a group of RNAseq models, optionally perform a default analysis with edgeR (assumes groups are tumor and normal) tumor is compared to normal.
EOS
}

######################################################################
########### set up main program that generates matrix ################
######################################################################

sub execute
{

    my $self = shift;

    my @models = $self->models;

    # declare/initalize various variables

    my @builds;
    my $annotation_build;
    my $reference_build;
    
    # loop through the models in the group and perform various quality checks

    foreach my $model (@models)
    {

        # die if a model does not have a succesfull build, otherwise build array with build id's

        my $build = $model -> last_succeeded_build;

        unless($build)
        {
            die('No successful builds were found for: '. $model->id);
        }

        push @builds, $build;

        # die if annotation builds for models do not match

        my $model_annotation_build = $model->annotation_build;

        if($annotation_build)
        {
            unless($annotation_build->id eq $model_annotation_build->id)
            {
                die('Mis-match annotation builds!');
            }

        } else {

                $annotation_build = $model_annotation_build;
    
        }

        # die if the reference sequence builds do not match

        my $model_reference_sequence_build = $model->reference_sequence_build;

        if($reference_build)
        {
            unless($reference_build->id eq $model_reference_sequence_build->id)
            {
                die("Mis-match reference sequence builds!")
            } 
        } else {

                $reference_build = $model_reference_sequence_build;

        }

    }


    ######################################################################
    ################ command to create the matrix ########################
    #######################################################################
    

    my %annotation_hash;

    # get the path to the .gtf (annotation file)
    
    my $gene_gtf_path = $annotation_build->annotation_file('gtf',$reference_build->id);

    $self->status_message('Loading known genes annotation file: '. $gene_gtf_path);

    # open the gtf file and read it creating a hash of ensemble_id => gene name

    open(GTF, "$gene_gtf_path") || die "Can not open path to file: $!";

    while(<GTF>)
    {
        chomp $_;
        my @tmp = split("\t", $_);
       
    # extract just the gene name from the annotation file

        $tmp[8] =~ /gene_name\s*"[A-Za-z0-9._-]*"/;
        my $gene_name = $&;
        $gene_name =~ s/gene_name\s*"//g;
        $gene_name =~ s/"//g;

    # extract just the ensemble id from the annotation file

        $tmp[8] =~ /gene_id\s*"\w+"/;
        my $ensemble_id = $&;
        $ensemble_id =~ s/gene_id\s*"//g;
        $ensemble_id =~ s/"//g;

    # build the hash
        
        $annotation_hash{$ensemble_id} = "$gene_name\t";
    }

    $self->status_message('There are '. scalar(keys %annotation_hash). ' genes in the annotation file: '. $gene_gtf_path);

    # for each rna-seq model open the raw counts file and take out ensemble_id and raw_counts
    # then if the ensemble_id exist in the annotation hash created above append the raw_counts
    # to the value of the hash, else kill the program as annotation file does not match data file

    my @subject;
    my $subject;
    my $subject_a;
    my $subject_b;

    foreach my $build (@builds)
    {

        # set up @subject array to hold subjects for each build ID to use as header later

         $subject_a = $build->subject->common_name;
         $subject_b = $build->subject->name;
         $subject = $subject_a . $subject_b;
         push(@subject, $subject);

         my $gene_count_tracking = $build->data_directory .'/results/digital_expression_result/gene-counts.tsv';

         $self->status_message('Loading raw count file: '. $gene_count_tracking);

         unless(-e $gene_count_tracking)
         {
             die("Failed to find: $gene_count_tracking\n");
         }

        open(DATA_FILE, "$gene_count_tracking") || die "Can not open $gene_count_tracking: $!";

        # read the DATA_FILE, ignore anything that does not have an ensemble id (i.e. no feature, etc.), elsif the key exist in the hash append the raw count
        # information to the value of the key, else give a warning that the key is not in the hash

        while(<DATA_FILE>)
        {
            chomp $_;
            my @tmp2 = split("\t", $_);
            my $gene_id = $tmp2[0];
            my $raw_count = $tmp2[1];

            if($gene_id eq 'no_feature' || $gene_id eq 'ambiguous' || $gene_id eq 'too_low_aQual' || $gene_id eq 'not_aligned' || $gene_id eq 'alignment_not_unique')
            {
                next;

            } elsif(exists $annotation_hash{$gene_id}) {

                $annotation_hash{$gene_id} .= "$raw_count\t";
            
            } else {

                print "\n!!!!ALERT!!!! The Following Ensemble id $gene_id is not found in the Annotation file: $annotation_hash{$gene_id}\n";
            }
        }

        close(DATA_FILE);
    }

    # write the headers and the hash to a file in a .tsv format

    my $output_file = $self->output;

    open(OUTPUT, ">$output_file\_tmp") || die "Can not open $output_file\_tmp";

    print OUTPUT "Ensemble_id\tGene_name";

    foreach(@subject)
    {
        print OUTPUT "\t$_";
    }

    print OUTPUT "\t\n";


    while( my ($key, $value) = each %annotation_hash)
    {
        print OUTPUT "$key\t$value\n";
    }

    close OUTPUT;

    # reformat the output file to get rid of the extra tabs at the end of the file

    $self->status_message('Reformating the output file');

    # read the tmp output file

    open(OUTPUT, "$output_file\_tmp") || die "Can not open $output_file\_tmp";

    # overwrite the original output file using a regexp to remove the unecessary tab on each line

    open(NEW_OUTPUT, ">$output_file") || die "Can not open $output_file";

    while(<OUTPUT>)
    {
        $_ =~ s/\t$//;
        print NEW_OUTPUT "$_";
    }

    unlink "$output_file\_tmp";

    close (NEW_OUTPUT);
    close (OUTPUT);

    # call R depending on optional parameter edgeR

    my $boolean = $self -> edgeR;

    if($boolean)
    {

        $self->status_message('Performing default differential expression analysis with edgeR');

        my $r_script_path = $self ->__meta__->module_path;
        $r_script_path =~ s/\.pm/\.R/;

        my $r_cmd = 'Rscript '. $r_script_path .' '. $output_file;
        Genome::Sys->shellcmd(cmd => $r_cmd, input_files => [], output_files => [],); 

    }

    $self->status_message('Finished!');

    return 1;

}

1;
