package Genome::Model::Tools::ManualReview::UploadVariantReviewList;

use strict;
use warnings;
use Genome::Utility::VariantReviewListReader;
use Data::Dumper;


use Genome;

UR::Object::Type->define(
    class_name => __PACKAGE__, is => 'Command',
    has => [
        list => {
            is => 'String',
            doc => 'Master csv to update from', 
        },
        #project_dir => {
        #    is => 'String',
        #    doc => 'Directory of review project',
        #},
    ],
    has_optional => [
        db_list_name =>{
            is => 'String',
            doc => 'Name of list to be updated',
        },
        db_list_id => { 
            is => 'String', 
            doc => 'ID of list to be updated',
        },
        separation_character => {
            is => 'String',
            doc => 'separation character in list of variants to be added',
        },
    ]
);

sub help_brief{
    return "Add new variants to an existing list";
}

sub help_synopsis{
    return "gmt manual-review add-variants-to-list --list <list> --db_list_id <id>";
}

sub help_detail{
"Adds variants to existing lists given the list name or list id. The list must be a character delimited list.  The default delimiter is the '|' char, but this can be specified on the command line.  The list columns must be in the following format:\n".join('|', Genome::Utility::VariantReviewListReader->list_columns)."\n Once uploaded, the variants can be viewed and edited online at this address: https://gscweb.gsc.wustl.edu/view/variant_review_list.html";
}

BEGIN {
    my @det_property_names;
    my @props =
        sort { $a->property_name cmp $b->property_name }
        grep { $_->column_name ne '' }
        grep { defined($_->column_name) }
        Genome::VariantReviewDetail->get_class_object->all_property_metas;
    @det_property_names = map { $_->property_name } @props;
    
    my @rev_property_names;
    my @rev_props =
        sort { $a->property_name cmp $b->property_name }
        grep { $_->column_name ne '' }
        grep { defined($_->column_name) }
        Genome::SNVManualReview->get_class_object->all_property_metas;
    @rev_property_names = map { $_->property_name } @rev_props;
    sub fix_hash_data
    {
        my ($line_hash) = @_;
        my %temp_line_hash;
        foreach my $col (@det_property_names)
        {
            $temp_line_hash{$col} = $line_hash->{$col} if exists $line_hash->{$col};
        }
        $temp_line_hash{start_position} = $line_hash->{begin_position} if exists $line_hash->{begin_position};
        $temp_line_hash{stop_position} = $line_hash->{end_position} if exists $line_hash->{end_position};
        $temp_line_hash{stop_position} = $temp_line_hash{start_position} if !defined $temp_line_hash{stop_position};
        $line_hash = \%temp_line_hash;

        return $line_hash;

    }
    sub get_review_data
    {
        my ($self, $line_hash) = @_;
        my %temp_line_hash;
        foreach my $col (@rev_property_names)
        {
            if( exists $line_hash->{$col})
            {
                if ($col eq 'pass_manual_review' or $col eq 'data_needed'){
                    my ($match) = $line_hash->{$col} =~ /^([PFYN])/;
                    $temp_line_hash{$col} = $match;                    
                }else{
                    $temp_line_hash{$col} = $line_hash->{$col};
                }
            }
            if((($col eq 'pass_manual_review') || ($col eq 'genotype_iub_code'))&& 
               !(defined ($temp_line_hash{$col}) && length ($temp_line_hash{$col}||'')))
            {
                $temp_line_hash{$col} = '-';
#                print "setting $col to $temp_line_hash{$col}\n";
            }
        }
        $line_hash = \%temp_line_hash;
        return $line_hash;
    }
}

sub execute{
    my $self = shift;
    my $list = Genome::Utility::VariantReviewListReader->new($self->list, $self->separation_character);
    my $db_list = Genome::VRList->get($self->db_list_name ? (name=>$self->db_list_name) : (id=>$self->db_list_id)); 

    unless ($db_list){
        $self->error_message("List doesn't exist");
        return 0;
    }
    my $list_id = $db_list->id;
    my ($detail) = $db_list->details;
    my $subject_name = $detail->subject_name;
    while (my $line_hash = $list->next_line_data()){
        last unless $line_hash;
        next if $line_hash->{header};

        my $det_hash = fix_hash_data($line_hash);
        my $current_member = Genome::VariantReviewDetail->get( start_position => $det_hash->{start_position}, chromosome => $det_hash->{chromosome}, subject_name => $subject_name );
        unless ($current_member){
            print "can't get member for pos=> ".$det_hash->{start_position}." chromosome => ".$det_hash->{chromosome}." subject_name => ".$subject_name;
        die;
        }
#        print Dumper $det_hash;
        my $rev_hash = $self->get_review_data($line_hash);
        next unless $rev_hash;
        my @file_name_split = split(/\./, $self->list);
        my $reviewer_name = $file_name_split[-1];
        my $review = Genome::SNVManualReview->get_or_create(detail_id => $current_member->id, dump_date => "2008-12-11 13:50:59", build_id => 93385989, reviewer => $reviewer_name,%{$rev_hash});
#        my $review = Genome::SNVManualReview->get_or_create(detail_id => $current_member->id, dump_date => "2008-12-11 13:50:59", build_id => 1, reviewer => $ENV{USERNAME});
#        $review->set(%$rev_hash);
#        print Dumper $rev_hash;
    }  
    return 1;
}

1;
