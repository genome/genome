package Genome::Model::Tools::ManualReview::CreateVariantReviewList;

use strict;
use warnings;
use Data::Dumper;
use IO::File;
use Genome::Utility::VariantReviewListReader;

use Genome;

class Genome::Model::Tools::ManualReview::CreateVariantReviewList
{
    is => 'Command',
    has => [
        list => {
            is => 'String',
            doc => 'Master csv to import', 
        },
        list_name => { 
            is => 'String', 
            doc => 'List name',
        },
        subject_name => {
            is => 'String',
            doc => 'Subject name',
            is_optional => 1,
        },
    ],
    has_optional => [
        separation_character =>{
            is => 'String',
            doc => 'character or string that separates the fields in the list',
            default => '|',
        },
        author => { 
            is => 'String', 
            doc => 'Author or authors of the list',
        },
        rt_ticket => { 
            is => 'String', 
            doc => 'RT ticket id number',
        },
        #we need to add subject name as a parameter
    ]
};

sub help_brief{
    return "Adds a Variant Review List to the Database";
}

sub help_synopsis{
    return "gmt manual-review create-variant-review-list --list <list> --name <list_name> --author <list_author> ";
}

sub help_detail{
    return "This command takes in a character delimited list of variants. The separation character is the '|' char by default, but can be specified on the command line. The list should be formatted as follows:\n".join('|', Genome::Utility::VariantReviewListReader->list_columns)."\n Once uploaded, the variants can be viewed and edited online at this address: $ENV{GENOME_SYS_SERVICES_FILES_URL}/view/variant_review_list.html
";
}

sub execute {
    my $self = shift;
    my $name = $self->list_name;
    my $author = $self->vtest($self->author);
    my $rt_ticket = $self->vtest($self->rt_ticket);
    my $separation_character = $self->separation_character;
    my $subject_name = $self->subject_name || $self->list_name;
    
    my @props =sort { 
        $a->property_name cmp $b->property_name
    } grep {
        $_->column_name ne ''
    } Genome::VariantReviewDetail->get_class_object->all_property_metas;
    my @det_property_names = map { $_->property_name } @props;
    eval{
        my $review_list = Genome::VRList->get(name => $name);
        if ($review_list){
            $self->error_message("A list with this name($name) already exists in the DB! id:".$review_list->id);
            return;
        }else{
            $review_list = Genome::VRList->create(name => $name, author => $author, rt_ticket => $rt_ticket);            
        }
        my $list_id = $review_list->id;
        my $list_reader = Genome::Utility::VariantReviewListReader->new($self->list, $separation_character);
        my $i=1;
        while (my $line_hash = $list_reader->next_line_data()){
            last unless $line_hash;
            next if $line_hash->{header};
            print "adding line: $i.\n";
            $i++;
            
            my $info;
            my %temp_line_hash;
            foreach my $col (@det_property_names)
            {
                $temp_line_hash{$col} = $line_hash->{$col} if exists $line_hash->{$col};
            }
            $temp_line_hash{start_position} = $line_hash->{begin_position} if exists $line_hash->{begin_position};
            $temp_line_hash{stop_position} = $line_hash->{end_position} if exists $line_hash->{end_position};
            $temp_line_hash{stop_position} = $temp_line_hash{start_position} if !defined $temp_line_hash{stop_position};
            $line_hash = \%temp_line_hash;
            
            my $current_member = Genome::VariantReviewDetail->get( start_position => $line_hash->{start_position}, chromosome => $line_hash->{chromosome}, subject_name => $subject_name );
            #if(!$current_member) {print "Adding variant to database.\n";}
            if ($current_member){

                #genes, this should maybe be a hangoff table
                my $current_genes = $current_member->genes;
                my @current_genes;
                @current_genes = split(/[,\/]/, $current_genes) if $current_genes;
                my $new_genes = $line_hash->{genes};
                my @new_genes;
                @new_genes = split(/[,\/]/, $new_genes) if $new_genes;
                my @genes_to_add;
                foreach my $new_gene (@new_genes){
                    my $found = 0;
                    foreach my $current_gene (@current_genes){
                        $found++ if $new_gene =~ /$current_gene/i;
                    }
                    push @genes_to_add, $new_gene unless $found;
                }
                if (@genes_to_add){
                    $info .= "adding genes ".join(" ", @genes_to_add)."\n";
                }
                my @genes = (@current_genes, @genes_to_add);
                $current_member->genes(join(",", @genes));
                #rest
                foreach my $col (keys %$line_hash){
                    next if $col =~ /genes/;
                    my $current_val = $current_member->$col;
                    my $new_val = $line_hash->{$col};
                    if ($current_val){
                        if ($new_val){
                            if ($current_val ne $new_val){
                                $current_member->$col("DISCREPANCY($current_val : $new_val");
                                $info .= "discrepancy at $col : $current_val || $new_val\n";
                            }
                        }
                    }elsif ($new_val){
                        $current_member->$col($new_val);
                    }
                }
                print $info if defined $info;
            }else{
                $current_member = Genome::VariantReviewDetail->create(
                    %$line_hash, subject_name => $subject_name   
                );
            }
            my $member_id = $current_member->id;
            my $review_list_member = Genome::VRListMember->get_or_create(
                list_id => $list_id,
                member_id => $member_id,
            );
        }  #while ( my $line = getline);
    };
    if ($@){
        $self->error_message("error in execution. $@");
        return 0;
    }
    return 1;
} 

sub vtest{
    my ($self, $v) = @_;
    return $v? $v : '';
}

1;
