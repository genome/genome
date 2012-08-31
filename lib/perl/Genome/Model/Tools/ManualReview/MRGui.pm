package Genome::Model::Tools::ManualReview::MRGui;

use strict;
use warnings;

our $initialized = 0;

sub init_gtk {
    return if $initialized;
    eval {
        require Gtk2;
        Gtk2->init;
        require Gtk2::GladeXML;
        require Glib;
    };
    die $@ if $@;
    $initialized = 1;
}

use IO::File;
use Genome::Utility::VariantReviewListReader;

use File::Basename ('fileparse','basename');
use base qw(Class::Accessor);
use Time::HiRes qw(usleep);
use Cwd 'abs_path';

Genome::Model::Tools::ManualReview::MRGui->mk_accessors(qw(current_file g_handle re_g_handle header));

my %iub_hash = ( A => 1,
                C => 2,
                G => 3,
                T => 4,
                M => 5,
                K => 6,
                Y => 7,
                R => 8,
                W => 9,
                S => 10,
                D => 11,
                B => 12,
                H => 13,
                V => 14,
                N => 15,
);
my %rev_iub_hash = map { $iub_hash{$_},$_; } keys %iub_hash;

my %somatic_status = (WT => 1,
                      O => 2,
                      LQ => 3,
                      A => 4,
                      S => 5,
                      G => 6,
                      V => 7,
                      NC => 8,
                      LOH => 9,
);

my %rev_somatic_status = map { $somatic_status{$_},$_; } keys %somatic_status;

my $last_chrom='';
my $last_pos='';
my @last_vals;

sub new 
{
    my ($caller, %params) = @_;   
    init_gtk() unless $initialized;

    croak("__PKG__:new:no class given, quitting") if @_ < 1;
	my $mr_dir = 'Genome/Model/Tools/ManualReview';
    foreach my $path (@INC) {
        my $fullpath = $path . "/" .$mr_dir;
        if( -e $fullpath) {
            $mr_dir = $fullpath;
            last;
        }
    }
    
    chomp $mr_dir;
    
    $mr_dir .= '/manual_review.glade';
    
    my $glade = new Gtk2::GladeXML($mr_dir,"manual_review");
    
    $params{g_handle} = $glade;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = $class->SUPER::new(\%params);    
    $glade->signal_autoconnect_from_package($self);

    my $mainWin = $glade->get_widget("manual_review");
    my $treeview = $glade->get_widget("review_list");
    $self->build_review_tree;

    $mainWin->signal_connect("destroy", sub { Gtk2->main_quit; nuke_consed_from_orbit(); });   

    my $project_file = $params{project_file};
    $self->open_file(abs_path($project_file)) if($project_file && -e $project_file);  
    return $self;
}

sub open_file_dialog
{
	my ($self,@parms) = @_;
	my $fc = Gtk2::FileChooserDialog->new("Open File",undef, 'open',
	'gtk-cancel' => 'cancel',
	'gtk-ok' => 'ok');
  	$fc->set_select_multiple(0);
 	my $response = $fc->run;
    if($response eq 'ok')
    {
 	    my $file = $fc->get_filename;
  	    $fc->destroy;
	    $self->open_file($file);
    }
    else
    {
        $fc->destroy;
    }
}

my %cmp_hash = (x => 30, X => 30, y => 31, Y=>31);
for(my $i=0;$i<30;$i++) { $cmp_hash{$i}=$i; }
sub chrom_sort
{
    my ($liststore, $itera, $iterb) = @_;
    my ($a1, $a2) = $liststore->get($itera,1,2);
    my ($b1, $b2) = $liststore->get($iterb,1,2);
    if(!defined $a1 || !defined $b1) {return 0;}
    if($a1 eq $b1) 
	{
		return $a2 <=> $b2;
	}
	else
	{
        if(exists $cmp_hash{$a1} && exists $cmp_hash{$b1})
        {
		    return $cmp_hash{$a1} <=> $cmp_hash{$b1};
        }
        else
        {
            return $a1 cmp $b1;
        }
	}
    
}

sub build_review_tree
{
    my ($self) = @_;
    my $handle = $self->g_handle;
    
    my $tree = $handle->get_widget("review_list");

    my @col = (Gtk2::TreeViewColumn->new_with_attributes
                        ('Sample Name', Gtk2::CellRendererText->new, text => 0),
         Gtk2::TreeViewColumn->new_with_attributes
                        ('Chromosome', Gtk2::CellRendererText->new, text => 1),
         Gtk2::TreeViewColumn->new_with_attributes
                        ('Position', Gtk2::CellRendererText->new, text => 2),
         Gtk2::TreeViewColumn->new_with_attributes
                        ('Delete Sequence', Gtk2::CellRendererText->new, text => 3),
         Gtk2::TreeViewColumn->new_with_attributes
                        ('Insert Sequence Allele 1', Gtk2::CellRendererText->new, text => 4),
         Gtk2::TreeViewColumn->new_with_attributes
                        ('Insert Sequence Allele 2', Gtk2::CellRendererText->new, text => 5),
         Gtk2::TreeViewColumn->new_with_attributes
                        ('Reference Genotype', Gtk2::CellRendererText->new, text => 6),
         Gtk2::TreeViewColumn->new_with_attributes
                        ('Pass Manual Review', Gtk2::CellRendererText->new, text => 7),
         Gtk2::TreeViewColumn->new_with_attributes
                        ('Manual Genotype Variant', Gtk2::CellRendererText->new, text => 8),
         Gtk2::TreeViewColumn->new_with_attributes
                        ('Notes', Gtk2::CellRendererText->new, text => 9),
         Gtk2::TreeViewColumn->new_with_attributes
                        ('Somatic Status', Gtk2::CellRendererText->new, text => 10),
         Gtk2::TreeViewColumn->new_with_attributes
                        ('Data Needed', Gtk2::CellRendererText->new, text => 11),
                                                
        );
    foreach (@col)
    {
        $tree->append_column($_);
    }
    
    $tree->signal_connect("row-activated", \&open_cb, $self);

    return $tree;


}

sub set_project_consedrc {
    my $edit_dir =shift;
    my $rc ="$edit_dir/.consedrc";

    my $wide = 'consed.alignedReadsWindowInitialCharsWide: 120';
    my $expand = 'consed.alignedReadsWindowAutomaticallyExpandRoomForReadNames: false';

    my %rc_hash;
    my $i;
    if(-e $rc ){
        open F, $rc or warn("Failed to open consed rc ($rc).\n") and return; 
        while(<F>){
            chomp;
            $i++;
            $rc_hash{$_}=$i;
        }
        close F;

    }

    my @to_append;
    unless(exists $rc_hash{$wide}){
        push @to_append, $wide;
    }
    unless(exists $rc_hash{$expand}){
        push @to_append, $expand;
    }

    return unless @to_append;

    open F, ">>$rc" or warn("Failed to open consed rc ($rc).\n") and return; 
    for(@to_append){
        print F $_ ."\n";
    }
    close F;

    return 1;
}

sub get_contig_start_pos
{
    my ($self,$ace_file_name) = @_;
    my $reader = Genome::Model::Tools::Consed::AceReader->create(file => $ace_file_name);
    die $self->error_message('Failed to create ace reader for '.$ace_file_name) if not $reader;
    while(my $line = $reader->_fh->getline)
    {
        if($line =~ /CT\{/)
        {
            $reader->_fh->seek(-length($line),1);
            my $tag = $reader->next;
            if($tag->{type} eq 'startNumberingConsensus')
            {
                my ($pos) = split (/\n/,$tag->{data});
                return $pos;
            }
            
        }
    }
}

sub open_consed
{
    my ($self, $proj_name, $relative_target_base_pos) = @_;
    
    my $consed= 'consed';
    my($file, $dir)= fileparse($self->current_file);
    my $suffix = '.1';
    my $edit_dir = "$dir/$proj_name/edit_dir";
    if(-e $edit_dir."/consedSocketLocalPortNumber")
    {
        unlink $edit_dir."/consedSocketLocalPortNumber";
    }
    my $pid = fork();
    if ($pid)
    {
        my $i=0;
        while(!-e $edit_dir."/consedSocketLocalPortNumber")
        {
            usleep 100000;
            $i++;
            if($i>50)
            {
                warn "Timed out waiting for consed to launch.\n";
                last;
            }
        }
        unlink $edit_dir."/consedSocketLocalPortNumber";
        
        return 1;
    }
    
    if(!defined $pid) {print "fork unsuccessful.\n"; }
    
    my $ace1 = "$proj_name.ace$suffix";
    $ace1 = "$proj_name.ace" unless (-e $edit_dir.'/'.$ace1);
    #print "opening ",$edit_dir.'/'.$ace1,"\n";
    my $contig_start_pos = $self->get_contig_start_pos($edit_dir.'/'.$ace1);
    if($contig_start_pos)#don't assume consed starts at 1, convert coords if necessary
    {
        $relative_target_base_pos = $relative_target_base_pos - $contig_start_pos + 1;
    }
    else
    {
        $relative_target_base_pos = 1001;# if(!defined $relative_target_base_pos);    
    }
    if( -d $edit_dir){
        chdir $edit_dir or die "can't cd to $edit_dir"; 
        #$ace1 = "$proj_name.ace" unless (-e $ace1);
        unless(-e $ace1){
            print "ERROR no: $ace1 ... skipping (report to apipe)\n"; 
            exit;
        }        
        set_project_consedrc($edit_dir);

        my $out = `grep SABBOTT $ace1`;# hack that Lynn added, may remove later, since we can just use the position
        #that is provided in the review editor
        chomp $out;
        if($out){

            (undef, undef, undef, $relative_target_base_pos) =split(/\s+/,$out);

            $relative_target_base_pos =&resolve_padded_pos($ace1, $relative_target_base_pos);
        }
        #end hack
        
        my $c_command= "$consed -socket 0 -ace $ace1 -mainContigPos $relative_target_base_pos";
        #print $c_command,"\n";
        open (STDOUT, '>/dev/null');
        open (STDERR, '>/dev/null');
#        *STDOUT = *LOG;
#        *STDERR = *LOG;
        my $rc = system($c_command);
        if($rc)
        {
            warn "Failed to launch consed.\n";
            #print "just finished opening consed\n";
            `touch consedSocketLocalPortNumber`;
        }
    }else{
        warn("ERROR  Can't find $edit_dir ... skipping (report to sabbott)\n");
    }
    exit;
    
}

sub resolve_padded_pos{
    my ($ace, $rel)=@_;

     my $obj =  AceFile->new($ace);
         die unless $obj;

    my $c;
    my @contigs =$obj->get_valid_contigs;
    if(@contigs == 1){
        $c = $contigs[0];
    }else{
        my @grep =grep{$_ == 1}@contigs;
        if(@grep){
            $c=1;
        }else{
            die "can't find valid contig in ace $ace";
        }
    }
    my %h =$obj->unpadded_positions_for_contig($c);
    unless(%h){
        warn "not able to get unpadded_positions for contig $c";
        return;
    }

    return  $h{$rel};
}

sub get_col_order
{
    my ($self) = @_;
    my $header = $self->header;
    my @vis_cols = (
        'sample_name',
        'chromosome',
        'start_position',
        'delete_sequence',
        'insert_sequence_allele1',
        'insert_sequence_allele2',
        'genotype_iub_code',
        'pass_manual_review',
        'manual_genotype_iub_variant',
        'notes',
        'somatic_status',
        'data_needed',
    );
    my %vis_cols = map { $_,1; } @vis_cols;
    
    foreach my $col (@$header)
    {
        push @vis_cols, $col unless(exists $vis_cols{$col});    
    }
    return @vis_cols;
}

sub open_file
{
	my ($self,$file) = @_;

	return unless (-e $file);
	$self->current_file($file);    

    my $handle = $self->g_handle;
    my $mainWin = $handle->get_widget("manual_review");
    my $tree = $handle->get_widget("review_list");

    my $list_reader = Genome::Utility::VariantReviewListReader->new($self->current_file, '|');
    $self->header($list_reader->{separated_value_reader}->headers);
     
    my @col_order = $self->get_col_order(@{$self->header});
    my $col_count = scalar @col_order;
    my $model = Gtk2::ListStore->new(('Glib::String')x$col_count);
    $tree->set_model($model);
    $model->set_sort_func (0, \&chrom_sort);
    while (my $line_hash = $list_reader->next_line_data())
    {
        last unless $line_hash;
        if ($line_hash->{header}) { print $line_hash->{header},"\n"; }
        next if $line_hash->{header};
        my $iter = $model->append;
        
        #check if we at least have a chromosome and position for the variant, if not,
        #then we are clearly looking at invalid data (or have a parser error)
        if(!defined $line_hash->{$col_order[1]} || !defined $line_hash->{$col_order[2]})
        {
            $tree->set_model(Gtk2::ListStore->new(('Glib::String')x$col_count));
            my $dialog = Gtk2::MessageDialog->new ($mainWin,
                                      'destroy-with-parent',
                                      'error', # message type
                                      'ok', # which set of buttons?
                                      "There was an error opening the csv file.'\n");
            my $response = $dialog->run;
  
            $dialog->destroy;
            return;        
        }
        for(my $i = 0;$i<@col_order;$i++)
        {
            $model->set($iter,
            $i => delete $line_hash->{$col_order[$i]});
        }            
    }
    $model->set_sort_column_id(0,'GTK_SORT_ASCENDING');
    
    my $project_file = $self->current_file;
    my $filename = basename($project_file);
    my $uid = $>;
    my $uname = getpwuid($uid);
    $project_file.=".$uname" unless($filename =~ /$uname/);
    $self->current_file($project_file);
    $filename = basename($project_file);
    $mainWin->set_title("Manual Review - $filename");
}

sub open_cb
{
	my ($tree, $mpath, $col, $self)=@_;
    
    $self->on_review_button_clicked;    
}

sub save_file
{
    my ($self) = @_;
    my $fh = IO::File->new(">".$self->current_file);
    my @col_order = $self->get_col_order;

    my $tree = $self->g_handle->get_widget("review_list");
    my $model = $tree->get_model;
    my $header = join '|',@col_order;
    $header .= "\n";
    print $fh $header;
    my $iter = $model->get_iter_first;
    no warnings;
    do
    {   
        my @cols = $model->get($iter);
        my $row = join '|',@cols;
        $row .= "\n";
        print $fh $row;    
    }
    while($iter = $model->iter_next($iter));
    use warnings;
    $fh->close;
}

sub on_save_file
{
	my ($self) = @_;
	
	unless(defined $self->current_file)
	{
		$self->save_file_as_dialog;
	}
	if(defined $self->current_file)
	{
		$self->save_file;
	}	
}

sub save_file_as_dialog
{
	my ($self) = @_;
	my $fc = Gtk2::FileChooserDialog->new("Save File As",undef, 'save',
	'gtk-cancel' => 'cancel',
	'gtk-ok' => 'ok');
  	$fc->set_select_multiple(0);
	if($self->current_file)
	{
  		$fc->set_current_folder ($self->current_file);
 	}
	my $response = $fc->run;
 	my $file = $fc->get_filename;
  	$fc->destroy;
    
	$self->current_file($file) if(defined $file);	
}

sub on_save_file_as
{
	my ($self) = @_;
	
	$self->save_file_as_dialog;
	if(defined $self->current_file)
	{
        my $handle = $self->g_handle;
        my $mainWin = $handle->get_widget("manual_review");        
        my $filename = basename($self->current_file);
        $mainWin->set_title("Manual Review - $filename");
		$self->save_file;
	}
}

sub on_re_ok
{
    my ($button, $data) = @_;
    my ($self,$glade,$review_editor,$model, $row) = @{$data};
    my %pf = (Pass => 1, Fail => 2);
    my %rpf = (1 => 'Pass', 2 => 'Fail');
    my %dn = (Yes => 1, No => 2);
    my %rdn = (1 => 'Yes', 2 => 'No');
    
    my $cb = $glade->get_widget('re_genotype');
    my $active = $cb->get_active();
    $model->set($row,6 => $rev_iub_hash{$active}) if exists $rev_iub_hash{$active};
    $cb = $glade->get_widget('re_passfail');
    $active = $cb->get_active();
    $model->set($row,7 => $rpf{$active}) if exists $rpf{$active};
    $cb = $glade->get_widget('re_genotype_variant');
    $active = $cb->get_active();
    $model->set($row,8 => $rev_iub_hash{$active}) if exists $rev_iub_hash{$active};
    $cb = $glade->get_widget('re_somatic_status');
    $active = $cb->get_active();
    $model->set($row,10 => $rev_somatic_status{$active}) if exists $rev_somatic_status{$active};
    $cb = $glade->get_widget('re_data_needed');
    $active = $cb->get_active();
    $model->set($row,11 => $rdn{$active}) if exists $rdn{$active};
    
    my $tb = $glade->get_widget('re_text_view');    
    my ($s,$e) = $tb->get_buffer->get_bounds;
    my $text = $tb->get_buffer->get_text($s,$e,0);
    $text =~ tr/\n\t/  /;
    $model->set($row, 9 => $text) if defined $text;
    #added later this destroys it too
    $self->on_review_editor_destroy($review_editor);    
    return 1;
}

sub on_review_editor_destroy
{
    my ($self, $review_editor) = @_;

    #saving window location
    my ($x, $y) = $review_editor->get_position;    
    $self->{last_x} = $x;
    $self->{last_y} = $y;
    $review_editor->destroy;
    $last_chrom = -1;
    $last_pos = -1;
    @last_vals = ();
    nuke_consed_from_orbit();
    return 1;
}

sub on_review_button_clicked
{
    my ($self) = @_;
    my $g_handle = $self->g_handle;
    my $review_list = $g_handle->get_widget("review_list");
    my $mainWin = $g_handle->get_widget("manual_review");
    if(!defined $review_list)
    {
        my $dialog = Gtk2::MessageDialog->new ($mainWin,
                                      'destroy-with-parent',
                                      'error', # message type
                                      'ok', # which set of buttons?
                                      "There are no reviews to edit, please open a review list.\n");
        my $response = $dialog->run;
        $dialog->destroy;
        return 1;
    }
    my ($path,$column) = $review_list->get_cursor;
    if(!defined $path || !defined $column)
    {
        my $dialog = Gtk2::MessageDialog->new ($mainWin,
                                      'destroy-with-parent',
                                      'error', # message type
                                      'ok', # which set of buttons?
                                      "Please select a review to edit first, then click Review.\n");
        my $response = $dialog->run;
        $dialog->destroy;
        return 1;    
    }
    my $mr_dir = 'Genome/Model/Tools/ManualReview';
    foreach my $path (@INC) {
        my $fullpath = $path . "/" .$mr_dir;
        if( -e $fullpath) {
            $mr_dir = $fullpath;
            last;
        }
    }    
    
    chomp $mr_dir;
    
    $mr_dir .= '/manual_review.glade';
    my $glade = new Gtk2::GladeXML($mr_dir,"review_editor");
    $self->re_g_handle($glade);
    #$glade->signal_autoconnect_from_package($self);
    my $review_editor = $glade->get_widget("review_editor");
    if(defined $self->{last_x} && defined $self->{last_y})
    {
        $review_editor->move($self->{last_x},$self->{last_y});
    }
    
    $review_editor->signal_connect("delete_event", sub { $self->on_review_editor_destroy($review_editor); });
    my $model = $review_list->get_model;
    
    my $row = $model->get_iter($path);
    my @val = $model->get($row,6,7,8,9,10,11);    
    my $re_hpaned = $glade->get_widget("re_hpaned");
    my ($width) = $review_editor->get_size_request;
    $re_hpaned->set_position($width/2.5);
    my $cancel = $glade->get_widget("re_cancel");
    $cancel->signal_connect("clicked", sub { $self->on_review_editor_destroy($review_editor); });
    my $ok = $glade->get_widget("re_ok");
    $ok->signal_connect("clicked", \&on_re_ok,[$self, $glade,$review_editor, $model, $row]);
    my $prev = $glade->get_widget("re_prev");
    $prev->signal_connect("clicked", sub {$self->on_prev_button_clicked($glade);});
    my $next = $glade->get_widget("re_next");
    $next->signal_connect("clicked", sub {$self->on_next_button_clicked($glade);});
    #set widgets

    my %pf = (Pass => 1, Fail => 2);
    my %dn = (Yes => 1, No => 2);
    my $cb = $glade->get_widget('re_genotype');
    $cb->set_active($iub_hash{$val[0]}) if(defined $val[0] && exists $iub_hash{$val[0]});
    $cb = $glade->get_widget('re_passfail');
    $cb->set_active($pf{$val[1]}) if(defined $val[1] && exists $pf{$val[1]});
    $cb = $glade->get_widget('re_genotype_variant');
    $cb->set_active($iub_hash{$val[2]}) if(defined $val[2] && exists $iub_hash{$val[2]});
    my $tb = $glade->get_widget('re_text_view');
    $tb->get_buffer->set_text($val[3]) if $val[3];
    $cb = $glade->get_widget('re_somatic_status');
    $cb->set_active($somatic_status{$val[5]}) if(defined $val[5] && exists $somatic_status{$val[5]}); 
    $cb = $glade->get_widget('re_data_needed');
    $cb->set_active($dn{$val[6]}) if(defined $val[6] && exists $dn{$val[6]});
    
    @val = $model->get($row,1,2);
    my $row_hash = $self->get_row_hash($model,$row);
    my $proj_dir;
    if(exists $row_hash->{project_type} &&
       $row_hash->{project_type} eq 'amplicon')
    {
        $proj_dir = $val[0];
    }
    else
    {
        $proj_dir = join '_',@val;
    }
    $self->open_consed($proj_dir,$val[1]);
    return $review_editor;
}

sub get_row_hash
{
    my ($self,$model, $row) = @_;
    my @col_order = $self->get_col_order;
    my @val = $model->get($row);
    my %row_hash;
    @row_hash{@col_order} = @val;
    return \%row_hash;

}

sub display_row
{
    my ($self,$model, $row, $glade) = @_;
    
    my @val = $model->get($row,6,7,8,9,10,11);
    #set widgets
    my ($chrom,$pos) = $model->get($row,1,2);

    if($last_chrom == $chrom && $last_pos == $pos)
    {
        
        if(@last_vals >1 && !scalar(grep {defined $_;} @val))
        {
            @val = @last_vals;
            #for(my$i = 0;$i<@last_vals;$i++)
            #{
            #    $val[$i] = $last_vals[$i];
            #}
        }    
    }

    my %pf = (Pass => 1, Fail => 2);
    my %dn = (Yes => 1, No => 2);
    my $cb = $glade->get_widget('re_genotype');
    if(defined $val[0] && exists $iub_hash{$val[0]})
    {
        $cb->set_active($iub_hash{$val[0]});        
    }
    else
    {
            $cb->set_active(0);        
    }
    $cb = $glade->get_widget('re_passfail');
    
    if(defined $val[1] && exists $pf{$val[1]})
    {
        $cb->set_active($pf{$val[1]});
    }
    else
    {
        $cb->set_active(0);
    }
    $cb = $glade->get_widget('re_genotype_variant');
    
    if(defined $val[2] && exists $iub_hash{$val[2]})
    {
        $cb->set_active($iub_hash{$val[2]});
    }
    else
    {
        $cb->set_active(0);
    }    
    
    $cb = $glade->get_widget('re_somatic_status');
    if(defined $val[4] && exists $somatic_status{$val[4]})
    {
        $cb->set_active($somatic_status{$val[4]});
    }
    else
    {
        $cb->set_active(0);
    }
    $cb = $glade->get_widget('re_data_needed');
    if(defined $val[5] && exists $dn{$val[5]})
    {
        $cb->set_active($dn{$val[5]});
    }
    else
    {
        $cb->set_active(0);
    }
    
    my $tb = $glade->get_widget('re_text_view');
    $tb->get_buffer->set_text($val[3]||'');

    @val = $model->get($row,1,2);
    my $proj_dir = join '_',@val;
    unless($last_chrom == $chrom && $last_pos == $pos)
    {
        nuke_consed_from_orbit();
        #system "killall -s QUIT -q consed";
        $self->open_consed($proj_dir,$val[1]);
    }
    return ;
}

sub save_row
{
    my ($self,$model, $row, $glade) = @_;
    my %pf = (Pass => 1, Fail => 2);
    my %rpf = (1 => 'Pass', 2 => 'Fail');
    my %dn = (Yes => 1, No => 2);
    my %rdn = (1 => 'Yes', 2 => 'No');
    my $cb = $glade->get_widget('re_genotype');
    my $active = $cb->get_active();
    $model->set($row,6 => $rev_iub_hash{$active}) if exists $rev_iub_hash{$active};
    $last_vals[0] = $rev_iub_hash{$active};
    $cb = $glade->get_widget('re_passfail');
    $active = $cb->get_active();
    $model->set($row,7 => $rpf{$active}) if exists $rpf{$active};
    $last_vals[1] = $rpf{$active};
    $cb = $glade->get_widget('re_genotype_variant');
    $active = $cb->get_active();   
    
    $model->set($row,8 => $rev_iub_hash{$active}) if exists $rev_iub_hash{$active};
    $last_vals[2] = $rev_iub_hash{$active};
    $cb = $glade->get_widget('re_somatic_status');
    $active = $cb->get_active();
    $model->set($row,10 => $rev_somatic_status{$active}) if exists $rev_somatic_status{$active};
    $last_vals[4] = $rev_somatic_status{$active};
    $cb = $glade->get_widget('re_data_needed');
    $active = $cb->get_active();
    $model->set($row,11 => $rdn{$active}) if exists $rdn{$active};
    $last_vals[5] = $rdn{$active};
    
    my $tb = $glade->get_widget('re_text_view');
    my ($s,$e) = $tb->get_buffer->get_bounds;
    my $text = $tb->get_buffer->get_text($s,$e,0);
    $text =~ tr/\n\t/  /;
    $model->set($row, 9 => $text) if defined $text;    
    $last_vals[3] = $text;
    return 1;

}

sub on_prev_button_clicked
{
    my ($self, $glade) = @_;
    my $g_handle = $self->g_handle;
    my $review_list = $g_handle->get_widget("review_list");
    my $mainWin = $g_handle->get_widget("manual_review");
    if(!defined $review_list)
    {
        my $dialog = Gtk2::MessageDialog->new ($mainWin,
                                      'destroy-with-parent',
                                      'error', # message type
                                      'ok', # which set of buttons?
                                      "There are no reviews to edit, please open a review list.\n");
        my $response = $dialog->run;
        $dialog->destroy;
        return 1;
    }
    my ($path,$column) = $review_list->get_cursor;
    if(!defined $path || !defined $column)
    {
        my $dialog = Gtk2::MessageDialog->new ($mainWin,
                                      'destroy-with-parent',
                                      'error', # message type
                                      'ok', # which set of buttons?
                                      "Please select a review to edit first, then click Review.\n");
        my $response = $dialog->run;
        $dialog->destroy;
        return 1;    
    }

    my $model = $review_list->get_model;
    my $iter = $model->get_iter($path);  
    ($last_chrom,$last_pos) = $model->get($iter,1,2);  
    $self->save_row($model, $iter, $glade) if $glade; 
    return unless $path;
    return unless $path->prev;
    $iter = $model->get_iter($path);
    $self->display_row($model, $iter, $glade) if ($glade && $path && $iter);
    
    $review_list->set_cursor($path);

}

sub on_next_button_clicked
{
    my ($self, $glade) = @_;
    my $g_handle = $self->g_handle;
    my $review_list = $g_handle->get_widget("review_list");
    my $mainWin = $g_handle->get_widget("manual_review");
    if(!defined $review_list)
    {
        my $dialog = Gtk2::MessageDialog->new ($mainWin,
                                      'destroy-with-parent',
                                      'error', # message type
                                      'ok', # which set of buttons?
                                      "There are no reviews to edit, please open a review list.\n");
        my $response = $dialog->run;
        $dialog->destroy;
        return 1;
    }
    my ($path,$column) = $review_list->get_cursor;
    if(!defined $path || !defined $column)
    {
        my $dialog = Gtk2::MessageDialog->new ($mainWin,
                                      'destroy-with-parent',
                                      'error', # message type
                                      'ok', # which set of buttons?
                                      "Please select a review to edit first, then click Review.\n");
        my $response = $dialog->run;
        $dialog->destroy;
        return 1;    
    }
    
    my $model = $review_list->get_model;
    my $iter = $model->get_iter($path);
    ($last_chrom,$last_pos) = $model->get($iter,1,2);   
    $self->save_row($model, $iter, $glade) if $glade; 
    return unless $path;
    $path->next;
    $iter = $model->get_iter($path);
    $path->prev if(!$iter);
    $self->display_row($model, $iter, $glade) if ($glade && $path && $iter);
    $review_list->set_cursor($path);
    #set widgets

    return ;
}

sub nuke_consed_from_orbit #killall doesn't seem to work since the default variant of consed keeps changing it's name, i.e. consed04032007 blah
{
    my @procs = `ps -o pid,command=`;
    chomp @procs;
    @procs = grep { /consed/; } @procs;    
    foreach (@procs)
    {
        $_ =~ s/^\s+//;#trim leading white space
        my ($pid) = split /\s+/,$_;
        `kill -s 9 $pid`;
    }    
}

sub gtk_main_quit
{
    Gtk2->main_quit;
    nuke_consed_from_orbit();
    #system "killall -s QUIT -q consed";
}


1;
