#:boberkfe this probably should use a templating engine rather than html mixed in

package Genome::Model::ReferenceAlignment::Report::Mapcheck;

use strict;
use warnings;

use Genome;

use CGI;
use IO::String;

class Genome::Model::ReferenceAlignment::Report::Mapcheck {
    is => 'Genome::Model::Report',
    has => [
            name => { default_value => 'Mapcheck' },
            description => {
                calculate => q|
                return '<div>maq mapcheck coverage for ' . $self->model->name . " as of " . UR::Context->current->now.'</div>';
                |,
            },
    ],
};


sub default_format {
    'html'
}

sub available_formats {
    return ('html')
}

sub _add_to_report_xml {
    my $self = shift;

    return {
        #description => $self->generate_report_brief,
        html => $self->generate_report_detail,
    };
}

sub generate_report_brief {
    my $self = shift;

    my $model = $self->model;
    return '<div>maq mapcheck coverage for ' . $model->name . " as of " . UR::Context->current->now.'</div>';
}

sub generate_report_detail {
    my $self       = shift;
    my $dedup_name = $self->model->duplication_handler_name;

    if ($dedup_name eq 'maq') { 
        return $self->get_maq_content;
    } 
    else {
        return $self->get_bam_content;
    }
}

sub get_bam_content {
    my $self = shift;
    my $build = $self->build;

    #my $pileup_file = $build->bam_pileup_file;
    #$self->debug_message("Using pileup file $pileup_file to generate Bam coverage.");
    #my $coverage = Genome::Model::Tools::Sam::Coverage->create(pileup_file=>$pileup_file);
    my $ref_build = $build->model->reference_sequence_build;
    my $reference_file = $ref_build->full_consensus_path('fa');
    my $aligned_reads = $build->whole_rmdup_bam_file;
    unless (-s $aligned_reads){
        $self->error_message("aligned reads file $aligned_reads doesn't exist!");
    }
    unless (-s $reference_file){
        $self->error_message("reference file $reference_file doesn't exist!");
    }

    my $flagstat_file = "$aligned_reads.flagstat";
    unless (-s $flagstat_file){
        $self->warning_message("No flagstat file found for merged file $aligned_reads, proceeding w/o checking aligned_reads count");
    }

    my $flagstat_data = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_file);
        
    unless($flagstat_data) {
        $self->error_message('No output from samtools flagstat');
        return;
    }
    
    my $skip;
    if(exists $flagstat_data->{errors}) {
        for my $error (@{ $flagstat_data->{errors} }) {
            if($error =~ m/Truncated file/) {
                $self->error_message('Flagstat output for ' . $aligned_reads . ' indicates possible truncation.');
            }
        }
    }else{
        if ($flagstat_data->{reads_mapped} == 0){
            $self->debug_message("0 mapped reads in merged bam, skipping bam-check");
            $skip = 1;
        }
    }

    if ($skip){
        $self->debug_message("No bam coverage report generated due to no aligned reads");
        return;
    }else{
        $self->debug_message("Using:  Reference File: $reference_file, Aligned Reads File: $aligned_reads");
        my $coverage = Genome::Model::Tools::Sam::Coverage->create( aligned_reads_file=>$aligned_reads,
            reference_file =>$reference_file,
            return_output =>1,
            coverage_command => $ENV{GENOME_SW} . '/samtools/bamcheck/bamcheck-v0.13/bam-check -q 1', 
            #use_version => 'r350wu1', #lift this once set default to r350wu1 in G::M::T::Sam
        );
        my $bam_coverage_report = $coverage->execute;
        if (defined($bam_coverage_report) ) {
            $self->debug_message("Bam coverage report successfully generated.");
            $self->debug_message("Bam coverage report string: \n".$bam_coverage_report);
            return $bam_coverage_report;
        }  else {
            $self->error_message("Could not generate Bam coverage report.");  
            return;
        }
    }
}

sub get_maq_content {
    my $self = shift;

    my $model = $self->model;
    my $build = $self->build;

    my $mapcheck_output = Genome::Sys->create_temp_file_path('mapcheck');
    #my $accumulated_alignments_file = $build->accumulate_maps;
    my $accumulated_alignments_file = $build->whole_rmdup_map_file;
    unless (Genome::Model::Tools::Maq::Mapcheck->execute(
            use_version => $build->maq_version_for_pp_parameter,
            bfa_file => $model->reference_sequence_build->full_consensus_path('bfa'),
            map_file => $accumulated_alignments_file,
            output_file => $mapcheck_output,
        ) ) {
        $self->error_message();
        die($self->error_message);
    }
    #my $rm_cmd = "rm $accumulated_alignments_file";
    #$self->shellcmd(cmd => $rm_cmd);

    my $mapcheck_output_fh = Genome::Sys->open_file_for_reading($mapcheck_output);
    my @maq = $mapcheck_output_fh->getlines;

    my $rpt = join('',@maq);
    $rpt = $self->format_maq_content($rpt);

    my $fh = IO::String->new();
    $fh->print($rpt);

    $fh->seek(0 ,0);
    return join('', $fh->getlines);
}

sub get_coverage_filename {
    my $self = shift;
    my $reports_dir = $self->build->resolve_reports_directory;
    #my $reports_dir = $self->model->resolve_reports_directory;
    my $model = $self->model;
    return $reports_dir . '/' .  $model->genome_model_id . '_coverage_detail.html';
}

sub format_maq_report
{
    #format plain text for html viewability
    my ($self,$content) = @_;

    my ($stats, $table);
    if ($content=~m/(.*)(\n\n)(.*)/sm)
    {
        ($stats, $table) = ($1, $3); 
        $stats=~s/\n/\<br>\n/g; 
        $stats = "<div id=\"stats\">$stats</div>";

        my @table = split("\n",$table);
        for (my $row = 0, my $cell, my $formatted_cell = ''; $row < scalar(@table); $row++, $formatted_cell = '')
        {   
            $cell = $table[$row];
            #trim leading & trailing
            $cell=~ s/^\s+//;
            $cell=~ s/\s+$//;

            #wrap cells
            if ($row == 0) #header
            {
                $cell=~s/\s*:\s/ /g;
                $cell=~s/\s+/<\/th><th>/g;
                $formatted_cell="<tr><td id=\"corner\"></td><th>$cell</th></tr>";
            }
            else
            {
                #color-code colon-delimited sections
                if ($cell=~m/(\s*)(\d+)(\s*)(.*?)(\s*:\s*)(.*?)(\s*:\s*)(.*?)(\s*:\s*)(.*)/)
                {
                    my (@sections) = ($2, $4, $6, $8, $10);
                    for (my $i = 0, my $sec; $i < scalar(@sections); $i++)
                    {
                        $sec = $sections[$i];
                        $sec=~s/\s+/<\/td><td class=\"sec$i\">/g;
                        $sec = "<td class=\"sec$i\">$sec</td>";
                        $formatted_cell .= $sec;
                    }
                }
                $formatted_cell = "<tr>$formatted_cell</tr>";
            } 
            $table[$row] = $formatted_cell;
        }

        $table = join('',@table);
        $table = "\n<table border=1 id=\"data\">$table</table>";

        return "<!--\n$content\n-->\n" . 
        "<div id=\"maq_report\">" .
        $stats . 
        $table . 
        "</div>" .   
        $self->get_style;
    }
    else
    {
        die("Expected format STATS \n\n TABLE");
    }
}

sub format_maq_content
{
    my ($self,$content) = @_;

    my ($stats, $table, @table);

    if ($content=~m/(.*)(\n\n)(.*)/sm)
    {
        ($stats, $table) = ($1, $3); 
        $stats=~s/\n/\<br>\n/g; 
        $stats = "<div id=\"stats\">$stats</div>";

        @table = split("\n",$table);
        for (my $row = 0, my $cell, my $formatted_cell = ''; $row < scalar(@table); $row++, $formatted_cell = '')
        {   
            $cell = $table[$row];
            #trim leading & trailing
            $cell=~ s/^\s+//;
            $cell=~ s/\s+$//;

            #wrap cells
            if ($row == 0) #header
            {
                $cell=~s/\s*:\s/ /g;
                $cell=~s/\s+/<\/th><th>/g;
                $formatted_cell="<tr><td id=\"corner\"></td><th>$cell</th></tr>";

            }
            else
            {
                #color-code colon-delimited sections
                if ($cell=~m/(\s*)(\d+)(\s*)(.*?)(\s*:\s*)(.*?)(\s*:\s*)(.*?)(\s*:\s*)(.*)/)
                {
                    my (@sections) = ($2, $4, $6, $8, $10);
                    for (my $i = 0, my $sec; $i < scalar(@sections); $i++)
                    {
                        $sec = $sections[$i];
                        $sec=~s/\s+/<\/td><td class=\"sec$i\">/g;
                        $sec = "<td class=\"sec$i\">$sec</td>";
                        $formatted_cell .= $sec;
                    }

                    $formatted_cell = "<tr>$formatted_cell</tr>";
                }
            } 

            $table[$row] = $formatted_cell;
        }
    }
    $table = join('',@table);
    $table = "\n<table border=1 id=\"data\">$table</table>";

    return "<!--\n$content\n-->\n" . 
    "<div id=\"maq_report\">" .
    $stats .
    $table . 
    "</div>" . 
    $self->get_css;
}

sub get_css
{    
    return
    "    <style>

    #maq_report #data #corner
    {
    border-bottom: 2px solid #6699CC;
    border-right: 1px solid #6699CC;
    border-top, border-left:#000000;
    background-color: #BEC8D1;
    }

    th
    { border-bottom: 2px solid #6699CC;
    border-left: 1px solid #6699CC;
    background-color: #BEC8D1;
    text-align: center;
    font-family: Verdana;
    font-weight: bold;
    font-size: 16px;
    color: #404040; }

    td
    { border-bottom: 1px solid #000;
    border-top: 0px;
    border-left: 1px solid #000;
    border-right: 0px;
    font-family: Verdana, sans-serif, Arial;
    font-weight: normal;
    font-size: 14px;
    padding: 0 0 0 0;

    background-color: #fafafa;
    border-spacing: 0px;
    margin-top: 0px;
    }


    table
    {
    text-align: center;
    font-family: Verdana;
    font-weight: normal;
    font-size: 14px;
    color: #404040;
    background-color: #fafafa;
    border: 1px #000 solid;
    border-collapse: collapse;
    border-spacing: 0px;
    }

    tr td{
    background-color: #fafafa;
    }

    tr.row0 td{
    background-color: #fafafa;
    }

    tr.row1 td{
    background-color: #eeeeee;
    }
    a img {
    border: medium none;
    border-collapse: collapse;
    }
    .drag_it {
    }
    .sec0
    { border-bottom: 2px solid #6699CC;
    border-right: 1px solid #6699CC;
    background-color: #BEC8D1;
    text-align: center;
    font-family: Verdana;
    font-weight: bold;
    font-size: 16px;
    color: #404040; }
    .sec1
    {
    background-color:#CCFFFF;
    }
    .sec2
    {
    background-color:#FFFF99;
    }
    .sec3
    {
    background-color:#FF9999;
    }
    .sec4
    {
    background-color:#CCFFCC;
    }

    </style>";
}
1;
