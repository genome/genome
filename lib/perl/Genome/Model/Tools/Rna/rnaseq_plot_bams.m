function rnaseq_plots(input_file)

%clear all;
%Read input file
input_data = read_input_file_bams(input_file);

fid = fopen(input_data.alignment_file_paths,'r'); 
i = 1;
a = fgets(fid);
while(a > 0 )
    alignments{i} = strtok(a);
    a = fgets(fid);
    i = i +1;
end
fclose(fid);

fid = fopen(input_data.junctions_file_paths,'r'); 
i = 1;
a = fgets(fid);
while(a > 0 )
    junctions{i} = strtok(a);
    a = fgets(fid);
    i = i +1;
end
fclose(fid);

fid = fopen(input_data.transcripts_file_paths,'r'); 
i = 1;
a = fgets(fid);
while(a > 0 )
    transcripts{i} = strtok(a);
    a = fgets(fid);
    i = i +1;
end

fid = fopen(input_data.sample_name,'r'); 
i = 1;
a = fgets(fid);
while(a > 0 )
    sample_name{i} = strtok(a);
    a = fgets(fid);
    i = i +1;
end

fclose(fid);
if(strcmp(input_data.chr,'X') ||strcmp(input_data.chr,'Y')) 
    chr = input_data.chr;
else
    chr = str2num(input_data.chr);
end
start = str2num(input_data.start);
stop = str2num(input_data.stop);

%Configure plots
plot_number = 1;
hf1 = figure('Position',[100 100 1500 500],'name',num2str(plot_number));
subplot_count = 1;
hold on;
cm = colormap;

row = 1;
col = 1;
counter = 1;
for ii = 1:length(alignments)
    if(counter > str2num(input_data.subplots_col))
        row = row+1;
        col = 1;
        counter = 1;
    end
    %Create temporary output dir if does not exist
    if(~exist(input_data.temp_directory))
        mkdir(input_data.temp_directory);
    else
        delete([input_data.temp_directory,'/*']);
    end
    build_out = [input_data.temp_directory,'out.bam'];
    
    roi = [input_data.temp_directory,input_data.chr,'-',num2str(start),'-',num2str(stop),'.txt'];
    fid = fopen(roi,'w'); %create locus file
    if(strcmp(input_data.chr,'X') ||strcmp(input_data.chr,'Y') )
    fprintf(fid, '%s\t%d\t%d',chr,start,stop);
    else
    fprintf(fid, '%d\t%d\t%d',chr,start,stop);
    end
    fclose(fid);
    
    if(plot_number > str2num(input_data.subplots_row) + str2num(input_data.subplots_col)+1)
      saveas(gcf,[input_data.output_directory,strcat(input_data.chr,':',input_data.start,'-',input_data.stop),'_',num2str(subplot_count)],'eps');
      saveas(gcf,[input_data.output_directory,strcat(input_data.chr,':',input_data.start,'-',input_data.stop),'_',num2str(subplot_count)],'fig');
      subplot_count = subplot_count + 1;
        plot_number = 1;
        row = 1;
        col = 1;
        hf1 = figure('Position',[100 100 1500 500],'name',num2str(plot_number));
    end
    ha1 = subplot(str2num(input_data.subplots_row),str2num(input_data.subplots_col),plot_number);
    hold on;
    range = [0 1 0 1];
    hight = range(4) - range(3);
    width = range(2) - range(1);
    ha1 = gca;
    axis(ha1,range);
    if(strcmp(input_data.chr,'X') ||strcmp(input_data.chr,'Y'))
        xlabel(ha1,['chr',chr,': ',num2str(start),'-',num2str(stop)],...
        'fontsize',8,'fontweight','bold');
    else
        xlabel(ha1,['chr',num2str(chr),': ',num2str(start),'-',num2str(stop)],...
        'fontsize',8,'fontweight','bold');
    end
    set(gca,'xtick',[]);
    ylabel(ha1,'read counts','fontsize',8,'fontweight','bold');
    %title(ha1,[sample_name{ii},' alignment: ',num2str(plot_number)],'fontweight','bold');
    title(ha1,sample_name{ii},'fontweight','bold');
    
    %Get chunk of bam file for the current loci to be annalyszed
    if(strcmp(input_data.chr,'X') ||strcmp(input_data.chr,'Y') )
        system(['samtools view -b ',alignments{ii},' ', chr,':',num2str(start),'-',num2str(stop),' > ',build_out]);
    else
        system(['samtools view -b ',alignments{ii},' ', num2str(chr),':',num2str(start),'-',num2str(stop),' > ',build_out]);    
    end

    %Make index file for bam chunk
    system(['samtools index ',build_out]);
    %If flagstat does not exist create it
    %     if(~exist([build_out,'.flagstat']))
    %         system(['samtools flagstat ',build,' | cut -d " " -f 1 > ',build_out,'.flagstat']);
    %     end
    
    %Perl scripts to pase files to plot
    gene_annotation_pl = '/gscmnt/gc2105/info/medseq/jbrea/rnaseq/scripts/gene_annotation.pl';
    trans_annotation_pl = '/gscmnt/gc2105/info/medseq/jbrea/rnaseq/scripts/trans_annotation.pl';
    repeatmask_pl = '/gscmnt/gc2105/info/medseq/jbrea/rnaseq/scripts/repeatmask.pl';
    junctions_pl = '/gscmnt/gc2105/info/medseq/jbrea/rnaseq/scripts/junctions.pl';
    cufflinks_pl = '/gscmnt/gc2105/info/medseq/jbrea/rnaseq/scripts/cufflinks.pl';
    
    if(strcmp(input_data.bamreadcounts,'yes'))
        %norm = load([build_out,'.flagstat']);
        norm = 1;
        system(['bam-readcount',' ','-q 1',' -l ',roi,' ',build_out,' ','| cut -f 2,4 >',input_data.temp_directory,'bamreadcount.txt']);
        r = load([input_data.temp_directory,'bamreadcount.txt']);
        if(size(r) > 0)
            r(:,2) = r(:,2)/norm;
            plot(r(:,1),r(:,2),'d','markersize',2);
            range = [start stop 0 max(r(:,2))+1];
            hight = range(4) - range(3);
            width = range(2) - range(1);
            axis(ha1,range);
        else
            range = [start stop 0 1];
            hight = range(4) - range(3);
            width = range(2) - range(1);
            axis(ha1,range);
        end
    end
    
    if(strcmp(input_data.gene_annotation,'yes'))
        vertical_loc = 0.32*hight;
        %perl(gene_annotation_pl,roi,input_data.gene_annotation_file,input_data.temp_directory);
        a = load(input_data.gene_annotation_file);
        if(size(a) > 0)
            for i = 1:length(a(:,1))
                if(a(i,4)==1)
                    plot(ha1,[a(i,2),a(i,3)],[vertical_loc,vertical_loc],...
                        'linewidth',5,'color',[0.7,0.7,0.7]);
                else
                    plot(ha1,[a(i,2),a(i,3)],[vertical_loc,vertical_loc],...
                        'linewidth',5,'color',[0.0,0.0,0.0]);
                end
                axis(ha1,range);
            end
        end
    end
    
    if(strcmp(input_data.transcript_annotation,'yes'))
        vertical_loc_min = 0.70*hight;
        vertical_loc_inc = 0.05*hight;
        %perl(trans_annotation_pl,roi,input_data.transcript_annotation_file,input_data.temp_directory);
        t = read_trans_annot_gtf(input_data.transcript_annotation_file);
        text(stop-0.5*width,range(4)-0.05*hight,t(1).gene_id);
        ht = 1;
        genes = {'null',0};
        while(length(t))
            exon_size_x = [str2double(t(1).start), str2double(t(1).stop)];
            exon_size_y = [(vertical_loc_min+vertical_loc_inc*ht) (vertical_loc_min+vertical_loc_inc*ht)];
            exon_links_y = [];
            exon_links_x = [];
            del = [];
            temp = {genes{1:2:end}};
            not_found = 1;
            for g = 0:length(genes)/2-1
                if(strcmp(genes{2*g+1},t(1).gene_id))
                    genes{2*g+2} = genes{2*g+2}+1;
                    ht = genes{2*g+2};
                    not_found = 0;
                    break
                end
            end
            if(not_found)
                genes = [genes,t(1).gene_id,1];
            end
            for k = 2:length(t)
                if(strcmp(t(1).transcript_id,t(k).transcript_id))
                    exon_size_x = [exon_size_x;[str2double(t(k).start), str2double(t(k).stop)]];
                    exon_size_y = [exon_size_y;[vertical_loc_min+vertical_loc_inc*ht vertical_loc_min+vertical_loc_inc*ht]];
                    exon_links_x = [exon_links_x;[str2double(t(1).stop), str2double(t(k).start)]];
                    exon_links_y = [exon_links_y;[vertical_loc_min+vertical_loc_inc*ht vertical_loc_min+vertical_loc_inc*ht]];
                    del = [del;k];
                end
            end
            t(del) = [];
            if(length(exon_links_y) > 0)
                if(str2double(t(1).exon_type_number) == 0)
                    line_width = 2;
                    color = [0.0,0.0,0.0];
                elseif(str2double(t(1).exon_type_number) == 1)
                    line_width = 4;
                    color = [0.0,0.0,0.0];
                elseif(str2double(t(1).exon_type_number) == 2)
                    line_width = 6;
                    color = [0.0,0.0,0.0];
                elseif(str2double(t(1).exon_type_number) == 3)
                    line_width = 6;
                    color = [0.5,0.5,0.5];
                end
                
                line(exon_size_x',exon_size_y','linewidth',line_width,'color',color);
                line(exon_links_x',exon_links_y','linewidth',2,'color',color);
                ht  = ht + 1;
            end
            t(1) = [];
        end
    end
    
    
    if(strcmp(input_data.repeat_mask,'yes'))
        vertical_loc = 0.40*hight;
       % if(strcmp(input_data.chr,'X') ||strcmp(input_data.chr,'Y') )
       %     perl(repeatmask_pl,roi,[input_data.mask_annotation_dir,'chr',chr,'_rmsk.txt'],input_data.temp_directory);
       % else
       %     perl(repeatmask_pl,roi,[input_data.mask_annotation_dir,'chr',num2str(chr),'_rmsk.txt'],input_data.temp_directory);
       % end
        al = load(input_data.mask_annotation_file);
        if(length(al(:,1)) ~= 0)
            for i = 1:length(al(:,1))
                plot(ha1,[al(i,1),al(i,2)],[vertical_loc vertical_loc],'linewidth',5,'color',[0.5,0.5,0.5]);
                axis(ha1,range);
            end
            
        end
    end
    
    if(strcmp(input_data.cufflinks_junctions,'yes'))
        %perl(junctions_pl,roi,junctions{ii},input_data.temp_directory);
j = load(junctions{ii});
        if(size(j))
            vertical_loc_min = 0.25*hight;
            vertical_loc_max = 0.30*hight;
            for i = 1:length(j(:,1))
                if(j(i,4) > 0)
                    tempx = [j(i,2) j(i,2)];
                    tempy = [vertical_loc_min vertical_loc_max];
                    line(tempx,tempy,'linewidth',2,'color',cm(round((63*j(i,4)/max(j(:,4))))+1,:));
                    tempx = [j(i,3) j(i,3)];
                    tempy = [vertical_loc_min vertical_loc_max];
                    line(tempx,tempy,'linewidth',2,'color',cm(round((63*j(i,4)/max(j(:,4))))+1,:));
                    axis(ha1,range);
                end
            end
            
            colorbar_x = 0.30+0.28*(col-1);
            colorbar_y = 0.86-0.47*(row-1);
            colobar_width = 0.01;
            colorbar_height = 0.05;
            caxis([0 max(j(:,4))]);
ch = colorbar('position',[colorbar_x colorbar_y colobar_width colorbar_height],'fontsize',8,'ytick',[0:max(j(:,4)):max(j(:,4))]);
set(get(ch,'ylabel'),'String', 'JUNC', 'Rotation', 270,'VerticalAlignment', 'top','fontsize',6,'FontWeight','bold')
        end
    end
    
    if(strcmp(input_data.cufflinks_transcript_asembly,'yes'))
        %perl(cufflinks_pl,roi,transcripts{ii},input_data.temp_directory);
c = load(transcripts{ii});
        if(size(c) > 1)
            vertical_loc_min = 0.45*hight;
            vertical_loc_inc = 0.04*hight;
           %fpkm_norm = zeros(length(c(:,1)),1);
           fpkm_max = 0;
	    if(c(1,2) == c(2,2))
	      if(c(1,3) == c(2,3))
              %fpkm_norm(c(1,2)) = max(fpkm_norm(c(1,2)),c(1,8));
              fpkm_max = max(fpkm_max,c(1,8));
              end
            end
            
            for i = 2:length(c(:,1))-1
		   if(c(i,2) == c(i-1,2) || c(i,2) == c(i+1,2))
		     if(c(i,3) == c(i-1,3) || c(i,3) == c(i+1,3))
                    %fpkm_norm(c(i,2)) = max(fpkm_norm(c(i,2)),c(i,8));
                    fpkm_max = max(fpkm_max,c(i,8));
                   end
                end
            end
	      if(c(end,2) == c(end-1,2))
		 if(c(end,3) == c(end-1,3))
               % fpkm_norm(c(end,2)) = max(fpkm_norm(c(end,2)),c(end,8));
                fpkm_max = max(fpkm_max,c(i,8));
              end
            end
            %fpkm_max = max(fpkm_norm);
            ht = 1;
            for i = 1:length(c(:,1))
                x = [];
                y = [];
                for j =i+1:length(c(:,1))
                    if((c(i,2)==c(j,2)) && c(i,3)==c(j,3))
                        if((i > 1) && (c(i,3) > c(i-1,3)))
                            ht = ht + 1;
                        end
                        if((i > 1) && (c(i,2) ~= c(i-1,2)))
                            ht = 1;
                        end
                        x = [x;[c(i,7),c(j,6)]];
                        y = [y;[vertical_loc_min+vertical_loc_inc*ht vertical_loc_min+vertical_loc_inc*ht]];
                        break
                    end
                end
                if(length(x) ~= 0)
                    tempx1 = [c(i,6),c(i,7)];
                    tempx2 = [c(j,6),c(j,7)];
                    tempy = [vertical_loc_min+vertical_loc_inc*ht vertical_loc_min+vertical_loc_inc*ht];
                    line(tempx1,tempy,'linewidth',5,'color','k');
                    line(tempx2,tempy,'linewidth',5,'color','k');
                    %line(x',y','linewidth',2,'color',cm(round(63*c(i,8)/fpkm_norm(c(i,2)))+1,:));
                    line(x',y','linewidth',2,'color',cm(round(63*c(i,8)/fpkm_max)+1,:));
                end
            end
            axis(ha1,range);
            %colorbar_x = 0.325+0.28*(col-1);
            colorbar_x = 0.15+0.28*(col-1);
            colorbar_y = 0.86-0.47*(row-1);
            colobar_width = 0.01;
            colorbar_height = 0.05;   
            caxis([0 fpkm_max]); 
           
            ch = colorbar('position',[colorbar_x colorbar_y colobar_width colorbar_height],'fontsize',8,'ytick',[0:fpkm_max:fpkm_max]);
            set(get(ch,'ylabel'),'String', 'FPKM', 'Rotation', 270,'VerticalAlignment', 'top','fontsize',6,'FontWeight','bold')
        end     
    end
    col = col+1;
    counter = counter + 1;
    plot_number = plot_number + 1;
end
saveas(gcf,[input_data.output_directory,strcat(input_data.chr,':',input_data.start,'-',input_data.stop),'_',num2str(subplot_count)],'eps');
saveas(gcf,[input_data.output_directory,strcat(input_data.chr,':',input_data.start,'-',input_data.stop),'_',num2str(subplot_count)],'fig');



function struct = read_input_file_bams(file)

fid = fopen(file);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.alignment_file_paths] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.junctions_file_paths] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.transcripts_file_paths] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.chr,x] = strtok(x);
[struct.start,x] = strtok(x);
[struct.stop,x] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.subplots_row] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.subplots_col] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.bamreadcounts] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.gene_annotation] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.transcript_annotation] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.repeat_mask] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.cufflinks_junctions] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.cufflinks_transcript_asembly] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.temp_directory] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.output_directory] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.gene_annotation_file] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.transcript_annotation_file] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.mask_annotation_file] = strtok(x);
x = fgets(fid);
[temp,x] = strtok(x);
[struct.sample_name] = strtok(x);
fclose(fid);

function struct_array = read_trans_annot_gtf(file)
    fid = fopen(file);
    x = fgets(fid);
    i = 1;
    while(x > 0)
       [struct.chr, x] = strtok(x);
       struct.chr = struct.chr(4:end);
       [struct.source,x] = strtok(x);
       [struct.exon_type, x] = strtok(x);
       [struct.start,x] = strtok(x);
       [struct.stop, x] = strtok(x);
       [struct.strand, x] = strtok(x);
       [struct.exon_type_number, x] = strtok(x);
       [gene_id_label,x] = strtok(x);
       [struct.gene_id, x] = strtok(x);
       struct.gene_id = struct.gene_id(2:end-2);
       [trascript_id_label,x] = strtok(x);
       [struct.transcript_id, x] = strtok(x);
       struct.transcript_id = struct.transcript_id(2:end-2);
       x = fgets(fid);
       struct_array(i) = struct;
       i = i + 1;
    end

