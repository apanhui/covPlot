#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2018-11-18 21:24:36    |
#-----------------------------------------------+
# name: asMap.pl
# func: 
# version: 1.0

use strict;
use warnings;

use SVG;
use Config::General;
use Getopt::Long;
use Math::Round;
use File::Basename qw/basename dirname/;
use FindBin qw($Bin);
use List::Util qw/min max sum/;

use lib "$Bin/lib";
use Font;
use CONF qw/load_conf/;
use General qw/timeLOG WARN ERROR/;

use lib "/home/aipeng/work/develepment/SBV/lib";
use SBV::STAT qw(dividing bezier3_xy);

use lib "/Bio/User/aipeng/bin/OSGO/lib";
use Pattern;

my $samtools = "/Bio/bin/samtools-1.9";

#-------------------------------------------------------------------------------
#  set the options
#-------------------------------------------------------------------------------
my %OPTS = (conf=>"$Bin/covPlot.conf");
GetOptions(\%OPTS,
    'conf:s','gene:s','gtf:s','beds:s','labels:s','colors:s','juncs:s',
    'upstream:i','downstream:i',
    'intron_fix:i','intron_scale:s','intron_background:s',
    'show_heads','show_tails',
    'strand','methy_type:s',
    'width:i','from_top:s','help');

&usage if ($OPTS{'help'});

unless ($OPTS{conf}){
    WARN("[FATA ERROR] you must defined the conf file!\n");
    &usage;
}

# load conf file 
my $conf = load_conf($OPTS{conf});
check_conf($conf);

# re defined some important options as var
my $legend = $conf->{legend};
my $gene   = $conf->{gene};

#-------------------------------------------------------------------------------
# load font 
my $font_path = "/Bio/User/aipeng/bin/OSGO/fonts/afm/";
load_font($font_path);

#-------------------------------------------------------------------------------
#  fetch the data
#-------------------------------------------------------------------------------
# fetch the coordinate
my %loci = fetch_gene_loci($conf->{gtf},$conf->{gene});
my ($chr,$start,$end,$strand) = ($loci{chr},$loci{start},$loci{end},$loci{strand});
timeLOG "fetch the gene coordinate done ... ";

# fetch the depth
my ($depths,$samples);
my %styles;

# read the bam files from ARGV
if (@ARGV){
    ($depths,$samples) = fetch_depths_from_ARGV(\@ARGV,$chr,$start,$end);
     %styles = init_styles($conf,$samples);
# read the bam files from conf file
}else{
    ($depths,$samples) = fetch_depths_from_conf($conf,$chr,$start,$end);
     %styles = parse_styles($conf,$samples);
}

my %depths  = %{$depths};
my @samples = sort { $samples->{$a}{order} <=> $samples->{$b}{order} } keys %$samples;

# reset the samples order which defined in conf file
@samples = split /,/ , $conf->{samples} if ($conf->{samples});

timeLOG "fetch the depths done ... ";

# fetch the junction
my %juncs = fetch_junctions($conf,$samples,$chr,$start,$end,$conf->{gene});
my $max_junc = fetch_max_junc(%juncs);
timeLOG "fetch the junctions done ... ";

# fetch marker info
my @marker = fetch_marker($conf,$chr,$start,$end);

#-------------------------------------------------------------------------------
# prepare to draw
#-------------------------------------------------------------------------------
my @flags   = split /,/ , $conf->{element_flags};
my %colors  = (upstream=>"#000000",five_prime_utr=>"#C01C30","CDS"=>"#29B473",
                three_prime_utr=>"#2A3890",exon=>"#000000",downstream=>"#000000",
                mRNA=>"#000000",gene=>"#000000");
my %heights = (upstream=>8,five_prime_utr=>14,CDS=>18,three_prime_utr=>14,exon=>18,downstream=>8,mRNA=>18,gene=>18);
%colors  = fetch_styles($conf->{colors},%colors);
%heights = fetch_styles($conf->{heights},%heights);

@flags = reverse @flags if ($strand eq "-");

my %flags = exon_or_intron(\%loci,\@flags,$conf);
my %coord = rebuild_gene_coord(\%loci,\@flags,$conf);
my $count = $coord{$end};
timeLOG "rebuild the region coordinate system done ... ";

# calc depths for each samples with window
my @regions = fetch_regions(\%flags,$start,$end);
my @exon_regions = @regions;

$exon_regions[0]  = $loci{gene_start};
$exon_regions[-1] = $loci{gene_end};

my %rdepths;
my $all_max_depth = 0;
for my $sample (@samples) {
    ERROR("${sample}'s info is not defined!") unless $styles{$sample}{window};
    my @temp = calc_depths($depths{$sample},$styles{$sample}{window},\@regions,\%coord);
    $rdepths{$sample}{depths} = \@temp;
    my @nums = map { $_->[1] } @temp;
    my $max = max(@nums);
    my $min = min(@nums);
    $rdepths{$sample}{max_depth} = max($max,abs($min));
    $rdepths{$sample}{min_depth} = $min;
    $all_max_depth = $rdepths{$sample}{max_depth} if ($all_max_depth < $rdepths{$sample}{max_depth});
}

#-------------------------------------------------------------------------------
#  start to draw AS plot 
#-------------------------------------------------------------------------------
# define some default options
my $font = Font->new("font-size:12");
my $font_height = $font->fetch_text_height;
my $font_style = $font->toStyle;

my $sr_font = Font->new("font-size:12");
my $sr_font_h = $sr_font->fetch_text_height;
my $sr_font_style = $sr_font->toStyle;

my $margin          = $conf->{size}->{margin}            // 20;
my $module_margin   = $conf->{size}->{module_margin}     // 20;
my $ctrl_height     = $conf->{size}->{ctrl_height}       // 20;
my $samples_spacing = $conf->{size}->{samples_spacing}   // 10;
my $lw              = $conf->{size}->{legend_item_width} // 40; # legend item width
my $spacing         = $conf->{size}->{spacing}           // 6;
my $width           = $conf->{width};
my $height          = 1000;

my $max_flag_width = $font->fetch_max_text_width(\@flags);
my $legend_width = $legend ?  $max_flag_width + $spacing*4 + $lw : 0;
my $sum_width = $width + $legend_width;

my $svg = SVG->new(width=>$sum_width,height=>$height,id=>"as_map");
my $ox = $margin;
my $oy = $margin;

## create the defs.
my $defs = $svg->defs(id=>"defs_1");

## create backgroud color for text with a filter
my $textbg = $defs->filter(x=>0,y=>0,width=>1,height=>1,id=>"textbg");
$textbg->fe(-type=>"flood","flood-color"=>"white");
$textbg->fe(-type=>"composite",in=>"SourceGraphic");

## groups for intron backgroud
my $intron_bg = $svg->g(id=>"intron_bg");

# init some coord
my @tsids = keys %{$loci{loci}};
my $gene_width = $font->fetch_max_text_width( [$gene,@tsids] );
my $max_sample_width = $font->fetch_max_text_width(\@samples);
my $max_tick_width = $all_max_depth != 1 ? $font->fetch_text_width( $all_max_depth ) : $font->fetch_text_width("0.2");
my $left_width = $spacing;
if ($conf->{samples_label_pos} eq "left"){
    my $tmp = $max_tick_width + $max_sample_width + $spacing;
    $left_width += $tmp > $gene_width ? $tmp : $gene_width;
}else{
    my $tmp = $max_tick_width;
    $left_width += $tmp > $gene_width ? $tmp : $gene_width;
}

my $chr_width = $font->fetch_text_width($chr);
my $end_width = $font->fetch_text_width($end);
my $gox = $ox + $left_width;
my $gex = $width - $margin;
my $tsh = $heights{transcript} // 24;

# init height 
$height = $oy;

if ($conf->{show_transcirpts}){

## draw ref genes
# draw chr
    $svg->text(x=>$gox-$chr_width-$spacing,y=>$oy+$font_height,style=>$font_style)->cdata($chr);
    $svg->text(x=>$gox,y=>$oy+$font_height,style=>$font_style)->cdata($start);
    $svg->text(x=>$gex-$end_width,y=>$oy+$font_height,style=>$font_style)->cdata($end);

# draw gene 
    my $geneid_width = $font->fetch_text_width($gene);
    $svg->text(x=>$gox-$spacing-$geneid_width,y=>$oy+$font_height+$spacing+$font_height,style=>$font_style)->cdata($gene);
    my $gene_arrowid = "gene_arrow";
    end_arrow($defs,$gene_arrowid,10,6,-zh=>0,-style=>"stroke-width:1;stroke:#000;fill:#000");

    my $goy = $oy + $font_height + $spacing + $font_height/2;

    my ($gx1,$gx2) = ($gox,$gex);
    ($gx1,$gx2) = ($gx2,$gx1) if ($strand eq "-");
    $svg->line(x1=>$gx1,x2=>$gx2,y1=>$goy,y2=>$goy,style=>"stroke-width:2;stroke:#000000",
            'marker-end'=>"url(#$gene_arrowid)");

# draw transcript
    my $toy = $goy + $font_height/2 + $spacing*3;
    foreach my $tsid (sort keys %{$loci{loci}}){
        my $ts_width = $font->fetch_text_width($tsid);
        $svg->text(x=>$gox-$spacing-$ts_width,y=>$toy+$tsh/2+$font_height/2,style=>$font_style)->cdata($tsid);
        $svg->line(x1=>$gox,x2=>$gex,y1=>$toy+$tsh/2,y2=>$toy+$tsh/2,style=>"stroke-width:1;stroke:#000000;stroke-dasharray:4 4 2 4");

        my $lastx = 0;
        foreach my $flag (@flags)
        {
            next unless $loci{loci}{$tsid}{$flag};
            my @recoreds = @{$loci{loci}{$tsid}{$flag}};
            for my $record (@recoreds)
            {
                my ($rs,$re) = @$record;
                $lastx = $rs if ($rs > $lastx);
            }
        }

        foreach my $flag (@flags)
        {
            next unless $loci{loci}{$tsid}{$flag};
            my @recoreds = @{$loci{loci}{$tsid}{$flag}};
            my $color = $colors{$flag};
            my $sw = $heights{$flag}/2;

            my @order_recoreds = sort {$a->[0] <=> $b->[0]} @recoreds;
            foreach my $record (@order_recoreds)
            {
                my ($rs,$re) = @$record;
                my $startx = $gox + fetch_axis_dis($coord{$rs},1,$count,$gex-$gox);
                my $endx = $gox + fetch_axis_dis($coord{$re},1,$count,$gex-$gox);
                $svg->rect(x=>$startx,y=>$toy+$tsh/2-$heights{$flag}/2,
                        width=>$endx-$startx,height=>$heights{$flag},style=>"stroke-width:0;fill:$color;");
            }
        }

        $toy += $tsh + $spacing;
    }

    timeLOG("the gene and transcript info was drew ...");

# draw marker 
    if (@marker){
        my $name   = $conf->{marker}->{name};
        my $color  = $conf->{marker}->{color};
        my $height = $conf->{marker}->{height};

        foreach my $mark (@marker) {
            my ($msta,$mend,$mid) = @$mark;
            my $startx = $gox + fetch_axis_dis($coord{$msta},1,$count,$gex-$gox);
            my $endx = $gox + fetch_axis_dis($coord{$mend},1,$count,$gex-$gox);
            $svg->rect(x=>$startx,y=>$toy+$tsh/2-$height/2,
                    width=>$endx-$startx,height=>$height,style=>"stroke-width:0;fill:$color;");
        }

        push @flags , $name;
        $heights{$name} = $height;
        $colors{$name}  = $color;
        $toy += $tsh + $spacing;

        timeLOG("the marker info was drew ...");
    }

# reset the height
    $height = $toy - $tsh/2 - $spacing + $module_margin;
}

#-------------------------------------------------------------------------------
#  draw reads depth distribution and the splice reads info
#-------------------------------------------------------------------------------
# draw connect from top
my $from_top = $conf->{from_top};

# unit height of depth 
my $uhd = $conf->{unit_depth_height};

# unit height of depth spacing
my $uhds = 20;

# ratio of depth height (with $uhd)
my $rdh = 0.8;

foreach my $sample (@samples){
    $height += $uhd;

    my @depths = @{$rdepths{$sample}{depths}};
    my $max_num  = $conf->{fix_axis} && $conf->{max_depth} ? $conf->{max_depth}  :
                                        $conf->{fix_axis}  ? $all_max_depth/$rdh : $rdepths{$sample}{max_depth}/$rdh;
    my $min_num  = $rdepths{$sample}{min_depth} < 0 ? -$max_num : 0;

    my $dividing = dividing($min_num,$max_num,-ntrue=>1,-xtrue=>1);
    my ($min,$max,$step) = split /\s/,$dividing;

    # draw depths distribution
    my $plot_oy = $min_num < 0 ? $height - $uhd/2 : $height;
    if ($styles{$sample}{plot} eq "area"){
        my @xv = ($gox);
        my @yv = ($plot_oy);
        foreach my $depth (@depths){
            my ($pos,$num) = @$depth;
            push @xv , $gox + fetch_axis_dis($pos,1,$count,$gex-$gox);
            my $ydis = $min_num < 0 ? fetch_axis_dis($num,0,$max,$uhd/2) : fetch_axis_dis($num,$min,$max,$uhd);
            push @yv , $plot_oy - $ydis;
        }
        push @xv , $gex;
        push @yv , $height;
        my $points = $svg->get_path(x=>\@xv,y=>\@yv,-type=>'polygon');
        $svg->polygon(%$points,style=>{fill=>$styles{$sample}{color},'stroke-width'=>0},id=>$sample);
    }elsif ($styles{$sample}{plot} eq "bar"){
        foreach my $depth (@depths){
            my ($pos,$num) = @$depth;
            next if ($num == 0);
            my $xv = $gox + fetch_axis_dis($pos,1,$count,$gex-$gox);
            my $ydis = $min_num < 0 ? fetch_axis_dis($num,0,$max,$uhd/2) : fetch_axis_dis($num,$min,$max,$uhd);
            my $yv = $plot_oy - $ydis;
            $svg->line(x1=>$xv,x2=>$xv,y1=>$plot_oy,y2=>$yv,style=>"stroke-width:$styles{$sample}{sr_size};stroke:$styles{$sample}{color}");
        }
    }

    # draw y axis 
    if ($conf->{show_y_axis}){
        $svg->line(x1=>$gox,x2=>$gox,y1=>$height,y2=>$height-$uhd,style=>"stroke-width:1;stroke:#000");
        my $i;
        for ($i=$min;$i<=$max;$i+=$step){
            my $ticky = $height - fetch_axis_dis($i,$min,$max,$uhd);
            $svg->line(x1=>$gox,x2=>$gox+4,y1=>$ticky,y2=>$ticky,style=>"stroke-width:1;stroke:#000");
            my $tick_width = $font->fetch_text_width($i);
            $svg->text(x=>$gox-$spacing-$tick_width,y=>$ticky+$font_height/2,style=>$font_style)->cdata($i)
        }
    }
    
    # show samples 
    my $sample_font = $font;
    $sample_font->setAttr("fill:$styles{$sample}{label_color};weight:bold;font-size:$styles{$sample}{label_size};");
    my $sample_label_width = $sample_font->fetch_text_width($sample);
    if ($conf->{samples_label_pos} eq "left"){
        $svg->text(x=>$gox-$sample_label_width-$max_tick_width-$spacing*2,y=>$height-$uhd/2+$font_height/2,
            style=>$sample_font->toStyle)->cdata($sample);
    }else{
        $svg->text(x=>$gex-$sample_label_width,y=>$height-$uhd+$font_height,style=>$sample_font->toStyle)->cdata($sample);
    }
    
    # draw junctions 
    if ($juncs{$sample}){
        my @juncs = @{$juncs{$sample}};
        my %sites_count;
        my $bottom_flag = 0;
        my $pre_sta = 0;
        my $pre_end = 0;
        foreach my $junc ( sort { $a->[0] <=> $b->[0] || $a->[1] <=> $a->[1] } @juncs){
            my ($jsta,$jend,$jnum,$jid) = @$junc;
            
            next if ($jnum < $conf->{min_junction_depth});
            my $jx1 = $gox + fetch_axis_dis($coord{$jsta},1,$count,$gex-$gox);
            my $jx2 = $gox + fetch_axis_dis($coord{$jend},1,$count,$gex-$gox);
            my $cx1 = $jx1;
            my $cx2 = $jx2;
            my ($jy1,$jy2,$cy);
            $from_top = isoverlap($jsta,$jend,$pre_sta,$pre_end,$from_top,$conf->{from_top});
            if ( $from_top ){
                my $jsta_depth = fetch_site_depth($coord{$jsta},@depths);
                my $jend_depth = fetch_site_depth($coord{$jend},@depths);
                $jy1 = $height - fetch_axis_dis($jsta_depth,$min,$max,$uhd);
                $jy2 = $height - fetch_axis_dis($jend_depth,$min,$max,$uhd);
                $cy  = min($jy1,$jy2) - $ctrl_height;
            } else {
                $jy1 = $height;
                $jy2 = $height;
                $cy  = $height + $ctrl_height;
                $bottom_flag = 1;
            }
            
            # add junction connect line 
            my $pathd = "M $jx1 $jy1 C $cx1 $cy, $cx2 $cy, $jx2 $jy2";
            my $sr_size = sprintf "%.2f" , $styles{$sample}{sr_size} * $jnum / $max_junc;
            $sr_size = 1 if ($sr_size < 1);
            $svg->path(d=>$pathd,style=>"stroke-width:$sr_size;fill:none;stroke:$styles{$sample}{sr_color}",class=>"bezier");
            
            # add the splice reads num [text]
            my $jlabel_w = $sr_font->fetch_text_width($jnum);
            my $jlabel_x = ($jx1 + $jx2)/2 - $jlabel_w/2;
            my $jlabel_y = bezier3_xy([$jx1,$jy1],[$cx1,$cy],[$cx2,$cy],[$jx2,$jy2],0.5) + $sr_font_h/2;
            
            if ($conf->{text_filter}){
                $svg->text(x=>$jlabel_x,y=>$jlabel_y,style=>$sr_font_style,filter=>"url(#textbg)")->cdata($jnum);
            } else {
                my $sr_font_w = $sr_font->fetch_text_width($jnum);
                $svg->rect(x=>$jlabel_x,y=>$jlabel_y-$sr_font_h,width=>$sr_font_w,height=>$sr_font_h,style=>"stroke-width:0;fill:#fff");
                $svg->text(x=>$jlabel_x,y=>$jlabel_y,style=>$sr_font_style)->cdata($jnum);
            }

            $pre_sta = $jsta;
            $pre_end = $jend;
        }

        $height += $ctrl_height if ($bottom_flag);
    }
    
    $height +=  $samples_spacing;
}

timeLOG("the depth and junction reads number was drew ...");

# add the intron_background
if ($conf->{intron_background}){
    for(my$i=1; $i<= $#regions-2; $i+=2){
        my $intron_sta = $regions[$i]+1;
        my $intron_end = $regions[$i+1]-1;
        my $intron_x1 = $gox + fetch_axis_dis($coord{$intron_sta},1,$count,$gex-$gox);
        my $intron_x2 = $gox + fetch_axis_dis($coord{$intron_end},1,$count,$gex-$gox);
        
        $intron_bg->rect(x=>$intron_x1,y=>$oy,width=>$intron_x2-$intron_x1,
            height=>$height-$samples_spacing,style=>"stroke-width:0;fill:#ddd;fill-opacity:0.2");
    }
}

#-------------------------------------------------------------------------------
#  draw the exons in the bottom
#-------------------------------------------------------------------------------
if ($conf->{show_tails} || $OPTS{show_tails}){
    my $bottom = $svg->g(id=>"bottom");

## draw ref genes
    my $goy = $height + $heights{exon}/2;

# draw backbone line
    my $gene_arrowid = "gene_arrow_tail";
    end_arrow($defs,$gene_arrowid,10,6,-zh=>0,-style=>"stroke-width:1;stroke:#000;fill:$colors{exon}");

    my ($gx1,$gx2) = ($gox,$gex);
    ($gx1,$gx2) = ($gx2,$gx1) if ($strand eq "-");
    #$bottom->line(x1=>$gx1,x2=>$gx2,y1=>$goy,y2=>$goy,style=>"stroke-width:2;stroke:$colors{backbone}",
    #        'marker-end'=>"url(#$gene_arrowid)");
    $bottom->line(x1=>$gx1,x2=>$gx2,y1=>$goy,y2=>$goy,style=>"stroke-width:2;stroke:$colors{backbone}");

## draw the backgroud of arrow for line plot
    my $alw = 60;
    my $alh = 10;
    my $alx1 = $alw - $alh * 0.7;
    my $aly1 = $alh/2;

    my $triangle = blank(parent=>$defs,id=>"alb",x=>$gox,y=>$goy-$alh/2,width=>$alw,height=>$alh);
    
    if ($strand eq "-"){
        $triangle->path(d=>"M$alx1 $aly1 L$alw 0 L$alw $alh Z",style=>"stroke-width:1;fill:$colors{backbone}");
    }else{
        $triangle->path(d=>"M$alw $aly1 L$alx1 0 L$alx1 $alh Z",style=>"stroke-width:1;fill:$colors{backbone}");
    }

    $bottom->rect(x=>$gox,y=>$goy-$alh/2,width=>$gex-$gox,height=>$alh,style=>"stroke-width:0;fill:url(#alb)");

# draw exons
    for(my$i=0; $i<= $#exon_regions-1; $i+=2){
        my $exon_sta = $exon_regions[$i];
        my $exon_end = $exon_regions[$i+1];

        my $exon_x1 = $gox + fetch_axis_dis($coord{$exon_sta},1,$count,$gex-$gox);
        my $exon_x2 = $gox + fetch_axis_dis($coord{$exon_end},1,$count,$gex-$gox);

        $bottom->rect(x=>$exon_x1,y=>$height,width=>$exon_x2-$exon_x1,
            height=>$heights{exon},style=>"stroke-width:0;fill:$colors{exon};");
    }

# draw gene label
    my $geneid_width = $font->fetch_text_width($gene);
    my $glx = $gox + ($gex-$gox)/2 - $gene_width/2;
    my $gly = $goy + $heights{exon}/2 + $spacing + $font_height;
    $bottom->text(x=>$glx,y=>$gly,style=>$font_style)->cdata($gene);

# draw chr
    $svg->text(x=>$gox-$chr_width-$spacing,y=>$gly,style=>$font_style)->cdata($chr);
    $svg->text(x=>$gox,y=>$gly,style=>$font_style)->cdata($start);
    $svg->text(x=>$gex-$end_width,y=>$gly,style=>$font_style)->cdata($end);

    $height += $heights{exon} + $spacing + $font_height;
}

# draw legend for TS
my $ly = $oy;
my $lx = $gex + $spacing * 2;
my $lh = $tsh;
if ($legend && $conf->{show_transcirpts})
{
    foreach my $item (@flags)
    {
        my $item_h = $heights{$item};
        $svg->rect(x=>$lx,y=>$ly+$lh/2-$item_h/2,width=>$lw,height=>$item_h,style=>"stroke-width:0;fill:$colors{$item};");
        $svg->text(x=>$lx+$lw+$spacing,y=>$ly+$lh/2+$font_height/2,style=>$font_style)->cdata($item);
        $ly += $lh + $spacing;
    }

    $ly += $lh + $spacing;
    timeLOG("the Legend was drew ...");
}

# draw unit scale 
if (defined $conf->{scale}){

    my $scale = $conf->{scale} == 0 ? auto_scale($legend_width*0.5,$count,$gex-$gox) : $conf->{scale};
    my $size  = fetch_axis_dis($scale,1,$count,$gex-$gox);
    my $scale_y = $ly + $font_height + $spacing;
    my $label_w = $font->fetch_text_width($scale);
    
    $svg->text(x=>$lx+$size/2-$label_w/2,y=>$ly + $font_height,style=>$font_style)->cdata($scale);
    $svg->line(x1=>$lx,x2=>$lx+$size,y1=>$scale_y,y2=>$scale_y,style=>"stroke-width:2;stroke:#000000");
    $svg->line(x1=>$lx,x2=>$lx,y1=>$scale_y-4,y2=>$scale_y+4,style=>"stroke-width:2;stroke:#000000");
    $svg->line(x1=>$lx+$size,x2=>$lx+$size,y1=>$scale_y-4,y2=>$scale_y+4,style=>"stroke-width:2;stroke:#000000");
    
    timeLOG("the Scale was drew ...");
}

##  save figure
# reset the svg heigt
my $root = $svg->getElementByID("as_map");
$root->setAttribute("height",$height+$margin);

my $outfile = $conf->{'methy_type'} ? "$gene.$conf->{methy_type}.CovPlot.svg" : "$gene.CovPlot.svg";

open OUT,">$outfile" or die $!;
print OUT $svg->xmlify;
close OUT;
timeLOG("the figure [$outfile] was finished :)");

#-------------------------------------------------------------------------------
#  sub functions
#-------------------------------------------------------------------------------
sub fetch_gene_loci {
    my ($gtf,$gene) = @_;

    my %flags = map { $_ => 1 } split /,/ , $conf->{element_flags};
    my @loci;
    my %loci;
    
    system("grep $gene $gtf > $gene.tmp.gtf");
    $gtf = "$gene.tmp.gtf";

    open GTF,$gtf or die "can't open gtf file, $gtf $!";
    while(<GTF>){
        next if -1 == index($_,"$gene\"");
        chomp;
        
        my ($chr,$flag,$start,$end,$strand,$attrs) = (split /\t/)[0,2,3,4,6,8];
        next unless $flags{$flag};
        
        my $tsid     = $attrs =~ /transcript_id \"([\w\-\.]+)\";/  ? $1 : die;
        my $geneid   = $attrs =~ /gene_id \"([\w\-\.]+)\";/       ? $1 : die;
        my $genename = $attrs =~ /gene_name \"([\w\-\.]+)\";/     ? $1 : "";
        next unless ($geneid eq $gene || ($genename && $genename eq $gene) );

        push @{$loci{loci}{$tsid}{$flag}} , [$start,$end];
        $loci{chr} = $chr unless $loci{chr};
        $loci{strand} = $strand unless $loci{strand};

        push @loci , $start;
        push @loci , $end;
    }
    close GTF;
    system("rm $gene.tmp.gtf");
    
    ERROR("gene_not_in_gtf, $gene is not exists in gtf file!") if ($#loci == -1);

    $loci{start} = min(@loci);
    $loci{end} = max(@loci);
    $loci{gene_start} = $loci{start};
    $loci{gene_end} = $loci{end};
    
    if ($loci{strand} eq "+"){
        $loci{start} -= $conf->{upstream}   if $conf->{upstream};
        $loci{end}   += $conf->{downstream} if $conf->{downstream};
    }elsif ($loci{strand} eq "-"){
        $loci{start} -= $conf->{downstream} if $conf->{downstream};
        $loci{end}   += $conf->{upstream}   if $conf->{upstream};
    }
   
    return %loci;
}

# check the conf and set some default values 
sub check_conf {
    my $conf = shift;
    
    default_OPTS($conf,"gtf");
    default_OPTS($conf,"loci");
    default_OPTS($conf,"gene");
    default_OPTS($conf,"intron_fix");
    default_OPTS($conf,"intron_scale");
    default_OPTS($conf,"width");
    default_OPTS($conf,"unit_depth_height");
    default_OPTS($conf,"intron_background");
    default_OPTS($conf,"upstream");
    default_OPTS($conf,"downstream");
    default_OPTS($conf,"methy_type");
    $conf->{legend} = 1 if ($OPTS{legend});

    ERROR("[no_gtf_ERROR] gtf file must be defined") unless $conf->{gtf};
    unless ($conf->{gene}){
        WARN("[no_gene_ERROR] the target gene must be defined\n");
        &usage;
    }
    
    default_set($conf,"width",1000);
    default_set($conf,"unit_depth_height",100);
    default_set($conf,"intron_background",1);
    default_set($conf,"intron_scale","auto");
    default_set($conf,"min_junction_depth",3);
    default_set($conf,"format","tophat2");
    default_set($conf,"samples_label_pos","left");
    default_set($conf,"upstream",0);
    default_set($conf,"downstream",0);
    default_set($conf,"element_flags","exon");
    default_set($conf,"show_y_axis",1);

    return $conf;
}

# the OPTS defined will cover the conf file defined
sub default_OPTS {
    my ($conf,$attr) = @_;
    $conf->{$attr} = $OPTS{$attr} if (defined $OPTS{$attr});
}

# set the default value which not defined 
sub default_set {
    my ($conf,$attr,$val) = @_;
    $conf->{$attr} = $val unless (defined $conf->{$attr} && $conf->{$attr} ne "");
}

sub inherit_opts {
    my $parent = shift;
    my $child = shift;

    foreach my $attr (keys %$parent){
        $child->{$attr} = $parent->{$attr} unless exists $child->{$attr};
    }
}

# init_styles
sub init_styles {
    my ($conf,$samples) = @_;
    
    my %styles = ();

    foreach my $sample (sort {$samples->{$a}->{order} <=> $samples->{$b}->{order}} keys %$samples){
        my $style = {};
        inherit_opts($conf->{ases},$style);
        default_set($style,"label_size",12);
        default_set($style,"sr_size",3);
        default_set($style,"label_color",$samples->{$sample}->{color});
        default_set($style,"sr_color",$samples->{$sample}->{color});
        default_set($style,"color",$samples->{$sample}->{color});
        default_set($style,"plot",$samples->{$sample}->{plot});
        $styles{$sample} = $style;
    }
    
    return %styles;
}

# parse styles 
sub parse_styles {
    my $conf = shift ;
    my $samples = shift;
    my %styles = ();
    
    my @ases = ref $conf->{ases}->{as} eq "HASH"  ? ( $conf->{ases}->{as}  ) : 
               ref $conf->{ases}->{as} eq "ARRAY" ? @{ $conf->{ases}->{as} } : die;
   
    foreach my $as (@ases){
        inherit_opts($conf->{ases},$as);
        default_set($as,"label_color",$as->{color});
        default_set($as,"label_size",12);
        default_set($as,"sr_color",$as->{color});
        default_set($as,"sr_size",3);
        
        my $label = $as->{label} ? $as->{label} : (split /\./ ,  basename($as->{bam}) )[0];
        $styles{$label} = $as;
    }

    return %styles;
}

# fetch the styles 
sub fetch_styles {
    my $conf = shift;
    my %hash = @_;

    return %hash unless $conf;

    foreach my $key (keys %$conf){
        $hash{$key} = $conf->{$key};
    }

    return %hash;
}

sub exon_or_intron {
    my ($loci,$flags,$conf) = @_;
    my ($start,$end) = ($loci->{start},$loci->{end});

    my %coord = map { $_ => 0 } $start .. $end;
    foreach my $tsid (keys %{$loci->{loci}}){
        foreach my $flag (@$flags){
            next unless $loci->{loci}{$tsid}{$flag};
            my @recoreds = @{$loci->{loci}{$tsid}{$flag}};
            foreach my $record (sort {$a->[0] <=> $b->[0]} @recoreds){
                map { $coord{$_} = 1 } $record->[0] .. $record->[1];
            }
        }
    }
    
    if ( ($conf->{upstream} && $loci->{strand} eq "+") || ($conf->{downstream} && $loci->{strand} eq "-") ){
        for ($start .. $loci->{gene_start}-1){ $coord{$_} = 1 }
    }

    if ( ($conf->{downstream} && $loci->{strand} eq "+") || ($conf->{upstream} && $loci->{strand} eq "-" ) ){
        for ($loci->{gene_end}+1 .. $end){ $coord{$_} =1 }
    }

    return %coord;
}

#-------------------------------------------------------------------------------
# rebuild the gene coordinate
# rebuild the CDS and UTR coord
# intron_fix, set all introns as length [intron_fix]
# intron_scale, zoom the intron size as [intron_scale] ratio to raw size
# intron_fix is prior to intron_scale
#-------------------------------------------------------------------------------
sub rebuild_gene_coord {
    my ($loci,$flags,$conf) = @_;
    my ($start,$end) = ($loci->{start},$loci->{end});

    my %coord = map { $_ => 0 } $start .. $end;
    foreach my $tsid (keys %{$loci->{loci}}){
        foreach my $flag (@$flags){
            next unless $loci->{loci}{$tsid}{$flag};
            my @recoreds = @{$loci->{loci}{$tsid}{$flag}};
            foreach my $record (sort {$a->[0] <=> $b->[0]} @recoreds){
                map { $coord{$_} = 1 } $record->[0] .. $record->[1];
            }
        }
    }
   
    if ($conf->{upstream}){
        for ($start .. $loci->{gene_start}-1){ $coord{$_} = 1 }
    }

    if ($conf->{downstream}){
        for ($loci->{gene_end}+1 .. $end){ $coord{$_} =1 }
    }

    my $count = 0;
    if ($conf->{intron_fix}){
        my $flag = 1;
        my $intron_scale = 1;
        
        foreach my $i ( $start .. $end ){
            if (1 == $coord{$i}){
                $count ++;
                $flag = 1;
            } elsif (0 == $coord{$i} && 1 == $flag){
                $intron_scale = fetch_intron_scale(\%coord,$i,$conf->{intron_fix});
                $count += $intron_scale;
                $flag = 0;
            } elsif (0 == $coord{$i} && 0 == $flag){
                $count += $intron_scale;
            }
            $coord{$i} = $count;
        }

    } elsif ($conf->{intron_scale} ne "auto") {
        map {
                $count = $coord{$_} == 1 ? $count + 1 : $count + $conf->{intron_scale};
                $coord{$_} = $count;
            } $start .. $end;
    } elsif ( $conf->{intron_scale} eq "auto" ) {
        my $exon_len = 0;
        map { $exon_len ++ if $coord{$_} == 1 } $start .. $end;
        $conf->{intron_scale} = $end-$start+1-$exon_len == 0 ? 1 : $exon_len/($end-$start+1-$exon_len);
        timeLOG("  the 'intron_scale' will be set as [$conf->{intron_scale}] auto ...");
        map {
                $count = $coord{$_} == 1 ? $count + 1 : $count + $conf->{intron_scale};
                $coord{$_} = $count;
            } $start .. $end;
        
    }

    return %coord;
}

sub fetch_intron_scale {
    my ($coord,$i,$intron_fix) = @_;
    my $intron_len = 0;
    while(exists $coord->{$i}){
        if (0 == $coord->{$i}){
            $intron_len ++;
            $i++
        } else {
            last;
        }
    }

    return  $intron_fix/$intron_len;
}

sub end_arrow
{
    my ($defs,$markid,$w,$h,%par) = @_;
    my $zh = defined $par{'-zh'} ? $par{'-zh'} : 1;
    $zh = 1 if ($zh < 0 || $zh > 1);
    my $zhx = $zh * $w;
    my $halfh = $h/2;
    my $marker = $defs->marker(
            refX => $w, refY => $h/2, # the basic point coord, is very important
            markerWidth=>$w,markerHeight=>$h,orient => 'auto',
            markerUnits => 'strokeWidth',id=>$markid);
    $marker->setAttribute('style',$par{-style}) if ($par{-style});
    $marker->setAttribute('class',$par{-class}) if ($par{-class});

    $marker->path(d => "M0,0 L$w,$halfh L0,$h L$zhx,$halfh z");
}

sub fetch_depths_from_conf {
    my ($conf,$chr,$sta,$end) = @_;
    my (%depths,%samples);
    
    my @ases = ref $conf->{ases}->{as} eq "HASH"  ? ( $conf->{ases}->{as}  ) : 
               ref $conf->{ases}->{as} eq "ARRAY" ? @{ $conf->{ases}->{as} } : die;
    
    my $order = 0;
    my @colors = $OPTS{colors} ? split /,/ , $OPTS{colors} : split /,/ , $conf->{colors}->{samples};
    @colors = map { "#000000" } 0 .. $#ases if ($#ases > $#colors);
    foreach my $as (@ases){
        my $bam_file = $as->{bam};
        my $label = $as->{label} ? $as->{label} : (split /\./ ,  basename($bam_file) )[0];
        $depths{$label} = fetch_depth($bam_file,$chr,$sta,$end);
        $samples{$label}{order} = $order;
        $samples{$label}{bam}   = $bam_file;
        $samples{$label}{juncs} = $as->{junction} if ($as->{junction});
        $samples{$label}{color} = $as->{color} // $colors[$order];

        $samples{$label}{plot}  = $as->{plot}            ? $as->{plot} : 
                                  $bam_file =~ /\.bam$/  ? "area"      : 
                                  $bam_file =~ /\.cout$/ ? "bar"       : 
                                  $bam_file =~ /\.bed/   ? "area"      : "area";

        $order ++;
    }

    return (\%depths,\%samples);
}

sub fetch_depths_from_ARGV {
    my ($bam_ref,$chr,$sta,$end) = @_;

    my (%depths,%samples);
    my @colors = $OPTS{colors} ? split /,/ , $OPTS{colors} : split /,/ , $conf->{colors}->{samples};
    @colors = map { "#000000" } 0 .. $#$bam_ref if ($#$bam_ref > $#colors);
    my @labels = split /,/ , $OPTS{labels} if ($OPTS{labels});
    my @juncs  = split /,/ , $OPTS{juncs}  if ($OPTS{juncs});

    my $order = 0;
    foreach my $bam_file (@$bam_ref){
        my $label = $OPTS{labels} ? $labels[$order] : (split /\./ ,  basename($bam_file) )[0];
        $depths{$label} = fetch_depth($bam_file,$chr,$sta,$end);
        $samples{$label}{order} = $order;
        $samples{$label}{bam}   = $bam_file;
        $samples{$label}{juncs} = $juncs[$order]  if ($OPTS{juncs});
        $samples{$label}{color} = $colors[$order];
        
        $samples{$label}{plot}  = $bam_file =~ /\.bam$/  ? "area"      : 
                                  $bam_file =~ /\.cout$/ ? "bar"       : 
                                  $bam_file =~ /\.bed/   ? "area"      : "area";

        $order ++;
    }

    return (\%depths,\%samples);
}

sub fetch_depth {
    my ($bam,$chr,$sta,$end) = @_;
    my %depth;
    
    ERROR("bam file is not exists, [$bam]") unless -e $bam;

    if ($bam =~ /\.bam$/){
        ERROR("cannot found the index for [$bam]") unless -e "$bam.bai";

        my $text = `$samtools depth -d 0 -r "$chr:$sta-$end" $bam`;
        my @lines = split /\n/ , $text;

        foreach my $line (@lines){
            my ($scf,$site,$num) = split /\t/ , $line;
            $depth{$site} = $num;
        }
    }elsif ($bam =~ /\.cout$/){
        open my $fh,$bam or die $!;
        while(<$fh>){
            chomp;
            my ($scf,$site,$strand,$type,$mC,$nmC) = (split /\t/)[0,1,2,3,6,7];
            next unless ($scf eq $chr && $site>=$sta && $site <=$end);
            next if ($conf->{methy_type} && $conf->{methy_type} ne $type);
            my $sum = $nmC + $mC;
            next if $sum == 0;
            my $ratio = $mC / $sum;
            $ratio = -$ratio if ($strand eq "-" && $OPTS{strand});
            $depth{$site} = $ratio;
        }
        close $fh;
    }elsif ($bam =~ /\.bed$/){
        
    }
    
    my @depths = map { $depth{$_} || 0 } $sta .. $end;
    return \@depths;
}

sub calc_depths {
    my ($depths,$window,$regions,$coord) = @_;
    my @region_depths = ();

    my $start = $regions->[0];

    # the first point
    my $index = 0;
    my $pos   = $index;
    my $depth = 0;
    
    for (my$i=0; $i<= $#regions-1; $i+=2){
        my $exon_sta = $regions->[$i] - $start;
        my $exon_end = $regions->[$i+1] - $start;
        
        # the exon start site
        push @region_depths , [ $exon_sta + $start , $depths->[$exon_sta] ];
        
        # the exon region with window
        my $j;
        for ($j=$exon_sta;$j<$exon_end-$window;$j+=$window){
            $depth = sum( @{$depths}[$j+1 .. $j+$window] )/$window;
            push @region_depths , [ $j + $window + $start , $depth ];
        }
        
        # besides exon region (the exon stop site)
        if ($j+$window > $exon_end){
            my @nums = @{$depths}[ $j+1 .. $exon_end ];
            $depth = sum (@nums)/($#nums+1);
            push @region_depths , [ $exon_end + $start , $depth ];
        }
        
        last if ($i == $#regions - 1); # the last exon

        ##  intron region 
        my $intron_sta = $exon_end + 1;
        my $intron_end = $regions->[$i+2] - 1 - $start;
        
        # the intron start site 
        push @region_depths , [ $intron_sta + $start , $depths->[$intron_sta] ];
        
        # the intron region with window
        for ($j=$intron_sta;$j<$intron_end-$window;$j+=$window){
            $depth = sum( @{$depths}[$j+1 .. $j+$window] )/$window;
            push @region_depths , [ $j + $window + $start , $depth ];
        }           

        # besides intron region (the intron stop site)
        if ($j+$window > $intron_end){
            my @nums = @{$depths}[ $j+1 .. $intron_end ];
            $depth = sum (@nums)/($#nums+1);
            push @region_depths , [ $intron_end + $start , $depth ];
        }
    }
    
    @region_depths = map { [ $coord->{$_->[0]} , $_->[1] ] } @region_depths;
    return @region_depths;
}

sub fetch_junctions {
    my ($conf,$samples,$chr,$sta,$end,$aimgene) = @_;
    my %junctions;
    
    foreach my $sample (keys %{$samples}){
        my $junc_file = $samples->{$sample}{juncs} || next;
        $junctions{$sample} = $conf->{format} eq "tophat2" ? fetch_junction($junc_file,$chr,$sta,$end) : fetch_bed($junc_file,$aimgene);
    }

    return %junctions;
}

sub fetch_junction {
    my ($junc_file,$aim_chr,$aim_sta,$aim_end) = @_;
    my @juncs;
    
    open my $fh_junc , $junc_file or die $!;
    <$fh_junc>;
    while(<$fh_junc>){
        chomp;
        my @arr = split /\t/;
        my ($chr,$sta,$end) = @arr[0,1,2];
        if ( $chr eq $aim_chr && $sta >= $aim_sta && $end <= $aim_end )
        {
            my ($junc_id,$num,$junc_len,$junc_sta) = @arr[3,4,10,11];
            my @junc_lens = split /,/ , $junc_len;
            my @junc_stas = split /,/ , $junc_sta;
            
            for my $i (1 .. $#junc_lens){
                my $jsta = $sta + $junc_stas[$i-1] + $junc_lens[$i-1];
                my $jend = $sta + $junc_stas[$i] + 1;
                push @juncs , [$jsta,$jend,$num,$junc_id];
            }
        }
    }
    close $fh_junc;

    return \@juncs;
}

sub fetch_bed {
    my ($file,$gene) = @_;
    my @juncs;

    open my $fh , $file or die $!;
    while(<$fh>){
        chomp;
        my ($chr,$sta,$end,$id,$geneid,$num) = split /\t/;
        next unless $geneid eq $gene;

        push @juncs , [$sta,$end,$num,$id];
    }
    close $fh;

    return \@juncs;
}

# fetch marker region
sub fetch_marker {
    my ($conf,$chr,$sta,$end) = @_;
    my @marker;
    
    return @marker unless ($conf->{marker}  && $conf->{marker}->{loci});

    open my $fh_marker , $conf->{marker}->{loci} or die $!;
    while(<$fh_marker>){
        chomp;
        next if (/^#/);
        my ($mchr,$msta,$mend,$id) = split /\t/;
        if ($mchr eq $chr && $msta >= $sta && $mend <= $end){
            push @marker , [$msta,$mend,$id];
        }
    }
    close $fh_marker;

    return @marker;
}

# fetch the regions with exon_or_intron flags hash
sub fetch_regions {
    my ($hash,$sta,$end) = @_;
    my @exons;
    
    my @regions = ($sta);
    for my $i ( $sta+1 .. $end ){
        if ($hash->{$i} != $hash->{$i-1}){
            push @regions , $i;
        }
    }

    push @regions , $end;
    return @regions;
}

sub fetch_axis_dis{
    my ($val,$min,$max,$axis_len) = @_;
    my $dis = sprintf ( "%.4f" , (($val-$min)/($max-$min))*$axis_len );
    return $dis;
}

sub fetch_max_junc {
    my %hash = @_;
    my @maxes = map {
        my @nums = map { $_->[2] } @$_;
        @nums ? max(@nums) : 0;
        } values %hash;
    return max(@maxes);
}

sub fetch_site_depth {
    my ($site,@depths) = @_;
    
    my $res;
    foreach (@depths) {
        my ($pos,$num) = @$_;
        if ($site == $pos){
            $res = $num;
            last;
        } elsif ($site > $pos) {
            $res = $num;
            next;
        } else {
            last;
        }
    }
    return $res;
}

sub isoverlap {
    my ($s1,$e1,$s2,$e2,$flag,$rawflag) = @_;

    my $isoverlap = ($s1 >= $s2 && $s1 <= $e2) || ($e1 >= $s2 && $e1 <= $e2) ? 1 : 0;
    if ($isoverlap) { $flag = not $flag; return $flag } else { return $rawflag }
}

sub auto_scale {
    my ( $size , $count , $total_size ) = @_;
    my $len = sprintf ("%e" , $size * $count / $total_size);
    my ($num,$unit) = split /e/ , $len;
    return (round($num) * (10**$unit));
}

sub usage {
    print <<HELP;
Usage:   

set bams in conf: perl $0 [options] --conf your.conf [ --gene target_gene ]
set bams in cmd:  perl $0 [options] --conf your.conf [ --gene target_gene ] <*.bam|cout|bed>[s]

Import Options:
         --help                      print the simple usage info
         --conf   <FILE>             the configuration file, demo: "$Bin/covPlot.conf"
         --gene   <STR>              set the gene name which you want to draw (must be in the gtf file), conlict with --loci
         --gtf    <FILE>             set the gtf file

Options for simple model ( bam files defined in cmd ):
         --juncs  <FILE>[s]          set the junction bed files, Separated by commas
         --labels <STR>[s]           set the labels for bam files, Separated by commas
         --colors <STR>[s]           set the colors for bam files, Separated by commas

Styles Options (it's better to be set in conf file):
         --intron_fix         <INT>      fix the intron sizse as [200] nt, recommended if you set upstream or downstream
         --intron_scale       <FLOAT>    zoom the intron size as [auto] ratio to raw size, [scale the introns total size == exons]
         --intron_background             display the intron backgroud or not
         --width              <INT>      the width size of the figure, [1000]
         --from_top                      the connect junction lines draw from top first, otherwise will from bottom
         --unit_depth_height  <INT>      the unit depth height of each sample, [100]
         --show_heads                    show the head part ( plot the gene and each transcript in the head part )
         --show_tails                    show the tail part ( plot the gene and exons in the tail part )

Note: for --intron_scale, the upstream and downstream are considered as exons

HELP
    exit 1;
}
