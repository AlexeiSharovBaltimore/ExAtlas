#!/usr/bin/perl
use strict;
use POSIX;
use File::Copy;

#*******************************************************
# Copyright @2014, The National Institute on Aging (NIA/NIH).
# All rights reserved.
# 
# This software is provided "AS IS".  NIA/NIH makes no warranties, express
# or implied, including no representation or warranty with respect to
# the performance of the software and derivatives or their safety,
# effectiveness, or commercial viability.  NIA does not warrant the
# merchantability or fitness of the software and derivatives for any
# particular purpose, or that they may be exploited without infringing
# the copyrights, patent rights or property rights of others. NIA shall
# not be liable for any claim, demand or action for any loss, harm,
# illness or other damage or injury arising from access to or use of the
# software or associated information, including without limitation any
# direct, indirect, incidental, exemplary, special or consequential
# damages.
# 
# This software program may not be sold, leased, transferred, exported
# or otherwise disclaimed to anyone, in whole or in part, without the
# prior written consent of NIA.
#
# Programmer: Alexei Sharov (sharoval@mail.nih.gov)
# National Institute on Aging, Genetics Lab
#*******************************************************

my %hashIni;
read_configuration("../../exatlas.ini",\%hashIni);
my $PATH_HOME=$hashIni{"PATH_HOME"};
my $PATH_PROG=$hashIni{"PATH_PROG"};
my $CGI_ADDRESS=$hashIni{"CGI_ADDRESS"};
my $HOME_ADDRESS=$hashIni{"HOME_ADDRESS"};
my $SECURITY=$hashIni{"SECURITY"};
my $UNIX=$hashIni{"UNIX"};
my $PATH_INFO = "$PATH_PROG/info";
my $PATH_DATA = "$PATH_PROG/data";
my $PATH_OUTPUT = "$PATH_HOME/output";
my $PATH_BIN = "$PATH_PROG/bin";
my $LINE = 2;
my $TEXT = 10;
my $BOX = 14;
my $CIRCLE = 1;
my $FILL = 8;
my $radius = 2;
my $black = 15;
my $gray = 14;
my $blue = 8;
my $magenta = 6;
my $red = 21;
my $green = 16;
my $ltgreen = 31;
my $orange = 17;
my $tinyFont = 0;
my $smallFont = 1;
my $largeFont = 2;
my $smallFontInbox = 4;
my $largeFontInbox = 5;
my $vertTinyFont = 6;
my $vertSmallFont = 7;
my $vertLargeFont = 8;
my $printNumbers = 1;
my $thin = 1;
my $thick2 = 2;
my $thick3 = 3;
my $solid = 0;
my $MISSING = -9999;
my $MAX_GENES = 3000;
my $EPSYLON = 1.0e-15;
my $ITMAX =200;
my $EPS =3.0e-7;
my $FPMIN =1.0e-30;
my $maxHeaderLength =   50;
my ($RUN_QUALITY, $RUN_OVERLAP, $RUN_SIGNIFICANT, $RUN_MATRIX_UPLOAD, $RUN_PAGE, 
	$RUN_CORRELATION, $RUN_ANOVA_FIRST, $RUN_ANOVA_PARAM, $RUN_ANOVA_NORM, $RUN_NORMALIZE,
	$RUN_METAANALYSIS, $RUN_PCA_PCA, $RUN_PCA_CLUSTER) = 1..13;
my $READ_ANOVA_ALL_ROW = 1;
my $READ_ANOVA_ALL_COL = 2;
my $READ_ANOVA_ONE_ROW = 3;
my $READ_ANOVA_ONE_COL = 4;
my $READ_ANOVA_NONE    = 5;
my $READ_BY_NUMBER     = 1;
my $READ_BY_NAME       = 2;
my $READ_GENESET_ALL   = 1;
my $READ_GENESET_ONE   = 2;
my $READ_GENESET_NONE  = 3;
my $READ_NAMES         = 1;
my $READ_NAMES_DESCRIP = 2;
my @correlation_list = ("Pearson","Spearman","Covariance");
my @colors=($blue,$magenta,$orange,$green,$ltgreen);
my @method_name=("Randrom effect","Fixed effect","Fisher's method","Z-score method","All methods");
my @unacceptable = ("asshole","balls","bastard","bitch","bloody","bollocks","bugger","bullshit","christ","cock","cocksucker","fuck","jesus","nigger","retard","shit","slut","whore");
my @dirNames = ("_up","_down");

my @organisms = ([9606,9000,"Homo sapiens","Human"],
[10090,9000,"Mus musculus","Mouse"],
[9544,9000,"Macaca mulatta","Rhesus monkey"],
[9541,9000,"Macaca fascicularis","Macaque"],
[9598,9000,"Pan troglodytes","Chimpanzee"],
[9600,9000,"Pongo pygmaeus","Orangutan"],
[9483,9000,"Callithrix jacchus","Marmoset","primate"],
[10116,9000,"Rattus norvegicus","Rat"],
[9615,9000,"Canis lupus familiaris","Dog"],
[9940,9000,"Ovis aries","Sheep"],
[9823,9000,"Sus scrofa","Pig"],
[9913,9000,"Bos taurus","Cow"],
[9796,9000,"Equus caballus","Horse"],
[9986,9000,"Oryctolagus cuniculus","Rabbit"],
[9031,9000,"Gallus gallus","Chicken"],
[9103,8000,"Meleagris gallopavo","Turkey"],
[8364,8000,"Xenopus tropicalis","Xenopus frog"],	
[8355,8000,"Xenopus laevis","Xenopus frog"],
[7955,8000,"Danio rerio","Zebrafish"],
[8022,8000,"Oncorhynchus mykiss","Rainbow trout"],
[8030,8000,"Salmo salar","Salmon"],
[7719,7000,"Ciona intestinalis","Acsidia"],
[7668,7000,"Strongylocentrotus purpuratus","Sea urchin"],
[7227,8000,"Drosophila melanogaster","Fruit fly"],
[7460,8000,"Apis mellifera","Honey bee"],
[7159,8000,"Aedes aegypti","Mosquito"],
[7165,8000,"Anopheles gambiae","Mosquito malaria"],
[7091,8000,"Bombyx mori","Silkworm"],
[6239,7000,"Caenorhabditis elegans","Nematode"],
[3702,7000,"Arabidopsis thaliana","Thale cress"],
[4530,7000,"Oryza sativa","Rice"],
[3847,7000,"Glycine max","Soybean"],
[4565,7000,"Triticum aestivum","Wheat"],
[4113,7000,"Solanum tuberosum","Potato"],
[4081,7000,"Solanum lycopersicum","Tomato"],
[4577,7000,"Zea mays","Maize"],
[3635,7000,"Gossypium hirsutum","Cotton"],
[927,4000,"Saccharomyces pombe","Yeast"],
[4932,4000,"Saccharomyces cerevisiae","Yeast"],
[44689,7000,"Dictyostelium discoideum","Slime mold"],
[562,4000,"Escherichia coli","Bacterium"],
[300852,2000,"Thermus thermophilus","Bacterium"],
[1280,2000,"Staphylococcus aureus", "Bacterium"]);

my %hashOrganism;
my %hashOrganismNum;
my %hashSpeciesNum;
for(my $i=0; $i<@organisms; $i++){
	my $ref = $organisms[$i];
	$hashOrganism{$ref->[0]}="$ref->[3] ($ref->[2])";
	$hashOrganismNum{$ref->[0]} = $i+1;
	$hashSpeciesNum{$ref->[2]} = $i+1;
}
# Read data from the input form
my %hashInput;
my $request_method = $ENV{'REQUEST_METHOD'};
my $size = $ENV{'CONTENT_LENGTH'};
my $remoteID = $ENV{'REMOTE_ADDR'};
my $form_info;
if($size){
	read (STDIN, $form_info, $size);
}else{
	$form_info = $ENV{'QUERY_STRING'};
}
if(!$form_info){
	$form_info = $ARGV[0];
}
my $loginname;
my $passwd;
my $sessionID;

if(!$form_info){ print "content-type: text/html","\n\n"; error_message("No data","continue"); }
my $date = `date \'+%y%m%d\'`;
chop($date);
my $date_record = `date \'+%Y-%m-%d\'`;
chop($date_record);
if($form_info =~ /\n/){
	print "content-type: text/html","\n\n";
	parse_data(\$form_info,\%hashInput);
	$sessionID = $hashInput{"sessionID"};
	foreach my $key (keys %hashInput){
		my $value = $hashInput{$key};
		if(check_form($key,$value)){
			$loginname = check_sessionID($sessionID);
			if(ref($value) eq 'ARRAY'){ $value = "$value->[0]"; }
			my $email1 = "sharoval\@mail.nih.gov";
			`echo "$loginname: $key - $value" | mailx -s "ExAtlas failed" $email1`;
			error_message("Operation failed!","continue");
		}
	}
	$loginname = check_sessionID($sessionID);
	my $changeOrganism = $hashInput{"changeOrganism"};
	if($changeOrganism && $hashInput{"action"} eq "file_upload"){
		file_upload();
	}
	my $rename = $hashInput{"rename"};
	my ($filename,$list);
	my $ref = $hashInput{"upload_filename"};
	if($ref){
		($filename,$list) = @$ref;
	}
	if($rename){
		$filename = $rename;
	}
	my $list1 = $hashInput{"pasted_text"};
	if($list1){
		$list = $list1;
	}
	my $file_type = $hashInput{"file_type"};
	my $action = $hashInput{"action"};
	if($action !~ /_expression_profile/){
		if(!$list){ error_message("No data found!","continue"); }
		if(!$filename && $file_type !~/^(genelist|expression)$/){ error_message("Filename is empty!","continue"); }
		if(!$file_type){ error_message("No file type!","continue"); }
		if($file_type ne "genelist" && (length($filename) > 150 || $filename=~/[\s\\\/\+\@<>\!\*\|]/)){
			error_message("Invalid filename (too long or includes special characters)","continue");
		}
	}
	my $organismID = $hashInput{"organismID"};
	if(!$organismID){ error_message("No organism ID!","continue"); }
	my @lines = split(/\n/,$list);
	my $description = $hashInput{"description"};
	my $file_info;
	my $platform;
	if($file_type eq "expression" || $action =~ /_expression_profile/){
		save_expression_profile(\@lines,$organismID,$filename,$description);
	}
	if($file_type eq "genelist"){
		my $genelist_fileID = save_list_of_genes($list,$organismID,$filename,$description);
		$hashInput{"genelist_fileID"} = $genelist_fileID;
	}
	if($file_type eq "legend"){
		save_legend(\@lines,$filename);
		terminal_window("<H3>Legend '$filename' is uploaded</H3>","continue");
	}
	elsif($file_type eq "matrix"){
		check_matrix_data(\@lines,$filename);
	}
	elsif($file_type eq "geneset"){
		$file_info = check_geneset_data(\@lines,$organismID,$filename,$description);
	}
	elsif($file_type eq "samples"){
		$file_info = check_samples_data(\@lines,$organismID,$filename,$description);
	}
	elsif($file_type eq "output"){
		$file_info = check_output_data(\@lines,$organismID,$filename,$description);
	}
	elsif($file_type eq "annotation"){
		$file_info = check_annotation_data(\@lines,$organismID,$filename,$description);
		$platform = $filename;
		$platform =~ s/_annot.txt$//;
		if($filename !~ /_annot.txt$/){ $filename .= "_annot.txt"; }
	}
	if($file_type !~ /^genelist$|^expression$/){
		finish_file_upload(\@lines,$filename,$file_type,$organismID,$platform);
		terminal_window("<H3>File '$filename' has been uploaded</H3>$file_info","continue");
	}
}else{
	print "content-type: text/html","\n\n";
	my @items = split(/[\&=]/,$form_info);
	while(@items){
		my $key = shift(@items);
		my $value = shift(@items);
		$key = substitute_chars($key);
		$value = substitute_chars($value);
		if($key eq ""){ error_message("Empty key!"); }
		my $x = length($value);
		if($key=~/^(new_)*passwd1*$/ && ($x<4 || $x>20) && $hashInput{"loginname"} !~ /^guest$/i){ error_message("Password should be from 5 to 20 characters!","register"); }
		$hashInput{$key} = $value;
		if(check_form($key,$value)){
			my $email1 = "sharoval\@mail.nih.gov";
			`echo "$key - $value" | mailx -s "ExAtlas failed" $email1`;
			error_message("Operation failed!");
		}
	}
	$loginname = $hashInput{"loginname"};
	$passwd = $hashInput{"passwd"};
	if($hashInput{"register"}){
		register_new_user();
	}elsif($hashInput{"login"}){
		if($loginname =~ /^guest$/i){
			my $xxx = int(rand 1000);
			$loginname = "guest$date$xxx";
			$hashInput{"loginname"} = $loginname;
			copy "$PATH_DATA/default-config.txt", "$PATH_INFO/$loginname-config.txt";
		}
		validate($loginname,$passwd);
		$sessionID = start_session($loginname);
	}elsif($hashInput{"reset_password"}){
		reset_password($loginname);
	}else{
		$sessionID = $hashInput{"sessionID"};
		$loginname = check_sessionID($sessionID);
	}
}
#print "FORM: $form_info\n";
#print %hashInput;
#print "$loginname - $action<p>\n";

$hashInput{"loginname"} = $loginname;
clean_up();
my $action = $hashInput{"action"};
if($action eq "change_password"){
	change_password_form($loginname,$sessionID);
}
elsif($action eq "update_password"){
	$passwd = $hashInput{"passwd"};
	validate($loginname,$passwd);
	update_password($loginname);
	$action = "continue";
}
if(!file_exist("$PATH_INFO/$loginname-config.txt")){
	if(file_exist("$PATH_INFO/$loginname-config.backup") && $loginname !~ /^guest/){
		print "<b>WARNING:</b> Configuration file was recovered from backup!<br>\n";
		copy "$PATH_INFO/$loginname-config.backup", "$PATH_INFO/$loginname-config.txt";
	}else{
		my @files = glob("$PATH_DATA/$loginname-*");
		if(@files>0){
			print "Configuration file not found!<br>Contact sharoval(at)mail.nih.gov to restore it.\n";
		}else{
			copy "$PATH_DATA/default-config.txt", "$PATH_INFO/$loginname-config.txt";
		}
	}
	#`chmod 666 $PATH_INFO/$loginname-config.txt`;
}elsif($hashInput{"login"} && $loginname !~ /^guest/){
	copy "$PATH_INFO/$loginname-config.txt", "$PATH_INFO/$loginname-config.backup";
}
if($hashInput{"mainPage"} eq "mainPage"){
	check_configuration();
}

######################### MAIN SUBROUTINE SWITCH ####################################
if($hashInput{"terminate_task"}){	terminate_task(); }
elsif($hashInput{"file_delete"}){	file_delete(); }
elsif($hashInput{"file_download"}){	file_download(); }
elsif($hashInput{"remove_samples"}){	remove_samples(); }
if($action eq "file_edit1"){		file_edit1(); $action="continue"; }
if($action eq "file_edit2"){		file_edit2(); $action="continue"; }

if($action eq "matrix_explore"){ 	run_anova($RUN_ANOVA_FIRST); matrix_explore(); }
elsif($action eq "matrix_explore1"){	matrix_explore1(); }
elsif($action eq "geneset_explore"){	geneset_explore(); }
elsif($action eq "geneset_explore1"){	geneset_explore1(); }
elsif($action eq "correlation"){	correlation(); }
elsif($action eq "correlation1"){	correlation1(); }
elsif($action eq "geneset_overlap"){	geneset_overlap(); }
elsif($action eq "geneset_overlap1"){	geneset_overlap1(); }
elsif($action eq "geneset"){		geneset_analysis(); }
elsif($action eq "geneset1"){		geneset_analysis1(); }
elsif($action eq "matrix2geneset1"){	matrix2geneset1(); }
elsif($action eq "search_GEO"){		search_GEO(); }
elsif($action eq "search_GEO1"){	search_GEO1(); }
elsif($action eq "GEO_search_results"){	GEO_search_results(); }
elsif($action eq "open_samples"){	open_samples(); }
elsif($action eq "output_explore"){	output_explore(); }
elsif($action eq "output_explore1"){	output_explore1(); }
elsif($action eq "file_upload"){	file_upload(); }
elsif($action eq "file_edit"){		file_edit(); }
elsif($action eq "plot_matrix"){	plot_matrix(); }
elsif($action eq "generate_matrix"){	generate_matrix(); }
elsif($action eq "interrupt_program"){	interrupt_program(); }
elsif($action eq "gene_list_explore"){	gene_list_explore(); }
elsif($action eq "meta-analysis"){	meta_analysis_combine(); }
elsif($action eq "add_geneset"){	add_geneset(); }
elsif($action eq "show_gene_list"){	show_gene_list(); }
elsif($action eq "admin_job" && $loginname eq "administrator"){ admin_job(); }
elsif($hashInput{"login"} || $action eq "continue"){ main_page(); }
else{ error_message("Unknown command $action!"); }
exit(0);

#************************************
sub  main_page
#************************************
{
my @matrix_list;
my @geneset_list;
my @samples_list;
my @annotation_list;
my @annotation_GPL;
my %hashGPL;
my @output_list;
my %hashDefault;
my $update_config=0;
my $update_text;
my $organismID = $hashInput{"organismID"};
my $orgID_last;
my $nLines;

open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!","register");
while(my $line = <INFO>){
	$line =~ s/\n$//;
	my %hash=();
	if(length($line) < 3) { next; }
	$nLines++;
	read_config_line($line,\%hash);
	my $descr = $hash{"description"};
	my $org = $hash{"organismID"};
	my $date1 = $hash{"date"};
	my %hashDuplicate;
	if($line=~/^type_/){
		my @items = split(/[\t=]/,$line);
		if($hashDuplicate{"$items[0]=$items[1]"}){
			$update_config=1;
			next;
		}
		$hashDuplicate{"$items[0]=$items[1]"}++;
		if(!$date1){
			my $file_name1 = "$loginname-$items[1]";
			if($line=~/^type_annot/){ $file_name1 .= "_annot.txt"; }
			if(!file_exist("$PATH_DATA/$file_name1")){
				$update_config=1;
				next;
			}
			my $resp = `find $PATH_DATA/$file_name1 -maxdepth 0 -printf "%TY-%Tm-%Td"`;
			if($resp =~ /^\d\d\d\d-\d\d-\d\d$/){
				$date1 = $resp;
				$line .= "\tdate=$date1";
				$update_config=1;
			}
		}
	}
	if($line =~ /^mainPage=mainPage/){
		$line =~ s/^mainPage=mainPage/mainPage=$org/;
		if(!$organismID){ $organismID=$org; }
		$update_config=1;
	}elsif($line =~ /^mainPage/){
		my $orgID = $hash{"mainPage"};
		if(!$organismID || $orgID==$organismID){
			%hashDefault = %hash;
			$orgID_last = $orgID;
		}
	}
	$update_text .= $line."\n";
	my $file_matrix = $hash{"type_matrix"};
	if($file_matrix){ push(@matrix_list,[$file_matrix,$descr,$org,$date1]); next; }
	my $file_geneset = $hash{"type_geneset"};
	if($file_geneset){ push(@geneset_list,[$file_geneset,$descr,$org,$date1]); next; }
	my $file_samples = $hash{"type_samples"};
	if($file_samples){ push(@samples_list,[$file_samples,$descr,$org,$date1]); next; }
	my $file_annotation = $hash{"type_annotation"};
	if($file_annotation){
		if($file_annotation !~ /^GPL\d+$/){
			push(@annotation_list,[$file_annotation,$descr,$org,$date1]);
		}else{
			my $num=$file_annotation;
			$num =~ s/^GPL//;
			push(@annotation_GPL,[$file_annotation,$descr,$org,$date1,$num]);
			$hashGPL{$num}=1;
		}
		next;
	}
	my $file_output = $hash{"type_output"};
	if($file_output){ push(@output_list,[$file_output,$descr,$org,$date1]); next; }
}
close INFO;
if(!$organismID){
	if($orgID_last){ $organismID = $orgID_last; }
	else{ $organismID=10090; }
}
if($update_config){
	open (OUT,">$PATH_INFO/$loginname-config1.txt");
	print OUT "$update_text";
	close OUT;
	my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
	$nLines1 =~ s/\s+$//;
	#print "$nLines1 $nLines<br>\n";
	if($nLines1 && $nLines1 > $nLines*0.9){
		copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
	}else{
		error_message("Failed to update configuration file! $nLines1 $nLines");
	}
}
@matrix_list = sort {lc($a->[0]) cmp lc($b->[0])} @matrix_list;
@geneset_list = sort {lc($a->[0]) cmp lc($b->[0])} @geneset_list;
@samples_list = sort {lc($a->[0]) cmp lc($b->[0])} @samples_list;
@annotation_list = sort {lc($a->[0]) cmp lc($b->[0])} @annotation_list;
@output_list = sort {lc($a->[0]) cmp lc($b->[0])} @output_list;
if($loginname eq "public"){
	push(@annotation_list, sort {$a->[4]<=>$b->[4]} @annotation_GPL);
}elsif(open(INFO,"<$PATH_INFO/public-config.txt")){
	my @matrix_list1=();
	my @geneset_list1=();
	my @annotation_list1=();
	while(my $line = <INFO>){
		$line =~ s/\n$//;
		my %hash=();
		read_config_line($line,\%hash);
		my $descr = $hash{"description"};
		my $org = $hash{"organismID"};
		my $date1 = $hash{"date"};
		my $file_matrix = $hash{"type_matrix"};
		if($file_matrix){ push(@matrix_list1,["public-$file_matrix",$descr,$org,$date1]); next; }
		my $file_geneset = $hash{"type_geneset"};
		if($file_geneset){ push(@geneset_list1,["public-$file_geneset",$descr,$org,$date1]); next; }
		my $file_annotation = $hash{"type_annotation"};
		if($file_annotation){
			if($file_annotation !~ /^GPL\d+$/){
				push(@annotation_list1,["public-$file_annotation",$descr,$org,$date1]);
			}else{
				my $num=$file_annotation;
				$num =~ s/^GPL//;
				push(@annotation_GPL,["public-$file_annotation",$descr,$org,$date1,$num]);
				$hashGPL{$num}=1;
			}
		}
	}
	push(@matrix_list, sort {lc($a->[0]) cmp lc($b->[0])} @matrix_list1);
	push(@geneset_list, sort {lc($a->[0]) cmp lc($b->[0])} @geneset_list1);
	push(@annotation_list, sort {lc($a->[0]) cmp lc($b->[0])} @annotation_list1);
	push(@annotation_list, sort {$a->[4]<=>$b->[4]} @annotation_GPL);
}
close INFO;
filter_list_by_organism(\@matrix_list, $organismID);
filter_list_by_organism(\@geneset_list, $organismID);
filter_list_by_organism(\@samples_list, $organismID);
filter_list_by_organism(\@output_list, $organismID);
filter_list_by_organism(\@annotation_list, $organismID);
print "<HTML><HEAD><TITLE>ExAtlas</TITLE>\n";
print_header("update_description();");
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
my ($items,$descriptions) = get_array_lists(\@samples_list);
print "samples_list = new Array($items);\n";
print "samples_description = new Array($descriptions);\n";
($items,$descriptions) = get_array_lists(\@matrix_list);
print "matrix_list = new Array($items);\n";
print "matrix_description = new Array($descriptions);\n";
($items,$descriptions) = get_array_lists(\@geneset_list);
print "geneset_list = new Array($items);\n";
print "geneset_description = new Array($descriptions);\n";
($items,$descriptions) = get_array_lists(\@annotation_list);
print "annotation_list = new Array($items);\n";
print "annotation_description = new Array($descriptions);\n";
($items,$descriptions) = get_array_lists(\@output_list);
print "output_list = new Array($items);\n";
print "output_description = new Array($descriptions);\n";

print "function update_description() {\n";
print "  var index;\n";
print "  index = document.exatlas.file_samples.selectedIndex;\n";
print "  document.exatlas.description_samples.value = samples_description[index];\n";
print "  index = document.exatlas.file_matrix.selectedIndex;\n";
print "  document.exatlas.description_matrix.value = matrix_description[index];\n";
print "  index = document.exatlas.file_geneset.selectedIndex;\n";
print "  document.exatlas.description_geneset.value = geneset_description[index];\n";
print "  index = document.exatlas.file_output.selectedIndex;\n";
print "  document.exatlas.description_output.value = output_description[index];\n";
print "}\n";
print "function clear_marks() {\n";
print "	document.exatlas.target = \"\"\n";
print "	document.exatlas.file_download.value = \"\";\n";
print "	document.exatlas.file_delete.value = \"\";\n";
print "	document.exatlas.mainPage.value = \"mainPage\";\n";
print "}\n";
print "function do_analysis(command) {\n";
print "	clear_marks();\n";
print "	if(command != \"file_upload\"){\n";
print "		var x = Math.round(Math.random()*10000);\n";
print "		document.exatlas.target = \"_BLANK\"+x;\n";
print "	}\n";
print "	document.exatlas.action.value = command;\n";
print "	document.exatlas.submit();\n";
print "}\n";
print "function change_organism() {\n";
print "	clear_marks();\n";
print "	document.exatlas.action.value = \"continue\";\n";
print "	document.exatlas.mainPage.value = \"changeOrganism\";\n";
print "	document.exatlas.submit();\n";
print "}\n";
print "function download_file(file_type) {\n";
print "	var filename;\n";
print "	if(file_type == \"geneset\"){ filename = document.exatlas.select_geneset_file.options[document.exatlas.select_geneset_file.selectedIndex].value; }\n";
print "	else if(file_type == \"matrix\"){ filename = document.exatlas.select_matrix_file.options[document.exatlas.select_matrix_file.selectedIndex].value; }\n";
print "	else if(file_type == \"output\"){ filename = document.exatlas.select_output_file.options[document.exatlas.select_output_file.selectedIndex].value; }\n";
print "	else if(file_type == \"samples\"){ filename = document.exatlas.select_samples_file.options[document.exatlas.select_samples_file.selectedIndex].value; }\n";
print "	else if(file_type == \"annotation\"){ filename = document.exatlas.select_annotation_file.options[document.exatlas.select_annotation_file.selectedIndex].value+\"_annot.txt\"; }\n";
print "	if(!filename){ alert(\"File not selected\"); return false; }\n";
print "	clear_marks();\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.exatlas.target = \"_BLANK\"+x;\n";
print "	document.exatlas.action.value = \"\";\n";
print "	document.exatlas.file_download.value = filename;\n";
print "	document.exatlas.submit();\n";
print "}\n";
print "function edit_refresh() {\n";
print "	clear_marks();\n";
print "	document.exatlas.action.value = \"continue\";\n";
print "	document.exatlas.submit();\n";
print "}\n";
print "function delete_file(file_type) {\n";
print "	var filename;\n";
print "	if(file_type == \"samples\"){ filename = document.exatlas.select_samples_file.options[document.exatlas.select_samples_file.selectedIndex].value; }\n";
print "	else if(file_type == \"annotation\"){ filename = document.exatlas.select_annotation_file.options[document.exatlas.select_annotation_file.selectedIndex].value; }\n";
print "	else if(file_type == \"geneset\"){ filename = document.exatlas.select_geneset_file.options[document.exatlas.select_geneset_file.selectedIndex].value; }\n";
print "	else if(file_type == \"matrix\"){ filename = document.exatlas.select_matrix_file.options[document.exatlas.select_matrix_file.selectedIndex].value; }\n";
print "	else if(file_type == \"output\"){ filename = document.exatlas.select_output_file.options[document.exatlas.select_output_file.selectedIndex].value; }\n";
print "	if(!filename){ alert(\"File not selected\"); return false; }\n";
if($loginname !~ /^public/){
	print "if(filename.search(/^public-/) >=0){ alert(\"Public files cannot be deleted\"); return false; }\n";
}
print "	clear_marks();\n";
print "	document.exatlas.file_delete.value = filename;\n";
print "	document.exatlas.file_type.value = file_type;\n";
print "	document.exatlas.action.value = \"continue\";\n";
print "	document.exatlas.submit();\n";
print "}\n";
print "function edit_file(file_type) {\n";
print "	var filename;\n";
print "	if(file_type == \"samples\"){ filename = document.exatlas.select_samples_file.options[document.exatlas.select_samples_file.selectedIndex].value; }\n";
print "	else if(file_type == \"annotation\"){ filename = document.exatlas.select_annotation_file.options[document.exatlas.select_annotation_file.selectedIndex].value; }\n";
print "	else if(file_type == \"geneset\"){ filename = document.exatlas.select_geneset_file.options[document.exatlas.select_geneset_file.selectedIndex].value; }\n";
print "	else if(file_type == \"matrix\"){ filename = document.exatlas.select_matrix_file.options[document.exatlas.select_matrix_file.selectedIndex].value; }\n";
print "	else if(file_type == \"output\"){ filename = document.exatlas.select_output_file.options[document.exatlas.select_output_file.selectedIndex].value; }\n";
print "	if(!filename){ alert(\"File not selected\"); return false; }\n";
print "	clear_marks();\n";
print "	document.exatlas.file_edit.value = filename;\n";
print "	document.exatlas.file_type.value = file_type;\n";
print "	document.exatlas.action.value = \"file_edit\";\n";
print "	document.exatlas.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
my $x = int(10000*rand());
print "<p>Use \"<font color=3000AA>Refresh</font>\" button if you don't see newly created files. &nbsp; HELP: <a href=../exatlas-help.html TARGET=\"_blank$x\">A guide to ExAtlas</a><br>\n";
print "<HR NOSHADE>\n";
print "<FORM NAME=exatlas ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST LANGUAGE=\"javascript\" onsubmit=\"return exatlas_onsubmit()\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=hidden VALUE=\"continue\">\n";
print "<INPUT NAME=\"file_download\" TYPE=hidden VALUE=\"\">\n";
print "<INPUT NAME=\"file_delete\" TYPE=hidden VALUE=\"\">\n";
print "<INPUT NAME=\"file_edit\" TYPE=hidden VALUE=\"\">\n";
print "<INPUT NAME=\"file_type\" TYPE=hidden VALUE=\"\">\n";
print "<INPUT NAME=\"mainPage\" TYPE=hidden VALUE=\"mainPage\">\n";
print "<TABLE border=0>\n";
print "<TR><TD><INPUT TYPE=\"button\" VALUE=\"Refresh\" onClick=\"edit_refresh();\">\n";
print " &nbsp; &nbsp; <img src=\"../images/arrow_right.gif\" border=0><TD><select name=organismID style=width:260px; onChange=change_organism();>\n";
foreach my $ref (@organisms){
	print "<option value=$ref->[0]";
	if($ref->[0] == $organismID){ print " selected"; }
	print ">$ref->[3] ($ref->[2])\n";
}
print "</select><TD COLSPAN=2>(page is reloaded after change)\n";
print "<TR><TD><center><FONT SIZE=+1><b>Select Data Files<TD><b><center>Description<TD><TD><b>File type\n";
print "<TR><TD><select name=\"file_matrix\" style=\"width: 300px;\" onChange=update_description();>\n";
for(my $i=0; $i<@matrix_list; ++$i){ 
	print "<option value=\"$matrix_list[$i]->[0]\"";
	if($matrix_list[$i]->[0] eq $hashDefault{"file_matrix"}){ print " selected"; }
	print "> $matrix_list[$i]->[0] $matrix_list[$i]->[3]\n";
}
print "</select><td><INPUT NAME=\"description_matrix\" style=width:260px;><TD>\n";
if(@matrix_list){
	print "<INPUT TYPE=\"button\"  VALUE=\" Open \" onClick=\"do_analysis('matrix_explore');\">\n";
}
print "<TD> Expression profiles\n";
print "<TR><TD><select name=\"file_geneset\" style=\"width: 300px;\" onChange=update_description();>\n";
for(my $i=0; $i<@geneset_list; ++$i){ 
	print "<option value=\"$geneset_list[$i]->[0]\"";
	if($geneset_list[$i]->[0] eq $hashDefault{"file_geneset"}){ print " selected"; }
	print "> $geneset_list[$i]->[0] $geneset_list[$i]->[3]\n";
}
print "</select><td><INPUT NAME=\"description_geneset\" style=width:260px;><TD>\n";
if(@geneset_list){
	print "<INPUT TYPE=\"button\" VALUE=\" Open \" onClick=\"do_analysis('geneset_explore');\">\n";
}
print "<TD> Gene set file\n";
print "<TR><TD><select name=\"file_output\" style=\"width: 300px;\" onChange=update_description();>\n";
for(my $i=0; $i<@output_list; ++$i){ 
	print "<option value=\"$output_list[$i]->[0]\"";
	if($output_list[$i]->[0] eq $hashDefault{"file_output"}){ print " selected"; }
	print "> $output_list[$i]->[0] $output_list[$i]->[3]\n";
}
print "</select><td><INPUT NAME=\"description_output\" style=width:260px;><TD>\n";
if(@output_list){
	print "<INPUT TYPE=\"button\"  VALUE=\" Open \" onClick=\"do_analysis('output_explore');\">\n";
}
print "<TD>Output file\n";
print "<TR><TD><select name=\"file_samples\" style=\"width: 300px;\" onChange=update_description();>\n";
for(my $i=0; $i<@samples_list; ++$i){ 
	print "<option value=\"$samples_list[$i]->[0]\"";
	if($samples_list[$i]->[0] eq $hashDefault{"file_samples"}){ print " selected"; }
	print "> $samples_list[$i]->[0] $samples_list[$i]->[3]\n";
}
print "</select><td><INPUT NAME=\"description_samples\" style=width:260px;><TD>\n";
if(@samples_list){
	print "<INPUT TYPE=\"button\"  VALUE=\" Open \" onClick=\"do_analysis('open_samples');\">\n";
}
print "<TD> Samples file\n";
print "<TR><TD><INPUT TYPE=\"button\"  VALUE=\"Find samples in GEO\" style=\"width: 300px;\" onClick=\"do_analysis('search_GEO');\">\n";
print "<TD COLSPAN=3>Extract samples from GEO database\n";
print "<TR><TD><INPUT TYPE=\"button\"  VALUE=\"Upload data file\" style=\"width: 300px;\" onClick=\"do_analysis('file_upload');\">\n";
print "<TD COLSPAN=3>Upload data files\n";
print "</TABLE><p>\n";

print "<TABLE border=0>\n";
print "<TR><TD><center><FONT SIZE=+1><b>File Management\n";
print "<TR><TD><select name=\"select_matrix_file\" style=\"width: 300px;\">\n";
print "<option value=\"\"> ----------------------- select -----------------------\n";
for(my $i=0; $i<@matrix_list; ++$i){ 
	print "<option value=\"$matrix_list[$i]->[0]\"> $matrix_list[$i]->[0] $matrix_list[$i]->[3]\n";
}
print "</select><TD>Expression profiles\n";
print "<td><INPUT TYPE=button VALUE=Download onClick=download_file('matrix');>\n";
print "<td><INPUT TYPE=button VALUE=Delete onClick=delete_file('matrix');>\n";
print "<td><INPUT TYPE=button VALUE=Edit onClick=edit_file('matrix');>\n";

print "<TR><TD><select name=select_geneset_file style=\"width: 300px;\">\n";
print "<option value=\"\"> ----------------------- select -----------------------\n";
for(my $i=0; $i<@geneset_list; ++$i){ 
	print "<option value=\"$geneset_list[$i]->[0]\"> $geneset_list[$i]->[0] $geneset_list[$i]->[3]\n";
}
print "</select><TD>Gene set files\n";
print "<td><INPUT TYPE=button VALUE=Download onClick=download_file('geneset');>\n";
print "<td><INPUT TYPE=button VALUE=Delete onClick=delete_file('geneset');>\n";
print "<td><INPUT TYPE=button VALUE=Edit onClick=edit_file('geneset');>\n";

print "<TR><TD><select name=select_output_file style=\"width: 300px;\">\n";
print "<option value=\"\"> ----------------------- select -----------------------\n";
for(my $i=0; $i<@output_list; ++$i){ 
	print "<option value=\"$output_list[$i]->[0]\"> $output_list[$i]->[0] $output_list[$i]->[3]\n";
}
print "</select><TD>Output files\n";
print "<td><INPUT TYPE=button VALUE=Download onClick=download_file('output');>\n";
print "<td><INPUT TYPE=button VALUE=Delete onClick=delete_file('output');>\n";
print "<td><INPUT TYPE=button VALUE=Edit onClick=edit_file('output');>\n";

print "<TR><TD><select name=select_samples_file style=\"width: 300px;\">\n";
print "<option value=\"\"> ----------------------- select -----------------------\n";
for(my $i=0; $i<@samples_list; ++$i){ 
	print "<option value=\"$samples_list[$i]->[0]\"> $samples_list[$i]->[0] $samples_list[$i]->[3]\n";
}
print "</select><TD>Samples files\n";
print "<td><INPUT TYPE=button VALUE=Download onClick=download_file('samples');>\n";
print "<td><INPUT TYPE=button VALUE=Delete onClick=delete_file('samples');>\n";
print "<td><INPUT TYPE=button VALUE=Edit onClick=edit_file('samples');>\n";

print "<TR><TD><select name=select_annotation_file style=\"width: 300px;\">\n";
print "<option value=\"\"> ----------------------- select -----------------------\n";
for(my $i=0; $i<@annotation_list; ++$i){ 
	print "<option value=\"$annotation_list[$i]->[0]\"> $annotation_list[$i]->[0] $annotation_list[$i]->[3]\n";
}
print "</select><TD>Annotation files\n";
print "<td><INPUT TYPE=button VALUE=Download onClick=download_file('annotation');>\n";
print "<td><INPUT TYPE=button VALUE=Delete onClick=delete_file('annotation');>\n";
print "</TABLE>\n";
if($loginname eq "administrator"){
	print "<INPUT TYPE=button  VALUE=\" Start administrator job \" onClick=\"do_analysis('admin_job');\">\n";
	print "<INPUT NAME=admin_param SIZE=40>\n";
}
print "</FORM>\n";
print "<p><HR NOSHADE><p>\n";
if($loginname !~ /^guest\d+/i){
	my $x = int(10000*rand());
	print "<p><FORM NAME=change_pas ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST TARGET=\"_blank$x\">\n";
	print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
	print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=change_password>\n";
	print "<INPUT NAME=\"organismID\" TYPE=\"hidden\" VALUE=\"$organismID\">\n";
	print "<INPUT TYPE=\"submit\" VALUE=\"Change password\">\n";
	print "</table></FORM>\n";
}
print "<i>Please report any problems to <a href=mailto:sharoval\@mail.nih.gov>webmaster</a><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  correlation
#**************************************
{
my @output_list;
my @matrix_list;
my @file_list;

my $file_matrix = $hashInput{"file_matrix"};
my $description_matrix = $hashInput{"description_matrix"};
my $organismID = $hashInput{"organismID"};
my $organismID1 = $hashInput{"organismID1"};
if(!$organismID1){ $organismID1=$organismID; }
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Config file not found!");
while(my $line = <INFO>){
	$line =~ s/\n$//;
	my %hash=();
	if(length($line) < 3) { next; }
	my @items = split(/[=\t]/,$line);
	read_config_line($line,\%hash);
	my $file_matrix1 = $hash{"type_matrix"};
	if($file_matrix1){
		push(@matrix_list,[$file_matrix1,$hash{"description"},$hash{"organismID"},$hash{"date"}]);
	}
	if($items[0] =~ /^type_output/){
		push(@output_list,$items[1]);
	}elsif($items[0] =~ /^type_/){
		push(@file_list,$items[1]);
	}
}
close INFO;
@matrix_list = sort {lc($a->[0]) cmp lc($b->[0])} @matrix_list;
if($loginname ne "public" && open(INFO,"<$PATH_INFO/public-config.txt")){
	my @matrix_list1;
	while(my $line = <INFO>){
		$line =~ s/\n$//;
		my %hash=();
		if(length($line) < 3) { next; }
		read_config_line($line,\%hash);
		my $file_matrix1 = $hash{"type_matrix"};
		if($file_matrix1){
			push(@matrix_list1,["public-$file_matrix1",$hash{"description"},$hash{"organismID"},$hash{"date"}]);
		}
	}
	push(@matrix_list, sort {lc($a->[0]) cmp lc($b->[0])} @matrix_list1);
}
close INFO;
my %hashOrganismID;
foreach my $ref (@matrix_list){
	$hashOrganismID{$ref->[2]}=1;
}
filter_list_by_organism(\@matrix_list, $organismID1);
my @anova_headers;
for(my $i=0; $i<@matrix_list; ++$i){
	my $file = $matrix_list[$i]->[0];
	my $file_anova1 = $file;
	$file_anova1 =~ s/^public-/public-anova-/;
	if($file_anova1 !~ /^public-/){
		$file_anova1 = "$loginname-anova-$file";
	}
	my @headers = get_anova_headers($file_anova1,1);
	if(!@headers){
		splice(@matrix_list,$i,1);
		$i--;
		next;
	}
	my $headers_list="\"Median profile\"";
	for(my $j=0; $j<@headers; ++$j){
		$headers_list .= ",\"".$headers[$j]."\""; 
	}
	push(@anova_headers,$headers_list);
}

my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_anova !~ /^public-/){
	$file_anova = "$loginname-anova-$file_matrix";
}
my @headers = get_anova_headers($file_anova,1);
if(!@headers){ error_message("Invalid ANOVA file"); }
my @anova_headers_file1 = ("Median profile",@headers);

my $file_list="";
my $output_list="";
foreach my $name (@output_list){
	if(!$output_list){ $output_list = "\"".$name."\""; }
	else{ $output_list .= ",\"".$name."\""; }
}
foreach my $name (@file_list){
	if(!$file_list){ $file_list = "\"".$name."\""; }
	else{ $file_list .= ",\"".$name."\""; }
}
my $logFileID = get_outputID(1);
my ($items,$descriptions) = get_array_lists(\@matrix_list);

# Print page header
print "<HTML><HEAD><TITLE>ExAtlas - correlation</TITLE>\n";
print_header("update_description();");
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "var output_list = new Array($output_list);\n";
print "var file_list = new Array($file_list);\n";
print "var matrix_list = new Array($items);\n";
print "var matrix_description = new Array($descriptions);\n";
print "var anova_headers = new Array();\n";
for(my $i=0; $i<@anova_headers; ++$i){
	print "anova_headers[$i] = new Array($anova_headers[$i]);\n";
}
print "function update_description() {\n";
print "	var index = document.form_correlation.file_matrix1.selectedIndex;\n";
print "	document.form_correlation.baseline1.options.length=0;\n";
print "	for(i=0; i<anova_headers[index].length; ++i){\n";
print "		var x = new Option(anova_headers[index][i],i);\n";
print "		document.form_correlation.baseline1.options[i] = x;\n";
print "	}\n";
print "	document.form_correlation.baseline1.selectedIndex=0;\n";
print "	document.form_correlation.description_matrix1.value = matrix_description[index];\n";
print "}\n";
print "function change_organism() {\n";
print "	document.form_correlation.action.value = \"correlation\";\n";
print "	document.form_correlation.submit();\n";
print "}\n";
print "function alert_onsubmit() {\n";
print "	document.form_correlation.action.value = \"correlation1\";\n";
print "	document.form_correlation.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";

my @FDR_list = (1,0.5,0.2,0.1,0.05,0.01,0.001,0.0001);
my @fold_change = (1,1.5,2,3,4,5,10);
my ($FDR,$fold_change) = (0.05,2);
my @filter_expr = (0,0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,100,300,1000,3000);
print "<H2>Estimate correlation between expression profiles</H2>\n";
print "<FORM NAME=form_correlation ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD WIDTH=60><TD><b>File name<TD><b>Description<TD><b>Organism\n";
print "<TR><TD><b>File 1:<TD>$file_matrix<TD WIDTH=350>$description_matrix<TD>$hashOrganism{$organismID}\n";
print "<TR><TD><TD><center><i>(Select second expression profile data)</i>\n";
print "<TR><TD><b>File 2:<TD><select name=file_matrix1 style=width:300px; onChange=update_description();>\n";
for(my $i=0; $i<@matrix_list; ++$i){
	print "<option value=\"$matrix_list[$i]->[0]\"> $matrix_list[$i]->[0] $matrix_list[$i]->[3]\n";
}
print "</select><TD><INPUT NAME=description_matrix1 SIZE=40><TD>$hashOrganism{$organismID1}\n";
print "</TABLE>\n";
if(keys %hashOrganismID > 1){
	print "<p>Change organism for 2nd file:\n";
	print "<TD><select name=organismID1 style=width:200px; onChange=change_organism();>\n";
	foreach my $ref (@organisms){
		if(!$hashOrganismID{$ref->[0]}){ next; }
		print "<option value=$ref->[0]";
		if($ref->[0] == $organismID1){ print " selected"; }
		print ">$ref->[3] ($ref->[2])\n";
	}
	print "</select> (page is reloaded after change)\n";
}else{
	print "<INPUT NAME=organismID1 TYPE=hidden VALUE=$organismID1>\n";
}
print "<p><TABLE BORDER=0>\n";
print "<TR><TD><b>Parameters<TD><center><b>File 1<TD WIDTH=15><TD><center><b>File 2\n";
print "<TR><TD>Use as baseline\n";
print "<TD><select name=baseline style=\"width: 150px;\">\n";
for(my $i=0; $i<@anova_headers_file1; ++$i){ 
	print "<option value=$i> $anova_headers_file1[$i]\n";
}
print "</select><TD>\n";
print "<TD><select name=baseline1 style=\"width: 150px;\"></select>\n";
print "<TR><TD><a href=../exatlas-help.html#fdr>FDR</a> threshold\n";
print "<TD><select name=FDR style=\"width: 150px;\">\n";
for(my $i=0; $i<@FDR_list; ++$i){ 
	print "<option value=$FDR_list[$i]"; if($FDR==$FDR_list[$i]){ print " selected"; } print "> $FDR_list[$i]\n";
}
print "</select><TD>\n";
print "<TD><select name=FDR1 style=\"width: 150px;\">\n";
for(my $i=0; $i<@FDR_list; ++$i){ 
	print "<option value=$FDR_list[$i]"; if($FDR==$FDR_list[$i]){ print " selected"; } print "> $FDR_list[$i]\n";
}
print "</select>\n";
print "<TR><TD>Fold change threshold\n";
print "<TD><select name=fold_change style=\"width: 150px;\">\n";
for(my $i=0; $i<@fold_change; ++$i){ 
	print "<option value=$fold_change[$i]"; if($fold_change==$fold_change[$i]){ print " selected"; } print "> $fold_change[$i]\n";
}
print "</select><TD>\n";
print "<TD><select name=fold_change1 style=\"width: 150px;\">\n";
for(my $i=0; $i<@fold_change; ++$i){ 
	print "<option value=$fold_change[$i]"; if($fold_change==$fold_change[$i]){ print " selected"; } print "> $fold_change[$i]\n";
}
print "</select>\n";
print "<TR><TD>Expression threshold\n";
print "<TD><select name=expr_thresh style=\"width: 150px;\">\n";
foreach my $x (@filter_expr){ print "<option value=$x> $x\n"; }
print "</select><TD>\n";
print "<TD><select name=expr_thresh1 style=\"width: 150px;\">\n";
foreach my $x (@filter_expr){ print "<option value=$x> $x\n"; }
print "</select>\n";

print "<TR><TD>&nbsp;<TR><TD>Correlation method\n";
print "<TD><select name=correlation_method style=\"width: 150px;\">\n";
for(my $i=0; $i<@correlation_list; ++$i){
	my $value = $i+1;
	print "<option value=$value> $correlation_list[$i]\n";
}
print "</select>\n";
print "<TD><TD><INPUT NAME=option_genes TYPE=checkbox> Identify coregulated genes\n";
my @EPFP_list = (0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6);
print "<TR><TD>Restrict values\n";
print "<TD><select name=restrict style=\"width: 150px;\">\n";
print "<option value=0> No restriction\n";
print "<option value=1> x > 0\n";
print "<option value=2> x > 0 and y > 0\n";
print "<option value=3> x < 0\n";
print "<option value=4> x < 0 and y < 0\n</select>\n";
print "<TD><TD><select name=coregulated_epfp style=\"width: 60px;\">\n";
for(my $i=0; $i<@EPFP_list; ++$i){
	my $x = $EPFP_list[$i];
	print "<option value=$x"; if($x==0.5){ print " selected"; } print "> $x\n";
}
print "</select> <a href=../exatlas-help.html#epfp>EPFP</a> (for coregulated genes)\n";
print "<TR><TD>Direction of change (file #2)\n";
print "<TD><select name=direction style=\"width: 150px;\">\n";
print "<option value=1> Original\n";
print "<option value=-1> Reversed\n</select>\n";
print "<TD><TD><INPUT NAME=option_subtract_baseline TYPE=checkbox checked> Subtract baseline for each gene\n";
print "</TABLE><p>\n";
print "<b>Notes:</b><br>(1) Use FDR threshold and/or fold change threshold to limit the number of genes for analysis.\n";
print "Lower values of FDR  and higher values of fold change correspond to more stringent filtering.<br>\n";
print "(2) If direction of change = \"Reversed\" then gene expression change for File #2 is inverted (multiplied by -1).<br>\n";
print "(3) In the option \"Restrict values\", x = values from file #1 and y = values from file #2.<br>\n";
print "(4) For algorithm see <a href=../exatlas-help.html#correlation>help file</a>.<br>\n";
print "(5) Option to subtract baseline expression of each gene (default) is used to compare gene expression change relative to baseline\n";
print "To compare absolute gene expression values, uncheck this option.<p>\n";
print "<INPUT TYPE=button value=\"Estimate correlation matrix\" style=\"width: 200px;\" onClick=\"alert_onsubmit();\">\n";
print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=action TYPE=hidden VALUE=correlation1>\n";
print "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
print "<INPUT NAME=description_matrix TYPE=hidden VALUE=\"$description_matrix\">\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=\"$organismID\">\n";
print "<INPUT NAME=logFileID TYPE=hidden VALUE=\"$logFileID\">\n";
print "</FORM><p>\n";
print "<HR NOSHADE></HR>\n";
print "<p><i>Please report any problems to <a href=mailto:sharoval\@mail.nih.gov>webmaster</a><p>\n";
print "<INPUT TYPE=button VALUE=\"Cancel (close window)\" style=\"width: 200px;\" LANGUAGE=javascript onClick=window.close();><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  correlation1
#**************************************
{
my $logFileID = $hashInput{"logFileID"};
my $file_matrix = $hashInput{"file_matrix"};
my $file_matrix1 = $hashInput{"file_matrix1"};
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_anova = "$loginname-anova-$file_matrix";
}
open(INFO,"<$PATH_DATA/$file_anova") or error_message("No ANOVA file #1");
my $line = <INFO>;
chop $line;
close INFO;
my ($junk,$junk1,@headers) =split(/\t/,$line);
if(!@headers){ error_message("Invalid ANOVA file #1"); }
my $nCol=0;
while($nCol<@headers && $headers[$nCol] !~ /^Var\(/){ ++$nCol; }
my $file_anova1 = $file_matrix1;
$file_anova1 =~ s/^public-/public-anova-/;
if($file_matrix1 !~ /^public-/){
	$file_anova1 = "$loginname-anova-$file_matrix1";
}
my $logFileID = get_outputID(3);  #0-log,1-web,2-output
$hashInput{"logFileID"} = $logFileID;
$hashInput{"runID"} = $RUN_CORRELATION;
if(!open(INFO,"<$PATH_DATA/$file_anova1")){
	interrupt_program(1);
	exit(0);
}
$line = <INFO>;
chop $line;
close INFO;
($junk,$junk1,@headers) =split(/\t/,$line);
if(!@headers){ error_message("Invalid ANOVA file #2"); }
my $nRow=0;
while($nRow<@headers && $headers[$nRow] !~ /^Var\(/){ ++$nRow; }
if($nCol<=0 || $nRow<=0){ error_message("Problem in input files: $nCol $nRow"); }
my $scale=10000;
if($hashInput{"option_genes"} eq "on"){ $scale *= 0.2; }
interrupt_program($nCol*$nRow/$scale);
exit(0);
}

#**************************************
sub  correlation2
#**************************************
{
my $logFileID = shift;
my $webPageID = $hashInput{"logFileID"}+1;
my $outputID = $webPageID+1;

my $organismID = $hashInput{"organismID"};
my $organismID1 = $hashInput{"organismID1"};
if(!$organismID1){ $organismID1=$organismID; }
my $file_matrix = $hashInput{"file_matrix"};
my $file_matrix1 = $hashInput{"file_matrix1"};
my $FDR_thresh = $hashInput{"FDR"};
my $FDR_thresh1 = $hashInput{"FDR1"};
my $expr_thresh = $hashInput{"expr_thresh"};
my $expr_thresh1 = $hashInput{"expr_thresh1"};
my $fold_change_thresh = $hashInput{"fold_change"};
my $fold_change_thresh1 = $hashInput{"fold_change1"};
my $correlation_method = $hashInput{"correlation_method"};
my $baseline = $hashInput{"baseline"};
my $baseline1 = $hashInput{"baseline1"};
my $file_description = $hashInput{"description"};
my $description_matrix = $hashInput{"description_matrix"};
my $description_matrix1 = $hashInput{"description_matrix1"};
my $direction = $hashInput{"direction"};
my $restrict = $hashInput{"restrict"};
my $option_subtract_baseline = 0;
my $option_genes = 0;
if($hashInput{"option_genes"} eq "on"){ $option_genes = 1; }
if($hashInput{"option_subtract_baseline"} eq "on"){ $option_subtract_baseline = 1; }
my $log10 = log(10);
my $logratio_thresh = 0;
if($fold_change_thresh>1){ $logratio_thresh=log($fold_change_thresh)/$log10; }
my $logratio_thresh1 = 0;
if($fold_change_thresh1>1){ $logratio_thresh1=log($fold_change_thresh1)/$log10; }
my $expr_log_thresh=$MISSING; if($expr_thresh>0){ $expr_log_thresh=log($expr_thresh)/$log10; }
my $expr_log_thresh1=$MISSING; if($expr_thresh1>0){ $expr_log_thresh1=log($expr_thresh1)/$log10; }
my $baseline_text;
my $baseline_text1;

my $file_matrix_full = $file_matrix;
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_matrix_full = "$loginname-$file_matrix";
	$file_anova = "$loginname-anova-$file_matrix";
}
my $file_anova1 = $file_matrix1;
$file_anova1 =~ s/^public-/public-anova-/;
my $file_matrix_full1 = $file_matrix1;
if($file_matrix1 !~ /^public-/){
	$file_matrix_full1 = "$loginname-$file_matrix1";
	$file_anova1 = "$loginname-anova-$file_matrix1";
}
if(!file_exist("$PATH_DATA/$file_anova1")){
	if($file_matrix1 =~ /^public-/){ error_message("Missing public ANOVA file #2",$logFileID); }
	run_anova1($file_matrix1,$logFileID);
}
if($logFileID){
	file_append("Processing ANOVA output for $file_matrix","$PATH_OUTPUT/$logFileID.txt");
}
open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file #1",$logFileID);
my $line = <INFO>;
chop $line;
my ($junk,$junk1,@headers) =split(/\t/,$line);
if(!@headers){ error_message("Headers not found in $file_anova",$logFileID); }
my $nCol=0;
my $nRep=0;
while($nCol<@headers && $headers[$nCol] !~ /^Var\(/){ ++$nCol; }
for(my $i=0; $i<$nCol; ++$i){
	if($headers[$i]=~/ \(\d+\)$/ && $headers[$i] !~ /^Mean/){
		my @items = split(/ \(/,$headers[$i]);
		my $n = pop(@items);
		my $name = join(" (",@items);
		$n =~ s/\)$//;
		$nRep += $n;
		$headers[$i] = $name;
	}else{
		$headers[$i] =~ s/^Mean\(//;
		$headers[$i] =~ s/\)$//;
		$nRep++;
	}
}
splice(@headers,$nCol);
my $baseline_value;
if($baseline==0){
	$baseline_value = "Median profile";
}else{
	$baseline_value = $headers[$baseline-1];
}
my @Data;
my %symbols;
my %symbolsUC;
my @score;
my $irow=0;
my $log_thresh = log($fold_change_thresh)/$log10;
my %hashGenes;
my %hashAlias;
get_official_symbols($organismID,\%hashGenes,\%hashAlias);
while(my $line = <INFO>){
	chop $line;
	my ($junk,$expr,@data) =split(/\t/,$line);
	my $symbol = $data[$nCol+8];
	if(!$symbol){ next; }
	if($expr_thresh>0 && $expr<$expr_log_thresh){ next; }
	if(!$hashGenes{$symbol}){
		my $symbol1 = $hashAlias{$symbol};
		if($symbol1){ $symbol=$symbol1; }
	}
	my $FDR = $data[$nCol+6];
	if($FDR_thresh>0 && $FDR > $FDR_thresh){ next; }
	my $F = $data[$nCol+4];
	splice(@data,$nCol);
	if($logratio_thresh>0){
		my @sorted = sort {$a<=>$b} @data;
		while($sorted[0]<=$MISSING){ shift(@sorted); }
		if(@sorted>1){
			my $logratio = $sorted[@sorted-1] - $sorted[0];
			if($logratio < $logratio_thresh){ next; }
		}
	}
	push(@score,$F);
	my $median;
	if($baseline==0){
		$median = median(\@data);
	}else{
		$median = $data[$baseline-1];
		if($median==$MISSING){ $median = median(\@data); }
	}
	my $ref = $symbols{$symbol};
	if($ref){
		if($F > $ref->[1]){
			my $iii = $ref->[0];
			for(my $i=0; $i<$nCol; ++$i){
				if($data[$i]==$MISSING){ $Data[$i]->[$iii] = $MISSING; }
				else{
					my $x = $data[$i];
					if($option_subtract_baseline){
						$x = floor(10000*($data[$i]-$median)+0.5)/10000;
					}
					$Data[$i]->[$iii] = $x;
				}
			}
			$ref->[1] = $F;
		}
	}else{
		for(my $i=0; $i<$nCol; ++$i){
			if($data[$i]==$MISSING){ push(@{$Data[$i]},$MISSING); }
			else{
				my $x = $data[$i];
				if($option_subtract_baseline){
					$x = floor(10000*($data[$i]-$median)+0.5)/10000;
				}
				push(@{$Data[$i]},$x);
			}
		}
		$symbols{$symbol} = [$irow,$F];
		my $symbol1 = uc($symbol);
		$symbolsUC{$symbol1} = $symbol;
		$irow++;
	}
}
close INFO;
my $Ngenes = $irow;
if($Ngenes==0){
	my $text = "<H3>No genes match selection criteia for file#1</H3>Consider using less stringent criteria. For example, set Fold change threshold = 1 and FDR threshold = 0.5";
	if($logFileID){ error_message($text,$logFileID); }
	terminal_window($text);
}
my %hashRows;
foreach my $symbol (keys %symbols){
	my ($ii,$expr) = @{$symbols{$symbol}};
	$hashRows{$symbol}=[$ii,-1];
}
if($logFileID){
	file_append("Processing ANOVA output for $file_matrix1","$PATH_OUTPUT/$logFileID.txt");
}
my %geneHomolog=();
my %hashGenes1;
my %hashAlias1;
if($organismID1 != $organismID){
	get_gene_homolog($organismID1,$organismID,\%geneHomolog,$logFileID);
	get_official_symbols($organismID1,\%hashGenes1,\%hashAlias1);
}
open(INFO,"<$PATH_DATA/$file_anova1") or error_message("Cannot open ANOVA file",$logFileID);
$line = <INFO>;
chop $line;
my ($junk2,$junk3,@headers1) =split(/\t/,$line);
if(!@headers1){ error_message("Headers not found in $file_anova",$logFileID); }
my $nCol1=0;
my $nRep1=0;
while($nCol1<@headers1 && $headers1[$nCol1] !~ /^Var\(/){ ++$nCol1; }
for(my $i=0; $i<$nCol1; ++$i){
	if($headers1[$i]=~/ \(\d+\)$/ && $headers1[$i] !~ /^Mean/){
		my @items = split(/ \(/,$headers1[$i]);
		my $n = pop(@items);
		my $name = join(" (",@items);
		$n =~ s/\)$//;
		$nRep1 += $n;
		$headers1[$i] = $name;
	}else{
		$headers1[$i] =~ s/^Mean\(//;
		$headers1[$i] =~ s/\)$//;
		$nRep1++;
	}
}
splice(@headers1,$nCol1);
my $baseline_value1;
if($baseline1==0){
	$baseline_value1 = "Median profile";
}else{
	$baseline_value1 = $headers1[$baseline1-1];
}
my @Data1;
my %symbols1;
@score=();
$irow=0;
my $log_thresh1 = log($fold_change_thresh1)/$log10;
my ($n1,$n2)=(0,0);
while(my $line = <INFO>){
	chop $line;
	my ($id,$expr,@data) =split(/\t/,$line);
	my $symbol = $data[$nCol1+8];
	if(!$symbol){ next; }
	if($expr_thresh1>0 && $expr<$expr_log_thresh1){ next; }
	if($organismID1==$organismID){
		if(!$hashGenes{$symbol}){
			my $symbol1 = $hashAlias{$symbol};
			if($symbol1){ $symbol=$symbol1; }
		}
	}else{
		if(!$hashGenes1{$symbol}){
			my $symbol1 = $hashAlias1{$symbol};
			if($symbol1){ $symbol=$symbol1; }
		}
		my $symbol1 = $geneHomolog{$symbol};
		if($symbol1){
			$symbol = $symbol1;
			$n1++;
		}else{
			$symbol1 = uc($symbol);
			if($symbolsUC{$symbol1}){
				$symbol = $symbolsUC{$symbol1};
				$n2++;
			}
		}
		if(!$symbol){ next; }
	}
	my $FDR = $data[$nCol1+6];
	if($FDR_thresh1>0 && $FDR > $FDR_thresh1){ next; }
	my $F = $data[$nCol+4];
	splice(@data,$nCol1);
	if($logratio_thresh1>0){
		my @sorted = sort {$a<=>$b} @data;
		while($sorted[0]<=$MISSING){ shift(@sorted); }
		if(@sorted>1){
			my $logratio = $sorted[@sorted-1] - $sorted[0];
			if($logratio < $logratio_thresh1){ next; }
		}
	}
	push(@score,$F);
	my $median;
	if($baseline1==0){ $median = median(\@data); }
	else{ $median = $data[$baseline1-1]; }
	my $ref = $symbols1{$symbol};
	if($ref){
		if($F > $ref->[1]){
			my $iii = $ref->[0];
			for(my $i=0; $i<$nCol1; ++$i){
				if($data[$i]==$MISSING){ $Data1[$i]->[$iii] = $MISSING; }
				else{		
					my $x = $data[$i];
					if($option_subtract_baseline){
						$x = floor(10000*($data[$i]-$median)+0.5)/10000;
					}
					if($direction<0){ $x=-$x; }
					$Data1[$i]->[$iii] = $x;
				}
			}
			$ref->[1] = $F;
		}
	}else{
		for(my $i=0; $i<$nCol1; ++$i){
			if($data[$i]==$MISSING){ push(@{$Data1[$i]},$MISSING); }
			else{
				my $x = $data[$i];
				if($option_subtract_baseline){
					$x = floor(10000*($data[$i]-$median)+0.5)/10000;
				}
				if($direction<0){ $x=-$x; }
				push(@{$Data1[$i]},$x);
			}
		}
		$symbols1{$symbol} = [$irow,$F];
		$irow++;
	}
}
close INFO;
my $Ngenes1 = $irow;
if($Ngenes1==0){
	my $text = "<H3>No genes match selection criteia for file#2</H3>Consider using less stringent criteria. For example, set Fold change threshold = 1 and FDR threshold = 0.5";
	if($logFileID){ error_message($text,$logFileID); }
	terminal_window($text);
}
my $Ngenes_common = 0;
foreach my $symbol (keys %symbols1){
	my ($ii1,$F1) = @{$symbols1{$symbol}};
	my $ref = $hashRows{$symbol};
	if($ref){
		my ($ii,$F) = @{$symbols{$symbol}};
		$ref->[1]=$ii1;
		$Ngenes_common++;
	}
}
if($Ngenes_common<30){
	my $text = "<H3>Not enough common significant genes between files#1 & $2 (N=$Ngenes_common)</H3>Consider using less stringent criteria. For example, set Fold change threshold = 1 and FDR threshold = 0.5";
	if($logFileID){ error_message($text,$logFileID); }
	terminal_window($text);
}
if($organismID1 != $organismID && $logFileID){
	file_append("Gene symbol conversion. From homologene: N1=$n1; Symbol match: N2=$n2","$PATH_OUTPUT/$logFileID.txt");
}
if($logFileID){
	file_append("Number of genes from expression profiles data #1=\t$Ngenes\nNumber of genes from expression profiles data #2=\t$Ngenes1\nNumber of common genes=\t$Ngenes_common","$PATH_OUTPUT/$logFileID.txt");
}
if($Ngenes_common < 20){
	error_message("Too few genes common between data sets N=$Ngenes_common",$logFileID);
}
my @genes_common=();

my $matrix_fileID = get_outputID(1);
open(OUT, ">$PATH_OUTPUT/$matrix_fileID.txt") or error_message("Cannot open $matrix_fileID.txt",$logFileID);
my $text="Symbols";
print OUT "Symbols";
for(my $i=0; $i<$nCol; ++$i){ print OUT "\t$headers[$i]"; } print OUT "\n";
for(my $i=0; $i<$nCol1; ++$i){ $text .= "\t$headers1[$i]"; } $text .= "\n";
foreach my $symbol (keys %hashRows){
	my ($ii,$ii1) = @{$hashRows{$symbol}};
	if($ii1<0){ next; }
	push(@genes_common,$symbol);
	print OUT $symbol;
	$text .= $symbol;
	for(my $i=0; $i<$nCol; ++$i){
		my $x = $Data[$i]->[$ii];
		print OUT "\t$x";
	}
	print OUT "\n";
	for(my $i=0; $i<$nCol1; ++$i){
		my $x = $Data1[$i]->[$ii1];
		$text .= "\t$x";
	}
	$text .= "\n";
}
print OUT "\n$text";
close OUT;

# Start writing the output file
if(!$outputID || !open(OUT,">$PATH_OUTPUT/$outputID.txt")){
	error_message("Cannot open output file",$logFileID);
}
my $temp_name = $file_matrix;
$temp_name =~ s/^public-//;
my $temp_name1 = $file_matrix1;
$temp_name1 =~ s/^public-//;
print OUT "!Output_name\tCorrelation of $temp_name and $temp_name1\n";
print OUT "!Output_data_set1\t$file_matrix\n";
print OUT "!Output_data_set2\t$file_matrix1\n";
print OUT "!Output_data_type1\tmatrix\n";
print OUT "!Output_data_type2\tmatrix\n";
print OUT "!Output_data_description1\t$description_matrix\n";
print OUT "!Output_data_description2\t$description_matrix1\n";
print OUT "!Output_taxid\t$organismID\n";
if($organismID1 != $organismID){
	print OUT "!Output_taxid2\t$organismID1\n";
}
print OUT "!Output_method\tCorrelation of expression profiles ($correlation_list[$correlation_method-1])\n";
print OUT "!Output_FDR1\t$FDR_thresh\n";
print OUT "!Output_FDR2\t$FDR_thresh1\n";
print OUT "!Output_fold_change1\t$fold_change_thresh\n";
print OUT "!Output_fold_change2\t$fold_change_thresh1\n";
print OUT "!Output_expr_thresh1\t$expr_thresh\n";
print OUT "!Output_expr_thresh2\t$expr_thresh1\n";
print OUT "!Output_correlation_method\t$correlation_list[$correlation_method-1]\n";
my $epfp = $hashInput{"coregulated_epfp"};
if($option_genes){
	print OUT "!Output_method_epfp\t$epfp\n";
}
print OUT "!Output_baseline1\t$baseline_value\n";
print OUT "!Output_baseline2\t$baseline_value1\n";
print OUT "!Output_direction_change\t$direction\n";
if($restrict){
	if($restrict==1){ print OUT "!Output_restrict_values\tx>0\n"; }
	elsif($restrict==2){ print OUT "!Output_restrict_values\tx>0 & y>0\n"; }
	elsif($restrict==3){ print OUT "!Output_restrict_values\tx<0\n"; }
	elsif($restrict==4){ print OUT "!Output_restrict_values\tx<0 & y<0\n"; }
}
my $autocorrelation=0;
if($file_matrix && $file_matrix eq $file_matrix1 && $FDR_thresh==$FDR_thresh1 && $fold_change_thresh==$fold_change_thresh1 && $baseline_value eq $baseline_value1 && $direction>0){
	if($restrict%2==0){
		$autocorrelation=1;
	}
}
print OUT "!Output_autocorrelation\t$autocorrelation\n";
print OUT "!Output_N_columns\t$nCol\n";
print OUT "!Output_N_rows\t$nCol1\n";
print OUT "!Output_Ngenes_file1\t$Ngenes\n";
print OUT "!Output_Ngenes_file2\t$Ngenes1\n";
print OUT "!Output_Ngenes_common\t$Ngenes_common\n";
print OUT "!Output_N_matrixes\t";
if($option_genes){ print OUT "4\n"; }
else{ print OUT "2\n"; }
print OUT "!Output_matrix1_name\tz-value\n";
if($correlation_method<3){
	print OUT "!Output_matrix2_name\tcorrelation\n";
}else{
	print OUT "!Output_matrix2_name\tcovariance\n";
}
if($option_genes){
	print OUT "!Output_matrix3_name\tgenes_upregulated\n";
	print OUT "!Output_matrix4_name\tgenes_downregulated\n";
}
print OUT "!Output_matrix1_type\tnumbers\n";
print OUT "!Output_matrix2_type\tnumbers\n";
if($option_genes){
	print OUT "!Output_matrix3_type\tgenes\tdatabase\n";
	print OUT "!Output_matrix4_type\tgenes\tdatabase\n";
	print OUT "!Output_database_headers\tgenes\tEPFP\n";
}
close OUT;
if($logFileID){
	file_append("Generating the correlation matrix..","$PATH_OUTPUT/$logFileID.txt");
}
my @correl_option = ("P","S","C");
my @command = ("$PATH_BIN/correlation_exatlas","$PATH_OUTPUT/$matrix_fileID.txt","$PATH_OUTPUT/$outputID.txt");
push(@command,$nCol,$nCol1,$Ngenes_common,"$correl_option[$correlation_method-1]");
if($logFileID){
	push(@command,"-l","$PATH_OUTPUT/$logFileID.txt");
}
if($autocorrelation){ push(@command,"-auto"); }
if($option_genes){ push(@command,"-genes","-epfp","$epfp"); }
#print "A1 @command<br>\n";
my $resp = system(@command);
if($resp){ error_message("correlation_exatlas crashed",$logFileID); }
if(!$logFileID){
	output_explore($webPageID);
}
return;
}

#**************************************
sub  geneset_analysis
#**************************************
{
my $organismID = $hashInput{"organismID"};
my $organismID1 = $hashInput{"organismID1"};
if(!$organismID1){ $organismID1=$organismID; }
my @output_list;
my @file_list;
my @geneset_list;
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Config file not found!");
while(my $line = <INFO>){
	chop $line;
	my @items = split(/[=\t]/,$line);
	if($items[0] =~ /^type_output/){
		push(@output_list,$items[1]);
	}elsif($items[0] =~ /^type_/){
		push(@file_list,$items[1]);
	}
	my %hash=();
	read_config_line($line,\%hash);
	my $file_geneset = $hash{"type_geneset"};
	if($file_geneset){
		push(@geneset_list,[$file_geneset,$hash{"description"},$hash{"organismID"},$hash{"date"}]);
	}
}
close INFO;
@geneset_list = sort {lc($a->[0]) cmp lc($b->[0])} @geneset_list;
if($loginname ne "public" && open(INFO,"<$PATH_INFO/public-config.txt")){
	my @geneset_list1=();
	while(my $line = <INFO>){
		chop $line;
		my %hash=();
		read_config_line($line,\%hash);
		my $file_geneset = $hash{"type_geneset"};
		if($file_geneset){
			push(@geneset_list1,["public-".$file_geneset,$hash{"description"},$hash{"organismID"}]);
		}
	}
	close INFO;
	push(@geneset_list, sort {lc($a->[0]) cmp lc($b->[0])} @geneset_list1);
}
my %hashOrganismID;
foreach my $ref (@geneset_list){
	$hashOrganismID{$ref->[2]}=1;
}
if($organismID1==$organismID && !$hashOrganismID{$organismID} && %hashOrganismID){
	if($hashOrganismID{9606}){ $organismID1=9606; }
	else{
		my @items = keys %hashOrganismID;
		$organismID1=$items[0];
	}
}
filter_list_by_organism(\@geneset_list, $organismID1);

my ($items,$descriptions) = get_array_lists(\@geneset_list);
my $file_list="";
my $output_list="";
foreach my $name (@output_list){
	if(!$output_list){ $output_list = "\"".$name."\""; }
	else{ $output_list .= ",\"".$name."\""; }
}
foreach my $name (@file_list){
	if(!$file_list){ $file_list = "\"".$name."\""; }
	else{ $file_list .= ",\"".$name."\""; }
}
my $file_matrix = $hashInput{"file_matrix"};
my $file_matrix_full = $file_matrix;
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_matrix_full = "$loginname-$file_matrix";
	$file_anova = "$loginname-anova-$file_matrix";
}
my $description_matrix = $hashInput{"description_matrix"};
my $file_geneset = $hashInput{"file_geneset"};
my $description_geneset = $hashInput{"description_geneset"};
my @headers = get_anova_headers($file_anova,1);
if(!@headers){ error_message("Headers not extracted"); }
my @anova_headers_file = ("Median profile",@headers);

my $logFileID = get_outputID(1);
# Print page header
print "<HTML><HEAD><TITLE>ExAtlas - geneset analysis</TITLE>\n";
print_header("update_description();");
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "output_list = new Array($output_list);\n";
print "file_list = new Array($file_list);\n";
print "geneset_list = new Array($items);\n";
print "geneset_description = new Array($descriptions);\n";
print "function alert_onsubmit() {\n";
print "	document.form_geneset.action.value = \"geneset1\";\n";
print "	document.form_geneset.submit();\n";
print "}\n";
print "function change_organism() {\n";
print "	document.form_geneset.action.value = \"geneset\";\n";
print "	document.form_geneset.submit();\n";
print "}\n";
print "function update_description() {\n";
print "	var index;\n";
print "	index = document.form_geneset.file_geneset.selectedIndex;\n";
print "	document.form_geneset.description_geneset.value = geneset_description[index];\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print "<FORM NAME=form_geneset ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<H2>Geneset enrichment analysis</H2>\n";
print "<b>Method:</b> Modified <a href=../exatlas_help.html#geneset_enrich>PAGE</a>. Modification: applied separately to upregulated and downregulated genes<p>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD><TD><b>File name<TD width=15><TD><center><b>File type<TD width=15><TD><center><b>Organism<TD width=15><TD WIDTH=300><b>Description\n";
print "<TR><TD WIDTH=60><b>File 1:<TD>$file_matrix<TD><TD><center><font size=-1>Expression profiles</font>\n";
print "<TD><TD><font size=-1>$hashOrganism{$organismID}</font><TD><TD><font size=-1>$description_matrix</font>\n";
print "<TR><TD><TD><b>Select geneset file</b>\n";
print "<TR><TD><b>File 2:\n";
print "<TD><select name=file_geneset style=\"width:250px;\" onChange=update_description();>\n";
for(my $i=0; $i<@geneset_list; ++$i){ 
	print "<option value=\"$geneset_list[$i]->[0]\"";
	if($geneset_list[$i]->[0] =~ /^public-GO/){ print " selected"; }
	print "> $geneset_list[$i]->[0]\n";
}
print "</select><TD><TD><center><font size=-1>Gene set file</font><TD><TD><center><font size=-1>$hashOrganism{$organismID1}</font>\n";
print "<TD><TD><INPUT NAME=description_geneset style=width:300px;>\n";
print "<TR><TD><TD COLSPAN=3><INPUT TYPE=CHECKBOX NAME=use_attribute> Use gene attributes (if available)\n";
print "</TABLE><br>\n";
if(keys %hashOrganismID > 1){
	print "Change organism (for file 2):\n";
	print "<select name=organismID1 style=width:200px; onChange=change_organism();>\n";
	foreach my $ref (@organisms){
		if(!$hashOrganismID{$ref->[0]}){ next; }
		print "<option value=$ref->[0]";
		if($ref->[0] == $organismID1){ print " selected"; }
		print ">$ref->[3] ($ref->[2])\n";
	}
	print "</select> <font size=-1>(page is reloaded after change)</font><p>\n";
}else{
	print "<INPUT NAME=organismID1 TYPE=hidden VALUE=$organismID1><p>\n";
}
print "<TABLE BORDER=0>\n";
print "<TR><TD>Select a baseline\n";
print "<TD><select name=baseline style=\"width: 200px;\">\n";
for(my $i=0; $i<@anova_headers_file; ++$i){ 
	print "<option value=$i> $anova_headers_file[$i]\n";
}
print "</select>\n";
print "<TD WIDTH=15><TD>Gene expression change is measured relative to baseline\n";
my @EPFP_list = (0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6);
my @fold_list = (1,1.05,1.1,1.25,1.5,2,3,4,5,7,10);
my @cutoff_list = (0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5);
my @FDR_list = (1,0.9,0.5,0.2,0.1,0.05,0.01,0.001,0.0001);
my @filter_expr = (0,0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,100,300,1000,3000);
print "<TR><TD>Cutoff proportion\n";
print "<TD><select name=target_cutoff style=\"width: 200px;\">\n";
for(my $i=0; $i<@cutoff_list; ++$i){
	my $x = $cutoff_list[$i];
	print "<option value=$x"; if($x==0.25){ print " selected"; } print "> $x\n";
}
print "</select>\n";
print "<TD><TD>Proportion of genes analyzed (from 0 to 0.5)\n";
print "<TR><TD>FDR threshold\n";
print "<TD><select name=FDR_threshold style=\"width: 200px;\">\n";
foreach my $x (@FDR_list){ print "<option value=$x> $x\n"; }
print "</select>\n";
print "<TD><TD>Used to filter genes for analysis\n";
print "<TR><TD>Fold change (global)\n";
print "<TD><select name=fold_change_threshold style=\"width: 200px;\">\n";
foreach my $x (@fold_list){ print "<option value=$x> $x\n"; }
print "</select>\n";
print "<TD><TD>Used to filter genes for analysis\n";
print "<TR><TD>Expression threshold\n";
print "<TD><select name=expr_threshold style=width:200px;>\n";
foreach my $x (@filter_expr){ print "<option value=$x> $x\n"; }
print "</select>\n";
print "<TD><TD>Used to filter genes for analysis\n";

print "<TR><TD>&nbsp;<TR><TD>Identify associated genes\n";
print "<TD><INPUT NAME=option_genes TYPE=checkbox>\n";
print "<TR><TD><a href=../exatlas-help.html#epfp>EPFP</a>\n";
print "<TD><select name=target_epfp style=\"width: 200px;\">\n";
for(my $i=0; $i<@EPFP_list; ++$i){
	my $x = $EPFP_list[$i];
	print "<option value=$x"; if($x==0.3){ print " selected"; } print "> $x\n";
}
print "</select>\n";
print "<TD><TD>EPFP threshold for associated genes\n";
print "<TR><TD>Fold change in associated genes\n";
print "<TD><select name=target_fold_change style=\"width: 200px;\">\n";
for(my $i=0; $i<@fold_list; ++$i){
	my $x = $fold_list[$i];
	print "<option value=$x"; if($x==1.5){ print " selected"; } print "> $x\n";
}
print "</select>\n";
print "<TD><TD>Minimum fold change in associated genes\n";
print "</TABLE><p>\n";
print "(*) Proportion of genes analyzed in each direction (upregulated and downregulated).<p>\n";
print "<INPUT TYPE=button value=\"Geneset analysis\" onClick=alert_onsubmit(); style=\"width: 200px;\">\n";
print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=action TYPE=hidden VALUE=\"geneset1\">\n";
print "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
print "<INPUT NAME=description_matrix TYPE=hidden VALUE=\"$description_matrix\">\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT NAME=logFileID TYPE=hidden VALUE=\"$logFileID\">\n";
print "</FORM><p>\n";
print "<HR NOSHADE></HR>\n";
print "<INPUT TYPE=button VALUE=\"Cancel (close window)\" LANGUAGE=\"javascript\" onClick=\"window.close();\" style=\"width: 200px;\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  geneset_analysis1
#**************************************
{
my $file_matrix = $hashInput{"file_matrix"};
my $file_geneset = $hashInput{"file_geneset"};
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_anova = "$loginname-anova-$file_matrix";
}
my $file_geneset_full = $file_geneset;
if($file_geneset !~ /^public-/){ $file_geneset_full = "$loginname-$file_geneset"; }
open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file");
my $line = <INFO>;
chop $line;
close INFO;
my ($junk,$junk1,@headers) =split(/\t/,$line);
if(!@headers){ error_message("Invalid ANOVA file"); }
my $nCol=0;
while($nCol<@headers && $headers[$nCol] !~ /^Var/){ ++$nCol; }
my %hashFile;
parse_geneset_file("$PATH_DATA/$file_geneset_full", \%hashFile);
my $ref = $hashFile{"geneset_names"};
my $nRow=0; 
if($ref && ref($ref) eq 'ARRAY'){
	$nRow = @$ref;
}
if($nCol<=0 || $nRow<=0){ error_message("Problem in input files: Ncol=$nCol, Nrow=$nRow"); }
my $logFileID = get_outputID(3);  #0-log,1-web,2-output
$hashInput{"logFileID"} = $logFileID;
$hashInput{"runID"} = $RUN_PAGE;
my $scale=10000;
if($hashInput{"use_attribute"} eq "on"){ $scale=100; }
interrupt_program($nCol*$nRow/$scale);
exit(0);
}

#**************************************
sub  geneset_analysis2
#**************************************
{
my $logFileID = shift;
my $webPageID = $hashInput{"logFileID"}+1;
my $outputID = $webPageID+1;

my $log10 = log(10);
my $organismID = $hashInput{"organismID"};
my $organismID1 = $hashInput{"organismID1"};
if(!$organismID1){ $organismID1=$organismID; }
my $file_matrix = $hashInput{"file_matrix"};
my $file_geneset = $hashInput{"file_geneset"};
my $file_matrix_full = $file_matrix;
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_matrix_full = "$loginname-$file_matrix";
	$file_anova = "$loginname-anova-$file_matrix";
}
my $file_geneset_full = $file_geneset;
if($file_geneset !~ /^public-/){ $file_geneset_full = "$loginname-$file_geneset"; }
my $baseline = $hashInput{"baseline"};
my $fold_change_threshold = $hashInput{"fold_change_threshold"};
my $logratio_threshold = 0;
if($fold_change_threshold>1){ $logratio_threshold =log($fold_change_threshold)/$log10; }
my $FDR_threshold = $hashInput{"FDR_threshold"};
my $expr_thresh = $hashInput{"expr_threshold"};
my $expr_log_thresh=$MISSING; if($expr_thresh>0){ $expr_log_thresh=log($expr_thresh)/$log10; }
my $use_attribute=0;
if($hashInput{"use_attribute"} eq "on"){
	$use_attribute=1;
}
open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file",$logFileID);
my $line = <INFO>;
chop $line;
my ($junk,$junk1,@headers) =split(/\t/,$line);
if(!@headers){ error_message("Invalid ANOVA file",$logFileID); }
my $nCol=0;
my $nRep=0;
while($nCol<@headers && $headers[$nCol] !~ /^Var\(/){ ++$nCol; }
for(my $i=0; $i<$nCol; ++$i){
	if($headers[$i]=~/ \(\d+\)$/ && $headers[$i] !~ /^Mean/){
		my @items = split(/ \(/,$headers[$i]);
		my $n = pop(@items);
		my $name = join(" (",@items);
		$n =~ s/\)$//;
		$nRep += $n;
		$headers[$i] = $name;
	}else{
		$headers[$i] =~ s/^Mean\(//;
		$headers[$i] =~ s/\)$//;
		$nRep++;
	}
}
splice(@headers,$nCol);
my $irow=0;
my @response;
my %symbols;
my %symbolsUC;
my @symbols;
my %hashGenes;
my %hashAlias;
get_official_symbols($organismID,\%hashGenes,\%hashAlias);
while(my $line = <INFO>){
	chop $line;
	my ($junk,$expr,@data1) =split(/\t/,$line);
	my $symbol = $data1[$nCol+8];
	if(!$symbol){ next; }
	if($expr_thresh>0 && $expr<$expr_log_thresh){ next; }
	if(!$hashGenes{$symbol}){
		my $symbol1 = $hashAlias{$symbol};
		if($symbol1){ $symbol=$symbol1; }
	}
	my $F = $data1[$nCol+4];
	my $FDR = $data1[$nCol+6];
	if($FDR_threshold>0 && $FDR > $FDR_threshold){
		next;
	}
	splice(@data1,$nCol);
	if($logratio_threshold>0){
		my @sorted = sort {$a<=>$b} @data1;
		while($sorted[0]<=$MISSING){ shift(@sorted); }
		if(@sorted>1){
			my $logratio = $sorted[@sorted-1] - $sorted[0];
			if($logratio < $logratio_threshold){ next; }
		}
	}
	my $median;
	if($baseline==0){
		$median = median(\@data1);
	}else{
		$median = $data1[$baseline-1];
		if($median==$MISSING){ $median = median(\@data1); }
	}
	my $ref = $symbols{$symbol};
	if($ref){
		if($F > $ref->[1]){
			my $iii = $ref->[0];
			for(my $i=0; $i<$nCol; ++$i){
				my $x = $data1[$i];
				if($x==$MISSING){ $response[$iii]->[$i]=$MISSING; }
				else{ $response[$iii]->[$i] = floor(10000*($data1[$i]-$median)+0.5)/10000; }
			}
			$ref->[1] = $F;
		}
	}else{
		for(my $i=0; $i<$nCol; ++$i){
			my $x = $data1[$i];
			if($x==$MISSING){ $response[$irow]->[$i]=$MISSING; }
			else{ $response[$irow]->[$i] = floor(10000*($data1[$i]-$median)+0.5)/10000; }
		}
		$symbols{$symbol} = [$irow,$F];
		my $symbol1 = uc($symbol);
		$symbolsUC{$symbol1} = $symbol;
		push(@symbols,$symbol);
		$irow++;
	}
}
close INFO;
my $nRow = $irow;
my $fileID = get_outputID(1);
open(OUT,">$PATH_OUTPUT/$fileID.txt") or error_message("Cannot open temp file",$logFileID);
print OUT "ProbeID\t".join("\t",@headers)."\n";
for(my $irow=0; $irow<$nRow; $irow++){
	print OUT "$symbols[$irow]\t".join("\t",@{$response[$irow]})."\n";
}
print OUT "\n";
my %geneHomolog=();
my %hashGenes1;
my %hashAlias1;
if($organismID1 != $organismID){
	get_gene_homolog($organismID1,$organismID,\%geneHomolog,$logFileID);
	get_official_symbols($organismID1,\%hashGenes1,\%hashAlias1);
}
my @geneSetList;
open(INFO,"<$PATH_DATA/$file_geneset_full") or error_message("Cannot open geneset file",$logFileID);
my @geneset_data;
my $name_old;
my @data;
while(my $line = <INFO>){
	if($line =~ /^!|^#/){ next; }
	chop $line;
	if(!$line){ next; }
	my($setName,$descrip,@symbol_list)=split(/\t/, $line);
	if($setName){
		my $name1 = $setName;
		my $n = @symbol_list;
		my $descripTerm = quotemeta($descrip);
		if($descrip && $name1 !~ /^$descripTerm$/i && $descrip !~ /^http|^html/){ $name1 .= ", $descrip"; }
		$descrip = $name1;
		if($name_old){
			push(@geneset_data,[@data]);
			@data=();
		}
		$name_old = $setName;
	}
	unshift(@symbol_list,$descrip);
	push(@data,\@symbol_list);
}
close INFO;
if($name_old){
	push(@geneset_data,[@data]);
}
my @original_order;
my $attrib_header;
my $reverse = 0;
for(my $igeneset=0; $igeneset<@geneset_data; $igeneset++){
	my $ref = $geneset_data[$igeneset];
	if(ref($ref) ne "ARRAY"){ error_message("Ref1 in geneset_analysis1"); }
	my ($symbols_ref,@attrib) = @{$geneset_data[$igeneset]};
	if(ref($symbols_ref) ne "ARRAY"){ error_message("Ref2 in geneset_analysis1"); }
	if(@$symbols_ref-1 < 5){ next; }
	my $sort_by = -1;
	for(my $i1=0; $i1<@attrib && $use_attribute; $i1++){
		my $header = $attrib[$i1]->[0];
		if($attrib_header){
			if(uc($header) eq uc($attrib_header)){
				$sort_by = $i1;
				last;
			}
		}elsif($header){
			if($i1==0){
				$sort_by=0;
				if($header =~ /^(FDR|EPFP|p|p-value)$/i){ $reverse = 1; }
			}
			if($sort_by==0 && $header =~ /^logratio|^log-ratio|^ClipRPM/i){
				$sort_by = $i1;
				$reverse = 0;
			}
			$attrib_header = $header;
		}
	}
	if($sort_by < 0 && $use_attribute){
		$use_attribute = 0;
	}
	if(!$attrib_header && $sort_by>=0){
		$attrib_header = $attrib[$sort_by]->[0];
	}
	my %hash;
	for(my $i=1; $i<@$symbols_ref; $i++){
		my $symbol = $symbols_ref->[$i];
		if(!$symbol){ next; }
		if($organismID1==$organismID){
			if(!$hashGenes{$symbol}){
				my $symbol1 = $hashAlias{$symbol};
				if($symbol1){ $symbol=$symbol1; }
			}
		}else{
			if(!$hashGenes1{$symbol}){
				my $symbol1 = $hashAlias1{$symbol};
				if($symbol1){ $symbol=$symbol1; }
			}
			my $symbol1 = $geneHomolog{$symbol};
			if($symbol1){
				$symbol = $symbol1;
			}else{
				$symbol1 = $symbolsUC{uc($symbol)};
				if($symbol1){
					$symbol = $symbol1;
				}else{
					next;
				}
			}
		}
		my $attrib1 = 1;
		if($sort_by>=0){ $attrib1 = $attrib[$sort_by]->[$i]; }
		my $attrib0 = $hash{$symbol};
		if(!defined($attrib0) || $reverse==0 && $attrib1>$attrib0 || $reverse==1 && $attrib1<$attrib0){
			$hash{$symbol}=$attrib1;
		}
	}
	my @symbol_list;
	if($reverse){ @symbol_list = sort {$hash{$a}<=>$hash{$b}} keys %hash; }
	else{ @symbol_list = sort {$hash{$b}<=>$hash{$a}} keys %hash; }
	my @row_numbers;
	for(my $i=0; $i<@symbol_list; ++$i){
		my $symbol = $symbol_list[$i];
		my $ref = $symbols{$symbol};
		if($ref){
			push(@row_numbers,$ref->[0]);
		}
	}
	if(@row_numbers >= 5){
my $nn1 = @row_numbers;
		print OUT "$symbols_ref->[0]\t".join("\t",@row_numbers)."\n";
		push(@original_order,$igeneset+1);
	}
}
my $nGeneSet = @original_order;
close OUT;
if(!$outputID || !open(OUT,">$PATH_OUTPUT/$outputID.txt")){
	error_message("Cannot open output file",$logFileID);
}
my $description_matrix = $hashInput{"description_matrix"};
my $description_geneset = $hashInput{"description_geneset"};
my $option_genes = 0;
if($hashInput{"option_genes"} eq "on"){ $option_genes = 1; }
my $EPFP_thresh = $hashInput{"target_epfp"};
my $fold_change_thresh = $hashInput{"target_fold_change"};
my $cutoff = $hashInput{"target_cutoff"};
my $baseline_value;
if($baseline==0){
	$baseline_value = "Median profile";
}else{
	$baseline_value = $headers[$baseline-1];
}
my $temp_name = $file_matrix;
$temp_name =~ s/^public-//;
my $temp_name1 = $file_geneset;
$temp_name1 =~ s/^public-//;
print OUT "!Output_name\tGeneset analysis of $temp_name vs. $temp_name1\n";
print OUT "!Output_data_set1\t$file_matrix\n";
print OUT "!Output_data_set2\t$file_geneset\n";
print OUT "!Output_data_type1\tmatrix\n";
print OUT "!Output_data_type2\tgeneset\n";
print OUT "!Output_data_description1\t$description_matrix\n";
print OUT "!Output_data_description2\t$description_geneset\n";
print OUT "!Output_Ngenes_file1\t$nRow\n";
print OUT "!Output_method\tParametric analysis of gene enrichment (PAGE)\n";
print OUT "!Output_method_cutoff\t$cutoff\n";
print OUT "!Output_fdr1\t$FDR_threshold\n";
print OUT "!output_fold_change1\t$fold_change_threshold\n";
print OUT "!Output_expr_thresh1\t$expr_thresh\n";
if($option_genes){
	print OUT "!Output_method_foldchange\t$fold_change_thresh\n";
	print OUT "!Output_method_epfp\t$EPFP_thresh\n";
}
print OUT "!Output_original_order2\t".join("\t",@original_order)."\n";
print OUT "!Output_taxid\t$organismID\n";
if($organismID1 != $organismID){
	print OUT "!Output_taxid2\t$organismID1\n";
}
my $Nmatrix = 3;
print OUT "!Output_N_rows\t$nGeneSet\n";
print OUT "!Output_N_columns\t$nCol\n";
print OUT "!Output_baseline\t$baseline_value\n";
if($option_genes){ $Nmatrix+=2; }
print OUT "!Output_N_matrixes\t$Nmatrix\n";
print OUT "!Output_matrix1_name\tz-value_upregulated\n";
print OUT "!Output_matrix2_name\tz-value_downregulated\n";
print OUT "!Output_matrix3_name\tz-value_combined\n";
if($option_genes){
	print OUT "!Output_matrix4_name\tgenes_upregulated\n";
	print OUT "!Output_matrix5_name\tgenes_downregulated\n";
}
print OUT "!Output_matrix1_type\tnumbers\n";
print OUT "!Output_matrix2_type\tnumbers\n";
print OUT "!Output_matrix3_type\tnumbers\n";
if($option_genes){
	print OUT "!Output_matrix4_type\tgenes\tdatabase\n";
	print OUT "!Output_matrix5_type\tgenes\tdatabase\n";
	print OUT "!Output_database_headers\tgenes\tEPFP\n";
}
close OUT;
if($logFileID){
	file_append("PAGE analysis started..","$PATH_OUTPUT/$logFileID.txt");
}
my @command = ("$PATH_BIN/page_exatlas","$PATH_OUTPUT/$fileID.txt","$PATH_OUTPUT/$outputID.txt");
push(@command,"$nCol","$nRow","$nGeneSet","-cut","$cutoff");
if($use_attribute){ push(@command,"-attrib"); }
if($logFileID){ push(@command,"-l","$PATH_OUTPUT/$logFileID.txt"); }
if($option_genes){ push(@command,"-genes","-epfp","$EPFP_thresh","-fold","$fold_change_thresh"); }
#print "@command<br>\n";
my $resp = system(@command);
if($resp){
	error_message("page_analysis crashed",$logFileID);
}
if(!$logFileID){
	output_explore($webPageID);
}
return;
}

#**************************************
sub  geneset_overlap
#**************************************
{
my $file_geneset = $hashInput{"file_geneset"};
my $file_geneset1 = $hashInput{"file_geneset1"};

if(!$file_geneset1){ error_message("file_geneset1 is empty!"); }
my $file_geneset_full1 = $file_geneset1;
if($file_geneset1 !~ /^public-/){ $file_geneset_full1 = "$loginname-$file_geneset1"; }
$file_geneset_full1 = "$PATH_DATA/$file_geneset_full1";
my $file_geneset_full;
if($file_geneset){
	$file_geneset_full = $file_geneset;
	if($file_geneset !~ /^public-/){ $file_geneset_full = "$loginname-$file_geneset"; }
	$file_geneset_full = "$PATH_DATA/$file_geneset_full";
}else{
	my $file_upload = $hashInput{"upload_geneset"};
	if(!$file_upload){ error_message("Geneset not specified in geneset_overlap"); }
	my $fileID = get_outputID(1);
	$file_geneset = "TEMPORARY!-$fileID.txt";
	$hashInput{"file_geneset"} = $file_geneset;
	open (INFO,"<$PATH_OUTPUT/$file_upload") or error_message("Cannot open $file_upload");
	my @headers;
	my $is=0;
	my @symbols;
	my ($title,$description);
	my $count=1;
	open (OUT,">$PATH_OUTPUT/$fileID.txt") or error_message("Cannot open temp file");
	while(my $line=<INFO>){
		chop $line;
		if(!$line){ next; }
		if($line=~/^!/){
			if(@headers){
				if(!$description && $hashInput{"description_geneset"}){
					$description = $hashInput{"description_geneset"};
				}
				if(!$hashInput{"description_geneset"}){
					$hashInput{"description_geneset"} = $description;
				}
				if(!$title){ $title = "Temporary-$count"; }
				print OUT "$title\t$description\t".join("\t",@symbols)."\n";
				@headers=();
				@symbols=();
				($title,$description)=("","");
				$count++;
			}
			my ($key,$value) = split(/\t/,$line);
			if($key =~ /title/){ $title = $value; }
			if($key =~ /description/){ $description = $value; }
			next;
		}
		if(!@headers){
			@headers=split(/\t/,$line);
			while($is<@headers && $headers[$is] !~ /^Gene symbol|^Symbol/i){ $is++; }
			if($is==@headers){ error_message("Gene symbols not fount in geneset_overlap"); }
			next;
		}
		my @items=split(/\t/,$line);
		if($items[$is]){
			push(@symbols,$items[$is]);
		}
	}
	close INFO;
	if(!$description && $hashInput{"description_geneset"}){
		$description = $hashInput{"description_geneset"};
	}
	if(!$hashInput{"description_geneset"}){
		$hashInput{"description_geneset"} = $description;
	}
	if(!$title){
		$title = "Temporary-$count";
	}
	print OUT "$title\t$description\t".join("\t",@symbols)."\n";
	close OUT;
	$file_geneset_full = "$PATH_OUTPUT/$fileID.txt";
}
my $nGeneset=0;
open(INFO,'<',$file_geneset_full) or error_message("Cannot open $file_geneset");
while(my $line=<INFO>){
	if($line=~/^!/){ next; }
	my @items = split(/\t/,$line);
	if($items[0]){ $nGeneset++; }
}
close INFO;
my $nGeneset1=0;
open(INFO,'<',$file_geneset_full1) or error_message("Cannot open $file_geneset1");
while(my $line=<INFO>){
	if($line=~/^!/){ next; }
	my @items = split(/\t/,$line);
	if($items[0]){ $nGeneset1++; }
}
close INFO;
my $logFileID = get_outputID(3); #0-log,1-web,2-output
$hashInput{"logFileID"} = $logFileID;
$hashInput{"runID"} = $RUN_OVERLAP;
my $scale=100000;
if($hashInput{"use_attribute"} eq "on"){ $scale=1000; }
interrupt_program($nGeneset*$nGeneset1/$scale);
exit(0);
}

#**************************************
sub  geneset_overlap1
#**************************************
{
my $logFileID = shift;
my $webPageID = $hashInput{"logFileID"}+1;
my $outputID = $webPageID+1;

my $organismID = $hashInput{"organismID"};
my $organismID1 = $hashInput{"organismID1"};
if(!$organismID1){ $organismID1=$organismID; }
my $FDR_thresh = $hashInput{"FDR"};
my $fold_enrichment_thresh = $hashInput{"fold_enrichment"};
my $minN = $hashInput{"minimum_genes"};
my $file_geneset1 = $hashInput{"file_geneset"};
my $file_geneset2 = $hashInput{"file_geneset1"};
my $description_geneset1 = $hashInput{"description_geneset"};
my $description_geneset2 = $hashInput{"description_geneset1"};
my $use_attribute=0;
if($hashInput{"use_attribute"} eq "on"){
	$use_attribute=1;
}
my $autocorrelation=0;
if($file_geneset1 && $file_geneset1 eq $file_geneset2){
	$autocorrelation=1;
}
my ($file_geneset_full1,$file_geneset_full2);
if($file_geneset1 =~ /^TEMPORARY!-/){
	my $file_temp = $file_geneset1;
	$file_temp =~ s/^TEMPORARY!-//;
	$file_geneset_full1 = "$PATH_OUTPUT/$file_temp";
}else{
	$file_geneset_full1 = "$PATH_DATA/$file_geneset1";
	if($file_geneset1 !~ /^public-/){ $file_geneset_full1 = "$PATH_DATA/$loginname-$file_geneset1"; }
}
if($file_geneset2 =~ /^TEMPORARY!-/){
	my $file_temp = $file_geneset2;
	$file_temp =~ s/^TEMPORARY!-//;
	$file_geneset_full2 = "$PATH_OUTPUT/$file_temp";
}else{
	$file_geneset_full2 = "$PATH_DATA/$file_geneset2";
	if($file_geneset2 !~ /^public-/){ $file_geneset_full2 = "$PATH_DATA/$loginname-$file_geneset2"; }
}
my %hashGenes;
my %hashAlias;
get_official_symbols($organismID,\%hashGenes,\%hashAlias);
my @geneset_data1;
my @geneset_name1;
my %allGenes;
my %symbolsUC;
my @original_order1;
my $count=0;
open(INFO,'<',$file_geneset_full1) or error_message("Cannot open $file_geneset1",$logFileID);
while(my $line=<INFO>){
	if($line =~ /^!|^#/){ next; }
	chop $line;
	my($name,$descrip,@symbols)=split(/\t/,$line);
	if(!$name){ next; }
	while(@symbols && !$symbols[0]){ shift(@symbols); }
	$count++;
	if(!$autocorrelation && @symbols < $minN){ next; }
	my %hash;
	foreach my $symbol (@symbols){
		if(!$hashGenes{$symbol}){
			my $symbol1 = $hashAlias{$symbol};
			if($symbol1){ $symbol=$symbol1; }
		}
		$hash{$symbol}=1;
		$allGenes{$symbol}=1;
		my $symbol1 = uc($symbol);
		$symbolsUC{$symbol1} = $symbol;
	}
	my $n = keys %hash;
	if(!$autocorrelation && $n < $minN){ next; }
	my $descripTerm = quotemeta($descrip);
	if($descrip && $name !~ /^$descripTerm$/i && $descrip !~ /^http|^html/){ $name .= ", $descrip"; }
	push(@geneset_name1,[$name,$n]);
	push(@geneset_data1,\%hash);
	push(@original_order1,$count);
}
close INFO;
if(!@geneset_data1){ error_message("Empty geneset file #1"); }
my %geneHomolog=();
my %hashGenes1;
my %hashAlias1;
if($organismID1 != $organismID){
	get_gene_homolog($organismID1,$organismID,\%geneHomolog,$logFileID);
	get_official_symbols($organismID1,\%hashGenes1,\%hashAlias1);
}
my @original_order2;
$count=1;
my @geneset_name2;
my @geneset_data2;
my @results;
my $name_old;
my @data;
open(INFO,'<',$file_geneset_full2) or error_message("Cannot open $file_geneset2",$logFileID);
while(my $line=<INFO>){
	if($line =~ /^!|^#/){ next; }
	chop $line;
	my($name,$descrip,@symbols)=split(/\t/,$line);
	if($name){
		my $name1 = $name;
		my $n = @symbols;
		my $descripTerm = quotemeta($descrip);
		if($descrip && $name1 !~ /^$descripTerm$/i && $descrip !~ /^http|^html/){ $name1 .= ", $descrip"; }
		push(@geneset_name2,[$name1,$n]);
		push(@original_order2,$count++);
		$descrip = "symbols";
		if($name_old){
			push(@geneset_data2,[@data]);
			@data=();
		}
		$name_old = $name;
	}
	#Load symbols and genes attributes if present: zero position=attribute title, then go values
	unshift(@symbols,$descrip);
	push(@data,\@symbols);
}
close INFO;
if($name_old){
	push(@geneset_data2,[@data]);
}
my $attrib_header;
my $reverse = 0;
my $genesetDelta = 100;
if(@geneset_data2 < 10){ $genesetDelta = 1; }
elsif(@geneset_data2 < 100){ $genesetDelta = 5; }
elsif(@geneset_data2 < 500){ $genesetDelta = 20; }
for(my $igeneset=0; $igeneset<@geneset_data2; $igeneset++){
	if($logFileID && $igeneset%$genesetDelta==0){
		file_append("Geneset #$igeneset","$PATH_OUTPUT/$logFileID.txt");
	}
	my $ref = $geneset_data2[$igeneset];
	if(ref($ref) ne "ARRAY"){ error_message("Ref1 in geneset_overlap1"); }
	my ($symbols_ref,@attrib) = @{$geneset_data2[$igeneset]};
	if(ref($symbols_ref) ne "ARRAY"){ error_message("Ref2 in geneset_overlap1"); }
	#Array @$symbols_ref starts with word "symbols" followed by gene symbols
	if(@$symbols_ref-1 < $minN){ next; }
	my $sort_by = -1;
	for(my $i1=0; $i1<@attrib && $use_attribute; $i1++){
		my $header = $attrib[$i1]->[0];
		if($attrib_header){
			if(uc($header) eq uc($attrib_header)){
				$sort_by = $i1;
				last;
			}
		}elsif($header){
			if($i1==0){
				$sort_by=0;
				if($header =~ /^(FDR|EPFP|p|p-value)$/i){ $reverse = 1; }
			}
			if($sort_by==0 && $header =~ /^logratio|^log-ratio|^ClipRPM/i){
				$sort_by = $i1;
				$reverse = 0;
			}
			$attrib_header = $header;
		}
	}
	if($sort_by < 0 && $use_attribute){
		$use_attribute = 0;
	}
	if(!$attrib_header && $sort_by>=0){
		$attrib_header = $attrib[$sort_by]->[0];
	}
	my %hash;
	for(my $i=1; $i<@$symbols_ref; $i++){
		my $symbol = $symbols_ref->[$i];
		if(!$symbol){ next; }
		if($organismID1==$organismID){
			if(!$hashGenes{$symbol}){
				my $symbol1 = $hashAlias{$symbol};
				if($symbol1){ $symbol=$symbol1; }
			}
		}else{
			if(!$hashGenes1{$symbol}){
				my $symbol1 = $hashAlias1{$symbol};
				if($symbol1){ $symbol=$symbol1; }
			}
			my $symbol1 = $geneHomolog{$symbol};
			if($symbol1){
				$symbol = $symbol1;
			}else{
				$symbol1 = uc($symbol);
				if($symbolsUC{$symbol1}){
					$symbol = $symbolsUC{$symbol1};
				}
			}
		}
		$allGenes{$symbol}=1;
		my $attrib1 = 1;
		if($sort_by>=0){ $attrib1 = $attrib[$sort_by]->[$i]; }
		my $attrib0 = $hash{$symbol};
		if(!defined($attrib0) || $reverse==0 && $attrib1>$attrib0 || $reverse==1 && $attrib1<$attrib0){
			$hash{$symbol}=$attrib1;
		}
	}
	my @symbols;
	if($reverse){ @symbols = sort {$hash{$a}<=>$hash{$b}} keys %hash; }
	else{ @symbols = sort {$hash{$b}<=>$hash{$a}} keys %hash; }
	my $n = @symbols;
	if($n < $minN){ next; }
	$geneset_name2[$igeneset]->[1] = $n;
	for(my $j=0; $j<@geneset_name1; $j++){
		if($autocorrelation && $igeneset==$j){ 
			next;
		}
		my @overlap=();
		my $attrib_thresh = $hash{$symbols[@symbols-1]};
		for(my $i=0; $i<@symbols; $i++){
			my $symbol = $symbols[$i];
			if($geneset_data1[$j]->{$symbol}){
				push(@overlap,$symbol,$i,$hash{$symbol});
			}
		}
		if(@overlap){
			push(@results,[0,$igeneset,$j,1,0,$attrib_thresh,\@overlap]);
		}
	}
}
close INFO;

my $Ngenes = keys %allGenes;
my $numOrganism = $hashOrganismNum{$organismID};
if($numOrganism && $Ngenes < $organisms[$numOrganism-1]->[1]){ $Ngenes = $organisms[$numOrganism-1]->[1]; }
for(my $ii=0; $ii<@results; $ii++){
	my $ref = $results[$ii];
	my ($z,$i2,$i1,$FDR,$ratio,$attrib_thresh,$overlap_ref) = @$ref;
	my $N1 = $geneset_name2[$i2]->[1];
	my $m = int(@$overlap_ref/3);
	my $p = $geneset_name1[$i1]->[1]/$Ngenes;
	if($p==1){ $p=0.999999; }
	if($p==0){ $p=0.000001; }
	my $q = 0;
	if($N1>0){ $q = $m/$N1; }
	$ratio = $q/$p;
	if($N1>=$minN && $m>=$minN){
		#use hypergeometric distribution test
		if($N1>0 && $Ngenes>1){
			$z = ($q-$p)/sqrt($p*(1-$p)*($Ngenes-$N1)/($Ngenes-1)/$N1);
		}
	}
	if($use_attribute && $attrib_header && $m>1 && $overlap_ref->[2] != $overlap_ref->[$m*3-1]){
		#print "B2 ". join("; ",@$overlap_ref)."<br>\n";
		for(my $i=$m-2; $i>=0; $i--){
			if($i+1<$minN){ last; }
			my $N1a = $overlap_ref->[$i*3+1];
			if($i<$m-2){ $N1a = ($N1a + $overlap_ref->[($i+1)*3+1])/2; }
			my $qa = 0;
			if($N1a>0){ $qa = ($i+1)/$N1a; }
			my $ratioa = (1,$qa/$p);
			if($ratioa < $fold_enrichment_thresh){ next; }
			my $za =0;
			if($N1a>0 && $Ngenes>1){
				$za = ($qa-$p)/sqrt($p*(1-$p)*($Ngenes-$N1a)/($Ngenes-1)/$N1a);
			}
			if($za > $z){
				#print "B1 $z $i $za $ratio $ratioa<br>\n";
				$z = $za;
				$ratio = $ratioa;
				$attrib_thresh = $overlap_ref->[$i*3+2];
			}
		}
	}
	$z = floor(10000*$z+0.5)/10000;
	if($z<0){ $z=0; }
	$ratio = floor(10000*$ratio+0.5)/10000;
	$ref->[0] = $z;
	$ref->[4] = $ratio;
	$ref->[5] = $attrib_thresh;
}
my @sorted = sort {$a->[0]<=>$b->[0]} @results;
my $N = @sorted;
my $rank = $N;
my $FDR1=1;
my @rows_output;
my @cols_output;
for(my $i=0; $i<@sorted; ++$i){
	my $ref = $sorted[$i];
	my $z = $ref->[0];
	my $p = 1;
	if($z>0){
		$p = 2*(1 - normal_distribution($z));
	}
	my $FDR = int(10000*$p*$N/$rank+0.5)/10000;
	if($FDR > $FDR1){ $FDR = $FDR1; }
	else{ $FDR1 = $FDR; }
	$ref->[3] = $FDR;
	my $ratio = $ref->[4];
	if($FDR <= $FDR_thresh && $ratio >= $fold_enrichment_thresh){
		$rows_output[$ref->[1]]=1;
		$cols_output[$ref->[2]]=1;
	}
}
my $nCol = 0;
my $nRow = 0;
my @original_order1a;
for(my $i=0; $i<@cols_output; ++$i){
	if($cols_output[$i]){
		push(@original_order1a,$original_order1[$i]);
		$nCol++;
	}
}
my @original_order2a;
for(my $i=0; $i<@rows_output; ++$i){
	if($rows_output[$i]){
		push(@original_order2a,$original_order2[$i]);
		$nRow++;
	}
}
if(!$nCol || !$nRow){
	my $text = "<H3>No significant operlap between gene sets</H3>Consider using less stringent criteria. For example, set Fold enrichment threshold = 0.001 and FDR = 1";
	if($logFileID){ error_message($text,$logFileID); }
	terminal_window($text);
}

# Write the output file
if(!$outputID || !open(OUT,">$PATH_OUTPUT/$outputID.txt")){
	error_message("Cannot open output file",$logFileID);
}
my $temp_name = $file_geneset1;
$temp_name =~ s/^public-//;
my $temp_name1 = $file_geneset2;
$temp_name1 =~ s/^public-//;
print OUT "!Output_name\tGeneset overlap of $temp_name and $temp_name1\n";
print OUT "!Output_description\tGeneset overlap of $temp_name and $temp_name1\n";
print OUT "!Output_data_set1\t$file_geneset1\n";
print OUT "!Output_data_set2\t$file_geneset2\n";
print OUT "!Output_data_description1\t$description_geneset1\n";
print OUT "!Output_data_description2\t$description_geneset2\n";
print OUT "!Output_data_type1\tgeneset\n";
print OUT "!Output_data_type2\tgeneset\n";
print OUT "!Output_taxid\t$organismID\n";
if($organismID1 && $organismID1!=$organismID){
	print OUT "!Output_taxid2\t$organismID1\n";
}
print OUT "!Output_method\tOverlap of gene sets, hypergeomentric test\n";
print OUT "!Output_autocorrelation\t$autocorrelation\n";
print OUT "!Output_original_order1\t".join("\t",@original_order1a)."\n";
print OUT "!Output_original_order2\t".join("\t",@original_order2a)."\n";
print OUT "!Output_N_rows\t$nRow\n";
print OUT "!Output_N_columns\t$nCol\n";
print OUT "!Output_FDR\t$FDR_thresh\n";
print OUT "!Output_fold_enrichment\t$fold_enrichment_thresh\n";
print OUT "!Output_minimum_overlap\t$minN\n";
my $nMatrix = 4;
if($use_attribute){
	print OUT "!Output_geneset_attribute_name\t$attrib_header\n";
	print OUT "!Output_geneset_attribute_reverse\t$reverse\n";
	$nMatrix=5;
}
print OUT "!Output_N_matrixes\t$nMatrix\n";
print OUT "!Output_matrix1_name\tz-value\n";
print OUT "!Output_matrix2_name\tfold_enrichment\n";
print OUT "!Output_matrix3_name\tFDR\n";
print OUT "!Output_matrix4_name\tgenes\n";
print OUT "!Output_matrix1_type\tnumbers\n";
print OUT "!Output_matrix2_type\tnumbers\n";
print OUT "!Output_matrix3_type\tnumbers\n";
print OUT "!Output_matrix4_type\tgenes\n";
if($use_attribute){
	print OUT "!Output_matrix5_name\tattribute_threshold\n";
	print OUT "!Output_matrix5_type\tnumbers\n";
}
my %hash;
for(my $ii=0; $ii<@results; ++$ii){
	my ($z,$i,$j) = @{$results[$ii]};
	$hash{"$i,$j"} = $ii;
}
print OUT "!Matrix1_start\n";
print OUT "GeneSetName";
for(my $i=0; $i<@cols_output; ++$i){
	if(!$cols_output[$i]){ next; }
	print OUT "\t$geneset_name1[$i]->[0]";
}
print OUT "\n";
for(my $j=0; $j<@rows_output; ++$j){
	if(!$rows_output[$j]){ next; }
	print OUT "$geneset_name2[$j]->[0]";
	for(my $i=0; $i<@cols_output; ++$i){
		if(!$cols_output[$i]){ next; }
		my $ii = $hash{"$j,$i"};
		my $x = 0;
		if(defined($ii)){ $x = $results[$ii]->[0]; }
		print OUT "\t$x";
	}
	print OUT "\n";
}
print OUT "!Matrix1_end\n";
print OUT "!Matrix2_start\n";
print OUT "GeneSetName";
for(my $i=0; $i<@cols_output; ++$i){
	if(!$cols_output[$i]){ next; }
	print OUT "\t$geneset_name1[$i]->[0]";
}
print OUT "\n";
for(my $j=0; $j<@rows_output; ++$j){
	if(!$rows_output[$j]){ next; }
	print OUT "$geneset_name2[$j]->[0]";
	for(my $i=0; $i<@cols_output; ++$i){
		if(!$cols_output[$i]){ next; }
		my $ii = $hash{"$j,$i"};
		my $x = 1;
		if(defined($ii)){ $x = $results[$ii]->[4]; }
		print OUT "\t$x";
	}
	print OUT "\n";
}
print OUT "!Matrix2_end\n";
print OUT "!Matrix3_start\n";
print OUT "GeneSetName";
for(my $i=0; $i<@cols_output; ++$i){
	if(!$cols_output[$i]){ next; }
	print OUT "\t$geneset_name1[$i]->[0]";
}
print OUT "\n";
for(my $j=0; $j<@rows_output; ++$j){
	if(!$rows_output[$j]){ next; }
	print OUT "$geneset_name2[$j]->[0]";
	for(my $i=0; $i<@cols_output; ++$i){
		if(!$cols_output[$i]){ next; }
		my $ii = $hash{"$j,$i"};
		my $x = 1;
		if(defined($ii)){ $x = $results[$ii]->[3]; }
		print OUT "\t$x";
	}
	print OUT "\n";
}
print OUT "!Matrix3_end\n";
print OUT "!Matrix4_start\n";
print OUT "GeneSetName";
for(my $i=0; $i<@cols_output; ++$i){
	if(!$cols_output[$i]){ next; }
	print OUT "\t$geneset_name1[$i]->[0]";
}
print OUT "\n";
for(my $j=0; $j<@rows_output; ++$j){
	if(!$rows_output[$j]){ next; }
	print OUT "$geneset_name2[$j]->[0]";
	for(my $i=0; $i<@cols_output; ++$i){
		if(!$cols_output[$i]){ next; }
		my $ii = $hash{"$j,$i"};
		my $x = "";
		if(defined($ii)){
			my @symbols;
			my $attrib_thresh = $results[$ii]->[5];
			my $ref = $results[$ii]->[6];
			my $n = 0;
			if(ref($ref) eq 'ARRAY'){ $n = @{$ref}; }
			for(my $j=0; $j<$n; $j+=3){
				if($use_attribute && (!$reverse && $ref->[$j+2]<$attrib_thresh || $reverse && $ref->[$j+2]>$attrib_thresh)){ last; }
				push(@symbols,$ref->[$j]);
			}
			$x = join(',',@symbols);
		}elsif($autocorrelation && $i==$j){
			$x = join(',',sort keys %{$geneset_data1[$i]});
		}
		print OUT "\t$x";
	}
	print OUT "\n";
}
print OUT "!Matrix4_end\n";
if($use_attribute){
	print OUT "!Matrix5_start\n";
	print OUT "GeneSetName";
	for(my $i=0; $i<@cols_output; ++$i){
		if(!$cols_output[$i]){ next; }
		print OUT "\t$geneset_name1[$i]->[0]";
	}
	print OUT "\n";
	for(my $j=0; $j<@rows_output; ++$j){
		if(!$rows_output[$j]){ next; }
		print OUT "$geneset_name2[$j]->[0]";
		for(my $i=0; $i<@cols_output; ++$i){
			if(!$cols_output[$i]){ next; }
			my $ii = $hash{"$j,$i"};
			my $x = 0;
			if(defined($ii)){
				$x = $results[$ii]->[5];
			}
			print OUT "\t$x";
		}
		print OUT "\n";
	}
	print OUT "!Matrix5_end\n";
}
close OUT;
if(!$logFileID){
	output_explore($webPageID);
}
return;
}

#**************************************
sub  matrix2geneset
#**************************************
{
my @geneset_list;
my @file_list;
my $file_matrix=$hashInput{"file_matrix"};
my $description_matrix=$hashInput{"description_matrix"};
my $file_geneset = $hashInput{"file_geneset"};
my $organismID = $hashInput{"organismID"};

open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Config file not found!");
while(my $line = <INFO>){
	if(length($line) < 3) { next; }
	my @items = split(/[=\t]/,$line);
	if($items[0] =~ /^type_geneset/){
		push(@geneset_list,$items[1]);
	}elsif($items[0] =~ /^type_/){
		push(@file_list,$items[1]);
	}
}
close INFO;

my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_anova = "$loginname-anova-$file_matrix";
}
open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file");
my $line = <INFO>;
close INFO;
chop $line;
my ($junk,$junk1,@header_col) =split(/\t/,$line);
if(!@header_col){ error_message("Headers not found in $file_anova"); }
my $nLevel=0;
while($nLevel<@header_col && $header_col[$nLevel] !~ /^Var\(/){ ++$nLevel; }
for(my $i=0; $i<$nLevel; ++$i){
	$header_col[$i] =~ s/ \(\d+\)$//;
}
splice(@header_col,$nLevel);

my $file_list="";
my $geneset_list="";
foreach my $name (@geneset_list){
	if(!$geneset_list){ $geneset_list = "\"".$name."\""; }
	else{ $geneset_list .= ",\"".$name."\""; }
}
foreach my $name (@file_list){
	if(!$file_list){ $file_list = "\"".$name."\""; }
	else{ $file_list .= ",\"".$name."\""; }
}
print "<HTML><HEAD><TITLE>ExAtlas - matrix2geneset</TITLE>\n";
print_header();
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "file_list = new Array($file_list);\n";
print "geneset_list = new Array($geneset_list);\n";
print "function alert_onsubmit() {\n";
print "	var file = document.form_matrix2geneset.file_geneset.value;\n";
print "	if(!file){\n";
print "		alert(\"The geneset file name is missing. Please type it in\"); return false;\n";
print "	}\n";
print "	var file1=file;\n";
print "	if(file.search(/\\.txt\$/)>=0){\n";
print "		if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "		file1=file.substring(0,file.length-4);\n";
print "	}\n";
print "	if(file1.search(/^[-\\w]+\$/)<0){\n";
print "		alert(\"File name should be one word with no special characters except underscore and dash.\\nPlease rename it!\");\n";
print "		return(false);\n";
print "	}\n";
print "	for(i=0; i<geneset_list.length; ++i){\n";
print "		if(file == geneset_list[i] && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	if(file.search(/^public-/i) >= 0){\n";
print "		alert(\"File name cannot start with 'public-'\"); return(false);\n";
print "	}\n";
print "	for(i=0; i<file_list.length; ++i){\n";
print "		if(file == file_list[i]){\n";
print "			alert(\"A file with this name already exists\"); return false;\n";
print "		}\n";
print "	}\n";
print "	var descrip = document.form_matrix2geneset.description_geneset.value;\n";
print "	if(descrip.search(/\\=|\\&/) >= 0){\n";
print "		alert(\"Description should not include character \'=\' or \'&\'\");\n";
print "		return false;\n";
print "	}\n";
print "	document.form_matrix2geneset.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";

my @FDR_list = (1,0.5,0.2,0.1,0.05,0.01,0.001,0.0001);
my @fold_enrichment = (1,1.5,2,3,4,5,7,10);
my @specific_list = (0,2,3,4,5,6,7,10,15,20);
my @filter_expr = (0,0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,100,300,1000,3000);
my ($FDR,$fold_enrichment) = (0.05,2);
print "<H2>Extract significant genes and save them in a new geneset file</H2>\n";
print "<b>Input file:</b> $file_matrix &nbsp; &nbsp; &nbsp; <b>Description:</b> $description_matrix<br>\n";
print "<b>Organism:</b> $hashOrganism{$organismID}<br>\n";
print "<b>Note:</b> Expression can be compared with the median value (default), or with any experiment (column in matrix)<br>\n";
print "Another option is to extract genes that are both significant and specific to cells/tissues<p>\n";
print "<FORM NAME=form_matrix2geneset ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD WIDTH=170><b>Parameters\n";
print "<TR><TD>Select baseline\n";
print "<TD><select name=select_column style=width:250px;>\n";
print "<option value=0> Median profile\n";
for(my $i=0; $i<@header_col; ++$i){
	my $x=$i+1;
	print "<option value=$x> $header_col[$i]\n";
}
print "</select><TD> All columns are compared to baseline\n";
print "<TR><TD>FDR threshold<TD><select name=FDR style=width:250px;>\n";
for(my $i=0; $i<@FDR_list; ++$i){ 
	print "<option value=$FDR_list[$i]"; if($FDR==$FDR_list[$i]){ print " selected"; } print "> $FDR_list[$i]\n";
}
print "</select>\n";
print "<TR><TD>Fold change threshold<TD><select name=fold_enrichment style=width:250px;>\n";
for(my $i=0; $i<@fold_enrichment; ++$i){ 
	print "<option value=$fold_enrichment[$i]"; if($fold_enrichment==$fold_enrichment[$i]){ print " selected"; } print "> $fold_enrichment[$i]\n";
}
print "</select>\n";
print "<TR><TD>Expression threshold<TD><select name=expr_threshold style=width:250px;>\n";
foreach my $x (@filter_expr){ print "<option value=$x> $x\n"; }
print "</select>\n";
print "<TR><TD>Specific gene filter<sup>*)</sup><TD><select name=specific style=width:250px;>\n";
for(my $i=0; $i<@specific_list; ++$i){ 
	print "<option value=$specific_list[$i]> $specific_list[$i]\n";
}
print "</select><TD>z-value (0 = no filtering,7 = stringently specific)\n";
print "<TR><TD>Genes without symbol<TD><INPUT NAME=include_all_genes TYPE=checkbox><TD>Include genes with no symbol\n";
print "<TR><TD>Make output tables<TD><INPUT NAME=output_significant TYPE=checkbox><TD>Make tables (output file) with z-value, fold change, FDR, etc.\n";
print "<TR><TD>&nbsp;<TR><TD>Output geneset file:<TD><INPUT NAME=file_geneset style=width:270px;>\n";
print "<TR><TD>Description:<TD><INPUT NAME=description_geneset style=width:270px;>\n";
print "<TR><TD><INPUT TYPE=button value=\" Save significant genes \" LANGUAGE=javascript onClick=\"alert_onsubmit();\">\n";
print "</TABLE><p>\n";
print "<b>Note:</b> <sup>*)</sup>Filtering of specific genes is based on z-value, which measures the significance of difference<br>\n";
print "between gene expression in one column compared to the average in other uncorrelated columns.<p>\n";
print "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=hidden VALUE=\"matrix2geneset1\">\n";
print "<INPUT NAME=\"file_matrix\" TYPE=hidden VALUE=\"$file_matrix\">\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "</FORM>\n";
print "<HR NOSHADE></HR>\n";
print "<p><INPUT TYPE=button VALUE=\" Cancel (close window) \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  matrix2geneset1
#**************************************
{
my $file_matrix = $hashInput{"file_matrix"};
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_anova = "$loginname-anova-$file_matrix";
}
open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file");
my $line = <INFO>;
close INFO;
chop $line;
my ($junk,$junk1,@headers) =split(/\t/,$line);
if(!@headers){ error_message("Headers not found in anova file"); }
my $nCol=0;
while($nCol<@headers && $headers[$nCol] !~ /^Var\(/){ ++$nCol; }
my $logFileID = get_outputID(1);
$hashInput{"logFileID"} = $logFileID;
$hashInput{"runID"} = $RUN_SIGNIFICANT;
interrupt_program($nCol/150);
exit(0);
}

#**************************************
sub  generate_significant_genesets
#**************************************
{
my $logFileID = shift;
my $file_geneset = $hashInput{"file_geneset"};
my $description = $hashInput{"description_geneset"};
my $file_geneset_full = "$loginname-$file_geneset";
my $file_matrix = $hashInput{"file_matrix"};
my $FDR_thresh = $hashInput{"FDR"};
my $fold_thresh = $hashInput{"fold_enrichment"};
my $expr_thresh = $hashInput{"expr_threshold"};
my $file_anova = $file_matrix;
my $organismID = $hashInput{"organismID"};
my $select_column = $hashInput{"select_column"};
my $specific = $hashInput{"specific"};
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_anova = "$loginname-anova-$file_matrix";
}
open(OUT,">$PATH_DATA/$file_geneset_full") or error_message("Cannot open $file_geneset,$logFileID");
my $file_description = $hashInput{"description"};
print OUT "!Geneset_name\t$file_geneset\n";
print OUT "!Geneset_description\t$description\n";
print OUT "!Geneset_source_matrix\t$file_matrix\n";
print OUT "!Geneset_taxid\t$organismID\n";
print OUT "!Geneset_FDR\t$FDR_thresh\n";
print OUT "!Geneset_fold_enrichment\t$fold_thresh\n";
if($expr_thresh>0){ print OUT "!Geneset_expr_threshold\t$expr_thresh\n"; }
if($specific){ print OUT "!Geneset_specific\t$specific\n"; }
my @command = ("$PATH_BIN/pairwise","$PATH_DATA/$file_anova","$PATH_DATA/$file_geneset_full","$select_column");
push(@command,"-fold","$fold_thresh","-fdr","$FDR_thresh","-spec","$specific","-add");
if($expr_thresh>0){ push(@command,"-expr","$expr_thresh"); }
if($hashInput{"include_all_genes"} ne "on"){ push(@command,"-symbol"); }
if($hashInput{"output_significant"} eq "on"){
	my $fileOutputID = get_outputID(1);
	print OUT "!Geneset_output_table\t$fileOutputID\n";
	push(@command,"-out","$PATH_OUTPUT/$fileOutputID.txt","-source","$file_matrix");
}
if($logFileID){ push(@command,"-log","$PATH_OUTPUT/$logFileID.txt"); }
close OUT;
my $resp = system(@command);
if($resp){ error_message("Pairwise comparison crashed",$logFileID); }
my %hashGeneset;
parse_geneset_file("$PATH_DATA/$file_geneset_full", \%hashGeneset,$logFileID);
my $ref = $hashGeneset{"geneset_names"};
if(!$ref || ref($ref) eq 'ARRAY' && @$ref==0){
	my $text = "<H3>No significant genes found!</H3>Consider using less stringent criteria (e.g., set Fold change threshold = 1.5 and FDR threshold = 0.5)";
	if($logFileID){ error_message($text,$logFileID); }
	unlink "$PATH_DATA/$file_geneset_full";
	terminal_window($text);
}
#`chmod 666 $PATH_DATA/$file_geneset_full`;
open(OUT, ">$PATH_INFO/$loginname-config1.txt");
open(INFO,"<$PATH_INFO/$loginname-config.txt");
my $nLines;
while(my $line = <INFO>){
	$nLines++;
	if($line !~ /^type_geneset=$file_geneset\s/ && length($line)>2){
		print OUT $line;
	}
}
close INFO;
print OUT "type_geneset=$file_geneset\torganismID=$organismID";
if($description){ print OUT "\tdescription=$description"; }
print OUT "\tdate=$date_record\n";
close OUT;
my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
if($nLines1 && $nLines1 > $nLines*0.9){
	copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
}else{
	error_message("Failed to update configuration file!");
}
unlink("$PATH_INFO/$loginname-config1.txt");
return;
}

#**************************************
sub  output_explore
#**************************************
{
my $webPageID = shift;

my $logFileID = $hashInput{"logFileID"};
my $file_output = $hashInput{"file_output"};
my $organismID = $hashInput{"organismID"};
my $save_file = $hashInput{"save_file"};
my $file_output_title;
my $from_main_page = $hashInput{"mainPage"};
if($from_main_page){
	$file_output_title=$file_output;
	my $file_copyID = get_outputID(1);
	copy "$PATH_DATA/$loginname-$file_output", "$PATH_OUTPUT/$file_copyID.txt";
	$file_output = "$file_copyID.txt";
}elsif($webPageID){
	my $outputID = $webPageID+1;
	$file_output = "$outputID.txt";
}
if($hashInput{"file_output_title"}){
	$file_output_title = $hashInput{"file_output_title"};
}
my $file_output_full = "$PATH_OUTPUT/$file_output";
if($save_file =~ /^output$|^geneset$/){
	save_output_file($file_output_full);
}
my %hashOutput;

parse_file_headers($file_output_full,\%hashOutput,$logFileID);
my $dataset1 = $hashOutput{"output_data_set1"};
my $dataset2 = $hashOutput{"output_data_set2"};
my $nMatrix = $hashOutput{"output_n_matrixes"};
my $file_list="";
my $organismID1;
if($dataset2){
	open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Config file not found!",$logFileID);
	while(my $line = <INFO>){
		chop $line;
		if(length($line) < 3) { next; }
		my @items = split(/[=\t]/,$line);
		if($items[0] =~ /^type_/){
			if(!$file_list){ $file_list = "\"".$items[1]."\""; }
			else{ $file_list .= ",\"".$items[1]."\""; }
			if($items[1] eq $dataset2){
				my %hash=();
				read_config_line($line,\%hash);
				$organismID1 = $hash{"organismID"};
			}
		}
	}
	close INFO;
	if($dataset2=~/^public-/ && open(INFO,"<$PATH_INFO/$loginname-config.txt")){
		while(my $line = <INFO>){
			chop $line;
			my @items = split(/[=\t]/,$line);
			if($items[0] =~ /^type_/ && $dataset2 eq "public-$items[1]"){
				my %hash=();
				read_config_line($line,\%hash);
				$organismID1 = $hash{"organismID"}; last;
			}
		}
		close INFO;
	}
}
my $description = $hashOutput{"output_description"};
my $nRow = $hashOutput{"output_n_rows"};
my $nCol = $hashOutput{"output_n_columns"};
my $headers = $hashOutput{"column_titles"};
my $rows = $hashOutput{"row_titles"};
my $description_data1 = $hashOutput{"output_data_description1"};
my $description_data2 = $hashOutput{"output_data_description2"};
my $data_type1 = $hashOutput{"output_data_type1"};
my $data_type2 = $hashOutput{"output_data_type2"};
my $FDR_thresh1 = $hashOutput{"output_fdr1"};
my $FDR_thresh2 = $hashOutput{"output_fdr2"};
my $Fold_thresh1 = $hashOutput{"output_fold_change1"};
my $Fold_thresh2 = $hashOutput{"output_fold_change2"};
my $expr_thresh1 = $hashOutput{"output_expr_thresh1"};
my $expr_thresh2 = $hashOutput{"output_expr_thresh2"};
my $baseline = $hashOutput{"output_baseline"};
my $baseline1 = $hashOutput{"output_baseline1"};
my $baseline2 = $hashOutput{"output_baseline2"};
my $direction = $hashOutput{"output_direction_change"};
my $Ngenes1 = $hashOutput{"output_ngenes_file1"};
my $Ngenes2 = $hashOutput{"output_ngenes_file2"};
my $Ngenes12 = $hashOutput{"output_ngenes_common"};
my $genes_foldchange = $hashOutput{"output_method_foldchange"};
my $EPFP_angle = $hashOutput{"output_method_angle"};
my $EPFP_threshold = $hashOutput{"output_method_epfp"};
my $PAGE_cutoff = $hashOutput{"output_method_cutoff"};
my $overlap_FDR = $hashOutput{"output_fdr"};
my $overlap_enrichment = $hashOutput{"output_fold_enrichment"};
my $overlap_ngenes = $hashOutput{"output_minimum_overlap"};
my $attribute_header = $hashOutput{"output_geneset_attribute_name"};
my $par_ref = $hashOutput{"output_epfp_parameters"};
if($par_ref && ref($par_ref) eq 'ARRAY'){
	foreach my $item (@$par_ref){
		my ($name,$value) = split(/=/,$item);
		if($name eq "EPFP_threshold"){ $EPFP_threshold=$value; }
		elsif($name eq "fold_change"){ $genes_foldchange=$value; }
		elsif($name eq "angle"){ $EPFP_angle=$value; }
		elsif($name eq "cutoff"){ $PAGE_cutoff=$value; }
	}
}
my $include_graphs = 0;
if($data_type1 =~ /^(matrix|geneset)$/ && $data_type2 =~ /^(matrix|geneset)$/){
	$include_graphs = 1;
}
my $geneset_type = "overlapping";
if($data_type1 eq "matrix"){
	$data_type1 = "expression profile";
	if($data_type2 eq "matrix"){ $geneset_type = "coregulated"; }
	else{ $geneset_type = "enriched"; }
}
if($data_type2 eq "matrix"){ $data_type2 = "expression profile"; }
my $taxid = $hashOutput{"output_taxid"};
if($taxid && $taxid != $organismID){ print "In explore output organism ID don't match.<br>\n"; }
my $title_dataset1 = $dataset1;
if($dataset1 =~ /^TEMPORARY!-/i){
	$title_dataset1 = "Temporary";
}

my $method = $hashOutput{"output_method"};
my $up_down_columns = detect_up_down($headers);
my $up_down_rows = detect_up_down($rows);
my @matrix_list=();
my @matrix_genes=();
for(my $i=1; $i<=$nMatrix; ++$i){
	my $matrix_type = $hashOutput{"output_matrix$i"."_type"};
	if(ref($matrix_type) eq 'ARRAY'){
		$matrix_type = $matrix_type->[0];
	}
	my $matrix_name = $hashOutput{"output_matrix$i"."_name"};
	if($matrix_type eq "numbers" && $matrix_name !~ /^FDR|^p-value/i || $matrix_type eq "genes" || $matrix_name =~ /^gene_symbols/i){
		push(@matrix_list,[$i,$matrix_name]);
	}
	if($matrix_type eq "genes" || $matrix_name =~ /^gene_symbols/i){
		push(@matrix_genes,[$i,$matrix_name]);
	}
}
my @filter_thresholds = (0,10,20,30,40,50,60,70);
my $fileID;
if($webPageID){
	open(OUT,">$PATH_OUTPUT/$webPageID.txt") or error_message("cannot open output file ID",$logFileID);
}else{
	$fileID = get_outputID(1);
	open(OUT,">$PATH_OUTPUT/$fileID.txt") or error_message("cannot open output file ID",$logFileID);
}
print OUT "<HTML><HEAD><TITLE>ExAtlas: Explore output</TITLE>\n";
print OUT get_header();
print OUT "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print OUT "<!--\n";
print OUT "file_list = new Array($file_list);\n";
print OUT "function output_process(key){\n";
if($nCol>1 && $nRow>1){
	print OUT "	if(key == \"profile\" && document.form_output.select_column.selectedIndex==0 && document.form_output.select_row.selectedIndex==0){\n";
	print OUT "		alert(\"You need to select column or row, or both\");\n";
	print OUT "		return(false);\n";
	print OUT "	}\n";
}
print OUT "	var x = Math.round(Math.random()*10000);\n";
print OUT "	document.form_output.target = \"_BLANK\"+x;\n";
print OUT "	document.form_output.output_type.value = key;\n";
print OUT "	document.form_output.file_output.value = \"$file_output\";\n";
print OUT "	document.form_output.action.value = \"output_explore1\";\n";
print OUT "	document.form_output.save_file.value = \"\";\n";
print OUT "	document.form_output.submit();\n";
print OUT "}\n";
print OUT "function save_output_file(){\n";
print OUT "	var file = document.form_output_save.file_output_new.value;\n";
print OUT "	if(!file){\n";
print OUT "		alert(\"The output file name is missing. Please type it in\"); return false;\n";
print OUT "	}\n";
print OUT "	var file1=file;\n";
print OUT "	if(file.search(/\\.txt\$/)>=0){\n";
print OUT "		if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print OUT "		file1=file.substring(0,file.length-4);\n";
print OUT "	}\n";
print OUT "	if(file1.search(/^[-\\w]+\$/)<0){\n";
print OUT "		alert(\"File name should be one word with no special characters except underscore and dash.\\nPlease rename it!\");\n";
print OUT "		return(false);\n";
print OUT "	}\n";
print OUT "	if(file.search(/^public-/i) >= 0){\n";
print OUT "		alert(\"File name cannot start with 'public-'\"); return(false);\n";
print OUT "	}\n";
print OUT "	for(i=0; i<file_list.length; ++i){\n";
print OUT "		if(file == file_list[i]){\n";
print OUT "			alert(\"A file with this name already exists\"); return false;\n";
print OUT "		}\n";
print OUT "	}\n";
print OUT "	var descrip = document.form_output_save.description_new.value;\n";
print OUT "	if(descrip.search(/\\=|\\&/) >= 0){\n";
print OUT "		alert(\"Description should not include character \'=\' or \'&\'\");\n";
print OUT "		return false;\n";
print OUT "	}\n";
print OUT "	var x = Math.round(Math.random()*10000);\n";
print OUT "	document.form_output_save.target = \"_BLANK\"+x;\n";
print OUT "	document.form_output_save.save_file.value = \"output\";\n";
print OUT "	document.form_output_save.submit();\n";
print OUT "}\n";
print OUT "function save_geneset_file(){\n";
print OUT "	var file = document.form_geneset_save.file_geneset_new.value;\n";
print OUT "	if(!file){\n";
print OUT "		alert(\"The geneset file name is missing. Please type it in\"); return false;\n";
print OUT "	}\n";
print OUT "	var file1=file;\n";
print OUT "	if(file.search(/\\.txt\$/)>=0){\n";
print OUT "		if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print OUT "		file1=file.substring(0,file.length-4);\n";
print OUT "	}\n";
print OUT "	if(file1.search(/^[-\\w]+\$/)<0){\n";
print OUT "		alert(\"File name should be one word with no special characters except underscore and dash.\\nPlease rename it!\");\n";
print OUT "		return(false);\n";
print OUT "	}\n";
print OUT "	if(file.search(/^public-/i) >= 0){\n";
print OUT "		alert(\"Geneset file name cannot start with 'public-'\"); return(false);\n";
print OUT "	}\n";
print OUT "	for(i=0; i<file_list.length; ++i){\n";
print OUT "		if(file == file_list[i]){\n";
print OUT "			alert(\"A file with this name already exists\"); return false;\n";
print OUT "		}\n";
print OUT "	}\n";
print OUT "	var descrip = document.form_geneset_save.description_geneset_new.value;\n";
print OUT "	if(descrip.search(/\\=|\\&/) >= 0){\n";
print OUT "		alert(\"Description should not include character \'=\' or \'&\'\");\n";
print OUT "		return false;\n";
print OUT "	}\n";
print OUT "	var x = Math.round(Math.random()*10000);\n";
print OUT "	document.form_geneset_save.target = \"_BLANK\"+x;\n";
print OUT "	document.form_geneset_save.save_file.value = \"geneset\";\n";
print OUT "	document.form_geneset_save.submit();\n";
print OUT "}\n";
print OUT "// -->\n";
print OUT "</SCRIPT>\n";
if(!$file_output_title){
	my $title1 = $hashOutput{"output_name"};
	if($title1=~/TEMPORARY!-\d+.txt/){ $title1=~s/TEMPORARY!-\d+.txt/Temporary/; }
	if($title1){ $file_output_title=$title1; }
}
if($description=~/TEMPORARY!-\d+.txt/){ $description=~s/TEMPORARY!-\d+.txt/Temporary/; }
print OUT "<H3>Explore output file '$file_output_title'</H3>\n";
if($description && $description ne $file_output_title){
	my $description1 = add_hyperlinks($description); print OUT "<b>Description:</b> $description1<br>";
}
print OUT "<b>N rows</b>= $nRow &nbsp; &nbsp; &nbsp; <b>N columns</b>= $nCol<p>\n";
my @data_table;
if($dataset1 || $dataset2){
	push(@{$data_table[0]},"File name");
	push(@{$data_table[1]},$dataset1);
	push(@{$data_table[2]},$dataset2);
}
if($data_type1 || $data_type2){
	push(@{$data_table[0]},"File type");
	push(@{$data_table[1]},$data_type1);
	push(@{$data_table[2]},$data_type2);
}
if($organismID && $organismID1 && $organismID1 != $organismID){
	push(@{$data_table[0]},"Organism");
	push(@{$data_table[1]},$hashOrganism{$organismID});
	push(@{$data_table[2]},$hashOrganism{$organismID1});
}
if($description_data1 || $description_data2){
	push(@{$data_table[0]},"Description");
	push(@{$data_table[1]},$description_data1);
	push(@{$data_table[2]},$description_data2);
}
if($FDR_thresh1>0 && $FDR_thresh1<1 || $FDR_thresh2>0 && $FDR_thresh2<1){
	push(@{$data_table[0]},"FDR");
	push(@{$data_table[1]},$FDR_thresh1);
	push(@{$data_table[2]},$FDR_thresh2);
}
if(!$Fold_thresh1 && !$baseline1 && $baseline){
	$Fold_thresh1=$genes_foldchange;
	$baseline1=$baseline;
}
if($Fold_thresh1>1 || $Fold_thresh2>1){
	push(@{$data_table[0]},"Fold change");
	push(@{$data_table[1]},$Fold_thresh1);
	push(@{$data_table[2]},$Fold_thresh2);
}
if($expr_thresh1>0 || $expr_thresh2>0){
	push(@{$data_table[0]},"Expr. thresh");
	push(@{$data_table[1]},$expr_thresh1);
	push(@{$data_table[2]},$expr_thresh2);
}
if($baseline1 || $baseline2){
	push(@{$data_table[0]},"Baseline");
	push(@{$data_table[1]},$baseline1);
	push(@{$data_table[2]},$baseline2);
}
if($Ngenes1 || $Ngenes2){
	push(@{$data_table[0]},"N genes");
	push(@{$data_table[1]},$Ngenes1);
	push(@{$data_table[2]},$Ngenes2);
}
if(@data_table){
	print OUT "<TABLE BORDER=0>\n";
	print OUT "<TR><TD><b>Parental data<TD WIDTH=10><TD><b>".join("<TD WIDTH=10><TD><b>",@{$data_table[0]})."\n";
	print OUT "<TR><TD WIDTH=60><b>File 1:<TD WIDTH=10><TD><font size=-1>".join("<TD WIDTH=10><TD><font size=-1>",@{$data_table[1]})."\n";
	if($dataset2){
		print OUT "<TR><TD><b>File 2:<TD WIDTH=10><TD><font size=-1>".join("<TD WIDTH=10><TD><font size=-1>",@{$data_table[2]})."\n";
	}
	print OUT "</TABLE>\n";
}
my $spacer = " &nbsp; &nbsp; &nbsp;";
if($method){ print OUT "<b>Method:</b> $method$spacer\n"; }
if($Ngenes12){ print OUT "<b>Common genes:</b> $Ngenes12$spacer\n"; }
if($direction<0){ print OUT "<b>Expression in file 2:</b>inverted$spacer\n"; }
if($PAGE_cutoff){ print OUT "<b>PAGE cutoff prop. genes:</b> $PAGE_cutoff$spacer\n"; }
if($EPFP_angle){ print OUT "<b>Angle:</b> $EPFP_angle$spacer\n"; }
if($EPFP_threshold){ print OUT "<b>EPFP:</b> $EPFP_threshold$spacer\n"; }
if($overlap_FDR){ print OUT "<b>FDR:</b> $overlap_FDR$spacer\n"; }
if($overlap_enrichment){ print OUT "<b>Fold enrichment:</b> $overlap_enrichment$spacer\n"; }
if($overlap_ngenes){ print OUT "<b>Minimum N genes:</b> $overlap_ngenes$spacer<br>\n"; }
if($attribute_header){ print OUT "<b>Optimized by attribute:</b> $attribute_header\n"; }
my $x = int(10000*rand());
print OUT "<br><b>Output file as text:</b> <a href=$HOME_ADDRESS/output/$file_output target=_blank$x>output file</a>\n";
print OUT "<p><FORM NAME=form_output ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print OUT "<H3>1. Plot output table</H3>\n";
if($up_down_columns || $up_down_rows){
	print OUT "Data on upregulated and downregulated genes are shown separately, select the data type you need<br>\n";
}
print OUT "<TABLE BORDER=0>\n";
if($up_down_columns){
	print OUT "<TR><TD><select name=updown_columns style=width:250px;><option value=updown> Both up- and down-regulation\n";
	print OUT "<option value=up selected> Up-regulation<option value=down> Down-regulation<option value=subtract> Subtract: Up minus Down</select>\n";
	print OUT "<TD>Select \'up-\' or \'down-regulation\' for data set $title_dataset1:\n";
}
if($up_down_rows){
	print OUT "<TR><TD><select name=updown_rows style=width:250px;><option value=updown> Both up- and down-regulation\n";
	print OUT "<option value=up selected> Up-regulation<option value=down> Down-regulation</select>\n";
	print OUT "<TD>Select \'up-\' or \'down-regulation\' for data set $dataset2:\n";
}
print OUT "<TR><TD><select name=matrix_plot style=width:250px;>\n";
foreach my $ref (@matrix_list){
	my $name = $ref->[1];
	if($name =~ /^genes/){ $name = "Number of ".$name; }
	print OUT "<option value=$ref->[0]> $name\n";
}
print OUT "</select><TD COLSPAN=2>Select data type to plot\n";
print OUT "<TR><TD><select name=matrix_cluster style=width:250px;>\n";
print OUT "<option value=3> Hierarchical clustering\n";
print OUT "<option value=4> Diagonal clustering\n";
print OUT "<option value=2> Cluster columns only\n";
print OUT "<option value=1> Cluster rows only\n";
print OUT "<option value=0> No clustering\n";
print OUT "</select><TD COLSPAN=2>Clustering options\n";
print OUT "<TR><TD><select name=filter_threshold style=width:250px;>\n";
foreach my $x (@filter_thresholds){
	print OUT "<option value=$x"; if($x==20){ print OUT " selected"; } print OUT "> $x\n";
}
print OUT "</select><TD> Filtering threshold, % of maximum (0 = no filtering)\n";
print OUT "<TR><TD><select name=matrix_center style=width:250px;>\n";
print OUT "<option value=0> No centering\n";
print OUT "<option value=1> Center columns\n";
print OUT "<option value=2> Center rows\n";
print OUT "<option value=3> Center columns and rows\n";
print OUT "</select><TD> Center data (subtract median)\n";
print OUT "<TR><TD><INPUT TYPE=submit VALUE=\"Plot output table\" LANGUAGE=javascript style=width:160px; onClick=\"output_process('plot');\">\n";
print OUT "<TD><font size=-1>NOTE: Rows and columns with values &lt; threshold will be removed</font>\n";
print OUT "</TABLE><p>\n";
print OUT "<H3>2. Profiles of rows, columns, and cells</H3>\n";
if($up_down_columns || $up_down_rows){
	print OUT "Data on upregulated and downregulated genes are shown separately, select the data type you need<br>\n";
}
if(@matrix_list>1){
	print OUT "<select name=matrix_profile style=width:250px;>\n";
	foreach my $ref (@matrix_list){
		my $name = $ref->[1];
		if($name =~ /^genes/){ $name = "Number of ".$name; }
		print OUT "<option value=$ref->[0]> $name\n";
	}
	print OUT "</select> &nbsp; Select data type\n";
}else{
 	print OUT "<INPUT TYPE=hidden NAME=matrix_profile VALUE=$matrix_list[0]->[0]>\n";
}
print OUT "<TABLE>\n";
if($nCol > 1 && $nRow > 1){
	print OUT "<TR><TD WIDTH=250>File 1: <font size=-1>$title_dataset1</font><TD WIDTH=350>File 2: <font size=-1>$dataset2</file>\n";
	print OUT "<TR><TD><select name=select_column style=width:250px;>\n";
	print OUT "<option value=0> ------Select column------\n";
	my @sorted = sort {lc($headers->[$a]) cmp lc($headers->[$b])} 0..(@$headers-1);
	for(my $i=1; $i<=@$headers; ++$i){
		my $j = $sorted[$i-1]+1;
		my $title1 = $headers->[$j-1];
		if(length($title1)> 70){ $title1=substr($title1,0,70); }
		print OUT "<option value=$j> $title1\n";
	}
	print OUT "</select><TD><select name=select_row style=width:350px;>\n";
	print OUT "<option value=0> ---------Select row--------\n";
	@sorted = sort {lc($rows->[$a]) cmp lc($rows->[$b])} 0..(@$rows-1);
	for(my $i=1; $i<=@$rows; ++$i){
		my $j = $sorted[$i-1]+1;
		my $title1 = $rows->[$j-1];
		if(length($title1)> 70){ $title1=substr($title1,0,70); }
		print OUT "<option value=$j> $title1\n";
	}
	print OUT "</select>\n";
}else{
	my ($col_no,$row_no)=(0,0);
	if($nCol==1){ $col_no=1; }
	if($nRow==1){ $row_no=1; }
 	print OUT "<INPUT TYPE=hidden NAME=select_column VALUE=$col_no>\n";
 	print OUT "<INPUT TYPE=hidden NAME=select_row VALUE=$row_no>\n";
}
print OUT "<TR><TD><INPUT TYPE=button VALUE=\"Get profile\" LANGUAGE=javascript style=width:160px; onClick=\"output_process('profile');\">\n";
print OUT "<TD><TABLE BORDER=0><TR>\n";
my $notes;
if($include_graphs){
	print OUT "<TD><INPUT TYPE=checkbox NAME=graphs> Add graphics<sup>*)</sup>\n";
	$notes = "Additional graphics includes either rank-plot, or scatter-plot, or venn diagram, depending on the parental data sets.";
}
print OUT "</TABLE></TABLE>\n";
print OUT "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
print OUT "<INPUT NAME=\"action\" TYPE=hidden VALUE=\"output_explore1\">\n";
print OUT "<INPUT NAME=\"save_file\" TYPE=hidden>\n";
print OUT "<INPUT NAME=\"output_type\" TYPE=hidden>\n";
print OUT "<INPUT NAME=\"file_output\" TYPE=hidden VALUE=\"$file_output\">\n";
print OUT "<INPUT NAME=\"file_output_title\" TYPE=hidden VALUE=\"$file_output_title\">\n";
print OUT "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print OUT "</FORM><p>\n";
if(!$from_main_page || @matrix_genes){
	print OUT "<b><H3>3. Save output file</H3><TABLE>\n";
}
if(!$from_main_page){
	my $description1;
	if($description){ $description1 = $description; }
	my $title1;
	print OUT "<p><FORM NAME=form_output_save ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
	print OUT "<TABLE>\n";
	print OUT "<TR><TD WIDTH=150>Save output file as:<TD><INPUT NAME=file_output_new SIZE=35>\n";
	print OUT "<TD><INPUT TYPE=button VALUE=\"Save output\" onClick=save_output_file(); style=width:160px;><p>\n";
	print OUT "<TR><TD>Description: <TD><INPUT NAME=description_new SIZE=35 VALUE=\"$description1\">\n";
	print OUT "</TABLE>\n";
	print OUT "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
	print OUT "<INPUT NAME=\"action\" TYPE=hidden VALUE=\"output_explore\">\n";
	print OUT "<INPUT NAME=\"save_file\" TYPE=hidden>\n";
	print OUT "<INPUT NAME=\"file_output\" TYPE=hidden VALUE=\"$file_output\">\n";
	print OUT "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID></FORM><p>\n";
}
if(@matrix_genes){
	print OUT "<p><FORM NAME=form_geneset_save ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
	print OUT "<TABLE>\n";
	print OUT "<b>Save lists of genes as a new geneset file</b><br><TABLE>\n";
	print OUT "<TR><TD WIDTH=150>Geneset file name:<TD><INPUT NAME=file_geneset_new SIZE=35>\n";
	print OUT "<TD><INPUT TYPE=button VALUE=\"Save lists of genes\" onClick=save_geneset_file(); style=width:160px;><p>\n";
	print OUT "<TR><TD>Description: <TD><INPUT NAME=description_geneset_new SIZE=35>\n";
	if(@matrix_genes>1){
		print OUT "<TR><TD>Select geneset data<TD><select name=matrix_genes style=width:250px;>\n";
		foreach my $ref (@matrix_genes){
			my $name = $ref->[1];
			print OUT "<option value=$ref->[0]> $name\n";
		}
		print OUT "</select>\n";
	}else{
 		print OUT "<INPUT TYPE=hidden NAME=matrix_genes VALUE=$matrix_genes[0]->[0]>\n";
	}
	print OUT "</TABLE>\n";
	print OUT "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
	print OUT "<INPUT NAME=\"action\" TYPE=hidden VALUE=\"output_explore\">\n";
	print OUT "<INPUT NAME=\"save_file\" TYPE=hidden>\n";
	print OUT "<INPUT NAME=\"file_output\" TYPE=hidden VALUE=\"$file_output\">\n";
	print OUT "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID></FORM><p>\n";
}
if($include_graphs){
	print OUT "<sup>*)</sup><b>Notes:</b> $notes<p>\n";
}
print OUT "<p><HR NOSHADE><p>\n";
print OUT "<p><i>Please report any problems to <a href=mailto:sharoval\@mail.nih.gov>webmaster</a><p>\n";
print OUT "<INPUT TYPE=button VALUE=\"Close window\" style=width:160px; LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print OUT "</BODY></HTML>\n";
close OUT;
if(!$webPageID){
	print_web_page($fileID);
	exit(0);
}
return;
}

#**************************************
sub  output_explore1
#**************************************
{
my $file_output = $hashInput{"file_output"};
my $matrix_plot = $hashInput{"matrix_plot"}; # Matrix number
my $filter_threshold = $hashInput{"filter_threshold"};
my $matrix_cluster = $hashInput{"matrix_cluster"};
my $matrix_center = $hashInput{"matrix_center"};
my $updown_columns = $hashInput{"updown_columns"};
my $updown_rows = $hashInput{"updown_rows"};
my $output_type = $hashInput{"output_type"};
my $organismID = $hashInput{"organismID"};
my $file_output_title = $hashInput{"file_output_title"};

my $file_output_full = "$PATH_OUTPUT/$file_output";
if(!$file_output){ error_message("Missing output file name"); }
my $diagonal = 0;
if($matrix_cluster==4){
	$matrix_cluster=3;
	$diagonal=1;
}
my %hashOutput;
parse_file_headers($file_output_full,\%hashOutput);
my $nMatrix = $hashOutput{"output_n_matrixes"};
my $nRow = $hashOutput{"output_n_rows"};
my $nCol = $hashOutput{"output_n_columns"};
my $headers = $hashOutput{"column_titles"};
my $rows = $hashOutput{"row_titles"};
my $taxid = $hashOutput{"output_taxid"};
if($taxid && !$organismID){
	$organismID = $taxid;
}
my $organismID1 = $hashOutput{"output_taxid2"};
if(!$organismID1){ $organismID1=$organismID; }
my $database_headers = $hashOutput{"output_database_headers"};
my $autocorrelation = $hashOutput{"output_autocorrelation"};
my @matrix_names;
my @matrix_types;
my $use_database=0;
for(my $im=1; $im<=$nMatrix; ++$im){
	my $matrix_name = $hashOutput{"output_matrix$im"."_name"};
	push(@matrix_names,$matrix_name);
	my $matrix_type = $hashOutput{"output_matrix$im"."_type"};
	if(ref($matrix_type) eq 'ARRAY'){
		if($matrix_type->[1] eq "database"){
			$use_database=1;
		}
		$matrix_type = $matrix_type->[0];
	}
	push(@matrix_types,$matrix_type);
}
my @database;
my @database_headers;
if($use_database){
	if(!$database_headers){ error_message("In output_explore1 - no database_headers"); }
	if(ref($database_headers) eq 'ARRAY'){
		@database_headers = @$database_headers;
	}else{
		@database_headers = ($database_headers);
	}
	open(INFO,'<',$file_output_full) or error_message("In output_explore1 - cannot open file");
	while(my $line = <INFO>){
		if($line =~ /^!database_start/i){ last; }
	}
	while(my $line = <INFO>){
		chop $line;
		if($line =~ /^!/){ last; }
		if($line !~ /^>/ || !$line){ next; }
		my $data_id = $line;
		$data_id =~ s/^>//;
		foreach my $header (@database_headers){
			my $line1 = <INFO>;
			chop $line1;
			my($key,$value)=split(/\t/,$line1);
			if($key ne $header){ error_message("Database header does not match"); }
			if($key eq $database_headers[0]){
				$database[$data_id] = $value;
			}
		}
		if(!$database[$data_id]){ error_message("\'genes\' not found in database ID=$data_id"); }
	}
	close INFO;
}
my @remove_columns;
my @remove_rows;
my ($nRow1,$nCol1)=($nRow,$nCol);
my @subtractCol;
my @subtractSign;
my %hashCol;
my @col_new;
for(my $i=0; $i<$nCol && $updown_columns && $updown_columns ne "updown"; $i++){
	if($updown_columns eq "subtract"){
		my $header1 = $headers->[$i];
		$header1 =~ s/_up$|_down$//i;
		my $ii = $hashCol{$header1};
		if(!$ii){
			push(@col_new,$header1);
			$hashCol{$header1} = @col_new;
			$ii = @col_new;
		}
		$subtractCol[$i] = $ii;
		$subtractSign[$i] = 1;
		if($headers->[$i] =~ /_down$/){ $subtractSign[$i] = -1; }
	}elsif($headers->[$i] !~ /_$updown_columns$/i){
		push(@remove_columns,1);
		$nCol1--;
	}else{ push(@remove_columns,0); }
}
for(my $i=0; $i<$nRow && $updown_rows && $updown_rows ne "updown"; $i++){
	if($rows->[$i] !~ /_$updown_rows$/i){
		push(@remove_rows,1);
		$nRow1--;
	}
	else{ push(@remove_rows,0); }
}
open(INFO,'<',$file_output_full) or error_message("In output_explore1 - cannot open file");
if($output_type eq "plot"){
	if($updown_columns eq "subtract"){
		$nCol1 = @col_new;
		@$headers = @col_new;
	}
	my $matrix_name = $matrix_names[$matrix_plot-1];
	my $matrix_type = $matrix_types[$matrix_plot-1];
	for(my $i=$nCol-1; $i>=0; $i--){
		if($remove_columns[$i]){
			splice(@$headers,$i,1);
		}
	}
	my $searchline = "!Matrix$matrix_plot"."_start";
	while(my $line = <INFO>){
		if($line =~ /^$searchline/){ last; }
	}
	my $line = <INFO>;
	chop $line;
	my @data;
	my @data_columns;
	@$rows=();
	my $irow=0;
	while(my $line = <INFO>){
		if($line =~ /^!/){ last; }
		chop $line;
		my ($row_header,@data1)=split(/\t/, $line);
		if($autocorrelation){
			if($matrix_name =~ /^z-value|^fold/){ $data1[$irow]=1000; }
			if($matrix_name =~ /^correlation/){ $data1[$irow]=1; }
		}
		if($remove_rows[$irow]){
			$irow++;
			next;
		}
		my @data2;
		for(my $i=$nCol-1; $i>=0; $i--){
			if($remove_columns[$i] && $i<@data1){
				splice(@data1,$i,1);
				next;
			}
			if($matrix_type eq "genes" || $matrix_name =~ /^gene_symbols/i){
				if(!$data1[$i]){
					$data1[$i] = 0;
					next;
				}
				if($use_database){
					$data1[$i] = $database[$data1[$i]];
				}
				$data1[$i] = split(/,/,$data1[$i]);
			}
			my $ii = $subtractCol[$i];
			if($ii){
				if($data1[$i]>=$MISSING){
					$data2[$ii-1] += $subtractSign[$i]*$data1[$i];
				}
			}
		}
		if(@subtractCol){ @data1 = @data2; }
		if(!@data1){ @data1=(); }
		push(@data,\@data1);
		push(@$rows,$row_header);
		$irow++;
	}
	close INFO;
	if($matrix_center==1 || $matrix_center==3){
		for(my $i=0; $i<$nRow; ++$i){
			my $ref = $data[$i];
			for(my $j=0; $j<$nCol; ++$j){
				push(@{$data_columns[$j]},$ref->[$j]);
			}
		}
		for(my $j=0; $j<$nCol; ++$j){
			my $median = median($data_columns[$j]);
			if($median>$MISSING){
				for(my $i=0; $i<$nRow; ++$i){
					my $x = $data[$i]->[$j];
					if($x > $MISSING){
						$x = floor(10000*($x-$median)+0.5)/10000;
						$data[$i]->[$j] = $x;
					}
				}
			}
		}
	}
	if($matrix_center==2 || $matrix_center==3){
		for(my $i=0; $i<$nRow; ++$i){
			my $median = median($data[$i]);
			if($median>$MISSING){
				for(my $j=0; $j<$nCol; ++$j){
					my $x = $data[$i]->[$j];
					if($x > $MISSING){
						$x = floor(10000*($x-$median)+0.5)/10000;
						$data[$i]->[$j] = $x;
					}
				}
			}
		}
	}
	if($filter_threshold){
		my $xmax = -100000;
		my @maxRow;
		my @maxCol;
		for(my $icol=0; $icol<$nCol1; ++$icol){
			$maxCol[$icol] = -100000;
		}
		for(my $irow=0; $irow<$nRow1; ++$irow){
			$maxRow[$irow] = -100000;
			my $ref = $data[$irow];
			for(my $icol=0; $icol<$nCol1; ++$icol){
				my $x = $ref->[$icol];
				if($x==$MISSING || $autocorrelation && $icol==$irow){ next; }
				if($xmax<$x){ $xmax=$x; }
				if($maxCol[$icol]<$x){ $maxCol[$icol]=$x; }
				if($maxRow[$irow]<$x){ $maxRow[$irow]=$x; }
			}
		}
		$filter_threshold = $xmax*$filter_threshold/100;
		for(my $irow=$nRow1-1; $irow>=0; $irow--){
			if($maxRow[$irow] < $filter_threshold){
				splice(@data,$irow,1);
				splice(@$rows,$irow,1);
				$nRow1--;
			}
		}
		for(my $icol=$nCol1-1; $icol>=0; $icol--){
			if($maxCol[$icol] < $filter_threshold){
				splice(@$headers,$icol,1);
				$nCol1--;
				foreach my $ref1 (@data){
					splice(@$ref1,$icol,1);
				}
			}
		}
	}
	if(!@data){ error_message("Empty matrix after filtering!"); }
	for(my $i=0; $i<$nCol; ++$i){
		my $len = length($headers->[$i]);
		if($len>$maxHeaderLength){ $headers->[$i]=substr($headers->[$i],0,$maxHeaderLength); }
	}
	for(my $i=0; $i<$nRow; ++$i){
		my $len = length($rows->[$i]);
		if($len>$maxHeaderLength){ $rows->[$i]=substr($rows->[$i],0,$maxHeaderLength); }
	}
	my $fileID = get_outputID(3);
	open(OUT, ">$PATH_OUTPUT/$fileID.txt") or error_message("In output_explore1 cannot open file $fileID.txt"); 
	print OUT "Headers";
	for(my $i=0; $i<$nCol1; ++$i){
		print OUT "\t$headers->[$i]";
	}
	print OUT "\n";
	my $logtransform = 0;
	if($matrix_type =~ /genes|gene_symbols/i || $matrix_name =~ /^gene_symbols/i){
		$logtransform = 1;
		$matrix_name =~ s/^genes_/genes /;
		$matrix_name = "Number of $matrix_name";
	}
	for(my $irow=0; $irow<$nRow1; ++$irow){
		print OUT $rows->[$irow];
		my $ref = $data[$irow];
		for(my $i=0; $i<$nCol1; ++$i){
			my $x=$ref->[$i];
			if($logtransform){
				if($x<0){ $x=0; }
				else{ $x = int(100000*log($x+1))/100000; }
			}
			print OUT "\t$x";
		}
		print OUT "\n";
	}
	close OUT;
	$hashInput{"runID"}=$RUN_PCA_CLUSTER;
	my $logFileID = get_outputID(3);
	$hashInput{"logFileID"}=$logFileID;
	my $data = [$fileID,$nCol1,$nRow1,$matrix_cluster,$diagonal,$logtransform,$matrix_name];
	my $evaluate = $nCol1*$nRow1/500000;
	if($evaluate < $nCol1/400){ $evaluate < $nCol1/400; }
	interrupt_program($evaluate,$data);
	exit(0);
}
if($output_type eq "profile"){
	my $select_row = $hashInput{"select_row"};
	my $select_col = $hashInput{"select_column"};
	my $FDR_thresh = $hashOutput{"output_fdr"};
	my $fold_thresh = $hashOutput{"output_fold_enrichment"};
	my $output_type_name = "correlation";
	if($hashOutput{"output_data_type2"} eq "geneset"){
		$output_type_name = "geneset enrichment";
		if($hashOutput{"output_data_type1"} eq "geneset"){
			$output_type_name = "geneset overlap";
		}
	}
	my $iMatrix = $hashInput{"matrix_profile"};
	my $matrix_profile_name = $matrix_names[$iMatrix-1];
	my $ending="";
	if($matrix_profile_name =~ /_upregulated$/i){ $ending="_upregulated"; }
	elsif($matrix_profile_name =~ /_downregulated$/i){ $ending="_downregulated"; }
	elsif($matrix_profile_name =~ /_pos$/i){ $ending="_pos"; }
	elsif($matrix_profile_name =~ /_neg$/i){ $ending="_neg"; }
	my @data;
	my ($iGenes,$iFDR,$iFold,$iZvalue,$iCorrel,$iAttrib,$iGenesNeg) = (-1,-1,-1,-1,-1,-1,-1);
	my $iMatrix1 = $iMatrix - 1;
	for(my $im=0; $im<$nMatrix; ++$im){
		my $matrix_name = $matrix_names[$im];
		my $matrix_type = $matrix_types[$im];
		if($iGenes<0 && ($matrix_type =~ /^gene/i || $matrix_name =~ /^gene/) || $ending && ($matrix_name eq "genes$ending" || $matrix_name eq "gene_symbols$ending")){
			$iGenes=$im;
		}elsif($iFDR<0 && $matrix_name =~ /^FDR/i || $ending && $matrix_name eq "FDR$ending"){
			$iFDR=$im;
		}elsif($iFold<0 && $matrix_name =~ /^fold_enrichment/i || $ending && $matrix_name eq "fold_enrichment$ending"){
			$iFold=$im;
		}elsif($iZvalue<0 && $matrix_name =~ /^z-value/i || $ending && $matrix_name eq "z-value$ending"){
			$iZvalue=$im;
		}elsif($iCorrel<0 && $matrix_name =~ /^correlation/i || $ending && $matrix_name eq "correlation$ending"){
			$iCorrel=$im;
		}elsif($iAttrib<0 && $matrix_name =~ /^attribute_threshold/i || $ending && $matrix_name eq "attribute_threshold$ending"){
			$iAttrib=$im;
		}elsif($iGenes>=0 && $iGenesNeg<0 && $ending !~ /_upregulated|_pos/ && $matrix_name =~ /genes_(downregulated|neg)/i){
			$iGenesNeg=$im;
		}
	}
	my $database_id;
	my $database_id_neg;
	for(my $im=1; $im<=$nMatrix; ++$im){
		my $matrix_name = $matrix_names[$im-1];
		my $matrix_type = $matrix_types[$im-1];
		my $searchline = "!Matrix$im"."_start";
		my $line;
		while($line = <INFO>){
			if($line =~ /^$searchline/){ last; }
		}
		if($line !~ /^$searchline/){ error_message("Matrix$im not found."); }
		$line = <INFO>;
		if(!$select_col){
			for(my $i=0; $i<$select_row; ++$i){ $line = <INFO>; }
			chop $line;
			my ($row_header,@data1)=split(/\t/, $line);
			for(my $j=0; $j<@data1; ++$j){
				if($iGenes==$im-1 && $use_database && $data1[$j]){
					$data1[$j] = $database[$data1[$j]];
				}elsif($iGenesNeg==$im-1 && $use_database && $data1[$j]){
					$data1[$j] = $database[$data1[$j]];
				}
			}
			$data[$im-1] = \@data1;
		}elsif(!$select_row){
			for(my $i=0; $i<$nRow; ++$i){
				$line = <INFO>;
				chop $line;
				my ($row_header,@data1)=split(/\t/, $line);
				my $x = $data1[$select_col-1];
				if($iGenes==$im-1 && $use_database && $x){
					$x = $database[$x];
				}elsif($iGenesNeg==$im-1 && $use_database && $x){
					$x = $database[$x];
				}
				push(@{$data[$im-1]},$x);
			}
		}else{
			for(my $i=0; $i<$select_row; ++$i){ $line = <INFO>; }
			chop $line;
			my ($row_header,@data1)=split(/\t/, $line);
			my $x = $data1[$select_col-1];
			if($iGenes==$im-1 && $use_database && $x){
				$database_id = $x;
				$x = $database[$x];
			}
			if($iGenesNeg==$im-1 && $use_database && $x){
				$database_id_neg = $x;
				$x = $database[$x];
			}
			$data[$im-1] = $x;
		}
	}
	my @attributes;	#Additional gene attributes
	my @attributes_neg;	#Additional gene attributes
	my $n_database = 0;
	if($database_id){ $n_database++; }
	if($database_id_neg){ $n_database++; }
	if(@database_headers>1 && $select_col && $select_row && $n_database){
		while(my $line = <INFO>){
			if(!$n_database){ last; }
			if($database_id && $line =~ /^>$database_id\n/ || $database_id_neg && $line =~ /^>$database_id_neg\n/){
				$n_database--;
				foreach my $header (@database_headers){
					my $line1 = <INFO>;
					chop $line1;
					my($key,$value)=split(/\t/,$line1);
					if($key =~ /^genes/){ next; }
					if($line =~ /^>$database_id\n/){ push(@attributes,$value); }
					elsif($line =~ /^>$database_id_neg\n/){ push(@attributes_neg,$value); }
				}
			}
		}
	}
	close INFO;
	
	my @output_order;
	if($iZvalue>=0){ push(@output_order,$iZvalue); }
	if($iCorrel>=0){ push(@output_order,$iCorrel); }
	if($iFDR>=0){ push(@output_order,$iFDR); }
	elsif($iZvalue>=0){ push(@output_order,$nMatrix); }
	if($iFold>=0){ push(@output_order,$iFold); }
	if($iAttrib>=0){ push(@output_order,$iAttrib); }
	if($iGenes>=0){ push(@output_order,$iGenes); }
	if($iGenesNeg>=0){ push(@output_order,$iGenesNeg); }

	print "<HTML><HEAD><TITLE>ExAtlas: Output profile</TITLE>\n";
	if($select_col && $select_row){
		print_header("update_description();");
	}else{
		print_header();
	}
	my $plot_type;
	if($select_col){ $plot_type = "Data#1 = $headers->[$select_col-1]"; }
	if($select_row){ $plot_type .= " Data#2 = $rows->[$select_row-1]"; }

	my $name = $matrix_profile_name;
	if($name =~ /^genes/){ $name = "Number of ".$name; }
	print "<SCRIPT language=JavaScript>\n";
	print "<!--\n";
	print "function plot(num) {\n";
	if($select_col){ print "	document.form_profile.select_row.value=num;\n"; }
	elsif($select_row){ print "	document.form_profile.select_column.value=num;\n"; }
	print "	var x = Math.round(Math.random()*10000);\n";
	print "	document.form_profile.target = \"_BLANK\"+x;\n";
	print "	document.form_profile.submit();\n";
	print "}\n";
	print "<!-- end script --></SCRIPT>\n";
	print "<H2>Output profile: $name</H2><p>\n";
	print "<b><u>Description:</u></b> $output_type_name: $plot_type<br>\n";
	print "<b><u>Organism:</u></b> $hashOrganism{$organismID}<br>\n";
	my $data_filtered=0;
	if(!$select_col || !$select_row){
		#O U T P U T    O F   M U L T I P L E    M A T R I X    E L E M E N T S

		my $plotID = get_outputID(2);
		my $fileID = get_outputID(1);
		my $x = int(10000*rand());
		print "<b><u>Table as text:</u></b> <a href=\"$HOME_ADDRESS/output/$fileID.txt\" target=_BLANK$x>table</a><p>\n";
		my @Ngenes;
		my @NgenesNeg;
		my $N = @{$data[0]};
		for(my $i=0; $i<$N && $iGenes>=0; ++$i){
			if($data[$iGenes]->[$i]){
				$Ngenes[$i] = split(/,/,$data[$iGenes]->[$i]);
			}else{
				$Ngenes[$i]=0;
			}
			if($iGenesNeg>=0){
				if($data[$iGenesNeg]->[$i]){
					$NgenesNeg[$i] = split(/,/,$data[$iGenesNeg]->[$i]);
				}else{
					$NgenesNeg[$i]=0;
				}
			}
		}
		my $ptr = $rows;
		if(!$select_col){
			$ptr = $headers;
		}
		my @sorted;
		my @sortedZ;
		my @pvalue;
		my @FDR;
		if($iZvalue>=0 && $iFDR<0){
			@sortedZ = sort {abs($data[$iZvalue]->[$b])<=>abs($data[$iZvalue]->[$a])} (0..($N-1));
			my $FDR1=1;
			for(my $i=@sortedZ-1; $i>=0; $i--){
				my $j = $sortedZ[$i];
				$pvalue[$j] = 2*(1 - normal_distribution(abs($data[$iZvalue]->[$j])));
				my $FDR = int(10000*$pvalue[$j]*$N/($i+1)+0.5)/10000;
				if($FDR > $FDR1){ $FDR = $FDR1; }
				else{ $FDR1 = $FDR; }
				$FDR[$j] = $FDR;
			}
		}
		if($iMatrix1 == $iGenes){
			@sorted = sort {$Ngenes[$b]<=>$Ngenes[$a]} (0..($N-1));
		}else{
			@sorted = sort {$data[$iMatrix1]->[$b]<=>$data[$iMatrix1]->[$a]} (0..($N-1));
		}
		for(my $i=@sorted-1; $i>=0 && @sorted>30; $i--){
			my $j = $sorted[$i];
			if($iZvalue>=0 && ($output_type_name=~/^geneset/ && $data[$iZvalue]->[$j] < 2 || $output_type_name !~/^geneset/ && abs($data[$iZvalue]->[$j]) < 2) ||
			   $iFDR>=0 && $FDR_thresh>0 && $data[$iFDR]->[$j] > $FDR_thresh ||
			   $iFold>=0 && $fold_thresh>0 && $data[$iFold]->[$j] < $fold_thresh ||
			   $iMatrix1==$iGenes && $Ngenes[$j] < 3){
				splice(@sorted,$i,1);
			}
		}
		for(my $im=0; $im<$nMatrix; ++$im){
			my @new;
			for(my $i=0; $i<@sorted; ++$i){
				$new[$i] = $data[$im]->[$sorted[$i]];
			}
			$data[$im] = \@new;
		}
		if($iZvalue>=0 && $iFDR<0){
			my @new;
			for(my $i=0; $i<@sorted; ++$i){
				push(@new,$FDR[$sorted[$i]]);
			}
			$data[$nMatrix] = \@new;
		}
		my @new1;
		for(my $i=0; $i<@sorted; ++$i){ $new1[$i] = $Ngenes[$sorted[$i]]; }
		@Ngenes = @new1;
		if($iGenesNeg){
			for(my $i=0; $i<@sorted; ++$i){ $new1[$i] = $NgenesNeg[$sorted[$i]]; }
			@NgenesNeg = @new1;
		}
		my @new;
		for(my $i=0; $i<@sorted; ++$i){
			$new[$i] = $ptr->[$sorted[$i]];
		}
		$ptr = \@new;
		$N = @sorted;
		my $thresh=0;
		if($iMatrix1 == $iGenes){
			plot_histogram_horiz("bars",$N,1,\@Ngenes,$ptr,$plotID,$name,"zero");
		}else{
			my $option_zero="";
			if($matrix_profile_name =~ /z-value/){ $option_zero="zero"; }
			plot_histogram_horiz("bars",$N,1,$data[$iMatrix1],$ptr,$plotID,$name,$option_zero);
		}
		print "<IMG SRC=../output/$plotID.gif BORDER=0><p>\n";
		if($data_filtered){
			print "<b>Note:</b> data filtered out if z < $thresh<p>\n";
		}
		print "<H3>Output data table</H3>\n";
		open(OUT,">$PATH_OUTPUT/$fileID.txt") or error_message("Cannot open fileID");
		print "<TABLE border=1><TR><TD><b>Title\n";
		print OUT "Title";
		my ($colWidth, $geneCount) = (500,50);
		my $geneDir = "";
		if($iGenesNeg>=0){ $geneDir = "<br>upreg."; ($colWidth, $geneCount) = (200,20); }
		foreach my $i (@output_order){
			if($i==$iGenes){
				print "<TD><b>N genes$geneDir<TD WIDTH=500><font size=+1><b>$matrix_names[$i]";
				$geneDir=~ s/<br>/ /;
				print OUT "\tN genes$geneDir\t$matrix_names[$i]";
			}elsif($i==$iGenesNeg){
				print "<TD><b>N genes<br>downreg.<TD WIDTH=500><font size=+1><b>$matrix_names[$i]";
				print OUT "\tN genes downreg.\t$matrix_names[$i]";
			}elsif($i==$nMatrix){
				print "<TD><b>FDR";
				print OUT "\tFDR";
			}else{
				my $title = $matrix_names[$i];
				$title =~ s/_/ /g;
				print "<TD><font size=+1><b>$title";
				print OUT "\t$title";
			}
		}
		for(my $i=0; $i<$N; ++$i){
			my $item = $i+1;
			if(@sorted){ $item=$sorted[$i]+1; }
			if(!$select_col){
				$item = $i+1;
				if(@sorted){ $item = $sorted[$i]+1; }
			}
			print "\n<TR><TD><a href=\"JavaScript:plot('$item');\">$ptr->[$i]<a>";
			print OUT "\n$ptr->[$i]";
			foreach my $j (@output_order){
				if($j==$iGenes){
					my $genes = $data[$j]->[$i];
					if($genes eq "0"){ $genes=""; }
					print OUT "\t$Ngenes[$i]\t$genes";
					if($Ngenes[$i] > $geneCount){
						my @symbols = split(/,/,$genes);
						splice(@symbols,$geneCount);
						$genes = join(", ",@symbols).",... (click on the name to see the full list)"; 
					}else{
						$genes =~ s/,/, /g;
					}
					print "<TD><center>$Ngenes[$i]<TD>$genes";
				}elsif($j==$iGenesNeg){
					my $genes = $data[$j]->[$i];
					if($genes eq "0"){ $genes=""; }
					print OUT "\t$NgenesNeg[$i]\t$genes";
					if($NgenesNeg[$i] > $geneCount){
						my @symbols = split(/,/,$genes);
						splice(@symbols,$geneCount);
						$genes = join(", ",@symbols).",... (click on the name to see the full list)"; 
					}else{
						$genes =~ s/,/, /g;
					}
					print "<TD><center>$NgenesNeg[$i]<TD>$genes";
				}else{
					print "<TD><center>$data[$j]->[$i]";
					print OUT "\t$data[$j]->[$i]";
				}
			}
		}
		print "\n</TABLE>\n";
		print OUT "\n";
		close OUT;
		print "<FORM NAME=form_profile ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
		print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
		print "<INPUT NAME=action TYPE=hidden VALUE=\"output_explore1\">\n";
		print "<INPUT NAME=file_output TYPE=hidden VALUE=\"$file_output\">\n";
		print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID><p>\n";
		print "<INPUT NAME=matrix_profile TYPE=hidden VALUE=$iMatrix><p>\n";
		print "<INPUT NAME=select_column TYPE=hidden VALUE=$select_col><p>\n";
		print "<INPUT NAME=select_row TYPE=hidden VALUE=$select_row><p>\n";
		print "<INPUT NAME=output_type TYPE=hidden VALUE=profile><p>\n";
		if($hashInput{"graphs"} eq "on"){
			print "<INPUT NAME=graphs TYPE=hidden VALUE=on><p>\n";
		}
		print "</FORM></HTML>\n";
	}else{
		#O U T P U T    O F   A   S I N G L E    M A T R I X    E L E M E N T
		my $data_type1 = $hashOutput{"output_data_type1"};
		my $data_type2 = $hashOutput{"output_data_type2"};
		my $dataset1 = $hashOutput{"output_data_set1"};
		my $dataset2 = $hashOutput{"output_data_set2"};
		my $baseline_name = $hashOutput{"output_baseline"};
		my ($col1,$row1) = ($select_col,$select_row);
		my $ref1 = $hashOutput{"output_original_order1"};
		if($ref1){ if(ref($ref1) eq 'ARRAY'){ $col1 = $ref1->[$col1-1]; } elsif($col1==1){ $col1=$ref1; }}
		my $ref2 = $hashOutput{"output_original_order2"};
		if($ref2){ if(ref($ref2) eq 'ARRAY'){ $row1 = $ref2->[$row1-1]; } elsif($row1==1){ $row1=$ref2; }}
		my $hyperlink=0;
		if($data_type1 eq "matrix"){ $hyperlink=1; }
		print "<H3>Output data table</H3>\n";
		print "<TABLE border=0><TR><TD><b>Title<TD WIDTH=20><TD WIDTH=800><b>Value\n";
		my @gene_list=();
		my @gene_list_neg=();
		my ($nGenes,$nGenesNeg)=(0,0);
		my %hashGenes=();
		my ($geneDir) = ("");
		if($iGenesNeg>=0){ $geneDir = " upreg."; }
		foreach my $j (@output_order){
			my $x = $data[$j];
			if($j==$iGenes){
				@gene_list = split(/,/,$x);
				$nGenes = @gene_list;
				print "<TR><TD><font size=+1>N genes$geneDir<TD><TD>$nGenes\n";
			}elsif($j==$iGenesNeg){
				@gene_list_neg = split(/,/,$x);
				$nGenesNeg = @gene_list_neg;
				if($nGenesNeg){
					print "<TR><TD><font size=+1>N genes downreg.<TD><TD>$nGenesNeg\n";
				}
			}else{
				my $title = $matrix_names[$j];
				$title =~ s/_/ /g;
				print "<TR><TD><font size=+1>$title<TD><TD>$x\n";
			}
		}
		print "\n</TABLE>\n";
		print "<SCRIPT language=JavaScript>\n";
		print "<!--\n";
		print "function plot(symbol) {\n";
		print "	document.form_symbols.search_term.value=symbol;\n";
		print "	var x = Math.round(Math.random()*10000);\n";
		print "	document.form_symbols.target = \"_BLANK\"+x;\n";
		print "	document.form_symbols.submit();\n";
		print "}\n";
		my $gene_attribute_min=0;
		if($iGenes>=0 && $nGenes>0){
			my $fileID = get_outputID(1);
			if($nGenes>=5 || $nGenesNeg>=5){
				my @geneset_list = get_geneset_list();
				filter_list_by_organism(\@geneset_list, $organismID);
				my ($items,$descriptions) = get_array_lists(\@geneset_list);
				print "geneset_list = new Array($items);\n";
				print "geneset_description = new Array($descriptions);\n";
				print "function geneset_overlap(file) {\n";
				print "	document.form_genelist.upload_geneset.value=file;\n";
				print "	document.form_genelist.action.value=\"geneset_overlap\";\n";
				print "	var x = Math.round(Math.random()*10000);\n";
				print "	document.form_genelist.target = \"_BLANK\"+x;\n";
				print "	document.form_genelist.submit();\n";
				print "}\n";
				print "function update_description() {\n";
				print "	var index;\n";
				print "	index = document.form_genelist.file_geneset1.selectedIndex;\n";
				print "	document.form_genelist.description_geneset1.value = geneset_description[index];\n";
				print "}\n";
				print "<!-- end script --></SCRIPT>\n";
				my $x = int(10000*rand());
				print "<b>Gene list as tab-delimited text:</b> <a href=\"$HOME_ADDRESS/output/$fileID.txt\" target=_BLANK$x>gene list</a><p>\n";
				print "<FORM NAME=form_genelist ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
				print "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
				print "<INPUT NAME=action TYPE=hidden>\n";
				print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
				print "<INPUT NAME=upload_geneset TYPE=hidden>\n";
				print "<TABLE BORDER=0>\n";
				my $menu_text = menu_geneset_overlap(\@geneset_list,$fileID);
				print $menu_text;
				print "</TABLE></FORM><p>\n";
			}else{
				print "<!-- end script --></SCRIPT>\n";
			}
			get_gene_annotation($organismID,\@gene_list,\%hashGenes,"symbol");
			my @attributes1=();
			my @attrib_names1=();
			if($data_type1 eq "geneset"){
				@attrib_names1 = get_gene_attributes($dataset1,$col1,\@gene_list,\@attributes1,$organismID,$organismID);
				if($attrib_names1[0]==NULL){ @attrib_names1=(); }
			}
			my @attributes2=();
			my @attrib_names2=();
			if($data_type2 eq "geneset"){
				@attrib_names2 = get_gene_attributes($dataset2,$row1,\@gene_list,\@attributes2,$organismID,$organismID1);
				if($attrib_names2[0]==NULL){ @attrib_names2=(); }
			}
			my $headers = $hashGenes{"HEADERS"};
			open(OUT,">$PATH_OUTPUT/$fileID.txt") or error_message("Cannot open fileID");
			if($matrix_profile_name=~/upregulated/ || $geneDir =~/upreg/){ $geneDir = "upregulated"; }
			if($matrix_profile_name=~/downregulated/){ $geneDir = "downregulated"; }
			print OUT "!Genelist_description\t$output_type_name $geneDir: $plot_type\n";
			print "<FORM NAME=form_symbols ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
			print "<FONT SIZE=+2>List of genes $geneDir</FONT><p>\n";
			my $header_line_html = "<TABLE><TR><TD><b>Symbol";
			my $header_line_text = "Symbol";
			for(my $i1=0; $i1<@attrib_names1; $i1++){
				my $header = $attrib_names1[$i1];
				$header =~ s/_/ /g;
				$header_line_text .= "\t$header";
				if($header=~/^gene[ _]*(title|name)$/){
					$header_line_html .= "<TD WIDTH=450><b>$header";
				}else{
					$header_line_html .= "<TD><b>$header";
				}
			}
			my $sort_by = -1;
			my $reverse = 0;
			for(my $i1=0; $i1<@attrib_names2; $i1++){
				my $header = $attrib_names2[$i1];
				$header =~ s/_/ /g;
				if($i1==0){ 
					$sort_by=0;
					if($header =~ /^(FDR|EPFP|p|p-value)$/i){ $reverse = 1; }
				}
				if($sort_by==0 && $header =~ /^logratio|^log-ratio|^ClipRPM/i){
					$sort_by = $i1;
					$reverse = 0;
				}
				$header_line_text .= "\t$header";
				$header_line_html .= "<TD><b>$header";
			}
			my @gene_attributes;
			for(my $ii=1; $ii<@database_headers; $ii++){
				my $x = $database_headers[$ii];
				$header_line_text .= "\t$x";
				$header_line_html .= "<TD><b>$x";
				my @values = split(/,/,$attributes[$ii-1]);
				push(@gene_attributes,\@values);
			}
			if($headers && ref($headers) eq 'ARRAY'){
				my ($junk,@items) = @$headers;
				$header_line_text .= "\t".join("\t",@items);
				$header_line_html .= "<TD><b>".join("<TD><b>",@items);
			}
			print OUT "$header_line_text\n";
			print "$header_line_html\n";

			my @sorted = (0..($nGenes-1));
			if($sort_by>=0){
				if($reverse){
					@sorted = sort {abs($attributes2[$sort_by]->[$a])<=>abs($attributes2[$sort_by]->[$b])} @sorted;
				}else{
					@sorted = sort {abs($attributes2[$sort_by]->[$b])<=>abs($attributes2[$sort_by]->[$a])} @sorted;
				}
			}
			for(my $i1=0; $i1<$nGenes; $i1++){
				my $i = $sorted[$i1];
				print OUT $gene_list[$i];
				if($hyperlink){
					print "<TR><TD><a href=\"JavaScript:plot('$gene_list[$i]');\">$gene_list[$i]</a>";
				}else{
					print "<TR><TD>$gene_list[$i]";
				}
				for(my $j=0; $j<@attributes1; $j++){
					my $x = $attributes1[$j]->[$i];
					print OUT "\t$x";
					print "<TD ALIGN=CENTER>$x";
				}
				for(my $j=0; $j<@attributes2; $j++){
					my $x = $attributes2[$j]->[$i];
					print OUT "\t$x";
					print "<TD ALIGN=CENTER>$x";
					if($j==0 && ($gene_attribute_min==0 || $gene_attribute_min>$x)){
						$gene_attribute_min=$x;
					}
				}
				for(my $j=0; $j<@gene_attributes; $j++){
					my $x = "";
					if($gene_attributes[$j]){
						$x = $gene_attributes[$j]->[$i];
					}
					print OUT "\t$x";
					print "<TD>$x";
				}
				my $ref = $hashGenes{$gene_list[$i]};
				if($ref && ref($ref) eq 'ARRAY'){
					for(my $j=1; $j<@$ref; $j++){
						my $x = $ref->[$j];
						print OUT "\t$x";
						if(length($x)>60){ print "<TD WIDTH=400>$x"; }
						else{ print "<TD>$x"; }
					}
				}
				print OUT "\n";
				print "\n";
			}
			print "</TABLE><p>\n";
			if($iGenesNeg>=0 && $nGenesNeg>0){
				my @gene_attributes_neg;
				for(my $ii=1; $ii<@database_headers; $ii++){
					my @values = split(/,/,$attributes_neg[$ii-1]);
					push(@gene_attributes_neg,\@values);
				}
				my %hashGenes1;
				get_gene_annotation($organismID,\@gene_list_neg,\%hashGenes1,"symbol");
				%hashGenes = (%hashGenes,%hashGenes1);
				print OUT "!Genelist_description\t$output_type_name downregulated: $plot_type\n";
				print "<FONT SIZE=+2>List of genes downregulated</FONT><p>\n";
				print OUT "$header_line_text\n";
				print "$header_line_html\n";
				for(my $i=0; $i<$nGenesNeg; $i++){
					print OUT $gene_list_neg[$i];
					print "<TR><TD><a href=\"JavaScript:plot('$gene_list_neg[$i]');\">$gene_list_neg[$i]</a>";
					for(my $j=0; $j<@gene_attributes_neg; $j++){
						my $x = "";
						if($gene_attributes_neg[$j]){
							$x = $gene_attributes_neg[$j]->[$i];
						}
						print OUT "\t$x";
						print "<TD>$x";
					}
					my $ref = $hashGenes1{$gene_list_neg[$i]};
					if($ref && ref($ref) eq 'ARRAY'){
						for(my $j=1; $j<@$ref; $j++){
							my $x = $ref->[$j];
							print OUT "\t$x";
							if(length($x)>60){ print "<TD WIDTH=400>$x"; }
							else{ print "<TD>$x"; }
						}
					}
					print OUT "\n";
					print "\n";
				}
				print "</TABLE><p>\n";
			}
			print "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
			print "<INPUT NAME=action TYPE=hidden VALUE=matrix_explore1>\n";
			print "<INPUT NAME=analysis TYPE=hidden VALUE=search>\n";
			print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
			print "<INPUT NAME=search_term TYPE=hidden>\n";
			print "<INPUT NAME=category TYPE=hidden VALUE=1>\n";
			print "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$dataset1\">\n";
			print "</FORM><p>\n";
			close OUT;
		}
		#GENERATE: RANKPLOT, SCATTERPLOT, VENN DIAGRAM IF APPLICABLE
		if($hashInput{"graphs"} eq "on"){
			if($data_type1 eq "matrix" && $data_type2 eq "geneset"){
				my $FDR = $hashOutput{"output_fdr1"};
				my $fold_change = $hashOutput{"output_fold_change1"};
				my $expr_thresh = $hashOutput{"output_expr_thresh1"};
				my $fileRankplotID = make_rankplot($dataset1,$dataset2,$col1,$row1,$baseline_name,$FDR,$fold_change,$expr_thresh,$organismID,$organismID1,$gene_attribute_min);
				if($fileRankplotID){
					print "<H3>Rankplot</H3>\n<IMG SRC=../output/$fileRankplotID.gif BORDER=0><p>\n";
				}
			}elsif($data_type1 eq "matrix" && $data_type2 eq "matrix"){
				my $baseline_name1 = $hashOutput{"output_baseline1"};
				my $baseline_name2 = $hashOutput{"output_baseline2"};
				my $FDR1 = $hashOutput{"output_FDR1"};
				my $FDR2 = $hashOutput{"output_FDR2"};
				my $fold_change1 = $hashOutput{"output_fold_change1"};
				my $fold_change2 = $hashOutput{"output_fold_change2"};
				my $expr_thresh1 = $hashOutput{"output_expr_thresh1"};
				my $expr_thresh2 = $hashOutput{"output_expr_thresh2"};
				my $fileScatterplotID = make_scatterplot_twofiles($dataset1,$dataset2,$col1,$row1,$baseline_name1,$baseline_name2,$FDR1,$FDR2,$fold_change1,$fold_change2,$expr_thresh1,$expr_thresh2,$organismID,$organismID1,\%hashGenes);
				if($fileScatterplotID){
					my $file2 = $fileScatterplotID+1;
					print "<H3>Scatterplots</H3>\n<TABLE BORDER=0><TR><TD align=center>Log-ratio of expression<TD align=center>Rank of expression change\n";
					print "<TR><TD><IMG SRC=../output/$fileScatterplotID.gif>\n";
					print "<TD><IMG SRC=../output/$file2.gif></TABLE><p>\n";
				}
			}elsif($data_type1 eq "geneset" && $data_type2 eq "geneset"){
				my $fileVennID = make_venn_diagram($dataset1,$dataset2,$col1,$row1,$organismID,$organismID1,$nGenes);
				if($fileVennID){
					print "<H3>Venn diagram</H3>\n<IMG SRC=../output/$fileVennID.gif BORDER=0><p>\n";
				}
			}
		}
	}
	print "<HR NOSHADE></HR><p>\n";
	print "<INPUT TYPE=button VALUE=\" Close window \" LANGUAGE=\"javascript\" onClick=\"window.close();\">\n";
	print "</BODY></HTML>\n";
	exit(0);
}
return;
}

#********************************
sub  save_output_file
#********************************
{
my $file_output_full=shift;
my $organismID = $hashInput{"organismID"};
my $save_file = $hashInput{"save_file"};
my %hashOutput;
parse_file_headers($file_output_full,\%hashOutput);
my $nMatrix = $hashOutput{"output_n_matrixes"};
if($save_file eq "output"){
	my $file_output_new = $hashInput{"file_output_new"};
	my $description_new = $hashInput{"description_new"};
	open (INFO, "<$file_output_full") or error_message("In save_output_file - cannot open file");
	open (OUT, ">$PATH_DATA/$loginname-$file_output_new") or error_message("In save_output_file - cannot open file");
	print OUT "!Output_name\t$file_output_new\n";
	if($description_new){
		print OUT "!Output_description\t$description_new\n";
	}
	while(my $line=<INFO>){
		if($line =~ /^\!Output_name\t/ || $description_new && $line =~ /^\!Output_description\t/){ next; }
		print OUT $line;
	}
	close OUT;
	my $textline = "type_output=$file_output_new\torganismID=$organismID";
	if($description_new){ $textline .= "\tdescription=$description_new"; }
	$textline .= "\tdate=$date_record";
	file_append("$textline","$PATH_INFO/$loginname-config.txt");
	terminal_window("<H3>Output saved in file '$file_output_new'</H3>\n");
	exit(0);
}
if($save_file eq "geneset"){
	my $matrix_genes = $hashInput{"matrix_genes"};
	if(!$matrix_genes){ error_message("Matrix with genes not passed"); }
	my $file_geneset_new = $hashInput{"file_geneset_new"};
	parse_file_headers($file_output_full,\%hashOutput);
	my $use_database=0;
	my $matrix_name = $hashOutput{"output_matrix$matrix_genes"."_name"};
	my $matrix_type = $hashOutput{"output_matrix$matrix_genes"."_type"};
	if(!$matrix_name){ error_message("Matrix not found"); }
	if(!$matrix_type){ error_message("No matrix type"); }
	if(ref($matrix_type) eq 'ARRAY' && $matrix_type->[1] eq "database"){
		$use_database=1;
	}elsif($matrix_type ne "genes" && $matrix_name !~ /^gene_symbols/i){
		error_message("Wrong matrix type");
	}
	my @database;
	my @database_headers;
	my $database_headers = $hashOutput{"output_database_headers"};
	if($use_database){
		if(!$database_headers){ error_message("In output_explore - no database_headers"); }
		if(ref($database_headers) eq 'ARRAY'){
			@database_headers = @$database_headers;
		}else{
			@database_headers = ($database_headers);
		}
		open(INFO,'<',$file_output_full) or error_message("In output_explore - cannot open file");
		while(my $line = <INFO>){
			if($line =~ /^!database_start/i){ last; }
		}
		while(my $line = <INFO>){
			chop $line;
			if($line =~ /^!/){ last; }
			if($line !~ /^>/ || !$line){ next; }
			my $data_id = $line;
			$data_id =~ s/^>//;
			for(my $ih=0; $ih<@database_headers; $ih++){
				my $header  = $database_headers[$ih];
				my $line1 = <INFO>;
				chop $line1;
				my($key,$value)=split(/\t/,$line1);
				if($key ne $header){ error_message("Database header does not match"); }
				$database[$data_id]->[$ih] = $value;
			}
			if(!$database[$data_id]){ error_message("\'genes\' not found in database ID=$data_id"); }
		}
		close INFO;
	}
	my $text = "!Output_name\t$file_geneset_new\n";
	my $description_new = $hashInput{"description_geneset_new"};
	if($description_new){
		$text .= "!Output_description\t$description_new\n";
	}
	open (OUT, ">$PATH_DATA/$loginname-$file_geneset_new") or error_message("In output_explore - cannot open file_geneset_new");
	print OUT $text;
	open(INFO,'<',$file_output_full) or error_message("In output_explore - cannot open file");
	my $searchline = "!Matrix$matrix_genes"."_start";
	my $found=0;
	while(my $line = <INFO>){
		if($line =~ /^$searchline/){ $found=1; last; }
	}
	if(!$found){ error_message("Matrix w. genes not found"); }
	my $line = <INFO>;
	chop $line;
	my ($junk,@headers)=split(/\t/, $line);
	while(my $line = <INFO>){
		if($line =~ /^!/){ last; }
		chop $line;
		my ($row_header,@data)=split(/\t/, $line);
		for(my $i=0; $i<@data; $i++){
			if(!$data[$i]){ next; }
			my $title = "$headers[$i].$row_header";
			if($matrix_name =~ /_up$|_upregulated$/i){ $title .= "_up"; }
			if($matrix_name =~ /_down$|_downregulated$/i){ $title .= "_down"; }
			if($use_database){
				my @list = split(/,/,$database[$data[$i]]->[0]);
				print OUT join("\t",$title,$title,@list)."\n";
				for(my $ih=1; $ih<@database_headers; $ih++){
					my @list = split(/,/,$database[$data[$i]]->[$ih]);
					print OUT join("\t","",$database_headers[$ih],@list)."\n";
				}
			}else{
				my @symbols = split(/,/,$data[$i]);
				print OUT join("\t",$title,$title,@symbols)."\n";
			}
		}
	}
	close INFO;
	close OUT;
	my $textline = "type_geneset=$file_geneset_new\torganismID=$organismID";
	if($description_new){ $textline .= "\tdescription=$description_new"; }
	$textline .= "\tdate=$date_record";
	file_append("$textline","$PATH_INFO/$loginname-config.txt");
	terminal_window("<H3>Lists of genes saved in file '$file_geneset_new'</H3>\n");
	exit(0);
}
return;
}

#**************************************
sub  matrix_explore
#**************************************
{
if($hashInput{"analysis"} eq "run_anova"){
	run_anova($RUN_ANOVA_PARAM);
}
my $organismID = $hashInput{"organismID"};
my $file_matrix = $hashInput{"file_matrix"};
my $description_matrix = $hashInput{"description_matrix"};
my $file_matrix_full = $file_matrix;
my $file_anova = $file_matrix;
my $legend_file = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
$legend_file =~ s/^public-/public-legend-/;
my $file_anova1 = $file_anova;
my $legendID;
if($file_matrix !~ /^public-/){
	$file_matrix_full = "$loginname-$file_matrix";
	$file_anova = "$loginname-anova-$file_matrix";
	$file_anova1 = "anova-".$file_anova1;
	$legend_file = "$loginname-legend-$file_matrix";
}
if(file_exist("$PATH_DATA/$legend_file")){
	$legendID = get_outputID(1);
	copy "$PATH_DATA/$legend_file", "$PATH_OUTPUT/$legendID.txt";
}
my %hashMatrix=();
my $file_platform = get_array_platform($file_matrix_full,\%hashMatrix);
if($file_platform=~/ /){ error_message($file_platform); }
my $ref = $hashMatrix{"sample_title"};
if(!$ref){ error_message("Matrix file has no sample titles!"); }
my @headers = @$ref;
my $nCol = @headers;
if(!file_exist("$PATH_DATA/$file_anova")){ error_message("ANOVA missing"); }
my @headers1 = ("Symbol/probe");
my ($count_FDR005,$count_FDR1,$nLevel,$nRow,$nSymbols,$nRedundant,$nMissing) = count_significant_genes($file_anova);
if($file_platform ne "None"){
	open(INFO,"<$PATH_DATA/$file_platform") or error_message("Cannot open annotation $file_platform");
	my $line;
	while($line=<INFO>){
		if($line !~ /^!/){ last; }
	}
	chop $line;
	close INFO;
	@headers1 = split(/\t/,$line);
}
my $seriesID = $hashMatrix{"series_geo_accession"};
if(!$seriesID && $file_matrix =~ /GSE/){
	$seriesID = $file_matrix;
	$seriesID =~ s/.*GSE/GSE/;
	$seriesID =~ s/\D+$//;
}
my $platform = $hashMatrix{"series_platform_id"};
#my $nRow=0;
#$ref = $hashMatrix{"sample_data_row_count"};
#if($ref){
#	$nRow = $ref->[0];
#}
my @header_col1 = get_anova_headers($file_anova);
my @header_col;
my $max_repl=1;
for(my $i=0; $i<@header_col1; ++$i){
	if($header_col1[$i]=~/ \(\d+\)$/ && $header_col1[$i] !~ /^Mean/){
		my @items = split(/ \(/,$header_col1[$i]);
		my $n = pop(@items);
		my $name = join(" (",@items);
		$n =~ s/\)$//;
		push(@header_col,[$name,$n]);
		if($max_repl < $n){ $max_repl = $n; }
	}else{
		my $name = $header_col1[$i];
		$name =~ s/^Mean\(//;
		$name =~ s/\)$//;
		push(@header_col,[$name,1]);
	}
}
if(!@header_col){ error_message("No headers in ANOVA file"); }
my @sorted = 0..(@header_col-1);
if(0){
	@sorted = sort {lc($header_col[$a]->[0]) cmp lc($header_col[$b]->[0])} @sorted;
}

print "<HTML><HEAD><TITLE>ExAtlas</TITLE>\n";
print "<STYLE>\n";
print "HH1 {font-size:120%; font-weight:bold}\n";
print "</STYLE>\n";
print_header("update_replications1(); update_replications2();");
print "<script type=text/JavaScript language=JavaScript>\n";
print "<!-- \n";
print "var replication_names = new Array();\n";
print "var replications = new Array();\n";
print "replication_names[0] = \"All\";\n";
for(my $i=1; $i<=$max_repl; ++$i){
	print "replication_names[$i] = $i;\n";
}
for(my $i=0; $i<@sorted; ++$i){
	print "replications[$i] = $header_col[$sorted[$i]]->[1];\n";
}
print "function update_replications1() {\n";
print "	var index = document.form_matrix.select_column.selectedIndex;\n";
print "	document.form_matrix.select_replication1.options.length=0;\n";
print "	for(i=0; i<=replications[index]; ++i){\n";
print "		var x = new Option(replication_names[i],i);\n";
print "		document.form_matrix.select_replication1.options[i] = x;\n";
print "	}\n";
print "}\n";
print "function update_replications2() {\n";
print "	var index = document.form_matrix.compare_column.selectedIndex;\n";
print "	document.form_matrix.select_replication2.options.length=0;\n";
print "	if(index==0){\n";
print "		var x = new Option(\"n/a\",0);\n";
print "		document.form_matrix.select_replication2.options[0] = x;\n";
print "	}else{\n";
print "		for(i=0; i<=replications[index-1]; ++i){\n";
print "			var x = new Option(replication_names[i],i);\n";
print "			document.form_matrix.select_replication2.options[i] = x;\n";
print "		}\n";
print "	}\n";
print "}\n";
print "function search_onsubmit() {\n";
print "	if(!document.form_matrix.search_term.value){\n";
print "		alert(\"You need to put a search term\");\n";
print "		return(false);\n";
print "	}\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_matrix.target = \"_BLANK\"+x;\n";
print "	document.form_matrix.analysis.value = \"search\";\n";
print "	document.form_matrix.file_download.value = \"\";\n";
print "	document.form_matrix.submit();\n";
print "}\n";
print "function do_analysis(type) {\n";
print "	if(type==\"\"){\n";
print "		alert(\"Empty analysis type! Contact webmaster.\");\n";
print "		return(false);\n";
print "	}\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	if(type==\"nonredundant\" && !confirm(\"Do you want to remove redundant probes?\\nThis operation cannot be reversed.\")){\n";
print "		return false;\n";
print "	}\n";
print "	if(type==\"normalize\" && !confirm(\"Do you want to normalize data?\\nThis operation cannot be reversed.\")){\n";
print "		return false;\n";
print "	}\n";
print "	if(type==\"nonredundant\" || type==\"normalize\" || type==\"run_anova\"){\n";
print "		document.form_matrix.target = \"\";\n";
print "	}else{\n";
print "		document.form_matrix.target = \"_BLANK\"+x;\n";
print "	}\n";
print "	if(type==\"anova\"){\n";
print "		document.form_matrix.file_download.value = \"$file_anova1\";\n";
print "	}else if(type==\"raw_data\"){\n";
print "		document.form_matrix.file_download.value = \"$file_matrix\";\n";
print "	}else{\n";
print "		document.form_matrix.file_download.value = \"\";\n";
print "	}\n";
print "	if(type==\"pairwise\" || type==\"meta-analysis\"){\n";
print "		if(document.form_matrix.select_column.selectedIndex==document.form_matrix.compare_column.selectedIndex-1 && document.form_matrix.select_replication1.selectedIndex==document.form_matrix.select_replication2.selectedIndex){\n";
print "			alert(\"You cannot compare expression data to itself\");\n";
print "			return(false);\n";
print "		}\n";
print "	}\n";
print "	if(type==\"meta-analysis\" && document.form_matrix.select_replication1.selectedIndex + document.form_matrix.select_replication2.selectedIndex>0){\n";
print "		alert(\"Select 'All' replications for meta-analysis\");\n";
print "		return(false);\n";
print "	}\n";
print "	document.form_matrix.analysis.value = type;\n";
print "	document.form_matrix.submit();\n";
print "}\n";
print "<!-- end script --></SCRIPT></HEAD>\n";
print "<H2>Expression profiles table: $file_matrix</H2>\n";
if(!$nSymbols){
	print "<b><font color=red>WARNING:</font></b> Gene symbols not found. See notes below<br>\n";
}
if($count_FDR005<5 || $count_FDR1 < 50){
	print "<b><font color=red>WARNING:</font></b>The number of significant probes (N=$count_FDR005) is too low! Possible cause: values are too noisy or no difference between cell types.<br>\n";
}
if($nMissing){
	print "<b><font color=red>WARNING:</font></b> Some probes (n=$nMissing) have missing names in the 1st column.<br>\n";
}
if($nRedundant){
	print "<b><font color=red>WARNING:</font></b> Some probes (n=$nRedundant) are redundant. Information can be lost<br>\n";
}
print "<b>Organism:</b> $hashOrganism{$organismID} &nbsp; &nbsp;\n";
if($seriesID =~ /GSE\d+$/ && length($seriesID) <=12){
	my $x = int(10000*rand());
	print "<b>GEO series:</b> <a href=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$seriesID target=_BLANK$x>$seriesID</a> &nbsp; &nbsp;\n";
}
print "<b>Platform:</b> \n";
if($platform =~/^GPL\d+$/){
	my $x = int(10000*rand());
	print "<a href=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$platform target=_BLANK$x>$platform</a>\n";
}elsif($platform == $organismID || $platform eq "symbol_$organismID"){
	print "Gene symbols\n";
}else{
	print "$platform\n";
}
my $nLevels = @header_col;
print " &nbsp; &nbsp; <b>N columns:</b> $nCol &nbsp; &nbsp;\n";
print "<b>N sample types:</b> $nLevels &nbsp; &nbsp;\n";
print "<b>N rows:</b> $nRow &nbsp; &nbsp; <b>Probes with gene symbols:</b> $nSymbols\n";
if($description_matrix){
	my $description1 = add_hyperlinks($description_matrix);
	print "<br><b>Description:</b> $description1\n";
}
if($legendID){
	my $x = int(10000*rand());
	print "<br><b>Legend:</b> <a href=$HOME_ADDRESS/output/$legendID.txt target=_blank$x>Legend file</a>\n";
}
print "<br><b>N significant genes</b> (FDR<0.05) = $count_FDR005\n";

print "<p><FORM NAME=form_matrix ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
my @filter_foldChange = (1,1.1,1.2,1.5,2,3,4,5,7,10,15,20);
my @filter_FDR = (1,0.9,0.5,0.2,0.1,0.05,0.01,0.001,0.0001);
my @cluster_correlation = (0.2,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9);
my @cluster_fold_thresh = (1,1.1,1.2,1.5,2,3,4,5,7,10,15,20);
my @filter_expr = (0,0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,100,300,1000,3000);
print "<TABLE BORDER=0><TR><TD WIDTH=250 VALIGN=center><HH1>1. Make heatmap or PCA</HH1>\n";
print "<TD VALIGN=center><INPUT TYPE=button VALUE=\"Make heatmap\" style=width:150px; onClick=\"do_analysis('heatmap');\">\n";
print "<TD VALIGN=center><INPUT TYPE=button VALUE=\"PCA\" style=width:150px; onClick=\"do_analysis('pca');\">\n";
print "<TD><INPUT TYPE=checkbox NAME=use_replications><br><INPUT TYPE=checkbox NAME=pca_cluster>\n";
print "<TD>Show replications<br>PC gene clusters\n";
print "</TABLE><TABLE BORDER=0>\n";
print "<TR><TD WIDTH=70><select name=filter_FDR style=width:60px;>\n";
foreach my $x (@filter_FDR){
	print "<option value=$x"; if($x==0.05){ print " selected"; } print "> $x\n";
}
print "</select><TD WIDTH=180>FDR threshold\n";
print "<TD><select name=matrix_filter style=width:240px;>\n";
print "<option value=none> No filtering\n";
print "<option value=rows> Filter rows only\n";
print "<option value=rows_and_columns selected> Filter rows and columns\n";
print "</select> Filtering\n";
print "<TR><TD><select name=filter_foldChange style=width:60px;>\n";
foreach my $x (@filter_foldChange){
	print "<option value=$x"; if($x==2){ print " selected"; } print "> $x\n";
}
print "</select><TD>Fold change threshold\n";
print "<TD><select name=matrix_cluster style=width:240px;>\n";
print "<option value=3> Hierarchical clustering\n";
print "<option value=4> Diagonal clustering\n";
print "<option value=2> Cluster columns only\n";
print "<option value=1> Cluster rows only\n";
print "<option value=0> No clustering\n";
print "</select> Heatmap clustering\n";
print "<TR><TD><select name=cluster_correlation style=width:60px;>\n";
foreach my $x (@cluster_correlation){
	print "<option value=$x"; if($x==0.7){ print " selected"; } print "> $x\n";
}
print "</select><TD>Correlation (PC clustring)\n";
print "<TD><select name=cluster_fold_thresh style=width:240px;>\n";
foreach my $x (@cluster_fold_thresh){
	print "<option value=$x"; if($x==2){ print " selected"; } print "> $x\n";
}
print "</select> Fold change (PC clustring)\n";
print "</TABLE><p>\n";

print "<TABLE BORDER=0><TR><TD WIDTH=250><HH1>2. Pairwise comparison</HH1>\n";
print "<TD><INPUT TYPE=button VALUE=\"Scatter-plot\" style=width:150px; onClick=\"do_analysis('pairwise');\">\n";
print "<TD><INPUT TYPE=button VALUE=\"Meta-analysis\" style=width:150px; onClick=\"do_analysis('meta-analysis');\">\n";
print "<TD><INPUT TYPE=checkbox NAME=useSymbols> Use non-redundant symbols\n";
print "</TABLE><TABLE BORDER=0>\n";
print "<TR><TD WIDTH=120>Select data<TD><select name=select_column style=width:250px; onChange=update_replications1();>\n";
my $nn = @header_col;
for(my $i=0; $i<$nn; ++$i){
	my $x=$i+1;
	print "<option value=$x> $header_col[$sorted[$i]]->[0]\n";
}
print "</select>\n";
print "<TD><select name=select_replication1 style=width:50px;></select>\n";

print "<TD WIDTH=50>repl.<TD><select name=FDR_pairwise style=width:120px;>\n";
foreach my $x (@filter_FDR){
	print "<option value=$x"; if($x==0.05){ print " selected"; } print "> $x\n";
}
#print "</select><TD> FDR threshold\n";
print "</select><TD> <select name=FDR_pvalue><option value=FDR>FDR<option value=p-value>p-value</select> threshold\n";
print "<TR><TD>Compare with:<TD><select name=compare_column style=width:250px; onChange=update_replications2();>\n";
print "<option value=0> Median profile\n";
for(my $i=0; $i<$nn; ++$i){
	my $x=$i+1;
	print "<option value=$x> $header_col[$sorted[$i]]->[0]\n";
}
print "</select>\n";
print "<TD><select name=select_replication2 style=width:50px;></select>\n";

print "<TD WIDTH=50>repl.<TD><select name=foldChange_pairwise style=width:120px;>\n";
foreach my $x (@filter_foldChange){
	print "<option value=$x"; if($x==2){ print " selected"; } print "> $x\n";
}
print "</select><TD>Fold change threshold\n";
print "<TR><TD><TD><TD><TD><TD><select name=exprThreshold_pairwise style=width:120px;>\n";
foreach my $x (@filter_expr){ print "<option value=$x> $x\n"; }
print "</select><TD>Expression threshold\n";
print "</TABLE><p>\n";

print "<TABLE BORDER=0><TR><TD WIDTH=250><HH1>3. Find genes</HH1>\n";
print "<TD><INPUT type=button value=\"Search\" style=width:150px; onClick=search_onsubmit();>\n";
print "</TABLE><TABLE BORDER=0>\n";
print "<TR><TD>Search term:<TD><TD><INPUT NAME=search_term style=width:200px;>\n";
print "<TD><select name=category style=width:150px;>\n";
print " <option value=-1>All categories\n";
for(my $i=1; $i<@headers1; $i++){
	print " <option value=$i> $headers1[$i]\n";
}
print " <option value=0> $headers1[0]\n";
print "</select>\n";
print "<TD><INPUT TYPE=checkbox NAME=sort_histogram> Sort expression values<TD>\n";
print "</TABLE><p>\n";

print "<HH1>4. Other functions</HH1><br>\n";
print "<INPUT TYPE=button VALUE=\"Correlation ...\" style=width:200px; onClick=\"do_analysis('correlation');\"";
if(!$nSymbols){ print "disabled"; }
print "> Correlation with other gene expression profiles<br>\n";
print "<INPUT TYPE=button VALUE=\"Geneset analysis ...\" style=width:200px; onClick=\"do_analysis('geneset_analysis');\"";
if(!$nSymbols){ print "disabled"; }
print "> Gene set enrichment analysis (e.g., Gene Ontology)<br>\n";
print "<INPUT TYPE=button VALUE=\"Signiticant genes ...\" style=width:200px; onClick=\"do_analysis('matrix2geneset');\"> Generate sets of significant/specific genes<br>\n";
print "<INPUT TYPE=button VALUE=\"Data quality\" style=width:200px; onClick=\"do_analysis('quality');\"> Evaluate data consistency<br>\n";
print "<INPUT TYPE=button VALUE=\"Get ANOVA output\" style=width:200px; onClick=\"do_analysis('anova');\"> Download ANOVA results<br>\n";
print "<INPUT TYPE=button VALUE=\"Run ANOVA again...\" style=width:200px; onClick=\"do_analysis('run_anova');\"";
if($file_matrix =~/^public-/){ print "disabled"; }
print "> Run ANOVA with different parameters<br>\n";
print "<INPUT TYPE=button VALUE=\"Get raw data\" style=width:200px; onClick=\"do_analysis('raw_data');\"> Download raw data (with symbols)<br>\n";
if($file_matrix !~ /^public-/ || $loginname eq "public"){
	if(!$hashMatrix{"series_normalized"}){
		print "<INPUT TYPE=button VALUE=\"Normalize\" style=width:200px; onClick=\"do_analysis('normalize');\"> Normalize using quantile method<br>\n";
	}
	if(!$hashMatrix{"series_nonredundant"}){
		print "<INPUT TYPE=button VALUE=\"Best probe\" style=width:200px; onClick=\"do_analysis('nonredundant');\"> Leave best probe for each gene (remove redundancy)<br>\n";
	}
}
print "<p><b>Notes:</b><br>\n";
if(!$nSymbols){ 
	print "Gene symbols are missing in the data set; thus many functions of ExAtlas are disabled. If you used gene symbols as probe identifiers (i.e., in the first column) then you need to select platform \"Gene symbols\".\n";
	print "Alternatively, use platform that connects probe identifiers with gene symbols. Platforn can be changed using \"Edit\" button in File management section of the main menu.\n";
}
print "FDR = False Discovery Rate. It is used instead of p-values to account for multiple hypotheses testing<br>\n";
print "PCA = Principal Component Analysis (spatial representation of data). PC clustering = identification of genes whose expression profile correlates with principal components.<p>\n";
print "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=hidden VALUE=\"matrix_explore1\">\n";
print "<INPUT NAME=\"analysis\" TYPE=hidden>\n";
print "<INPUT NAME=\"file_matrix\" TYPE=hidden VALUE=\"$file_matrix\">\n";
print "<INPUT NAME=\"description_matrix\" TYPE=hidden VALUE=\"$description_matrix\">\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT NAME=file_download TYPE=hidden VALUE=\"\">\n";
if($legendID){ print "<INPUT NAME=legendID TYPE=hidden VALUE=$legendID>\n"; }
print "</FORM>\n";
print "<HR NOSHADE></HR>\n";
print "<INPUT TYPE=button VALUE=\" Cancel (close window) \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub   get_anova_headers
#**************************************
{
my $file_anova = shift;
my $remove_repl = shift;

open(INFO,"<$PATH_DATA/$file_anova") or return;
my $line = <INFO>;
chop $line;
my ($junk,$junk1,@header_col) =split(/\t/,$line);
if(!@header_col){ return; }
my $nLevel=0;
while($nLevel<@header_col && $header_col[$nLevel] !~ /^Var\(/){ ++$nLevel; }
if($remove_repl){
	for(my $i=0; $i<$nLevel; ++$i){
		$header_col[$i] =~ s/ \(\d+\)$//;
	}
}
splice(@header_col,$nLevel);
close INFO;
return (@header_col);
}

#**************************************
sub   count_significant_genes
#**************************************
{
my $file_anova = shift;
open(INFO,"<$PATH_DATA/$file_anova") or error_message("ANOVA file is missing");
my $line = <INFO>;
chop $line;
my ($junk,$junk1,@header_col) =split(/\t/,$line);
if(!@header_col){ error_message("Invalid ANOVA file!"); }
my $nLevel=0;
while($nLevel<@header_col && $header_col[$nLevel] !~ /^Var\(/){ ++$nLevel; }
my ($count_FDR005,$count_FDR1,$nrows,$nSymbols,$nRedundant,$nMissing)=(0,0,0,0,0,0);
my %hashProbe;
while(my $line = <INFO>){
	chop $line;
	my ($id,$expr,@data) =split(/\t/, $line);
	if(!$id){ $nMissing++; }
	else{ $hashProbe{"$id"}=1; }
	$nrows++;
	my $FDR = $data[$nLevel+6];
	if($FDR<=0.05){ $count_FDR005++; }
	if($FDR<0.95){ $count_FDR1++; }
	if($data[$nLevel+8]){ $nSymbols++; }
}
close INFO;
my $n1 = keys %hashProbe;
my $nRedundant = $nrows-$nMissing-$n1;
return ($count_FDR005,$count_FDR1,$nLevel,$nrows,$nSymbols,$nRedundant,$nMissing);
}

#**************************************
sub   run_anova
#**************************************
{
my $runID = shift;
$hashInput{"runID"}=$runID;
my $file_matrix = $hashInput{"file_matrix"};
if($file_matrix =~ /^public-/){ return; }
my $file_matrix_full = "$loginname-$file_matrix";
my $file_anova = "$loginname-anova-$file_matrix";
if($runID==$RUN_ANOVA_FIRST && file_exist("$PATH_DATA/$file_anova")){ return; }
my %hashMatrix=();
my $file_platform = get_array_platform($file_matrix_full,\%hashMatrix);
my $ref = $hashMatrix{"sample_title"};
if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Matrix file has no sample titles!"); }
my $nCol = @$ref;
if($nCol < 200 || $runID==$RUN_ANOVA_NORM && $nCol < 80){
	run_anova1($file_matrix);
	return;
}
my $logFileID = get_outputID(1);
$hashInput{"logFileID"}=$logFileID;
interrupt_program(1);
}

#**************************************
sub   run_anova1
#**************************************
{
my $file_matrix = shift;
my $logFileID = shift;
my $err_model = shift;

my $file_matrix_full = $file_matrix;
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_matrix_full = "$loginname-$file_matrix";
	$file_anova = "$loginname-anova-$file_matrix";
}
my %hashMatrix=();
my $file_platform = get_array_platform($file_matrix_full,\%hashMatrix);
if($file_platform=~/ /){ error_message($file_platform); }
my $nRow=0;
my $ref = $hashMatrix{"sample_data_row_count"};
if($ref){ $nRow = $ref->[0]; }
my $runID = $hashInput{"runID"};
if($runID==$RUN_ANOVA_NORM){
	my $fileID = get_outputID(1);
	my $response = system("$PATH_BIN/norm_new","$PATH_DATA/$file_matrix_full","$PATH_OUTPUT/$fileID.txt");
	if($response){ error_message("Normalization crashed!",$logFileID); }
	copy "$PATH_OUTPUT/$fileID.txt", "$PATH_DATA/$file_matrix_full";
}
my @command = ("$PATH_BIN/anova_oneway","-i","$PATH_DATA/$file_matrix_full","-o","$PATH_DATA/$file_anova");
if($file_platform eq "None"){ push(@command,"-probe","2"); }
else{ push(@command,"-a","$PATH_DATA/$file_platform"); }
if($runID == $RUN_ANOVA_PARAM){
	my $cutoff = $hashInput{"cutoff"}+0; if($cutoff<=0){ $cutoff=0; }
	my $z_outliers = $hashInput{"z_outliers"}+0; if($z_outliers<=0){ $z_outliers=0; }
	my $prop_var = $hashInput{"prop_var"}+0; if($prop_var<=0){ $prop_var=0; }
	my $window = int($hashInput{"window_width"}); if($window > $nRow){ $window = $nRow; }
	my $error_model = int($hashInput{"error_model"}); if($error_model<=0||$error_model>5){ $error_model=4; }
	my $df_bayesian = int($hashInput{"df_bayesian"}); if($df_bayesian<=0||$df_bayesian>10000){ $df_bayesian=10; }
	my $use_ID = int($hashInput{"use_probeID"}); if($use_ID>2){ $use_ID=2; } if($use_ID<0){ $use_ID=0; }
	if($use_ID){ push(@command,"-probe","$use_ID"); }
	push(@command,"-cutoff","$cutoff","-z","$z_outliers","-prop","$prop_var","-win","$window","-err","$error_model","-bayes","$df_bayesian");
}elsif($err_model){
	push(@command,"-err",$err_model);
}elsif($nRow < 500 || @$ref > 30){
	push(@command,"-err",1);
}
#print "A1 @command<br>\n";
my $response = system(@command);
my $text = check_anova_output($file_matrix);
if($text || $response){
	error_message("ANOVA crashed!\n$text",$logFileID);
}
return;
}

#**************************************
sub  check_anova_output
#**************************************
{
my $file_matrix = shift;
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_anova = "$loginname-anova-$file_matrix";
}
if(!open(INFO,"<$PATH_DATA/$file_anova")){
	return("No output");
}
my $line = <INFO>;
chop $line;
my $text=$line;
if($line=~ /error/i){
	while($line = <INFO>){ $text.=$line; }
	close INFO;
	return($text);
}
close INFO;
return 0;
}

#**************************************
sub   interrupt_program
#**************************************
{
my $condition = shift;
my $data = shift;

my $logFileID = $hashInput{"logFileID"};
my $runID = $hashInput{"runID"};
my $organismID = $hashInput{"organismID"};
my $file_info;
#print "interrupt_program: LogFileID = $logFileID<br>\n";
#print "interrupt_program: runID = $runID<br>\n";
if($logFileID && open(INFO,"<$PATH_OUTPUT/$logFileID.txt")){
	my $response="";
	my %hash;
	while(my $line=<INFO>){
		if($runID==$RUN_MATRIX_UPLOAD && $line eq "File_info\n"){
			$file_info = 99;
			next;
		}
		if($file_info && $runID==$RUN_MATRIX_UPLOAD){
			if($file_info==99){ $file_info=$line; }
			else{ $file_info .= $line; }
		}
		$response .= $line;
		if($line =~ /\t/){
			chop $line;
			my @items = split(/\t/,$line);
			$hash{$items[0]}=$items[1];
		}
	}
	close INFO;
	if($response =~ /error/i){
		if($runID == $RUN_MATRIX_UPLOAD || $runID == $RUN_NORMALIZE){
			error_message("Task stopped!\n<b>Log information:</b>\n$response","continue");
		}else{
			terminal_window("<h3>Errors in task. Task stopped!</h3><b>Log information:</b><br><pre>\n$response\n</pre>","continue");
		}
	}
	if($response =~ /Task completed/i){
		if($runID == $RUN_SIGNIFICANT){
			geneset_explore();
		}elsif($runID==$RUN_OVERLAP || $runID==$RUN_PAGE || $runID==$RUN_CORRELATION){
			output_explore($logFileID+1);
		}elsif($runID == $RUN_NORMALIZE){
			terminal_window("<H3>Expression profile data are combined</H3>","continue");
		}elsif($runID == $RUN_MATRIX_UPLOAD){
			terminal_window("<H3>File is uploaded</H3>$file_info","continue");
		}elsif($runID == $RUN_ANOVA_FIRST || $runID == $RUN_ANOVA_PARAM || $runID == $RUN_ANOVA_NORM){
			matrix_explore();
		}
		print_web_page($logFileID+1);
	}
	my $pid = $hashInput{"process_id"};
	if(!$pid){ $pid = $hash{"process_id"}; }
	print "<HTML><HEAD><TITLE>ExAtlas - interrupt</TITLE>\n";
	print_header();
	print "<h3>The task is not finished yet..</h3>\n";
	print "To see if it is finished, click on \"Check your task\" button<p>\n";
	my $x = int(10000*rand());
	print "The status of your task can be checked here: <a href=$HOME_ADDRESS/output/$logFileID.txt target=_blank$x>Log file</a><p>\n";
	print "<FORM ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
	print "<TABLE BORDER=0>\n";
	print "<TR><TD><INPUT TYPE=submit VALUE=\"  Check your task  \">\n";
	print "<TD WIDTH=100><TD><INPUT NAME=terminate_task TYPE=submit VALUE=\"  Cancel the task  \">\n";
		print "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
	print "<INPUT NAME=\"process_id\" TYPE=hidden VALUE=\"$pid\">\n";
	print "<INPUT NAME=logFileID TYPE=hidden VALUE=\"$logFileID\">\n";
	print "<INPUT NAME=runID TYPE=hidden VALUE=\"$runID\">\n";
	print "<INPUT NAME=action TYPE=hidden VALUE=interrupt_program>\n";
	print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
	print "<INPUT NAME=runID TYPE=hidden VALUE=$runID>\n";
	if($runID == $RUN_SIGNIFICANT){
		my $file_geneset = $hashInput{"file_geneset"};
		print "<INPUT NAME=file_geneset TYPE=hidden VALUE=\"$file_geneset\">\n";
	}elsif($runID==$RUN_ANOVA_FIRST || $runID==$RUN_ANOVA_PARAM || $runID==$RUN_ANOVA_NORM){
		my $file_matrix = $hashInput{"file_matrix"};
		print "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
	}
	print "</TABLE></FORM><p><HR NOSHADE>\n";
	print "</BODY>\n";
	print "</HTML>\n";
	exit(0);
}
if($logFileID){ file_append("Task started","$PATH_OUTPUT/$logFileID.txt",1); }
my $logFileID1;
if($condition >= 1){
	$logFileID1 = $logFileID;
	my $pid = fork();
	if($pid){
		$|++;
		$pid = pid_record($pid);
		print "<HTML><HEAD><TITLE>ExAtlas - interrupt</TITLE>\n";
		print_header();
		print "<h3>The task has started!</h3>\n";
		print "To see if it is finished, click on \"Check your task\" button<p>\n";
		my $x = int(10000*rand());
		print "The status of your task can be checked here: <a href=$HOME_ADDRESS/output/$logFileID.txt target=_blank$x>Log file</a><p>\n";
		print "<FORM ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
		print "<TABLE BORDER=0>\n";
		print "<TR><TD><INPUT TYPE=submit VALUE=\"  Check your task  \">\n";
		print "<TD WIDTH=100><TD><INPUT NAME=terminate_task TYPE=submit VALUE=\"  Cancel the task  \">\n";
		print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
		print "<INPUT NAME=process_id TYPE=hidden VALUE=\"$pid\">\n";
		print "<INPUT NAME=logFileID TYPE=hidden VALUE=\"$logFileID\">\n";
		print "<INPUT NAME=action TYPE=hidden VALUE=interrupt_program>\n";
		print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
		print "<INPUT NAME=runID TYPE=hidden VALUE=$runID>\n";
		if($runID == $RUN_SIGNIFICANT){
			my $file_geneset = $hashInput{"file_geneset"};
			print "<INPUT NAME=file_geneset TYPE=hidden VALUE=\"$file_geneset\">\n";
		}elsif($runID==$RUN_ANOVA_FIRST || $runID==$RUN_ANOVA_PARAM || $runID==$RUN_ANOVA_NORM){
			my $file_matrix = $hashInput{"file_matrix"};
			print "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
		}
		print "</TABLE></FORM><p><HR NOSHADE>\n";
		print "</BODY>\n";
		print "</HTML>\n";
		exit(0);
	}
}else{
	$logFileID1 = "";
}
# Child process
if($runID == $RUN_QUALITY){
	matrix_quality1($logFileID1);
}elsif($runID == $RUN_OVERLAP){
	geneset_overlap1($logFileID1);
}elsif($runID == $RUN_SIGNIFICANT){
	generate_significant_genesets($logFileID1);
}elsif($runID == $RUN_NORMALIZE){
	file_edit3($logFileID1);
}elsif($runID == $RUN_MATRIX_UPLOAD){
	$file_info = check_matrix_data1($logFileID1,$data);
	finish_file_upload($data->[0],$data->[1],"matrix",$organismID);
}elsif($runID == $RUN_PAGE){
	geneset_analysis2($logFileID1);
}elsif($runID == $RUN_CORRELATION){
	correlation2($logFileID1);
}elsif($runID==$RUN_ANOVA_FIRST || $runID==$RUN_ANOVA_PARAM || $runID==$RUN_ANOVA_NORM){
	my $file_matrix = $hashInput{"file_matrix"};
	run_anova1($file_matrix,$logFileID1);
}elsif($runID == $RUN_METAANALYSIS){
	meta_compute1($logFileID1);
}elsif($runID == $RUN_PCA_PCA){
	make_PCA_output_page($logFileID1,$data);
}elsif($runID == $RUN_PCA_CLUSTER){
	my ($fileID,$nCol,$nRow,$matrix_cluster,$diagonal,$logtransform,$matrix_name)=@$data;
	cluster_matrix($fileID,$nCol,$nRow,$matrix_cluster,$diagonal,$logFileID1);
	if($logtransform){ exp_transform($fileID,1); }
	plot_matrix($fileID,$matrix_name,$logFileID1);
}else{
	error_message("Wrong runID in interrupt_program: $runID");
}
if($logFileID1){ file_append("Task completed","$PATH_OUTPUT/$logFileID1.txt"); }
else{
	if($runID == $RUN_SIGNIFICANT){
		geneset_explore();
	}elsif($runID == $RUN_NORMALIZE){
		terminal_window("<H3>Expression profile data are combined</H3>","continue");
	}elsif($runID == $RUN_MATRIX_UPLOAD){
		terminal_window("<H3>File $data->[1] is uploaded</H3>$file_info","continue");
	}
	print_web_page($logFileID+1);
}
exit(0);
}

#**************************************
sub  matrix_explore1
#**************************************
{
my $analysis_type = $hashInput{"analysis"};
my $search_term = $hashInput{"search_term"};
if($analysis_type eq "matrix2geneset"){
	matrix2geneset();
}elsif($analysis_type eq "correlation"){
	correlation();
}elsif($analysis_type eq "run_anova"){
	anova_parameters();
}elsif($analysis_type eq "geneset_analysis"){
	geneset_analysis();
}elsif($analysis_type eq "normalize"){
	run_anova($RUN_ANOVA_NORM);
	matrix_explore();
}elsif($analysis_type eq "nonredundant"){
	nonredundant_matrix();
}elsif($analysis_type eq "quality"){
	matrix_quality();
}elsif($analysis_type eq "meta-analysis"){
	meta_analysis_combine();
}elsif($analysis_type eq "search" && $search_term =~ /^replications_/){
	show_replications();
}

my $replication1 = $hashInput{"select_replication1"};
my $replication2 = $hashInput{"select_replication2"};
if($analysis_type eq "pairwise" && ($replication1>0 || $replication2>0)){
	pairwise_replications();
}
if($analysis_type eq "pairwise"){
	pairwise_compare();
}
my $file_matrix = $hashInput{"file_matrix"};
my $category = $hashInput{"category"};
my $matrix_filter = $hashInput{"matrix_filter"};
my $matrix_cluster = $hashInput{"matrix_cluster"};
my $organismID = $hashInput{"organismID"};
my $legendID = $hashInput{"legendID"};
my $diagonal = 0;
if($matrix_cluster==4){
	$matrix_cluster=3;
	$diagonal=1;
}
my $file_matrix_full = $file_matrix;
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_matrix_full = "$loginname-$file_matrix";
	$file_anova = "$loginname-anova-$file_matrix";
}
my $clusterPC = $hashInput{"clusterPC"};
my $clusterDir = $hashInput{"clusterDir"};
my @cluster_list;
my %hashClusterRow;
if($analysis_type eq "pca_cluster"){
	my $cluster_file = $hashInput{"cluster_file"};
	open(INFO,"<$PATH_OUTPUT/$cluster_file") or error_message("Cannot open $cluster_file");
	my $line = <INFO>;
	while(my $line = <INFO>){
		chop $line;
		my ($id,$logChange,$r,$PC,$dir) =split(/\t/, $line);
		if($id =~ / \(/){
			$id =~ s/ \(.+$//;
		}
		if($clusterPC eq "all" || $PC eq "PC$clusterPC" && $dir eq $clusterDir){
			push(@cluster_list,[$id,$logChange,$r,$PC,$dir]);
		}
	}
	close INFO;
	if(!@cluster_list){
		error_message("No genes found in this cluster!");
	}
	if($clusterDir eq "positive"){
		@cluster_list = sort {$b->[1]<=>$a->[1]} @cluster_list;
	}elsif($clusterDir eq "negative"){
		@cluster_list = sort {$a->[1]<=>$b->[1]} @cluster_list;
	}elsif($clusterDir ne "all"){
		error_message("Wrong direction: $clusterDir");
	}
	for(my $i=0; $i<@cluster_list; ++$i){
		my $ref = $cluster_list[$i];
		$hashClusterRow{$ref->[0]} = $i+1;
	}
}

my $similar_genes_id;
my @search_pattern;
if($analysis_type eq "search" && $search_term =~ /^similar_genes_/){
	$similar_genes_id = $search_term;
	$similar_genes_id =~ s/^similar_genes_//;
	open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file");
	while(my $line = <INFO>){
		chop $line;
		my ($id,$expr,@data1) =split(/\t/, $line);
		if($id eq $similar_genes_id){
			@search_pattern = @data1;
			last;
		}
	}
	close INFO;
}

#THIS SECTION IS USED BOTH FOR GENE SEARCH, HEATMAP
open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file");
my $line = <INFO>;
chop $line;
my ($junk,$junk1,@headers) =split(/\t/,$line);
if(!@headers){ error_message("Headers not found in $file_anova"); }
my @n_repl;
my $nCol=0;
my $nRep=0;
while($nCol<@headers && $headers[$nCol] !~ /^Var\(/){ ++$nCol; }
for(my $i=0; $i<$nCol; ++$i){
	$n_repl[$i] = 1;
	if($headers[$i]=~/ \(\d+\)$/ && $headers[$i] !~ /^Mean/){
		my @items = split(/ \(/,$headers[$i]);
		my $n = pop(@items);
		my $name = join(" (",@items);
		$n =~ s/\)$//;
		$nRep += $n;
		$n_repl[$i] = $n;
		if($n_repl[$i] < 1){ $n_repl[$i]=1; }
		$headers[$i] = $name;
	}else{
		$headers[$i] =~ s/^Mean\(//;
		$headers[$i] =~ s/\)$//;
		$nRep++;
	}
}

my $pattern_scope;
if(@search_pattern){
	splice(@search_pattern,$nCol);
	my ($min,$max)=(1000000,-1000000);
	for(my $i=0; $i<$nCol; ++$i){
		my $x = $search_pattern[$i];
		if($x==$MISSING){ next; }
		if($min>$x){ $min=$x; }
		if($max<$x){ $max=$x; }
	}
	$pattern_scope = $max-$min;
}
my @annotations = @headers;
splice(@headers,$nCol);
splice(@annotations,0,$nCol+8);
my $searchTerm;
my @aliasTerm;
my $annot_search;
if($analysis_type eq "search"){
	$search_term =~ s/<//g;	#destroy hyperlinks
	$search_term =~ s/^\s+//;
	$search_term =~ s/\s+$//;
	if($category>0){ $annot_search = $annotations[$category-1]; }
	if(($category==0 || $annot_search=~/probeid/i) && $search_term =~ /^Z\d\d\d\d\d\d\d\d$/){
		$search_term .= "-1";
	}
	$searchTerm = quotemeta($search_term);
	#print "$search_term $searchTerm<br>\n";
	if($category >0 && $annot_search =~ /^genbank|^refseq/i){
		$search_term =~ s/\.\d+$//;
		open (INFO1, "$PATH_DATA/genbank_$organismID"."_annot.txt");
		while(my $line = <INFO1>){
			chop $line;
			my ($acc,$symbol,$geneName) = split(/\t/,$line);
			$acc =~ s/\.\d+$//;
			if($acc =~ /^$searchTerm$/i){
				push(@aliasTerm,quotemeta($symbol));
				last;
			}
		}
		close INFO1;
	}elsif($category==-1 || $category==1){
		#$search_term =~ s/ .+$|\(.+$//;
		my $found=0;
		open (INFO1, "$PATH_DATA/gene_info_$organismID.txt");
		while(my $line = <INFO1>){
			chop $line;
			my ($org_id,$entrez_id,$symbol,$junk1,$alias_list,$altName,$junk3,$junk4,$geneName) = split(/\t/,$line);
			if($symbol =~ /^$searchTerm$/i){
				$found = 1;
				@aliasTerm=();
				last;
			}
			if($alias_list !~ /$searchTerm/i && $alias_list =~ /-/){
				$alias_list =~ s/-//g;
			}
			if($alias_list =~ /$searchTerm/i){
				foreach my $alias (split(/\|/,$alias_list)){
					if($alias =~ /^$searchTerm$/i){
						push(@aliasTerm,quotemeta($symbol)); last;
					}
				}
			}
			elsif($altName =~ /$searchTerm/i){
				foreach my $expr (split(/\|/,$altName)){
					if($expr !~ /$searchTerm/i){ next; }
					$expr =~ s/^.+://;
					if($expr =~ /^$searchTerm$/i){
						push(@aliasTerm,quotemeta($symbol)); last;
					}
				}
			}
		}
		close INFO1;
	}
}
my @data_list;
my @optional_list;
my @alias_list;
my @data;
my @rows;
my @symbols;
my $irow=0;
my $log10 = log(10);
my $logratio_threshold = $hashInput{"filter_foldChange"};
if($logratio_threshold>1){
	$logratio_threshold = log($logratio_threshold)/$log10;
}else{
	$logratio_threshold = 0;
}
my $correlation = $hashInput{"correlation"};
my $fold_change_search = $hashInput{"fold_change_search"};
my $log_fold_change_search = 0;
if($fold_change_search){
	$log_fold_change_search = log($fold_change_search)/$log10;
}
my @logratio_col;
my $FDRthresh = $hashInput{"filter_FDR"};
my $filter_rows = 0;
if($matrix_filter =~ /rows/){ $filter_rows = 1; }
my $filter_cols = 0;
if($matrix_filter =~ /columns/){ $filter_cols = 1; }
my $symbol_found = 0;
my @average;

while(my $line = <INFO>){
	chop $line;
	$line =~ s/\s+$//;
	my ($id,$expr,@data1) =split(/\t/, $line);
	my $FDR = $data1[$nCol+6];
	my $symbol = $data1[$nCol+8];
	my $gene_name = $data1[$nCol+9];
	my $F = $data1[$nCol+4];
	my $MSE = $data1[$nCol+3];
	if($analysis_type eq "pca_cluster"){  ######## USED FOR PC CLUSTER
		my $ii = $hashClusterRow{$id};
		if($ii){
			push(@{$cluster_list[$ii-1]},$symbol,$gene_name);
		}
	}
	elsif($analysis_type =~ /^heatmap$|^pca$/){  ######## USED FOR HEATMAP, PCA
		if($filter_rows && $FDR > $FDRthresh){ next; }
		splice(@data1,$nCol);
		if($filter_rows && $logratio_threshold>0){
			my @sorted = sort {$a<=>$b} @data1;
			while($sorted[0]<=$MISSING){ shift(@sorted); }
			if(@sorted>1){
				my $logratio = $sorted[@sorted-1] - $sorted[0];
				if($logratio < $logratio_threshold){ next; }
			}
		}
		for(my $i=0; $i<$nCol; ++$i){
			my $x = $data1[$i];
			if($x > $MISSING){
				$data1[$i] = floor(10000*($data1[$i]-$expr)+0.5)/10000;
				my $absx = abs($data1[$i]); 
				if($logratio_col[$i] < $absx){
					$logratio_col[$i] = $absx;
				}
			}
		}
		push(@data,\@data1);
		push(@rows,$id);
		push(@symbols,$symbol);
		push(@average,$expr);
	}elsif($analysis_type eq "search"){  ######## USED FOR GENE SEARCH
		my @annot_list = @data1;
		splice(@annot_list,0,$nCol+8);
		splice(@data1,$nCol);
		if($similar_genes_id){
			if($FDR > 0.05){ next; }
			my ($r,$n,$bb,$aa) = pearson_correlation(\@search_pattern,\@data1);
			if(abs($r) < $correlation){ next; }
			my $scope = $pattern_scope*abs($bb);
			if($scope >= $log_fold_change_search){
				$r = floor(10000*$r+0.5)/10000;
				my $fold = floor(100*exp($scope*$log10)+0.5)/100;
				push(@data_list,[$id,$MSE,\@data1,\@annot_list,$expr,$F,$r,$fold]);
			}
			next;
		}
		if($category<0){   	#search for ALL CATEGORIES
			if($line !~ /$searchTerm/i){ next; }
			if(!$symbol_found && uc($symbol) eq uc($searchTerm)){ $symbol_found=1; @data_list=(); }
			if($symbol_found && uc($symbol) eq uc($searchTerm)){
				push(@data_list,[$id,$MSE,\@data1,\@annot_list,$expr,$F]);
			}else{
				push(@optional_list,[$id,$MSE,\@data1,\@annot_list,$expr,$F]);
			}
		}elsif($category==0){   #search for probe ID
			if($id =~ /^$searchTerm$/i){	#search for exact ID
				push(@data_list,[$id,$MSE,\@data1,\@annot_list,$expr,$F]);
				last;
			}
		}elsif($category>0 && $annot_search =~ /Entrez|GenBank|Refseq|Ensembl/i){
			if($annot_list[$category-1]=~ /^$searchTerm$/i){
				push(@data_list,[$id,$MSE,\@data1,\@annot_list,$expr,$F]);
			}
		}elsif($category==1){   #search for symbol
			my $term = $annot_list[$category-1];
			if($term=~ /$searchTerm/i){ push(@optional_list,[$id,$MSE,\@data1,\@annot_list,$expr,$F]); }
			#Test exact matching:
			if(!$symbol_found && uc($symbol) eq uc($searchTerm)){ $symbol_found=1; @data_list=(); }
			if($symbol_found && uc($symbol) eq uc($searchTerm)){
				push(@data_list,[$id,$MSE,\@data1,\@annot_list,$expr,$F]);
			}
		}elsif($annot_list[$category-1] =~ /$searchTerm/i){
			push(@data_list,[$id,$MSE,\@data1,\@annot_list,$expr,$F]);
		}
		if(@aliasTerm){
			for(my $j=0; $j<@aliasTerm; $j++){
				if($symbol =~ /^$aliasTerm[$j]$/i || $id =~ /^$aliasTerm[$j]$/i){
					$annot_list[0] .= "($search_term)";
					push(@alias_list,[$id,$MSE,\@data1,\@annot_list,$expr,$F]);
				}
			}
		}
	}else{
		error_message("Wrong analysis type: $analysis_type");
	}
}
close INFO;

if($analysis_type eq "pca_cluster"){
	print_cluster_table(\@cluster_list);
}
my $fileID = get_outputID(1);
my $nRow;
if($analysis_type =~ /^heatmap$|^pca$/){
	for(my $i=$nCol-1; $i>=0 && $filter_cols; $i--){
		if($logratio_col[$i] < $logratio_threshold){
			foreach my $ref (@data){
				splice(@$ref,$i,1);
			}
			splice(@headers,$i,1);
			$nCol--;
		}
	}
	$nRow = @rows;
	if(!$nRow){ error_message("Empty matrix after filtering!"); }

	# If 'use_replications' then extract data from matrix file
	if($hashInput{"use_replications"} eq "on"){
		my %hashMatrix;
		parse_matrix_file("$PATH_DATA/$file_matrix_full", \%hashMatrix);
		my $ref = $hashMatrix{"sample_title"};
		my @headers_raw = @$ref;
		my @headers1;
		my @replic_col;
		my %hash;
		my $ii=0;
		for(my $i=0; $i<@headers_raw; $i++){
			if(!$hash{$headers_raw[$i]}){
				$headers1[$ii] = $headers_raw[$i];
				$hash{$headers_raw[$i]} = $ii+1;
				$replic_col[$ii] = [$i];
				$ii++;
			}else{
				my $i1 = $hash{$headers_raw[$i]}-1;
				push(@{$replic_col[$i1]},$i);
			}
		}
		for(my $i=0; $i<@headers; $i++){
			while($i<@headers1 && $headers[$i] ne $headers1[$i]){
				splice(@headers1,$i,1);
				splice(@replic_col,$i,1);
				$i--;
			}
			if($i==@headers1){ error_message("Headers don't match: $headers[$i]"); }
		}
		@headers=();
		for(my $i=0; $i<@headers1; $i++){
			my $len = length($headers1[$i]);
			if($len>$maxHeaderLength-6){ $headers1[$i]=substr($headers1[$i],0,$maxHeaderLength-6); }
			for(my $j=0; $j<@{$replic_col[$i]}; $j++){
				my $irep = $j+1;
				push(@headers,"$headers1[$i] rep$irep");
			}
		}
		$nCol = @headers;
		open(INFO,"<$PATH_DATA/$file_matrix_full") or error_message("Cannot open matrix file");
		while(my $line = <INFO>){
			chop $line;
			if($line =~ /^!/ || !$line){ next; }
			last;
		}
		@data = ();
		my $log10 = log(10);
		my $iii = 0;
		my $xmin = 0.001*$hashMatrix{"series_xmin"};
		while(my $line = <INFO>){
			chop $line;
			if($line =~ /^!/ || $iii > @rows){ last; }
			my ($id,@data1) = split(/\t/, $line);
			$id =~ s/\"//g;
			if($id ne $rows[$iii]){ next; }
			my @data2;
			for(my $i=0; $i<@headers1; $i++){
				for(my $j=0; $j<@{$replic_col[$i]}; $j++){
					my $xx = $data1[$replic_col[$i]->[$j]];
					if($xx > $MISSING){
						if($xx < $xmin){ $xx=$xmin; }
						$xx = log($xx)/$log10;
					}
					push(@data2,$xx);
				}
			}
			my $median = median(\@data2);
			for(my $i=0; $i<@data2; $i++){
				if($data2[$i] > $MISSING){
					$data2[$i] = floor(10000*($data2[$i]-$median)+0.5)/10000;
				}else{
					$data2[$i] = 0;
				}
			}
			$data[$iii++] = \@data2;
		}
		close INFO;
	}
	for(my $i=0; $i<$nRow; $i++){
		my $id = $rows[$i];
		my $symbol = $symbols[$i];
		my $len = length($symbol);
		if($symbol && $id ne $symbol && $len<25){
			$rows[$i] .= " ($symbol)";
		}
	}
	for(my $i=0; $i<$nCol; ++$i){
		my $len = length($headers[$i]);
		if($len>$maxHeaderLength){ $headers[$i]=substr($headers[$i],0,$maxHeaderLength); }
	}
	open(OUT, ">$PATH_OUTPUT/$fileID.txt") or error_message("In matrix_explore1 cannot open file $fileID.txt"); 
	print OUT "Genes";
	for(my $i=0; $i<$nCol; ++$i){
		print OUT "\t$headers[$i]";
	}
	print OUT "\n";
	for(my $irow=0; $irow<$nRow; ++$irow){
		print OUT $rows[$irow];
		my $ref = $data[$irow];
		for(my $i=0; $i<$nCol; ++$i){
			print OUT "\t$ref->[$i]";
		}
		print OUT "\n";
	}
	close OUT;
}
#my $fileID1 = get_outputID(1);
#`cp $PATH_OUTPUT/$fileID.txt $PATH_OUTPUT/$fileID1.txt`;
#print "Original <a href=../output/$fileID1.txt>file</a><p>\n";
if($analysis_type =~ /^heatmap$|^pca$/){
	if($nCol*$nRow > 30000000){
		my $text = "<h3>The data set is too large</h3>\n";
		$text .= "N columns = $nCol; N rows = $nRow<br>\n";
		$text .= "Link to the data file (tab-delimited text): <a href=$HOME_ADDRESS/output/$fileID.txt>DATA</a><br>\n";
		$text .= "Use filtering to reduce the size of matrix or use alternative plotting software<br>";
		terminal_window($text);
	}elsif($nRow>0 && $nRow <= $nCol && $analysis_type =~ /^pca$/){
		terminal_window("<H3>Not enough genes match selection criteia</H3>N genes = $nRow should be greater than the number of columns.<br>Consider using less stringent criteria. For example, Fold change = 1 and FDR = 0.5");
	}elsif($nCol*$nRow==0){
		terminal_window("<H3>No genes match selection criteia</H3>Consider using less stringent criteria. For example, set Fold change threshold = 1 and FDR threshold = 0.5");
	}
}
my $evaluate = $nCol*$nRow/500000;
if($evaluate < $nCol/400){ $evaluate < $nCol/400; }
if($analysis_type =~ /^heatmap$/){
	$hashInput{"runID"}=$RUN_PCA_CLUSTER;
	my $logFileID = get_outputID(3);
	$hashInput{"logFileID"}=$logFileID;
	my $data = [$fileID,$nCol,$nRow,$matrix_cluster,$diagonal,0,"Gene expression heatmap"];
	interrupt_program($evaluate,$data);
	exit(0);
}
if($analysis_type =~ /^pca$/){
	my $npc = int(1.7*sqrt($nCol));
	if($npc > 7){ $npc=8; }
	$hashInput{"runID"}=$RUN_PCA_PCA;
	$hashInput{"fileID"}=$fileID;
	my $logFileID = get_outputID(3);
	$hashInput{"logFileID"}=$logFileID;
	my $data = [$fileID,$nCol,$nRow,$matrix_cluster,\@average];
	interrupt_program($evaluate,$data);
	exit(0);
}
########### THIS SECTION IS USED FOR GENE SEARCH ONLY!
if(!@data_list && @optional_list){
	@data_list = @optional_list;
}
if(@alias_list){
	push(@data_list,@alias_list);
}
if($search_term =~ /^similar_genes_/ && (!@data_list || @data_list==1 && $data_list[0]->[0] eq $similar_genes_id)){
	terminal_window("<H2>No genes found with similar expression profile to $similar_genes_id</H2>");
}
my @table2;
my $plotID = -1;
my $sort_histogram = $hashInput{"sort_histogram"};
if(@data_list == 1){
	my($id,$MSE,$row_ref) = @{$data_list[0]};
	my @SE=();
	for(my $i=0; $i<$nCol; ++$i){
		$SE[$i] = sqrt($MSE/$n_repl[$i]);
	}
	$plotID = get_outputID(2);
	@table2 = plot_histogram_horiz("deviation",$nCol,1,$row_ref,\@headers,$plotID,"Expression (Log10)",0,\@SE,\@n_repl);
}
print "<HTML><HEAD><TITLE>ExAtlas: Search Results</TITLE>\n";
my $parameter="";
if(@data_list>=5){ $parameter = "update_description();"; }
print_header($parameter);
if($plotID > 0){	#Search results for a single gene
	my($id,$MSE,$row_ref,$annot_ref,$expr,$F) = @{$data_list[0]};
	print "<SCRIPT language=JavaScript>\n";
	print "<!--\n";
	print "function find_similar_genes() {\n";
	print "	document.form_search.search_term.value=\"similar_genes_$id\";\n";
	print "	var x = Math.round(Math.random()*10000);\n";
	print "	document.form_search.target = \"_BLANK\"+x;\n";
	print "	document.form_search.submit();\n";
	print "}\n";
	print "function show_replications() {\n";
	print "	document.form_search.search_term.value=\"replications_$id\";\n";
	print "	var x = Math.round(Math.random()*10000);\n";
	print "	document.form_search.target = \"_BLANK\"+x;\n";
	print "	document.form_search.submit();\n";
	print "}\n";
	print "<!-- end script --></SCRIPT>\n";
	print "<FORM NAME=form_search ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
	print "<font size=+3><b>Search Results for $annot_search \"$search_term\"</b></font><p>\n";
	print "<b>ProbeID</b>: $id &nbsp; <b>Average log-expression</b>: $expr &nbsp; <b>F-statistics</b>: $F<br>\n";
	for(my $i=0; $i<@annotations; ++$i){
		print "<b>$annotations[$i]</b>: $annot_ref->[$i] &nbsp; ";
	}
	if($legendID){
		my $x = int(10000*rand());
		print "<br><b>Legend:</b> <a href=$HOME_ADDRESS/output/$legendID.txt target=_blank$x>Legend file</a>\n";
	}
	print "<br><INPUT TYPE=button value=\"   Show replications   \" onclick=show_replications();><p>\n";
	my @correlation_list = (0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9);
	my @fold_list = (1,1.5,2,3,4,5,7);
	print "<b>Find genes with same profile:</b> Correlation = <select name=correlation>\n";
	for(my $i=0; $i<@correlation_list; ++$i){ 
		print "<option value=$correlation_list[$i]"; if($correlation_list[$i]==0.7){ print " selected"; } print "> $correlation_list[$i]\n";
	}
	print "</select>\n";
	print " &nbsp; Fold change = <select name=fold_change_search>\n";
	for(my $i=0; $i<@fold_list; ++$i){ 
		print "<option value=$fold_list[$i]"; if($fold_list[$i]==2){ print " selected"; } print "> $fold_list[$i]\n";
	}
	print "</select>\n";
	print "<INPUT TYPE=button value=\" Find genes \" onclick=find_similar_genes();><p>\n";
		print "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
	print "<INPUT NAME=\"action\" TYPE=hidden VALUE=\"matrix_explore1\">\n";
	print "<INPUT NAME=\"analysis\" TYPE=hidden VALUE=search>\n";
	print "<INPUT NAME=\"file_matrix\" TYPE=hidden VALUE=\"$file_matrix\">\n";
	print "<INPUT NAME=\"organismID\" TYPE=hidden VALUE=$organismID>\n";
	print "<INPUT NAME=search_term TYPE=hidden VALUE=\"similar_genes_$id\">\n";
	if($legendID){ print "<INPUT NAME=legendID TYPE=hidden VALUE=$legendID>\n"; }
	print "</FORM>\n";
	print "<IMG SRC=../output/$plotID.gif BORDER=0></p>\n";
	print "<TABLE>\n";
	print "<TR><TD><b>Sample</b><TD><b>Expression</b><TD><b>N repl.</b>\n";
	for(my $i=0; $i<@table2; ++$i){
		print $table2[$i]."\n";
	}
	print "</TABLE>\n";
}
elsif(!@data_list){
	print "<H2>Search term \"$search_term\" not found</H2>\n";
	print "Search for gene symbols is based on exact match. To find partial match try to\n";
	print "select 'All' as a search type.\n";
	print "Also you can try other search targets (GenBank acc#, etc)<p>\n";
}
else{
	my $fileID = get_outputID(1);
	my @geneset_list;
	if(@data_list>=5){
		@geneset_list = get_geneset_list();
		filter_list_by_organism(\@geneset_list, $organismID);
		my ($items,$descriptions) = get_array_lists(\@geneset_list);
		print "<SCRIPT language=JavaScript>\n";
		print "<!--\n";
		print "geneset_list = new Array($items);\n";
		print "geneset_description = new Array($descriptions);\n";
		print "function geneset_overlap(file) {\n";
		print "	document.form_genelist.upload_geneset.value=file;\n";
		print "	document.form_genelist.action.value=\"geneset_overlap\";\n";
		print "	var x = Math.round(Math.random()*10000);\n";
		print "	document.form_genelist.target = \"_BLANK\"+x;\n";
		print "	document.form_genelist.submit();\n";
		print "}\n";
		print "function update_description() {\n";
		print "	var index;\n";
		print "	index = document.form_genelist.file_geneset1.selectedIndex;\n";
		print "	document.form_genelist.description_geneset1.value = geneset_description[index];\n";
		print "}\n";
		print "<!-- end script --></SCRIPT>\n";
	}
	open(OUT,">$PATH_OUTPUT/$fileID.txt");
	my @sorted;
	my @negative;
	my ($Npos,$Nneg);
	if($similar_genes_id){
		@sorted = sort {$b->[6]<=>$a->[6]} @data_list;
		my $ii=0;
		while($ii<@sorted && $sorted[$ii]->[6] > 0){ ++$ii; }
		if($ii<@sorted){
			@negative = splice(@sorted,$ii);
			@negative = sort {$a->[6]<=>$b->[6]} @negative;
		}
		$Npos = @sorted;
		$Nneg = @negative ;
		print "<font size=+3><b>Genes with pattern similar to \"$similar_genes_id\"</b></font><p>\n";
		print "<b>Number of genes:</b> Positive correlation = $Npos &nbsp; &nbsp; Negative correlation = $Nneg<br>\n";
		print "<b>Correlation:</b> $correlation<br>\n";
		print "<b>Fold change:</b> $fold_change_search<br>\n";
		print OUT "!Genelist_description\tGenes positively correlated with $similar_genes_id\n";
	}else{
		@sorted = sort {$b->[5]<=>$a->[5]} @data_list;
		print "<font size=+3><b>Search Results for \"$search_term\"</b></font><p>\n";
		print OUT "!Genelist_description\tSearch results for \'$search_term\'\n";
	}
	my $x = int(10000*rand());
	print "<b>Gene list as tab-delimited text:</b> <a href=$HOME_ADDRESS/output/$fileID.txt target=_BLANK$x>gene list</a>\n";
	if($legendID){
		my $x = int(10000*rand());
		print "<br><b>Legend:</b> <a href=$HOME_ADDRESS/output/$legendID.txt target=_blank$x>Legend file</a>\n";
	}
	print "<p><FORM NAME=form_genelist ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
	print "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
	print "<INPUT NAME=action TYPE=hidden>\n";
	print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
	print "<INPUT NAME=upload_geneset TYPE=hidden>\n";
	if($legendID){ print "<INPUT NAME=legendID TYPE=hidden VALUE=$legendID>\n"; }
	if(@data_list>=5){
		print "<TABLE BORDER=0>\n";
		my $menu_text = menu_geneset_overlap(\@geneset_list,$fileID);
		print $menu_text;
		print "</TABLE><p>\n";
	}
	print "</FORM><TABLE>\n";
	print "<TR><TD><b>ProbeID<TD WIDTH=120><b>Aver expression<TD><b>F-statistics";
	print OUT "ProbeID\tAver expression\tF-statistics";
	if($similar_genes_id){
		print "<TD><b>Correlation<TD WIDTH=100><b>Fold change";
		print OUT "\tCorrelation\tFold change";
	}
	for(my $i=0; $i<2; ++$i){
		print "<TD";
		if($i==1){ print " WIDTH=500"; }
		print "><b>$annotations[$i]";
		print OUT "\t$annotations[$i]";
	}
	print "\n";
	print OUT "\n";
	for(my $i=0; $i<@sorted; ++$i){
		my($id,$MSE,$row_ref,$annot_ref,$expr,$F,$r,$fold) = @{$sorted[$i]};
		my $x = int(10000*rand());
		print "<TR><TD><a href=\"exatlas.cgi?category=0&search_term=$id&action=matrix_explore1&analysis=search&file_matrix=$file_matrix&sessionID=$sessionID&organismID=$organismID&sort_histogram=$sort_histogram\" target=_BLANK$x>$id<a>";
		print "<TD><center>$expr<TD><center>$F";
		print OUT "$id\t$expr\t$F";
		if($similar_genes_id){
			print "<TD><center>$r<TD><center>$fold";
			print OUT "\t$r\t$fold";
		}
		for(my $i=0; $i<2; ++$i){
			print "<TD>$annot_ref->[$i]"; 
			print OUT "\t$annot_ref->[$i]";
		}
		print "\n";
		print OUT "\n";
	}
	print "</TABLE><p>\n";
	if(@negative){
		print "<font size=+2><b>Negatively correlated genes</b></font><p>\n";
		print OUT "!Genelist_description\tGenes negatively correlated with $similar_genes_id\n";
		print "<TABLE>\n";
		print "<TR><TD><b>ProbeID<TD WIDTH=120><b>Aver. expression<TD><b>F-statistics<TD><b>Correlation<TD WIDTH=100><b>Fold change";
		print OUT "ProbeID\tAver expression\tF-statistics\tCorrelation\tFold change";
		for(my $i=0; $i<2; ++$i){
			print "<TD";
			if($i==1){ print " WIDTH=500"; }
			print "><b>$annotations[$i]";
			print OUT "\t$annotations[$i]";
		}
		print "\n";
		print OUT "\n";
		for(my $i=0; $i<@negative; ++$i){
			my($id,$MSE,$row_ref,$annot_ref,$expr,$F,$r,$fold) = @{$negative[$i]};
			my $x = int(10000*rand());
			print "<TR><TD><a href=\"exatlas.cgi?category=0&search_term=$id&action=matrix_explore1&analysis=search&file_matrix=$file_matrix&sessionID=$sessionID&organismID=$organismID\" target=_BLANK$x>$id<a>";
			print "<TD><center>$expr<TD><center>$F<TD><center>$r<TD><center>$fold";
			print OUT "$id\t$expr\t$F\t$r\t$fold";
			for(my $i=0; $i<2; ++$i){
				print "<TD>$annot_ref->[$i]"; 
				print OUT "\t$annot_ref->[$i]";
			}
			print "\n";
			print OUT "\n";
		}
		print "</TABLE><p>\n";
	}
}
print "<p><i>Please report any problems to <a href=mailto:sharoval\@mail.nih.gov>webmaster</a><p>\n";
print "<INPUT TYPE=button VALUE=\" Close window \" LANGUAGE=\"javascript\" onClick=\"window.close();\">\n";
print "</BODY></HTML>\n";
exit(0);
}

#**************************************
sub  read_anova_output
#**************************************
{
my $headers_ref = shift;
my $nReplications_ref = shift;
my $columns_ref = shift;
my $rowHeaders_ref = shift;
my $symbols_ref = shift;
my $geneNames_ref = shift;
my $MSE_ref = shift;
my $FDR_ref = shift;
my $median_ref = shift;
my $useSymbols = shift;
my $subtractMedian = shift;

my $file_matrix = $hashInput{"file_matrix"};
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_anova = "$loginname-anova-$file_matrix";
}
open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file");
my $line = <INFO>;
chop $line;
my ($junk,$junk1,@headers) =split(/\t/,$line);
if(!@headers){ error_message("Headers not found in anova file"); }
my $nCol=0;
while($nCol<@headers && $headers[$nCol] !~ /^Var\(/){ ++$nCol; }
my $nRep=0;
for(my $i=0; $i<$nCol; ++$i){
	$nReplications_ref->[$i] = 1;
	if($headers[$i]=~/ \(\d+\)$/ && $headers[$i] !~ /^Mean/){
		my @items = split(/ \(/,$headers[$i]);
		my $n = pop(@items);
		my $name = join(" (",@items);
		$n =~ s/\)$//;
		if($n < 1){ $n=1; }
		$nRep += $n;
		$nReplications_ref->[$i] = $n;
		$headers[$i] = $name;
	}else{
		$headers[$i] =~ s/^Mean\(//;
		$headers[$i] =~ s/\)$//;
		$nRep++;
	}
}
@$headers_ref = @headers;
my %symbols;
my $count=0;
while(my $line = <INFO>){
	chop $line;
	my ($probe_id,$expr,@data1) =split(/\t/, $line);
	my $FDR = $data1[$nCol+6];
	my $F = $data1[$nCol+4];
	my $MSE = $data1[$nCol+3];
	my $symbol = $data1[$nCol+8];
	my $gene_name = $data1[$nCol+9];
	#if($nRep<=$nCol){ $F=$expr; }
	splice(@data1,$nCol);
	my $median = median(\@data1);
	my $row = $count;
	if($useSymbols){
		if(!$symbol){ next; }
		my $ref = $symbols{$symbol};
		if(!$ref){
			$symbols{$symbol} = [$probe_id,$F,$count];
			$row=$count++;
		}elsif($F > $ref->[1]){
			$ref->[0] = $probe_id;
			$ref->[1] = $F;
			$row = $ref->[2];
		}
	}else{
		$row=$count++;
	}
	for(my $j=0; $j<$nCol; $j++){
		my $x = $data1[$j];
		if($subtractMedian){
			if($x>$MISSING){ $x -= $median; }
		}
		$columns_ref->[$j]->[$row] = $x;
	}
	$rowHeaders_ref->[$row] = $probe_id;
	$symbols_ref->[$row] = $symbol;
	$geneNames_ref->[$row] = $gene_name;
	$MSE_ref->[$row] = $MSE;
	$FDR_ref->[$row] = $FDR;
	$median_ref->[$row] = $median;
}
close INFO;
return ($nCol,$count);
}

#**************************************
sub  pairwise_compare
#**************************************
{
my @headers=();
my @n_repl=();
my @columns=();
my @probeID=();
my @symbols=();
my @geneNames=();
my @MSE=();
my @FDR=();
my @median=();
my $useSymbols=0;
if($hashInput{"useSymbols"} eq "on"){ $useSymbols=1; }
my $column1 = $hashInput{"select_column"};
my $column2 = $hashInput{"compare_column"};
my ($nCol,$nRow) = read_anova_output(\@headers,\@n_repl,\@columns,\@probeID,\@symbols,\@geneNames,\@MSE,\@FDR,\@median,$useSymbols,0);

if(!$nRow){ error_message($nCol); }
my @sorted_index;
my @DIST;
my $icol = $column1-1;
for(my $i=0; $i<$nCol; $i++){
	my $dist=0;
	for(my $k=0; $k<$nRow; $k++){
		my $d = $columns[$icol]->[$k]-$columns[$i]->[$k];
		$dist += $d*$d;
	}
	$dist = sqrt($dist);
	$DIST[$i] = $dist;
}
@sorted_index = sort {$DIST[$b]<=>$DIST[$a]} 0..($nCol-1);
my $maxdist = $DIST[$sorted_index[0]];
my $i = int($nCol/2);
while($i<$nCol-1){
	if( $i>=3 && $DIST[$sorted_index[$i]] < $maxdist/3){ last; }
	$i++;
}
splice(@sorted_index,$i); #Remove correlated columns, used for z_spec
my @logratio_col;
my @points_pairwise;
my $coefficient = 1.0/$n_repl[$column1-1];
my ($minx,$miny) = (1000000,1000000);
my ($maxx,$maxy) = (-1000000,-1000000);
if($column2>0){ $coefficient += 1.0/$n_repl[$column2-1]; }
else{ $coefficient += 1.0/$nCol; }
my $FDRthresh = $hashInput{"filter_FDR"};
my $exprThreshold = $hashInput{"exprThreshold_pairwise"};
if($exprThreshold>0){ $exprThreshold=log($exprThreshold)/log(10); }
else{ $exprThreshold = -1000; }
for(my $irow=0; $irow<$nRow; $irow++){
	my $y = $columns[$icol]->[$irow];
	my $x;
	if($column2==0){ $x = $median[$irow]; }
	else{ $x = $columns[$column2-1]->[$irow]; }
	if($x==$MISSING || $y==$MISSING){ next; }
	$x = floor(10000*$x+0.5)/10000;
	$y = floor(10000*$y+0.5)/10000;
	if($minx > $x){ $minx=$x; }
	if($maxx < $x){ $maxx=$x; }
	if($miny > $y){ $miny=$y; }
	if($maxy < $y){ $maxy=$y; }
	if($MSE[$irow]<=1.0E-10){ $MSE[$irow]=1.0E-10; }
	my $z = abs($y-$x)/sqrt($MSE[$irow]*$coefficient);
	$z = floor(10000*$z+0.5)/10000;
	if($x<$exprThreshold && $y<$exprThreshold){ $z=0; }
	my $sqrtMSE = sqrt($MSE[$irow]);
	my $z_spec = 0;
	my ($m,$sd,$n) = (0,0,0);
	foreach my $j (@sorted_index){
		my $x1 = $columns[$j]->[$irow];
		if($x1==$MISSING){ next; }
		$m += $x1;
		$sd += $x1*$x1;
		$n++;
	}
	if($n>1 && ($y-$x)*($y*$n-$m) > 0){
		$sd -= $m*$m/$n;
		if($sd<0.000001){ $sd=0.000001; }
		$sd = sqrt($sd/($n-1));
		if($sd<$sqrtMSE){ $sd=$sqrtMSE; }
		$m /= $n;
		if($sd>0){ $z_spec = floor(1000*abs($y-$m)/$sd+0.5)/1000; }
	}
	push(@points_pairwise,[$z,$x,$y,1,1,1,$symbols[$irow],$geneNames[$irow],$z_spec,$probeID[$irow]]);
}
my $header2 = "median";
if($column2>0){ $header2 = $headers[$column2-1]; }
make_scatterplot(\@points_pairwise,$headers[$column1-1],$header2,[$minx,$miny],[$maxx,$maxy]);
exit(0);
}

#**************************************
sub  pairwise_replications
#**************************************
{
my $file_matrix = $hashInput{"file_matrix"};
my $column1 = $hashInput{"select_column"};
my $column2 = $hashInput{"compare_column"};
my $replication1 = $hashInput{"select_replication1"};
my $replication2 = $hashInput{"select_replication2"};
my $organismID = $hashInput{"organismID"};
my $file_matrix_full = $file_matrix;
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_matrix_full = "$loginname-$file_matrix";
	$file_anova = "$loginname-anova-$file_matrix";
}
my %hashMatrix;
parse_matrix_file("$PATH_DATA/$file_matrix_full", \%hashMatrix);
my $organismID1 = $hashMatrix{"series_sample_taxid"};
if($organismID1 && $organismID != $organismID1){ error_message("Organism ID do not match"); }
my $ref = $hashMatrix{"sample_title"};
my @headers_raw = @$ref;
my $nCol = @headers_raw;
my @headers;
my @replic_col;
my %hash;
my $ii=0;
for(my $i=0; $i<$nCol; $i++){
	if(!$hash{$headers_raw[$i]}){
		$headers[$ii] = $headers_raw[$i];
		$hash{$headers_raw[$i]} = $ii+1;
		$replic_col[$ii] = [$i];
		$ii++;
	}else{
		my $i1 = $hash{$headers_raw[$i]}-1;
		push(@{$replic_col[$i1]},$i);
	}
}
open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file");
my $line = <INFO>;
chop $line;
my ($junk,$junk1,@headers_anova) =split(/\t/,$line);
if(!@headers_anova){ error_message("Headers not found in file_anova"); }
my $nCol1=0;
while($nCol1<@headers_anova && $headers_anova[$nCol1] !~ /^Var\(/){ ++$nCol1; }
my @annotations;
while(my $line = <INFO>){
	chop $line;
	my ($id,$expr,@data1) =split(/\t/, $line);
	my $MSE = $data1[$nCol1+3];
	my $symbol = $data1[$nCol1+8];
	my $gene_name = $data1[$nCol1+9];
	push(@annotations,[$id,$MSE,$symbol,$gene_name]);
}
close INFO;
open(INFO,"<$PATH_DATA/$file_matrix_full") or error_message("Cannot open matrix file");
while(my $line = <INFO>){
	chop $line;
	if($line =~ /^!/ || !$line){ next; }
	last;
}
my @points_pairwise;
my $log10 = log(10);
my ($minx,$miny) = (1000000,1000000);
my ($maxx,$maxy) = (-1000000,-1000000);
my $xmin = 0.001*$hashMatrix{"series_xmin"};
while(my $line = <INFO>){
	chop $line;
	if($line =~ /^!/){ last; }
	my ($id,@data) =split(/\t/, $line);
	$id =~ s/\"//g;
	my $ref = shift(@annotations);
	if($ref->[0] ne $id){ error_message("Probe IDs don't match: $ref->[0] ne $id"); }
	my ($junk,$MSE,$symbol,$gene_name) = @$ref;
	my ($x,$y);
	my @data1;
	my @nn1;
	my $coefficient=0;
	for(my $i=0; $i<@headers; $i++){
		my ($nn,$ss)=(0,0);
		for(my $j=0; $j<@{$replic_col[$i]}; $j++){
			my $xx = $data[$replic_col[$i]->[$j]];
			if($xx==$MISSING){ next; }
			if($xx < $xmin){ $xx=$xmin; }
			$xx = log($xx)/$log10;
			$ss += $xx;
			$nn++;
		}
		if($nn){ $data1[$i] = $ss/$nn; }
		else{ $data1[$i] = $MISSING; }
		$nn1[$i] = $nn;
	}
	if($replication1==0){
		$y = $data1[$column1-1];
		if($y==$MISSING){ next; }
		$coefficient = 1/$nn1[$column1-1];
	}else{
		my $xx = $data[$replic_col[$column1-1]->[$replication1-1]];
		if($xx==$MISSING){ next; }
		if($xx < $xmin){ $xx=$xmin*(0.5+0.5*rand()); }
		$y = log($xx)/$log10;
		$coefficient = 1;
	}
	if($column2==0){
		$x = median(\@data1);
		$coefficient += 1/@headers_raw;
	}elsif($replication2==0){
		$x = $data1[$column2-1];
		if($x==$MISSING){ next; }
		$coefficient += 1/$nn1[$column2-1];
	}else{
		my $xx = $data[$replic_col[$column2-1]->[$replication2-1]];
		if($xx==$MISSING){ next; }
		if($xx < $xmin){ $xx=$xmin*(0.5+0.5*rand()); }
		$x = log($xx)/$log10;
		$coefficient += 1;
	}
	if($x==$MISSING){ next; }
	$x = floor(10000*$x+0.5)/10000;
	$y = floor(10000*$y+0.5)/10000;
	if($minx > $x){ $minx=$x; }
	if($maxx < $x){ $maxx=$x; }
	if($miny > $y){ $miny=$y; }
	if($maxy < $y){ $maxy=$y; }
	my $z = 0;
	if($MSE*$coefficient > 0){
		$z = abs($y-$x)/sqrt($MSE*$coefficient);
		$z = int(1000*$z+0.5)/1000;
	}
	push(@points_pairwise,[$z,$x,$y,1,1,1,$symbol,$gene_name,0,$id]);
}
close INFO;
my $header1 = $headers[$column1-1];
if($replication1 ne "All"){ $header1 .= "_rep$replication1"; }
my $header2 = "median";
if($column2>0){
	$header2 = $headers[$column2-1];
	if($replication2 ne "All"){ $header2 .= "_rep$replication2"; }
}
make_scatterplot(\@points_pairwise,$header1,$header2,[$minx,$miny],[$maxx,$maxy]);
exit(0);
}

#**************************************
sub  show_replications
#**************************************
{
my $file_matrix = $hashInput{"file_matrix"};
my $organismID = $hashInput{"organismID"};
my $file_matrix_full = $file_matrix;
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_matrix_full = "$loginname-$file_matrix";
	$file_anova = "$loginname-anova-$file_matrix";
}
my $replications_id = $hashInput{"search_term"};
$replications_id =~ s/^replications_//;
my $legendID = $hashInput{"legendID"};

my %hashMatrix;
parse_matrix_file("$PATH_DATA/$file_matrix_full", \%hashMatrix);
my $organismID1 = $hashMatrix{"series_sample_taxid"};
if($organismID1 && $organismID != $organismID1){ error_message("Organism ID do not match"); }
my $ref = $hashMatrix{"sample_title"};
my @headers_raw = @$ref;
my $nCol = @headers_raw;
my @headers;
my @replic_col;
my %hash;
my $ii=0;
for(my $i=0; $i<$nCol; $i++){
	if(!$hash{$headers_raw[$i]}){
		$headers[$ii] = $headers_raw[$i];
		$hash{$headers_raw[$i]} = $ii+1;
		$replic_col[$ii] = [$i];
		$ii++;
	}else{
		my $i1 = $hash{$headers_raw[$i]}-1;
		push(@{$replic_col[$i1]},$i);
	}
}
open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file");
my $line = <INFO>;
chop $line;
my ($junk,$junk1,@headers_anova) =split(/\t/,$line);
if(!@headers_anova){ error_message("Headers not found in file_anova"); }
my $nCol1=0;
while($nCol1<@headers_anova && $headers_anova[$nCol1] !~ /^Var\(/){ ++$nCol1; }
my @annotation;
while(my $line = <INFO>){
	chop $line;
	my ($id,$expr,@data1) =split(/\t/, $line);
	if($id eq $replications_id){
		push(@annotation,$data1[$nCol1+8],$data1[$nCol1+9]);
		last;
	}
}
close INFO;
open(INFO,"<$PATH_DATA/$file_matrix_full") or error_message("Cannot open matrix file");
while(my $line = <INFO>){
	chop $line;
	if($line =~ /^!/ || !$line){ next; }
	last;
}
my @data1;
my $log10 = log(10);
while(my $line = <INFO>){
	chop $line;
	if($line =~ /^!/){ last; }
	my ($id,@data) =split(/\t/, $line);
	$id =~ s/\"//g;
	if($id eq $replications_id){
		@data1 = @data;
		last;
	}
}
close INFO;
my $median = median(\@data1);
my @data2;
my @headers2;
my $xmin = 0.001*$hashMatrix{"series_xmin"};
for(my $i=0; $i<@headers; $i++){
	for(my $j=0; $j<@{$replic_col[$i]}; $j++){
		my $xx = $data1[$replic_col[$i]->[$j]];
		if($xx==$MISSING){ $xx=$median; }
		if($xx < $xmin){ $xx=$xmin; }
		push(@data2,log($xx)/$log10);
		my $rep = $j+1;
		push(@headers2,"$headers[$i] rep$rep");
	}
}
close INFO;
my $plotID = get_outputID(2);
my $n = @headers2;
my @table2 = plot_histogram_horiz("deviation",$n,1,\@data2,\@headers2,$plotID,"Expression (Log10)",0);

print "<HTML><HEAD><TITLE>ExAtlas: Replications Plot</TITLE>\n";
print_header();
my($symbol,$gene_name) = @annotation;
print "<font size=+3><b>Replication Plot for \"$replications_id\"</b></font><p>\n";
print "<b>Symbol</b>: $symbol &nbsp; <b>Gene name</b>: $gene_name\n";
if($legendID){
	my $x = int(10000*rand());
	print "<br><b>Legend:</b> <a href=$HOME_ADDRESS/output/$legendID.txt target=_blank$x>Legend file</a>\n";
}
print "<p><IMG SRC=../output/$plotID.gif BORDER=0></p>\n";
print "<TABLE>\n";
print "<TR><TD><b>Sample</b><TD><b>Expression</b>\n";
for(my $i=0; $i<@table2; ++$i){
	print $table2[$i]."\n";
}
print "</TABLE>\n";
exit(0);
}

#*************************************
sub   print_cluster_table
#*************************************
{
my $cluster_list = shift;

my $file_matrix = $hashInput{"file_matrix"};
my $description_matrix = $hashInput{"description_matrix"};
my $clusterPC = $hashInput{"clusterPC"};
my $clusterDir = $hashInput{"clusterDir"};
my $organismID = $hashInput{"organismID"};
my $legendID = $hashInput{"legendID"};
if($clusterPC eq "all"){  #Save PC clustes to a geneset
	my @geneset;
	foreach my $ref (@$cluster_list){
		my ($id,$logChange,$r,$PC,$dir,$symbol,$title) = @$ref;
		if(!$symbol){ next; }
		$PC =~ s/^PC//;
		my $ii = $PC*2;
		if($dir eq "negative"){ $ii++; }
		$geneset[$ii]->{$symbol}=1;
	}		
	my $file_geneset = $hashInput{"file_geneset"};
	my $file_geneset_full = "$loginname-$file_geneset";

	my $description_geneset = $hashInput{"description_geneset"};
	open (OUT, ">$PATH_DATA/$file_geneset_full") or error_message("Cannot write to geneset");
	for(my $i=2; $i<@geneset; $i+=2){
		for(my $j=0; $j<2; $j++){
			my $ii=$i+$j;
			my $refHash = $geneset[$ii];
			if(!$refHash){ next; }
			my @symbols = sort keys %$refHash;
			my $iPC = $i/2;
			my $setName = "PC$iPC$dirNames[$j]";
			print OUT "$setName\t$setName\t".join("\t",@symbols)."\n";
		}
	}
	close OUT;
	open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!");
	open(OUT, ">$PATH_INFO/$loginname-config1.txt") or error_message("Cannot update configuration file!");
	my $nLines;
	while(my $line = <INFO>){
		$nLines++;
		chop $line;
		if(length($line)<3){ next; }
		my @items = split(/[=\t]/,$line);
		if($items[0] ne "type_geneset" || $items[1] ne $file_geneset){
			print OUT $line."\n";
		}
	}
	close INFO;
	print OUT "type_geneset=$file_geneset\torganismID=$organismID";
	if($description_geneset){ print OUT "\tdescription=$description_geneset"; }
	print OUT "\tdate=$date_record\n";
	close OUT;
	my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
	if($nLines1 && $nLines1 > $nLines*0.9){
		copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
	}else{
		error_message("Failed to update configuration file!");
	}
	unlink "$PATH_INFO/$loginname-config1.txt";
	geneset_explore();
	exit(0);
}
my $log10 = log(10);
my @geneset_list = get_geneset_list();
filter_list_by_organism(\@geneset_list, $organismID);
my ($items,$descriptions) = get_array_lists(\@geneset_list);
my $fileID = get_outputID(1);
my $form_text .= "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
$form_text .= "<INPUT NAME=action TYPE=hidden VALUE=matrix_explore1>\n";
$form_text .= "<INPUT NAME=analysis TYPE=hidden VALUE=search>\n";
$form_text .= "<INPUT NAME=category TYPE=hidden VALUE=0>\n";
$form_text .= "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
$form_text .= "<INPUT NAME=description_matrix TYPE=hidden VALUE=\"$description_matrix\">\n";
$form_text .= "<INPUT NAME=description_geneset1 TYPE=hidden VALUE=\"Clustrer PC$clusterPC-$clusterDir\">\n";
$form_text .= "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
$form_text .= "<INPUT NAME=search_term TYPE=hidden>\n";
$form_text .= "<INPUT NAME=upload_geneset TYPE=hidden>\n";
if($legendID){ $form_text .= "<INPUT NAME=legendID TYPE=hidden VALUE=$legendID>\n"; }
$form_text .= "</FORM><p>\n";
$form_text .= "<HR NOSHADE></HR>\n";
$form_text .= "<INPUT TYPE=button VALUE=\"    Close window   \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
$form_text .= "</BODY>\n";
$form_text .= "</HTML>\n";
my $geneset_text = "<TABLE BORDER=0>\n";
my $menu_text = menu_geneset_overlap(\@geneset_list,$fileID);
$geneset_text .= "$menu_text</TABLE><p>\n";
my $nGenes = @$cluster_list;
open (OUT1, ">$PATH_OUTPUT/$fileID.txt") or error_message("Cannot write to $fileID.txt");
print OUT1 "!Number of genes = $nGenes\n";
print OUT1 "ProbeID\tFold change\tCorrelation\tGene symbol\tGene name\n";
print "<HTML><HEAD><TITLE>ExAtlas: Table of overexpressed genes</TITLE>\n";
if(@geneset_list && $nGenes>=5){
	print_header("update_description();");
}else{
	print_header();
}
print "<SCRIPT language=JavaScript>\n";
print "<!--\n";
print "geneset_list = new Array($items);\n";
print "geneset_description = new Array($descriptions);\n";
print "function plot_histogram(probeid) {\n";
print "	document.form_scatterplot.search_term.value=probeid;\n";
print "	document.form_scatterplot.action.value=\"matrix_explore1\";\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_scatterplot.target = \"_BLANK\"+x;\n";
print "	document.form_scatterplot.submit();\n";
print "}\n";
if(@geneset_list && $nGenes>=5){
	print "function geneset_overlap(file) {\n";
	print "	document.form_scatterplot.upload_geneset.value=file;\n";
	print "	document.form_scatterplot.action.value=\"geneset_overlap\";\n";
	print "	var x = Math.round(Math.random()*10000);\n";
	print "	document.form_scatterplot.target = \"_BLANK\"+x;\n";
	print "	document.form_scatterplot.submit();\n";
	print "}\n";
	print "function update_description() {\n";
	print "	var index;\n";
	print "	index = document.form_scatterplot.file_geneset1.selectedIndex;\n";
	print "	document.form_scatterplot.description_geneset.value = geneset_description[index];\n";
	print "}\n";
}
print "<!-- end script --></SCRIPT></HEAD>\n";
print "<H2>Cluster of genes for PC$clusterPC ($clusterDir direction)</H2>\n";
print "Number of genes = $nGenes<br>\n";
my $x = int(10000*rand());
print "<a href=$HOME_ADDRESS/output/$fileID.txt target=_BLANK$x>Table as tab-delimited text</a><p>\n";
print "<FORM NAME=form_scatterplot ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
if(@geneset_list && $nGenes>=5){
	print $geneset_text;
}
print "<TABLE BORDER=0><TR><TD><b>ProbeID<TD WIDTH=100><b>Fold change<TD><b>Correlation<TD><b>Gene symbol<TD WIDTH=500><b>Gene name\n";
foreach my $line (@$cluster_list){
	my ($id,$logChange,$r,$PC,$dir,$symbol,$title) = @$line;
	if($logChange>10){ $logChange=10; }
	elsif($logChange<-10){ $logChange=-10; }
	my $fold = int(1000*exp(abs($logChange)*$log10)+0.5)/1000;
	print "<TR><TD><a href=\"JavaScript:plot_histogram('$id');\">$id</a>\n";
	print "<TD><center>$fold<TD><center>$r<TD>$symbol<TD>$title\n";
	print OUT1 "$id\t$fold\t$r\t$symbol\t$title\n";
}
print "</TABLE>\n".$form_text;
close OUT1;
exit(0);
}

#*************************************
sub  make_scatterplot
#*************************************
{
my $points_pairwise = shift;
my $header1 = shift;
my $header2 = shift;
my $minx = shift;
my $maxx = shift;

my $organismID = $hashInput{"organismID"};
my $file_matrix = $hashInput{"file_matrix"};
my $description_matrix = $hashInput{"description_matrix"};
my $FDR_thresh = $hashInput{"FDR_pairwise"};
my $FDR_pvalue = $hashInput{"FDR_pvalue"};
my $foldChange_thresh = $hashInput{"foldChange_pairwise"};
my $legendID = $hashInput{"legendID"};
if(!$FDR_pvalue){ $FDR_pvalue="FDR"; }

my $FDR1 = 1;
my $log10 = log(10);
my $totalProbes = @$points_pairwise;
my $rank = $totalProbes;
my $thresh = log($foldChange_thresh)/$log10;
my @table_headers = ("!Matrix_name=$file_matrix");
if($description_matrix){ push(@table_headers,"!Matrix_description=$description_matrix"); }
push(@table_headers,"!Column_name=$header1","!Comparing_with=$header2");
push(@table_headers,"!Fold_change_threshold=$foldChange_thresh");
if($FDR_pvalue eq "FDR"){ push(@table_headers,"!FDR_threshold=$FDR_thresh"); }
else{ push(@table_headers,"!pvalue_threshold=$FDR_thresh"); }
if($legendID){ push(@table_headers,"!legendID=$legendID"); }
my @genes_under;
my @genes_over;
foreach my $ref (sort {$a->[0]<=>$b->[0]} @$points_pairwise){
	my $z = $ref->[0];
	my $p =  2*(1 - normal_distribution($z));
	$ref->[3] = $rank;
	my $FDR = int(10000*$p*$totalProbes/($rank--)+0.5)/10000;
	if($FDR > $FDR1){ $FDR = $FDR1; }
	$FDR1 = $FDR;
	my $logratio = $ref->[2]-$ref->[1];
	if($logratio>10){ $logratio=10; }
	if($logratio<-10){ $logratio=-10; }
	my $fold = int(1000*exp(abs($logratio)*$log10)+0.5)/1000;
	$p = int(10000*$p+0.5)/10000;
	my $lineOut = "$ref->[9]\t$ref->[2]\t$ref->[1]\t$fold\t$ref->[0]\t$p\t$FDR";
	if($header2 eq "median"){ $lineOut .= "\t$ref->[8]"; }
	$lineOut .= "\t$ref->[6]\t$ref->[7]";
	#if($z>0 && abs($logratio) >= $thresh && $FDR <= $FDR_thresh){
	if($z>0 && abs($logratio) >= $thresh && ($FDR_pvalue eq "FDR" && $FDR <= $FDR_thresh || $FDR_pvalue ne "FDR" && $p <= $FDR_thresh)){
		if($logratio>0){
			push(@genes_over,$lineOut);
		}
		else{
			push(@genes_under,$lineOut);
		}
	}
	$ref->[4] = $FDR;
	if($p>0.000001){
		$p = int($p*100000+0.5)/100000;
	}else{
		$p = sprintf("%e.1",$p);
	}
	$ref->[5] = $p;
}

@genes_over = reverse(@genes_over);
@genes_under = reverse(@genes_under);
my $nOverexpressed=@genes_over;
my $nUnderexpressed=@genes_under;
my $fileID = get_outputID(4);
my $fileTableOver = $fileID+2;
my $fileTableUnder = $fileID+3;

open (OUT1, ">$PATH_OUTPUT/$fileTableOver.txt") or error_message("Cannot write to fileTableOver");
print OUT1 join("\n",@table_headers)."\n";
print OUT1 "!Number_of_over-expressed_genes=$nOverexpressed\n";
print OUT1 "ProbeID\t$header1\t$header2\tFold change\tz-value\tp\tFDR";
if($header2 eq "median"){ print OUT1 "\tSpecificity"; }
print OUT1 "\tGene symbol\tGene name\n";
print OUT1 join("\n",@genes_over)."\n";
close OUT1;

open (OUT1, ">$PATH_OUTPUT/$fileTableUnder.txt") or error_message("Cannot write to fileTableUnder");
print OUT1 join("\n",@table_headers)."\n";
print OUT1 "!Number_of_under-expressed_genes=$nUnderexpressed\n";
print OUT1 "ProbeID\t$header1\t$header2\tFold change\tz-value\tp\tFDR";
if($header2 eq "median"){ print OUT1 "\tSpecificity"; }
print OUT1 "\tGene symbol\tGene name\n";
print OUT1 join("\n",@genes_under)."\n";
close OUT1;

print "<HTML><HEAD><TITLE>ExAtlas: Pairwise comparison</TITLE>\n";
print_header();
print "<SCRIPT language=JavaScript>\n";
print "<!--\n";
print "function plot(probeid) {\n";
print "	document.form_scatterplot.search_term.value=probeid;\n";
print "	document.form_scatterplot.action.value=\"matrix_explore1\";\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_scatterplot.target = \"_BLANK\"+x;\n";
print "	document.form_scatterplot.submit();\n";
print "}\n";
print "function show_genes(direction) {\n";
print "	document.form_scatterplot.action.value=\"show_gene_list\";\n";
print "	document.form_scatterplot.file_genelist.value=$fileTableOver;\n";
print "	if(direction==2){ document.form_scatterplot.file_genelist.value=$fileTableUnder; }\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_scatterplot.target = \"_BLANK\"+x;\n";
print "	document.form_scatterplot.submit();\n";
print "}\n";

print "<!-- end script --></SCRIPT></HEAD>\n";
print "<h2>Scatter-plot of '$header1' vs. '$header2'</h2>\n";
print "<TABLE BORDER=0>\n";
if($nOverexpressed){
	my $x = int(10000*rand());
	print "<TR><TD><b>Overexpressed genes</b><TD>(N = $nOverexpressed)\n";
	print "<TD><INPUT TYPE=button value=\" Show list of genes \" onclick=show_genes(1);>\n";
}
if($nUnderexpressed){
	if($nOverexpressed){ print "<br>"; }
	my $x = int(10000*rand());
	print "<TR><TD><b>Underexpressed genes</b><TD>(N = $nUnderexpressed)\n";
	print "<TD><INPUT TYPE=button value=\" Show list of genes \" onclick=show_genes(2);>\n";
}
if($nOverexpressed+$nUnderexpressed==0){
	print "<TR><TD><b>No significant genes found!</b>\n";
}
print "</TABLE><p>\n";
print "Click on points in the graph to get information<p>\n";
print "<FORM NAME=form_scatterplot ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
my @span = (0,0);
my @delta = (0,0);
for(my $i=0; $i<2; ++$i){
	$span[$i] = $maxx->[$i] - $minx->[$i];
	if ($span[$i] < 3) { $span[$i] = 3; }
	$delta[$i] = 0.5;
	if ($span[$i] > 5) { $delta[$i] = 1; }
}
open (OUT1, ">$PATH_OUTPUT/$fileID.txt") or error_message("Cannot write to $fileID.txt");
my $WID = 600;
my $HGT = ($WID-90)/$span[0]*$span[1]+90;
if($HGT > 2000){ $HGT = 2000; }
my $scale = $span[0]/($WID - 90);
my $x_center = ($maxx->[0] + $minx->[0])/2 - 25*$scale;
my $y_center = ($maxx->[1] + $minx->[1])/2 - 15*$scale;
my $nobjects = @$points_pairwise + 10000;
print OUT1 "$scale\n$x_center\n$y_center\n$nobjects\n";
# Axis horiz & vertical
my $border = 6*$scale;
my $x1 = $minx->[0] - $border;
my $x2 = $maxx->[0] + $border;
my $y1 = $minx->[1] - $border;
my $y2 = $maxx->[1] + $border;
print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x2 $y1\n";
print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x1 $y2\n";
# Make ticks on the horizontal axis
my $y3 = $y1-5*$scale;
my $y4 = $y1-15*$scale;
for (my $i = 0; $i <= $maxx->[0]/$delta[0]; ++$i) {
	my $x3 = $i*$delta[0];
	if($x3 >= $minx->[0]){
		print OUT1 "$LINE $black $thin $solid  $x3 $y1 $x3 $y3\n";
		print OUT1 "$TEXT $black $smallFont $x3 $y4 $x3\n";
	}
}
for (my $i = 1; $i < -$minx->[0]/$delta[0]; ++$i) {
	my $x3 = -$i*$delta[0];
	print OUT1 "$LINE $black $thin $solid  $x3 $y1 $x3 $y3\n";
	print OUT1 "$TEXT $black $smallFont $x3 $y4 $x3\n";
}
my $x5 = ($maxx->[0] + $minx->[0])/2;
my $y5 = $minx->[1]-40*$scale;
print OUT1 "$TEXT $black $largeFont $x5 $y5 $header2\n";

# Make ticks on the vertical axis
my $x3 = $x1-5*$scale;
my $x4 = $x1-20*$scale;
for (my $i = 0; $i <= $maxx->[1]/$delta[1]; ++$i) {
	$y3 = $i*$delta[1];
	if($y3 < $minx->[1]){ next; }
	$y4 = $y3+2*$scale;
	print OUT1 "$LINE $black $thin $solid  $x1 $y3 $x3 $y3\n";
	print OUT1 "$TEXT $black $smallFont $x4 $y4 $y3\n";
}
for (my $i = 1; $i < -$minx->[1]/$delta[1]; ++$i) {
	$y3 = -$i*$delta[1];
	if($y3 > $maxx->[1]){ next; }
	$y4 = $y3+4*$scale;
	print OUT1 "$LINE $black $thin $solid  $x1 $y3 $x3 $y3\n";
	print OUT1 "$TEXT $black $smallFont $x4 $y4 $y3\n";
}
$x5 = $minx->[0]-47*$scale;
$y5 = ($maxx->[1] + $minx->[1])/2-length($header1)/2*6*$scale;
print OUT1 "$TEXT $black $vertLargeFont $x5 $y5 $header1\n";

# Print points
my @listPlot=();
foreach my $ref (sort {$a->[0]<=>$b->[0]} @$points_pairwise){
	my ($z,$x1,$y1,$rank,$FDR,$p,$symbol,$name,$zspec,$id) = @$ref;
	my $logratio = $y1-$x1;
	my $color = $gray;
	#if($z>0 && abs($logratio) >= $thresh && ($FDR_pvalue eq "FDR" && $FDR <= $FDR_thresh || $FDR_pvalue ne "FDR" && $p <= $FDR_thresh)){

	if($z==0 || ($FDR_pvalue eq "FDR" && $FDR > $FDR_thresh || $FDR_pvalue ne "FDR" && $p > $FDR_thresh) || abs($logratio) < $thresh){
		if(abs($logratio) >= $thresh && $rank < $MAX_GENES){ push(@listPlot,[$id,$x1,$y1]); } 
		print OUT1 "$CIRCLE $color $radius $color $x1 $y1\n";
	}else{
		if($logratio>0){ $color = $red; }
		else{ $color = $green; }
		print OUT1 "$CIRCLE $black $radius $color $x1 $y1\n";
		if($rank < $MAX_GENES){ push(@listPlot,[$id,$x1,$y1]); } 
	}
}
close (OUT1);

# Create gif file
my $gif_file = $fileID+1;
$gif_file .= ".gif";
system("$PATH_BIN/togif","$PATH_OUTPUT/$fileID.txt","$PATH_OUTPUT/$gif_file","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");
print "<img src=\"$HOME_ADDRESS/output/$gif_file\" border=0  useMap=#myMap align=center><map NAME='myMap'>\n";
foreach my $record (@listPlot){
	my ($id,$x,$y) = @$record;
	$x = int($WID/2 + ($x-$x_center)/$scale);
	$y = int($HGT/2 - ($y-$y_center)/$scale);
	my $x1 = $x-3;
	my $x2 = $x+3;
	my $y1 = $y-3;
	my $y2 = $y+3;
	print "<area shape=\"rect\" href=\"JavaScript: plot('$id');\" coords=\"$x1,$y1,$x2,$y2\">\n";
}
print "</map><img src=\"$HOME_ADDRESS/images/legend.gif\" align=center><p>\n";
print "<B>Significance criteria: </B> FDR &le; $FDR_thresh and fold-change &le; $foldChange_thresh<p>";
print "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
print "<INPUT NAME=action TYPE=hidden VALUE=matrix_explore1>\n";
print "<INPUT NAME=analysis TYPE=hidden VALUE=search>\n";
print "<INPUT NAME=category TYPE=hidden VALUE=0>\n";
print "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
print "<INPUT NAME=description_matrix TYPE=hidden VALUE=\"$description_matrix\">\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT NAME=search_term TYPE=hidden>\n";
print "<INPUT NAME=file_genelist TYPE=hidden>\n";
if($legendID){ print "<INPUT NAME=legendID TYPE=hidden VALUE=$legendID>\n"; }
print "</FORM><p>\n";
print "<HR NOSHADE></HR>\n";
print "<INPUT TYPE=button VALUE=\"    Close window   \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#*************************************
sub  show_gene_list
#*************************************
{
my $organismID = $hashInput{"organismID"};
my $file_genelist = $hashInput{"file_genelist"};

my ($file_matrix,$description_matrix,$header1,$header2,$overunder,$FDR_thresh,$foldChange_thresh,$nGenes,$legendID);
my @gene_list;
open(INFO,"<$PATH_OUTPUT/$file_genelist.txt") or error_message("File genelist not found!");
while(my $line = <INFO>){
	chop $line;
	if($line !~ /^!/){ last; }
	my ($key,$value) = split("=",$line);
	if($key eq "!Matrix_name"){ $file_matrix=$value; }
	elsif($key eq "!Matrix_description"){ $description_matrix=$value; }
	elsif($key eq "!Column_name"){ $header1=$value; }
	elsif($key eq "!Comparing_with"){ $header2=$value; }
	elsif($key eq "!FDR_threshold"){ $FDR_thresh=$value; }
	elsif($key eq "!Fold_change_threshold"){ $foldChange_thresh=$value; }
	elsif($key eq "!legendID"){ $legendID=$value; }
	elsif($key =~ /^!Number_of_/){
		$nGenes=$value;
		$key =~ s/^!Number_of_//; $key =~ s/_genes$//; $overunder=$key;
	}
}
while(my $line = <INFO>){
	chop $line;
	if($line){ push(@gene_list,$line); }
}

my @file_list;
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!");
while(my $line = <INFO>){
	chop $line;
	my %hash;
	read_config_line($line,\%hash);
	my $descr = $hash{"description"};
	my $org = $hash{"organismID"};
	my @items = split(/[=\t]/,$line);
	if($items[0] =~ /^type_/){
		push(@file_list,$items[1]);
	}
}
close INFO;
my @geneset_list = get_geneset_list();
filter_list_by_organism(\@geneset_list, $organismID);
my ($items,$descriptions) = get_array_lists(\@geneset_list);
my @copy_file_list = @geneset_list;
splice(@copy_file_list,0,0,["--- New file ---",""]);
my $file_list="";
foreach my $name (@file_list){
	if(!$file_list){ $file_list = "\"".$name."\""; }
	else{ $file_list .= ",\"".$name."\""; }
}
print "<HTML><HEAD><TITLE>ExAtlas: Table of $overunder genes</TITLE>\n";
if(@geneset_list && $nGenes>=5){
	print get_header("update_description();");
}else{
	print get_header()."\n";
}
print "<SCRIPT language=JavaScript>\n";
print "<!--\n";
print "geneset_list = new Array($items);\n";
print "geneset_description = new Array($descriptions);\n";
print "function plot(probeid) {\n";
print "	document.form_show_genes.search_term.value=probeid;\n";
print "	document.form_show_genes.action.value=\"matrix_explore1\";\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_show_genes.target = \"_BLANK\"+x;\n";
print "	document.form_show_genes.submit();\n";
print "}\n";
print "function save_gene_set() {\n";
print "	if(!document.form_save_geneset.geneset_name_new.value){\n";
print "		alert(\"Please, enter geneset name\"); return false;\n";
print "	}\n";
print "	var file = document.form_save_geneset.copy_file.options[document.form_save_geneset.copy_file.selectedIndex].value;\n";
print "	if(document.form_save_geneset.copy_file.selectedIndex==0){\n";
print "		file = prompt(\"Provide name of a new file where to copy selected items\");\n";
print "		if(!file){ return(false); }\n";
print "		var file1=file;\n";
print "		if(file.search(/\\.txt\$/)>=0){\n";
print "			if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "			file1=file.substring(0,file.length-4);\n";
print "		}\n";
print "		if(file1.search(/^[-\\w]+\$/)<0){\n";
print "			alert(\"File name should have no spaces or special characters\");\n";
print "			return(false);\n";
print "		}\n";
print "		if(file.search(/^public/i) >= 0){\n";
print "			alert(\"File name cannot start with 'public'\");\n";
print "			return(false);\n";
print "		}\n";
print "		for(i=0; i<file_list.length; ++i){\n";
print "			if(file == file_list[i]){\n";
print "				alert(\"File with this name already exists\"); return false;\n";
print "			}\n";
print "		}\n";
print "		var descrip = document.form_save_geneset.description_copy_file.value;\n";
print "		if(descrip.search(/\\=|\\&/) >= 0){\n";
print "			alert(\"Description should not include character \'=\' or \'&\'\");\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	if(document.form_save_geneset.geneset_name_new.value.search(/\\=|\\&/) >= 0 || document.form_save_geneset.geneset_description_new.value.search(/\\=|\\&/) >= 0){\n";
print "		alert(\"Geneset name and description should not include character \'=\' or \'&\'\");\n";
print "		return false;\n";
print "	}\n";
print "	document.form_save_geneset.copy_to_geneset.value = file;\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_save_geneset.target = \"_BLANK\"+x;\n";
print "	document.form_save_geneset.submit();\n";
print "}\n";
if(@geneset_list){
	print "file_list = new Array($file_list);\n";
	my ($items1,$descriptions1) = get_array_lists(\@copy_file_list);
	print "copy_file_list = new Array($items1);\n";
	print "copy_file_description = new Array($descriptions1);\n";
	print "function update_description() {\n";
	print "	var index;\n";
	print "	index = document.form_save_geneset.copy_file.selectedIndex;\n";
	print "	document.form_save_geneset.description_copy_file.value = copy_file_description[index];\n";
	print "	index = document.form_show_genes.file_geneset1.selectedIndex;\n";
	print "	document.form_show_genes.description_geneset1.value = geneset_description[index];\n";
	print "}\n";
	print "function geneset_overlap(file) {\n";
	print "	document.form_show_genes.upload_geneset.value=file;\n";
	print "	document.form_show_genes.action.value=\"geneset_overlap\";\n";
	print "	var x = Math.round(Math.random()*10000);\n";
	print "	document.form_show_genes.target = \"_BLANK\"+x;\n";
	print "	document.form_show_genes.submit();\n";
	print "}\n";
}
print "<!-- end script --></SCRIPT></HEAD>\n";
print "<H2>Table of $overunder genes in '$header1' vs. '$header2'</H2>\n";
print "Number of $overunder genes = $nGenes<br>\n";
my $x = int(10000*rand());
print "<a href=\"$HOME_ADDRESS/output/$file_genelist.txt\" target=_BLANK$x>Table as tab-delimited text</a><p>\n";
print "<FORM NAME=form_show_genes ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
if(@geneset_list && $nGenes>=5){
	my $menu_text = menu_geneset_overlap(\@geneset_list,$file_genelist);
	print "<TABLE BORDER=0>$menu_text</TABLE><p>\n";
}
print "<TABLE BORDER=0><TR><TD><b>ProbeID<TD><b>$header1<TD><b>$header2<TD WIDTH=100><b>Fold change<TD><b>z-value<TD><b>p<TD><b>FDR\n";
if($header2 eq "median"){ print "<TD><b>Specificity"; }
print "<TD><b>Gene symbol<TD WIDTH=500><b>Gene name\n";
my $count=0;
foreach my $line (@gene_list){
	my @item = split(/\t/,$line);
	$item[0] = "<a href=\"JavaScript:plot('$item[0]');\">$item[0]</a>\n";
	print "<TR><TD>".join("<TD>",@item)."\n";
	if(++$count > $MAX_GENES){ last; }
}
print "</TABLE><p>\n";
if($count > $MAX_GENES){ print "<b>Warning:</b>Table is limited to $MAX_GENES genes. The full list of genes is available <a href=$HOME_ADDRESS/output/$file_genelist.txt target=_BLANK$x>here</a> as text.<p>\n"; }
print "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
print "<INPUT NAME=action TYPE=hidden VALUE=matrix_explore1>\n";
print "<INPUT NAME=analysis TYPE=hidden VALUE=search>\n";
print "<INPUT NAME=category TYPE=hidden VALUE=0>\n";
print "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
print "<INPUT NAME=description_matrix TYPE=hidden VALUE=\"$description_matrix\">\n";
print "<INPUT NAME=description_geneset TYPE=hidden VALUE=\"Overexpressed in $header1 vs. $header2\">\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT NAME=search_term TYPE=hidden>\n";
print "<INPUT NAME=upload_geneset TYPE=hidden VALUE=$file_genelist>\n";
if($legendID){ print "<INPUT NAME=legendID TYPE=hidden VALUE=$legendID>\n"; }
print "</FORM><p>\n";
print "<HR NOSHADE></HR>\n";

print "<p><font size=+2><b>Save genes as geneset</b></font> &nbsp; &nbsp; &nbsp; &nbsp;\n";
print "<FORM NAME=form_save_geneset ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<INPUT TYPE=button VALUE=\"Save Genes\" onClick=save_gene_set(); style=width:250px;>\n";
print "<p><TABLE BORDER=0><TR><TD><b>Geneset name</b><TD><b>Geneset description</b><TD><b>Select geneset file</b><TD><b>File description</b>\n";
print "<TR><TD><INPUT NAME=geneset_name_new style=width:220px;>\n";
print "<TD><INPUT NAME=geneset_description_new style=width:220px;>\n";
print "<TD><select name=copy_file onChange=update_description(); style=width:220px;>\n";
print "<option value=0> $copy_file_list[0]->[0]\n";
for(my $i=1; $i<@copy_file_list; ++$i){ 
	print "<option value=\"$copy_file_list[$i]->[0]\"> $copy_file_list[$i]->[0]\n";
}
print "</select>\n";
print "<TD><INPUT NAME=description_copy_file style=width:220px;>\n";
print "</TABLE><br>\n";
print "If you select \"New file\" you will be prompted for file name, which shound be one-word without special characters (underscore allowed)<p>\n";
print "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
print "<INPUT NAME=action TYPE=hidden VALUE=add_geneset>\n";
print "<INPUT NAME=copy_to_geneset TYPE=hidden>\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT NAME=fileID TYPE=hidden VALUE=$file_genelist>\n";
print "</FORM><p>\n";
print "<HR NOSHADE></HR>\n";
print "<INPUT TYPE=button VALUE=\" Close window \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  make_rankplot
#**************************************
{
my $file_matrix = shift;
my $file_geneset = shift;
my $icol = shift;
my $irow = shift;
my $baseline_name = shift;
my $FDR_thresh = shift;
my $fold_change_thresh = shift;
my $expr_thresh = shift;
my $organismID = shift;
my $organismID1 = shift;
my $gene_attribute_min = shift;

my %symbols=();
my $header1 = get_column_from_matrix($file_matrix,$icol,$baseline_name,$FDR_thresh,$fold_change_thresh,$expr_thresh,\%symbols,$organismID);
$header1 =~ s/ \(\d+\)$//;
my %symbols1=();
get_geneset_from_file($file_geneset,$irow,\%symbols1,$organismID,$organismID1,$gene_attribute_min);
my @sorted = sort {$symbols{$a}->[0]<=>$symbols{$b}->[0]} keys %symbols;
my $NN = @sorted;
my $sum=0;
my @rankplot;
my @response;
my ($max,$xmin,$xmax) = (-1,100000000,-100000000);
my $iii=0;
while($symbols{$sorted[$iii]}->[0] <= $MISSING){ $iii++; }
if($iii > @sorted-1000){ return 0; }
my $window = 500;
if(@sorted - $iii<10000){ $window = int((@sorted-$iii)/20); }
my $iii1 = $iii;
my $imax = $iii+$window-1;
if($imax>@sorted){ $imax=@sorted; }
for(; $iii1<$imax; $iii1++){
	if($symbols1{$sorted[$iii1]}){ $sum++; }
}
my $step = 20;
for(my $i=0; $i<@sorted-$imax; $i++){
	if($symbols1{$sorted[$iii1+$i]}){ $sum++; }
	if($i%$step==0 || $i==@sorted-$imax-1){
		push(@rankplot,$sum/$window);
		if($max < $sum/$window){ $max=$sum/$window; }
		my $x = $symbols{$sorted[$i+int(($iii+$iii1)/2)]}->[0];
		push(@response,$x);
		if($xmax < $x){ $xmax = $x; }
		if($xmin > $x){ $xmin = $x; }
	}
	if($symbols1{$sorted[$iii+$i]}){ $sum--; }
}
if($max < 15/$window){ $max=15/$window; }
my $fileID = get_outputID(2);
my $scriptID = $fileID+1;
my $spany = $max;
my $deltay = 0.005;
if ($spany > 0.05){ $deltay = 0.01; }
if ($spany > 0.1){ $deltay = 0.02; }
if ($spany > 0.2){ $deltay = 0.05; }
if ($spany > 0.5){ $deltay = 0.1; }
my $spany1 = $xmax-$xmin;
if($spany1 < 0.2){ $spany1=0.2; }
my $deltay1 = 0.05;
if ($spany1 > 1){ $deltay1 = 0.2; }
if ($spany1 > 2){ $deltay1 = 0.5; }
if ($spany1 > 5){ $deltay1 = 1; }
if ($spany1 >10){ $deltay1 = 2; }
my $spanx = @rankplot*$step;
my $deltax = 500;
if ($spanx > 5000){ $deltax = 1000; }
if ($spanx > 10000){ $deltax = 2000; }
if ($spanx > 20000){ $deltax = 5000; }
open (OUT1, ">$PATH_OUTPUT/$scriptID.txt") or error_message("Cannot write to file");
my $WID = 700;
my $HGT = 700;
my $scaley = $spany/($HGT/2 - 60);
my $scaley1 = $spany1/($HGT/2 - 60);
my $scalex = $spanx/($WID - 80);
my $x_center = $spanx/2/$scalex - 25;
my $y_center = -30;
my $nobjects = @rankplot*2 + 10000;
print OUT1 "1\n$x_center\n$y_center\n$nobjects\n";
# Axis horiz & vertical
my $x2 = $spanx/$scalex+3;
my $y2 = $spany/$scaley+5;
print OUT1 "$LINE $black $thick2 $solid  0 0 $x2 0\n";
print OUT1 "$LINE $black $thick2 $solid  0 0 0 $y2\n";
# Make ticks on the horizontal axis
my $y3 = -5;
my $y4 = -15;
for (my $i=0; $i<=$spanx/$deltax; ++$i) {
	my $x3 = $i*$deltax;
	my $x4 = $x3/$scalex;
	print OUT1 "$LINE $black $thin $solid  $x4 0 $x4 $y3\n";
	print OUT1 "$TEXT $black $smallFont $x4 $y4 $x3\n";
}
my $x5 = ($WID-120)/2;
my $y5 = -33;
print OUT1 "$TEXT $black $smallFont $x5 $y5 Gene rank by expression in $header1\n";
# Make ticks on the vertical axis
my $x3 = -5;
my $x4 = -20;
for (my $i=0; $i<=$spany/$deltay; ++$i) {
	my $y3 = $i*$deltay;
	my $y4 = $y3/$scaley;
	print OUT1 "$LINE $black $thin $solid  $x3 $y4 0 $y4\n";
	print OUT1 "$TEXT $black $smallFont $x4 $y4 $y3\n";
}
$x5 = -5;
$y5 = $spany/$scaley+15;
print OUT1 "$TEXT $black $largeFont $x5 $y5 Proportion\n";
# Make line for rankplot
for(my $i=1; $i<@rankplot; $i++){
	my $x1 = ($i-1)*$step/$scalex;
	my $x2 = $i*$step/$scalex;
	my $y1 = $rankplot[$i-1]/$scaley;
	my $y2 = $rankplot[$i]/$scaley;
	print OUT1 "$LINE $blue $thick2 $solid  $x1 $y1 $x2 $y2\n";
}

# Axis horiz & vertical for response
my $shift = $HGT/2-10;
$x2 = $spanx/$scalex+3;
$y2 = $spany1/$scaley1+5-$shift;
print OUT1 "$LINE $black $thick2 $solid  0 -$shift $x2 -$shift\n";
print OUT1 "$LINE $black $thick2 $solid  0 -$shift 0 $y2\n";
$y3 = -5-$shift;
$y4 = -15-$shift;
for (my $i=0; $i<=$spanx/$deltax; ++$i) {
	my $x3 = $i*$deltax;
	my $x4 = $x3/$scalex;
	print OUT1 "$LINE $black $thin $solid  $x4 -$shift $x4 $y3\n";
	print OUT1 "$TEXT $black $smallFont $x4 $y4 $x3\n";
}
$x5 = ($WID-120)/2;
$y5 = -33-$shift;
print OUT1 "$TEXT $black $largeFont $x5 $y5 Gene rank\n";
# Make ticks on the vertical axis
$x3 = -5;
$x4 = -20;
my $i1 = 0;
if($xmin<0){
	$i1 = int(-$xmin/$deltay1);
}
for (my $i=0; $i<=$spany1/$deltay1; ++$i) {
	my $y3 = ($i-$i1)*$deltay1;
	my $y4 = ($y3 - $xmin)/$scaley1-$shift;
	if($y3 > $xmax){ last; }
	print OUT1 "$LINE $black $thin $solid  0 $y4 $x3 $y4\n";
	print OUT1 "$TEXT $black $smallFont $x4 $y4 $y3\n";
}
$x5 = -5;
$y5 = $spany/$scaley+15-$shift;
print OUT1 "$TEXT $black $largeFont $x5 $y5 Response\n";
# Make line for response
for(my $i=1; $i<@response; $i++){
	my $x1 = ($i-1)*$step/$scalex;
	my $x2 = $i*$step/$scalex;
	my $y1 = ($response[$i-1] - $xmin)/$scaley1-$shift;
	my $y2 = ($response[$i] - $xmin)/$scaley1-$shift;
	print OUT1 "$LINE $blue $thick2 $solid  $x1 $y1 $x2 $y2\n";
}
close (OUT1);
system("$PATH_BIN/togif","$PATH_OUTPUT/$scriptID.txt","$PATH_OUTPUT/$fileID.gif","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");
return $fileID;
}

#**************************************
sub  get_column_from_matrix
#**************************************
{
my $file_matrix = shift;
my $icol = shift;
my $baseline_name = shift;
my $FDR_thresh = shift;
my $fold_change_thresh = shift;
my $expr_thresh = shift;
my $symbols = shift;
my $organismID = shift;
my $organismID1 = shift;	#optional

my $log10 = log(10);
if($fold_change_thresh > 1){ $fold_change_thresh = log($fold_change_thresh)/$log10; }
my $expr_log_thresh=$MISSING; if($expr_thresh>0){ $expr_log_thresh=log($expr_thresh)/$log10; }
else{ $fold_change_thresh = 0; }
if(!$organismID1){$organismID1=$organismID; }
if($organismID*$organismID1==0 && $organismID+$organismID1 != 0){
	error_message("Bad organism ID in get_column_from_matrix");
}
if(!$symbols || ref($symbols) ne 'HASH'){ error_message("No symbols pointer in get_column_from_matrix"); }
%$symbols = ();
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_anova = "$loginname-anova-$file_matrix";
}
my %hashGenes;
my %hashAlias;
my %geneHomolog=();
if($organismID1 != $organismID){
	get_gene_homolog($organismID1,$organismID,\%geneHomolog);
	get_official_symbols($organismID1,\%hashGenes,\%hashAlias);
}else{
	get_official_symbols($organismID,\%hashGenes,\%hashAlias);
}
open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file");
my $line = <INFO>;
chop $line;
my ($junk,$junk1,@headers) =split(/\t/,$line);
if(!@headers){ error_message("Headers not found in anova file"); }
my $nCol=0;
while($nCol<@headers && $headers[$nCol] !~ /^Var\(/){ ++$nCol; }
my $baseline=0;
if($baseline_name ne "Median profile"){
	for(my $i=0; $i<$nCol; ++$i){
		if($headers[$i]=~/ \(\d+\)$/ && $headers[$i] !~ /^Mean/){
			my @items = split(/ \(/,$headers[$i]);
			my $n = pop(@items);
			my $name = join(" (",@items);
			$n =~ s/\)$//;
			$headers[$i] = $name;
		}else{
			$headers[$i] =~ s/^Mean\(//;
			$headers[$i] =~ s/\)$//;
		}
		if($baseline_name eq $headers[$i]){ $baseline=$i+1; }
	}
	if(!$baseline){ error_message("Baseline header not found"); }
}
#printf "$baseline_name $baseline<br>\n";
while(my $line = <INFO>){
	chop $line;
	my ($id,$expr,@data1) =split(/\t/,$line);
	my $FDR = $data1[$nCol+6];
	my $symbol = $data1[$nCol+8];
	if(!$symbol || $FDR_thresh>0 && $FDR>$FDR_thresh){
		next;
	}
	if($expr_thresh>0 && $expr<$expr_log_thresh){ next; }
	if(!$hashGenes{$symbol}){
		my $symbol1 = $hashAlias{$symbol};
		if($symbol1){ $symbol=$symbol1; }
	}
	if($organismID1 != $organismID){
		my $symbol1 = $geneHomolog{$symbol};
		if($symbol1){
			$symbol = $symbol1;
		}else{
			next;
		}
	}
	my $F = $data1[$nCol+4];
	my $x = $data1[$icol-1];
	splice(@data1,$nCol);
	if($fold_change_thresh>0){
		my @sorted = sort {$a<=>$b} @data1;
		while($sorted[0]<=$MISSING){ shift(@sorted); }
		if(@sorted>1){
			my $logratio = $sorted[@sorted-1] - $sorted[0];
			if($logratio < $fold_change_thresh){ next; }
		}
	}
	my $median;
	if($baseline==0){
		$median = median(\@data1);
	}else{
		$median = $data1[$baseline-1];
	}
	if($x>$MISSING && $median>$MISSING){
		$x -= $median;
	}
	my $ref = $symbols->{$symbol};
	if($ref){
		if($F > $ref->[1] && $x>$MISSING || $ref->[0]==$MISSING){
			$ref->[1] = $F;
			$ref->[0] = $x;
		}
	}else{
		$symbols->{$symbol} = [$x,$F];
	}
}
close INFO;
return $headers[$icol-1];
}

#**************************************
sub  get_geneset_from_file
#**************************************
{
my $file_geneset = shift;
my $irow = shift;
my $symbols = shift;
my $organismID = shift;
my $organismID1 = shift;
if(!$organismID1){ $organismID1=$organismID; }
my $gene_attribute_min = shift;

if(!$symbols || ref($symbols) ne 'HASH'){ error_message("No symbols pointer in get_geneset_from_file"); }
%$symbols = ();
my $file_geneset_full = "$PATH_DATA/$file_geneset";
if($file_geneset =~ /^TEMPORARY!-/i){
	$file_geneset =~ s/^TEMPORARY!-//i;
	$file_geneset_full = "$PATH_OUTPUT/$file_geneset";
}
elsif($file_geneset !~ /^public-/){
	$file_geneset_full = "$PATH_DATA/$loginname-$file_geneset";
}
my @geneSetList;
my %hashGenes;
my %hashAlias;
my %geneHomolog=();
if($organismID1 != $organismID){
	get_gene_homolog($organismID1,$organismID,\%geneHomolog);
	get_official_symbols($organismID1,\%hashGenes,\%hashAlias);
}else{
	get_official_symbols($organismID,\%hashGenes,\%hashAlias);
}
open(INFO,'<',$file_geneset_full) or error_message("Cannot open geneset file");
my $count=0;
my $header;
while(my $line = <INFO>){
	if($line =~ /^!|^#/){ next; }
	chop $line;
	if(!$line){ next; }
	my($setName,$descrip,@symbol_list)=split(/\t/, $line);
	if(!$setName){ next; }
	if(++$count == $irow){
		$header = $setName;
		$line = <INFO>;
		my($junk1,$descrip1,@attribute)=split(/\t/, $line);
		for(my $is=0; $is<@symbol_list; $is++){
			if($gene_attribute_min>0 && $attribute[$is]>0 && $attribute[$is]<$gene_attribute_min){ next; }
			my $symbol = $symbol_list[$is];
			if(!$hashGenes{$symbol}){
				my $symbol1 = $hashAlias{$symbol};
				if($symbol1){ $symbol=$symbol1; }
			}
			if($organismID1 != $organismID){
				my $symbol1 = $geneHomolog{$symbol};
				if($symbol1){
					$symbol = $symbol1;
				}else{
					next;
				}
			}
			$symbols->{$symbol}=1;
		}
		last;
	}
}
return $header;
}

#**************************************
sub  get_gene_attributes
#**************************************
{
my $file_geneset = shift;
my $irow = shift;
my $symbol_ref = shift;
my $attrib_ref = shift;
my $organismID = shift;
my $organismID1 = shift;
if(!$organismID1){ $organismID1=$organismID; }

if(!$symbol_ref || ref($symbol_ref) ne 'ARRAY'){ error_message("No symbol pointer in get_gene_attributes"); }
if(!$attrib_ref || ref($attrib_ref) ne 'ARRAY'){ error_message("No attrib pointer in get_gene_attributes"); }
@$attrib_ref = ();
my %hashGeneOrder;
for(my $i=0; $i<@$symbol_ref; $i++){
	my $symbol = $symbol_ref->[$i];
	$hashGeneOrder{uc($symbol)}=$i+1;
}
my @attrib_names;
my $file_geneset_full = "$PATH_DATA/$file_geneset";
if($file_geneset =~ /^TEMPORARY!-/i){
	$file_geneset =~ s/^TEMPORARY!-//i;
	$file_geneset_full = "$PATH_OUTPUT/$file_geneset";
}
elsif($file_geneset !~ /^public-/){
	$file_geneset_full = "$PATH_DATA/$loginname-$file_geneset";
}
my @geneSetList;
my %hashGenes;
my %hashAlias;
my %geneHomolog=();
if($organismID1 != $organismID){
	get_gene_homolog($organismID1,$organismID,\%geneHomolog);
	get_official_symbols($organismID1,\%hashGenes,\%hashAlias);
}else{
	get_official_symbols($organismID,\%hashGenes,\%hashAlias);
}
open(INFO,'<',$file_geneset_full) or return NULL;
my $count=0;
while(my $line = <INFO>){
	if($line =~ /^!|^#/){ next; }
	chop $line;
	if(!$line){ next; }
	my($setName,$descrip,@symbol_list)=split(/\t/, $line);
	if(!$setName){ next; }
	if(++$count == $irow){
		my @gene_order;
		for(my $i=0; $i<@symbol_list; $i++){
			my $symbol = $symbol_list[$i];
			if(!$hashGenes{$symbol}){
				my $symbol1 = $hashAlias{$symbol};
				if($symbol1){ $symbol=$symbol1; }
			}
			my $ii;
			if($organismID1 != $organismID){
				my $symbol1 = $geneHomolog{$symbol};
				if($symbol1){
					$ii = $hashGeneOrder{uc($symbol1)};
					if($ii){ $gene_order[$ii-1] = $i+1; }
				}
			}
			if(!$ii){
				$ii = $hashGeneOrder{uc($symbol)};
				if($ii){ $gene_order[$ii-1] = $i+1; }
			}
		}
		while($line = <INFO>){
			chop $line;
			my($blank,$attribName,@attributes)=split(/\t/, $line);
			if($blank){ last; }
			my @attrib1;
			for(my $j=0; $j<@$symbol_ref; $j++){
				if($gene_order[$j]){
					$attrib1[$j] = $attributes[$gene_order[$j]-1];
				}
			}
			push(@attrib_names,$attribName);
			push(@$attrib_ref,\@attrib1);
		}
		last;
	}
}
return @attrib_names;
}

#**************************************
sub  make_scatterplot_twofiles
#**************************************
{
my $file_matrix1 = shift;
my $file_matrix2 = shift;
my $icol = shift;
my $irow = shift;
my $baseline_name1 = shift;
my $baseline_name2 = shift;
my $FDR1 = shift;
my $FDR2 = shift;
my $fold_change1 = shift;
my $fold_change2 = shift;
my $expr_thresh1 = shift;
my $expr_thresh2 = shift;
my $organismID = shift;
my $organismID1 = shift;
my $hashGenes = shift;

if(!$FDR1){ $FDR1=0.05; }
if(!$FDR2){ $FDR2=0.05; }
if(!$fold_change1){ $fold_change1=1.5; }
if(!$fold_change2){ $fold_change2=1.5; }
my %symbols1=();
my $header1 = get_column_from_matrix($file_matrix1,$icol,$baseline_name1,$FDR1,$fold_change1,$expr_thresh1,\%symbols1,$organismID);
my %symbols2=();
my $header2 = get_column_from_matrix($file_matrix2,$irow,$baseline_name2,$FDR2,$fold_change2,$expr_thresh2,\%symbols2,$organismID,$organismID1);
$header1 =~ s/ \(\d+\)$//;
$header2 =~ s/ \(\d+\)$//;
my @data_points;
my @rank_points;
my @minx=(100000000,100000000);
my @maxx=(-100000000,-100000000);
foreach my $symbol (keys %symbols1){
	my $ref1 = $symbols1{$symbol};
	my $ref2 = $symbols2{$symbol};
	if(!$ref2 || $ref1->[0]==$MISSING || $ref2->[0]==$MISSING){
		next;
	}
	my $significant=0;
	if($hashGenes->{$symbol}){
		$significant=1;
	}
	my @x = ($ref1->[0],$ref2->[0],$significant);
	for(my $i=0; $i<2; $i++){
		if($minx[$i]>$x[$i]){ $minx[$i]=$x[$i]; }
		if($maxx[$i]<$x[$i]){ $maxx[$i]=$x[$i]; }
	}
	push(@data_points,\@x);
	my @y = @x;
	push(@rank_points,\@y);
}
my $N = @rank_points;
for(my $i=0; $i<2; $i++){
	my @rank = sort {$rank_points[$a]->[$i]<=>$rank_points[$b]->[$i]} (0..($N-1));
	for(my $j=0; $j<$N; $j++){
		$rank_points[$rank[$j]]->[$i] = $j;
	}
}
my $fileID = get_outputID(4);
my $scriptID = $fileID+2;
my @span = (0,0);
my @delta = (0,0);
for(my $i=0; $i<2; ++$i){
	$span[$i] = $maxx[$i] - $minx[$i];
	if ($span[$i] < 1) { $span[$i] = 1; }
	$delta[$i] = 0.2;
	if ($span[$i] > 2) { $delta[$i] = 0.5; }
	if ($span[$i] > 4) { $delta[$i] = 1; }
}
open (OUT1, ">$PATH_OUTPUT/$scriptID.txt") or error_message("Cannot write to file");
my $WID = 500;
my $HGT = 500;
my $scalex = $span[0]/($WID - 90);
my $scaley = $span[1]/($HGT - 90);
my $x_center = $span[0]/2/$scalex - 25;
my $y_center = $span[1]/2/$scaley - 15;
my $nobjects = @data_points*2 + 10000;
print OUT1 "1\n$x_center\n$y_center\n$nobjects\n";
# Axis horiz & vertical
my $x2 = $span[0]/$scalex+5;
my $y2 = $span[1]/$scaley+5;
print OUT1 "$LINE $black $thin $solid  0 0 $x2 0\n";
print OUT1 "$LINE $black $thin $solid  0 0 0 $y2\n";
# Make ticks on the horizontal axis
my $y3 = -5;
my $y4 = -15;
my $i1 = 0;
if($minx[0]<0){
	$i1 = int(-$minx[0]/$delta[0]);
}
for (my $i = 0; $i <= $span[0]/$delta[0]; ++$i) {
	my $x3 = int(($i-$i1)*$delta[0]);
	my $x4 = int(($x3-$minx[0])/$scalex);
	if($x4<=$x2){
		print OUT1 "$LINE $black $thin $solid  $x4 0 $x4 $y3\n";
		print OUT1 "$TEXT $black $smallFont $x4 $y4 $x3\n";
	}
}
my $x5 = $span[0]/2/$scalex;
my $y5 = -33;
print OUT1 "$TEXT $black $largeFont $x5 $y5 $header1\n";
# Make ticks on the vertical axis
my $x3 = -5;
my $x4 = -20;
$i1 = 0;
if($minx[1]<0){
	$i1 = int(-$minx[1]/$delta[1]);
}
for (my $i = 0; $i <= $span[1]/$delta[1]; ++$i) {
	$y3 = ($i-$i1)*$delta[1];
	$y4 = ($y3-$minx[1])/$scaley;
	if($y4<=$y2){
		print OUT1 "$LINE $black $thin $solid  $x3 $y4 0 $y4\n";
		print OUT1 "$TEXT $black $smallFont $x4 $y4 $y3\n";
	}
}
$x5 = -40;
$y5 = $span[1]/2/$scaley;
print OUT1 "$TEXT $black $vertLargeFont $x5 $y5 $header2\n";
# Print points
foreach my $ref (@data_points){
	my ($x1,$y1,$sig) = @$ref;
	$x1 = int(($x1-$minx[0])/$scalex);
	$y1 = int(($y1-$minx[1])/$scaley);
	if($sig){ print OUT1 "$LINE $magenta $thick3 $solid  $x1 $y1 $x1 $y1\n"; }
	else{ print OUT1 "$LINE $blue $thick2 $solid  $x1 $y1 $x1 $y1\n"; }
}
close (OUT1);
# Create gif file
system("$PATH_BIN/togif","$PATH_OUTPUT/$scriptID.txt","$PATH_OUTPUT/$fileID.gif","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");

$scriptID++;
open (OUT1, ">$PATH_OUTPUT/$scriptID.txt") or error_message("Cannot write to file");
$scalex = $N/($WID - 90);
$x_center = $N/2/$scalex - 25;
$y_center = $N/2/$scalex - 15;
my $nobjects = $N*2 + 10000;
print OUT1 "1\n$x_center\n$y_center\n$nobjects\n";
# Axis horiz & vertical
$x2 = $N/$scalex+5;
$y2 = $N/$scaley+5;
print OUT1 "$LINE $black $thin $solid  0 0 $x2 0\n";
print OUT1 "$LINE $black $thin $solid  0 0 0 $y2\n";
$x5 = $N/2/$scalex;
$y5 = -33;
print OUT1 "$TEXT $black $largeFont $x5 $y5 $header1\n";
$x5 = -40;
$y5 = $span[1]/2/$scaley;
print OUT1 "$TEXT $black $vertLargeFont $x5 $y5 $header2\n";
# Print points
foreach my $ref (@rank_points){
	my ($x1,$y1,$sig) = @$ref;
	$x1 = int($x1/$scalex);
	$y1 = int($y1/$scalex);
	if($sig){ print OUT1 "$LINE $magenta $thick2 $solid  $x1 $y1 $x1 $y1\n"; }
	else{ print OUT1 "$LINE $blue 1 $solid  $x1 $y1 $x1 $y1\n"; }
}
close (OUT1);
# Create gif file
my $fileID1 = $fileID+1;
system("$PATH_BIN/togif","$PATH_OUTPUT/$scriptID.txt","$PATH_OUTPUT/$fileID1.gif","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");
return $fileID;
}

#**************************************
sub  make_venn_diagram
#**************************************
{
my $file_geneset1 = shift;
my $file_geneset2 = shift;
my $icol = shift;
my $irow = shift;
my $organismID = shift;
my $organismID1 = shift;
my $N12 = shift;

my %symbols1=();
my $header1 = get_geneset_from_file($file_geneset1,$icol,\%symbols1,$organismID);
my %symbols2=();
my $header2 = get_geneset_from_file($file_geneset2,$irow,\%symbols2,$organismID,$organismID1);
my $N1 = keys %symbols1;
my $N2 = keys %symbols2;

my ($Nmax,$Nmin,$N11,$N22) = ($N1,$N2,$N1-$N12,$N2-$N12);
my ($r1,$r2) = (150, int(150*sqrt($N2/$N1)));
my ($Dmin,$Dmax) = ($r1-$r2,$r1+$r2);
if($Nmax<$N2){
	($Nmax,$Nmin) = ($N2,$N1);
	($r1,$r2) = (int(150*sqrt($N1/$N2)),150);
	($Dmin,$Dmax) = ($r2-$r1,$r1+$r2);
}
my $prop = 1-acos(($Nmin-2*$N12)/$Nmin)/3.141593;
my $d = $Dmin+int(($Dmax-$Dmin)*$prop);

my $delta1 = int($d*$r1/$Dmax);
my $delta2 = int($d*$r2/$Dmax);

my $fileID = get_outputID(2);
my $scriptID = $fileID+1;
open (OUT1, ">$PATH_OUTPUT/$scriptID.txt") or error_message("Cannot write to file");
my $WID = 700;
my $HGT = 500;
my $x_center = 0;
my $y_center = 0;
my $nobjects = 10000;
print OUT1 "1\n$x_center\n$y_center\n$nobjects\n";
print OUT1 "$CIRCLE $black $r1 -1 -$delta1 0\n";
print OUT1 "$CIRCLE $black $r2 -1 $delta2 0\n";
my $x1 = -150;
my $y1 = 200;
my $y2 = 180;
print OUT1 "$TEXT $black $largeFont $x1 $y1 $header1\n";
print OUT1 "$TEXT $black $smallFont $x1 $y2 N = $N1\n";
$x1 = 150;
print OUT1 "$TEXT $black $largeFont $x1 $y1 $header2\n";
print OUT1 "$TEXT $black $smallFont $x1 $y2 N = $N2\n";
$x1 = ($delta2-$delta1-$r2-$r1)/2;
print OUT1 "$FILL 52 $x1 0\n";
print OUT1 "$TEXT $black $smallFont $x1 0 $N11\n";
$x1 = ($delta2-$delta1+$r2+$r1)/2;
print OUT1 "$FILL 59 $x1 0\n";
print OUT1 "$TEXT $black $smallFont $x1 0 $N22\n";
if($prop < 0.98){
	$x1 = ($delta2-$delta1+$r1-$r2)/2;
	print OUT1 "$FILL $ltgreen $x1 0\n";
	print OUT1 "$TEXT $black $smallFont $x1 0 $N12\n";
}
close (OUT1);
# Create gif file
system("$PATH_BIN/togif","$PATH_OUTPUT/$scriptID.txt","$PATH_OUTPUT/$fileID.gif","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");
return $fileID;
}

#**************************************
sub   make_PCA_output_page
#**************************************
{
my $logFileID = shift;
my $data_sent = shift;
my $webPageID = $hashInput{"logFileID"}+1;
my $PCAfileID = $webPageID+1;

my $matrixFileID = $data_sent->[0];
my $nCol = $data_sent->[1];
my $nRow = $data_sent->[2];
my $cluster_type = $data_sent->[3];   # 2=cols, 3=both

my $MAXCOLUMNS = 12;
my @COLOR_GRADE = (36,37,38,39,40,41,42,43);

my $pca_cluster = $hashInput{"pca_cluster"};
my $legendID = $hashInput{"legendID"};
my $npc = int(1.7*sqrt($nCol));
if($npc > 7){ $npc=8; }
my $return = system("$PATH_BIN/pca","$PATH_OUTPUT/$matrixFileID.txt","$PATH_OUTPUT/$PCAfileID.txt","$npc","0","V","1");
if($return){ error_message("PCA crashed",$logFileID); }
open(INFO,"<$PATH_OUTPUT/$PCAfileID.txt") or error_message("Cannot open file $PCAfileID.txt",$logFileID);
# Read eigenvalues
while(my $line = <INFO>){
	if($line =~ /^Eigenvalues/){ last; }
}
my @PClist;
my @PCstdev;
my $sumVar=0;
while(my $line = <INFO>){
	chop($line);
	if(!$line){ last; }
	if($line>0){
		push(@PClist, $line);
		$sumVar += $line;
		push(@PCstdev,sqrt($line));
	}

}
if(!@PClist){ error_message("PCA failed! Use filtering to remove zero columns and rows",$logFileID); }
my @minx;
my @maxx;
my @miny;
my @maxy;
my $pc_toplot=$npc;
my $sum = 0;
for(my $i=0; $i<$npc; ++$i){
	push(@minx,1000000);
	push(@maxx,-1000000);
	push(@miny,1000000);
	push(@maxy,-1000000);
}
# Read column projections
my @lines_col;
my @pcColumns;
my @columnidList;
my $maxColumnHeaderLength=0;
my $line = <INFO>;	# skip header line
while(my $line = <INFO>){
	chop $line;
	if(!$line){ last; }
	my ($columnid, @coord)=split(/\t/,$line);
	my $len = length($columnid);
	if($maxColumnHeaderLength<$len){ $maxColumnHeaderLength=$len; }
	push(@columnidList,$columnid); 
	for(my $i=0;$i<$pc_toplot; ++$i){
		if($miny[$i] > $coord[$i]) { $miny[$i] = $coord[$i]; }
		if($maxy[$i] < $coord[$i]) { $maxy[$i] = $coord[$i]; }
		push(@{$pcColumns[$i]}, $coord[$i]); 
	}
	push(@lines_col, $line);
}
my $fontLabel = $vertSmallFont;
my $charSpacing = 6;
if(@columnidList > $MAXCOLUMNS){
	$fontLabel = $vertTinyFont;
	$charSpacing = 5;
}
# Read row-point projections
my @sumPC;
my @varPC;
my @sumproduct;
my @lines_row;
$line = <INFO>; #skip headers
while(my $line = <INFO>){
	chop $line;
	if(!$line){ last; }
	my ($rowid, $mean1, @coord)=split(/\t/,$line);
	for(my $i=0;$i<$pc_toplot; ++$i){
		if($minx[$i] > $coord[$i]) { $minx[$i] = $coord[$i]; }
		if($maxx[$i] < $coord[$i]) { $maxx[$i] = $coord[$i]; }
	}
	push(@lines_row, $line); 
}
close INFO;

open(INFO,"<$PATH_OUTPUT/$matrixFileID.txt") or error_message("Cannot open file $matrixFileID.txt",$logFileID);
$line = <INFO>;
chop $line;
my ($junk,@headers) = split(/\t/,$line);
my %hashCoord;
my @average;
while(my $line = <INFO>){
	chop $line;
	if(!$line){ last; }
	my ($id,@data1) = split(/\t/,$line);
	$id =~ s/\"//g;
	$hashCoord{$id} = \@data1;
	push(@average,average(\@data1));
}
close INFO;
my @spany;
my @spanx;
for(my $i=0;$i<$pc_toplot; ++$i){
	$spany[$i] = $maxy[$i]-$miny[$i];
	$spanx[$i] = $maxx[$i]-$minx[$i];
}

my $log10=log(10);
my %hashCluster;
my $outputPCcluster = get_outputID(1);
$outputPCcluster .= ".txt";

my $cluster_correlation = $hashInput{"cluster_correlation"};
my $cluster_fold_thresh = $hashInput{"cluster_fold_thresh"};
my $logthreshold = log($cluster_fold_thresh)/$log10;
if($pca_cluster eq "on"){
	open(OUT1, ">$PATH_OUTPUT/$outputPCcluster");
	print OUT1 "GeneID\tLog10Change\tCorrelation\tPCnumber\tDirection\n";
	foreach my $featureid (keys %hashCoord){
		my @coord = @{$hashCoord{$featureid}};
		my @coordCentered=();
		my $median = median(\@coord);
		for(my $i=0; $i<$nCol; ++$i){
			if($coord[$i]>$MISSING){ $coord[$i] -= $median; }
			push(@coordCentered,$coord[$i]);
		}
		my @results;
		for(my $ipc=0; $ipc < $pc_toplot; ++$ipc){
			my ($change,$correl) = get_change($pcColumns[$ipc],\@coordCentered);
			if(abs($change) < $logthreshold || abs($correl) < $cluster_correlation){ next; }
			my $ref = $hashCluster{$featureid};
			if(!$ref){
				$hashCluster{$featureid} = [$ipc,$change,$correl];
			}elsif(abs($correl) > abs($ref->[2])){
				@$ref = ($ipc,$change,$correl);
			}
		}
		my $ref = $hashCluster{$featureid};
		if(!$ref){ next; }
		my ($ipc,$change,$correl) = @$ref;
		$change = floor($change*1000+0.5)/1000;
		$correl = floor($correl*1000+0.5)/1000;
		my $PCno = $ipc+1;
		my $direction = "positive";
		if($change < 0){ $direction = "negative"; }
		print OUT1 "$featureid\t$change\t$correl\tPC$PCno\t$direction\n";
	}
	close OUT1;
}
my $minMean = 10000000;
my $maxMean = -10000000;
for(my $i=0; $i<@average; $i++){
	my $x = $average[$i];
	if($minMean > $x){ $minMean = $x; }
	if($maxMean < $x){ $maxMean = $x; }
}
my $spanMean = $maxMean - $minMean;
if(@lines_row > 3000){
	my @array;
	my @new_aver;
	my $ii=0;
	foreach my $line (@lines_row){
		my ($rowid, $rowmean, @coord)=split(/\t/,$line);
		my $dist=0;
		for(my $i=0;$i<$pc_toplot; ++$i){
			if($i==3){ last; }
			$dist += $coord[$i]*$coord[$i]*$PClist[$i];
		}
		push(@array,$dist);
	}
	my $nRow = @lines_row;
	my @new;
	foreach my $ii (sort {$array[$b]<=>$array[$a]} 0..($nRow-1)){
		push(@new, $lines_row[$ii]);
		push(@new_aver, $average[$ii]);
		if(@new==3000){ last; }
	}
	@lines_row = @new;
	@average = @new_aver;
}

my $fileID = get_outputID(4);
open (OUT1, ">$PATH_OUTPUT/$fileID.txt");
my $WID = 660;
my $HGT = 360;
my @scale;
$scale[0] = $spany[0]/250;
$scale[1] = $spany[1]/250;
my $x_center = 330;
my $y_center = 180;
my $nobjects = (@lines_col + @lines_row)*2 + 1000;
print OUT1 "1\n$x_center\n$y_center\n$nobjects\n";
# Axis horiz & vertical
my $SHIFTX = 40;
my $SHIFTY = 40;
my $border = 20;
my $x1 = $SHIFTX-$border;
my $x2 = $x1 + 300;
my $y1 = $SHIFTY-$border;
my $y2 = $y1 + 300;
my $x0 = -$miny[0]/$scale[0]+$SHIFTX;
my $y0 = -$miny[1]/$scale[1]+$SHIFTY;
print OUT1 "$LINE $black 2 $solid  $x1 $y1 $x2 $y1\n";
print OUT1 "$LINE $black 2 $solid  $x1 $y1 $x1 $y2\n";
print OUT1 "$LINE $black $thin $solid  $x0 $y1 $x0 $y2\n";
print OUT1 "$LINE $black $thin $solid  $x1 $y0 $x2 $y0\n";
my $x3 = $x2-20;
my $y3 = $y1-8;
print OUT1 "$TEXT $black $smallFont $x3 $y3 PC1\n";
$x3 = $x1+5;
$y3 = $y2+10;
print OUT1 "$TEXT $black $smallFont $x3 $y3 PC2\n";
# Plot column-points
my $radius=3;
my $iline = 0;
foreach my $line (@lines_col){
	my ($columnid,@coord) = split(/\t/, $line);
	$x1 = ($coord[0]-$miny[0])/$scale[0]+$SHIFTX;
	$y1 = ($coord[1]-$miny[1])/$scale[1]+$SHIFTY;
	$x2 = $x1 + length($columnid)/2*6 - 5;
	if($x2 + length($columnid)/2*6 > ($maxy[0]-$miny[0])/$scale[0]+$SHIFTX){
		$x2 = ($maxy[0]-$miny[0])/$scale[0]+$SHIFTX-length($columnid)/2*6;
	}
	$y2 = $y1-10;
	my $color = $green;
	print OUT1 "$CIRCLE $black $radius $color $x1 $y1\n";
	print OUT1 "$TEXT $black $smallFont $x2 $y2 $columnid\n";
	++$iline;
}
# Plot row-points
$scale[0] = $spanx[0]/260;
$scale[1] = $spanx[1]/260;
# Axis horiz & vertical
$SHIFTX = 350;
$border = 20;
$x1 = $SHIFTX-$border;
$x2 = $x1 + 300;
$y1 = $SHIFTY-$border;
$y2 = $y1 + 300;
$x0 = -$minx[0]/$scale[0]+$SHIFTX;
$y0 = -$minx[1]/$scale[1]+$SHIFTY;
print OUT1 "$LINE $black 2 $solid  $x1 $y1 $x2 $y1\n";
print OUT1 "$LINE $black 2 $solid  $x1 $y1 $x1 $y2\n";
print OUT1 "$LINE $black $thin $solid  $x0 $y1 $x0 $y2\n";
print OUT1 "$LINE $black $thin $solid  $x1 $y0 $x2 $y0\n";
$x3 = $x2-20; $y3 = $y1-8;
print OUT1 "$TEXT $black $smallFont $x3 $y3 PC1\n";
$x3 = $x1+5; $y3 = $y2+10;
print OUT1 "$TEXT $black $smallFont $x3 $y3 PC2\n";
# Print row points
$radius=2;
my $npoints = @lines_row;
my @distance;
for(my $i=$npoints-1; $i>=0; --$i){
	my ($rowid,$meanrow,$x,$y) = split(/\t/, $lines_row[$i]);
	$distance[$i] = sqrt($x*$x+$y*$y);
	$x1 = ($x-$minx[0])/$scale[0] + $SHIFTX;
	$y1 = ($y-$minx[1])/$scale[1] + $SHIFTY;
	my $aver = $average[$i];
	my $color = int(($aver-$minMean)*8/($spanMean+0.0001));
	if($color > 7){ $color=7; }
	if($color < 0){ $color=0; }
	print OUT1 "$CIRCLE $black $radius $COLOR_GRADE[$color] $x1 $y1\n";
}
close (OUT1);

# Create gif file
my $filename1 = "$fileID.gif";
system("$PATH_BIN/togif","$PATH_OUTPUT/$fileID.txt","$PATH_OUTPUT/$filename1","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");
my $image_map = "<map NAME='myMap'>\n";

my $count=0;
foreach my $i (sort {$distance[$b]<=>$distance[$a]} 0..($npoints-1)){
	my ($rowid,$meanrow,@coord) = split(/\t/, $lines_row[$i]);
	$x0 = int(($coord[0]-$minx[0])/$scale[0]+$SHIFTX);
	$y0 = int($HGT - (($coord[1]-$minx[1])/$scale[1]+$SHIFTY));
	if($x0 < 1 || $x0 > $WID-2 || $y0 < 1 || $y0 > $HGT-2) { next; }
	$x1 = $x0-3;
	$x2 = $x0+3;
	$y1 = $y0-3;
	$y2 = $y0+3;
	if($rowid =~ / \(/){
		$rowid =~ s/ \(.+$//;
	}
	$image_map .= "<area shape=rect href=\"JavaScript: plot_histogram('$rowid');\" coords=\"$x1,$y1,$x2,$y2\">\n";
	if($count++ >=2000){ last; }
}
$image_map .= "</map>\n";

my $file_matrix = $hashInput{"file_matrix"};
my $description_matrix = $hashInput{"description_matrix"};

my $vrml1 = make_vrml_file(1,\@minx,\@miny,\@spanx,\@spany,\@lines_col,\@lines_row,$minMean,$maxMean,\@average);
my $vrml2 = make_vrml_file(2,\@minx,\@miny,\@spanx,\@spany,\@lines_col,\@lines_row,$minMean,$maxMean,\@average);

# Plot graphs for each PC
$WID = 270;
my $delta;
if(@columnidList >= $MAXCOLUMNS){
	$WID = int(sqrt(@columnidList/$MAXCOLUMNS)*$WID);
	$delta = int(10*($WID-50)/(@columnidList-1))/10;
	if($delta < 8){
		$delta=8;
		$WID = $delta*(@columnidList-1)+50;
	}
}
$delta = int(10*($WID-50)/(@columnidList-1))/10;
$HGT = 120*$pc_toplot + 20 + $maxColumnHeaderLength*$charSpacing;
$radius = 2;
my $fileID1 = $fileID+1;
open (OUT1, ">$PATH_OUTPUT/$fileID1.txt");
$x_center = $WID/2;
$y_center = $HGT/2;
$nobjects = @lines_row*@lines_col*$pc_toplot + 1000;
print OUT1 "1\n$x_center\n$y_center\n$nobjects\n";
$SHIFTX = 35;
for(my $ipc=0; $ipc < $pc_toplot; ++$ipc){
	my $scale = $spany[$ipc]/110;
	if(!$scale){ $scale=0.00001; }
	$SHIFTY = ($pc_toplot-$ipc-1)*120+25+$maxColumnHeaderLength*$charSpacing;
	
	# Axis horiz & vertical
	$x1 = $SHIFTX-10;
	$x2 = $WID - 5;
	$y1 = $SHIFTY-5;
	$y2 = $y1 + 110;
	print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x2 $y1\n";
	print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x1 $y2\n";
	$y3 = $y1+60;
	my $PCnum = $ipc+1;
	print OUT1 "$TEXT $black $fontLabel 10 $y3 PC$PCnum\n";
	if($ipc == $pc_toplot-1){
		for(my $i=0; $i<@columnidList; ++$i){
			my $len = length($headers[$i]);
			my $y4 = $y1-4;
			my $y5 = $y1-$charSpacing*$len-10;
			$x1 = $i*$delta+$SHIFTX;
			print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x1 $y4\n";
			print OUT1 "$TEXT $black $fontLabel $x1 $y5 $headers[$i]\n";
		}
	}
	# Plot PC graphs
	for(my $i=0; $i<@columnidList; ++$i){
		my $coord = $pcColumns[$ipc]->[$i];
		$x1 = $i*$delta+$SHIFTX;
		$y1 = ($coord-$miny[$ipc])/$scale+$SHIFTY;
		if($i){
			print OUT1 "$LINE $black $thin $solid $x1 $y1 $x2 $y2\n";
		}
		$x2 = $x1;
		$y2 = $y1;
	}
	for(my $i=0; $i<@columnidList; ++$i){
		my $coord = $pcColumns[$ipc]->[$i];
		$x1 = $i*$delta+$SHIFTX;
		$y1 = ($coord-$miny[$ipc])/$scale+$SHIFTY;
		print OUT1 "$CIRCLE $black $radius $green $x1 $y1\n";
	}
}
close (OUT1);
my $filename2 = $fileID1.".gif";
system("$PATH_BIN/togif","$PATH_OUTPUT/$fileID1.txt","$PATH_OUTPUT/$filename2","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");

my ($filename3,$filename4);
my $file_list="";
if($pca_cluster eq "on"){
	my @file_list;
	open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Config file not found!",$logFileID);
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my @items = split(/[=\t]/,$line);
		if($items[0] =~ /^type_/){
			if(!$file_list){ $file_list = "\"".$items[1]."\""; }
			else{ $file_list .= ",\"".$items[1]."\""; }
		}
	}
	close INFO;
	# Create gif file for PC-clustes
	my $fileID2 = $fileID+2;
	open (OUT1, ">$PATH_OUTPUT/$fileID2.txt");
	$x_center = $WID/2;
	$y_center = $HGT/2;
	$nobjects = @lines_row*@lines_col*$pc_toplot + 1000;
	my $height = 110;
	print OUT1 "1\n$x_center\n$y_center\n$nobjects\n";
	for(my $ipc=0; $ipc < $pc_toplot; ++$ipc){
		$SHIFTY = ($pc_toplot-$ipc-1)*120+25+$maxColumnHeaderLength*$charSpacing;
		# Plot genes that are clustered by this PC
		plot_cluster($ipc, \%hashCluster, \%hashCoord, \@headers, $SHIFTX, $SHIFTY, "positive",$height,$delta,$WID,$pc_toplot);
	}
	close (OUT1);
	$filename3 = $fileID2.".gif";
	system("$PATH_BIN/togif","$PATH_OUTPUT/$fileID2.txt","$PATH_OUTPUT/$filename3","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");

	my $fileID3 = $fileID+3;
	open (OUT1, ">$PATH_OUTPUT/$fileID3.txt");
	$x_center = $WID/2;
	$y_center = $HGT/2;
	$nobjects = @lines_row*@lines_col*$pc_toplot + 1000;
	print OUT1 "1\n$x_center\n$y_center\n$nobjects\n";
	for(my $ipc=0; $ipc < $pc_toplot; ++$ipc){
		$SHIFTY = ($pc_toplot-$ipc-1)*120+25+$maxColumnHeaderLength*$charSpacing;
		# Plot genes that are clustered by this PC
		plot_cluster($ipc, \%hashCluster, \%hashCoord, \@headers, $SHIFTX, $SHIFTY, "negative",$height,$delta,$WID,$pc_toplot);
	}
	close (OUT1);
	$filename4 = $fileID3.".gif";
	system("$PATH_BIN/togif","$PATH_OUTPUT/$fileID3.txt","$PATH_OUTPUT/$filename4","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");

	for(my $idir=0; $idir < 2; ++$idir){
		my $direction = "positive";
		if($idir){ $direction = "negative"; }
		$image_map .= "<map NAME='mapCluster$direction'>\n";
		$x1 = $SHIFTX-5;
		$x2 = $WID-5;
		for(my $pc=1; $pc <= $pc_toplot; ++$pc){
			$SHIFTY = ($pc_toplot-$pc)*120+25+$maxColumnHeaderLength*$charSpacing;
			$y1 = $HGT - $SHIFTY;
			$y2 = $y1 - 110;
			$image_map .= "<area shape=\"rect\" href=\"JavaScript: get_cluster_list('$pc','$direction');\" coords=\"$x1,$y2,$x2,$y1\">\n";
		}
		$image_map .= "</map>\n";
	}
}
if(!$webPageID || !open(OUT,">$PATH_OUTPUT/$webPageID.txt")){
	error_message("cannot open web page ID",$logFileID);
}
print OUT "<HTML><HEAD><TITLE>ExAtlas: PCA output</TITLE>\n";
print OUT get_header();
print OUT "<SCRIPT language=JavaScript>\n";
print OUT "<!--\n";
if($pca_cluster eq "on"){
	print OUT "file_list = new Array($file_list);\n";
}
print OUT "function plot_histogram (rowname) {\n";
print OUT "	document.pca.analysis.value = \"search\";\n";
print OUT "	document.pca.search_term.value = rowname;\n";
print OUT "	var x = Math.round(Math.random()*10000);\n";
print OUT "	document.pca.target = \"_blank\"+x;\n";
print OUT "	document.pca.submit();\n";
print OUT "}\n";
print OUT "function view_vrml(file) {\n";
print OUT "	var source=\"$HOME_ADDRESS/output/\" + file;\n";
print OUT "	var x = Math.round(Math.random()*10000);\n";
print OUT "	open(source,\"_BLANK\"+x);\n";
print OUT "}\n";
if($pca_cluster eq "on"){
	print OUT "function get_cluster_list (pc, dir) {\n";
	print OUT "	document.pca.analysis.value = \"pca_cluster\";\n";
	print OUT "	document.pca.clusterPC.value = pc;\n";
	print OUT "	document.pca.clusterDir.value = dir;\n";
	print OUT "	document.pca.cluster_file.value = \"$outputPCcluster\";\n";
	print OUT "	var x = Math.round(Math.random()*10000);\n";
	print OUT "	document.pca.target = \"_blank\"+x;\n";
	print OUT "	document.pca.submit();\n";
	print OUT "}\n";
	print OUT "function save_PCA () {\n";
	print OUT "	var file = document.pca.file_geneset.value;\n";
	print OUT "	if(!file){\n";
	print OUT "		file = prompt(\"Enter file name\"); return false;\n";
	print OUT "		if(!file){ alert(\"You have to enter file name\"); return false; }\n";
	print OUT "	}\n";
	print OUT "	var file1=file;\n";
	print OUT "	if(file.search(/\\.txt\$/)>=0){\n";
	print OUT "		if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
	print OUT "		file1=file.substring(0,file.length-4);\n";
	print OUT "	}\n";
	print OUT "	if(file1.search(/^[-\\w]+\$/)<0){\n";
	print OUT "		alert(\"File name should have neither spaces nor special characters\");\n";
	print OUT "		return(false);\n";
	print OUT "	}\n";
	print OUT "	if(file.search(/^public-/i) >= 0){\n";
	print OUT "		alert(\"File name cannot start with 'public-'\"); return(false);\n";
	print OUT "	}\n";
	print OUT "	for(i=0; i<file_list.length; ++i){\n";
	print OUT "		if(file == file_list[i]){\n";
	print OUT "			alert(\"A file with this name already exists\"); return false;\n";
	print OUT "		}\n";
	print OUT "	}\n";
	print OUT "	var descrip = document.pca.description_geneset.value;\n";
	print OUT "	if(descrip.search(/\\=|\\&/) >= 0){\n";
	print OUT "		alert(\"Description should not include character \'=\' or \'&\'\");\n";
	print OUT "		return false;\n";
	print OUT "	}\n";
	print OUT "	document.pca.analysis.value = \"pca_cluster\";\n";
	print OUT "	document.pca.clusterPC.value = \"all\";\n";
	print OUT "	document.pca.clusterDir.value = \"all\";\n";
	print OUT "	document.pca.cluster_file.value = \"$outputPCcluster\";\n";
	print OUT "	document.pca.target = \"\";\n";
	print OUT "	document.pca.submit();\n";
	print OUT "}\n";
}
print OUT "<!-- end script --></SCRIPT></HEAD><p>\n";
print OUT "<table border=0><tr><td width=250 bgcolor=#DDDDDD rowspan=2><center>\n";
print OUT "To view PCA in 3D<br>you need a VRML viewer:<br>\n";
print OUT "<a href=http://freewrl.sourceforge.net>FreeWRL</a> or\n";
print OUT "<a href=http://www.cortona3d.com>Cortona3d</a><br>\n";
print OUT "<td width=600><center><b><font size=+2>Principal Component Analysis (PCA)</font><br>";
my $description = $description_matrix;
if(!$description){ $description ="None"; }
my $file_matrix1 = $file_matrix;
$file_matrix1 =~ s/^public-//;
print OUT "<font size=+1>for file $file_matrix1</b></font>";
print OUT "</table><p>\n";
print OUT "<b>Description:</b> $description &nbsp; &nbsp; <b>N columns:</b> $nCol &nbsp; &nbsp; <b>N rows:</b> $nRow<br>\n";
my $x = int(10000*rand());
print OUT "Link to download PCA output as text: <a href=\"$HOME_ADDRESS/output/$PCAfileID.txt\" target=_blank$x>PCA_file</a>\n";
if($legendID){
	my $x = int(10000*rand());
	print OUT "<br><b>Legend:</b> <a href=$HOME_ADDRESS/output/$legendID.txt target=_blank$x>Legend file</a>\n";
}
print OUT "<FORM NAME=pca ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print OUT "<p><table border=0><tr>\n";
print OUT "<td width=400><center>Experiments (tissues) &nbsp; &nbsp; ";
print OUT "<INPUT NAME=\"vrml1\" TYPE=\"button\" VALUE=\"3D View\" onClick=\"view_vrml(\'$vrml1\');\"><br>\n";
$x1 = int(10000*rand());
print OUT "Alternative <a href=$HOME_ADDRESS/output/$vrml1 target=_blank$x1>3d View</a>";
print OUT "<td><center>Genes &nbsp; &nbsp; &nbsp; ";
print OUT "<INPUT NAME=\"vrml2\" TYPE=\"button\" VALUE=\"3D View\" onClick=\"view_vrml(\'$vrml2\');\"><br>\n";
$x1 = int(10000*rand());
print OUT "Alternative <a href=$HOME_ADDRESS/output/$vrml2 target=_blank$x1>3d View</a>";
print OUT "</table>\n";
print OUT "<img src=\"$HOME_ADDRESS/output/$fileID.gif\" border=0 align=center useMap=#myMap><br>\n";
print OUT $image_map;
print OUT "This pair of graphs is a biplot, i.e. genes are more expressed in those tissues\n";
print OUT "which are located in the same area of the graph.<br>Genes are colored according to average expression (blue=low, red=high)<br>\n";
print OUT "Click on genes to view a histogram<p>\n";
if($pca_cluster eq "on"){
	print OUT "<h3>PC-based clustering of genes</h3>\n";
	print OUT "<i>Correlation threshold</i> = $cluster_correlation &nbsp; &nbsp; &nbsp; <i>Fold change threshold</i> = $cluster_fold_thresh<br>\n";
	$x1 = int(10000*rand());
	print OUT "Download tab-separated <a href=$HOME_ADDRESS/output/$outputPCcluster target=_blank$x1>table of gene clusters</a><br>\n";
	print OUT "Click on the graph with gene clusters to see the table of genes<br>\n";
	if(@columnidList < $MAXCOLUMNS){
		print OUT "<table border=0>\n";
		print OUT "<tr><td><center><b>Principal Components\n";
		print OUT "<td><center><b>Positive direction<td><center><b>Negative direction\n";
		print OUT "<tr><td><img src=\"$HOME_ADDRESS/output/$filename2\" border=0>\n";
		print OUT "<td><img src=\"$HOME_ADDRESS/output/$filename3\" border=0 useMap=#mapClusterpositive><br>\n";
		print OUT "<td><img src=\"$HOME_ADDRESS/output/$filename4\" border=0 useMap=#mapClusternegative><br>\n";
		print OUT "</table>\n";
	}else{
		print OUT "<table border=0>\n";
		print OUT "<tr><td><center><b>Principal Components\n";
		print OUT "<tr><td><img src=\"$HOME_ADDRESS/output/$filename2\" border=0>\n";
		print OUT "<tr><td><center><b>Positive direction\n";
		print OUT "<tr><td>Click on the graph with gene clusters to see the table of genes\n";
		print OUT "<tr><td><img src=\"$HOME_ADDRESS/output/$filename3\" border=0 useMap=#mapClusterpositive>\n";
		print OUT "<tr><td><center><b>Negative direction\n";
		print OUT "<tr><td>Click on the graph with gene clusters to see the table of genes\n";
		print OUT "<tr><td><img src=\"$HOME_ADDRESS/output/$filename4\" border=0 useMap=#mapClusternegative>\n";
		print OUT "</table>\n";
	}
	print OUT "<b>Legend:</b><br>Gray lines = centered gene intensity<br>Red line = average<p>\n";
	print OUT "<INPUT TYPE=button VALUE=\"Save PC clusters\" onClick=save_PCA(); style=\"width:250px;\"> Save PC clusters as geneset file</b>\n";
	print OUT "<INPUT NAME=file_geneset VALUE=$file_matrix1-PCAclusters style=\"width:250px;\"> file name<br>\n";
	print OUT "<INPUT NAME=description_geneset VALUE=$file_matrix1-PCAclusters style=width:250px;> description<p>\n";
}else{
	print OUT "<h3>Principal Components</h3>\n";
	print OUT "<img src=\"$HOME_ADDRESS/output/$filename2\" border=0><p>\n";
}
# Print the table of eigenvalues
my $cumPercent = 0;
print OUT "<h3>Table of eigenvalues</h3>\n";
print OUT "<table border=1><tr><td width=60><center>PC#<td width=60><center>Value<td width=60><center>Percent<td width=110><center>Cumulative percent\n";
for(my $i=1; $i<=@PClist; ++$i){
	if($i > 10 && $i > $pc_toplot){ last; }
	my $percent = int($PClist[$i-1]/$sumVar*100000)/1000;
	$cumPercent += $percent;
	print OUT "<tr><td><center>$i<td><center>$PClist[$i-1]<td><center>$percent<td><center>$cumPercent\n";
}
print OUT "</table>\n";
my $organismID=$hashInput{"organismID"};

print OUT "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
print OUT "<INPUT NAME=action TYPE=hidden VALUE=matrix_explore1>\n";
print OUT "<INPUT NAME=analysis TYPE=hidden VALUE=search>\n";
print OUT "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
print OUT "<INPUT NAME=description_matrix TYPE=hidden VALUE=\"$description_matrix\">\n";
if($pca_cluster eq "on"){
	print OUT "<INPUT NAME=clusterPC TYPE=hidden>\n";
	print OUT "<INPUT NAME=clusterDir TYPE=hidden>\n";
	print OUT "<INPUT NAME=cluster_file TYPE=hidden>\n";
}
print OUT "<INPUT NAME=search_term TYPE=hidden>\n";
print OUT "<INPUT NAME=category TYPE=hidden VALUE=0>\n";
print OUT "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
if($legendID){ print OUT "<INPUT NAME=legendID TYPE=hidden VALUE=$legendID>\n"; }
print OUT "</FORM><p>\n";
print OUT "<HR NOSHADE></HR>\n";
print OUT "<INPUT TYPE=button VALUE=\"    Close window   \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
close OUT;
return;
}

#******************************
sub make_vrml_file
#******************************
{
my $option = shift;
my $minx = shift;
my $miny = shift;
my $spanx = shift;
my $spany = shift;
my $lines_col = shift;
my $lines_row = shift;
my $minMean = shift;
my $maxMean = shift;
my $average = shift;

my $organismID = $hashInput{"organismID"}; 
my $file_matrix = $hashInput{"file_matrix"}; 
my $file = get_outputID();
$file .= ".wrl";
my @RGB=("0.4\t0.4\t0.8","0.6\t0.6\t1","0\t1\t1","0.6\t1\t0","0.8\t0.8\t0.2","1\t0.6\t0.2","1\t0.4\t0","1\t0\t0");
open (OUT2, ">$PATH_OUTPUT/$file");
print OUT2 "#VRML V2.0 utf8\n";
print OUT2 "NavigationInfo {speed 2 type \"EXAMINE\"}\n";
print OUT2 "DEF Camera01 Viewpoint {\n";
print OUT2 "  position 0 0 8\n";
print OUT2 "  orientation 0 1 0 0\n";
print OUT2 "  fieldOfView 0.9\n";
print OUT2 "  description \"Entry\"\n";
print OUT2 "}\n";
print OUT2 "Background {\n";
print OUT2 "  groundAngle [ 0.9, 1.5, 1.57]\n";
print OUT2 "  groundColor [ 0.77 0.8 0.82, 0.77 0.8 0.82, 0.77 0.8 0.82, 0.77 0.8 0.82]\n";
print OUT2 "  skyAngle [ 0.9, 1.5, 1.57 ]\n";
print OUT2 "  skyColor [ 0.21 0.18 0.66, 0.2 0.44 0.85, 0.51 0.81 0.95, 0.77 0.8 0.82]\n";
print OUT2 "}\n";
print OUT2 "DEF Omni01 PointLight {\n";
print OUT2 "  intensity 0.7\n";
print OUT2 "  color 1 1 1\n";
print OUT2 "  location 50 -60 -30\n";
print OUT2 "  on TRUE\n";
print OUT2 "}\n";
print OUT2 "DEF Omni02 PointLight {\n";
print OUT2 "  intensity 0.6\n";
print OUT2 "  color 1 1 1\n";
print OUT2 "  location 20 30 80\n";
print OUT2 "  on TRUE\n";
print OUT2 "}\n";
print OUT2 "Transform{translation 0 -2.5 -2.5 rotation 0 0 1 1.57 children[\n";
print OUT2 "Shape{appearance Appearance{material Material{diffuseColor 0 0 0.2}}\n";
print OUT2 "geometry Cylinder{radius 0.01 height 5}}]}\n";
print OUT2 "Transform{translation -2.5 0 -2.5 rotation 0 0 1 0 children[\n";
print OUT2 "Shape{appearance Appearance{material Material{diffuseColor 0 0 0.2}}\n";
print OUT2 "geometry Cylinder{radius 0.01 height 5}}]}\n";
print OUT2 "Transform{translation -2.5 -2.5 0 rotation 1 0 0 1.57 children[\n";
print OUT2 "Shape{appearance Appearance{material Material{diffuseColor 0 0 0.2}}\n";
print OUT2 "geometry Cylinder{radius 0.01 height 5}}]}\n";
print OUT2 "Transform{translation 2.5 -2.5 -2.5 children[Billboard{axisOfRotation 0 0 0 children[Shape{appearance Appearance{material Material{diffuseColor 0 0 0}}geometry Text{ string \"PC1\" fontStyle FontStyle{size 0.3}}}]}]}\n";
print OUT2 "Transform{translation -2.5 2.5 -2.5 children[Billboard{axisOfRotation 0 0 0 children[Shape{appearance Appearance{material Material{diffuseColor 0 0 0}}geometry Text{ string \"PC2\" fontStyle FontStyle{size 0.3}}}]}]}\n";
print OUT2 "Transform{translation -2.5 -2.5 2.5 children[Billboard{axisOfRotation 0 0 0 children[Shape{appearance Appearance{material Material{diffuseColor 0 0 0}}geometry Text{ string \"PC3\" fontStyle FontStyle{size 0.3}}}]}]}\n";
my $n3 = @$minx;
my $n5 = @$miny;
my $n4 = @$spanx;
my $n6 = @$spany;
my $npoints = @$lines_row;

if($option==2){
	my $spanMean = $maxMean-$minMean;
	for(my $i=0; $i<$npoints; ++$i){
		my ($rowid,$meanrow,@coord) = split(/\t/, $lines_row->[$i]);
		my ($x1,$y1,$z1,$color) = (0,0,0,1);
		if($spanx->[0]){
			$x1 = 5*($coord[0]-$minx->[0])/$spanx->[0]-2.5;
		}
		if($spanx->[1]){
			$y1 = 5*($coord[1]-$minx->[1])/$spanx->[1]-2.5;
		}
		if($spanx->[2]){
			$z1 = 5*($coord[2]-$minx->[2])/$spanx->[2]-2.5;
		}
		if($spanMean){
			$color = int(($average->[$i]-$minMean)*8/$spanMean);
		}
		if($color > 7){ $color=7; }
		if($color < 0){ $color=0; }
		if($rowid =~ / \(/){
			$rowid =~ s/ \(.+$//;
		}
		print OUT2 "Anchor{url \"javascript: var a = open('../bin/exatlas.cgi?category=0&search_term=$rowid&action=matrix_explore1&analysis=search&file_matrix=$file_matrix&sessionID=$sessionID&organismID=$organismID');\"\n";
		print OUT2 "children[Transform{translation $x1 $y1 $z1 children[Shape{appearance Appearance{material Material{diffuseColor $RGB[$color]}}geometry Box{size 0.05 0.05 0.05}}]}]}\n";
	}
}
foreach my $line (@$lines_col){
	my $color = "0.15 0.5 0";
	my ($columnid,@coord) = split(/\t/, $line);
	my $x1 = 5*($coord[0]-$miny->[0])/($spany->[0]+0.000001)-2.5;
	my $y1 = 5*($coord[1]-$miny->[1])/($spany->[1]+0.000001)-2.5;
	my $z1 = 5*($coord[2]-$miny->[2])/($spany->[2]+0.000001)-2.5;
	print OUT2 "Transform{translation $x1 $y1 $z1 children[Shape{appearance Appearance{material Material{diffuseColor $color}}geometry Sphere{radius 0.08}}]}\n";
	print OUT2 "Transform{translation $x1 $y1 $z1 children[Billboard{axisOfRotation 0 0 0 children[Shape{appearance Appearance{material Material{diffuseColor 0 0 0}}geometry Text{ string \"   $columnid\" fontStyle FontStyle{size 0.18}}}]}]}\n";
}
close OUT2;
return $file;
}

#***********************************
sub get_change
#***********************************
{
my $x_ref = shift;
my $y_ref = shift;

my $sx=0;
my $sy=0;
my $sxy=0;
my $sxx=0;
my $syy=0;
my $n=0;
my $minx=1000000;
my $maxx=-1000000;
for(my $i=0; $i<@$x_ref; ++$i){
	if($x_ref->[$i] <= $MISSING || $y_ref->[$i] <= $MISSING){ next; }
	$sx += $x_ref->[$i];
	$sy += $y_ref->[$i];
	$n++;
}
if($n < 2){
	print "No pattern defined\n";
	exit(0);
}
$sx /= $n;
$sy /= $n;
for(my $i=0; $i<@$x_ref; ++$i){
	if($x_ref->[$i] <= $MISSING || $y_ref->[$i] <= $MISSING){ next; }
	my $x = $x_ref->[$i]-$sx;
	my $y = $y_ref->[$i]-$sy;
	$sxy += $x*$y;
	$sxx += $x*$x;
	$syy += $y*$y;
	if($minx > $x){ $minx = $x; }
	if($maxx < $x){ $maxx = $x; }
}
my $correl = $sxy;
my $slope = $sxy;
if($sxx){ $slope /= $sxx; }
if($sxx && $syy){ $correl /= sqrt($sxx*$syy); }
$correl = floor(1000*$correl+0.5)/1000;
return ($slope*($maxx-$minx),$correl);
}

#***********************************
sub  plot_cluster
#***********************************
{
my $ipc = shift;
my $hashCluster_ref = shift;
my $hashCoord_ref = shift;
my $columnidList_ref = shift;
my $SHIFTX = shift;
my $SHIFTY = shift;
my $direction = shift;
my $height = shift;
my $delta = shift;
my $WID = shift;
my $pc_toplot = shift;

my $MAXCOLUMNS = 12;
my $fontLabel = $vertSmallFont;
my $charSpacing = 6;
if(@$columnidList_ref > $MAXCOLUMNS){
	$fontLabel = $vertTinyFont;
	$charSpacing = 5;
}
my $y1 = $SHIFTY-5;
my $y3 = $y1+60;
my $PCnum = $ipc+1;
print OUT1 "$TEXT $black $fontLabel 10 $y3 PC$PCnum\n";
if($ipc == $pc_toplot-1){
	for(my $i=0; $i<@$columnidList_ref; ++$i){
		my $len = length($columnidList_ref->[$i]);
		my $y4 = $y1-4;
		my $y5 = $y1-$charSpacing*$len-10;
		my $x1 = $i*$delta+$SHIFTX;
		print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x1 $y4\n";
		print OUT1 "$TEXT $black $fontLabel $x1 $y5 $columnidList_ref->[$i]\n";
	}
}
# Find min, max, means
my $miny = 1000000;
my $maxy = -1000000;
my @means = ();
my @nnn;
my $NNN;
foreach my $id (keys %$hashCluster_ref){
	if($hashCluster_ref->{$id}->[0] != $ipc){ next; }
	my $change = $hashCluster_ref->{$id}->[1];
	if($direction eq "positive" && $change < 0 || $direction eq "negative" && $change >= 0){ next; }
	$NNN++;
	for(my $i=0; $i<@$columnidList_ref; ++$i){
		my $y = $hashCoord_ref->{$id}->[$i];
		if($y<=$MISSING){ next; }
		if($miny > $y){ $miny = $y; }
		if($maxy < $y){ $maxy = $y; }
		$means[$i] += $y;
		$nnn[$i]++;
	}
}
# Axis horiz & vertical
my $x1 = $SHIFTX-10;
my $x2 = $WID - 5;
$y1 = $SHIFTY-5;
my $y2 = $y1 + $height;
print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x2 $y1\n";
print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x1 $y2\n";
if(!@nnn){
	my $x1 = $SHIFTX+100;
	my $y1 = $SHIFTY+50;
	print OUT1 "$TEXT $black $largeFont $x1 $y1 No genes found\n";
	return;
}
my $sumM = 0;
my $sumMM = 0;
my $n1 = @$columnidList_ref;
my $nn1=0;
for(my $i=0; $i<$n1; ++$i){
	if($nnn[$i]){
		$means[$i] /= $nnn[$i];
		$sumM += $means[$i];	
		$sumMM += $means[$i]*$means[$i];
		$nn1++;
	}else{
		$means[$i]=$MISSING;
	}	
}
if($nn1>0){
	$sumMM = ($sumMM - $sumM*$sumM/$n1)/($nn1-1);
}
my $scale = ($maxy-$miny)/$height;
# Plot PC graphs
foreach my $id (keys %$hashCluster_ref){
	if($hashCluster_ref->{$id}->[0] != $ipc){ next; }
	my $change = $hashCluster_ref->{$id}->[1];
	if($direction eq "positive" && $change < 0 || $direction eq "negative" && $change >= 0){ next; }
	my $last_missing=1;
	for(my $i=0; $i<@$columnidList_ref; ++$i){
		my $y = $hashCoord_ref->{$id}->[$i];
		if($y<=$MISSING){
			$last_missing=1;
			next;
		}
		my $x1 = $i*$delta+$SHIFTX;
		my $y1 = floor(($y-$miny)/$scale+0.5)+$SHIFTY;
		if(!$last_missing){
			print OUT1 "$LINE $gray $thin $solid $x1 $y1 $x2 $y2\n";
		}
		$last_missing=0;
		$x2 = $x1;
		$y2 = $y1;
	}
}
my $last_missing=1;
for(my $i=0; $i<@$columnidList_ref; ++$i){
	my $y = $means[$i];
	if($y<=$MISSING){
		$last_missing=1;
		next;
	}
	my $x1 = $i*$delta+$SHIFTX;
	my $y1 = ($y-$miny)/$scale+$SHIFTY;
	if(!$last_missing){
		print OUT1 "$LINE $red $thin $solid $x1 $y1 $x2 $y2\n";
	}
	$last_missing=0;
	$x2 = $x1;
	$y2 = $y1;
}
$y3 = $SHIFTY+$height-5;
my $x3 = $SHIFTX+15;
print OUT1 "$TEXT $black $smallFont $x3 $y3 N=$NNN\n";
return;
}

#**************************************
sub   check_matrix_headers
#**************************************
{
my $file_matrix = shift;
my $organismID = shift;
my $errorMessage = shift;

my %hashMatrix;
#print "Parsing matrix<br>\n";
parse_matrix_file($file_matrix, \%hashMatrix);
my $ref = $hashMatrix{"sample_title"};
if(!$ref){ $$errorMessage .= "<b>Error:</b> Matrix file has no sample titles!"; }
my $taxid = $hashMatrix{"series_sample_taxid"};
if(!$taxid || $taxid=~ /[,;]/){
	my $ref = $hashMatrix{"sample_taxid_ch1"};
	if($ref){
		my $taxid1 = $ref->[0];
		if($taxid1){ $taxid=$taxid1; }
	}
	$ref = $hashMatrix{"sample_taxid_ch2"};
	if(!$taxid && $ref){
		my $taxid1 = $ref->[0];
		if($taxid1){ $taxid=$taxid1; }
	}
}
if(!$taxid){ $$errorMessage .= "<b>Error:</b> Matrix file has no organismID!"; } 
elsif($taxid != $organismID){ $$errorMessage .= "<b>Error:</b> Wrong organism ID!"; } 
if($$errorMessage){ return 0; }
my @headers = @$ref;
my $nCol = @headers;
#print "nCol $nCol<br>\n";
my $changed = 0;
for(my $i=0; $i<$nCol; ++$i){
	if($headers[$i] =~ /(replicate|replication|rep)[\s_]*\d+$/i){
		$changed = 1;
		$headers[$i] =~ s/,*[\s_]+(biological[\s_]+|)(replicate|replication|rep)[\s_]*\d+$//i;
	}
}
if(!$changed){ return 1; }
my $fileID = get_outputID(1);
open(INFO,'<',$file_matrix);
open(OUT, ">$PATH_OUTPUT/$fileID.txt");
while(my $line=<INFO>){
	if($line =~ "!Sample_title"){
		print OUT "!Sample_title\t\"".join("\"\t\"",@headers)."\"\n";
	}else{
		print OUT $line;
	}
}
close OUT;
close INFO;
copy "$PATH_OUTPUT/$fileID.txt", "$file_matrix";
return 1;
}

#**************************************
sub   matrix_quality
#**************************************
{
my $file_matrix = $hashInput{"file_matrix"};
my $organismID = $hashInput{"organismID"};

my %hashMatrix;
my $file_matrix_full = $file_matrix;
if($file_matrix !~ /^public-/){
	$file_matrix_full = "$loginname-$file_matrix";
}
parse_matrix_file("$PATH_DATA/$file_matrix_full", \%hashMatrix);
my $organismID1 = $hashMatrix{"series_sample_taxid"};
if($organismID1 && $organismID != $organismID1){ error_message("Organism ID do not match"); }
my $ref = $hashMatrix{"sample_title"};
if(!$ref){ error_message("Matrix file has no sample titles!"); }
my @headers = @$ref;
my $nCol = @headers;
my $logFileID = get_outputID(7); #0-log,1-web,2-output,3-4 plot1,5-6 plot2
$hashInput{"logFileID"} = $logFileID;
$hashInput{"runID"} = $RUN_QUALITY;
interrupt_program($nCol/100);
exit(0);
}

#**************************************
sub   matrix_quality1
#**************************************
{
my $logFileID = shift;
my $webPageID = $hashInput{"logFileID"}+1;
my $outputID = $webPageID+1;

my $file_matrix = $hashInput{"file_matrix"};
my %hashMatrix;
my $file_matrix_full = $file_matrix;
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_matrix_full = "$loginname-$file_matrix";
	$file_anova = "$loginname-anova-$file_matrix";
}
parse_matrix_file("$PATH_DATA/$file_matrix_full", \%hashMatrix);
my $organismID = $hashMatrix{"series_sample_taxid"};
my $ref = $hashMatrix{"sample_title"};
my @headers = @$ref;
my $nCol = @headers;
my @replications;
my @used;
my $k=-1;
for(my $i=0; $i<$nCol; ++$i){
	if($used[$i]){ next; }
	$k++;
	push(@{$replications[$k]}, $i);
	for(my $j=$i+1; $j<$nCol; ++$j){
		if($used[$j]){ next; }
		if($headers[$j] eq $headers[$i]){
			push(@{$replications[$k]},$j);
			$used[$j] = 1;
		}
	}
}
$ref = $hashMatrix{"sample_geo_accession"};
my @sampleID;
if($ref){
	@sampleID = @$ref;
}else{
	@sampleID = @headers;
}
my %expression;
my $log10 = log(10);
if(open(INFO,"<$PATH_DATA/$organismID"."_expression.txt")){
	my $line = <INFO>;
	my $count=0;
	while($line = <INFO>){
		chop $line;
		my ($symbol,$expr,$stable)=split(/\t/,$line);
		if($stable>0){
			$expression{$symbol} = $expr;
			$count++;
		}
	}
	close INFO;
	if($count < 10){ %expression=(); }
}

my %bestProbe;
my %annotation;
open (INFO,"<$PATH_DATA/$file_anova") or error_message("No anova file!",$logFileID);
my $line = <INFO>;
chop $line;
my ($junk,$junk1,@headers1) =split(/\t/,$line);
if(!@headers1){ error_message("Headers not found in $file_anova",$logFileID); }
my $nCol1=0;
while($nCol1<@headers1 && $headers1[$nCol1] !~ /^Var\(/){ ++$nCol1; }
my $nRep1=0;
for(my $i=0; $i<$nCol1; ++$i){
	if($headers1[$i]=~/ \(\d+\)$/ && $headers1[$i] !~ /^Mean/){
		my @items = split(/ \(/,$headers1[$i]);
		my $n = pop(@items);
		my $name = join(" (",@items);
		$n =~ s/\)$//;
		$nRep1 += $n;
	}else{
		$nRep1++;
	}
}
while(my $line = <INFO>){
	chop $line;
	my ($id,$expr,@data1) =split(/\t/, $line);
	my $MSE = $data1[$nCol1+3];
	my $symbol = $data1[$nCol1+8];
	my $F = $data1[$nCol1+4];
	if($nRep1<=$nCol1){ $F=$expr; }
	if($expression{$symbol}){
		my $ref = $bestProbe{$symbol};
		if(!$ref || $ref->[1] < $F){
			$bestProbe{$symbol} = [$id,$F];
		}
	}
	$annotation{$id} = [$symbol,$MSE];
}
close INFO;
if($logFileID){ file_append("Correlation analysis started ..","$PATH_OUTPUT/$logFileID.txt"); }
my $count;
open(INFO,"<$PATH_DATA/$file_matrix_full") or error_message("Cannot open $file_matrix",$logFileID);
while(my $line=<INFO>){
	if($line =~ /^!series_matrix_table_begin/i){ last; }
}
my @missing_values;
my @Data;
my @stable_genes;
my @average;
my @MSE;
my $xmin = 0.001*$hashMatrix{"series_xmin"};
$line=<INFO>;
while($line=<INFO>){
	if($line =~ /^!/){ last; }
	chop $line;
	$line =~ s/\s+$//;
	my($probe,@data1) = split(/\t/,$line);
	$probe =~ s/^\"//;
	$probe =~ s/\"$//;
	my $symbol = $annotation{$probe}->[0];
	my $stable=0;
	my $expr = 0;
	if(%expression && $symbol && $bestProbe{$symbol} && $bestProbe{$symbol}->[0] eq $probe){
		$stable=1;
		push(@{$stable_genes[0]},log($expression{$symbol})/$log10);
	}
	my $i1 = 0;
	for(my $ir=0; $ir<@replications; ++$ir){
		foreach my $i (@{$replications[$ir]}){
			my $x = $data1[$i];
			if($x<=$MISSING){ $missing_values[$i1]++; }
			else{
				if($x < $xmin){ $x=$xmin; }
				$x = log($x)/$log10;
				$data1[$i] = $x;
				if(%expression && $stable){
					push(@{$stable_genes[$i1+1]},$x);
				}
			}
			$i1++;
		}
	}
	push(@average,average(\@data1));
	push(@Data,\@data1);
	push(@MSE,$annotation{$probe}->[1]);
}
close INFO;
my @correlation;
for(my $i=0; $i<$nCol; ++$i){
	if(%expression){
		my ($r,$n) = pearson_correlation($stable_genes[0],$stable_genes[$i+1]);
		$correlation[$i] = floor(10000*$r+0.5)/10000;
	}else{
		$correlation[$i] = 0;
	}
	if(!$missing_values[$i]){ $missing_values[$i]=0; }
}
if($logFileID){
	file_append("Correlation analysis done","$PATH_OUTPUT/$logFileID.txt");
	file_append("Analysis of replications started ..","$PATH_OUTPUT/$logFileID.txt");
}
my @sorted = sort {$a<=>$b} @average;
my @thresh;
my $nThresh = 10;
my $nn = @sorted;
for(my $i=0; $i<$nThresh; $i++){
	$thresh[$i] = $sorted[int($nn*($i+0.5)/($nThresh+1)+0.5)];
}
while($thresh[0]==$thresh[1] || $thresh[0]<=$MISSING){
	shift(@thresh);
	$nThresh--;
}
my @ss;
my @ss1;
my @nn;
my @nn1;
for(my $i=0; $i<@Data; ++$i){
	my $ref = $Data[$i];
	my $avr = $average[$i];
	my $SD = sqrt($MSE[$i]);
	my $t = 0;
	while($t<$nThresh && $avr > $thresh[$t]){ ++$t; }
	my $nn=0;
	for(my $j=0; $j<@$ref; ++$j){
		my $x = $ref->[$j];
		if($x>$MISSING){ $nn++; }
	}
	my $wgt0 = 1-1/$nn;
	my $i1=0;
	for(my $ir=0; $ir<@replications; ++$ir){
		my @array;
		my $i2 = 0;
		foreach my $k (@{$replications[$ir]}){
			my $x = $ref->[$k];
			if($x>$MISSING){
				$ss[$t]->[$i1+$i2] += ($x-$avr)*($x-$avr);
				$nn[$t]->[$i1+$i2] += $wgt0;
				push(@array,$x);
			}
			$i2++;
		}
		if(@array > 1){
			my $med = median(\@array);
			my $nn = @array;
			my $wgt = 1-(1/$nn);
			my @outliers;
			for(my $k=0; $k<@array; $k++){
				my $x = $array[$k];
				my $d = $x-$med;
				my $st = 1;
				if($SD>0){ $st = abs($x-$med)/$SD; }
				if($st>2.5 && abs($x-$avr) > abs($med-$avr)){
					push(@outliers,splice(@array,$k,1));
					$k--;
				}
			}
			if(@outliers){
				$med = median(\@array);
			}
			my $i3 = $i1;
			foreach my $k (@{$replications[$ir]}){
				my $x = $ref->[$k];
				if($x>$MISSING){
					my $d = $x-$med;
					$ss1[$t]->[$i3] += $d*$d;
					$nn1[$t]->[$i3] += $wgt;
				}
				$i3++;
			}
		}
		$i1 += $i2;
	}
}
my @midThresh;
for(my $t=0; $t<=$nThresh; ++$t){
	if($t>0 && $t<$nThresh){
		$midThresh[$t] = int(100*0.5*($thresh[$t-1]+$thresh[$t]))/100;
	}elsif($t==0){
		$midThresh[$t] = int(100*($thresh[0] - 0.5*($thresh[1]-$thresh[0])))/100;
	}else{
		$midThresh[$t] = int(100*($thresh[$nThresh-1] + 0.5*($thresh[$nThresh-1]-$thresh[$nThresh-2])))/100;
	}
	for(my $i=0; $i<$nCol; ++$i){
		if($nn[$t]->[$i] > 0){
			$ss[$t]->[$i] = int(10000*sqrt($ss[$t]->[$i]/$nn[$t]->[$i])+0.5)/10000;
		}else{
			$ss[$t]->[$i] = 0;
		}
		if($nn1[$t]->[$i] > 0){
			$ss1[$t]->[$i] = int(10000*sqrt($ss1[$t]->[$i]/$nn1[$t]->[$i])+0.5)/10000;
		}else{
			$ss1[$t]->[$i] = 0;
		}
	}
}
if($logFileID){
	file_append("Analysis of replications done","$PATH_OUTPUT/$logFileID.txt");
	file_append("Printing output results ..","$PATH_OUTPUT/$logFileID.txt");
}
if(!$outputID || !open(OUT,">$PATH_OUTPUT/$outputID.txt")){
	error_message("cannot open output file ID",$logFileID);
}
print OUT "Quality report for matrix $file_matrix\n\n";
print OUT "Col.#\tSample ID\tSeries ID\tPlatform ID\tN missing\tCorrelation\tS.D.(glob)\tS.D.(replic)\t Sample title\n";
my $ref1 = $hashMatrix{"sample_geo_accession"};
my $ref2 = $hashMatrix{"sample_platform_id"};
my $ref3 = $hashMatrix{"sample_series_geo_accession"};
my $seriesID = $hashMatrix{"series_geo_accession"};
my $text;
my @sampleID1;
my $icol = 0;
my $platform = $hashMatrix{"series_platform_id"};
for(my $ir=0; $ir<@replications; ++$ir){
	foreach my $i (@{$replications[$ir]}){
		my $sampleID = $headers[$i];
		if($ref1 && $ref1->[$i]){ $sampleID = $ref1->[$i]; }
		push(@sampleID1,$sampleID);
		my $platformID = $platform;
		if($ref2 && $i<@$ref2){ $platformID = $ref2->[$i]; }
		if($ref3 && $i<@$ref3){ $seriesID = $ref3->[$i]; }
		if(!$seriesID){ $seriesID="None"; }
		my $col = $icol+1;
		my $line = "$col\t$sampleID\t$seriesID\t$platformID\t$missing_values[$icol]\t$correlation[$icol]\t$ss[$nThresh-1]->[$icol]\t$ss1[$nThresh-1]->[$icol]\t$headers[$i]\n";
		print OUT $line;
		$text .= "<TR><TD>";
		if($file_matrix !~ /public-/ || $loginname eq "public"){
			$text .= "<INPUT TYPE=checkbox NAME=column$i> ";
		}
		$text .= $line;
		$icol++;
	}
}
my $nClass = $nThresh+1;
print OUT "\nTable of standard deviation from global mean for each set of genes grouped by average expression\n";
my $text1 = "\nTable of standard deviation for replications for each set of genes grouped by average expression\n";
$line = "Col.#\tSample ID";
for(my $t=0; $t<$nClass; ++$t){
	$line .= "\tExpr=$midThresh[$t]";
}
$line .= "\n";
print OUT $line;
$text1 .= $line;
for(my $icol=0; $icol<$nCol; ++$icol){
	my $col = $icol+1;
	print OUT "$col\t$sampleID1[$icol]";
	$text1 .= "$col\t$sampleID1[$icol]";
	for(my $t=0; $t<$nClass; ++$t){
		print OUT "\t$ss[$t]->[$icol]";
		$text1 .= "\t$ss1[$t]->[$icol]";
	}
	print OUT "\n";
	$text1 .= "\n";
}
print OUT $text1;
close OUT;
my $plotID1 = $outputID+1;
$hashInput{"sort_histogram"}="";
if(%expression){
	plot_histogram_horiz("bars",$nCol,1,\@correlation,\@sampleID1,$plotID1,"Correlation",0);
}
my $plotID2 = $outputID+3;
if(@replications<$nCol){
	plot_histogram_horiz("bars",$nCol,1,$ss1[$nThresh-1],\@sampleID1,$plotID2,"S.D. of replicates, log10","zero");
}
# Print page header
if(!$webPageID || !open(OUT,">$PATH_OUTPUT/$webPageID.txt")){
	error_message("cannot open web page ID",$logFileID);
}
print OUT "<HTML><HEAD><TITLE>ExAtlas - matrix quality</TITLE>\n";
print OUT get_header();
print OUT "<script type=text/JavaScript language=JavaScript>\n";
print OUT "<!-- \n";
print OUT "function count_checked_samples(){\n";
print OUT "	var nchecked = 0;\n";
print OUT "	for(i=0; i<document.quality_control.elements.length; ++i){\n";
print OUT "		if(document.quality_control.elements[i].type==\"checkbox\"){\n";
print OUT "			if(document.quality_control.elements[i].checked){ nchecked++; };\n";
print OUT "		}\n";
print OUT "	}\n";
print OUT "	return nchecked;\n";
print OUT "}\n";
print OUT "function get_warning() {\n";
print OUT "	var nchecked = count_checked_samples();\n";
print OUT "	if(!nchecked){\n";
print OUT "		alert(\"You have not selected any samples. Operation cancelled.\");\n";
print OUT "		return false;\n";
print OUT "	}\n";
print OUT "	if(!confirm(\"You have selected to delete \"+nchecked+\" samples.\\nIf you proceed, you will need to close this data set and opet it again to redo ANOVA.\\nDo you want to delete these samples?\")){\n";
print OUT "		return false;\n";
print OUT "	}\n";
print OUT "	document.quality_control.remove_samples.value=1;\n";
print OUT "	document.quality_control.submit();\n";
print OUT "}\n";
print OUT "<!-- end script --></SCRIPT></HEAD>\n";
print OUT "<p><font size=+2><b>Quality of matrix \"$file_matrix\"</b></font>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;";
print OUT "<INPUT TYPE=button VALUE=\" Cancel (close window) \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
my $x = int(10000*rand());
print OUT "Text file with full information on quality: <a href=$HOME_ADDRESS/output/$outputID.txt target=_BLANK$x>TABLE</a>.<p>\n";
print OUT "<TABLE BORDER=0>\n";
if(%expression){
	print OUT "<TR><TD><center><b>Correlation of gene expression for housekeeping genes</b>\n";
	print OUT "<TR><TD><IMG SRC=../output/$plotID1.gif BORDER=0></A><p>\n";
}
if(@replications<$nCol){
	print OUT "<TR><TD><center><b>Standard deviation (SD) for replications</b>\n";
	print OUT "<TR><TD><IMG SRC=../output/$plotID2.gif BORDER=0></A>\n";
}
print OUT "</TABLE><p>\n";
print OUT "<FORM NAME=quality_control ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
if($file_matrix !~ /public-/ || $loginname eq "public"){
	print OUT "<INPUT TYPE=button VALUE=\"Remove checked probes\" style=width:200px; onClick=\"get_warning();\"><p>\n";
}
print OUT "<TABLE BORDER=0>\n";
print OUT "<TR><TD>Col.#<TD>Sample ID<TD>Series ID<TD>Platform<TD>N missing<TD>Correlation<TD>S.D.(glob)<TD>S.D.(replic)<TD>Sample title\n";
$text =~ s/\t /<TD>/g;
$text =~ s/\t/<TD><center>/g;
print OUT $text;
print OUT "</TABLE>\n";
print OUT "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
print OUT "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
print OUT "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print OUT "<INPUT NAME=remove_samples TYPE=hidden>\n";
print OUT "</FORM><HR NOSHADE></HR>\n";
print OUT "<INPUT TYPE=button VALUE=\" Cancel (close window) \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print OUT "<p><i>Please report any problems to <a href=mailto:sharoval\@mail.nih.gov>webmaster</a><p>\n";
print OUT "</BODY>\n";
print OUT "</HTML>\n";
close OUT;
return;
}

#**************************************
sub   remove_samples
#**************************************
{
my $file_matrix = $hashInput{"file_matrix"};
my %hashMatrix;
my $file_matrix_full = $file_matrix;
my $file_anova = $file_matrix;
$file_anova =~ s/^public-/public-anova-/;
if($file_matrix !~ /^public-/){
	$file_matrix_full = "$loginname-$file_matrix";
	$file_anova = "$loginname-anova-$file_matrix";
}
parse_matrix_file("$PATH_DATA/$file_matrix_full", \%hashMatrix);
my $ref = $hashMatrix{"sample_title"};
if(!$ref || ref($ref) ne 'ARRAY'){ error_message("No column headres"); }
my $n = @$ref;
my @remove;
for(my $i=0; $i<$n; $i++){
	if($hashInput{"column$i"} eq "on"){ push(@remove,$i); }
}
if(!@remove){ error_message("No samples were selected for removal. Task cancelled"); }
@remove = sort {$b<=>$a} @remove;
my $fileID = get_outputID(1);
open(INFO,"<$PATH_DATA/$file_matrix_full") or error_message("Cannot open $file_matrix");
open(OUT,">$PATH_OUTPUT/$fileID.txt") or error_message("Cannot open temp file");
my $line;
my $done=0;
while($line=<INFO>){
	if($line =~ /^!sample/i){
		chop $line;
		$line =~ s/\s+$//;
		my ($key,@data) = split(/\t/,$line);
		if(@data>1){
			foreach my $k (@remove){ splice(@data,$k,1); }
			if(!@data){
				error_message("You cannot select all samples for deletion!");
			}
			$line = join("\t",$key,@data)."\n";
		}
	}
	print OUT $line;
	if($line =~ /^!series_matrix_table_begin/i){ last; }
}
while($line=<INFO>){
	chop $line;
	$line =~ s/\s+$//;
	my($probe,@data) = split(/\t/,$line);
	if(@data>1){
		foreach my $k (@remove){ splice(@data,$k,1); }
		$probe =~ s/^\"//;
		$probe =~ s/\"$//;
		$line = join("\t",$probe,@data)."\n";
	}
	print OUT $line;
}
close INFO;
close OUT;
copy "$PATH_OUTPUT/$fileID.txt", "$PATH_DATA/$file_matrix_full";
unlink "$PATH_OUTPUT/$fileID.txt";
unlink "$PATH_DATA/$file_anova";

print "<HTML><HEAD><TITLE>ExAtlas samples removed</TITLE>\n";
print_header();
print "<H3>Selected samples were removed from file '$file_matrix'</H3>\n";
print "When you open this file next time, ANOVA will run again to reflect the changes<p>\n";
print "<HR NOSHADE></HR>\n";
print "<INPUT TYPE=button VALUE=\" Close window \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub   nonredundant_matrix
#**************************************
{
my $file_matrix = $hashInput{"file_matrix"};

remove_redundant_genes($file_matrix);
print "<HTML><HEAD><TITLE>ExAtlas</TITLE>\n";
print_header();
print "<H3>Best probe is selected for each gene in '$file_matrix'</H3>\n";
print "Method: probes with highest F-statistics are selected for each gene. Non-annotated probes removed<p>\n";
print "When you open this file next time, ANOVA will run again to reflect the changes<p>\n";
print "<HR NOSHADE></HR>\n";
print "<INPUT TYPE=button VALUE=\" Close window \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub   remove_redundant_genes
#**************************************
{
my $file_matrix = shift;

my $file_matrix_full = "$loginname-$file_matrix";
my $file_anova = "$loginname-anova-$file_matrix";
my @sorted;
my %symbols;
open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file");
my $line = <INFO>;
chop $line;
my ($junk,$junk1,@headers) =split(/\t/,$line);
if(!@headers){ error_message("Headers not found in $file_anova"); }
my $nCol=0;
while($nCol<@headers && $headers[$nCol] !~ /^Var\(/){ ++$nCol; }
my $nRep=0;
for(my $i=0; $i<$nCol; ++$i){
	if($headers[$i]=~/ \(\d+\)$/ && $headers[$i] !~ /^Mean/){
		my @items = split(/ \(/,$headers[$i]);
		my $n = pop(@items);
		my $name = join(" (",@items);
		$n =~ s/\)$//;
		$nRep += $n;
	}else{
		$nRep++;
	}
}
my $count;
while($line=<INFO>){
	chop $line;
	my ($probe_id,$expr,@data) = split(/\t/,$line);
	my $F = $data[$nCol+4];
	#if($nRep<=$nCol){ $F=$expr; }
	my $symbol = $data[$nCol+8];
	my $ref = $symbols{$symbol};
	if(!$symbol){
	}elsif(!$ref){
		$symbols{$symbol} = [$probe_id,$F,$count];
	}elsif($ref->[1] < $F){
		$ref->[0] = $probe_id;
		$ref->[1] = $F;
		$ref->[2] = $count;
	}
	$count++;
}
my $nRow = keys %symbols;
@sorted=();
foreach my $symbol (keys %symbols){
	push(@sorted,$symbols{$symbol}->[2]);
}
@sorted = sort {$a<=>$b} @sorted;

my $fileID = get_outputID();
open(INFO,"<$PATH_DATA/$file_matrix_full") or error_message("Cannot open $file_matrix");
open(OUT,">$PATH_OUTPUT/$fileID.txt") or error_message("Cannot open temp.txt");
my $done=0;
while(my $line=<INFO>){
	if($line =~ /^!Sample_data_row_count\t/){
		my ($key,@items)=split(/\t/,$line);
		print OUT "!Sample_data_row_count";
		for(my $i=0; $i<@items; ++$i){
			print OUT "\t\"$nRow\"";
		}
		print OUT "\n";
		next;
	}
	if($line =~ /^!Series_nonredundant/i){ $done=1; }
	if(!$done && ($line !~ /^!series/i || $line =~ /^!series_matrix_table_begin/i)){
		print OUT "!Series_nonredundant\t\"true\"\n";
		$done=1;
	}
	print OUT $line;
	if($line =~ /^!series_matrix_table_begin/i){ last; }
}
$line = <INFO>;
print OUT $line;
$count = 0;
my $ii=0;
while($line=<INFO>){
	if($line =~ /^!/){ last; }
	if($count==$sorted[$ii]){
		print OUT $line;
		$ii++;
	}
	++$count;
}
close INFO;
if($line =~ /^!series_matrix_table_end/i){
	print OUT $line;
}else{
	print OUT "!series_matrix_table_end\n";
}
close OUT;
copy "$PATH_OUTPUT/$fileID.txt", "$PATH_DATA/$file_matrix_full";
unlink "$PATH_OUTPUT/$fileID.txt";
unlink "$PATH_DATA/$file_anova";
return;
}

#**************************************
sub  cluster_matrix
#**************************************
{
my $fileID = shift;
my $nCol = shift;
my $nRow = shift;
my $cluster_type = shift;   # 1=rows 2=cols, 3=both
my $diagonal = shift;
my $logFileID = shift;

if($cluster_type==0){ return; }
my @columnReorder;
my @rowReorder;
my $PCAoutput = 0;
my $npc=0;

my $PCA_lines=0;	
if($nCol==1){
	if($cluster_type==2){ return; }
	my @data;
	open (INFO,"<$PATH_OUTPUT/$fileID.txt");
	my $line = <INFO>;
	my $count=0;
	while(my $line = <INFO>){
		chop $line;
		my ($row_header,$data1)=split(/\t/, $line);
		push(@data,$data1);
	}
	close INFO;
	foreach my $x (sort {$data[$b]<=>$data[$a]} 0..($nRow-1)){
		push(@rowReorder,$x+1);
	}
}
elsif($nCol==2){
	if($cluster_type==2){ return; }
	my @data;
	open (INFO,"<$PATH_OUTPUT/$fileID.txt");
	my $line = <INFO>;
	while(my $line = <INFO>){
		chop $line;
		my ($row_header,@data1)=split(/\t/, $line);
		if($data1[0]==$MISSING){ $data1[0]=0; }
		if($data1[1]==$MISSING){ $data1[1]=0; }
		my $x = $data1[0]+$data1[1];
		my $y = $data1[0]*$data1[1];
		my $z = abs($data1[0])-abs($data1[1]);
		push(@data1,$x,$y,$z);
		push(@data,\@data1);
	}
	close INFO;
	@rowReorder = sort {$data[$b]->[3]<=>$data[$a]->[3]} 0..($nRow-1);
	my $threshold = $data[$rowReorder[0]]->[3]*0.25;
	my $ii=0;
	while($ii<@rowReorder && $data[$rowReorder[$ii]]->[3] > $threshold){ $ii++; }
	my @residual = splice(@rowReorder,$ii);
	@rowReorder = sort {$data[$b]->[2]<=>$data[$a]->[2]} @rowReorder;
	$ii=0;
	while($ii<@residual && $data[$residual[$ii]]->[3] > -$threshold){ $ii++; }
	my @residual1 = splice(@residual,0,$ii);
	@residual = sort {$data[$b]->[0]<=>$data[$a]->[0]} @residual;
	push(@rowReorder,@residual);   #Add opposite
	@residual1 = sort {$data[$b]->[4]<=>$data[$a]->[4]} @residual1;
	$ii=0;
	while($ii<@residual && $data[$residual[$ii]]->[4] >0){ $ii++; }
	@residual = splice(@residual1,0,$ii);
	@residual = sort {$data[$b]->[0]<=>$data[$a]->[0]} @residual;
	push(@rowReorder,@residual);
	@residual1 = sort {$data[$b]->[1]<=>$data[$a]->[1]} @residual1;
	push(@rowReorder,@residual1);
	for(my $i=0; $i<@rowReorder; ++$i){ $rowReorder[$i]++; }
}
elsif($nRow==1){
	if($cluster_type==1){ return; }
	my @data;
	open (INFO,"<$PATH_OUTPUT/$fileID.txt");
	my $line = <INFO>;
	my $count=0;
	while(my $line = <INFO>){
		chop $line;
		my ($row_header,@data1)=split(/\t/, $line);
		@data = @data1;
	}
	close INFO;
	foreach my $x (sort {$data[$b]<=>$data[$a]} 0..($nCol-1)){
		push(@columnReorder,$x+1);
	}
}
elsif($nRow==2){
	if($cluster_type==1){ return; }
	my @data1;
	my @data;
	open (INFO,"<$PATH_OUTPUT/$fileID.txt");
	my $line = <INFO>;
	while(my $line = <INFO>){
		chop $line;
		$line =~ s/$MISSING/0/g;
		my ($row_header,@data2)=split(/\t/, $line);
		push(@data1,\@data2);
	}
	close INFO;
	for(my $i=0; $i<$nCol; $i++){
		my $x1 = $data1[0]->[$i];
		my $x2 = $data1[1]->[$i];
		my $x = $x1+$x2;
		my $y = $x1*$x2;
		my $z = abs($x1)-abs($x2);
		push(@data,[$x1,$x2,$x,$y,$z]);
	}
	close INFO;
	@columnReorder = sort {$data[$b]->[3]<=>$data[$a]->[3]} 0..($nCol-1);
	my $threshold = $data[$columnReorder[0]]->[3]*0.25;
	my $ii=0;
	while($ii<@columnReorder && $data[$columnReorder[$ii]]->[3] > $threshold){ $ii++; }
	my @residual = splice(@columnReorder,$ii);
	@columnReorder = sort {$data[$b]->[2]<=>$data[$a]->[2]} @columnReorder;
	$ii=0;
	while($ii<@residual && $data[$residual[$ii]]->[3] > -$threshold){ $ii++; }
	my @residual1 = splice(@residual,0,$ii);
	@residual = sort {$data[$b]->[0]<=>$data[$a]->[0]} @residual;
	push(@columnReorder,@residual);   #Add opposite
	@residual1 = sort {$data[$b]->[4]<=>$data[$a]->[4]} @residual1;
	$ii=0;
	while($ii<@residual && $data[$residual[$ii]]->[4] >0){ $ii++; }
	@residual = splice(@residual1,0,$ii);
	@residual = sort {$data[$b]->[0]<=>$data[$a]->[0]} @residual;
	push(@columnReorder,@residual);
	@residual1 = sort {$data[$b]->[1]<=>$data[$a]->[1]} @residual1;
	push(@columnReorder,@residual1);
	for(my $i=0; $i<@columnReorder; ++$i){ $columnReorder[$i]++; }
}
elsif(!$diagonal){
	$PCAoutput = get_outputID(1);
	$npc = int(2*sqrt($nCol));
	if($npc > $nCol){ $npc=$nCol; }

	my $return = system("$PATH_BIN/pca","$PATH_OUTPUT/$fileID.txt","$PATH_OUTPUT/$PCAoutput.txt","$npc","0","V","1");
	if($return){ error_message("PCA crashed in heatmap", $logFileID); }
	my @eigenval;
	open (INFO,"<$PATH_OUTPUT/$PCAoutput.txt") or error_message("In cluster_matrix $PCAoutput.txt", $logFileID);
	while(my $line = <INFO>){ $PCA_lines++; }
	close INFO;
	if($PCA_lines>10){
		open (INFO,"<$PATH_OUTPUT/$PCAoutput.txt");
		if(($cluster_type==1 || $cluster_type==3) && $nRow > 2000){
			while(my $line = <INFO>){
				if($line=~ /^Eigenvalues/){ last; }
			}
			while(my $line = <INFO>){
				chop $line;
				if(!$line){ last; }
				if($line >0){
					push(@eigenval,sqrt($line));
				}
			}
			my @data;
			while(my $line = <INFO>){
				if($line =~ /^Row means/){ last; }
			}
			my $count=0;
			while(my $line = <INFO>){
				chop $line;
				if(!$line){ last; }
				my($header,$mean,@pc) = split(/\t/,$line);
				for(my $ipc=0; $ipc<$npc; ++$ipc){
					$pc[$ipc] = int(10000*$pc[$ipc]*$eigenval[$ipc]+0.5)/10000;
				}
				push(@data,[$count++,$header,\@pc]);
			}
			@rowReorder = pca_based_clustering(\@data);
		}
		if($cluster_type==2 || $cluster_type==3){
			my $line;
			while($line = <INFO>){
				if($line =~ /^Sorted columns/){ last; }
			}
			if($line !~ /^Sorted columns/){ error_message("Sorted columns not found", $logFileID); }
			$line = <INFO>;
			chop $line;
			@columnReorder = split(/,/,$line);
		}
		if(($cluster_type==1 || $cluster_type==3) && $nRow <= 2000){
			my $line;
			while($line = <INFO>){
				if($line =~ /^Sorted rows/){ last; }
			}
			if($line !~ /^Sorted rows/){ error_message("Sorted rows not found", $logFileID); }
			$line = <INFO>;
			chop $line;
			@rowReorder = split(/,/,$line);
		}
		close INFO;
	}
}
if($diagonal==1 || $nCol>2 && $nRow>2 && $PCA_lines<=1){
	my @data;
	open (INFO,"<$PATH_OUTPUT/$fileID.txt");
	my $line = <INFO>;
	while(my $line = <INFO>){
		chop $line;
		my ($row_header,@data1)=split(/\t/, $line);
		push(@data,\@data1);
	}
	close INFO;
	my @M;
	my @usedRow;
	my @usedCol;
	for(my $ir=0; $ir<$nRow; $ir++){
		$usedRow[$ir]=0;
	}
	for(my $ic=0; $ic<$nCol; $ic++){
		$usedCol[$ic]=0;
		for(my $ir=0; $ir<$nRow; $ir++){
			my $x = $data[$ir]->[$ic];
			if($x == $MISSING){ $x=0; }
			$M[$ir]->[$ic] = $x;
		}
	}
	foreach my $pair (sort_diagonal(\@M)){
		if(!$pair){ next; }
		my ($ic,$ir) = @$pair;
		#print "$ic $ir $M[$ir]->[$ic]<br>\n";
		if(!$usedCol[$ic]){
			push(@columnReorder,$ic+1);
		}
		if(!$usedRow[$ir]){
			push(@rowReorder,$ir+1);
		}
		$usedCol[$ic]=1;
		$usedRow[$ir]=1;
	}
}
open (INFO,"<$PATH_OUTPUT/$fileID.txt");
my $line = <INFO>;
chop $line;
my ($headerTop,@headers) = split(/\t/,$line);
my $nCol1 = @headers;
if($nCol1 != $nCol){ error_message("In cluster_matrix, nCol1 != nCol $nCol1 != $nCol", $logFileID); }
my @data;
my @rows;
while(my $line = <INFO>){
	chop $line;
	my ($row_header,@data1)=split(/\t/, $line);
	push(@data,\@data1);
	push(@rows,$row_header);
}
close INFO;
my $nRow1 = @rows;
if($nRow1 != $nRow){ error_message("In cluster_matrix, nRow1 != nRow", $logFileID); }
if($cluster_type==0){
	return;
}

#REORDER THE MATRIX AND SAVE IT
open (OUT, ">$PATH_OUTPUT/$fileID.txt");
print OUT "$headerTop";
if(!@columnReorder){
	for(my $i=0; $i<$nCol; ++$i){
		push(@columnReorder,$i+1);
	}
}
if(!@rowReorder){
	for(my $j=0; $j<$nRow; ++$j){
		push(@rowReorder,$j+1);
	}
}
my $replic=1;
my %hashRep;
for(my $i=0; $i<$nCol; ++$i){
	if($headers[$i] !~ /rep\d+$/){ $replic=0; last; }
	$hashRep{$headers[$i]}=$i+1;
}
my @columnReorder1;
if($replic){
	my %hash;
	for(my $i=0; $i<$nCol; ++$i){
		my @items = split(/\s+rep/,$headers[$columnReorder[$i]-1]);
		my $irep = pop(@items);
		my $name = join(" rep",@items);
		if($hash{$name}){ next; }
		$hash{$name}=1;
		for(my $j=1; $j<1000; ++$j){
			my $i1 = $hashRep{"$name rep$j"};
			if(!$i1){ last; }
			push(@columnReorder1,$i1);
		}
	}
	@columnReorder = @columnReorder1;
}
for(my $i=0; $i<$nCol; ++$i){
	print OUT "\t$headers[$columnReorder[$i]-1]";
}
print OUT "\n";
for(my $j=0; $j<$nRow; ++$j){
	my $jj = $rowReorder[$j]-1;
	print OUT $rows[$jj];
	for(my $i=0; $i<$nCol; ++$i){
		print OUT "\t$data[$jj]->[$columnReorder[$i]-1]";
	}
	print OUT "\n";
}
close OUT;
return;
}

#*****************************************************
sub  pca_based_clustering
#*****************************************************
{
my $data = shift;

my $pcMax;
my $N = @$data;
my $npc = @{$data->[0]->[2]};
for(my $ipc=0; $ipc<$npc; ++$ipc){
	my ($xmax,$xmin) = (-100000,100000);
	for(my $i=0; $i<$N; ++$i){
		my $x = $data->[$i]->[2]->[$ipc];
		if($xmax < $x){ $xmax=$x; }
		if($xmin > $x){ $xmin=$x; }
	}
	if(abs($xmin) > $xmax){
		for(my $i=0; $i<$N; ++$i){
			$data->[$i]->[2]->[$ipc] *= (-1);
		}
	}
}
my @sort1 = sort {$a->[2]->[0]<=>$b->[2]->[0]} @$data;
if($N<30){
	$pcMax = $sort1[$N-1]->[2]->[0];
}else{
	$pcMax = $sort1[int(0.99*($N-0.5))]->[2]->[0];
}
my $thresh = $pcMax/4;
#print "Thresh = $thresh<br>\n";
for(my $ipc=0; $ipc<$npc; ++$ipc){
	foreach my $ref (@$data){
		$ref->[2]->[$ipc] = int(10000*$ref->[2]->[$ipc])/10000;
	}
}
my @allClusters;
my @residuals = @$data;
for(my $ipc0=0; $ipc0<$npc; ++$ipc0){
	my @clusters=();
	my @sorted = sort {$b->[2]->[$ipc0]<=>$a->[2]->[$ipc0]} @residuals;
	my @positive;
	my $i=0;
	while($i<@sorted && $sorted[$i]->[2]->[$ipc0] > $thresh){
		my $max = 0;
		for(my $j=$ipc0+1; $j<$npc; ++$j){
			if($max < abs($sorted[$i]->[2]->[$j])){ $max=abs($sorted[$i]->[2]->[$j]); }
		}
		if($sorted[$i]->[2]->[$ipc0] > $max/1.5){
			push(@positive,splice(@sorted,$i,1));
			$i--;
		}
		$i++;
	}
	my @residuals1 = @positive;
	for(my $ipc1=$ipc0+1; $ipc1<$npc; ++$ipc1){
		my @sortPair = sort {$b->[2]->[$ipc0]*$b->[2]->[$ipc1]<=>$a->[2]->[$ipc0]*$a->[2]->[$ipc1]} @residuals1;
		my @cluster1;
		my @cluster2;
		for(my $i=0; $i<@sortPair; ++$i){
			my $ratio = $sortPair[$i]->[2]->[$ipc1] / $sortPair[$i]->[2]->[$ipc0];
			my $max = 0;
			for(my $j=$ipc0+1; $j<$npc; ++$j){
				if($j==$ipc1){ next; }
				if($max < abs($sortPair[$i]->[2]->[$j])){ $max=abs($sortPair[$i]->[2]->[$j]); }
			}
			if($sortPair[$i]->[2]->[$ipc1] > $thresh && $ratio<3 && $ratio>0.333){
				if($sortPair[$i]->[2]->[$ipc1] >= $max){
					push(@cluster1, splice(@sortPair,$i,1));
					--$i;
				}
			}
			elsif($sortPair[$i]->[2]->[$ipc1] < -$thresh && $ratio>-3 && $ratio<-0.333){
				if($sortPair[$i]->[2]->[$ipc1] <= -$max){
					splice(@cluster2,0,0,splice(@sortPair,$i,1));
					--$i;
				}
			}
		}
		@residuals1 = @sortPair;
		my $n = @cluster1;
		if($n >= 5){ push(@clusters,[$ipc0,$ipc1,1,1,$n,\@cluster1]); }
		else{ push(@residuals1,@cluster1); }
		$n = @cluster2;
		if($n >= 5){ push(@clusters,[$ipc0,$ipc1,1,-1,$n,\@cluster2]); }
		else{ push(@residuals1,@cluster2); }
	}
	my @sorted1 = sort {$b->[2]->[$ipc0]<=>$a->[2]->[$ipc0]} @residuals1;
	my $n = @sorted1;
	push(@clusters,[$ipc0,-1,1,0,$n,\@sorted1]);
	@clusters = sort {$b->[4]<=>$a->[4]} @clusters;
	push(@allClusters,@clusters);

	@clusters=();
	my @negative;
	$i = @sorted-1;
	while($i>=0 && $sorted[$i]->[2]->[$ipc0] < -$thresh){
		my $max = 0;
		for(my $j=$ipc0+1; $j<$npc; ++$j){
			if($max < abs($sorted[$i]->[2]->[$j])){ $max=abs($sorted[$i]->[2]->[$j]); }
		}
		if($sorted[$i]->[2]->[$ipc0] < -$max/1.5){
			push(@negative,splice(@sorted,$i,1));
		}
		$i--;
	}
	my @residuals2 = @negative;
	for(my $ipc1=$ipc0+1; $ipc1<$npc; ++$ipc1){
		my @sortPair = sort {$b->[2]->[$ipc0]*$b->[2]->[$ipc1]<=>$a->[2]->[$ipc0]*$a->[2]->[$ipc1]} @residuals2;
		my @cluster1;
		my @cluster2;
		for(my $i=0; $i<@sortPair; ++$i){
			my $ratio = $sortPair[$i]->[2]->[$ipc1] / $sortPair[$i]->[2]->[$ipc0];
			if($sortPair[$i]->[2]->[$ipc1] < -$thresh && $ratio<3 && $ratio>0.333){
				push(@cluster1, splice(@sortPair,$i,1));
				$i--;
			}
			elsif($sortPair[$i]->[2]->[$ipc1] > $thresh && $ratio>-3 && $ratio<-0.333){
				splice(@cluster2,0,0,splice(@sortPair,$i,1));
				$i--;
			}
		}
		@residuals2 = @sortPair;
		my $n = @cluster1;
		if($n >= 5){ push(@clusters,[$ipc0,$ipc1,-1,-1,$n,\@cluster1]); }
		else{ push(@residuals2,@cluster1); }
		$n = @cluster2;
		if($n >= 5){ push(@clusters,[$ipc0,$ipc1,-1,1,$n,\@cluster2]); }
		else{ push(@residuals2,@cluster2); }
	}
	my @sorted3 = sort {$a->[2]->[$ipc0]<=>$b->[2]->[$ipc0]} @residuals2;
	$n = @sorted3;
	push(@clusters,[$ipc0,-1,-1,0,$n,\@sorted3]);
	@clusters = sort {$b->[4]<=>$a->[4]} @clusters;
	push(@allClusters,@clusters);
	@residuals = @sorted;
}
@residuals = sort {$b->[2]->[0]<=>$a->[2]->[0]} @residuals;
my $n = @residuals;
push(@allClusters,[-1,-1,0,0,$n,\@residuals]);
my @reorder;
my $ic = $N+1;
my $root = $N + @allClusters + 1;
foreach my $cluster (@allClusters){
	my ($ipc0,$ipc1,$dir1,$dir2,$n,$ref) = @$cluster;
	#print "$ipc0, $ipc1, $dir1, $dir2, $n<br>\n";
	for(my $i=0; $i<$n; ++$i){
		push(@reorder,$ref->[$i]->[0]+1); 
	}
}
return @reorder;
}

#******************************************
sub   combine_clusters
#******************************************
{
my $dist = shift;
my $c = shift;
my $ic = shift;
my $c1 = shift;
my $c2 = shift;
my $n1 = shift;
my $n = shift;

my ($w1,$w2) = (1,1);
my $ic1 = $c->[$ic]->[0];
if($ic1 >= $n){
	$ic1 -= $n;
	$w1 = $c->[$ic1]->[3];
	if($c->[$ic1]->[5]->[$c2] < 0){ $c->[$ic1]->[4] = -1; }
	else{ $c->[$ic1]->[4] = 1; }
}
my $ic2 = $c->[$ic]->[1];
if($ic2 >= $n){
	$ic2 -= $n;
	$w2 = $c->[$ic2]->[3];
	if($c->[$ic2]->[5]->[$c1] < 0){ $c->[$ic2]->[4] = -1; }
	else{ $c->[$ic2]->[4] = 1; }
}
$c->[$ic]->[3] = $w1+$w2;
$c->[$ic]->[4] = 1;
for(my $i=0; $i<$n1; $i++){
	$c->[$ic]->[5]->[$i] = $dist->[$c1]->[$i] - $dist->[$c2]->[$i];
	$dist->[$i]->[$n1] = ($w1*$dist->[$c1]->[$i] + $w2*$dist->[$c2]->[$i])/($w1+$w2);
	$dist->[$n1]->[$i] = $dist->[$i]->[$n1]
}
for(my $i=0; $i<$n1+1; ++$i){
	for(my $j=$c1+1; $j<$n1+1; ++$j){
		my $j1 = $j-1;
		if($j>$c2){ $j1--; }
		$dist->[$i]->[$j1]=$dist->[$i]->[$j];
	}
}
my ($swap1,$swap2)=($dist->[$c1],$dist->[$c2]);
for(my $i=$c1; $i<=$c2-2; ++$i){
	$dist->[$i] = $dist->[$i+1];
}
for(my $i=$c2-1; $i<$n1-1; ++$i){
	$dist->[$i] = $dist->[$i+2];
}
($dist->[$n1-1],$dist->[$n1]) = ($swap1,$swap2);
for(my $i=0; $i<=$ic; $i++){
	my $x = ($w1*$c->[$i]->[5]->[$c1] + $w2*$c->[$i]->[5]->[$c2])/($w1+$w2);
	for(my $j=$c1+1; $j<$n1; ++$j){
		my $j1 = $j-1;
		if($j>$c2){ $j1--; }
		$c->[$i]->[5]->[$j1] = $c->[$i]->[5]->[$j];
	}
	$c->[$i]->[5]->[$n1-2] = $x;
}
return;
}

#*****************************************************
sub  sort_diagonal
#*****************************************************
{
my $matrix = shift;
my $fixedCol = shift;
my $fixedRow = shift;

my @pairs;
my $nRow = @$matrix;
my $nCol = @{$matrix->[0]};
my @usedCol;
my @usedRow;
if($fixedCol){
	for(my $j=0; $j<$nCol; $j++){
		my $i1=-1;
		my $max = -1000000000;
		for(my $i=0; $i<$nRow; $i++){
			if($max < $matrix->[$i]->[$j]){
				$max = $matrix->[$i]->[$j];
				$i1 = $i;
			}
		}
		push(@pairs,[$j,$i1]);
		$usedCol[$j]=1;
		$usedRow[$i1]=1;
	}
}	
if($fixedRow){
	for(my $i=0; $i<$nRow; ++$i){
		my $j1=-1;
		my $max = -1000000000;
		for(my $j=0; $j<$nCol;++$j){
			if($max < $matrix->[$i]->[$j]){
				$max = $matrix->[$i]->[$j];
				$j1 = $j;
			}
		}
		push(@pairs,[$j1,$i]);
		$usedCol[$j1]=1;
		$usedRow[$i]=1;
	}
}	
my %hash;
for(my $i=0; $i<$nRow; ++$i){
	for(my $j=0; $j<$nCol; ++$j){
		$hash{"$j,$i"} = $matrix->[$i]->[$j];
	}
}
#foreach my $ref (@pairs){ print "N $ref->[0] $ref->[1]<br>\n"; }
foreach my $key (sort {$hash{$b}<=>$hash{$a}} keys %hash){
	my ($j,$i) = split(/,/,$key);
	my $score = $hash{$key};
	if($score <= 0){ last; }
	if($usedCol[$j] && $usedRow[$i]){ next; }
	if(!$usedCol[$j] && !$usedRow[$i]){
		push(@pairs,[$j,$i,$score]);
		$usedCol[$j]=1;
		$usedRow[$i]=1;
		next;
	}
	my $ip=@pairs-1;
	if($usedCol[$j]){
		while($ip>0 && $pairs[$ip]->[0] != $j){ $ip--; }
		if($score > $pairs[$ip]->[2]/1.5){
			my $score1 = $pairs[$ip]->[2];
			$ip++;
			#if($ip<@pairs && ($ip>0 || $pairs[$ip-1]->[0]==$pairs[$ip]->[0] || $pairs[$ip-1]->[1]==$pairs[$ip]->[1])){ $ip++; }
			if($ip<@pairs){
				splice(@pairs,$ip,0,[$j,$i,$score1]);
			}else{
				push(@pairs,[$j,$i,$score]);
			}
		}else{
			push(@pairs,[$j,$i,$score]);
		}
		$usedRow[$i]=1;
	}else{
		while($ip>0 && $pairs[$ip]->[1] != $i){ $ip--; }
		if($score > $pairs[$ip]->[2]/1.5){
			my $score1 = $pairs[$ip]->[2];
			$ip++;
			#if($ip<@pairs && ($ip>0 || $pairs[$ip-1]->[0]==$pairs[$ip]->[0] || $pairs[$ip-1]->[1]==$pairs[$ip]->[1])){ $ip++; }
			if($ip<@pairs){
				splice(@pairs,$ip,0,[$j,$i,$score1]);
			}else{
				push(@pairs,[$j,$i,$score]);
			}
		}else{
			push(@pairs,[$j,$i,$score]);
		}
		$usedCol[$j]=1;
	}
}

my $ip=@pairs-1;
my $i = $pairs[$ip]->[1];
for(my $j=0; $j<$nCol; ++$j){
	if(!$usedCol[$j]){
		push(@pairs,[$j,$i,$matrix->[$i]->[$j]]);
	}
}
my $j = $pairs[$ip]->[0];
for(my $i=0;$i<$nRow; ++$i){
	if(!$usedRow[$i]){
		push(@pairs,[$j,$i,$matrix->[$i]->[$j]]);
	}
}
return @pairs;
}

#**************************************
sub  geneset_explore
#**************************************
{
my $file_geneset = $hashInput{"file_geneset"};
my $organismID = $hashInput{"organismID"};
my $organismID1 = $hashInput{"organismID1"};
if(!$organismID1){ $organismID1=$organismID; }
if($hashInput{"action_geneset"} =~ /^output_table_/){
	my $fileOutputID = $hashInput{"action_geneset"};
	$fileOutputID =~ s/^output_table_//;
	$hashInput{"file_output"} = "$fileOutputID.txt";
	output_explore();
	print "ERR2";
	exit(0);
}
my $sort_genesets = 0;
if($hashInput{"sorting_genesets"} eq "on"){ $sort_genesets = 1; }
my @geneset_list = get_geneset_list();
my $description_geneset = $hashInput{"description_geneset"};
if(!$description_geneset){
	foreach my $ref (@geneset_list){
		if($ref->[0] eq $file_geneset){ $description_geneset=$ref->[1]; last; }
	}
	if(!$description_geneset){
		$description_geneset="No description";
	}
}
my %hashOrganismID;
foreach my $ref (@geneset_list){
	$hashOrganismID{$ref->[2]}=1;
}
my @geneset_list_conspecific = @geneset_list;
filter_list_by_organism(\@geneset_list_conspecific, $organismID);
for(my $i=0; $i<@geneset_list_conspecific; $i++){
	if($geneset_list_conspecific[$i]->[0] =~ /^public-/ || $geneset_list_conspecific[$i]->[0] eq $file_geneset){
		splice(@geneset_list_conspecific,$i--,1);
	}
}
filter_list_by_organism(\@geneset_list, $organismID1);
my ($items,$descriptions) = get_array_lists(\@geneset_list);
my $file_geneset_full = $file_geneset;
if($file_geneset !~ /^public-/){ $file_geneset_full = "$loginname-$file_geneset"; }
my @geneset;
my %titles_short;
my $count=0;
my $table;
my $delete_geneset = -1;
my $copy_geneset = -1;
my $file_updateID;
my $message;
if($hashInput{"action_geneset"} eq "delete"){
	$delete_geneset = $hashInput{"geneset_item"};
}elsif($hashInput{"action_geneset"} eq "copy"){
	$copy_geneset = $hashInput{"geneset_item"};
}
if($delete_geneset >=0 || $copy_geneset >=0){
	$file_updateID = get_outputID(1);
	open(OUT,">$PATH_OUTPUT/$file_updateID.txt") or error_message("Cannot open temporary file");
}
my $count_non_updown=0;
my @attribNames;
my %hashAttrib;
my %hashGeneset;
open(INFO,"<$PATH_DATA/$file_geneset_full") or error_message("Cannot open file_geneset");
while(my $line=<INFO>){
	chop $line;
	if($line =~ /^[#!]/ || !$line){
		if($line =~ /^!/){
			my ($key,$value) = split(/\t/,$line);
			$key =~ s/^!G/g/;
			if($value){ $hashGeneset{$key}=$value; }
		}
		if($delete_geneset >=0){ print OUT $line."\n"; }
		next;
	}
	my ($title,$descrip,@genes) = split(/\t/,$line);
	my $nGenes = @genes;
	if(!$title){
		if($copy_geneset>=0 && $copy_geneset == $count-1){ print OUT $line."\n"; }
		if($delete_geneset>=0 && $delete_geneset != $count-1){ print OUT $line."\n"; }
		if($delete_geneset<0 || $delete_geneset != $count-1){
			my ($min,$max) = get_bounds(\@genes);
			if($min == $MISSING){ next; }
			my $ref = $hashAttrib{$descrip};
			if(!$ref){
				$hashAttrib{$descrip}=[$min,$max];
				push(@attribNames,$descrip);
			}else{
				my($min1,$max1) = @$ref;
				if($min<$min1){ $ref->[0] = $min; }
				if($max>$max1){ $ref->[1] = $max; }
			}
		}
		next;
	}
	if($copy_geneset>=0 && $copy_geneset == $count){
		my $file = $hashInput{"file_geneset_cons"};
		if($file){
			$message = "Geneset: $title - $descrip (N = $nGenes) copied to file $file";
		}else{
			$message = "Destination file $file not found. Copying geneset aborted.";
		}
		print OUT $line."\n";
	}
	if($delete_geneset>=0){
		if($delete_geneset == $count){
			$count++;			
			$message = "Geneset: $title - $descrip (N = $nGenes) deleted";
			next;
		}
		print OUT $line."\n";
	}
	$table .= "$title\t$descrip\t$nGenes\t".join(',',@genes)."\n";
	if($title =~ /_up$/i){
		my $short = $title;
		$short =~ s/_up$//i;
		$titles_short{$short}->[0]=$nGenes;
	}elsif($title =~ /_down$/i){
		my $short = $title;
		$short =~ s/_down$//i;
		$titles_short{$short}->[1]=$nGenes;
	}else{ $count_non_updown++; }	
	my $descripTerm = quotemeta($descrip);
	if($descrip && $title !~ /^$descripTerm$/i && $descrip !~ /^http|^html/){
		$title .= " - $descrip";
	}
	my $count1 = $count;
	if($delete_geneset>=0 && $count1 > $delete_geneset){
		$count1--;
	}
	push(@geneset,[$title,$count1,$nGenes]);
	$count++;
}
close INFO;
my $source_matrix = $hashGeneset{"geneset_source_matrix"};
my $FDR_matrix = $hashGeneset{"geneset_FDR"};
my $fold_enrichment = $hashGeneset{"geneset_fold_enrichment"};
my $expr_threshold = $hashGeneset{"geneset_expr_threshold"};
my $specific = $hashGeneset{"geneset_specific"};
my $legendID;
if($source_matrix){
	my $legend_file = $source_matrix; 
	$legend_file =~ s/^public-/public-legend-/;
	if($source_matrix !~ /^public-/){
		$legend_file = "$loginname-legend-$source_matrix";
	}
	if(file_exist("$PATH_DATA/$legend_file")){
		$legendID = get_outputID(1);
		copy "$PATH_DATA/$legend_file", "$PATH_OUTPUT/$legendID.txt";
	}
}
if($file_updateID){
	close OUT;
	if($delete_geneset>=0){
		copy "$PATH_OUTPUT/$file_updateID.txt", "$PATH_DATA/$file_geneset_full";
	}elsif($copy_geneset>=0){
		my $file_destination = $hashInput{"file_geneset_cons"};
		my $file_full = $loginname."-".$file_destination;
		open(INFO,"<$PATH_DATA/$file_full") or error_message("Cannot open file_full");
		my $found;
		my %hashAttrib;
		my $line_no;
		while(my $line=<INFO>){
			chop $line;
			if($line =~ /^!/ || !$line){ next; }
			my ($title,$descrip,@genes) = split(/\t/,$line);
			if($title){
				$line_no = 1;
				my $title1 = quotemeta($title);
				if($message =~ /^Geneset\: $title1 \-/){
					$found = $title;
				}
			}elsif($descrip && !$hashAttrib{$descrip}){
				$hashAttrib{$descrip} = $line_no++;
			}
		}
		close INFO;
		if($found){
			$message = "Geneset $found already exists in file $file_destination. Copying aborted";
		}elsif($file_destination){
			my @attributes;
			open(INFO,"<$PATH_OUTPUT/$file_updateID.txt") or error_message("Cannot open file_update");
			my $first_line=<INFO>;
			my ($title,$descrip,@genes) = split(/\t/,$first_line);
			while(my $line=<INFO>){
				my @items = split(/\t/,$line);
				if($items[0] || !$items[1]){ error_message("First item is not empty or second is empty."); }
				my $ii = $hashAttrib{$items[1]};
				if($ii){
					$attributes[$ii-1] = $line;
				}
			}
			close INFO;

			my @keys = sort {$hashAttrib{$a}<=>$hashAttrib{$b}} keys %hashAttrib;
			open(OUT, ">>$PATH_DATA/$file_full") or error_message("Cannot write to file_full");
			print OUT $first_line;
			for(my $i=0; $i<@keys; $i++){
				my $xx = $attributes[$i];
				if($xx){ print OUT $xx; }
				else{
					print OUT "\t$keys[$i]";
					foreach my $i1 (1..@genes){ print OUT "\t"; }
					print OUT "\n";
				}
			}
			close OUT;
		}
	}
	unlink "$PATH_OUTPUT/$file_updateID.txt";
}
my $up_down = 0;
my $n1 = @geneset;
my $n2 = keys %titles_short;
if($count_non_updown==0){ $up_down = 1; }
my $plotID1 = get_outputID(3);
my $tableID = $plotID1+2;
open(OUT,">$PATH_OUTPUT/$tableID.txt") or error_message("Cannot open file");
if(!$up_down){
	my @sorted = sort {$b->[2]<=>$a->[2]} @geneset;
	my $N = @sorted;
	if($N>200){ $N=200; }
	my @data;
	my @titles;
	for(my $i=0; $i<$N; ++$i){
		my $ref = $sorted[$i];
		push(@data,$ref->[2]);
		push(@titles,$ref->[0]);
	}
	plot_histogram_horiz("bars",$N,1,\@data,\@titles,$plotID1,"N genes","zero");
}else{
	print OUT "Gene set\tN genes upregulated\tN genes downregulated\tN total\n";
	my $N = keys %titles_short;
	if($N>200){ $N=200; }
	my @data;
	my @titles;
	foreach my $title (sort {$titles_short{$b}->[0]+$titles_short{$b}->[1]<=>$titles_short{$a}->[0]+$titles_short{$a}->[1]} keys %titles_short){
		my $ref = $titles_short{$title};
		if(!$ref->[0]){ $ref->[0]=0; }
		if(!$ref->[1]){ $ref->[1]=0; }
		my $sum = $ref->[0]+$ref->[1];
		push(@data,$ref);
		push(@titles,$title);
		print OUT "$title\t$ref->[0]\t$ref->[1]\t$sum\n";

	}
	plot_histogram_horiz("stacked",$N,2,\@data,\@titles,$plotID1,"N genes","zero");
}
my $nGeneset=@geneset;
print OUT "\n$table";
close OUT;
print "<HTML><HEAD><TITLE>ExAtlas - geneset</TITLE>\n";
print_header("update_description();");
print "<script type=text/JavaScript language=JavaScript>\n";
print "<!-- \n";
print "geneset_description = new Array($descriptions);\n";
if(@geneset_list_conspecific){
	my ($items_cons,$descriptions_cons) = get_array_lists(\@geneset_list_conspecific);
	print "geneset_cons_description = new Array($descriptions_cons);\n";
}
print "function update_description() {\n";
print "	var index;\n";
print "	index = document.form_geneset.file_geneset1.selectedIndex;\n";
print "	document.form_geneset.description_geneset1.value = geneset_description[index];\n";
if(@geneset_list_conspecific){
	print "	index = document.form_geneset.file_geneset_cons.selectedIndex;\n";
	print "	document.form_geneset.description_geneset_cons.value = geneset_cons_description[index];\n";
}
print "}\n";
print "function alert_onsubmit() {\n";
print "	document.form_geneset.action.value = \"geneset_overlap\";\n";
print "	document.form_geneset.action_geneset.value = \"\";\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_geneset.target = \"_blank\"+x;\n";
print "	document.form_geneset.submit();\n";
print "}\n";
print "function geneset_search(){\n";
print "	if(!document.form_geneset.search_term.value){\n";
print "		alert(\"You need to put a search term\");\n";
print "		return(false);\n";
print "	}\n";
print "	document.form_geneset.action.value = \"geneset_explore1\";\n";
print "	document.form_geneset.action_geneset.value = \"search\";\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_geneset.target = \"_blank\"+x;\n";
print "	document.form_geneset.submit();\n";
print "}\n";
print "function change_organism() {\n";
print "	document.form_geneset.action_geneset.value = \"\";\n";
print "	document.form_geneset.action.value = \"geneset_explore\";\n";
print "	document.form_geneset.target = \"\";\n";
print "	document.form_geneset.submit();\n";
print "}\n";

print "function sort_genesets() {\n";
print "	document.form_geneset.action_geneset.value = \"\";\n";
print "	document.form_geneset.action.value = \"geneset_explore\";\n";
print "	document.form_geneset.target = \"\";\n";
print "	document.form_geneset.submit();\n";
print "}\n";

print "function select_geneset(){\n";
print "	document.form_geneset.action.value = \"geneset_explore1\";\n";
print "	document.form_geneset.action_geneset.value = \"number\";\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_geneset.target = \"_blank\"+x;\n";
print "	document.form_geneset.submit();\n";
print "}\n";
print "function delete_geneset(){\n";
print "	if(!confirm(\"Do you want to delete this geneset?\")){\n";
print "		return false;\n";
print "	}\n";
print "	document.form_geneset.action.value = \"geneset_explore\";\n";
print "	document.form_geneset.action_geneset.value = \"delete\";\n";
print "	document.form_geneset.target = \"\";\n";
print "	document.form_geneset.submit();\n";
print "}\n";
print "function copy_geneset(){\n";
print "	if(!confirm(\"Do you want to copy this geneset?\")){\n";
print "		return false;\n";
print "	}\n";
print "	document.form_geneset.action.value = \"geneset_explore\";\n";
print "	document.form_geneset.action_geneset.value = \"copy\";\n";
print "	document.form_geneset.target = \"\";\n";
print "	document.form_geneset.submit();\n";
print "}\n";
print "function open_output(output_file){\n";
print "	document.form_geneset.action.value = \"geneset_explore\";\n";
print "	document.form_geneset.action_geneset.value = \"output_table_\"+output_file;\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_geneset.target = \"_blank\"+x;\n";
print "	document.form_geneset.submit();\n";
print "}\n";
print "<!-- end script --></SCRIPT>\n";
print "<H2>Geneset file '$file_geneset'</H2>\n";
if($message){
	print "<font color=blue><b>Message:</b> $message</font><p>\n";
}
print "<TABLE BORDER=0><TR><TD VALIGN=top>\n";
print "<b>Organism:</b> $hashOrganism{$organismID}<br>\n";
if($description_geneset){
	my $description1 = add_hyperlinks($description_geneset);
	print "<b>Geneset file description:</b> $description1<br>\n";
}
print "<b>Number of genesets:</b> $nGeneset<br>\n";
my $x = int(10000*rand());
print "<b>Genesets: tab-delimited text:</b> <a href=$HOME_ADDRESS/output/$tableID.txt target=_blank$x>text file</a><br>\n";
if($source_matrix){ print "<b>Source matrix:</b> $source_matrix<br>\n"; }
if($FDR_matrix){ print "<b>FDR_matrix:</b> $FDR_matrix<br>\n"; }
if($fold_enrichment){ print "<b>Fold_enrichment:</b> $fold_enrichment<br>\n"; }
if($expr_threshold){ print "<b>Expression threshold:</b> $expr_threshold<br>\n"; }
if($specific){ print "<b>Specific:</b> $specific<br>\n"; }
if($legendID){
	my $x = int(10000*rand());
	print "<b>Legend:</b> <a href=$HOME_ADDRESS/output/$legendID.txt target=_blank$x>Legend file</a>\n";
}
my $fileOutputID = $hashGeneset{"geneset_output_table"};
if($fileOutputID){
	print "<b>Tables:</b> <INPUT type=button value=\"Open tables\" onClick=open_output($fileOutputID); style=width:200px;>\n";
}
if(@attribNames){
	print "<TD WIDTH=60><TD VALIGN=top>\n";
	print "<b>Geneset attributes</b>\n";
	print "<TABLE BORDER=0>\n";
	print "<TR><TD>No.<TD>Attribute<TD ALIGN=center>Minimum<TD ALIGN=center>Maximum\n";
	for(my $i=0; $i<@attribNames; $i++){
		my $i1=$i+1;
		my $ref = $hashAttrib{$attribNames[$i]};
		my($min,$max);
		if($ref){ ($min,$max) = @$ref; }
		print "<TR><TD>$i1.<TD>$attribNames[$i]<TD>$min<TD>$max\n";
	}
	print "\n</TABLE>";
}
print "\n</TABLE><p>";
print "<b>Number of genes in gene sets</b><br>\n";
if($up_down){
	print "Orange = upregulated genes, dark blue = downregulated genes\n";
}
print "<p><IMG SRC=../output/$plotID1.gif BORDER=0><p>\n";
print "<FORM NAME=form_geneset ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<FONT SIZE=+1><b>1. Find and display a gene set (or gene)</b></FONT>\n";
print "&nbsp; &nbsp; &nbsp; <INPUT NAME=sorting_genesets type=checkbox onChange=sort_genesets();";
if($sort_genesets){ print " checked"; }
print "> Sort geneset names<br>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD WIDTH=200>Search term(s)<br>(comma-separated)<TD><INPUT NAME=search_term style=width:280px;>\n";
print "<TD><INPUT type=button value=Search onClick=geneset_search(); style=width:200px;>\n";
print "<TD WIDTH=15><TD><SELECT NAME=search_type style=width:200px;>\n";
print "<option value=geneset>Find genesets \n<option value=gene>Gene symbols</SELECT>\n";
print "<TR><TD WIDTH=200>Or select from list:\n";
print "<TD><select name=geneset_item style=width:280px;>\n";
my @sorted;
if($sort_genesets){ @sorted = sort {lc($a->[0]) cmp lc($b->[0])} @geneset; }
else{ @sorted = @geneset; }
for(my $i=0; $i<@sorted; ++$i){
	my $ref = $sorted[$i];
	my $name = $ref->[0];
	if(length($name) > 70){ $name = substr($name,0,70); }
	my $num = $ref->[1];
	print "<option value=$num>$name\n";
}
print "</select>\n";
print "<TD><INPUT type=button value=\"Display genes\" onClick=select_geneset(); style=width:200px;><p>\n";
if($file_geneset !~ /^public-/ || $loginname eq "public"){
	print "<TD WIDTH=15><TD><INPUT type=button value=\"Delete geneset\" onClick=delete_geneset(); style=width:200px;><p>\n";
}
if(@geneset_list_conspecific){
	print "<TR><TD>Copy geneset to file:\n";
	print "<TD><select name=file_geneset_cons onChange=update_description(); style=\"width: 280px;\">\n";
	for(my $i=0; $i<@geneset_list_conspecific; ++$i){ 
		print "<option value=\"$geneset_list_conspecific[$i]->[0]\"> $geneset_list_conspecific[$i]->[0]\n";
	}
	print "</select><TD><INPUT NAME=\"description_geneset_cons\" style=\"width: 200px;\">\n";
	print "<TD WIDTH=15><TD><INPUT type=button value=\"Copy geneset\" onClick=copy_geneset(); style=width:200px;><p>\n";
}
print "</TABLE><p>\n";

print "<TABLE BORDER=0><TR><TD WIDTH=483>\n";
print "<FONT SIZE=+1><b>2. Analyze overlap with another geneset file</b></FONT>\n";
print "<TD><INPUT TYPE=button value=\"Overlap analysis\" LANGUAGE=javascript onClick=\"alert_onsubmit();\" style=\"width: 200px;\">\n";
print "</TABLE>\n";
print "<TABLE BORDER=0>\n";
if(keys %hashOrganismID > 1){
	print "<TR><TD WIDTH=200>Change organism for 2nd file:\n";
	print "<TD><select name=organismID1 style=width:280px; onChange=change_organism();>\n";
	foreach my $ref (@organisms){
		if(!$hashOrganismID{$ref->[0]}){ next; }
		print "<option value=$ref->[0]";
		if($ref->[0] == $organismID1){ print " selected"; }
		print ">$ref->[3] ($ref->[2])\n";
	}
	print "</select><TD>(page is reloaded after change)\n";
}
print "<TR><TD WIDTH=200>Select file:\n";
print "<TD><select name=\"file_geneset1\" onChange=update_description(); style=\"width: 280px;\">\n";
for(my $i=0; $i<@geneset_list; ++$i){ 
	print "<option value=\"$geneset_list[$i]->[0]\"";
	if($geneset_list[$i]->[0] =~ /^public-GO/){ print " selected"; }
	print "> $geneset_list[$i]->[0]\n";
}
print "</select><TD><INPUT NAME=\"description_geneset1\" style=\"width: 280px;\">\n";
my @FDR_list = (1,0.5,0.2,0.1,0.05,0.01,0.001,0.0001);
my @fold_enrichment = (0.0001,1,1.5,2,3,4,5,10);
my @minN = (1,2,3,4,5,7,10);
my ($FDR,$fold_enrichment,$minN) = (0.05,2,5);
print "<TR><TD>Parameter: FDR threshold<TD><select name=FDR style=\"width: 140px;\">\n";
for(my $i=0; $i<@FDR_list; ++$i){ 
	print "<option value=$FDR_list[$i]"; if($FDR==$FDR_list[$i]){ print " selected"; } print "> $FDR_list[$i]\n";
}
print "</select><TD COLSPAN=2>Siginifance of overlap (hypergeometric test)\n";
print "<TR><TD>Parameter: Fold enrichment<TD><select name=fold_enrichment style=\"width: 140px;\">\n";
for(my $i=0; $i<@fold_enrichment; ++$i){ 
	print "<option value=$fold_enrichment[$i]"; if($fold_enrichment==$fold_enrichment[$i]){ print " selected"; } print "> $fold_enrichment[$i]\n";
}
print "</select><TD COLSPAN=2>Fold enrichment threshold\n";

print "<TR><TD>Parameter: N genes <TD><select name=minimum_genes style=\"width: 140px;\">\n";
for(my $i=0; $i<@minN; ++$i){ 
	print "<option value=$minN[$i]"; if($minN==$minN[$i]){ print " selected"; } print "> $minN[$i]\n";
}
print "</select><TD COLSPAN=2>N genes (min)\n";
print "<TR><TD><TD COLSPAN=3><INPUT TYPE=CHECKBOX NAME=use_attribute> Use gene attributes (if available)\n";
print "</TABLE><p>\n";
print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=action TYPE=hidden VALUE=\"geneset_explore1\">\n";
print "<INPUT NAME=file_geneset TYPE=hidden VALUE=\"$file_geneset\">\n";
print "<INPUT NAME=description_geneset TYPE=hidden VALUE=\"$description_geneset\">\n";
print "<INPUT NAME=file_type TYPE=hidden>\n";
print "<INPUT NAME=action_geneset TYPE=hidden>\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
if($legendID){ print "<INPUT NAME=legendID TYPE=hidden VALUE=$legendID>\n"; }
if($source_matrix){ print "<INPUT NAME=source_matrix TYPE=hidden VALUE=\"$source_matrix\">\n"; }
if(keys %hashOrganismID <= 1){
	print "<INPUT NAME=organismID1 TYPE=hidden VALUE=$organismID1>\n";
}
print "</FORM><p>\n";
print "<HR NOSHADE>\n";
print "<INPUT TYPE=button VALUE=\" Cancel (close window) \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  geneset_explore1
#**************************************
{
my $organismID = $hashInput{"organismID"};
my $file_geneset = $hashInput{"file_geneset"};
my $description_geneset = $hashInput{"description_geneset"};
my $legendID = $hashInput{"legendID"};
my $source_matrix = $hashInput{"source_matrix"};
my $file_geneset_full = $file_geneset;
if($file_geneset !~ /^public-/){ $file_geneset_full = "$loginname-$file_geneset"; }
my $geneset_item = $hashInput{"geneset_item"};
if(!$geneset_item){ $geneset_item=0; }
my $action_geneset = $hashInput{"action_geneset"};
if(!$action_geneset){ $action_geneset = "number"; }
my $search_type = $hashInput{"search_type"};
my $search_text = $hashInput{"search_term"};
$search_text =~ s/<//g;  #destroy hyperlinks
my @search_term = ($search_text);
if($search_text =~ /,/){
	@search_term = split(/,\s*/,$search_text);
}elsif($search_text =~ /\s/ && $search_type eq "gene"){
	@search_term = split(/\s+/,$search_text);
}
my %hashTerms;
for(my $i=0; $i<@search_term; $i++){
	$search_term[$i] =~ s/^\s+//;
	$search_term[$i] =~ s/\s+$//;
	$search_term[$i] =~ s/\s+/ /;
	$hashTerms{uc($search_term[$i])}=1;
}
if($source_matrix){
	if($source_matrix =~ /^public-/){
		if(!file_exist("$PATH_DATA/$source_matrix")){ $source_matrix=""; }
	}else{
		if(!file_exist("$PATH_DATA/$loginname-$source_matrix")){ $source_matrix=""; }
	}
}
my %hashTermsOrig = %hashTerms;
my $aliasList;
if($search_type eq "gene"){
	my %hashAlias;
	open (INFO1, "$PATH_DATA/gene_info_$organismID.txt");
	while(my $line = <INFO1>){
		chop $line;
		my ($org_id,$entrez_id,$symbol,$junk1,$alias_list,$altName,$junk3,$junk4,$geneName) = split(/\t/,$line);
		$alias_list =~ s/^-//;
		foreach my $alias (split(/\|/,$alias_list)){
			push(@{$hashAlias{uc($alias)}},$symbol);
		}
		if($hashTerms{uc($symbol)}){ $hashTerms{uc($symbol)}++; }
	}
	close INFO1;
	for(my $i=0; $i<@search_term; $i++){
		if($hashTerms{uc($search_term[$i])}==1){
			my $ref = $hashAlias{uc($search_term[$i])};
			if($ref){
				foreach my $symbol (@$ref){
					$hashTerms{uc($symbol)}++;
				}
			}
		}
	}
	@search_term = keys %hashTerms;
	for(my $i=0; $i<@search_term; $i++){
		if(!$hashTermsOrig{$search_term[$i]}){
			if(!$aliasList){ $aliasList=$search_term[$i]; }
			else{ $aliasList .= ", $search_term[$i]"; }
		}
		$search_term[$i] = quotemeta($search_term[$i]);
	}
}
if($aliasList){ $search_text .= " (Alias for: $aliasList)"; }
my @geneset=();
my $nGenes=0;
my @output;
my $count = 0;
open(INFO,"<$PATH_DATA/$file_geneset_full") or error_message("Cannot open file geneset");
while(my $line=<INFO>){
	chop $line;
	if($line =~ /^[!#]/ || !$line){ next; }
	my ($title,$descrip,@genes) = split(/\t/,$line);
	if(!$title){ next; }
	$nGenes = @genes;
	if($action_geneset eq "number" && $geneset_item==$count){
		@output = ($title,$descrip,$nGenes,\@genes);
		my @attribute;
		while($line=<INFO>){
			chop $line;
			my ($blank,@attr) = split(/\t/,$line);
			if($blank){ last; }
			push(@attribute,\@attr);
		}
		if(@attribute){ push(@output,\@attribute); }
	}
	if($action_geneset eq "search"){
		my $found=0;
		for(my $i=0; $i<@search_term; $i++){
			if($search_type eq "geneset"){
				if($title =~ /$search_term[$i]/i || $descrip =~ /$search_term[$i]/i){ $found=1; }
			}elsif($search_type eq "gene"){
				foreach my $symbol (@genes){
					if($symbol =~ /^$search_term[$i]$/i){ $found++; }
				}
			}
		}
		if($found){
			push(@output,[$title,$descrip,$count,$found,$nGenes]);
		}
	}
	$count++;
}
close INFO;
my @geneset_list=();
filter_list_by_organism(\@geneset_list, $organismID);
print "<!DOCTYPE html><HEAD><TITLE>ExAtlas - geneset</TITLE>\n";
if($nGenes>=5){
	print_header("update_description();");
	print "<script type=text/JavaScript language=JavaScript><!-- \n";
	@geneset_list = get_geneset_list();
	filter_list_by_organism(\@geneset_list, $organismID);
	my ($items,$descriptions) = get_array_lists(\@geneset_list);
	print "geneset_list = new Array($items);\n";
	print "geneset_description = new Array($descriptions);\n";
	print "function geneset_overlap(file) {\n";
	print "	document.form_geneset.upload_geneset.value=file;\n";
	print "	document.form_geneset.action.value=\"geneset_overlap\";\n";
	print "	document.form_geneset.file_geneset.value=\"\";\n";
	print "	var x = Math.round(Math.random()*10000);\n";
	print "	document.form_geneset.target = \"_BLANK\"+x;\n";
	print "	document.form_geneset.submit();\n";
	print "}\n";
	print "function update_description() {\n";
	print "	var index;\n";
	print "	index = document.form_geneset.file_geneset1.selectedIndex;\n";
	print "	document.form_geneset.description_geneset1.value = geneset_description[index];\n";
	print "}\n";
}else{
	print_header();
	print "<script type=text/JavaScript language=JavaScript><!-- \n";
}
print "function select_geneset(num){\n";
print "	document.form_geneset.action.value = \"geneset_explore1\";\n";
print "	document.form_geneset.action_geneset.value = \"number\";\n";
print "	document.form_geneset.geneset_item.value = num;\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_geneset.target = \"_blank\"+x;\n";
print "	document.form_geneset.submit();\n";
print "}\n";
print "<!-- end script --></SCRIPT>\n";
print "<FORM NAME=form_geneset ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
if($action_geneset eq "number"){
	my ($title,$description,$nGenes,$gene_ref,$attrib_ref) = @output;
	my %hashGenes=();
	get_gene_annotation($organismID,$gene_ref,\%hashGenes,"symbol");
	my $fileID = get_outputID(1);
	print "<H2>Gene set '$title'</H2>\n";
	print "<TABLE BORDER=0>\n";
	print "<TR><TD><b>Organism:<TD>$hashOrganism{$organismID}\n";
	if($description && $description ne $title){
		print "<TR><TD><b>Geneset description<TD>$description\n";
	}
	print "<TR><TD WIDTH=180><b>In file:<TD>$file_geneset\n";
	if($description_geneset && $description_geneset ne $file_geneset){
		my $description1 = add_hyperlinks($description_geneset);
		print "<TR><TD><b>File description:<TD>$description1\n";
	}
	print "<TR><TD><b>Number of genes:<TD>$nGenes\n";
	my $x = int(10000*rand());
	print "<TR><TD><b>Table of genes:<TD><a href=$HOME_ADDRESS/output/$fileID.txt target=_blank$x>text file</a>\n";
	if($legendID){
		$x+=182;
		print "<TR><TD><b>Legend:<TD><a href=$HOME_ADDRESS/output/$legendID.txt target=_blank$x>Legend file</a>\n";
	}
	print "</TABLE><p>\n";
	if($nGenes>=5){
		print "<TABLE BORDER=0>\n";
		my $menu_text = menu_geneset_overlap(\@geneset_list,$fileID);
		print $menu_text;
		print "</TABLE><p>\n";
	}
	print "<H3>List of genes</H3>\n";
	open(OUT,">$PATH_OUTPUT/$fileID.txt") or error_message("Cannot open fileID");
	print OUT "!Genelist_title\t$title\n";
	print OUT "!Genelist_description\t$description\n";
	print OUT "Gene symbol";
	print "<TABLE border=1><TR><TD><b>Gene symbol";
	my $sort_by = -1;
	my $reverse = 0;
	if($attrib_ref && ref($attrib_ref) eq "ARRAY"){
		$sort_by=0;
		if($attrib_ref->[0]->[0] =~ /^(FDR|EPFP|p|p-value)$/i){ $reverse = 1; }
		for(my $i=0; $i<@$attrib_ref; $i++){
			my $header = shift(@{$attrib_ref->[$i]});
			if($sort_by==0 && $header =~ /^logratio|^log-ratio|^ClipRPM/i){
				$sort_by = $i;
				$reverse = 0;
			}
			print OUT "\t$header";
			print "<TD><b>$header";
		}
	}
	my @sorted = (0..($nGenes-1));
	if($sort_by>=0){
		if($reverse){
			@sorted = sort {abs($attrib_ref->[$sort_by]->[$a])<=>abs($attrib_ref->[$sort_by]->[$b])} @sorted;
		}else{
			@sorted = sort {abs($attrib_ref->[$sort_by]->[$b])<=>abs($attrib_ref->[$sort_by]->[$a])} @sorted;
		}
	}
	print OUT "\tGene name\n";
	print "<TD WIDTH=350><b>Gene name\n";
	for(my $ii1=0; $ii1<$nGenes; $ii1++){
		my $ii = $sorted[$ii1];
		my $symbol = $gene_ref->[$ii];
		print OUT $symbol;
		print "<TR><TD>";
		if($source_matrix){
			my $x = int(10000*rand());
			print "<a href=\"exatlas.cgi?category=0&search_term=$symbol&category=1&action=matrix_explore1&analysis=search&file_matrix=$source_matrix&sessionID=$sessionID&organismID=$organismID\" target=_BLANK$x>";
		}
		print $symbol;
		if($source_matrix){ print "</a>"; }
		if($attrib_ref){
			for(my $i=0; $i<@$attrib_ref; $i++){
				my $x = $attrib_ref->[$i]->[$ii];
				print OUT "\t$x";
				print "<TD>$x";
			}
		}
		my $gene_name="";
		my $ref = $hashGenes{$symbol};
		if($ref && ref($ref) eq 'ARRAY'){ $gene_name=$ref->[1]; }
		print OUT "\t$gene_name\n";
		print "<TD>$gene_name\n";
	}
	print "</TABLE><p>\n";
	close OUT;
}else{
	print "<H2>Search results in file '$file_geneset'</H2>\n";
	print "<TABLE BORDER=0>\n";
	if($description_geneset && $description_geneset ne $file_geneset){ 
		print "<TR><TD><b>File description:<TD>$description_geneset\n";
	}
	print "<TR><TD><b>Organism:<TD WIDTH=800>$hashOrganism{$organismID}<br>\n";
	print "<TR><TD><b>Search term(s):<TD>$search_text\n";
	print "</TABLE>\n";
	if(!@output){
		print "<b>No items found!</b><p>\n";
	}else{
		my @sorted = (0..(@output-1));
		my $header1 = "";
		if($search_type eq "gene" && @search_term>1){
			$header1 = "<TD><b>N genes";
			@sorted = sort {$output[$b]->[3]<=>$output[$a]->[3]} @sorted;
		}
		print "<TABLE BORDER=0>\n";
		print "<TR><TD>&nbsp;<TR><TD><b>Geneset name$header1<TD><b>N-genes<TD><b>Description\n";
		for(my $i1=0; $i1<@output; $i1++){
			my $i = $sorted[$i1];
			my ($title,$description,$count,$Nhits,$nGenes) = @{$output[$i]};
			my $numHits = "";
			if($search_type eq "gene" && @search_term>1){ $numHits = "<TD ALIGN=CENTER>$Nhits"; }
			print "<TR><TD><a href=\"JavaScript: select_geneset('$count');\">$title</a>$numHits<TD>$nGenes<TD>$description\n";
		}
		print "</TABLE>\n";
	}
}
print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT NAME=action TYPE=hidden VALUE=\"geneset_explore1\">\n";
print "<INPUT NAME=file_geneset TYPE=hidden VALUE=\"$file_geneset\">\n";
print "<INPUT NAME=description_geneset TYPE=hidden VALUE=\"$description_geneset\">\n";
print "<INPUT NAME=action_geneset TYPE=hidden>\n";
print "<INPUT NAME=geneset_item TYPE=hidden>\n";
print "<INPUT NAME=upload_geneset TYPE=hidden>\n";
if($legendID){ print "<INPUT NAME=legendID TYPE=hidden VALUE=$legendID>\n"; }
print "</FORM>\n";
print "<HR NOSHADE></HR>\n";
print "<INPUT TYPE=button VALUE=\" Close the window \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub   search_GEO
#**************************************
{
my $organismID = $hashInput{"organismID"};
my @platform_list=(["All_platforms",""]);
open(INFO, "$PATH_DATA/array_platforms.txt") or error_message("Cannot open array_platforms.txt");
my $line = <INFO>;
while(my $line = <INFO>){
	chop $line;
	my($ID,$title,$technology,$taxonomy,$rows,$taxid,$type) = split(/\t/,$line);
	my $found=0;
	foreach my $id (split(/;/,$taxid)){
		if($id==$organismID){ $found=1; }
	}
	if(!$found){ next; }
	push(@platform_list,[$ID,"$title - $technology - $type - N rows=$rows"]);
}
close INFO;
my ($items,$descriptions) = get_array_lists(\@platform_list);

# Print page header
print "<HTML><HEAD><TITLE>ExAtlas</TITLE>\n";
print_header("update_description();");
print "<script type=text/JavaScript language=JavaScript>\n";
print "<!-- \n";
print "platform_list = new Array($items);\n";
print "platform_description = new Array($descriptions);\n";
print "function search_onsubmit() {\n";
print "	if(!document.form_search.search_term.value){\n";
print "		if(!confirm(\"There is no search term! Do you want to proceed?\")){;\n";
print "			return(false);\n";
print "		}\n";
print "	}\n";
print "	document.form_search.submit();\n";
print "}\n";
print "function update_description() {\n";
print "  var index;\n";
print "  index = document.form_search.platform.selectedIndex;\n";
print "  document.form_search.description_platform.value = platform_description[index];\n";
print "}\n";
print "<!-- end script --></SCRIPT></HEAD>\n";

print "<H2>Search GEO database for gene expression data</H2>\n";
print "<b>Organism:</b> $hashOrganism{$organismID}<p>\n";
print "<FORM NAME=form_search ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD><INPUT NAME=search_term  style=\"width: 300px;\"><TD WIDTH=10><TD>Search terms (comma separated!)\n";
print "<TR><TD><INPUT NAME=avoid_term  style=\"width: 300px;\"><TD WIDTH=10><TD>Avoid terms\n";

print "<TR><TD><select name=\"platform\" style=\"width: 300px;\" onChange=update_description();>\n";
for(my $i=0; $i<@platform_list; ++$i){ 
	print "<option value=\"$platform_list[$i]->[0]\"> $platform_list[$i]->[0]\n";
}
print "</select><td><td><INPUT NAME=\"description_platform\" style=width:400px;>\n";
print "<TD>Platform\n";
print "<TR><TD><select name=series_samples>\n";
print "	<option value=samples> Show individual samples\n";
print "	<option value=series> Show series of samples\n";
print "</select><TD WIDTH=10><TD>Show individual samples or series?\n";
print "<TR><TD><INPUT type=button value=\" Search \" onClick=search_onsubmit();>\n";
print "<TR><TD>&nbsp;\n";
print "</TABLE>\n";
print "<b>Notes:</b> Use comma-separated keywords (e.g., kidney, muscle, brain, ES cells) for search.\n";
print "Examples of \"Avoid\" terms: cancer,tumor,carcinoma,biopsy,patient,hepatitis,HIV.\n";
print "Search results can be filtered by platform ID (if selected).\n";
print "You can bypass filtering by putting GEO accession number(s) for data series (e.g., GSE22151)\n";
print "into the feild \'Search terms\'.<p>\n";
print "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=hidden VALUE=\"search_GEO1\">\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "</FORM>\n";
print "<HR NOSHADE></HR>\n";
print "<INPUT TYPE=button VALUE=\" Cancel (close window) \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "<p><i>Please report any problems to <a href=mailto:sharoval\@mail.nih.gov>webmaster</a><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub   search_GEO1
#**************************************
{
my $search_term = $hashInput{"search_term"};
my $contributors = $hashInput{"contributors"};
my $avoid_term = $hashInput{"avoid_term"};
my $series_samples = $hashInput{"series_samples"};
my $platform = $hashInput{"platform"};
my $organismID = $hashInput{"organismID"};

$search_term =~ s/^\s+//;
$search_term =~ s/\s+$//;
$search_term =~ s/\s+/ /g;
$search_term =~ s/,\s+/,/g;
$avoid_term =~ s/^\s+//;
$avoid_term =~ s/\s+$//;
$avoid_term =~ s/\s+/ /g;
$avoid_term =~ s/,\s+/,/g;
$avoid_term =~ s/\<//g;
$platform =~ s/^\s+//;
$platform =~ s/\s+$//;
$platform =~ s/\s+/ /g;
my @search_term = split(/,/,$search_term);
my @avoid_term = split(/,/,$avoid_term);
if($platform eq "All_platforms"){ $platform=""; }

my $N_words = 1;
for(my $i=0; $i<@search_term; ++$i){
	$search_term[$i] = quotemeta($search_term[$i]);
}
for(my $i=0; $i<@avoid_term; ++$i){
	$avoid_term[$i] = quotemeta($avoid_term[$i]);
}
my %series_info=();
my %hashSamples=();
open(INFO,"<$PATH_DATA/GEO_series_summary.txt") or error_message("Cannot open GEO_series_summary.txt");
my ($seriesID,$platformID,$title,$n_samples,$taxid,$text);
my @referenceKey = ("reference","refpool","control");
my %seriesMultiple;
my $seriesID_old;
my $search_adjusted=0;
my $Nseries_total=0;
my $minimum_score=1;
my $comments;
my $Nseries;
my $first_record=1;
while(my $line = <INFO>){
	chop $line;
	$line =~ s/\r$//;
	$line =~ s/\n$//;
	if($line =~ /^>/){
		my @matches;
		if($first_record){
			$line =~ s/^>//;
			($seriesID,$platformID,$title)=split(/\t/,$line);
			$text = $line;
			$first_record=0;
			next;
		}
		my $found = 1;
		my $avoid = 0;
		if($organismID != $taxid){ $avoid=1; }
		for(my $i=0; $i<@avoid_term && !$avoid; ++$i){
			if($text =~ /$avoid_term[$i]/i){
				$avoid=1;
			}
		}
		for(my $i=0; $i<@search_term; ++$i){
			if($seriesID eq $search_term[$i] && (!$platform || $platformID eq $platform)){
				$matches[$i]+=10000;
				$avoid=0;
				$found=1;
				if($organismID != $taxid){ error_message("Series $search_term[$i] comes from a different species.\nSelect correct species before searching GEO."); }
			}else{
				while($text =~ /$search_term[$i]/ig){
					$matches[$i]++;
				}
			}
		}
		if(!@search_term){
			$matches[0] = 1;
		}
		if(@matches && $platform && $found){
			$found=0;
			if($platformID eq $platform){ $found=1; }
		}
		if(!@matches || !$found || $avoid){
			$Nseries_total++;
			$line =~ s/^>//;
			($seriesID,$platformID,$title)=split(/\t/,$line);
			$text = $line;
			next;
		}
		my $score=0;
		my $total_match=0;
		for(my $i=0; $i<@matches; ++$i){
			if($matches[$i]){
				$score++;
				$total_match += int(log($matches[$i]));
			}
		}
		if($total_match>29){ $total_match=29; }
		$score += $total_match/30;
		if(!$n_samples){ $score=0; }
		if($score >= $minimum_score){
			my @lines = split(/\n/,$text);
			while(@lines && $lines[0] !~ /^GSM\d+\t/i){ shift(@lines); }
			my $seriesIDext = $seriesID;
			my $ref = $series_info{$seriesID};
			if($ref){
				my $num = 1;
				while($series_info{"$seriesID-$num"}){ $num++; }
				$seriesIDext = "$seriesID-$num";
			}
			my $nDye = @lines/$n_samples;
			my @countRef=(0,0,0);
			for(my $i=0; $i<@lines; ++$i){
				my($sampleID,$descr,@items)=split(/\t/,$lines[$i]);
				if($nDye==2){
					for(my $i1=0; $i1<@referenceKey; ++$i1){
						if($lines[$i] =~ /$referenceKey[$i1]/i){ $countRef[$i1]++; }
					}
				}
				my $seriesID1 = $hashSamples{$sampleID};
				if($seriesID1 && $seriesID1 ne $seriesIDext){
					my $refSamples = $series_info{$seriesID1}->[4];
					for(my $j=0; $j<@$refSamples; ++$j){
						if($refSamples->[$j] =~ /^$sampleID\t/){
							splice(@$refSamples,$j,1);
							$j--;
						}
					}
					$series_info{$seriesID1}->[2]--;
				}
				$hashSamples{$sampleID} = $seriesIDext;
			}
			my $control="";
			for(my $i1=0; $i1<@referenceKey; ++$i1){
				if($countRef[$i1]==$n_samples){ $control = $referenceKey[$i1]; }
			}
			for(my $i=0; $i<@lines; ++$i){
				if($control && $lines[$i] =~ /$control/){
					splice(@lines,$i,1);
					$i--;
					next;
				}
				my($sampleID,$descr)=split(/\t/,$lines[$i]);
				$lines[$i] = "$sampleID\t$descr";
			}
			$nDye=1;
			if($n_samples && @lines){
				$nDye = @lines/$n_samples;
			}
			if($nDye==2){   # Handle 2-dye arrays if still present
				my @ch1;
				my @ch2;
				my ($n1,$n2)=(0,0);
				for(my $i=0; $i<$n_samples; $i++){
					push(@ch1,[split(/\t/,$lines[2*$i])]);
					push(@ch2,[split(/\t/,$lines[2*$i+1])]);
					if($i>0){
						if($ch1[$i]->[1] eq $ch1[$i-1]->[1]){ $n1++; }
						if($ch2[$i]->[1] eq $ch2[$i-1]->[1]){ $n2++; }
					}
				}
				if($n1>$n2){
					@ch1 = @ch2;
				}
				@lines=();
				for(my $i=0; $i<$n_samples; ++$i){
					push(@lines,"$ch1[$i]->[0]\t$ch1[$i]->[1]");
				}
			}
			my $n1 = @lines;
			if(@lines){
				$series_info{$seriesIDext} = [$title,$score,$n_samples,$platformID,\@lines];
			}
			my @keys = keys %series_info;
			$Nseries = @keys;
			if($Nseries > 100 && !$search_adjusted){
				my @scores;
				foreach my $key (@keys){
					push(@scores, [$series_info{$key}->[1],$key]);
				}
				@scores = sort {$b->[0]<=>$a->[0]} @scores;
				my $ind = int($Nseries_total/30000*100);
				if($ind>=$Nseries){ $ind=$Nseries-1; }
				$minimum_score = $scores[$ind]->[0] + 0.001;
				$search_adjusted=1;
				if(int($scores[0])==1 || $scores[0]<$minimum_score){
					for(my $i1 = 100; $i1<$Nseries; $i1++){
						delete $series_info{$scores[$i1]->[1]};
					}
					$comments = "<b>Warning:</b> Too many data series match your search criteria! Search stopped after scanning $Nseries_total series\n";
					$comments .= "To make your search more specific try using multiple search terms (comma separated), or use \"Avoid\" terms<p>\n";
					last;
				}
				foreach my $key (@keys){
					if($series_info{$key}->[1] < $minimum_score){
						delete $series_info{$key};
					}
				}
			}
		}
		$Nseries_total++;
		$line =~ s/^>//;
		($seriesID,$platformID,$title)=split(/\t/,$line);
		if($seriesID eq $seriesID_old){
			$seriesMultiple{$seriesID}=1;
		}
		$seriesID_old = $seriesID;
		$text = $line;
	}elsif($line=~/^taxid:/i){
		$taxid = $line;
		$taxid =~ s/^taxid://i;
		foreach my $id (split(/,/,$taxid)){
			if($id eq $organismID){ $taxid = $organismID; }
		}
		$text .= "\n".$line;
	}elsif($line=~/^samples:/i){
		$n_samples = $line;
		$n_samples =~ s/^samples://i;
		$text .= "\n".$line;
	}elsif($line=~/^GSM\d/){
		$text .= "\n".$line;
	}
}
my @matches;
for(my $i=0; $i<@search_term; ++$i){
	while($text =~ /$search_term[$i]/ig){
		$matches[$i]++;
	}
}
my $avoid = 0;
if($organismID != $taxid){ $avoid=1; }
for(my $i=0; $i<@avoid_term && !$avoid; ++$i){
	if($text =~ /$avoid_term[$i]/i){
		$avoid=1;
	}
}
my $found = 1;
if(@matches && !$avoid && $platform){
	$found=0;
	if($platformID eq $platform){ $found=1; }
}
my $score=0;
if(@matches && !$avoid && $found){
	my $total_match=0;
	for(my $i=0; $i<@matches; ++$i){
		if($matches[$i]){
			$score++;
			$total_match += int(log($matches[$i]));
		}
	}
	if($total_match>29){ $total_match=29; }
	$score += $total_match/30;
}
if($score >= $minimum_score){
	my @lines = split(/\n/,$text);
	while(@lines && $lines[0] !~ /^GSM\d+\t/i){ shift(@lines); }
	my $seriesIDext = $seriesID;
	my $ref = $series_info{$seriesID};
	if($ref){
		my $num = 1;
		while($series_info{"$seriesID-$num"}){ $num++; }
		$seriesIDext = "$seriesID-$num";
	}
	my $nDye = @lines/$n_samples;
	my @countRef=(0,0,0);
	for(my $i=0; $i<@lines; ++$i){
		my($sampleID,$descr,@items)=split(/\t/,$lines[$i]);
		if($nDye==2){
			for(my $i1=0; $i1<@referenceKey; ++$i1){
				if($lines[$i] =~ /$referenceKey[$i1]/i){ $countRef[$i1]++; }
			}
		}
		my $seriesID1 = $hashSamples{$sampleID};
		if($seriesID1 && $seriesID1 ne $seriesIDext){
			my $refSamples = $series_info{$seriesID1}->[4];
			for(my $j=0; $j<@$refSamples; ++$j){
				if($refSamples->[$j] =~ /^$sampleID\t/){
					splice(@$refSamples,$j,1);
					$j--;
				}
			}
			$series_info{$seriesID1}->[2]--;
		}
		$hashSamples{$sampleID} = $seriesIDext;
	}
	my $control="";
	for(my $i1=0; $i1<@referenceKey; ++$i1){
		if($countRef[$i1]==$n_samples){ $control = $referenceKey[$i1]; }
	}
	for(my $i=0; $i<@lines; ++$i){
		if($control && $lines[$i] =~ /$control/){
			splice(@lines,$i,1);
			$i--;
			next;
		}
		my($sampleID,$descr)=split(/\t/,$lines[$i]);
		$lines[$i] = "$sampleID\t$descr";
	}
	$nDye=1;
	if($n_samples && @lines){
		$nDye = @lines/$n_samples;
	}
	if($nDye==2){   # Handle 2-dye arrays if still present
		my @ch1;
		my @ch2;
		my ($n1,$n2)=(0,0);
		for(my $i=0; $i<$n_samples; $i++){
			push(@ch1,[split(/\t/,$lines[2*$i])]);
			push(@ch2,[split(/\t/,$lines[2*$i+1])]);
			if($i>0){
				if($ch1[$i]->[1] eq $ch1[$i-1]->[1]){ $n1++; }
				if($ch2[$i]->[1] eq $ch2[$i-1]->[1]){ $n2++; }
			}
		}
		if($n1>$n2){
			@ch1 = @ch2;
		}
		@lines=();
		for(my $i=0; $i<$n_samples; ++$i){
			push(@lines,"$ch1[$i]->[0]\t$ch1[$i]->[1]");
		}
	}
	if(@lines){
		$series_info{$seriesIDext} = [$title,$score,$n_samples,$platformID,\@lines];
	}
	my @keys = keys %series_info;
}
$Nseries_total++;
close INFO;

if(!%series_info){
	if($search_term !~ /GSE\d/){
		error_message("No data found in GEO that match your keywords");
	}
	for(my $i=0; $i<@search_term; $i++){
		if($search_term[$i] !~ /^GSE\d+$/){
			splice(@search_term,$i--,1);
		}
	}
	if(@search_term){
		$comments = find_specific_series(\@search_term,\%series_info);
		if(!%series_info){
			terminal_window($comments);
		}
	}
}
#Write down search results
my $fileID = get_outputID(1);
open(OUT, ">$PATH_OUTPUT/$fileID.txt") or error_message("In search_GEO1 cannot open file $fileID.txt"); 
my $pageNum=0;
my $samplesOnPage=0;
foreach my $seriesID (sort {$series_info{$b}->[1]<=>$series_info{$a}->[1]} keys %series_info){
	if($pageNum==0 || $samplesOnPage > 100){
		$pageNum++;
		$samplesOnPage=0;
	}
	my ($title,$score,$n_samples,$platformID,$refSample) = @{$series_info{$seriesID}};
	if($n_samples<=0){ next; }
	$seriesID =~ s/-\d+$//;
	if($seriesMultiple{$seriesID}){
		$seriesID .= "-$platformID";
	}
	print OUT ">$seriesID\t$title\t$score\t$n_samples\t$platformID\t$pageNum\n";
	if($series_samples eq "series"){
		$samplesOnPage++;
	}
	for(my $i=0; $i<$n_samples; ++$i){
		$refSample->[$i] =~ s/,*[\s_]+(biological[\s_]+|)(rep|replicate|replication)[\s_]*\d+$//i;
		print OUT "$refSample->[$i]\n";
		if($series_samples eq "samples"){
			$samplesOnPage++;
		}
	}
}
close OUT;
GEO_search_results($fileID,$comments);
exit(0);
}

#**************************************
sub GEO_search_results
#**************************************
{
my $fileID = shift;
my $comments = shift;

if(!$fileID){ $fileID=$hashInput{"fileID"}; }
my $organismID = $hashInput{"organismID"};
my $page_number=$hashInput{"page_number"};
my $gotopage=$hashInput{"gotopage"};
if(!$page_number){ $page_number=1; $gotopage=1; }
my $n_pages=0;

my %hashSamples;
my $file_samples;
my $description_samples;
my $series_samples = $hashInput{"series_samples"};
my $search_term = $hashInput{"search_term"};
$search_term =~ s/,/_/g;
my $savesamples = $hashInput{"savesamples"};
if($savesamples eq "save_samples"){
	$file_samples = $hashInput{"new_file"};
	$description_samples = $hashInput{"description_new_file"};
	if(!$file_samples){ error_message("In GEO search: no new sample file name"); } 
	open(OUT,">$PATH_DATA/$loginname-$file_samples");
	$hashInput{"file_samples"} = $file_samples;
	$hashInput{"description_samples"} = $description_samples;
}
elsif($savesamples eq "add_samples"){
	$file_samples = $hashInput{"file_samples"};
	open(INFO,"<$PATH_DATA/$loginname-$file_samples") or error_message("In save_sample no file $file_samples");
	while(my $line = <INFO>){
		chop $line;
		if($line =~ /^!|^#/){ next; }
		my ($seriesID,$platformID,$sampleID,$title) = split(/\t/,$line);
		$hashSamples{$sampleID}=1;
	}
	close INFO;
	open(OUT,">>$PATH_DATA/$loginname-$file_samples");
}

#Read search results
my @series_info;
open(INFO,"<$PATH_OUTPUT/$fileID.txt") or error_message("In GEO_search_results no file $fileID.txt");
my @samples;
my $changed=0;
my ($seriesID,$title,$score,$n_samples,$platformID,$pageNum,$seriesChecked);
my $nSamples = 0;
my $n_checked = 0;
my $n_added = 0;
while(my $line = <INFO>){
	chop $line;
	$line =~ s/[^[:ascii:]]//g;
	if($line =~ /^>/){
		if($seriesID){
			my @samples1 = @samples;
			push(@series_info,[$seriesID,$title,$score,$n_samples,$platformID,$pageNum,\@samples1,$seriesChecked]);
		}
		$line =~ s/^>//;
		($seriesID,$title,$score,$n_samples,$platformID,$pageNum,$seriesChecked) = split(/\t/,$line);
		if($series_samples eq "series"){
			my $checked1 = $hashInput{$seriesID};
			if($checked1 eq "on" && !$seriesChecked){ $seriesChecked=1; $changed=1; }
			elsif(!$checked1 && $page_number==$pageNum && $seriesChecked==1){ $seriesChecked=0; $changed=1; }
		}
		$n_pages = $pageNum;
		@samples=();
	}else{
		$nSamples++;
		my ($sampleID,$sampleTitle,$checked) = split(/\t/,$line);
		if($series_samples eq "samples"){
			my $checked1 = $hashInput{$sampleID};
			if($checked1 eq "on" && !$checked){ $checked=1; $changed=1; }
			elsif(!$checked1 && $page_number==$pageNum && $checked==1){ $checked=0; $changed=1; }
		}
		push(@samples,[$sampleID,$sampleTitle,$checked]);
		if($savesamples && ($checked || $series_samples eq "series" && $seriesChecked)){
			if($savesamples eq "add_samples" && $hashSamples{$sampleID}){ next; } #Prevent duplication of samples
			if($savesamples eq "save_samples"){ $hashInput{$sampleID}=""; }
			print OUT "$seriesID\t$platformID\t$sampleID\t$sampleTitle\n";
			$n_added++;
		}
	}
}
close OUT;
if($savesamples eq "add_samples"){
	if(!$n_added){ error_message("Samples are already present in the destination file"); }
	terminal_window("<H3>Samples saved in '$file_samples' (N = $n_added).</H3>");
}
push(@series_info,[$seriesID,$title,$score,$n_samples,$platformID,$pageNum,\@samples,$seriesChecked]);
close INFO;
my $nSeries = @series_info;
# Replace file if changed
if($changed){
	open(OUT1, ">$PATH_OUTPUT/$fileID.txt") or error_message("In search_GEO1 cannot open file $fileID.txt"); 
	foreach my $ref (@series_info){
		my ($seriesID,$title,$score,$n_samples,$platformID,$pageNum,$refSample,$seriesChecked) = @$ref;
		print OUT1 ">$seriesID\t$title\t$score\t$n_samples\t$platformID\t$pageNum\t$seriesChecked\n";
		foreach my $ref1 (@$refSample){
			if($seriesChecked){ $ref1->[2]=1; }
			if($ref1->[2]){ $n_checked++; }
			print OUT1 join("\t",@$ref1)."\n";
		}
	}
	close OUT1;
}
if($savesamples eq "save_samples"){
	#Update configuration file
	open(OUT, ">$PATH_INFO/$loginname-config1.txt");
	open(INFO,"<$PATH_INFO/$loginname-config.txt");
	my $nLines;
	while(my $line = <INFO>){
		$nLines++;
		if($line !~ /^type_samples=$file_samples\s/ && length($line)>2){
			print OUT $line;
		}
	}
	close INFO;
	print OUT "type_samples=$file_samples\torganismID=$organismID";
	if($description_samples){ print OUT "\tdescription=$description_samples"; }
	print OUT "\tdate=$date_record\n";
	close OUT;
	my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
	if($nLines1 && $nLines1 > $nLines*0.9){
		copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
	}else{
		error_message("Failed to update configuration file!");
	}
	unlink("$PATH_INFO/$loginname-config1.txt");
	$hashInput{"page_number"}=0;
	open_samples();
}
$page_number = $gotopage;
my @samples_list;
my @file_list;
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Config file not found!");
while(my $line = <INFO>){
	chop $line;
	my %hash=();
	if(length($line) < 3) { next; }
	my @items = split(/[=\t]/,$line);
	read_config_line($line,\%hash);
	my $file_samples = $hash{"type_samples"};
	if($file_samples){
		my $descr = $hash{"description"};
		my $org = $hash{"organismID"};
		push(@samples_list,[$file_samples,$descr,$org]);
	}elsif($items[0] =~ /^type_/){
		push(@file_list,$items[1]);
	}
}
close INFO;
filter_list_by_organism(\@samples_list, $organismID);
@samples_list = sort {lc($a->[0]) cmp lc($b->[0])} @samples_list;
my ($samples_list,$sample_descriptions) = get_array_lists(\@samples_list);
if($sample_descriptions){ $sample_descriptions = "\"\",".$sample_descriptions; }
else{ $sample_descriptions = "\"\""; }
my $file_list="";
foreach my $name (@file_list){
	if(!$file_list){ $file_list = "\"".$name."\""; }
	else{ $file_list .= ",\"".$name."\""; }
}
print "<HTML><HEAD><TITLE>ExAtlas - GEO search</TITLE>\n";
print_header("update_description();");
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "samples_list = new Array($samples_list);\n";
print "samples_description = new Array($sample_descriptions);\n";
print "file_list = new Array($file_list);\n";
print "function goto_page(page){\n";
print "	clean_fields();\n";
print "	document.GEO_search.action.value=\"GEO_search_results\";\n";
print "	if(page<0){\n";
print "		page=document.GEO_search.select_page.options[document.GEO_search.select_page.selectedIndex].value;\n";
print "	}\n";
print "	document.GEO_search.gotopage.value=page;\n";
print "	document.GEO_search.submit();\n";
print "}\n";
print "function clean_fields(){\n";
print "	document.GEO_search.target = \"\";\n";
print "	document.GEO_search.action.value=\"GEO_search_results\";\n";
print "	document.GEO_search.savesamples.value=\"\";\n";
print "}\n";
print "function update_description() {\n";
print "	clean_fields();\n";
print " var index;\n";
print " index = document.GEO_search.file_samples.selectedIndex;\n";
print " document.GEO_search.description_samples.value = samples_description[index];\n";
print "}\n";
print "function add_samples() {\n";
print "	if(document.GEO_search.file_samples.selectedIndex==0){\n";
print "		alert(\"Select existing sample file to append.\"); return false;\n";
print "	}\n";
print "	clean_fields();\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.GEO_search.target = \"_BLANK\"+x;\n";
print "	document.GEO_search.savesamples.value = \"add_samples\";\n";
print "	document.GEO_search.submit();\n";
print "}\n";
print "function save_samples() {\n";
print "	var file = document.GEO_search.new_file.value;\n";
print "	if(!file){\n";
print "		alert(\"Enter new file name (and description) for samples\"); return false;\n";
print "	}\n";
print "	var file1=file;\n";
print "	if(file.search(/\\.txt\$/)>=0){\n";
print "		if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "		file1=file.substring(0,file.length-4);\n";
print "	}\n";
print "	if(file1.search(/^[-\\w]+\$/)<0){\n";
print "		alert(\"File name should be one word with no special characters except underscore and dash.\\nRename it in provided field\");\n";
print "		return(false);\n";
print "	}\n";
print "	for(i=0; i<samples_list.length; ++i){\n";
print "		if(file == samples_list[i] && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	for(i=0; i<file_list.length; ++i){\n";
print "		if(file == file_list[i]){\n";
print "			alert(\"A file with this name already exists\"); return false;\n";
print "		}\n";
print "	}\n";
print "	if(file.search(/^public-/i) >= 0){\n";
print "		alert(\"File name cannot start with 'public-'\"); return(false);\n";
print "	}\n";
print "	var descrip = document.GEO_search.description_new_file.value;\n";
print "	if(descrip.search(/\\=|\\&/) >= 0){\n";
print "		alert(\"Description should not include character \'=\' or \'&\'\");\n";
print "		return false;\n";
print "	}\n";
print "	clean_fields();\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.GEO_search.target = \"_BLANK\"+x;\n";
print "	document.GEO_search.gotopage.value = 0;\n";
print "	document.GEO_search.savesamples.value = \"save_samples\";\n";
print "	document.GEO_search.submit();\n";
print "}\n";
print "function check_all(){\n";
print "	for(i=0; i<document.GEO_search.elements.length; ++i){\n";
print "		if(document.GEO_search.elements[i].type==\"checkbox\"){\n";
print "			document.GEO_search.elements[i].checked=true;\n";
print "		}\n";
print "	}\n";
print "}\n";
print "function uncheck_all(){\n";
print "	for(i=0; i<document.GEO_search.elements.length; ++i){\n";
print "		if(document.GEO_search.elements[i].type==\"checkbox\"){\n";
print "			document.GEO_search.elements[i].checked=false;\n";
print "		}\n";
print "	}\n";
print "}\n";
print "<!-- end script --></SCRIPT></HEAD>\n";

print "<FORM NAME=GEO_search ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<H2>Results of GEO database search ($nSeries data series)</H2>\n";
print $comments;
if($savesamples eq "add_samples"){
	print "<b>Note:</b> Selected samples (N = $n_added) were added to file '$file_samples'<p>\n";
}
print "<b>Organism:</b> $hashOrganism{$organismID} &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;\n";
print "<u>Text file with search results:</u> <a href=$HOME_ADDRESS/output/$fileID.txt target=_blank329>$fileID.txt</a>\n";
print "<TABLE border=0>\n";
print "<TR><TD>Page $page_number. (N=$n_pages)<TD>\n";
if($page_number>1){
	print "<INPUT TYPE=button VALUE=\"Previous page\" LANGUAGE=\"javascript\" onClick=goto_page($page_number-1);>\n";
}
if($page_number<$n_pages){
	print "<INPUT TYPE=button VALUE=\"Next page\" LANGUAGE=\"javascript\" onClick=goto_page($page_number+1);>\n";
}
if($n_pages>2){
	print "Go to page<SELECT NAME=select_page LANGUAGE=\"javascript\" onChange=goto_page(-1);>\n";
	for(my $i=1; $i<=$n_pages; $i++){
		if($page_number==$i){ print "<option value=$i selected> $i\n"; }
		else{ print "<option value=$i> $i\n"; }
	}
	print "</SELECT>\n";
}
print "<TD><center><b>File name<TD><center><b>Description\n";
print "<TR><TD>N series = $nSeries<TD><center>Add samples to existing file:<TD><select name=\"file_samples\" onChange=update_description(); style=width:200px;>\n";
print "<option value=\"\"> ------------select------------\n";
for(my $i=0; $i<@samples_list; ++$i){ 
	print "<option value=\"$samples_list[$i]->[0]\"> $samples_list[$i]->[0]\n";
}
print "</select>\n";
my $new_file_name1 = $search_term;
$new_file_name1 =~ s/[\s\W]/_/g;
print "<TD><INPUT NAME=\"description_samples\" SIZE=30 VALUE=\"\">\n";
print "<TD><INPUT TYPE=button VALUE=\"Add samples\" LANGUAGE=\"javascript\" onClick=add_samples(); style=width:100px;>\n";
print "<TR><TD>N samples = $nSamples<TD><center>Or create a new samples file:\n";
print "<TD><INPUT NAME=new_file VALUE=\"$new_file_name1\" style=width:200px;><TD><INPUT NAME=\"description_new_file\" SIZE=30 VALUE=\"$search_term\">\n";
print "<TD><INPUT TYPE=button VALUE=\"Save samples\" LANGUAGE=\"javascript\" onClick=save_samples(); style=width:100px;>\n";
print "</TABLE>\n";

print "<H3>Table of GEO data sets and samples</H3>\n";
print "Click on the links to data sets or samples to see details &nbsp; &nbsp; &nbsp; &nbsp; \n";
print "<INPUT TYPE=button VALUE=\"Uncheck all samples\" LANGUAGE=\"javascript\" onClick=uncheck_all();>\n";
print "<INPUT TYPE=button VALUE=\"Check all samples\" LANGUAGE=\"javascript\" onClick=check_all();>\n";

my $count=0;
print "<TABLE border=0>\n";
foreach my $ref (@series_info){
	my ($seriesID,$title,$score,$n_samples,$platformID,$page,$refSample,$seriesChecked) = @$ref;
	$count++;
	if($page < $page_number){ next; }
	if($page > $page_number){ last; }
	my $seriesID1=$seriesID;
	$seriesID1 =~ s/-.+$//;
	print "<TR><TD>$count.\n";
	if($series_samples eq "series"){
		print "<TD><INPUT TYPE=checkbox NAME=$seriesID";
		if($seriesChecked){ print " checked"; }
		print ">\n";
	}
	my $x = int(10000*rand());
	print "<TD><b><a href=\"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$seriesID1\" target=_BLANK$x>$seriesID1</a>\n";
	$x = int(10000*rand());
	print "<TD><a href=\"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$platformID\" target=_BLANK$x>$platformID</a>\n";
	print "<TD>$n_samples<TD>$title\n";
	if($series_samples eq "samples"){
		foreach my $ref1 (@$refSample){
			my ($sampleID,$sampleTitle,$checked) = @$ref1;
			print "<TR><TD><TD><INPUT TYPE=checkbox NAME=$sampleID";
			if($checked){ print " checked"; }
			print ">\n";
			$x = int(10000*rand());
			print "<TD><a href=\"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$sampleID\" target=_BLANK$x>$sampleID</a><TD COLSPAN=4>$sampleTitle\n";
		}
	}
}
print "</TABLE>\n";

print "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
print "<INPUT NAME=action TYPE=hidden VALUE=GEO_search_results>\n";
print "<INPUT NAME=gotopage TYPE=hidden VALUE=$page_number>\n";
print "<INPUT NAME=page_number TYPE=hidden VALUE=$page_number>\n";
print "<INPUT NAME=fileID TYPE=hidden VALUE=$fileID>\n";
print "<INPUT NAME=series_samples TYPE=hidden VALUE=$series_samples>\n";
print "<INPUT NAME=savesamples TYPE=hidden>\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT NAME=search_term TYPE=hidden VALUE=$search_term>\n";
print "</FORM>\n";
print "<p><HR NOSHADE></HR>\n";
print "<INPUT TYPE=button VALUE=\" Cancel (close window) \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "<p><i>Please report any problems to <a href=mailto:sharoval\@mail.nih.gov>webmaster</a><p>\n";
print "</BODY>\n";
print "</HTML>\n";
}

#**************************************
sub  open_samples 
#**************************************
{
my %hashSamples;
my $file_samples = $hashInput{"file_samples"};
my $description_samples = $hashInput{"description_samples"};
my $organismID = $hashInput{"organismID"};
my $page_number=$hashInput{"page_number"};
my $gotopage=$hashInput{"gotopage"};
my $deletesamples=$hashInput{"deletesamples"};
my $copysamples=$hashInput{"copysamples"};
my $file_tempID=$hashInput{"file_tempID"};
my $LINES_PER_PAGE = 50;

my $error_notes;
if($file_tempID){
	open(INFO,"<$PATH_OUTPUT/$file_tempID.txt") or error_message("In open_samples no temp fileID");
}else{
	open(INFO,"<$PATH_DATA/$loginname-$file_samples") or error_message("In open_samples no file_samples");
}
my @series_info;
my @copy_info;
my %hashSeries;
my $count=0;
my $n_checked=0;
my $changed=0;
while(my $line = <INFO>){
	chop $line;
	if($line =~ /^!|^#/){ next; }
	my ($seriesID,$platformID,$sampleID,$title,$selected) = split(/\t/,$line);
	if(!$page_number){ $selected=0; }
	$hashSamples{$sampleID}=1;
	my $pageNum = int($count/$LINES_PER_PAGE)+1;
	if($page_number==$pageNum && $selected && $hashInput{$sampleID} ne "on"){ $selected = 0; }
	if(!$selected && $hashInput{$sampleID} eq "on"){ $selected = 1; }
	if($deletesamples && $selected){ $changed=1; next; }
	if($copysamples && $selected){
		push(@copy_info,[$seriesID,$platformID,$sampleID,$title]);
		$selected=0;
	}
	my $title1 = $hashInput{"$sampleID-title"}; #Rename sample if specified
	if($title1=~/\w/ && $title ne $title1){
		$title1=~ s/^\s+//; $title1=~ s/\s+$//;
		$title = $title1;
		$changed=1;
	}
	$n_checked += $selected;
	my $ii = $hashSeries{$seriesID};
	if(defined $ii){
		push(@{$series_info[$ii]},[$seriesID,$platformID,$sampleID,$title,$selected]);
		if($ii < @series_info-1){
			$changed=1;
		}
	}else{
		push(@series_info,[[$seriesID,$platformID,$sampleID,$title,$selected]]);
		$hashSeries{$seriesID} = @series_info-1;
	}
	$count++;
}
close INFO;
if(!$page_number){ $page_number=1; $gotopage=1; }

if($copysamples && !@copy_info){
	$error_notes = "No samples were selected for copying. ";
	$copysamples="";
}
my %hashSamplesCopy;
my $file_samples_copy = $hashInput{"file_samples_copy"};
my $add_to_config = 0;
if($file_samples_copy eq "new_file" && $copysamples){
	$file_samples_copy = $copysamples;
	$add_to_config = 1;
}
my $count_samples_copied=0;
if($copysamples && open(INFO,"<$PATH_DATA/$loginname-$file_samples_copy")){
	while(my $line = <INFO>){
		if($line =~ /^!|^#/){ next; }
		my ($seriesID,$platformID,$sampleID,$title) = split(/\t/,$line);
		$hashSamplesCopy{$sampleID}=1;
	}
	close INFO;
}
if($copysamples){
	my $text_samples;
	foreach my $ref (@copy_info){
		my($seriesID,$platformID,$sampleID,$title) = @$ref;
		if(!$hashSamplesCopy{$sampleID}){
			$text_samples .= "$seriesID\t$platformID\t$sampleID\t$title\n";
			$count_samples_copied++;
		}
	}
	if($count_samples_copied){
		if(open(OUT,">>$PATH_DATA/$loginname-$file_samples_copy")){
			print OUT $text_samples;
			close OUT;
		}else{
			$count_samples_copied=0;
			$error_notes = "Failure writing samples to file. ";
		}
	}else{
		$error_notes = "Samples already exist in the destination file! ";
	}
}
#Reorder and update samples file if necessary
if(!$file_tempID){
	$file_tempID = get_outputID(1);
}
open(OUT,">$PATH_OUTPUT/$file_tempID.txt");
$count=0;
for(my $i=0; $i<@series_info; ++$i){
	my $ref = $series_info[$i];
	my %hashDone=();
	my @reorder=();
	for(my $j=0; $j<@$ref; ++$j){
		if($hashDone{$j}){ next; }
		my $title = $ref->[$j]->[3];
		push(@reorder,$ref->[$j]);
		my $gaps=0;
		for(my $j1=$j+1; $j1<@$ref; ++$j1){
			my $title1 = $ref->[$j1]->[3];
			if($title1 eq $title){
				push(@reorder,$ref->[$j1]);
				$hashDone{$j1} = 1;
				if($gaps){ $changed=1; }
			}else{
				$gaps=1;
			}
		}
	}
	$count += @reorder;
	$series_info[$i] = \@reorder;
	foreach my $ref1 (@reorder){
		print OUT join("\t",@$ref1)."\n";
	}
}
close OUT;
if($changed){
	open(OUT,">$PATH_DATA/$loginname-$file_samples");
	foreach my $ref (@series_info){
		for(my $j=0; $j<@$ref; ++$j){
			my($seriesID,$platformID,$sampleID,$title) = @{$ref->[$j]};
			print OUT "$seriesID\t$platformID\t$sampleID\t$title\n";
		}
	}
	close OUT;
}
if($hashInput{"generatematrix"}){
	my $logFileID = get_outputID(1);
	$hashInput{"logFileID"} = $logFileID;
	generate_matrix(\@series_info);
	exit(0);
}
my $n_pages = int($count/$LINES_PER_PAGE)+1;
$page_number = $gotopage;
my @file_list;
my @matrix_list;
my @samples_list;
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Config file not found!");
while(my $line = <INFO>){
	chop $line;
	if(length($line) < 3) { next; }
	my @items = split(/[=\t]/,$line);
	if($items[0] =~ /^type_matrix/){
		push(@matrix_list,$items[1]);
	}elsif($items[0] =~ /^type_/){
		push(@file_list,$items[1]);
	}
	my %hash=();
	read_config_line($line,\%hash);
	my $file_samples1 = $hash{"type_samples"};
	if($file_samples1 && $file_samples1 ne $file_samples){
		my $descr = $hash{"description"};
		my $org = $hash{"organismID"};
		push(@samples_list,[$file_samples1,$descr,$org]);
	}
}
close INFO;
filter_list_by_organism(\@samples_list, $organismID);
@samples_list = sort {lc($a->[0]) cmp lc($b->[0])} @samples_list;
if($count_samples_copied && $add_to_config){
	open(OUT,">>$PATH_INFO/$loginname-config.txt") or error_message("Cannot open file config");
	my $description_samples_copy = $hashInput{"description_samples_copy"};
	print OUT "type_samples=$file_samples_copy\torganismID=$organismID";
	if($description_samples_copy){ print OUT "\tdescription=$description_samples_copy"; }
	print OUT "\tdate=$date_record\n";
	close OUT;
}
my ($samples_list,$sample_descriptions) = get_array_lists(\@samples_list);
if($sample_descriptions){ $sample_descriptions = "\"\",".$sample_descriptions; }
else{ $sample_descriptions = "\"\""; }
my $file_list="";
my $matrix_list="";
foreach my $name (@matrix_list){
	if(!$matrix_list){ $matrix_list = "\"".$name."\""; }
	else{ $matrix_list .= ",\"".$name."\""; }
}
foreach my $name (@file_list){
	if(!$file_list){ $file_list = "\"".$name."\""; }
	else{ $file_list .= ",\"".$name."\""; }
}
my $file_matrix = $hashInput{"file_new_matrix"};
my $description_matrix = $hashInput{"description_new_matrix"};
if(!$file_matrix){
	$file_matrix = $file_samples;
	$file_matrix =~ s/\..+$//;
	$file_matrix .= "-matrix";
}
if(!$description_matrix){ $description_matrix = $description_samples; }

print "<HTML><HEAD><TITLE>ExAtlas - open_samples</TITLE>\n";
print_header("update_description();");
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "file_list = new Array($file_list);\n";
print "matrix_list = new Array($matrix_list);\n";
print "samples_description = new Array($sample_descriptions);\n";
print "function count_checked_samples(){\n";
print "	var nchecked = $n_checked;\n";
print "	for(i=0; i<document.open_samples.elements.length; ++i){\n";
print "		if(document.open_samples.elements[i].type==\"checkbox\"){\n";
print "			if(document.open_samples.elements[i].checked){ nchecked++; };\n";
print "		}\n";
print "	}\n";
print "	return nchecked;\n";
print "}\n";
print "function update_description(){\n";
print " var index;\n";
print " index = document.open_samples.file_samples_copy.selectedIndex;\n";
print " document.open_samples.description_samples_copy.value = samples_description[index];\n";
print "}\n";
print "function delete_samples(){\n";
print "	if(count_checked_samples()==0){\n";
print "		alert(\"No samples selected! Cannot proceed.\"); return false;\n";
print "	}\n";
print "	clean_fields();\n";
print "	document.open_samples.deletesamples.value=\"delete_samples\";\n";
print "	document.open_samples.submit();\n";
print "}\n";
print "function copy_samples(){\n";
print "	if(count_checked_samples()==0){\n";
print "		alert(\"No samples selected! Cannot proceed.\"); return false;\n";
print "	}\n";
print "	clean_fields();\n";
print "	if(document.open_samples.file_samples_copy.selectedIndex==0){\n";
print " 	var file;\n";
print "		file = prompt(\"Provide name of a new file where to copy selected samples\");\n";
print "		if(!file){ return(false); }\n";
print "		var file1=file;\n";
print "		if(file.search(/\\.txt\$/)>=0){\n";
print "			if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "			file1=file.substring(0,file.length-4);\n";
print "		}\n";
print "		if(file1.search(/^[-\\w]+\$/)<0){\n";
print "			alert(\"File name should have neither spaces nor special characters\");\n";
print "			return(false);\n";
print "		}\n";
print "		if(file.search(/^public-/i) >= 0){\n";
print "			alert(\"File name cannot start with 'public-'\"); return(false);\n";
print "		}\n";
print "		var descrip = document.open_samples.description_samples_copy.value;\n";
print "		if(descrip.search(/\\=|\\&/) >= 0){\n";
print "			alert(\"Description should not include character \'=\' or \'&\'\");\n";
print "			return false;\n";
print "		}\n";
print "		for(i=0; i<file_list.length; ++i){\n";
print "			if(file == file_list[i] || file==\"$file_samples\"){\n";
print "				alert(\"File with this name already exists\"); return false;\n";
print "			}\n";
print "		}\n";
print "		document.open_samples.copysamples.value=file;\n";
print "	}else{\n";
print "		document.open_samples.copysamples.value=\"copy_samples\";\n";
print "	}\n";
print "	document.open_samples.submit();\n";
print "}\n";
print "function clean_fields(){\n";
print "	document.open_samples.deletesamples.value=\"\";\n";
print "	document.open_samples.copysamples.value=\"\";\n";
print "	document.open_samples.generatematrix.value=\"\";\n";
print "	document.open_samples.target=\"\";\n";
print "}\n";
print "function goto_page(page){\n";
print "	clean_fields();\n";
print "	if(page<0){\n";
print "		page=document.open_samples.select_page.options[document.open_samples.select_page.selectedIndex].value;\n";
print "	}\n";
print "	document.open_samples.gotopage.value=page;\n";
print "	document.open_samples.submit();\n";
print "}\n";
print "function save_matrix() {\n";
print "	clean_fields();\n";
print "	var file = document.open_samples.file_new_matrix.value;\n";
print "	if(!file){\n";
print "		alert(\"You need to enter new file name.\"); return false;\n";
print "	}\n";
print "	var file1=file;\n";
print "	if(file.search(/\\.txt\$/)>=0){\n";
print "		if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "		file1=file.substring(0,file.length-4);\n";
print "	}\n";
print "	if(file1.search(/^[-\\w]+\$/)<0){\n";
print "		alert(\"File name should be one word with no special characters except underscore and dash.\\nRename it in provided field\");\n";
print "		return(false);\n";
print "	}\n";
print "	var descrip = document.open_samples.description_new_matrix.value;\n";
print "	if(descrip.search(/\\=|\\&/) >= 0){\n";
print "		alert(\"Description should not include character \'=\' or \'&\'\");\n";
print "		return false;\n";
print "	}\n";
print "	if(file.search(/^public/i) >= 0){\n";
print "		alert(\"File name cannot start with 'public'\"); return(false);\n";
print "	}\n";
print "	for(i=0; i<matrix_list.length; ++i){\n";
print "		if(file == matrix_list[i] && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	for(i=0; i<file_list.length; ++i){\n";
print "		if(file == file_list[i]){\n";
print "			alert(\"File with this name already exists\"); return false;\n";
print "		}\n";
print "	}\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.open_samples.target = \"_BLANK\"+x;\n";
print "	document.open_samples.generatematrix.value = \"generate_matrix\";\n";
print "	document.open_samples.submit();\n";
print "}\n";
print "function check_all(){\n";
print "	for(i=0; i<document.open_samples.elements.length; ++i){\n";
print "		if(document.open_samples.elements[i].type==\"checkbox\"){\n";
print "			document.open_samples.elements[i].checked=true;\n";
print "		}\n";
print "	}\n";
print "}\n";
print "function uncheck_all(){\n";
print "	for(i=0; i<document.open_samples.elements.length; ++i){\n";
print "		if(document.open_samples.elements[i].type==\"checkbox\"){\n";
print "			document.open_samples.elements[i].checked=false;\n";
print "		}\n";
print "	}\n";
print "}\n";
print "<!-- end script --></SCRIPT></HEAD>\n";
print "<FORM NAME=open_samples ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<H2>Samples file '$file_samples'</H2>\n";
if($n_pages>1){
	print "<b>Page $page_number. (N=$n_pages)</b>\n";
}
if($page_number>1){
	print "<INPUT TYPE=button VALUE=\"Previous page\" LANGUAGE=\"javascript\" onClick=goto_page($page_number-1);>\n";
}
if($page_number<$n_pages){
	print "<INPUT TYPE=button VALUE=\"Next page\" LANGUAGE=\"javascript\" onClick=goto_page($page_number+1);>\n";
}
if($n_pages>2){
	print "Go to page<SELECT NAME=select_page LANGUAGE=\"javascript\" onChange=goto_page(-1);>\n";
	for(my $i=1; $i<=$n_pages; $i++){
		if($page_number==$i){ print "<option value=$i selected> $i\n"; }
		else{ print "<option value=$i> $i\n"; }
	}
	print "</SELECT>\n";
}
print "<p>\n";
if($copysamples){
	if($count_samples_copied){
		print "<b>Samples copied:</b> N = $count_samples_copied to file $file_samples_copy\n";
	}else{
		print "<b>No samples copied:</b> they are already in the target file\n";
	}
	print "<p>\n";
}
if($error_notes){
	print "<b>Errors:</b> $error_notes<p>\n";
}
if($description_samples){
	print "<b>Description:</b> $description_samples<br>\n";
}
print "<b>Organism:</b> $hashOrganism{$organismID}<p>\n";
print "<p>You can delete selected samples or edit their descriptions. Samples with identical descriptions in the same data set (series) are<br>\n";
print "considered as replications and will be placed together. If descriptions are cryptic, click on sample ID and<br>\n";
print "check descriptions at GEO database. After the list is finalized, generate a combined matrix of all samples.<p>\n";
print "<TABLE border=0>\n";
print "<TR><TD>Edit buttons<TD><TD><TD><center>File name<TD><center>Description<TD><center>Action buttons\n";
print "<TR><TD><INPUT TYPE=button VALUE=\"Refresh list\" LANGUAGE=\"javascript\" onClick=goto_page($page_number); style=width:200px;>\n";
print "<TD WIDTH=50><TD>Provide matrix name:\n";
print "<TD><INPUT NAME=file_new_matrix VALUE=\"$file_matrix\" style=width:150px;>\n";
print "<TD><INPUT NAME=description_new_matrix VALUE=\"$description_matrix\" style=width:150px;>\n";
print "<TD><INPUT TYPE=button VALUE=\"Generate matrix\" LANGUAGE=\"javascript\" onClick=save_matrix(); style=width:200px;>\n";
print "<TR><TD><INPUT TYPE=button VALUE=\"Delete selected samples\" LANGUAGE=\"javascript\" onClick=delete_samples(); style=width:200px;>\n";
print "<TD><TD>Copy selected sample to:\n";
print "<TD><select name=file_samples_copy onChange=update_description(); style=width:150px;>\n";
print "<option value=new_file> -------- New file -------\n";
for(my $i=0; $i<@samples_list; ++$i){ 
	print "<option value=\"$samples_list[$i]->[0]\"> $samples_list[$i]->[0]\n";
}
print "</select>\n";
print "<TD><INPUT NAME=description_samples_copy style=width:150px;>\n";
print "<TD><INPUT TYPE=button VALUE=\"Copy selected samples\" LANGUAGE=\"javascript\" onClick=copy_samples(); style=width:200px;>\n";
print "</TABLE><br>\n";
print "If you select \"New file\" you will be prompted for file name, which shound be one-word without special characters (underscore allowed)<br>\n";
print "<INPUT TYPE=button VALUE=\"Uncheck all samples\" LANGUAGE=\"javascript\" onClick=uncheck_all();>\n";
print "<INPUT TYPE=button VALUE=\"Check all samples\" LANGUAGE=\"javascript\" onClick=check_all();>\n";

print "<TABLE border=0>\n";
print "<TR><TD>No.<TD><TD><TD>Sample ID<TD WIDTH=10><TD>Series ID<TD WIDTH=10><TD>Platform ID<TD WIDTH=10><TD>Sample description (edit if needed)\n";
$count = 0;
foreach my $ref (@series_info){
	foreach my $ref1 (@$ref){
		my $page = int($count/$LINES_PER_PAGE)+1;
		$count++;
		if($page < $page_number){ next; }
		if($page > $page_number){ last; }
		my ($seriesID,$platformID,$sampleID,$title,$selected) = @$ref1;
		my $seriesID1=$seriesID;
		$seriesID1 =~ s/-.+$//;
		print "<TR><TD>$count.";
		print "<TD><INPUT TYPE=checkbox NAME=$sampleID";
		if($selected){ print " checked"; }
		print ">\n";
		my $x = int(10000*rand());
		print "<TD><TD><a href=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$sampleID target=_BLANK$x>$sampleID</a>\n";
		$x++;
		print "<TD><TD><a href=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$seriesID1 target=_BLANK$x>$seriesID1</a>\n";
		$x++;
		print "<TD><TD><a href=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$platformID target=_BLANK$x>$platformID</a>\n";
		print "<TD><TD WIDTH=500><INPUT NAME=$sampleID-title style=\"width: 100%;\" VALUE=\"$title\">\n";
	}
}
print "</TABLE><p>\n";

print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=action TYPE=hidden VALUE=open_samples>\n";
print "<INPUT NAME=file_samples TYPE=hidden VALUE=\"$file_samples\">\n";
print "<INPUT NAME=description_samples TYPE=hidden VALUE=\"$description_samples\">\n";
print "<INPUT NAME=generatematrix TYPE=hidden>\n";
print "<INPUT NAME=deletesamples TYPE=hidden>\n";
print "<INPUT NAME=copysamples TYPE=hidden>\n";
print "<INPUT NAME=gotopage TYPE=hidden VALUE=$page_number>\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT NAME=page_number TYPE=hidden VALUE=$page_number>\n";
print "<INPUT NAME=file_tempID TYPE=hidden VALUE=$file_tempID>\n";
print "</FORM><p>\n";

print "<p><HR NOSHADE></HR>\n";
print "<INPUT TYPE=button VALUE=\"     Close window     \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "<p><i>Please report any problems to <a href=mailto:sharoval\@mail.nih.gov>webmaster</a><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub    generate_matrix 
#**************************************
{
my $series_info = shift;
my $organismID = $hashInput{"organismID"};
my $logFileID = $hashInput{"logFileID"};
my $file_matrix = $hashInput{"file_new_matrix"};
my $description_matrix = $hashInput{"description_new_matrix"};
if($logFileID && open(INFO,"<$PATH_OUTPUT/$logFileID.txt")){
	my $response="";
	my %hash;
	while(my $line=<INFO>){
		$response .= $line;
		if($line =~ /\t/){
			chop $line;
			my @items = split(/\t/,$line);
			$hash{$items[0]}=$items[1];
		}
	}
	close INFO;
	if($response =~ /error/i){
		print "<HTML><HEAD><TITLE>ExAtlas - generate_matrix</TITLE>\n";
		print_header();
		print "<h3>Errors in processing, task stopped</h3>\n";
		print "<b>Log information</b><br><pre>\n";
		print "$response\n</pre>";
		print "<INPUT TYPE=button VALUE=\"  Close window  \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
		exit(0);
	}
	if($response =~ /Task completed/i){
		$hashInput{"file_matrix"} = $file_matrix;
		$hashInput{"description_matrix"} = $description_matrix;
		matrix_explore();
		exit(0);
	}
	my $pid = $hashInput{"process_id"};
	if(!$pid){ $pid = $hash{"process_id"}; }
	print "<HTML><HEAD><TITLE>ExAtlas - generate_matrix</TITLE>\n";
	print_header();
	print "<h3>Your task is not finished!</h3>\n";
	print "Check your task later ...<br>\n";
	print "Alternatively you can close this window and later access results from main menu: click the \"Refresh\" button, select file '$file_matrix' in the pull-down menu \"matrix files\"\n";
	print "and click the button \"Explore matrix\".<p>\n";
	my $x = int(10000*rand());
	print "The status of your task can be checked here: <a href=$HOME_ADDRESS/output/$logFileID.txt target=_blank$x>Log file</a><p>\n";
	print "<FORM ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
	print "<TABLE BORDER=0>\n";
	print "<TR><TD><INPUT TYPE=submit VALUE=\"  Check your task  \">\n";
	print "<TD WIDTH=100><TD><INPUT NAME=terminate_task TYPE=submit VALUE=\"  Cancel the task  \">\n";
	print "<TD WIDTH=100><TD><INPUT TYPE=button VALUE=\"  Close window  \" LANGUAGE=\"javascript\" onClick=\"window.close();\">\n";
	print "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
	print "<INPUT NAME=\"process_id\" TYPE=hidden VALUE=\"$pid\">\n";
	print "<INPUT NAME=logFileID TYPE=hidden VALUE=\"$logFileID\">\n";
	print "<INPUT NAME=file_new_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
	print "<INPUT NAME=description_new_matrix TYPE=hidden VALUE=\"$description_matrix\">\n";
	print "<INPUT NAME=action TYPE=hidden VALUE=\"generate_matrix\">\n";
	print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
	print "</TABLE></FORM><p><HR NOSHADE>\n";
	print "</BODY>\n";
	print "</HTML>\n";
	exit(0);
}
my $file_description = $hashInput{"description_new_matrix"};
if(!$file_matrix){ error_message("Blank name for matrix file"); } 
if(!$organismID){ error_message("No organism ID"); } 
if($logFileID){ file_append("Task started","$PATH_OUTPUT/$logFileID.txt",1); }

my $pid = fork();
if($pid){
	$|++;
	$pid = pid_record($pid);
	print "<HTML><HEAD><TITLE>ExAtlas - generate_matrix</TITLE>\n";
	print_header();
	print "<h2>Your task is submitted</h2>\n";
	print "Expression profiles will be saved in file \"$file_matrix\"\n";
	print "It may take some time to generate the file. To see the results, click on \"Check your task\" button<p>\n";
	my $x = int(10000*rand());
	print "The status of your task can be checked here: <a href=$HOME_ADDRESS/output/$logFileID.txt target=_blank$x>Log file</a><p>\n";
	print "<FORM ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
	print "<TABLE BORDER=0>\n";
	print "<TR><TD><INPUT TYPE=submit VALUE=\"  Check your task  \">\n";
	print "<TD WIDTH=100><TD><INPUT NAME=terminate_task TYPE=submit VALUE=\"  Cancel the task  \">\n";
	print "<TD WIDTH=100><TD><INPUT TYPE=button VALUE=\"  Close window  \" LANGUAGE=\"javascript\" onClick=\"window.close();\">\n";
	print "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
	print "<INPUT NAME=\"process_id\" TYPE=hidden VALUE=\"$pid\">\n";
	print "<INPUT NAME=logFileID TYPE=hidden VALUE=\"$logFileID\">\n";
	print "<INPUT NAME=file_new_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
	print "<INPUT NAME=description_new_matrix TYPE=hidden VALUE=\"$description_matrix\">\n";
	print "<INPUT NAME=action TYPE=hidden VALUE=\"generate_matrix\">\n";
	print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
	print "</TABLE></FORM><p><HR NOSHADE>\n";
	print "</BODY>\n";
	print "</HTML>\n";
	exit(0);
}
# Child process
my %hashPlatform;
foreach my $ref (@$series_info){
	my ($seriesID,$platformID,$sampleID,$title) = @{$ref->[0]};
	$hashPlatform{$platformID}=1;
}

#Generate annotation files:
my %noAnnotation;
foreach my $platform (keys %hashPlatform){
	my $error_message = "";
	my $warning_message = "";
	my $annotFilename = "public-$platform"."_annot.txt";
	if(!file_exist("$PATH_DATA/$annotFilename")){
		my $return = download_platform_GEO($platform,$organismID,\$error_message,\$warning_message);
		if($error_message || $warning_message){
			my $message = "$error_message. $warning_message.";
			$message =~ s/^. //;
			if($logFileID){ file_append("$message","$PATH_OUTPUT/$logFileID.txt"); }
		}
		if($return==0){
			$noAnnotation{$platform}=1;
			next;
		}
	}
}
my $log2 = log(2);
my $log10 = log(10);
my %Matrix;
my $logratio_in_matrix=0;
my %hashProbeID;
my @samples;
my @platforms = keys %hashPlatform;
%hashPlatform=();
for(my $i=0; $i<@$series_info; ++$i){
	my $ref = $series_info->[$i];
	my ($seriesID,$platformID,$sampleID,$title) = @{$ref->[0]};
	if($noAnnotation{$platformID}){
		splice(@$series_info,$i,1);
		$i--;
	}
	my $error_message = "";
	my $warning_message = "";
	my $return = download_matrix_file($seriesID,$platformID,\$error_message,\$warning_message,$organismID);
	if($error_message || $warning_message){
		my $message = "$error_message. $warning_message.";
		$message =~ s/^. //;
		if($logFileID){ file_append("$message","$PATH_OUTPUT/$logFileID.txt"); }
	}
	my $matrix_file = $seriesID."_series_matrix.txt";
	if(!file_exist("$PATH_DATA/$matrix_file")){
		$matrix_file = "$seriesID-$platformID"."_series_matrix.txt";
	}
	if($logFileID){ file_append("Reading matrix $matrix_file","$PATH_OUTPUT/$logFileID.txt"); }
	my $annot_file = $platformID;
	if(@platforms==1){ $annot_file = "native"; }
	$return = read_matrix_file($matrix_file,$annot_file,$series_info->[$i],\%Matrix,\%hashProbeID,$logFileID,$organismID);
	if(!$return){
		splice(@$series_info,$i,1);
		$i--;
		next;
	}
	if($return==2 && $i==0){
		$logratio_in_matrix=1;
	}
	foreach my $ref1 (@$ref){
		push(@$ref1,$return);
		push(@samples,$ref1);
		#print join("\t",@$ref1)."<br>\n";
	}
	$hashPlatform{$platformID}=1;
}
if(!%hashPlatform){
	error_message("Matrix assembly failed: no samples extracted",$logFileID);
}
#Save the combined matrix file
open(OUT,">$PATH_DATA/$loginname-$file_matrix") or error_message("Cannot open $file_matrix",$logFileID);
print OUT "!Series_title\t\"$file_matrix\"\n";
if($file_description){
	print OUT "!Series_summary\t\"$file_description\"\n";
}
print OUT "!Series_type\t\"Expression profiling by array\"\n";
if(@platforms > 1){
	print OUT "!Series_platform_id\t\"symbol_$organismID\"\n";
	print OUT "!Series_nonredundant\t\"true\"\n";
}else{
	print OUT "!Series_platform_id\t\"$platforms[0]\"\n";
}
print OUT "!Series_platform_taxid\t\"$organismID\"\n";
print OUT "!Series_sample_taxid\t\"$organismID\"\n";
my $nRow = keys %Matrix;
print OUT "\n!Sample_title";
for(my $i=0; $i<@samples; ++$i){
	my ($seriesID,$platformID,$sampleID,$name) = @{$samples[$i]};
	$name =~ s/^\s*|\s*$//g;
	$name =~ s/,*[\s_]+(biological[\s_]+|)(rep|replicate|replication)[\s_]*\d+$//i;
	#$name .= " $platformID $seriesID";
	print OUT "\t\"$name\"";
}
print OUT "\n!Sample_geo_accession";
for(my $i=0; $i<@samples; ++$i){
	print OUT "\t\"$samples[$i]->[2]\"";
}
print OUT "\n!Sample_series_geo_accession";
for(my $i=0; $i<@samples; ++$i){
	my $seriesIDshort = $samples[$i]->[0];
	$seriesIDshort =~ s/-.+$//;
	print OUT "\t\"$seriesIDshort\"";
}
print OUT "\n!Sample_taxid_ch1";
for(my $i=0; $i<@samples; ++$i){
	print OUT "\t\"$organismID\"";
}
print OUT "\n!Sample_platform_id";
for(my $i=0; $i<@samples; ++$i){
	print OUT "\t\"$samples[$i]->[1]\"";
}
print OUT "\n!Sample_logratio_used";
for(my $i=0; $i<@samples; ++$i){
	my $x = $samples[$i]->[5] - 1;
	print OUT "\t\"$x\"";
}
print OUT "\n!Sample_data_row_count";
for(my $i=0; $i<@samples; ++$i){
	print OUT "\t\"$nRow\"";
}
print OUT "\n!series_matrix_table_begin\n";
print OUT "\"ID_REF\"";
for(my $i=0; $i<@samples; ++$i){
	print OUT "\t\"$samples[$i]->[2]\"";
}
print OUT "\n";
foreach my $symbol (sort keys %Matrix){
	if(!$symbol){ next; }
	my $ref = $Matrix{$symbol};
	my $nMissing=0;
	foreach my $x (@$ref){ if($x==$MISSING){ $nMissing++; }}
	if($nMissing <= @$ref/2){
		print OUT "$symbol\t".join("\t",@$ref),"\n";
	}
}
print OUT "!series_matrix_table_end\n";
close OUT;
#if($loginname eq "public"){
#	`chmod 666 $PATH_DATA/$loginname-$file_matrix`;
#}
open(OUT, ">$PATH_INFO/$loginname-config1.txt");
open(INFO,"<$PATH_INFO/$loginname-config.txt");
my $nLines;
while(my $line = <INFO>){
	$nLines++;
	if($line !~ /^type_matrix=$file_matrix\s/ && length($line)>2){
		print OUT $line;
	}
}
close INFO;
print OUT "type_matrix=$file_matrix\torganismID=$organismID\tdescription=$file_description\tdate=$date_record\n";
close OUT;
my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
if($nLines1 && $nLines1 > $nLines*0.9){
	copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
}else{
	error_message("Failed to update configuration file!");
}
unlink("$PATH_INFO/$loginname-config1.txt");
if(file_exist("$PATH_DATA/$loginname-anova-$file_matrix")){
	unlink "$PATH_DATA/$loginname-anova-$file_matrix";
}
#Normalize matrix
if($logFileID){ file_append("Normalizing matrix and running ANOVA","$PATH_OUTPUT/$logFileID.txt"); }
$hashInput{"runID"}=$RUN_ANOVA_NORM;
my @param = ($file_matrix,$logFileID);
if($logratio_in_matrix){ push(@param,1); }
run_anova1(@param);
$hashInput{"runID"}=-1;
if($logFileID){ file_append("Task completed","$PATH_OUTPUT/$logFileID.txt"); }
exit(0);
}

#**************************************
sub   read_platform_annotation
#**************************************
{
my $platform = shift;
my $annotation = shift;

my $organismID = $hashInput{"organismID"};
my $count=0;
if($annotation){
	if(ref($annotation) eq 'HASH'){ %$annotation=(); }
	else{ return $count; }
}
my $annot_file = "$PATH_DATA/$loginname-$platform"."_annot.txt";
if(!open(INFO,'<',$annot_file)){
	$annot_file = "$PATH_DATA/$platform"."_annot.txt";
	if(!open(INFO,'<',$annot_file)){
		$annot_file = "$PATH_DATA/public-$platform"."_annot.txt";
		if(!open(INFO,'<',$annot_file)){
			if($platform =~ /^GPL\d+/){
				my ($warning_message,$error_message);
				if(!download_platform_GEO($platform,$organismID,\$error_message,\$warning_message)){
					return $count;
				}
				if(!open(INFO,'<',$annot_file)){ return $count; }
			}else{
				return $count;
			}
		}
	}
}
my $line;
while($line = <INFO>){
	if($line !~ /^!/){ last; }
}
my @items = split(/\t/, $line);
my $altColumn=0;
if($items[3] =~ /^ProbeName$/i){ $altColumn=1; }
while(my $line = <INFO>){
	chop $line;
	my ($probe_id,$symbol,$geneName,$probeName) = split(/\t/, $line);
	if($symbol){ $count++; }
	if($probe_id){
		if($symbol){ $annotation->{$probe_id} = $symbol; }
		else{ $annotation->{$probe_id} = "none"; }
	}
	if($altColumn && $probeName && !$annotation->{$probeName}){
		if($symbol){ $annotation->{$probeName} = $symbol; }
		else{ $annotation->{$probeName} = "none"; }
	}
}
close INFO;
return $count;
}

#**************************************
sub  find_specific_series
#**************************************
{
my $seriesID_ref = shift;
my $series_info = shift;
my $organismID = $hashInput{"organismID"};
my $comments="";

if(!$seriesID_ref || ref($seriesID_ref) ne 'ARRAY'){
	error_message("Arguments in find_specific_series<br>");
}
my $host = "ftp.ncbi.nlm.nih.gov";
my $user = "anonymous";
my $password = "user\@comcast.net";
use Net::FTP;
#my $f = Net::FTP->new($host, Port => 21, Passive => 0) or return("Can't open $host<br>\n");
my $f = Net::FTP->new($host) or return("Can't open $host<br>\n");
if(!$f){ return("FTP to GEO: Can't open $host"); }
if(!$f->login($user, $password)){ return("FTP to GEO: Can't log user in<br>"); }
$f->binary();
foreach my $seriesID (@$seriesID_ref){
	my $len = length($seriesID);
	if($len<6){ $len=6; }
	my $part = substr($seriesID,0,$len-3)."nnn";
	if(!$f->cwd("/geo/series/$part/$seriesID/matrix")){ $comments .= "FTP: cannot cwd to $seriesID/matrix<br>"; next; }
	my @matrix = $f->ls;
	foreach my $file_matrix (sort @matrix){
		if(file_exist("$PATH_DATA/$file_matrix")){
			unlink "$PATH_DATA/$file_matrix";
		}
		if(!$f->get($file_matrix,"$PATH_DATA/$file_matrix")){ $comments .= "Can't get $file_matrix<br>\n"; next; }
		system("gzip","-d","$PATH_DATA/$file_matrix");
		$file_matrix =~ s/\.gz$//;
		my %hashMatrix=();
		parse_matrix_file("$PATH_DATA/$file_matrix",\%hashMatrix);
		my $series_type = $hashMatrix{"series_type"};
		my $series_title = $hashMatrix{"series_title"};
		my $platform = $hashMatrix{"series_platform_id"};
		$platform = uc($platform);
		my $nrows = $hashMatrix{"sample_data_row_count"};
		if($nrows && ref($nrows) eq 'ARRAY'){ $nrows = $nrows->[0]; }
		else{ $nrows=0; }
		my $sample_type = $hashMatrix{"sample_type"};
		if($sample_type && ref($sample_type) eq 'ARRAY'){ $sample_type = $sample_type->[0]; }
		my $taxid = $hashMatrix{"series_sample_taxid"};
		my $numOrganism = $hashOrganismNum{$taxid};
		my @sampleID;
		my @sample_title;
		my $reject = 0;
		my $ref = $hashMatrix{"sample_title"};
		if(!$ref || ref($ref) ne 'ARRAY'){
			$reject = 1;
			$comments .= "File $file_matrix has no sample titles<br>\n";
		}else{
			@sample_title = @$ref;
		}
		$ref = $hashMatrix{"sample_geo_accession"};
		if(!$ref || ref($ref) ne 'ARRAY'){
			$reject = 1;
			$comments .= "File $file_matrix has no sample IDs<br>\n";
		}else{
			@sampleID = @$ref;
		}
		if(!$platform){
			$reject = 1;
			$comments .= "File $file_matrix has no_platform<br>\n";
		}
		if(!$taxid){
			$reject = 1;
			$comments .= "File $file_matrix has no taxid<br>\n";
		}
		if(!$numOrganism){
			$reject = 1;
			$comments .= "File $file_matrix is for organismID $taxid which is not supported by ExAtlas<br>\n";
		}
		if($sample_type !~ /^RNA/i){
			$reject = 1;
			$comments .= "Sample type in $file_matrix is $sample_type (not RNA)<br>\n";
		}
		if($numOrganism && ($nrows < $organisms[$numOrganism-1]->[1] || $nrows > 120000)){
			my $many = "few (<$organisms[$numOrganism-1]->[1])";
			if($nrows > 120000){ $many = "many (>120000)"; }
			$comments .= "Warning: file $file_matrix has too $many rows: $nrows<br>\n";
		}
		if($series_type =~ /expression profiling by array/i){
			if($reject){ next; } 
			my $seriesIDext = $seriesID;
			if($series_info->{$seriesIDext}){
				my $num = 1;
				while($series_info->{"$seriesID-$num"}){ $num++; }
				$seriesIDext = "$seriesID-$num";
			}
			my @lines;
			for(my $i=0; $i<@sampleID; ++$i){
				push(@lines,"$sampleID[$i]\t$sample_title[$i]");
			}
			if(@lines){
				my $n_samples = @sampleID;
				$series_info->{$seriesIDext} = [$series_title,100,$n_samples,$platform,\@lines];
			}
		}else{
			$comments .= "Series type \"$series_type\" is not supported by ExAtlas<br>\n";
			$comments .= "If this is RNA-seq data, you can download it from GEO, process manually and then upload to ExAtlas<br>\n";
			$comments .= "Here is information about data series $file_matrix<br><b>Title:</b>$series_title<b>Platform:</b>$platform<p>\n";
		}
	}
}
$f->quit();
if(!%$series_info){
	$comments = "<h2>Specified data cannot be processed</h2>\n".$comments;
}
return $comments;
}

#************************************
sub  download_platform_GEO
#***********************************
{
my $platform = shift;
my $organismID = shift;
my $error_ref = shift;
my $warning_ref = shift;
my $comment = shift;

my $len = length($platform);
if($platform !~ /^GPL\d+$/ && $len > 13){
	$$error_ref = "Unknown platform $platform"; return 0;
}
use LWP;
my $browser = LWP::UserAgent->new;
my $content = $browser->get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=$platform");
if(!$content){
	$$error_ref = "Platform not found in GEO"; return 0;
}
$content =~ s/[^[:ascii:]]//g;
my @lines = split(/\n/,$content);
my ($title,$technology,$taxonomy,$rows,$taxid,$type,$manufacturer,$description,$Nsamples);
my $found = 0;
my $min_rows = 0;
my $iline = 200;
while($iline<@lines){
	my $line = $lines[$iline++];
	if($line =~ /^<tr valign="top"><td nowrap>Title/i){
		$title = $lines[$iline++];
		$title =~ s/<[^>]+>//g;
	}
	elsif($line =~ /^<tr valign="top"><td nowrap>Technology type/i){
		$technology = $lines[$iline++];
		$technology =~ s/<[^>]+>//g;
	}
	elsif($line =~ /^<tr valign="top"><td nowrap>Description/i){
		my @items = split(/<[^>]+>/,$lines[$iline++]);
		$description = $items[1];
	}
	elsif($line =~ /^<tr valign="top"><td>Samples/i){
		$description = $lines[$iline++];
		my @items = split(/<[^>]+>/,$lines[$iline++]);
		$Nsamples = $items[1];
		$Nsamples =~ s/^Samples \(//;
		$Nsamples =~ s/\)$//;
		#print "N sample $platform = $Nsamples\n";
	}
	elsif($line =~ /^<tr valign="top"><td nowrap>Organism/i){
		my $line = $lines[$iline++];
		$line =~ s/<[^>]+>//g;
		my @species = split(/; /,$line);
		for(my $i=0; $i<@species; $i++){
			my @words = split(/ /,$species[$i]);
			if(@words>2){ $species[$i] = "$words[0] $words[1]"; }
			my $iii = $hashSpeciesNum{$species[$i]};
			if($iii){
				if($taxid){ $taxid .= ";"; $taxonomy .= ";"; }
				$taxid .= $organisms[$iii-1]->[0];
				$taxonomy .= $species[$i];
				if($min_rows < $organisms[$iii-1]->[1]){ $min_rows=$organisms[$iii-1]->[1]; }
			}
		}
		#print "Organism = $taxid\n";
	}
	elsif($line =~ /^<br>Total number of rows/i){
		my @items = split(/<[^>]+>/,$line);
		$rows = $items[2];
		last;
	}
}
if(!$taxonomy){ $$error_ref = "Species not supported"; return 0; }
if($title =~ /CGH|SNP|HELP|ROMA|CNV|HiSeq/ || $title =~ /promoter|tiling|chip-on-chip|chr[123456789XY _]/i ||
 $description =~ /CGH|SNP|HELP|CNV|MIP/ || $description =~ /promoter|tiling|chip-on-chip|chr[123456789XY _]/i){
	$$error_ref = "Wrong platform type (e.g., genomic, tiling)"; return 0;
}
$type = "mRNA";
if($title =~ /miR|microR|piRNA|PIWI|piwi/){ $type = "miRNA"; }
$title = s/^\[[^\]]*\]\s*//;
my $file1 = "public-$platform"."_annot.txt";
my $taxid1 = $taxid;
$taxid1 =~ s/;/,/g;
if($type eq "miRNA"){
	system("$PATH_BIN/download_miRNA_platform.pl","$platform","$file1","$taxid1","-data","$PATH_DATA");
}else{
	system("$PATH_BIN/download_platform.pl","$platform","$file1","$taxid1","-data","$PATH_DATA");
}
`chmod 666 $PATH_DATA/$file1`;
#print "$platform\t$rows\t$type\t$taxid\t$technology\n";
add_public_annotation($platform,$title,$taxid,$type);
return 1;
}

#***********************************
sub composite_string
#***********************************
{
my $word = shift;
my $n=0;
$n++ while ($word =~ m/[,=\/\|]/g);
if($n>4){ return 1; }
return 0;
}

#************************************
sub  download_matrix_file
#***********************************
{
my $seriesID = shift;
my $platformID = shift;
my $error_ref = shift;
my $warning_ref = shift;
my $organismID = shift;

my $len = length($seriesID);
if($seriesID !~ /^GSE\d+$|^GSE\d+-GPL\d+$/ || $len > 20){
	$$error_ref = "Unknown series ID $seriesID"; return 0;
}
# Check if matrix file exists
my $file_matrix = $seriesID."_series_matrix.txt";
if(file_exist("$PATH_DATA/$file_matrix")){
	$$warning_ref = "Using old matrix file $file_matrix";
	return 1;
}
my $file_matrix1 = $file_matrix.".gz";
my $response;
if(file_exist("$PATH_DATA/$file_matrix1")){
	system("gzip","-d","$PATH_DATA/$file_matrix1");
	#unlink "$PATH_DATA/$file_matrix1";
	my $return = check_matrix_headers("$PATH_DATA/$file_matrix",$organismID,$error_ref);
	if($return==0){ return 0; }
	$$warning_ref = "Old matrix file unzipped";
	return 1;
}
my $file_matrix_long = "$seriesID-$platformID"."_series_matrix.txt";
if(file_exist("$PATH_DATA/$file_matrix_long")){
	$$warning_ref = "Using old matrix file $file_matrix_long";
	return 1;
}
my $file_matrix_long1 = $file_matrix_long.".gz";
if(file_exist("$PATH_DATA/$file_matrix_long1")){
	system("gzip","-d","$PATH_DATA/$file_matrix_long1");
	my $return = check_matrix_headers("$PATH_DATA/$file_matrix_long",$organismID,$error_ref);
	if($return==0){ return 0; }
	$$warning_ref = "Old matrix file unzipped";
	return 1;
}
#Initiate FTP session
my $host = "ftp.ncbi.nlm.nih.gov";
my $user = "anonymous";
my $password = "user\@comcast.net";
use Net::FTP;
#my $f = Net::FTP->new($host, Port => 21, Passive => 0);
my $f = Net::FTP->new($host);
if(!$f){ $$error_ref="FTP to GEO: Can't open $host"; return 0; }
if(!$f->login($user, $password)){ $$error_ref="FTP to GEO: Can't log $user in"; return 0; }
$f->binary();
my $seriesIDshort = $seriesID;
$seriesIDshort =~ s/-.+$//;
$len = length($seriesIDshort);
if($len < 6){ $len=6 };
my $part = substr($seriesIDshort,0,$len-3)."nnn";
if(!$f->cwd("/geo/series/$part/$seriesIDshort/matrix")){ $$error_ref="FTP: cannot cwd to $seriesIDshort/matrix"; return 0; };
my @matrix = $f->ls;
for(my $i=0; $i<@matrix; $i++){
	if($matrix[$i] eq $file_matrix_long1){
		$file_matrix = $file_matrix_long;
		$file_matrix1 = $file_matrix_long1;
	}
	if($matrix[$i] eq $file_matrix1){
		if(file_exist("$PATH_DATA/$file_matrix")){ unlink "$PATH_DATA/$file_matrix"; }
		if(file_exist("$PATH_DATA/$file_matrix1")){ unlink "$PATH_DATA/$file_matrix1"; }
	}
}
if(!$f->get($file_matrix1,"$PATH_DATA/$file_matrix1")){ $$error_ref="FTP: cannot get $file_matrix1\n"; return 0; }
#if(!$f->get($file_matrix1)){ $$error_ref="FTP: cannot get $file_matrix1 in $PATH_DATA\n"; return 0; }
my $count=0;
while(!$f && $count < 3){
	$f = Net::FTP->new($host, Port => 21, Passive => 0);
	$f->login($user, $password);
	$f->binary();
	$f->cwd("/geo/series/$part/$seriesIDshort/matrix");
	$f->get($file_matrix1,"$PATH_DATA/$file_matrix1");
	$count++;
}
system("gzip","-d","$PATH_DATA/$file_matrix1");
$f->quit();
my $return = check_matrix_headers("$PATH_DATA/$file_matrix",$organismID,$error_ref);
if($return==0){ return 0; }
#unlink "$PATH_DATA/$file_matrix1";
return 1;
}

#**********************
sub read_matrix_file
#**********************
{
my $matrix_file = shift;
my $annot_file = shift;
my $series_info = shift;
my $Matrix = shift;
my $hashProbeID = shift;
my $logFileID= shift;
my $organismID= shift;

my $iorg = 0;
while($iorg<@organisms && $organismID != $organisms[$iorg]->[0]){ $iorg++; }
if($iorg==@organisms){
	if($logFileID){ file_append("Organism ID not found","$PATH_OUTPUT/$logFileID.txt"); }
	return 0;
}
my %annotation;
my $count;
if($annot_file ne "native"){
	$count = read_platform_annotation($annot_file,\%annotation);
	if(!$count){
		if($logFileID){ file_append("No file $annot_file; skipping $series_info->[0]->[0]","$PATH_OUTPUT/$logFileID.txt"); }
		return 0;
	}
}
my %exprData;
if(!open (INFO,"<$PATH_DATA/$matrix_file")){
	if($logFileID){ file_append("Cannot open $matrix_file; skipping $series_info->[0]->[0]","$PATH_OUTPUT/$logFileID.txt"); }
	return 0;
}
my @sample_index;
while(my $line = <INFO>){
	$line =~ s/[^[:ascii:]]//g;
	if($line =~ /^!Series_sample_taxid/){
		my $found=0;
		while($line =~ /^!Series_sample_taxid/){
			chop $line;
			$line =~ s/[\n\r]$//;
			$line =~ s/\"//g;
			my @items = split(/\t/,$line);
			foreach my $organismID1 (split(/,/,$items[1])){
				if($organismID1==$organismID){ $found=1; last; }
			}
			$line = <INFO>;
		}
		if(!$found){
			if($logFileID){ file_append("Organism ID do not match; skipping $series_info->[0]->[0]","$PATH_OUTPUT/$logFileID.txt"); }
			return 0;
		}
	}
	if($line =~ /^!Sample_geo_accession/){
		chop $line;
		$line =~ s/[\n\r]$//;
		$line =~ s/\"//g;
		my ($junk,@GSM) = split(/\t/,$line);
		for(my $j=0; $j<@$series_info; $j++){
			for(my $i=0; $i<@GSM; $i++){
				if($series_info->[$j]->[2] eq $GSM[$i]){
					$sample_index[$j] = $i;
					last;
				}
			}
		}
		last;
	}
}
while(my $line = <INFO>){
	if($line =~ /^\!series_matrix_table_begin/i){ last; }
}
my $line = <INFO>; # skip_headers
my @allData;
while(my $line = <INFO>){
	chop $line;
	$line =~ s/[^[:ascii:]]//g;
	$line =~ s/\s+$//;
	if($line =~ /^\!series_matrix_table_end/i){ last; }
	my ($probe_id,@data) = split(/\t/, $line);
	foreach my $x (@data){
		if($x =~ /nd|none|na|n\/a|null|nan|missing/i){ next; }
		push(@allData,$x);
	}
}
close INFO;
# Find the range of data values to determine transformations:
@allData = sort {$a<=>$b} @allData;
my $nn = @allData;
my $dataMedian = $allData[int($nn*0.5)];
my $dataMin = $allData[int($nn*0.10)];
if($dataMin > 0){ $dataMin /= 10; }
my $xx = $allData[int($nn*0.05)];
if($xx > 0){
	$xx /= 3;
	if($dataMin < $xx){ $dataMin=$xx; }
}
$xx = $allData[int($nn*0.01)];
if($xx > 0 && $dataMin < $xx){ $dataMin=$xx; }
my $dataMax = $allData[int($nn*0.99)];
my $logBase = 0;
my $adjustment = 0;
my $logratio_option=0;
my $ratio = 0;
if($dataMax > 0){
	$ratio = $dataMin/$dataMax;
}
if($dataMedian >= 0.3){
	if($dataMax > 8 && $dataMax< 20){ $logBase = 2; }
	elsif($dataMax >= 2.0 && $dataMax<= 8){ $logBase = 10; }
	elsif($dataMax < 2.0){
		if($logFileID){ file_append("Data max = $dataMax; skipping $matrix_file","$PATH_OUTPUT/$logFileID.txt"); }
		return 0;
	}
	if($logFileID){ file_append("Min = $dataMin\tMax = $dataMax\tmedian = $dataMedian\tLog=$logBase","$PATH_OUTPUT/$logFileID.txt"); }
	if($logBase){
		$dataMin = exp($dataMin*log($logBase));
		$dataMax = exp($dataMax*log($logBase));
	}
	if($dataMax > 1000){
		if($dataMin > 35 && $dataMin < 75){ $adjustment = 50; }
		elsif($dataMin >= 75 && $dataMin < 140){ $adjustment = 100; }
		if($adjustment){
			$dataMin -= $adjustment;
			$dataMax -= $adjustment;
			if($dataMin <= $dataMax*0.000001){ $dataMin = $dataMax*0.000001; }
		}
	}
}elsif(abs($dataMedian) < 0.6 && $dataMax<15 && $dataMin>-15){
	$logBase = 2;
	$adjustment = 0;
	$logratio_option=1;
	#if(!file_exist("$PATH_DATA/$organismID"."_expression.txt")){
	#	if($logFileID){ file_append("Logratio & no expression file; skipping $matrix_file","$PATH_OUTPUT/$logFileID.txt"); }
	#	return 0;
	#}
}else{
	if($logFileID){ file_append("Transformation unknown (median=$dataMedian, max=$dataMax, min=$dataMin); skipping $matrix_file","$PATH_OUTPUT/$logFileID.txt"); }
	return 0;
}

#file_append("$series_info->[0]->[0]\t$dataMin\t$dataMax\t$logBase\t$adjustment\t$logratio_option","$PATH_OUTPUT/$logFileID.txt");
open (INFO,"<$PATH_DATA/$matrix_file");
while(my $line = <INFO>){
	if($line =~ /^\!series_matrix_table_begin/i){ last; }
}
$line = <INFO>; # skip_headers
my $log2 = log(2);
my $log10 = log(10);
my $refPlaftorm = $hashProbeID->{$series_info->[0]->[1]};
my $iline=0;
while(my $line = <INFO>){
	$iline++;
	#if($logFileID){ file_append("Line=$iline","$PATH_OUTPUT/$logFileID.txt"); }
	chop $line;
	$line =~ s/[^[:ascii:]]//g;
	$line =~ s/\s+$//;
	$line =~ s/\"//g;
	if($line =~ /^\!series_matrix_table_end/i){ last; }
	my ($probe_id,@data) = split(/\t/, $line);
	$probe_id =~ s/\"//g;
	if($probe_id =~ /^AFFX-.+_at$|^\(\+\)E1A_r60_|^3xSLv1|^DarkCorner|^DCP_\d+_\d|^ERCC-\d+_\d|^ETG\d+_\d|^GE_Bright|^Bright/){
		next;
	}
	my $symbol;
	if($annot_file eq "native"){
		$symbol = $probe_id;
	}else{
		$symbol = $annotation{$probe_id};
		if($symbol =~ /^none$/i){ $symbol =""; }
		if(!$symbol){ next; }
	}
	my ($sd,$sum,$nn) = (0,0);
	for(my $i=0; $i<@data; $i++){
		my $x = $data[$i];
		if($x =~ /nd|none|na|n\/a|null|nan|missing|-9999/i){ $x=$MISSING; }
		else{
			if($logBase){
				my $x1 = $x;
				if($logBase==2){ $x1 *= $log2; }
				elsif($logBase==10){$x1 *= $log10; }
				if($x1 < -30){ $x1=-30; }
				elsif($x1 > 30){ $x1=30; }
				$x = exp($x1);
				if($logratio_option){ $x *= 1000; }
			}
			if($adjustment > 0){
				$x -= $adjustment;
			}
			if($logratio_option==0){
				if($x < $dataMin/2){ $x = $dataMin/2*(0.2+0.8*rand()); }
				$x = floor($x*10000+0.5)/10000;
			}elsif($x < 1.0e-7){
				if($x > -$dataMax*0.000001){ $x = $dataMax*0.000001*(0.2+0.8*rand()); }
				else{
					$data[$i] = $MISSING;
					next;
				}
			}
			if($annot_file ne "native" && $x>1.0E-20){
				my $y = log($x);
				$sum += $y;
				$sd += $y*$y;
				$nn++;
			}
		}
		$data[$i] = $x;
	}
	if($nn > 2 && $sd*$nn > $sum*$sum){
		$sd = sqrt(($sd-$sum*$sum/$nn)/($nn-1));
	}else{
		$sd = 0;
	}
	my @data1;
	for(my $i=0; $i<@$series_info; $i++){
		my $ii = $sample_index[$i];
		if(!defined($ii)){ error_message("Sample $series_info->[$i]->[2] not found",$logFileID); }
		push(@data1,$data[$ii]);
	}
	my $ref = $exprData{$symbol};
	if(!$ref){
		$exprData{$symbol} = [$sd,$probe_id,@data1];
	}elsif($refPlaftorm && $probe_id eq $refPlaftorm->{$symbol}){
		@$ref = ($sd,$probe_id,@data1);
	}elsif($ref->[0] < $sd){
		@$ref = ($sd,$probe_id,@data1);
	}
}
close INFO;
if(!$refPlaftorm){
	foreach my $symbol (keys %exprData){
		$refPlaftorm->{$symbol} = $exprData{$symbol}->[1];
	}
}
my $nCol = 0;
foreach my $symbol (keys %$Matrix){
	$nCol = @{$Matrix->{$symbol}};
	last;
}
foreach my $symbol (keys %exprData){
	my $ref = $Matrix->{$symbol};
	my ($sd,$probe_id,@data) = @{$exprData{$symbol}};
	if(!$ref){
		$Matrix->{$symbol} = [];
		$ref = $Matrix->{$symbol};
		for(my $i=0; $i<$nCol; ++$i){
			push(@$ref,$MISSING);
		}
	}
	my $ref1 = $exprData{$symbol};
	for(my $i=0; $i<@$series_info; ++$i){
		my $x = $ref1->[$i+2];
		push(@$ref,$x);
	}
}
foreach my $symbol (keys %$Matrix){
	my $ref = $Matrix->{$symbol};
	if(@$ref == $nCol){
		for(my $i=0; $i<@$series_info; ++$i){
			push(@$ref,$MISSING);
		}
	}
	elsif(@$ref < $nCol){
		if($logFileID){ file_append("Data is missing for $symbol in $series_info->[0]->[0]","$PATH_OUTPUT/$logFileID.txt"); }
	}
}
if($logFileID){ file_append("Finished $matrix_file Logratio option=$logratio_option","$PATH_OUTPUT/$logFileID.txt"); }
if($logratio_option){ return 2; }
return 1;
}

#**********************
sub  split_items
#**********************
{
my $list = shift;
my $single = shift;
my $comma = shift;
my @list;
if($list =~ /\/\/| \/ /){ @list = split(/\s*\/\/+\s*| \/ /,$list); }
elsif($comma && $list =~ /[,;]\s*/){ @list = split(/[,;]\s*/,$list); }
else{ return $list; }
if($single==1){ return $list[0]; }
elsif($single>1 && @list>$single){ splice(@list,$single); }
return join(',',@list);
}

#***********************************
sub remove_redundant
#***********************************
{
my $ref = shift;
if(!$ref || $ref && @$ref==0){ return; }
my %hash;
my $changed=0;
for(my $i=0; $i<@$ref; $i++){
	if($hash{$ref->[$i]}){
		splice(@$ref,$i,1); $i--;
		$changed=1;
	}else{
		$hash{$ref->[$i]}=1;
	}
}
return($changed);
}

#**********************
sub add_public_annotation
#**********************
{
my $platform=shift;
my $title=shift;
my $taxid=shift;
my $type=shift;

my $organismID = $hashInput{"organismID"};
open(OUT, ">$PATH_INFO/public-config1.txt");
open(INFO,"<$PATH_INFO/public-config.txt");
my $found=0;
my $nLines;
while(my $line = <INFO>){
	$nLines++;
	if($line =~ /^type_annotation=$platform\s/){
		$found=1;
	}else{
		print OUT $line;
	}
}
close INFO;
if($found){
	close OUT;
	unlink "$PATH_INFO/public-config1.txt";
	return;
}
print OUT "type_annotation=$platform";
if($title){ print OUT "\tdescription=$title"; }
print OUT "\torganismID=$taxid\tarray_type=$type\tdate=$date_record\n";
close OUT;
my $nLines1 = get_line_counts("$PATH_INFO/public-config1.txt");
if($nLines1 > $nLines*0.9){
	copy "$PATH_INFO/public-config1.txt", "$PATH_INFO/public-config.txt";
}else{
	error_message("Failed to update configuration file!");
}
unlink "$PATH_INFO/public-config1.txt";
return;
}

#*****************************
sub get_official_symbols
#*****************************
{
my $organismID= shift;
my $symbolRef= shift;
my $aliasRef= shift;

my %hashAlias;
open (INFO_SYMB,"<$PATH_DATA/gene_info_$organismID.txt");
while(my $line = <INFO_SYMB>){
	chop $line;
	my ($org_id,$entrez_id,$symbol,$junk1,$alias_list,$altName,$junk3,$junk4,$geneName) = split(/\t/,$line);
	if(!$symbol){ next; }
	if(!$geneName){ $geneName="none"; }
	$symbolRef->{$symbol} = $geneName;
	foreach my $alias (split(/\|/,$alias_list)){
		push(@{$hashAlias{$alias}},$symbol);
		if($alias=~/-/){
			$alias=~s/-//g;
			push(@{$hashAlias{$alias}},$symbol);
		}
	}
	if(uc($symbol) ne $symbol){ $aliasRef->{uc($symbol)}=$symbol; }
}
close INFO_SYMB;
foreach my $alias (keys %hashAlias){
	my @symbols = @{$hashAlias{$alias}};
	if(@symbols==1){ $aliasRef->{$alias} = $symbols[0]; next; }
	if(@symbols>10 && length($alias)<3){ next; }
	@symbols = sort @symbols;
	my @similarity;
	my @charAlias = split(//,uc($alias));
	foreach my $symb (@symbols){
		my $s = 0;
		my @charSymb = split(//,uc($symb));
		for(my $i=0; $i<@charAlias && $i<@charSymb; $i++){
			if($charSymb[$i] eq $charAlias[$i]){ $s++; }
			else{ last; }
		}
		push(@similarity,$s);
	}
	my @sorted = sort {$similarity[$b]<=>$similarity[$a]} 0..(@symbols-1);
	if($similarity[$sorted[0]]>1){
		#print "$alias $similarity[$sorted[0]] - $symbols[$sorted[0]]\n";
		$aliasRef->{$alias} = $symbols[$sorted[0]];
		next;
	}
	while(@symbols>1 && uc(substr($symbols[0],0,2)) ne uc(substr($symbols[1],0,2))){
		shift(@symbols);
	}
	if(@symbols==1){ next; }
	my $count=2;
	for(my $i=2; $i<@symbols; $i++){ 
		if(uc(substr($symbols[$i],0,2)) eq uc(substr($symbols[0],0,2))){ $count++; }
		else{ last; }
	}
	#print "DD $alias @symbols $count\n";
	if($count >= @symbols/2){ 
		$aliasRef->{$alias} = $symbols[0];
	}
}
open(INFO_SYMB,"<$PATH_DATA/symbol_$organismID"."_annot.txt");
my $line=<INFO_SYMB>;
while(my $line=<INFO_SYMB>){
	chop $line;
	my ($key,$symbol,$geneName) = split(/\t/,$line);
	if(!$symbolRef->{$symbol}){
		$symbolRef->{$symbol} = $geneName;
		if(uc($symbol) ne $symbol){ $aliasRef->{uc($symbol)}=$symbol; }
	}
}
close INFO_SYMB;
return;
}

#**********************
sub save_list_of_genes
#**********************
{
my $text = shift;
my $organismID = shift;
my $title = shift;
my $description = shift;

my $update_symbols=0; 
if($hashInput{"update_symbols"} eq "on"){ $update_symbols=1; }
if(!$title){ $title = "Uploaded"; }
$text =~ s/\n/,/g;
$text =~ s/\"//g;
$text =~ s/\t/,/g;
$text =~ s/;/,/g;
$text =~ s/\s//g;
$text =~ s/,+/,/g;
my @items = split(/,/,$text);
my $n = @items;
my %hashFound;
for(my $i=0; $i<@items; $i++){
	#if($hashFound{$items[$i]}){ splice(@items,$i--,1); next; }
	$hashFound{$items[$i]} = 1;
}

my %hashGenes=();
my $genes_found = get_gene_annotation($organismID,\@items,\%hashGenes,"all");
my $headers = $hashGenes{"HEADERS"};
my $genelist_fileID = get_outputID(1);
open(OUT, ">$PATH_OUTPUT/$genelist_fileID.txt");
print OUT "!Genelist_title\t$title\n";
if($description){
	print OUT "!Genelist_description\t$description\n";
}
print OUT "!Genelist_N_genes\t$n\n";
if($update_symbols){
	print OUT join("\t",@$headers)."\tOrig_name\n";
}else{
	print OUT join("\t",@$headers)."\n";
}

for(my $i=0; $i<@items; $i++){
	my $orig = $items[$i];
	my $ref = $hashGenes{$items[$i]};
	if(ref($ref) eq 'ARRAY'){
		if($update_symbols){
			print OUT join("\t",@$ref,$items[$i])."\n";
		}else{
			if($items[$i] eq $ref->[0]){
				print OUT join("\t",@$ref)."\n";
			}else{
				print OUT "$items[$i]\n";
			}
		}
	}else{
		print OUT "$items[$i]\n";
	}
}
close OUT;
return $genelist_fileID;
}

#**********************
sub  save_legend
#**********************
{
my $lines = shift;
my $filename = shift;
open(OUT, ">$PATH_DATA/$loginname-$filename");
foreach my $line (@$lines){
	print OUT "$line\n";
}
close OUT;
return();
}

#**********************
sub get_gene_annotation
#**********************
{
my $organismID = shift;
my $gene_list = shift;
my $hashGenes = shift;
my $geneID_type = shift;
my $comment = shift;

my $nfound=0;
if(!$gene_list || !$hashGenes){ return; }
if(!$organismID){ error_message("Missing organism ID",$comment); }
if(!$geneID_type){ error_message("Missing geneID type",$comment); }
my %hashTemp;
$hashGenes->{"HEADERS"} = ["Symbol","Gene title"];
for(my $i=0; $i<@$gene_list; $i++){
	my $x = $gene_list->[$i];
	$hashTemp{$x} = $x;
	if($x =~ /[a-z]/){ $hashTemp{uc($x)} = $x; }
	if($x =~ /\.\d+$/){ $x =~ s/\.\d+$//; $hashTemp{$x} = $gene_list->[$i]; }
}
if(($geneID_type eq "ensembl" || $geneID_type eq "all") &&
 open(INFO,"<$PATH_DATA/ensembl_$organismID"."_annot.txt")){
	my $line=<INFO>;
	while(my $line=<INFO>){
		chop $line;
		my ($key,$symbol,$geneName) = split(/\t/,$line);
		my $orig = $hashTemp{$key};
		if(!$orig){ $orig = $hashTemp{uc($key)}; }
		if($orig && !$hashGenes->{$orig}){
			$hashGenes->{$orig} = [$symbol,$geneName]; $nfound++;
		}
	}
	close INFO;
}
if(($geneID_type =~ /^refseq|^genbank|^all/) && open (INFO,"<$PATH_DATA/symbol2refseq_$organismID.txt")){
	while(my $line = <INFO>){
		chop $line;
		my ($symbol,$name,$acc_list)=split(/\t/,$line);
		foreach my $acc (split(/,/,$acc_list)){
			my $orig = $hashTemp{uc($acc)};
			if(!$orig){ $orig = $hashTemp{uc($acc)}; }
			if($orig && !$hashGenes->{$orig}){
				$hashGenes->{$orig} = [$symbol,$name]; $nfound++;
			}
		}
	}
	close INFO;
}
if(($geneID_type eq "symbol" || $geneID_type eq "all") &&
 open(INFO,"<$PATH_DATA/symbol_$organismID"."_annot.txt")){
	my $line=<INFO>;
	while(my $line=<INFO>){
		chop $line;
		my ($key,$symbol,$geneName) = split(/\t/,$line);
		my $orig = $hashTemp{uc($symbol)};
		if(!$orig){ $orig = $hashTemp{uc($key)}; }
		if($orig && !$hashGenes->{$orig}){
			$hashGenes->{$orig} = [$symbol,$geneName]; $nfound++;
		}
	}
	close INFO;
}
if(($geneID_type eq "entrez" || $geneID_type eq "all") &&
 open(INFO,"<$PATH_DATA/entrez_$organismID"."_annot.txt")){
	while(my $line = <INFO>){
		chop $line;
		my ($key,$symbol,$geneName) = split(/\t/,$line);
		my $orig = $hashTemp{uc($symbol)};
		if(!$orig){ $orig = $hashTemp{uc($key)}; }
		if($orig && !$hashGenes->{$orig}){
			$hashGenes->{$orig} = [$symbol,$geneName]; $nfound++;
		}
	}
	close INFO;
}
if($geneID_type !~ /^symbol|^refseq|^genbank|^all|^entrez|^ensembl/ || !$nfound){
	print "WARNING: Unknown geneID type $geneID_type or annotation file not found!<br>\n";
}
return $nfound;
}

#**********************
sub get_geneset_list
#**********************
{
my @geneset_list;
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!","register");
while(my $line = <INFO>){
	chop $line;
	my %hash=();
	read_config_line($line,\%hash);
	my $file_geneset = $hash{"type_geneset"};
	if($file_geneset){
		push(@geneset_list,[$file_geneset,$hash{"description"},$hash{"organismID"}]);
	}
}
close INFO;
@geneset_list = sort {lc($a->[0]) cmp lc($b->[0])} @geneset_list;
if($loginname ne "public" && open(INFO,"<$PATH_INFO/public-config.txt")){
	my @geneset_list1=();
	while(my $line = <INFO>){
		chop $line;
		my %hash=();
		read_config_line($line,\%hash);
		my $file_geneset = $hash{"type_geneset"};
		if($file_geneset){
			push(@geneset_list1,["public-".$file_geneset,$hash{"description"},$hash{"organismID"}]);
		}
	}
	close INFO;
	push(@geneset_list, sort {lc($a->[0]) cmp lc($b->[0])} @geneset_list1);
}
return @geneset_list;
}

#**************************************
sub  gene_list_explore
#**************************************
{
my $genelist_fileID = $hashInput{"genelist_fileID"};
my $organismID = $hashInput{"organismID"};

my @file_list;
my @copy_file_list;
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!","continue");
while(my $line = <INFO>){
	chop $line;
	my %hash;
	read_config_line($line,\%hash);
	my $descr = $hash{"description"};
	my $org = $hash{"organismID"};
	my $date1 = $hash{"date"};
	if($line=~/^type_geneset/){
		push(@copy_file_list,[$hash{"type_geneset"},$descr,$org,$date1]);
	}
	my @items = split(/[=\t]/,$line);
	if($items[0] =~ /^type_/){
		push(@file_list,$items[1]);
	}
}
close INFO;
@copy_file_list = sort {lc($a->[0]) cmp lc($b->[0])} @copy_file_list;
filter_list_by_organism(\@copy_file_list, $organismID);
splice(@copy_file_list,0,0,["--- New file ---",""]);
my $file_list="";
foreach my $name (@file_list){
	if(!$file_list){ $file_list = "\"".$name."\""; }
	else{ $file_list .= ",\"".$name."\""; }
}
open(INFO,"<$PATH_OUTPUT/$genelist_fileID.txt") or error_message("Cannot open the gene list","continue");
my @gene_list;
my @header_lines;
my $description;
my $genelist_title;
my $line;
while($line=<INFO>){
	$line =~ s/\s+$//;
	if(!$line){ next; }
	if($line =~ /^!/){ 
		push(@header_lines,$line);
		my($key,$value) = split(/\t/,$line);
		if($key =~ /^!Genelist_description/i){ $description = $value; }
		elsif($key =~ /^!Genelist_title/i){ $genelist_title = $value; }
		next;
	}
	last;
}
my $isymbol=0;
my @headers = split(/\t/,$line);
while($isymbol<@headers && $headers[$isymbol] !~ /symbol/i){ $isymbol++; }
if($isymbol>=@headers){ error_message("Symbols missing","continue"); }
push(@gene_list,$line);
my %hash;
while(my $line=<INFO>){
	$line =~ s/\s+$//;
	my @items = split(/\t/,$line);
	my $symbol = $items[$isymbol];
	#if($hash{$symbol}){ next; }
	$hash{$symbol}=1;
	push(@gene_list,$line);
}
close INFO;
my $fileID1 = get_outputID(1);
my $ortholog_organismID = $hashInput{"ortholog_organismID"};
my %geneHomolog;
if($ortholog_organismID && $organismID != $ortholog_organismID){
	get_gene_homolog($organismID,$ortholog_organismID,\%geneHomolog);
}
open(OUT, ">$PATH_OUTPUT/$fileID1.txt") or error_message("Cannot open file","continue");
foreach my $line (@header_lines){ print OUT $line."\n"; }
for(my $i=0; $i<@gene_list; $i++){
	my @items = split(/\t/,$gene_list[$i]);
	my $line1 = $gene_list[$i];
	if(%geneHomolog){
		if($i==0){
			$line1 .= "\tOrtholog in $hashOrganism{$ortholog_organismID}";
		}else{
			my $symbol = $items[$isymbol];
			my $symbol1 = $geneHomolog{$items[$isymbol]};
			if($symbol1){
				$line1 .= "\t$symbol1";
			}
		}
		$gene_list[$i] = $line1;
	}		
	print OUT $line1."\n";
}
close OUT;
my @geneset_list = get_geneset_list();
filter_list_by_organism(\@geneset_list, $organismID);
my ($items,$descriptions) = get_array_lists(\@geneset_list);
print "<HTML><HEAD><TITLE>ExAtlas: Table of overexpressed genes</TITLE>\n";
print_header("update_description();");
print "<SCRIPT language=JavaScript>\n";
print "<!--\n";
print "geneset_list = new Array($items);\n";
print "geneset_description = new Array($descriptions);\n";
print "file_list = new Array($file_list);\n";
my ($items1,$descriptions1) = get_array_lists(\@copy_file_list);
print "copy_file_list = new Array($items1);\n";
print "copy_file_description = new Array($descriptions1);\n";
print "function goto_main_menu() {\n";
#print "	document.form_genelist.upload_geneset.value=\"\";\n";
print "	document.form_genelist.action.value=\"continue\";\n";
print "	document.form_genelist.target = \"\";\n";
print "	document.form_genelist.submit();\n";
print "}\n";
print "function geneset_overlap(file) {\n";
print "	document.form_genelist.upload_geneset.value=file;\n";
print "	document.form_genelist.action.value=\"geneset_overlap\";\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_genelist.target = \"_BLANK\"+x;\n";
print "	document.form_genelist.submit();\n";
print "}\n";
print "function update_description() {\n";
print "	var index;\n";
print "	index = document.form_genelist.file_geneset1.selectedIndex;\n";
print "	document.form_genelist.description_geneset1.value = geneset_description[index];\n";
print "	index = document.form_genelist.copy_file.selectedIndex;\n";
print "	document.form_genelist.description_copy_file.value = copy_file_description[index];\n";
print "}\n";
print "function get_orthologs(){\n";
print "	document.form_genelist.upload_geneset.value=\"\";\n";
print "	document.form_genelist.action.value=\"gene_list_explore\";\n";
print "	document.form_genelist.target = \"\";\n";
print "	document.form_genelist.submit();\n";
print "}\n";
print "function save_gene_set() {\n";
print "	if(!document.form_genelist.geneset_name_new.value){\n";
print "		alert(\"Please, enter geneset name\"); return false;\n";
print "	}\n";
print "	var file = document.form_genelist.copy_file.options[document.form_genelist.copy_file.selectedIndex].value;\n";
print "	if(document.form_genelist.copy_file.selectedIndex==0){\n";
print "		file = prompt(\"Provide name of a new file where to copy selected items\");\n";
print "		if(!file){ return(false); }\n";
print "		var file1=file;\n";
print "		if(file.search(/\\.txt\$/)>=0){\n";
print "			if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "			file1=file.substring(0,file.length-4);\n";
print "		}\n";
print "		if(file1.search(/^[-\\w]+\$/)<0){\n";
print "			alert(\"File name should have neither spaces nor special characters\");\n";
print "			return(false);\n";
print "		}\n";
print "		if(file.search(/^public/i) >= 0){\n";
print "			alert(\"File name cannot start with 'public'\");\n";
print "			return(false);\n";
print "		}\n";
print "		for(i=0; i<file_list.length; ++i){\n";
print "			if(file == file_list[i]){\n";
print "				alert(\"File with this name already exists\"); return false;\n";
print "			}\n";
print "		}\n";
print "		var descrip = document.form_genelist.description_copy_file.value;\n";
print "		if(descrip.search(/\\=|\\&/) >= 0){\n";
print "			alert(\"Description should not include character \'=\' or \'&\'\");\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	if(document.form_genelist.geneset_name_new.value.search(/\\=|\\&/) >= 0 || document.form_genelist.geneset_description_new.value.search(/\\=|\\&/) >= 0){\n";
print "		alert(\"Geneset name and description should not include character \'=\' or \'&\'\");\n";
print "		return false;\n";
print "	}\n";
print "	document.form_genelist.copy_to_geneset.value = file;\n";
print "	document.form_genelist.action.value = \"add_geneset\";\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.form_genelist.target = \"_BLANK\"+x;\n";
print "	document.form_genelist.submit();\n";
print "}\n";
print "<!-- end script --></SCRIPT></HEAD>\n";

print "<FONT SIZE=+2><b>List of genes '$genelist_title'</b></FONT> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; \n";
print "<INPUT TYPE=button VALUE=\" Return to main menu  \" LANGUAGE=\"javascript\" onClick=goto_main_menu();><p>\n";
my $header_line = shift(@gene_list);
my $n = @gene_list;
if($description){
	print "<b>Description</b> = $description<br>\n";
}
print "<b>Organism:</b> $hashOrganism{$organismID}<br>\n";
print "<b>Number of genes</b> = $n<br>\n";
my $x = int(10000*rand());
print "<b>Table as tab-delimited text:</b> <a href=\"$HOME_ADDRESS/output/$fileID1.txt\" target=_BLANK$x>table</a><p>\n";
print "<FORM NAME=form_genelist ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<TABLE BORDER=0>\n";
my $menu_text = menu_geneset_overlap(\@geneset_list,$genelist_fileID);
print $menu_text;
print "<TR><TD><b>Get orthologs in:<TD><select name=ortholog_organismID style=width:250px;>\n";
print "<option value=0> ---- Select species -----\n";
foreach my $ref (@organisms){
	if($organismID==$ref->[0]){ next; }
	print "<option value=$ref->[0]>$ref->[3] ($ref->[2])\n";
}
print "</select><TD COLSPAN=2>\n";
print "<INPUT TYPE=button VALUE=\"Show gene orthologs\" onClick=\"get_orthologs();\" style=width:250px;>\n";
print "</TABLE><p>\n";
$header_line =~ s/\t/<TD><b>/g;
print "<TABLE BORDER=0><TR><TD><b>$header_line\n";
foreach my $line (@gene_list){
	$line =~ s/\t/<TD>/g;
	print "<TR><TD>$line\n";
}
print "</TABLE>\n";
print "<p><font size=+2><b>Save genes as geneset</b></font> &nbsp; &nbsp; &nbsp; &nbsp;\n";
print "<INPUT TYPE=button VALUE=\"Save Genes\" onClick=save_gene_set(); style=width:250px;>\n";
print "<p><TABLE BORDER=0><TR><TD><b>Geneset name</b><TD><b>Geneset description</b><TD><b>Select geneset file</b><TD><b>File description</b>\n";
print "<TR><TD><INPUT NAME=geneset_name_new VALUE=\"$genelist_title\" style=width:220px;>\n";
print "<TD><INPUT NAME=geneset_description_new VALUE=\"$description\" style=width:220px;>\n";
print "<TD><select name=copy_file onChange=update_description(); style=width:220px;>\n";
print "<option value=0> $copy_file_list[0]->[0]\n";
for(my $i=1; $i<@copy_file_list; ++$i){
	print "<option value=\"$copy_file_list[$i]->[0]\"> $copy_file_list[$i]->[0]\n";
}
print "</select>\n";
print "<TD><INPUT NAME=description_copy_file style=width:220px;>\n";
print "</TABLE><br>\n";
print "If you select \"New file\" you will be prompted for file name, which shound be one-word without special characters (underscore allowed)<p>\n";

print "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
print "<INPUT NAME=action TYPE=hidden VALUE=continue>\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT NAME=upload_geneset TYPE=hidden>\n";
print "<INPUT NAME=fileID TYPE=hidden VALUE=$fileID1>\n";
print "<INPUT NAME=copy_to_geneset TYPE=hidden>\n";
print "<INPUT NAME=genelist_fileID TYPE=hidden VALUE=$fileID1>\n";
print "<INPUT TYPE=button VALUE=\" Return to main menu  \" LANGUAGE=\"javascript\" onClick=goto_main_menu();><p>\n";
print "</FORM><p>\n";
print "<HR NOSHADE></HR>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#********************************
sub save_expression_profile
#********************************
{
my $lines = shift;
my $organismID = shift;
my $filename = shift;
my $description = shift;

my @geneID_types = ("refseq","genbank","ensembl","entrez","symbol");
my $expression_fileID = $hashInput{"expression_fileID"};
my $add_fileID = $hashInput{"add_fileID"};
my $action = $hashInput{"action"};
my $platform = $hashInput{"platform"};
if($platform eq "None"){ $platform = ""; }
my @expression_headers;
my %hashExpression;
my @expression_data;		#previously uploaded expression data
my %hashExpressionData;		#row number for each probeID
my $error_message;
if($expression_fileID){		#Read previously uploaded expression data
	if(!open(INFO,"<$PATH_OUTPUT/$expression_fileID.txt")){
		$error_message .= "Cannot open expression_file";
		$expression_fileID = "";
	}
	while(my $line = <INFO>){
		chop $line;
		if(!@expression_headers){
			my ($junk,$junk1,@headers)=split(/\t/,$line);
			@expression_headers = @headers;
			next;
		}
		my @items = split(/\t/,$line);
		push(@expression_data,\@items);
		$hashExpressionData{$items[0]} = @expression_data; #save row number for each probeID
	} 
	close INFO;
	#Update expression headers if changed by user
	for(my $i=0; $i<@expression_headers; $i++){
		my $x = $hashInput{"expression_header$i"};
		if($x && $x ne $expression_headers[$i]){
			$expression_headers[$i] = $x;
		}
	}
}
my $nRows = $hashExpression{"n_rows"};
my $nRowsAll = 0;
my $nColsAll = 0;
my $expression_title = $hashInput{"expression_title"};
my $expression_description = $hashInput{"expression_description"};
my $geneID_type = $hashInput{"geneID_type"};

if(!$expression_title && $filename){ $expression_title = $filename; }
if(!$expression_description && $description){ $expression_description = $description; }
my $nlines=0;
if(ref($lines) eq 'ARRAY'){
	$nlines = @$lines;
}
my $startLine=0;
my $nCol=0;
my $nHeaderLines=-1;
my @headers;
my @examples;
if($nlines){	#Load a new data file and save it into $add_fileID
	my @items_old;
	my $alt_headers;
	while($startLine<@$lines){
		if($startLine == 50){ last; }
		my $line = $lines->[$startLine++];
		if($line =~ /^#|^!/){
			next;
		}
		$line =~ s/\"|\.//g;
		my ($name,@items) = split(/\t/,$line);
		if($platform =~ /^genbank_/ && $name=~ /\.\d+/){
			$name=~ s/\.\d+$//;
			$lines->[$startLine-1] = join("\t",$name,@items);
		}
		if(@items && !@items_old && !$alt_headers){
			$alt_headers = $startLine;
		}
		if(@items && @items_old >= @items){
			my $nNumCol=0;
			for(my $i=0; $i<@items; $i++){
				my $x = $items[$i];
				my $y = $items_old[$i];
				if($x =~ /^\d+$|^\d+e\-\d+$|^n\/*d|^none$|^n\/*a$|^null$|^nan$|^missing$/i &&
				 $y !~ /^\d+$|^\d+e\-\d+$|^n\/*d$|^none$|^n\/*a$|^null$|^nan$|^missing$/i){
					$nNumCol++;
				}
			}
			if(@items && $nNumCol >= 0.5*@items){
				$nHeaderLines = $startLine-1;
			}
		}
		@items_old = @items;
	}
	if($nHeaderLines<=0){
		if($alt_headers){ $nHeaderLines=$alt_headers; }
		else{ $nHeaderLines=1; }
	}
	if($startLine < 50){
		$error_message .= "Too few rows N=$nlines.<br>";
		$nlines = 0;
		$nCol = 0;
	}else{
		@headers = split(/\t/,$lines->[$nHeaderLines-1]);
		@examples = split(/\t/,$lines->[$nHeaderLines+int(($nlines-$nHeaderLines)/2)]);
		for(my $i=0; $i<@examples; $i++){
			$examples[$i] =~ s/[<>]//g;
			my $len = length($examples[$i]);
			if($len > 20){ $examples[$i] = substr($examples[$i],0,20).".."; }
		}
		$nCol = @headers;
		if(@examples < 0.5*$nCol || @examples >$nCol){ $error_message .= "Number of columns not matching<br>\n"; }
		$add_fileID = get_outputID(1);
		open(OUT, ">$PATH_OUTPUT/$add_fileID.txt");
		for(my $i=$nHeaderLines-1; $i<$nlines; $i++){
			print OUT "$lines->[$i]\n";
		}
		close OUT;
	}
}
if($action eq "add_expression_profile" && $add_fileID && !$error_message){
	open(INFO,"<$PATH_OUTPUT/$add_fileID.txt") or $error_message .= "Cannot read expression_file<br>\n";
	my $line = <INFO>;
	my $nColAdd = split(/\t/,$line);
	my ($columnProbeID,$columnGeneID)=(-1,-1);
	my ($headerProbeID,$headerGeneID);
	my @columnExpr;
	for(my $i=0; $i<$nColAdd; $i++){
		my $x = $hashInput{"column$i"};
		my $y = $hashInput{"additional_header$i"};
		if($x==1){ $columnProbeID=$i; $headerProbeID=$y; }
		elsif($x==2){ $columnGeneID=$i; $headerGeneID=$y; }
		elsif($x==3){
			push(@columnExpr,$i);
			push(@expression_headers,$y);
		}
	}
	my @data_probeID;
	my @data_geneID;
	my @data_expr;
	my @sum_expression;
	if($platform =~ /symbol_$organismID|genbank_$organismID/){
		if($columnProbeID<0 && $columnGeneID>=0){ $columnProbeID = $columnGeneID; }
		elsif($columnProbeID>=0 && $columnGeneID<0){ $columnGeneID = $columnProbeID; }
		elsif($columnProbeID<0 && $columnGeneID<0){
			$error_message .= "No column selected as gene ID for platform $platform; uploading cancelled."; 
		}
	}	
	elsif($platform && $columnProbeID<0){
		$error_message .= "No column selected as probe ID for platform $platform; uploading cancelled."; 
	}
	elsif($columnProbeID <0 && $columnGeneID>=0 && !$platform){
		$columnProbeID = $columnGeneID;
	}
	if(!@columnExpr){
		$error_message .= "No column selected as expression; uploading cancelled."; 
	}
	if($geneID_type && $headerGeneID && $geneID_type ne $headerGeneID){
		$error_message .= "GeneID name does not match to saved data; uploading cancelled."; 
	}
	my %hashPos;
	while(defined($line=<INFO>) && !$error_message){
		chop $line;
		$line =~ s/\"//g;
		my @items = split(/\t/,$line);
		my $probeID = $items[$columnProbeID];
		my @data;
		my $sum = 0;
		foreach my $i (@columnExpr){
			my $x = $items[$i];
			if($x =~ /^(n\/*d|none|n\/*a|null|nan|missing|-)$/i){ $x=$MISSING; }
			else{
				$x = floor(10000*$x+0.5)/10000;
				if($x>0.001){ $sum += log($x*1000); }
			}
			push(@data,$x);
		}
		my $pos = $hashPos{$probeID};
		my $geneID;
		if($columnGeneID>=0){
			$geneID = $items[$columnGeneID];
			if($headerGeneID =~ /^(genbank|refseq)$/i){ $geneID =~ s/\.\d+$//; }
		}
		if(!defined($pos)){
			push(@data_probeID,$probeID);
			if($columnGeneID>=0){ push(@data_geneID,$geneID); }
			push(@data_expr,\@data);
			push(@sum_expression,$sum);
			$hashPos{$probeID} = @data_expr;
		}elsif($sum > $sum_expression[$pos]){
			if($columnGeneID>=0 && !$data_geneID[$pos]){
				$data_geneID[$pos] = $geneID;
			}
			@{$data_expr[$pos]} = @data;
			$sum_expression[$pos] = $sum;
		}
	}
	close INFO;
	if(@data_probeID<10){
		$error_message .= "Too few data loaded (N<10); uploading cancelled."; 
	}
	if($platform && !$error_message){  #Check if probeID is correct
		my %hashPlatform;
		if(!read_platform_annotation($platform,\%hashPlatform)){
			$error_message .= "Failed to read platform annotation"; 
		}else{
			my @counts=(0,0);
			foreach my $probeID (keys %hashPos){
				$counts[0]++;
				if($hashPlatform{$probeID}){ $counts[1]++; }
			}
			if($counts[1]/$counts[0] < 0.5 && $columnGeneID>=0){
				$error_message .= "Probe ID check failed: Too few probe IDs recognized ($counts[1] out of $counts[0])<br>\n";
			}
		}
	}
	if(!$platform && $columnGeneID>=0 && !$error_message){  #Check if geneID is correct
		if($geneID_type){ $headerGeneID=$geneID_type; }
		my @counts;
		my %hashSymbol;
		if($headerGeneID =~ /symbol/i){
			if(!read_platform_annotation("symbol_$organismID",\%hashSymbol)){
				$error_message .= "Cannot read symbols<br>";
			}
		}
		for(my $j=0; $j<@data_geneID; $j++){
			my $x = $data_geneID[$j];
			if($x =~ /\/\/\//){ my @item = split(/\s*\/\/\/\s*/,$x); $x = $item[0]; $data_geneID[$j]=$x; }
			if($x =~ /,/){ my @item = split(/,\s*/,$x); $x = $item[0]; $data_geneID[$j]=$x; }
			if($x =~ /\|/){ my @item = split(/\s*\|+\s*/,$x); $x = $item[0]; $data_geneID[$j]=$x; }
		}
		for(my $j=0; $j<1000 && $j<@data_geneID; $j++){
			my $x = $data_geneID[$j];
			my $ucx = uc($x);
			if(!$x || $ucx =~ /^(NULL|NA|N\/A|NONE)$/){ next; }
			if($x=~ /^NM_|^XM_/){ $counts[1]++; }
			elsif($x=~ /^ENS\D*T\d\d\d\d\d\d\d/){ $counts[3]++; }
			elsif($x=~ /^\D\D\d\d+$|^\D\d\d+$|^NM_|^XM_/){ $counts[2]++; }
			elsif($headerGeneID=~/entrez/i && $x=~ /^\d+$/){ $counts[4]++; }
			elsif(%hashSymbol && $hashSymbol{$x} && $hashSymbol{$x}!~ /^none$/i){ $counts[5]++; }
			$counts[0]++;
		}
		if($counts[0]<100){ $error_message .= "Gene ID check failed: too few data entried (N=$counts[0] out of 1000)<br>\n"; }
		my @sorted = sort {$counts[$b]<=>$counts[$a]} (1..5);
		my $ig=1;
		if($geneID_type){
			while($ig <= @geneID_types && $geneID_types[$ig-1] ne $geneID_type){ $ig++; }
		}elsif($headerGeneID =~ /^(symbol|refseq|genbank|entrez|ensembl)/i){
			while($ig <= @geneID_types && $geneID_types[$ig-1] !~ /$headerGeneID/i){ $ig++; }
		}else{
			$ig = $sorted[0];
		}
		if($ig==1 && $counts[2]>0){ $ig=2; }
		if($ig==2){ $counts[2] += $counts[1]; }
		if($counts[$ig]/$counts[0] < 0.2 && $columnGeneID>=0){
			$error_message .= "Gene ID check failed: Too few gene IDs recognized ($counts[$sorted[0]] out of $counts[0])<br>\n";
		}
		if(!$geneID_type && !$error_message){
			$geneID_type = $geneID_types[$ig-1];
		}
	}
	$nRowsAll = @expression_data;
	if($nRowsAll){ $nColsAll = @{$expression_data[0]}-2; }
	if(!$error_message){
		for(my $j=0; $j<@data_probeID; $j++){
			my $probeID = $data_probeID[$j];
			my $geneID="";
			if(@data_geneID){ $geneID = $data_geneID[$j]; }
			my $irow = $hashExpressionData{$probeID};
			if(!$irow){
				push(@{$expression_data[$nRowsAll]},$probeID,$geneID);
				for(my $i=0; $i<$nColsAll; $i++){
					push(@{$expression_data[$nRowsAll]},$MISSING);
				}
				$nRowsAll++;
				$irow = $nRowsAll;
			}
			if($irow && !$expression_data[$irow-1]->[1] && $geneID){
				$expression_data[$irow-1]->[1] = $geneID;
			}
			push(@{$expression_data[$irow-1]},@{$data_expr[$j]});
		}
		# Fill rows with no new entries with missing values
		for(my $j=0; $j<$nRowsAll; $j++){
			if(@{$expression_data[$j]} == $nColsAll){
				for(my $i=0; $i<@columnExpr; $i++){
					push(@{$expression_data[$nRowsAll]},$MISSING);
				}
			}
		}
		#Save combined file
		my $fileID = get_outputID(1);
		$nColsAll += @columnExpr;
		open(OUT, ">$PATH_OUTPUT/$fileID.txt") or $error_message .= "Cannot write to expression_file<br>\n";
		print OUT "probeID\tgeneID\t".join("\t",@expression_headers)."\n";
		for(my $j=0; $j<$nRowsAll; $j++){
			print OUT join("\t",@{$expression_data[$j]})."\n";
		}
		close OUT;
		if(!$error_message){
			if(!$expression_fileID){ $expression_fileID = $fileID; }
			else{
				copy "$PATH_OUTPUT/$fileID.txt", "$PATH_OUTPUT/$expression_fileID.txt";
			}
		}
	}
}
if($action eq "save_expression_profile" && !$expression_fileID){
	$error_message .= "The file is empty, nothing to save. To cancel click on \'Return to main menu\'.<br>\n";
}
if($action eq "save_expression_profile" && !$platform && !$geneID_type){
	$error_message .= "No gene ID information! Add annotation file that associates probe ID with gene ID.<br>\n";
}
if($action eq "save_expression_profile" && !$error_message){
	open(OUT, ">$PATH_DATA/$loginname-$expression_title") or $error_message .= "Cannot write to the output file<br>\n";
	print OUT "!Series_title\t\"$expression_title\"\n";
	if($expression_description){
		print OUT "!Series_summary\t\"$expression_description\"\n";
	}
	my $make_platform_annotation=0;
	if(!$platform){
		if($geneID_type eq "symbol"){ $platform = "symbol_$organismID"; }
		elsif($geneID_type eq "refseq" || $geneID_type eq "genbank"){ $platform = "genbank_$organismID"; }
		else{
			$make_platform_annotation = 1;
			$platform = $expression_title;
			$platform =~ s/\.txt$//;
		}
	}
	my @allData;
	print OUT "!Series_platform_id\t\"$platform\"\n";
	print OUT "!Series_platform_taxid\t\"$organismID\"\n";
	print OUT "!Series_sample_taxid\t\"$organismID\"\n";
	#Remove genes with all zeroes
	for(my $j=0; $j<@expression_data; $j++){
		my ($probeID,$geneID,@data) = @{$expression_data[$j]};
		my $found=0;
		for(my $i=0; $i<@data; $i++){
			my $x = $data[$i];
			if($x =~ /e-\d\d/i){
				while($x =~ /e-0\d/i){ $x =~ s/e-0/e-/i; }
				if($x =~ /e-\d\d/i){ $x=0; }
				$expression_data[$j]->[$i+2]=$x;
			}
			elsif($x =~ /^(nd|none|na|n\/a|null|nan|missing)$/i){ $x=$MISSING; $expression_data[$j]->[$i+2]=$MISSING; }
			if($x != $MISSING && $x > 0.000001){
				$found=1;
				push(@allData, $x);
			}
		}
		#if(!$found){
		#	splice(@expression_data,$j--,1);
		#}
	}
	@allData = sort {$a<=>$b} @allData;
	my $nn = @allData;
	my $dataMin = $allData[int($nn*0.10)]/10;
	my $xx = $allData[int($nn*0.05)]/3;
	if($dataMin < $xx){ $dataMin=$xx; }
	$xx = $allData[int($nn*0.01)];
	if($dataMin < $xx){ $dataMin=$xx; }
	my $dataMax = $allData[int($nn*0.99)];
	$nRowsAll = @expression_data;
	# Check if data are 2-color arrays with reference as one of the channels
	my $reference=0;
	if(@expression_headers%2==0 && $expression_headers[0] =~ /^reference$/i || $expression_headers[1] =~ /^reference$/i){
		$reference = 1;
		for(my $i=0; $i<@expression_headers; $i+=2){
			if($expression_headers[$i] !~ /^reference$/i && $expression_headers[$i+1] !~ /^reference$/i){
				$reference = 0;
			}
			if($expression_headers[$i] =~ /^reference$/i && $expression_headers[$i+1] =~ /^reference$/i){
				$reference = 0;
			}
		}
	}
	if($reference){  #Normalize data by reference
		my $log10 = log(10);
		my @headers1;
		my @swap;
		for(my $i=0; $i<@expression_headers; $i+=2){
			my ($sw,$hd)=(0,$expression_headers[$i]);
			if($expression_headers[$i] =~ /^reference$/i){
				$sw =1;
				$hd = $expression_headers[$i+1];
			}
			push(@headers1,$hd);
			push(@swap,$sw);
		}
		@expression_headers = @headers1;
		for(my $j=0; $j<@expression_data; $j++){
			my ($probeID,$geneID,@data) = @{$expression_data[$j]};
			my @grn1=();
			my @red1=();
			my ($mr,$mg,$nr,$ng)=(0,0);
			for(my $i=0; $i<@data; $i+=2){
				my ($x,$y) = ($data[$i],$data[$i+1]);
				if($swap[$i/2]){ ($x,$y) = ($y,$x); } 
				if($x>$MISSING){ push(@grn1,$x); }
				else{
					if($x < $dataMin/2){ $x = $dataMin/2; }
					my $z = log($x);
					$mg += $z;
					$ng++;
				}
				if($y>$MISSING){ push(@red1,$y); }
				else{
					if($y < $dataMin/2){ $y = $dataMin/2; }
					my $z = log($y);
					$mr += $z;
					$nr++
				}
			}
			if($nr){ $mr /= $nr; }
			if($ng){ $mg /= $ng; }
			#my ($r,$nn,$slope) = pearson_correlation(\@grn1,\@red1);
			splice(@{$expression_data[$j]},2);
			for(my $i=0; $i<@grn1; $i++){
				my $x = $grn1[$i];
				my $y = $red1[$i];
				if($x > $MISSING){
					#if($y > $MISSING && !($mr-$mg < -1.1 && $slope > 0.3)){
						$x = exp($x-$y+$mr);
					#}else{
					#	$x = exp($x);
					#}
				}
				if($x > $MISSING && $x < $dataMin/2){
					$x = $dataMin/2*(0.2+0.8*rand());
				}
				if($x > 10){ $x = floor(1000*$x+0.5)/1000; }
				elsif($x > 0.1){ $x = floor(100000*$x+0.5)/100000; }
				$expression_data[$j]->[$i+2] = $x;
			}
		}
	}
	my $nCol = @expression_headers;
	print OUT "!Sample_title";
	for(my $i=0; $i<$nCol; ++$i){
		my $x = $expression_headers[$i];
		if($x =~ /(replicate|replication|rep)[\s_]*\d+$/i){
			$x =~ s/,*[\s_]+(biological[\s_]+|)(replicate|replication|rep)[\s_]*\d+$//i;
			$expression_headers[$i] = $x;
		}
		print OUT "\t\"$x\"";
	}
	print OUT "\n!Sample_data_row_count";
	for(my $i=0; $i<@expression_headers; ++$i){
		print OUT "\t\"$nRowsAll\"";
	}
	print OUT "\n!series_matrix_table_begin\n";
	print OUT "ID_REF";
	for(my $i=0; $i<@expression_headers; ++$i){
		my $x = $expression_headers[$i];
		print OUT "\t$x";
	}
	print OUT "\n";
	my @gene_list;
	for(my $j=0; $j<$nRowsAll; $j++){
		my @data = @{$expression_data[$j]};
		my $geneID = splice(@data,1,1);
		if($geneID){ push(@gene_list,$geneID); }
		print OUT join("\t",@data)."\n";
	}
	print OUT "!series_matrix_table_end\n";
	close OUT;

	#Update configuration file
	my $fileID1 = get_outputID(1);
	open(OUT, ">$PATH_OUTPUT/$fileID1.txt");
	open(INFO,"<$PATH_INFO/$loginname-config.txt");
	my $nLines;
	while(my $line = <INFO>){
		$nLines++;
		if($line =~ /^type_matrix=$expression_title\s/){
			if(file_exist("$PATH_DATA/$loginname-anova-$expression_title")){
				unlink "$PATH_DATA/$loginname-anova-$expression_title";
			}
		}elsif($line !~ /^type_annotation=$expression_title\_annot\.txt\s/ && length($line)>2){
			print OUT $line;
		}
	}
	close INFO;
	print OUT "type_matrix=$expression_title\torganismID=$organismID";
	if($expression_description){ print OUT "\tdescription=$expression_description"; }
	print OUT "\tdate=$date_record\n";
	if($make_platform_annotation){
		print OUT "type_annotation=$expression_title"."_annot.txt\torganismID=$organismID\tdescription=For $expression_title\tdate=$date_record\n";
	}
	close OUT;
	my $nLines1 = get_line_counts("$PATH_OUTPUT/$fileID1.txt");
	if($nLines1 && $nLines1 > $nLines*0.9){
		copy "$PATH_OUTPUT/$fileID1.txt", "$PATH_INFO/$loginname-config.txt";
	}else{
		error_message("Failed to update configuration file!");
	}
	if($make_platform_annotation){	#SAVE PLATFORM ANNOTATION
		my %hashGenes;
		get_gene_annotation($organismID,\@gene_list,\%hashGenes,$geneID_type,"continue");
		open(OUT, ">$PATH_DATA/$loginname-$platform"."_annot.txt") or $error_message .= "Cannot write to the annotation file<br>\n";
		print OUT "ProbeID\tSymbol\tGene title\n";
		for(my $j=0; $j<$nRowsAll; $j++){
			my ($probeID,$geneID) = @{$expression_data[$j]};
			print OUT $probeID;
			if($geneID){
				my $ref = $hashGenes{$geneID};
				if(ref($ref) eq 'ARRAY'){ print OUT "\t$ref->[0]\t$ref->[1]"; }
			}
			print OUT "\n";
		}
		close OUT;
	}
	if(!$error_message){ terminal_window("<H2>Expression profile data is saved</H2>","continue"); }
}
print "<HTML><HEAD><TITLE>ExAtlas: Compile gene expression profile</TITLE>\n";
print_header();
print "<SCRIPT language=JavaScript>\n";
print "<!--\n";
print "function goto_main_menu(file) {\n";
print "	if(!confirm(\"The expression profile is not saved and will be lost. Do you want to proceed?\")){ return false; }\n";
print "	document.cancel.submit();\n";
print "}\n";
print "function save_expression(){\n";
print "	document.form_expression.action.value = \"save_expression_profile\";\n";
print "	document.form_expression.submit();\n";
print "}\n";
if($nlines){
	print "function add_data(){\n";
	print "	var nSelected = new Array(0,0,0,0);\n";
	print "	for(i=0; i<document.form_expression.elements.length; ++i){\n";
	print "		var elementName = document.form_expression.elements[i].name;\n";
	print "		if(elementName.search(/^column\\d+/)>=0){\n";
	print "			var x = document.form_expression.elements[i].selectedIndex;\n";
	print "			nSelected[x]++;\n";
	print "			if(x==2 && document.form_expression.elements[i-1].value.search(/^(symbol|refseq|genbank|entrez|ensembl)/) <0){\n";
	print "				alert(\"Rename column header as either \'symbol\', \'refseq\', \'genbank\', \'entrez\', or \'ensembl\'\"); return false;\n";
	print "			}else if(x==3 && !document.form_expression.elements[i-1].value){\n";
	print "				alert(\"One of the columns with expression values has empty name! Fill it in.\"); return false;\n";
	print "			}\n";
	print "		}\n";
	print "	}\n";
	if(!$platform){
		print "	if(nSelected[2]==0){ alert(\"No column selected as Gene ID/name; file upload failed.\\nIf Gene ID is not available, click 'back' & specify platform annotation\\nbefore uploading files.\"); return false; }\n";
	}else{
		print "	if(nSelected[1]==0){ alert(\"No column selected as Probe/tracking ID that matches platform $platform!\"); return false; }\n";
	}
	print "	if(nSelected[1] > 1){ alert(\"Select no more than one column as Probe/tracking ID\"); return false; }\n";
	print "	if(nSelected[2] > 1){ alert(\"Select no more than one column as Gene ID/name\"); return false; }\n";
	print "	if(nSelected[3]==0){ alert(\"Select at least one column as expression\"); return false; }\n";
	print "	document.form_expression.action.value =\"add_expression_profile\";\n";
	print "	document.form_expression.submit();\n";
	print "}\n";
	print "function cancel_file(){\n";
	print "	for(i=0; i<document.form_expression.elements.length; ++i){\n";
	print "		var elementName = document.form_expression.elements[i].name;\n";
	print "		if(elementName.search(/^column\\d+/)>=0){\n";
	print "			document.form_expression.elements[i].selectedIndex=0;\n";
	print "		}\n";
	print "	}\n";
	print "	document.form_expression.action.value =\"add_expression_profile\";\n";
	print "	document.form_expression.submit();\n";
	print "}\n";
}
if($expression_fileID){
	print "function upload_file(){\n";
	print "	document.form_expression.action.value =\"upload_expression_profile\";\n";
	print "	document.form_expression.submit();\n";
	print "}\n";
	print "function finish_and_save(){\n";
	print "	document.form_expression.action.value =\"save_expression_profile\";\n";
	print "	document.form_expression.submit();\n";
	print "}\n";
}
print "<!-- end script --></SCRIPT></HEAD>\n";
print "<p><FONT SIZE=+2><b>Compile gene expression profile</b></FONT> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; \n";
print "<INPUT TYPE=button VALUE=\"  Cancel and return to main menu  \" LANGUAGE=\"javascript\" onClick=goto_main_menu();><p>\n";
print "<b>File name</b> = $expression_title<br>\n";
if($expression_description){
	print "<b>Description</b> = $expression_description<br>\n";
}
if($platform){
	print "<b>Platform</b> = $platform<br>\n";
}else{
	my $geneID_type1 = "Not specified";
	if($geneID_type){ $geneID_type1=$geneID_type; }
	print "<b>GeneID type</b> = $geneID_type1<br>\n";
}
if($expression_fileID){
	my $x = int(10000*rand());
	print "View data table as tab-delimited text: <a href=$HOME_ADDRESS/output/$expression_fileID.txt target=_blank$x>Data table</a><p>\n";
}
if($error_message){
	print "<p><FONT SIZE=+1 COLOR=red><b>Errors in data processing:</b></FONT><br>$error_message\n";
}
print "<p><FORM NAME=form_expression ENCTYPE=multipart/form-data ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
my $text_saved;
if($expression_fileID){
	$text_saved .= "<p><FONT SIZE=+1><b>List of gene expression data that have been uploaded</b></FONT><br>\n";
	$text_saved .= "<b>Note:</b> You can edit sample names: samples with identical names will be treated as replications<br>\n";
	$text_saved .= "<TABLE BORDER=0><TR><TD WIDTH=30><b>No</b><TD><center><b>Column name</b></center>\n";
	my $ref = $hashExpression{"source_file"};
	for(my $icol=0; $icol<@expression_headers; $icol++){
		my $colNo = $icol+1;
		$text_saved .= "<TR><TD>$colNo.<TD><INPUT NAME=expression_header$icol SIZE=30 VALUE=\"$expression_headers[$icol]\">\n";
	}
	$text_saved .= "</TABLE><p>\n";
}
if($nlines && $nCol){
	print "<FONT SIZE=+1><b>Select columns with gene ID, gene name, and expression profile</b></FONT><p>\n";
	print "<TABLE BORDER=0><TR><TD><b>No</b><TD><b>Column name</b><TD><b>Example</b><TD><center><b>Select</b></center>\n";
	for(my $icol=0; $icol<$nCol; $icol++){
		my $colNo = $icol+1;
		my $len = length($headers[$icol]);
		if($len>40){ $headers[$icol] = substr($headers[$icol],0,40); }
		print "<TR><TD>$colNo.<TD><INPUT NAME=additional_header$icol SIZE=20 VALUE=$headers[$icol]>\n";
		print "<TD>$examples[$icol]\n";
		print "<TD><select name=column$icol style=\"width: 200px;\">\n";
		print "<option value=0> Ignore\n";
		print "<option value=1> Probe/tracking ID\n";
		print "<option value=2> Gene ID/name\n";
		print "<option value=3> Gene expression\n";
	}
	print "</TABLE><p>\n";
	print "<b>Note1:</b> Edit sample names (column names): samples with identical names will be treated as replications.\n";
	print "If you have 2-dye arrays and one channel is used for reference RNA, edit column name as 'reference'. In this case reference will be used for normalization.<br>\n";
	print "<b>Note2:</b> Probe/tracking ID is platform-specific. Platform ID should be specified in 'Upload file' screen before uploading files<br>\n";
	print "<b>Note3:</b> If platform is not specified, then use gene ID/name, e.g. symbol, GenBank, Entrez, or Ensembl ID.\n";
	print "Please, edit column header as \'symbol\', \'refseq\', \'genbank\', \'entrez\', or \'ensembl\' respectively.<br>\n";
	print "<b>Note4:</b> Probe/tracking ID or Gene ID/name should be the same for all data files that are assembled together.<p>\n";
	print "<INPUT TYPE=button VALUE=\"Add selected data\" style=\"width: 200px;\" onClick=add_data();>\n";
	print "<INPUT TYPE=button VALUE=\"Do not use this file\" style=\"width: 200px;\" onClick=cancel_file();><p>\n";
	print $text_saved;
}elsif(!$nlines){
	print $text_saved;
	if($expression_fileID){
		print "<INPUT TYPE=button VALUE=\"Finish and save data\" style=\"width: 200px;\" onClick=finish_and_save();><p>\n";
		print "<b>Save the compiled gene expression data</b><p>\n";
	}
	print "<b>Add another set of gene expression data:</b><br>\n";
	print "Paste text or browse files.<p>";
	print "<TEXTAREA NAME=pasted_text ROWS=5 COLS=75></TEXTAREA><br>\n";
	print "Or upload a file: <INPUT NAME=upload_filename TYPE=file SIZE=30>\n";
	print "<INPUT TYPE=button VALUE=\"Upload data\" style=\"width: 200px;\" onClick=upload_file();><p>\n";
	print "<b>Note:</b> The input file can be a Cufflinks output file (genes.fpkm_tracking) or\n";
	print "a tab-delimited text with column headers that specifies one or several gene identifiers.<p>\n";
}

print "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
print "<INPUT NAME=action TYPE=hidden VALUE=continue>\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT NAME=expression_title TYPE=hidden VALUE=\"$expression_title\">\n";
if($platform){ print "<INPUT NAME=platform TYPE=hidden VALUE=$platform>\n"; }
if($expression_description){
	print "<INPUT NAME=expression_description TYPE=hidden VALUE=\"$expression_description\">\n";
}
if($geneID_type){
	print "<INPUT NAME=geneID_type TYPE=hidden VALUE=$geneID_type>\n";
}
if($expression_fileID){
	print "<INPUT NAME=expression_fileID TYPE=hidden VALUE=$expression_fileID>\n";
}
if($add_fileID){
	print "<INPUT NAME=add_fileID TYPE=hidden VALUE=$add_fileID>\n";
}
print "</FORM><p>\n";
print "<HR NOSHADE></HR>\n";
print "<FORM NAME=cancel ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";

print "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
print "<INPUT NAME=action TYPE=hidden VALUE=continue>\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT TYPE=button VALUE=\"Cancel and return to main menu\" style=\"width: 250px;\" onClick=goto_main_menu();><p>\n";
print "</FORM>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#****************************
sub  get_coregulated_genes
#****************************
{
my $symbol_ref = shift;
my $xr = shift;
my $yr = shift;
my $sign = shift;

my @data;
my @genelist=();
my @EPFPlist=();
my $n = @$xr;
for(my $i=0; $i<$n; ++$i){
	my ($x,$y) = ($xr->[$i],$yr->[$i]);
	if(!defined($y)){ last; }
	if($x==$MISSING || $y==$MISSING){ next; }
	if($sign<0){
		$x = -$x;
		$y = -$y;
	}
	if($x<=0 || $y<=0){ next; }
	push(@data,[$x,$y,$symbol_ref->[$i],0,0,0,0,0]);
}
$n = @data;
if($n<5){ return ("",""); }
# Estimate rank and write it into @data, positions 3 & 4:
my @sorted = sort{$data[$b]->[0]<=>$data[$a]->[0]} (0..($n-1));
for(my $i=0; $i<$n; ++$i){
	$data[$sorted[$i]]->[3]=$i+1;
}
@sorted = sort{$data[$b]->[1]<=>$data[$a]->[1]} (0..($n-1));
for(my $i=0; $i<$n; ++$i){
	my $j = $sorted[$i];
	$data[$j]->[4]=$i+1;
	my ($x1,$y1) = ($data[$j]->[3],$data[$j]->[4]);
	$data[$j]->[5] = sqrt($x1*$x1+$y1*$y1);  # circular boundary for rank x and y
	$data[$j]->[6]=$x1*$y1/$n/$n;  # p-value
}

# Smooth p-values using sliding-window average
@sorted = sort{$data[$a]->[5]<=>$data[$b]->[5]} (0..($n-1));
my $sum=$data[$sorted[0]]->[6];
my $nsum=1;
$data[$sorted[0]]->[7] = ($sum+$data[$sorted[1]]->[6])/2;
my $half_window = 9;
for(my $ii=1; $ii<$n; $ii++){
	if($ii<=$half_window){
		if($nsum<$n){ $sum += $data[$sorted[$nsum++]]->[6]; }
		if($nsum<$n){ $sum += $data[$sorted[$nsum++]]->[6]; }
	}else{
		if($ii+$half_window<$n){
			$sum += $data[$sorted[$ii+$half_window]]->[6];
			$nsum++;
		}
		$sum -= $data[$sorted[$ii-$half_window-1]]->[6];
		$nsum--;
	}
	my $j = $sorted[$ii];
	$data[$j]->[7] = $sum/$nsum;
	if($data[$j]->[7] < $data[$sorted[$ii-1]]->[7]){
		$data[$j]->[7] = $data[$sorted[$ii-1]]->[7];
	}
}
my $EPFP1 = 1;
my $rank = $n;
for(my $ii=$n-1; $ii>=0; $ii--){
	my $j = $sorted[$ii];
	my $EPFP = $data[$j]->[7]*$n/$rank--;
	if($EPFP > $EPFP1){ $EPFP = $EPFP1; }
	else{ $EPFP1 = $EPFP; }
	if($EPFP <= 0.25){
		push(@genelist,$data[$j]->[2]);
		$EPFP = int(10000*$EPFP+0.5)/10000;
		push(@EPFPlist,$EPFP);
	}
}
@genelist = reverse(@genelist);
@EPFPlist = reverse(@EPFPlist);
my $genelist = join(',',@genelist);
my $EPFPlist = join(',',@EPFPlist);
return ($genelist,$EPFPlist);
}

#**********************
sub check_matrix_data
#**********************
{
my $lines = shift;
my $filename = shift;

my $nRow = @$lines;
my $nCol = 0;
my $iline = 0;
my $count = 0;
while($iline<@$lines){
	if($count>100){ last; }
	if($lines->[$iline] =~ /^[\!\#]/){ $iline++; next; }
	my @items = split(/\t/,$lines->[$iline]);
	if($nCol < @items){
		$nCol = @items;
	}
	$iline++;
	$count++
}
my $logFileID = get_outputID(1);
$hashInput{"logFileID"} = $logFileID;
$hashInput{"runID"} = $RUN_MATRIX_UPLOAD;
my @data = ($lines,$filename);
interrupt_program($nCol*$nRow/5000,\@data);
return;
}

#**********************
sub finish_file_upload
#**********************
{
my $lines = shift;
my $filename = shift;
my $file_type = shift;
my $organismID = shift;
my $platform = shift;

#Update configuration file
my $description = $hashInput{"description"};
my %default;
my $nLines;
open(OUT, ">$PATH_INFO/$loginname-config1.txt");
open(INFO,"<$PATH_INFO/$loginname-config.txt");
while(my $line = <INFO>){
	chop $line;
	if(!$line){ next; }
	$nLines++;
	if($line =~ /^type_matrix=\t/){ next; }
	if($line =~ /^type_$file_type=$filename\s/ || $file_type eq "annotation" && $line =~ /^type_annotation=$platform/){
		if($file_type eq "matrix" && file_exist("$PATH_DATA/$loginname-anova-$filename")){
			unlink "$PATH_DATA/$loginname-anova-$filename";
		}
	}elsif($line =~ /^mainPage=$organismID\t/){
		read_config_line($line,\%default);
	}else{
		print OUT $line."\n";
	}
}
close INFO;
if($file_type eq "annotation"){
	print OUT "type_$file_type=$platform\torganismID=$organismID";
}else{
	print OUT "type_$file_type=$filename\torganismID=$organismID";
}
if($description){ print OUT "\tdescription=$description"; }
print OUT "\tdate=$date_record\n";
$default{"file_$file_type"}=$filename;
print OUT "mainPage=$organismID";
foreach my $key (keys %default){
	if($key && $key ne "mainPage"){
		print OUT "\t$key=$default{$key}";
	}
}
print OUT "\n";
close OUT;
my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
if($nLines1 && $nLines1 > $nLines*0.9){
	copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
}else{
	error_message("Failed to update configuration file!");
}
unlink("$PATH_INFO/$loginname-config1.txt");
#Save the uploaded file
open(OUT, ">$PATH_DATA/$loginname-$filename");
foreach my $line (@$lines){
	print OUT "$line\n";
}
close OUT;
if($file_type eq "annotation" && $loginname eq "public"){
	`chmod 666 $PATH_DATA/$loginname-$filename`;
}
return;
}

#**********************
sub check_matrix_data1
#**********************
{
my $logFileID = shift;
my $data_sent = shift;

my $organismID = $hashInput{"organismID"};
my $lines = $data_sent->[0];
my $filename = $data_sent->[1];
my $description = $hashInput{"description"};

my $update_symbols=0; 
if($hashInput{"update_symbols"} eq "on"){ $update_symbols=1; }
my $startLine=0;
my $nHeaderLines=-1;
my $nRow=0;
my @nNumbers;
my %hashMatrix;
my $transform = $hashInput{"transform"};
my $adjustment = $hashInput{"adjustment"};
while($startLine<@$lines && $lines->[$startLine] =~ /^\!|^\s*$/){
	my $line = $lines->[$startLine];
	if($line =~ /\s+$/){
		$line =~ s/\s+$//;
		$lines->[$startLine]=$line;
	}
	$line =~ s/^\!//;
	$line =~ s/\"//g;
	if($line =~ /^series/i){
		my ($key,$value) = split(/\t/,$line);
		$key=lc($key);
		$hashMatrix{$key}=$value;
	}elsif($line =~ /^sample/i){
		my ($key,@values) = split(/\t/,$line);
		$key=lc($key);
		$hashMatrix{$key}=\@values;
	}
	$startLine++;
}
my $nCommentLines=$startLine;
my $platform = $hashMatrix{"series_platform_id"};
if(!$platform){
	$platform = $hashInput{"platform"};
}
my $error_message = "";
my $warning_message = "";
my %hash_annotation;
if($logFileID){
	file_append("Processing matrix file $filename...","$PATH_OUTPUT/$logFileID.txt");
}
if($platform ne "None" && !read_platform_annotation($platform,\%hash_annotation)){
	error_message("Cannot read platform annotation $platform.",$logFileID);
}
while($startLine<@$lines && $lines->[$startLine] !~ /^\!/){
	if($nHeaderLines==-1 && $nRow<30){
		my ($name,@items) = split(/ *\t */,$lines->[$startLine]);
		foreach my $x (@items){
			if($x =~ /^[-+]*\d*\.*\d+$|^[-+]*\d*\.*\d+e[-+]*\d+$|^nd$|^none$|^na$|^n\/a$|^null$|^nan$|^missing$/i){ $nNumbers[$nRow]++; }
		}
		if($nRow == 29){
			if(!$nNumbers[$nRow]){ error_message("No numbers found in file!",$logFileID); }
			for(my $i=30-2; $i>=0; $i--){
				if($nNumbers[$i] < $nNumbers[30-1]-1){
					$nHeaderLines = $i+1;
					last;
				}
			}
			if($nHeaderLines==-1){ $nHeaderLines=0; }
			if($nHeaderLines>0 && $nNumbers[$nHeaderLines] > 0 && $nNumbers[$nHeaderLines+1] - $nNumbers[$nHeaderLines] < 2){
				$nHeaderLines--;
			}
		}
	}else{
		last;
	}
	$nRow++;
	$startLine++;
}
$nRow=0;
if($nHeaderLines>1){
	splice(@$lines,$nCommentLines,$nHeaderLines-1);
	$nHeaderLines=1;
}
$startLine = $nCommentLines;
my $line = $lines->[$startLine++];
if($line =~ /\s+$/){
	$line =~ s/\s+$//;
	$lines->[$startLine-1]=$line;
}
$line =~ s/\"//g;
my ($first,@headers) = split(/ *\t */,$line);
my $data_start = $startLine;
my @probeID;
my @allData;
my $logTrans = 0;
my $nNullProbe=0;
my %hashProbe;
if($transform){ $logTrans = log($transform); }
while($startLine<@$lines && $lines->[$startLine] !~ /^\!/){
	my ($probe_id,@data) = split(/ *\t */, $lines->[$startLine]);
	$probe_id =~ s/\"//g;
	if(!$probe_id || $probe_id =~ /^AFFX-.+_at$|^\(\+\)E1A_r60_|^3xSLv1|^DarkCorner|^DCP_\d+_\d|^ERCC-\d+_\d|^ETG\d+_\d|^GE_Bright|^Bright/){
		if(!$probe_id){ $nNullProbe++; }
		splice(@$lines,$startLine,1);
		next;
	}
	if($platform =~ /^genbank_/ && $probe_id=~ /\.\d+/){
		$probe_id=~ s/\.\d+$//;
	}
	$hashProbe{"$probe_id"} = 1;
	my $n_positive = 0;
	for(my $i=0; $i<@data; $i++){
		my $x = $data[$i];
		while($x =~ /e-0\d/i){ $x =~ s/e-0/e-/i; }
		if($x =~ /e-\d\d/i){ $x=0; }
		elsif($x eq ""){ $x=$MISSING; }
		elsif($x =~ /^(n\/*d|none|n\/*a|null|nan|missing|-)$/i){ $x=$MISSING; }
		if($transform && $x>$MISSING){
			if($x*$logTrans > 40){ $x=$MISSING; }
			else{
				my $x1 = $x*$logTrans;
				if($x1 > 30){ $x1 = 30; }
				elsif($x1 < -30){ $x1 = -30; }
				$x1 = exp($x1);
				if($x1 > 10){ $x = int(1000*$x1+0.5)/1000; }
				elsif($x1 > 0.1){ $x = int(100000*$x1+0.5)/100000; }
				else{ $x = $x1; }
			}
		}
		if($adjustment && $x>$MISSING){
			$x -= $adjustment;
			if($x<$adjustment/100){
				$x = int(1000*(0.2+0.8*rand())*$adjustment/100+0.5)/1000;
			}
		}
		if($x>0){
			push(@allData,$x);
			$n_positive++;
		}
		$data[$i] = $x;
	}
	if($n_positive > 0){
		$lines->[$startLine] = join("\t",$probe_id,@data);
		push(@probeID,$probe_id);
		$nRow++;
	}else{
		splice(@$lines,$startLine,1);
		next;
	}
	if($startLine%1000==0 && $logFileID){
		file_append("Reading line $startLine","$PATH_OUTPUT/$logFileID.txt");
	}
	$startLine++;
}
if($nRow < 5){
	error_message("Too few rows N=$nRow",$logFileID);
}
@allData = sort {$a<=>$b} @allData;
my $nn = @allData;
my $dataMedian = $allData[int($nn*0.5)];
my $dataMin = $allData[int($nn*0.10)]/10;
my $xx = $allData[int($nn*0.05)]/3;
if($dataMin < $xx){ $dataMin=$xx; }
$xx = $allData[int($nn*0.01)];
if($dataMin < $xx){ $dataMin=$xx; }
my $dataMax = $allData[int($nn*0.99)];
if($dataMedian >= 0.4 || $dataMax > 2*abs($dataMin)){
	if($dataMax<20 && !$transform){
		$warning_message .= "Possibly data is logtransformed. In this case select transformation type and reload the file<br>\n";
	}
}else{
	error_message("Data max = $dataMax; Data min = $dataMin.\nApparently matrix uses logratio values, which cannot be processed.",$logFileID);
}
if($dataMin < -100 || $dataMax < 0){
	error_message("Too low negative values ($dataMin)! Check the data or select transformation.",$logFileID);
}
if($dataMax > 1000 && !$adjustment){
	my $adjust1 = 0;
	if($dataMin > 35 && $dataMin <= 75){ $adjust1=50; }
	elsif($dataMin > 75 && $dataMin <= 140){ $adjust1=100; }
	if($adjust1){
		$warning_message .= "Possibly data requires adjustment. In this case select adjustment=$adjust1 and reload the file<br>\n";
	}
}
if($dataMin <= $dataMax*0.0000001){ $dataMin = $dataMax*0.0000001; }
if(!$platform){
	error_message("Array annotation ($platform) is not found.<br>Upload annotations first, then use it for uploading expression profiles.",$logFileID);
}
if($nNullProbe){
	$warning_message .= "Some probes (N=$nNullProbe) are empty, i.e., not labeled. Data for these probes is ignored<br>\n";
}
my $nUnique = keys %hashProbe;
if($nUnique < $nRow){
	my $redundant = $nRow-$nUnique;
	$warning_message .= "Some probes (N=$redundant) are redundant, and thus, information can be lost.<br>\n";
}
my @add_lines_series;
my $add_sample_title;
if(!$hashMatrix{"series_title"}){
	if($description){
		push(@add_lines_series,"!Series_title\t\"$description\"");
	}else{
		push(@add_lines_series,"!Series_title\t\"$filename\"");
	}
}elsif(!$description){
	$hashInput{"description"} = $hashMatrix{"series_title"};
}
my $organismID1=$hashMatrix{"series_platform_taxid"};
if(!$organismID1){ $organismID1=$hashMatrix{"series_sample_taxid"}; }
if($organismID1 && $organismID1 != $organismID){
	error_message("Organism ID do not match! First select a correct organism species in the main menu, then upload your file.",$logFileID);
}
if(!$hashMatrix{"series_platform_id"}){
	push(@add_lines_series,"!Series_platform_id\t\"$platform\"");
}
if(!$hashMatrix{"series_platform_taxid"}){
	push(@add_lines_series,"!Series_platform_taxid\t\"$organismID\"");
}
if(!$hashMatrix{"series_sample_taxid"}){
	push(@add_lines_series,"!Series_sample_taxid\t\"$organismID\"");
}
if(@add_lines_series){
	splice(@$lines,0,0,@add_lines_series);
}
my $ref = $hashMatrix{"sample_title"};
if($ref && ref($ref) eq 'ARRAY'){
	@headers = @$ref;
}
my $nCol = @headers;
# Check if data are 2-color arrays with reference as one of the channels
my $reference=0;
if(@headers%2==0 && $headers[0] =~ /^reference$/i || $headers[1] =~ /^reference$/i){
	$reference = 1;
	for(my $i=0; $i<$nCol; $i+=2){
		if($headers[$i] !~ /^reference$/i && $headers[$i+1] !~ /^reference$/i){
			$reference = 0;
			$warning_message .= "Possibly missing reference column for 2-dye arrays<br>\n";
		}
		if($headers[$i] =~ /^reference$/i && $headers[$i+1] =~ /^reference$/i){
			$reference = 0;
			$warning_message .= "Possibly duplicated reference column for 2-dye arrays<br>\n";
		}
		if(!$reference){ last; }
	}
}
my @swap;
my $changed = 0;
if($reference){
	$nCol /= 2;
	my @headers1;
	for(my $i=0; $i<@headers; $i+=2){
		my ($sw,$hd)=(0,$headers[$i]);
		if($headers[$i] =~ /^reference$/i){
			$sw =1;
			$hd = $headers[$i+1];
		}
		push(@headers1,$hd);
		push(@swap,$sw);
	}
	@headers = @headers1;
	$changed = 1;
}
for(my $i=0; $i<$nCol; ++$i){
	if($headers[$i] =~ /(replicate|replication|rep)[\s_]*\d+$/i){
		$changed = 1;
		$headers[$i] =~ s/,*[\s_]+(biological[\s_]+|)(replicate|replication|rep)[\s_]*\d+$//i;
	}
}
$add_sample_title = "!Sample_title\t\"".join("\"\t\"",@headers)."\"";
$startLine=0;
while($startLine<@$lines){
	if(!$ref){
		if($lines->[$startLine] !~ /^\!/){
			$startLine--;
			last;
		}
		if($lines->[$startLine] =~ /^\!sample/i){ last; }
	}elsif($lines->[$startLine] =~ /^\!sample_title/i){ last; }
	$startLine++;
}
if(!$ref){
	splice(@$lines,++$startLine,0,$add_sample_title);
}elsif($changed){
	$lines->[$startLine] = $add_sample_title;
}
$ref = $hashMatrix{"sample_data_row_count"};
if(!$ref){
	my $row_counts = "!Sample_data_row_count";
	for(my $i=0; $i<$nCol; ++$i){
		$row_counts .= "\t\"$nRow\"";
	}
	while($startLine<@$lines && $lines->[$startLine] =~ /^\!sample/i){
		$startLine++;
	}
	splice(@$lines,$startLine,0,$row_counts);
}
while($startLine<@$lines && $lines->[$startLine] =~ /^\!|^\s*$/i){
	$startLine++;
}
if($lines->[$startLine-1] !~ /^\!series_matrix_table_begin/i){
	splice(@$lines,$startLine,0,"!series_matrix_table_begin");
	$startLine++;
}
if($changed){
	$lines->[$startLine] = join("\t",$first,@headers);
}
$startLine++;
my ($count_probes,$count_symbols)=(0,0);
my %officialSymbol;
my %alias;
my $probeid_symbol = 0;
if($platform == $organismID){ $platform ="symbol_$organismID"; }
if($platform =~ /^symbol_$organismID/){
	get_official_symbols($organismID,\%officialSymbol,\%alias);
	$probeid_symbol = 1;
}
while($startLine<@$lines){
	if($startLine%1000==0 && $logFileID){
		file_append("Processing line $startLine","$PATH_OUTPUT/$logFileID.txt");
	}
	if($lines->[$startLine] =~ /^\!/){ last; }
	while($startLine<@$lines && !$lines->[$startLine]){
		splice(@$lines,$startLine,1);
		if($startLine==@$lines){ last; }
	}
	my($probeID,@data) = split(/\t/,$lines->[$startLine]);
	if($probeid_symbol && $update_symbols && !$officialSymbol{$probeID}){
		my $symbol1 = $alias{$probeID};
		if($symbol1){ $probeID=$symbol1; }
	}
	$probeID =~ s/\"//g;
	if($platform =~ /^(ensembl|genbank)_\d/){
		$probeID =~ s/\.\d+$//;
	}
	if($platform eq "None"){
		$count_probes++;
	}else{
		my $symbol = $hash_annotation{$probeID};
		if($symbol){
			$count_probes++;
			if($symbol !~ /^none$/i){ $count_symbols++; }
		}
	}
	if($reference){
		my @grn1=();
		my @red1=();
		my ($mr,$mg,$nr,$ng)=(0,0);
		for(my $i=0; $i<@data; $i+=2){
			my ($x,$y) = ($data[$i],$data[$i+1]);
			if($swap[$i/2]){ ($x,$y) = ($y,$x); } 
			if($x==$MISSING){ push(@grn1,$x); }
			else{
				if($x < $dataMin/2){ $x = $dataMin/2; }
				my $z = log($x);
				push(@grn1,$z);
				$mg += $z;
				$ng++;
			}
			if($y==$MISSING){ push(@red1,$y); }
			else{
				if($y < $dataMin/2){ $y = $dataMin/2; }
				my $z = log($y);
				push(@red1,$z);
				$mr += $z;
				$nr++
			}
		}
		if($nr){ $mr /= $nr; }
		if($ng){ $mg /= $ng; }
		my ($r,$nn,$slope) = pearson_correlation(\@grn1,\@red1);
		@data=();
		for(my $i=0; $i<@grn1; $i++){
			my $x = $grn1[$i];
			my $y = $red1[$i];
			if($x > $MISSING){
				if($y > $MISSING && !($mr-$mg < -1.1 && $slope > 0.3)){
					$x = exp($x-$y+$mr);
				}else{
					$x = exp($x);
				}
			}
			if($x > $MISSING && $x < $dataMin/2){
				$x = $dataMin/2*(0.2+0.8*rand());
			}
			if($x > 10){ $data[$i] = floor(1000*$x+0.5)/1000; }
			elsif($x > 0.1){ $data[$i] = floor(100000*$x+0.5)/100000; }
			else{ $data[$i] = $x; }
		}
	}else{
		for(my $i=0; $i<@data; $i++){
			my $x = $data[$i];
			if($x > $MISSING && $x < $dataMin/2){
				$x = $dataMin/2*(0.2+0.8*rand());
				if($x > 10){ $data[$i] = floor(1000*$x+0.5)/1000; }
				elsif($x > 0.1){ $data[$i] = floor(100000*$x+0.5)/100000; }
				else{ $data[$i] = $x; }
			}
		}

	}
	$lines->[$startLine++] = join("\t",$probeID,@data);
}
if($startLine==@$lines || $lines->[$startLine] !~ /^\!series_matrix_table_end/i){
	$lines->[$startLine] = "!series_matrix_table_end";
}
my $prop_missing = int(1000*(1-$count_symbols/$nRow))/10;
if($prop_missing > 50){
	$warning_message .= "<br>Many probes ($prop_missing %) are not annotated! Check platform annotation.";
}
my $file_info = "<b>Organism:</b> $hashOrganism{$organismID}<br>\n";
$file_info .= "<b>N columns:</b> $nCol<br><b>N rows:</b> $nRow<br>\n";
$file_info .= "<b>Percent probes not annotated:</b> $prop_missing %<br>\n";
$file_info .= "<b>N probes with symbols:</b> $count_symbols<br>\n";
if($transform){
	$file_info .= "<b>Log base:</b> $transform<br>\n";
}
if($adjustment){
	$file_info .= "<b>Adjustment:</b> $adjustment<br>\n";
}
if($warning_message){
	$file_info .= "<b>Warnings:</b> $warning_message<br>\n";
}
if($logFileID){
	file_append("File_info\n$file_info","$PATH_OUTPUT/$logFileID.txt");
}
return $file_info;
}

#**************************************
sub   check_geneset_data
#**************************************
{
my $lines = shift;
my $organismID = shift;
my $filename = shift;
my $description = shift;

my $update_symbols=0; 
if($hashInput{"update_symbols"} eq "on"){ $update_symbols=1; }
my $startLine=0;
my %hashGeneset;
while($startLine<@$lines && $lines->[$startLine] =~ /^\!/){
	$lines->[$startLine] =~ s/\"//g;
	my $line = $lines->[$startLine];
	$line =~ s/^\!//;
	if($line =~ /^geneset/i){
		my ($key,$value) = split(/\t/,$line);
		$key=lc($key);
		$hashGeneset{$key}=$value;
	}
	$startLine++;
}
my $organismID1 = $hashGeneset{"geneset_taxid"};
my $description1 = $hashGeneset{"geneset_description"};
my $filename1 = $hashGeneset{"geneset_filename"};
if(!$filename1){
	splice(@$lines,$startLine++,0,"!Geneset_filename\t$filename");
}
if(!$description1){
	splice(@$lines,$startLine++,0,"!Geneset_description\t$description");
}
if($organismID1 && $organismID1 != $organismID){
	error_message("Organism ID do not match! First select a correct organism species in the main menu, then upload your file.","continue");
}
if(!$organismID1){
	splice(@$lines,$startLine++,0,"!Geneset_taxid\t$organismID");
}
my %hashGenes;
my %alias;
get_official_symbols($organismID,\%hashGenes,\%alias);
my ($n_genes,$n_found,$count,$count_all)=(0,0,0,0);
my $comma_separated=0;
my @removed;
while($startLine<@$lines){
	$lines->[$startLine] =~ s/\"//g;
	my $line = $lines->[$startLine];
	if($line =~ /^\!/){
		splice(@$lines,$startLine,1);
		next;
	}
	my ($name,$descrip,@genes) = split(/\t/,$line);
	if(!$comma_separated && @genes==1 && $genes[0]=~/,/){
		$comma_separated = 1;
	}
	if($comma_separated){
		@genes = split(/,/,$genes[0]);
	}
	my %hash=();
	if($name){
		@removed=(); 
		for(my $i=@genes-1; $i>=0; $i--){
			my $symbol=$genes[$i];
			if($update_symbols && !$hashGenes{$symbol}){
				my $symb1 = $alias{$symbol};
				if($symb1){
					$symbol=$symb1;
					$genes[$i]=$symbol;
				}else{
					splice(@genes,$i,1);
					$removed[$i]=1;
					next;
				}
			}
			if($hash{$symbol}){
				splice(@genes,$i,1);
				$removed[$i]=1;
			}else{
				$n_genes++;
				if(!$update_symbols || $hashGenes{$symbol}){ $n_found++; }
				$hash{$symbol}=1;
			}
		}
		$count_all++;
		if(@removed || $comma_separated){
			$lines->[$startLine] = "$name\t$descrip\t".join("\t",@genes);
		}
		$count++;
	}
	else{
		for(my $i=@removed-1; $i>=0 && $name; $i--){
			if($removed[$i]){
				splice(@genes,$i,1);
			}
		}
		if(@removed || $comma_separated){
			$lines->[$startLine] = "\t$descrip\t".join("\t",@genes);
		}
	}
	$startLine++;
}
if(!$count){ error_message("No genesets found in the file.","continue"); }
my $prop = floor(1000*$n_found/$n_genes+0.5)/1000;
my $file_info = "<b>Organism:</b>$hashOrganism{$organismID}<br>\n";
$file_info .= "<b>N geneset submitted:</b>$count_all<br>\n<b>N geneset loaded:</b>$count<br>\n";
$file_info .= "<b>Proportion of genes with official gene symbols:</b>$prop<p>\n";
if($prop < 0.5){
	$file_info .= "<font size=+1 color=red><b>WARNING:</b></font>The proportion of genes official gene symbol is too low!<br>\n";
}
return $file_info;
}

#**************************************
sub   check_samples_data
#**************************************
{
my $lines = shift;
my $organismID = shift;
my $filename = shift;
my $description = shift;

my $startLine=0;
my %hashSamples;
while($startLine<@$lines && $lines->[$startLine] =~ /^\!/){
	my $line = $lines->[$startLine];
	$line =~ s/^\!//;
	$line =~ s/\"//g;
	if($line =~ /^samples/i){
		my ($key,$value) = split(/\t/,$line);
		$key=lc($key);
		$hashSamples{$key}=$value;
	}
	$startLine++;
}
my $organismID1 = $hashSamples{"samples_taxid"};
my $description1 = $hashSamples{"samples_description"};
my $filename1 = $hashSamples{"samples_filename"};
if(!$filename1){
	splice(@$lines,$startLine++,0,"!Samples_filename\t$filename");
}
if(!$description1){
	splice(@$lines,$startLine++,0,"!Samples_description\t$description");
}
if($organismID1 && $organismID1 != $organismID){
	error_message("Organism ID do not match! First select a correct organism species in the main menu, then upload your file.","continue");
}
if(!$organismID1){
	splice(@$lines,$startLine++,0,"!Samples_taxid\t$organismID");
}
my $count=0;
while($startLine<@$lines){
	my $line = $lines->[$startLine];
	if($line =~ /^\!/){ error_message("Line starts with '!'"); }
	my ($seriesID,$platformID,$sampleID,$title) = split(/\t/,$line);
	if(!$title){ error_message("No title for sample no. $count."); }
	if($sampleID !~ /^GSM\d+$/){ error_message("Bad sample ID ($sampleID) in sample no. $count.","continue"); }
	if($seriesID !~ /^GSE\d+$/){ error_message("Bad series ID ($seriesID) in sample no. $count.","continue"); }
	if($platformID !~ /^GPL\d+$/){ error_message("Bad platform ID ($platformID) in sample no. $count.","continue"); }
	$startLine++;
	$count++;
}
if(!$count){ error_message("No samples found in the file.","continue"); }
my $file_info = "<b>Organism:</b>$hashOrganism{$organismID}<br>\n";
$file_info .= "<b>N samples:</b>$count<br>\n";
return $file_info;
}

#**************************************
sub   check_annotation_data
#**************************************
{
my $lines = shift;
my $organismID = shift;
my $filename = shift;
my $description = shift;

if($filename !~ /^[-\w]+$|^[-\w]+_annot\.txt$/){
	error_message("Annotation file name should have no spaces, dots, or special characters"); 
} 
my $update_symbols=0; 
if($hashInput{"update_symbols"} eq "on"){ $update_symbols=1; }
my $startLine=0;
my %hashAnnotation;
while($startLine<@$lines && $lines->[$startLine] =~ /^\!/){
	my $line = $lines->[$startLine];
	$line =~ s/^\!//;
	$line =~ s/\"//g;
	if($line =~ /^annotation/i){
		my ($key,$value) = split(/\t/,$line);
		$key=lc($key);
		$hashAnnotation{$key}=$value;
	}
	$startLine++;
}
my $organismID1 = $hashAnnotation{"annotation_taxid"};
my $description1 = $hashAnnotation{"annotation_description"};
my $filename1 = $hashAnnotation{"annotation_filename"};
if(!$filename1){
	splice(@$lines,$startLine++,0,"!Annotation_filename\t$filename");
}
if(!$description1){
	splice(@$lines,$startLine++,0,"!Annotation_description\t$description");
}
if($organismID1 && $organismID1 != $organismID){
	error_message("Organism ID do not match! First select a correct organism species in the main menu, then upload your file.","continue");
}
if(!$organismID1){
	splice(@$lines,$startLine++,0,"!Annotation_taxid\t$organismID");
}
$lines->[$startLine] =~ s/\"//g;
my $line = $lines->[$startLine++];
my @items = split(/\t/,$line);
if(@items < 2){
	error_message("Annotation file should have 3 columns: probeID, gene symbol, gene name.","continue");
}
if($items[1] !~ /^symbol|^gene[ _]*symbol/i){
	error_message("The header of 2nd column should be \"gene symbols\".","continue");
}
my %hashGenes;
my %alias;
get_official_symbols($organismID,\%hashGenes,\%alias);
my %hashAll;
my %hashFound;
while($startLine<@$lines){
	$lines->[$startLine] =~ s/\"//g;
	my $line = $lines->[$startLine];
	if($line =~ /^\!/){ error_message("Line starts with '!'","continue"); }
	my ($probeID,$symbol,@items) = split(/\t/,$line);
	if($symbol){
		$hashAll{$symbol}=1;
		my $changed=0;
		if($symbol =~ /^\s|\s$/){
			$symbol =~ s/^\s+|\s+$//g; $changed=1;
		}
		if($symbol =~ /\/\//){
			$symbol =~ s/\s*\/\/.*$//;
			if($items[0] =~ /\/\//){ $items[0] =~ s/\s*\/\/.*$//; }
			$changed=1;
		}
		if($update_symbols && !$hashGenes{$symbol}){
			my $symb1 = $alias{$symbol};
			if($symb1){
				$symbol=$symb1;
				
			}
		}
		if($update_symbols && $hashGenes{$symbol}){
			$hashFound{$symbol}=1;
		}
		if($changed){ $lines->[$startLine] = join("\t",$probeID,$symbol,@items); }
	}
	$startLine++;
}
my $countAll = keys %hashAll;
my $n_found = keys %hashFound;
if(!$countAll){ error_message("No genes found in the data file."); }
my $prop = floor(1000*$n_found/$countAll+0.5)/1000;
my $file_info = "<b>Organism:</b>$hashOrganism{$organismID}<br>\n";
$file_info .= "<b>Total lines:</b> $startLine<br>\n";
$file_info .= "<b>N gene symbols submitted:</b> $countAll<br>\n";
if($update_symbols){
	$file_info .= "<b>N gene symbols accepted:</b> $n_found<br>\n";
	if($prop < 0.3){
		$file_info .= "<font size=+2><b>WARNING:</b></font>The proportion of genes accepted is too low! Verify gene symbols with NCBI database.<br>\n";
	}
}
return $file_info;
}

#**************************************
sub   check_output_data
#**************************************
{
my $lines = shift;
my $organismID = shift;
my $filename = shift;
my $description = shift;

my $startLine=0;
my %hashOutput;
while($startLine<@$lines && $lines->[$startLine] =~ /^\!/){
	my $line = $lines->[$startLine];
	$line =~ s/^\!//;
	$line =~ s/\"//g;
	if($line =~ /^Output/i){
		my ($key,$value) = split(/\t/,$line);
		$key=lc($key);
		$hashOutput{$key}=$value;
	}
	$startLine++;
}
if($lines->[$startLine-1] ne "!Matrix1_start"){
	splice(@$lines,$startLine++,0,"!Matrix1_start");
}
my ($junk,@headers) = split(/\t/,$lines->[$startLine++]);
my $nCol = @headers;
my $nRow = 0;
while($startLine<@$lines && $lines->[$startLine] !~ /^\!/){
	if(!$lines->[$startLine]){ splice(@$lines,$startLine--,1); }
	my($rowHeader,@rowValue) = split(/\t/,$lines->[$startLine]);
	if(@rowValue ne $nCol){ error_message("Wrong file format","continue"); }
	foreach my $x (@rowValue){
		unless($x=~/^-*[.0-9]+$/ || $x=~/^[-.0-9]+e[-+]*\d+$/i){ error_message("Wrong file format: $x","continue"); }
	}
	$nRow++;
	$startLine++;
}
if($nCol<1 || $nRow<1){ error_message("Wrong file format","continue"); }

if($lines->[$startLine-1] ne "!Matrix1_end"){
	splice(@$lines,$startLine++,0,"!Matrix1_end");
}
my $organismID1 = $hashOutput{"output_taxid"};
my $description1 = $hashOutput{"output_description"};
my $filename1 = $hashOutput{"output_filename"};
my $nCol1 = $hashOutput{"output_n_columns"};
my $nRow1 = $hashOutput{"output_n_rows"};
my $nMatrix1 = $hashOutput{"output_n_matrixes"};

$startLine=0;
if(!$filename1){
	splice(@$lines,$startLine++,0,"!Output_filename\t$filename");
}
if(!$description1 && $description){
	splice(@$lines,$startLine++,0,"!Output_description\t$description");
}
if($organismID1 && $organismID1 != $organismID){
	error_message("Organism ID do not match! First select a correct organism species in the main menu, then upload your file.","continue");
}
if(!$organismID1){
	splice(@$lines,$startLine++,0,"!Output_taxid\t$organismID");
}
if(!$nCol1){
	splice(@$lines,$startLine++,0,"!Output_N_columns\t$nCol");
}
if(!$nRow1){
	splice(@$lines,$startLine++,0,"!Output_N_rows\t$nRow");
}
if(!$nMatrix1){
	splice(@$lines,$startLine++,0,"!Output_N_matrixes\t1");
	splice(@$lines,$startLine++,0,"!Output_matrix1_name\ttemp");
	splice(@$lines,$startLine++,0,"!Output_matrix1_type\tnumbers");
}
return;
}

#**************************************
sub file_upload 
#**************************************
{
my $organismID = $hashInput{"organismID"};
my @matrix_list;
my @geneset_list;
my @samples_list;
my @annotation_list;
my @output_list;
my @annotation_list1;
my @annotation_GPL;
my %hashGPL;
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!","continue");
while(my $line = <INFO>){
	chop $line;
	my @items = split(/[=\t]/,$line);
	if($items[0] =~ /^type_matrix/){
		push(@matrix_list,$items[1]);
	}elsif($items[0] =~ /^type_geneset/){
		push(@geneset_list,$items[1]);
	}elsif($items[0] =~ /^type_samples/){
		push(@samples_list,$items[1]);
	}elsif($items[0] =~ /^type_annotation/){
		push(@annotation_list,$items[1]);
		my %hash=();
		read_config_line($line,\%hash);
		my $file_annotation = $hash{"type_annotation"};
		if($file_annotation !~ /^GPL\d+$/){
			push(@annotation_list1,[$file_annotation,$hash{"description"},$hash{"organismID"}]);
		}else{
			my $num=$file_annotation;
			$num =~ s/^GPL//;
			push(@annotation_GPL,[$file_annotation,$hash{"description"},$hash{"organismID"},$num]);
			$hashGPL{$num}=1;
		}
	}elsif($items[0] =~ /^type_output/){
		push(@output_list,$items[1]);
	}
}
close INFO;
@annotation_list1 = sort {lc($a->[0]) cmp lc($b->[0])} @annotation_list1;
push(@annotation_list1,["Gene symbols","Official gene symbols",$organismID]);
if($organismID < 20000 && $organismID > 3700 && $organismID != 4530 && $organismID != 4932 && $organismID != 5476 || $organismID==562 || $organismID==1280){
	push(@annotation_list1,["GenBank","GenBank or RefSeq accession",$organismID]);
	push(@annotation_list1,["Entrez","Entrez gene ID",$organismID]);
}
if($organismID < 10100 && $organismID > 9000 || $organismID==4932 || $organismID==6239 || $organismID==7227 || $organismID==7955){
	push(@annotation_list1,["Ensembl","Ensembl ID (gene or transcript)",$organismID]);
}
splice(@annotation_list1,0,0,["None","",$organismID]);
if($loginname ne "public" && open(INFO,"<$PATH_INFO/public-config.txt")){
	my @annotation_list2=();
	while(my $line = <INFO>){
		chop $line;
		my %hash=();
		read_config_line($line,\%hash);
		my $file_annotation = $hash{"type_annotation"};
		if($file_annotation){
			if($file_annotation !~ /^GPL\d+$/){
				push(@annotation_list2,["public-$file_annotation",$hash{"description"},$hash{"organismID"}]);
			}else{
				my $num=$file_annotation;
				$num =~ s/^GPL//;
				push(@annotation_GPL,["public-$file_annotation",$hash{"description"},$hash{"organismID"},$num]);
				$hashGPL{$num}=1;
			}
		}
	}
	push(@annotation_list1, sort {lc($a->[0]) cmp lc($b->[0])} @annotation_list2);
}
close INFO;
if(open(INFO,"<$PATH_DATA/array_platforms.txt")){
	while(my $line = <INFO>){
		chop $line;
		my($ID,$title,$technology,$taxonomy,$rows,$taxid) = split(/\t/,$line);
		my $found=0;
		foreach my $id (split(/;/,$taxid)){
			if($id==$organismID){ $found=1; }
		}
		if(!$found){ next; }
		my $num=$ID;
		$num =~ s/^GPL//;
		if($hashGPL{$num}){ next; }
		$hashGPL{$num} = 1;
		$title =~ s/\[[^\]]+\]//g;
		$title =~ s/\([^\)]+\)//g;
		if(length($title) > 30){ $title = substr($title,0,30); }
		push(@annotation_GPL,[$ID,$title,$organismID,$num]);
	}
	push(@annotation_list1, sort {$a->[3]<=>$b->[3]} @annotation_GPL);
}
filter_list_by_organism(\@annotation_list1, $organismID);
my $matrix_list="";
my $geneset_list="";
my $samples_list="";
my $output_list="";
my $annotation_list="";
foreach my $name (@matrix_list){
	if(!$matrix_list){ $matrix_list = "\"".$name."\""; }
	else{ $matrix_list .= ",\"".$name."\""; }
}
foreach my $name (@geneset_list){
	if(!$geneset_list){ $geneset_list = "\"".$name."\""; }
	else{ $geneset_list .= ",\"".$name."\""; }
}
foreach my $name (@samples_list){
	if(!$samples_list){ $samples_list = "\"".$name."\""; }
	else{ $samples_list .= ",\"".$name."\""; }
}
foreach my $name (@output_list){
	if(!$output_list){ $output_list = "\"".$name."\""; }
	else{ $output_list .= ",\"".$name."\""; }
}
foreach my $name (@annotation_list){
	if(!$annotation_list){ $annotation_list = "\"".$name."\""; }
	else{ $annotation_list .= ",\"".$name."\""; }
}
print "<HTML><HEAD><TITLE>ExAtlas</TITLE>\n";
print_header();
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "matrix_list = new Array($matrix_list);\n";
print "annotation_list = new Array($annotation_list);\n";
print "geneset_list = new Array($geneset_list);\n";
print "samples_list = new Array($samples_list);\n";
print "output_list = new Array($output_list);\n";
my ($items,$descriptions) = get_array_lists(\@annotation_list1);
print "platform_list = new Array($items);\n";
print "platform_description = new Array($descriptions);\n";
print "function return_to_main_menu(){\n";
#print "	document.upload.changeOrganism.value = 0;\n";
#print "	document.cancel.action.value = \"continue\";\n";
print "	document.cancel.submit();\n";
print "}\n";
print "function update_description() {\n";
print " var index;\n";
print "	document.upload.changeOrganism.value = 0;\n";
print " index = document.upload.platform.selectedIndex;\n";
print " document.upload.description_platform.value = platform_description[index];\n";
print "}\n";
print "function change_organism() {\n";
print "	document.upload.action.value = \"file_upload\";\n";
print "	document.upload.changeOrganism.value = 1;\n";
print "	document.upload.submit();\n";
print "}\n";
print "function upload_onsubmit() {\n";
print "	document.upload.changeOrganism.value = 0;\n";
print "	document.upload.action.value = \"continue\";\n";
print "	if(!document.upload.upload_filename.value && !document.upload.pasted_text.value){\n";
print "		alert(\"Nothing to upload! Select a file (Browse.. button) or paste text\"); return false;\n";
print "	}\n";
print "	if(document.upload.upload_filename.value && document.upload.pasted_text.value){\n";
print "		alert(\"You cannot use both text and file for submission\\nRemove one of them\"); return false;\n";
print "	}\n";
print "	if(document.upload.file_type.selectedIndex==3){\n";
print "		document.upload.action.value = \"gene_list_explore\";\n";
print "		document.upload.rename.value = \"\";\n";
print "		return true;\n";
print "	}\n";
print "	document.upload.action.value = \"continue\";\n";
print "	var file = document.upload.upload_filename.value;\n";
print "	file = file.replace(/.+[/\\\\]/,\"\");\n";
print "	if(document.upload.rename.value){\n";
print "		file=document.upload.rename.value;\n";
print "		var file1=file;\n";
print "		if(file.search(/\\.txt\$/)>=0){\n";
print "			if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "			file1=file.substring(0,file.length-4);\n";
print "		}\n";
print "		else if(file.search(/\\./)>=0){\n";
print "			alert(\"File extension is not valid. Only tab-delimited text files can be uploaded\"); return(false);\n";
print "		}\n";
print "		if(file1.search(/^[-\\w]+\$/)<0){\n";
print "			alert(\"File name should be one word with no special characters except underscore and dash.\\nPlease rename it!\");\n";
print "			return(false);\n";
print "		}\n";
print "	}\n";
print "	else if(document.upload.pasted_text.value){\n";
print "		alert(\"Fill in the file name in the field \'Rename file as:\'\"); return false;\n";
print "	}\n";
print "	if(!file){\n";
print "		alert(\"Nothing to upload! Select a file of paste a text first.\"); return false;\n";
print "	}\n";
print "	var i;\n";
print "	var filetype = document.upload.file_type.value;\n";
print "	if(filetype==\"annotation\" && file.search(/^[-\\w]+\$/)<0){\n";
print "		alert(\"File name should be one word with no dots and no special characters\\nRename it in the provided field\");\n";
print "		return(false);\n";
print "	}\n";
print "	if(filetype==\"matrix\" && file.search(/^GSE\\d/) >= 0){\n";
print "		alert(\"Do not upload matrix files from GEO; instead use button \'Find samples in GEO\' in main menu.\\nIf you need this specific file, rename it so that it does not start with \'GSE\'.\"); return false;\n";
print "	}\n";
print "	var legend_found = 0;\n";
print "	for(i=0; i<matrix_list.length; ++i){\n";
print "		if(filetype==\"matrix\" || filetype==\"expression\"){\n";
print "			if(file == matrix_list[i] && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "				return false;\n";
print "			}\n";
print "		}else if(file == matrix_list[i]){\n";
print "			alert(\"File with this name already exists\"); return false;\n";
print "		}else if(file == \"legend-\" + matrix_list[i]){\n";
print "			legend_found = 1;\n";
print "		}\n";
print "	}\n";
print "	if(filetype==\"legend\" && legend_found==0){\n";
print "		alert(\"Legend file for specific ExpressionProfile should be named as \'legend-ExpressionProfile\'\"); return false;\n";
print "	}\n";
print "	for(i=0; i<geneset_list.length; ++i){\n";
print "		if(filetype==\"geneset\"){\n";
print "			if(file == geneset_list[i] && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "				return false;\n";
print "			}\n";
print "		}else if(file == geneset_list[i]){\n";
print "			alert(\"File with this name already exists\"); return false;\n";
print "		}\n";
print "	}\n";
print "	for(i=0; i<samples_list.length; ++i){\n";
print "		if(filetype==\"samples\"){\n";
print "			if(file == samples_list[i] && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "				return false;\n";
print "			}\n";
print "		}else if(file == samples_list[i]){\n";
print "			alert(\"File with this name already exists\"); return false;\n";
print "		}\n";
print "	}\n";
print "	for(i=0; i<annotation_list.length; ++i){\n";
print "		if(filetype==\"annotation\"){\n";
print "			if((file==annotation_list[i] || file==annotation_list[i]+\"_annot.txt\") && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "				return false;\n";
print "			}\n";
print "		}else if(file==annotation_list[i]+\"_annot.txt\"){\n";
print "			alert(\"File with this name already exists\"); return false;\n";
print "		}\n";
print "	}\n";
print "	for(i=0; i<output_list.length; ++i){\n";
print "		if(filetype==\"output\"){\n";
print "			if(file == output_list[i] && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "				return false;\n";
print "			}\n";
print "		}else if(file == output_list[i]){\n";
print "			alert(\"File with this name already exists\"); return false;\n";
print "		}\n";
print "	}\n";
print "	var file1=file;\n";
print "	if(file.search(/\\.txt\$/)>=0){\n";
print "		if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "		file1=file.substring(0,file.length-4);\n";
print "	}\n";
print "	else if(file.search(/\\./)>=0){\n";
print "		alert(\"File extension is not valid. Only tab-delimited text files can be uploaded\"); return(false);\n";
print "	}\n";
print "	if(file1.search(/^[-\\w]+\$/)<0){\n";
print "		alert(\"File name should be one word with no special characters except underscore and dash.\\nRename it in provided field\");\n";
print "		return(false);\n";
print "	}\n";
print "	if(file.search(/^public-/) != -1){\n";
print "		alert(\"File name should not start with \'public\'. Rename it\");\n";
print "		return false;\n";
print "	}\n";
print "	if(file.search(/\\.(php|cgi|htm|html)\$/) != -1){\n";
print "		alert(\"Improper file name extension. Rename your data\");\n";
print "		return false;\n";
print "	}\n";
print "	if(file.search(/\\.xl/) != -1){\n";
print "		alert(\"Excel files cannot be uploaded. Convert to tab-delimited text and use txt extension.\");\n";
print "		return false;\n";
print "	}\n";
print "	var descrip = documect.upload.description.value;\n";
print "	if(descrip.search(/\\=|\\&/) >= 0){\n";
print "		alert(\"Description should not include character \'=\' or \'&\'\");\n";
print "		return false;\n";
print "	}\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print "<p><font size=+3>Upload a new file</font> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;\n";
print "<INPUT TYPE=button VALUE=\"Cancel\" style=width:120px; onClick=return_to_main_menu();><p>\n";
print "<FORM NAME=upload ENCTYPE=\"multipart/form-data\" ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST onSubmit=\"return upload_onsubmit();\">\n";
print "<b>Organism:</b> <select name=organismID style=width:260px; onChange=change_organism();>\n";
foreach my $ref (@organisms){
	print "<option value=$ref->[0]";
	if($ref->[0] == $organismID){ print " selected"; }
	print ">$ref->[3] ($ref->[2])\n";
}
print "</select><TD COLSPAN=2>(page is reloaded after change)<p>\n";
print "Paste text or browse files; check <a href=../exatlas-help.html#upload>file format</a> before uploading!<p>";
print "<TEXTAREA NAME=pasted_text ROWS=5 COLS=80></TEXTAREA><br>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD WIDTH=120>Or upload a file:<TD><INPUT NAME=upload_filename TYPE=file SIZE=30>\n";
print "<TD WIDTH=120 ALIGN=RIGHT>Select file type:<TD><SELECT NAME=file_type style=width:350px;>";
print "<OPTION VALUE=matrix>Gene expression profile matrix\n";
print "<OPTION VALUE=geneset>Gene set file\n";
print "<OPTION VALUE=samples>Samples file\n";
print "<OPTION VALUE=genelist>List of genes (symbols or GenBank acc#)\n";
print "<OPTION VALUE=output>Output file\n";
print "<OPTION VALUE=annotation>Annotation file\n";
print "<OPTION VALUE=expression>Compile expression profile\n";
if($loginname eq "public"){
	print "<OPTION VALUE=legend>Legend text file\n";
}
print "</SELECT>\n";
print "<TR><TD>&nbsp;\n";
print "<TR><TD COLSPAN=2><b><font size=+1>Optional feilds:</font></b>\n";
print "<TR><TD>Rename file as:<TD><INPUT NAME=\"rename\" style=width:220px;>\n";
print "<TD ALIGN=right>Description:<TD><INPUT NAME=\"description\" style=width:350px;>\n";
print "<TR><TD>Select platform:<TD><SELECT NAME=platform style=width:220px; onChange=update_description();>\n";
for(my $i=0; $i<@annotation_list1; ++$i){ 
	my $name1 = $annotation_list1[$i]->[0];
	if($name1 =~ /^public-GPL/){ $name1 =~ s/^public-//; }
	if($name1 =~ /^Gene symbols/){
		print "<option value=symbol_$organismID> $name1\n";
	}elsif($name1 =~ /^GenBank$/){
		print "<option value=genbank_$organismID> $name1\n";
	}elsif($name1 =~ /^Ensembl$/){
		print "<option value=ensembl_$organismID> $name1\n";
	}elsif($name1 =~ /^Entrez$/){
		print "<option value=entrez_$organismID> $name1\n";
	}else{
		print "<option value=\"$name1\"> $name1\n";
	}
}
print "</select><TD ALIGN=right>Description:<TD ALIGN=right><INPUT NAME=\"description_platform\" style=width:350px;><TD>\n";
print "<TR><TD>Transformation:<TD><SELECT NAME=transform style=width:220px;>";
print "<OPTION VALUE=0>None";
print "<OPTION VALUE=2>Log2";
print "<OPTION VALUE=10>Log10";
print "<OPTION VALUE=2.71828>Loge";
print "<OPTION VALUE=10>Z-value";
print "</SELECT>\n";
print "<TD ALIGN=right>Adjustment(subtract):<TD><SELECT NAME=adjustment style=width:80px;>";
print "<OPTION VALUE=0>None";
print "<OPTION VALUE=50>50";
print "<OPTION VALUE=100>100";
print "</SELECT>\n";
print "&nbsp; &nbsp; <INPUT TYPE=CHECKBOX NAME=update_symbols CHECKED> Update gene symbols\n";
print "</TABLE><br>\n";
print "Look-up table of GEO array platforms: <a href=../download/array_platforms.txt>platforms</a><br>\n";
print "Select platform if expression profile matrix file has no platform information<p>\n";
print "<INPUT TYPE=SUBMIT value=\"Upload file\" style=width:120px;><br>\n";
print "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=hidden VALUE=\"continue\">\n";
print "<INPUT NAME=\"changeOrganism\" TYPE=hidden>\n";
#print "<INPUT NAME=\"organismID\" TYPE=hidden VALUE=\"$organismID\">\n";
print "</FORM><p>\n";
print "<b>Note 1:</b> If microarray platform is not available, then upload custom array annotation before uploading gene expression data.<br>\n";
print "<b>Note 2:</b> If microarray platform from GEO is not in the list, you can try to specify it in the header of expression profile matrix file (see <a href=../exatlas-help.html#upload>help</a>).<br>\n";
print "<b>Note 3:</b> Gene expression profile matrix may use raw data, log10, log2, or Z-value transformed; log-ratio is not supported<br>\n";
print "<b>Note 4:</b> To upload a list of genes, use gene symbols (one at each line or comma-separated).<br>\n";
print "<b>Note 5:</b> If the file is too big, you can upload it from Helix FTP server (helix.nih.gov). First load the file to the server (make your sub-directory within pub directory).\n";
print "Then, specify the address on the paste box: ftp://helix.nih.gov/pub/directory/filename<br>\n";
print "<b>Note 6:</b> Option 'Compile expression profile' is used for compiling gene expression data from\n";
print "multiple input files generated by scanning microarrays or processing RNA-seq data (e.g. Cufflinks file genes.fpkm_tracking).\n";
print "The file should be a tab-delimited text with column headers (which can be renamed later). For microarrays,\n";
print "select platform from pull-down menu or upload a custom array annotation which then will appear in the\n";
print "platform selection menu. If a data file includes a column with gene identifiers (e.g., gene symbol, RefSeq or\n";
print "GenBank acc#, Ensembl gene, Entrez ID) then it can be used for analysis without specifying platform.\n";
print "When the data is uploaded you can save it as Gene Expression Profile matrix.<p>\n";
print "<FORM NAME=cancel ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
print "<INPUT NAME=action TYPE=hidden VALUE=continue>\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print "<INPUT TYPE=submit VALUE=\"Cancel\" style=width:120px;>\n";
print "</FORM>\n";
return;
}

#**************************************
sub file_edit 
#**************************************
{
my $organismID = $hashInput{"organismID"};
my $filename = $hashInput{"file_edit"};
my $file_type = $hashInput{"file_type"};
if(!$filename){
	$filename = $hashInput{"filename"};
}
my $description;
my @file_list;
my @copy_file_list;
my @geneset_list;
my @annotation_list1;
my $organismID_old;
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!","continue");
while(my $line = <INFO>){
	chop $line;
	my %hash;
	read_config_line($line,\%hash);
	my $descr = $hash{"description"};
	my $org = $hash{"organismID"};
	my $date1 = $hash{"date"};
	if($file_type eq "matrix" && $line=~/^type_matrix/ && $filename ne $hash{"type_matrix"}){
		push(@copy_file_list,[$hash{"type_matrix"},$descr,$org,$date1]);
	}elsif($file_type eq "geneset" && $line=~/^type_geneset/ && $filename ne $hash{"type_geneset"}){
		push(@copy_file_list,[$hash{"type_geneset"},$descr,$org,$date1]);
	}elsif($file_type eq "matrix" && $hash{"type_annotation"}){
		push(@annotation_list1,[$hash{"type_annotation"},$hash{"description"},$hash{"organismID"}]);
	}
	if($hash{"type_$file_type"} eq $filename){
		$description = $descr;
		$organismID_old = $org;
	}
	my @items = split(/[=\t]/,$line);
	if($items[0] =~ /^type_/){
		push(@file_list,$items[1]);
	}
}
close INFO;
@annotation_list1 = sort {lc($a->[0]) cmp lc($b->[0])} @annotation_list1;
push(@annotation_list1,["Gene symbols","Official gene symbols",$organismID]);
if($organismID < 20000 && $organismID > 3700 && $organismID != 4530 && $organismID != 4932 && $organismID != 5476 || $organismID==562 || $organismID==1280){
	push(@annotation_list1,["GenBank","GenBank or RefSeq accession",$organismID]);
	push(@annotation_list1,["Entrez","Entrez gene ID",$organismID]);
}
if($organismID < 10100 && $organismID > 9000 || $organismID==4932 || $organismID==6239 || $organismID==7227 || $organismID==7955){
	push(@annotation_list1,["Ensembl","Ensembl ID (gene or transcript)",$organismID]);
}
splice(@annotation_list1,0,0,["None","",$organismID]);
if($loginname ne "public" && open(INFO,"<$PATH_INFO/public-config.txt")){
	my @annotation_list2=();
	while(my $line = <INFO>){
		chop $line;
		my %hash=();
		read_config_line($line,\%hash);
		my $file_annotation = $hash{"type_annotation"};
		if($file_annotation && $file_annotation !~ /^GPL\d/){
			push(@annotation_list2,["public-$file_annotation",$hash{"description"},$hash{"organismID"}]);
		}
	}
	push(@annotation_list1, sort {lc($a->[0]) cmp lc($b->[0])} @annotation_list2);
}
close INFO;
if(open(INFO,"<$PATH_DATA/array_platforms.txt")){
	my @annotation_list2=();
	while(my $line = <INFO>){
		chop $line;
		my($ID,$title,$technology,$taxonomy,$rows,$taxid) = split(/\t/,$line);
		my $found=0;
		foreach my $id (split(/;/,$taxid)){
			if($id==$organismID){ $found=1; }
		}
		if(!$found){ next; }
		my $num=$ID;
		$num =~ s/^GPL//;
		$title =~ s/\[[^\]]+\]//g;
		$title =~ s/\([^\)]+\)//g;
		if(length($title) > 30){ $title = substr($title,0,30); }
		push(@annotation_list2,[$ID,$title,$organismID,$num]);
	}
	push(@annotation_list1, sort {$a->[3]<=>$b->[3]} @annotation_list2);
}
filter_list_by_organism(\@annotation_list1, $organismID);
@copy_file_list = sort {lc($a->[0]) cmp lc($b->[0])} @copy_file_list;
filter_list_by_organism(\@copy_file_list, $organismID);
splice(@copy_file_list,0,0,["--- New file ---",""]);
my $file_list="";
foreach my $name (@file_list){
	if(!$file_list){ $file_list = "\"".$name."\""; }
	else{ $file_list .= ",\"".$name."\""; }
}
my $filename_full = "$loginname-$filename";
if($filename=~ /^public-/){
	$filename_full = $filename;
	my $filename1 = $filename;
	$filename1 =~ s/^public-//;
	open(INFO,"<$PATH_INFO/public-config.txt") or error_message("Configuration file not found!","continue");
	while(my $line = <INFO>){
		chop $line;
		my %hash;
		read_config_line($line,\%hash);
		if($hash{"type_$file_type"} eq $filename1){
			$description = $hash{"description"};
			last;
		}
	}
	close INFO;
}
my %hashFile;
my @item_names;
my @item_descriptions;
my @geneset_attributes;
my @geneset_attrib_min;
my @geneset_attrib_max;
my $item_name;
if($file_type eq "matrix"){
	parse_matrix_file("$PATH_DATA/$filename_full", \%hashFile, "continue");
	my $ref = $hashFile{"sample_title"};
	if($ref && ref($ref) eq 'ARRAY'){
		@item_names = @$ref;
	}
	$item_name = "Sample";
}
elsif($file_type eq "geneset"){
	parse_geneset_file("$PATH_DATA/$filename_full", \%hashFile, "continue");
	my $ref = $hashFile{"geneset_names"};
	if($ref && ref($ref) eq 'ARRAY'){
		@item_names = @$ref;
	}
	$ref = $hashFile{"geneset_descriptions"};
	if($ref && ref($ref) eq 'ARRAY'){
		@item_descriptions = @$ref;
	}
	$ref = $hashFile{"geneset_attributes"};
	if($ref && ref($ref) eq 'ARRAY'){
		@geneset_attributes = @$ref;
	}
	$ref = $hashFile{"geneset_attribute_min"};
	if($ref && ref($ref) eq 'ARRAY'){
		@geneset_attrib_min = @$ref;
	}
	$ref = $hashFile{"geneset_attribute_max"};
	if($ref && ref($ref) eq 'ARRAY'){
		@geneset_attrib_max = @$ref;
	}
	$item_name = "Gene set";
}
print "<HTML><HEAD><TITLE>ExAtlas</TITLE>\n";
my $command;
if(@item_names){ $command="update_description();" }
print_header($command);
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "file_list = new Array($file_list);\n";
my ($items,$descriptions) = get_array_lists(\@copy_file_list);
print "copy_file_list = new Array($items);\n";
print "copy_file_description = new Array($descriptions);\n";
print "function save_changes() {\n";
print "	if(!document.file_edit.new_filename.value){\n";
print "		alert(\"Empty file name! Type it in or cancel.\"); return false;\n";
print "	}\n";
print "	var file = file_edit.new_filename.value;\n";
print "	var file1=file;\n";
print "	if(file.search(/\\.txt\$/)>=0){\n";
print "		if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "		file1=file.substring(0,file.length-4);\n";
print "	}\n";
print "	if(file1.search(/^[-\\w]+\$/)<0){\n";
print "		alert(\"File name should have neither spaces nor special characters!\");\n";
print "		return false;\n";
print "	}\n";
print "	if(file.search(/^public-/i) >= 0){\n";
print "		alert(\"File name cannot start with \'public-\'!\");\n";
print "		return false;\n";
print "	}\n";
print "	var i;\n";
print "	for(i=0; i<file_list.length; ++i){\n";
print "		if(file != \"$filename\" && file == file_list[i]){\n";
print "			alert(\"This file name already exists! Rename it or cancel.\"); return false;\n";
print "		}\n";
print "	}\n";
print "	var descrip = document.file_edit.new_description.value;\n";
print "	if(descrip.search(/\\=|\\&/) >= 0){\n";
print "		alert(\"Description should not include character \'=\' or \'&\'\");\n";
print "		return false;\n";
print "	}\n";
print "	if(!check_attribute_limits()){ return false; }\n";
print "	document.file_edit.selected_items.value = \"\";\n";
print "	document.file_edit.action.value = \"file_edit1\";\n";
print "	document.file_edit.submit();\n";
print "}\n";
if(@item_names){
	print "function update_description() {\n";
	print "  var index;\n";
	print "  index = document.file_edit.copy_file.selectedIndex;\n";
	print "  document.file_edit.description_copy_file.value = copy_file_description[index];\n";
	print "}\n";
}
print "function count_checked_samples(){\n";
print "	var nchecked = 0;\n";
print "	for(i=0; i<document.file_edit.elements.length; ++i){\n";
print "		if(document.file_edit.elements[i].type==\"checkbox\"){\n";
print "			if(document.file_edit.elements[i].checked){ nchecked++; };\n";
print "		}\n";
print "	}\n";
print "	return nchecked;\n";
print "}\n";
print "function check_all(){\n";
print "	for(i=0; i<document.file_edit.elements.length; ++i){\n";
print "		if(document.file_edit.elements[i].type==\"checkbox\"){\n";
print "			document.file_edit.elements[i].checked=true;\n";
print "		}\n";
print "	}\n";
print "}\n";
print "function uncheck_all(){\n";
print "	for(i=0; i<document.file_edit.elements.length; ++i){\n";
print "		if(document.file_edit.elements[i].type==\"checkbox\"){\n";
print "			document.file_edit.elements[i].checked=false;\n";
print "		}\n";
print "	}\n";
print "}\n";
print "function check_attribute_limits(){\n";
for(my $i=0; $i<@geneset_attributes; $i++){
	print "	var xmin = new Number(document.file_edit.geneset_attribute_min$i.value);\n";
	print "	var xmax = new Number(document.file_edit.geneset_attribute_max$i.value);\n";
	print "	if(xmin > xmax){\n";
	print "		alert(\"Wrong limits for $geneset_attributes[$i] (max < min)! Cannot proceed.\"); return false;\n";
	print "	}\n";
}
print "	return true;\n";
print "}\n";
print "function copy_items() {\n";
print "	if(count_checked_samples()==0){\n";
print "		alert(\"No items selected! Cannot proceed.\"); return false;\n";
print "	}\n";
print " var file = document.file_edit.copy_file.options[document.file_edit.copy_file.selectedIndex].value;\n";
print " if(document.file_edit.copy_file.selectedIndex==0){\n";
print "		file = prompt(\"Provide name of a new file where to copy selected items\");\n";
print "		if(!file){ return(false); }\n";
print "		var file1=file;\n";
print "		if(file.search(/\\.txt\$/)>=0){\n";
print "			if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "			file1=file.substring(0,file.length-4);\n";
print "		}\n";
print "		if(file1.search(/^[-\\w]+\$/)<0){\n";
print "			alert(\"File name should have neither spaces nor special characters\");\n";
print "			return(false);\n";
print "		}\n";
print "		if(file.search(/^public/i) >= 0){\n";
print "			alert(\"File name cannot start with 'public'\");\n";
print "			return(false);\n";
print "		}\n";
print "		for(i=0; i<file_list.length; ++i){\n";
print "			if(file == file_list[i]){\n";
print "				alert(\"File with this name already exists\"); return false;\n";
print "			}\n";
print "		}\n";
print "		var descrip = document.file_edit.description_copy_file.value;\n";
print "		if(descrip.search(/\\=|\\&/) >= 0){\n";
print "			alert(\"Description should not include character \'=\' or \'&\'\");\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	if(!check_attribute_limits()){ return false; }\n";
print "	document.file_edit.selected_items.value = \"copy-\"+file;\n";
print "	document.file_edit.action.value = \"file_edit1\";\n";
print "	document.file_edit.submit();\n";
print "}\n";
print "function delete_items() {\n";
print "	if(count_checked_samples()==0){\n";
print "		alert(\"No samples selected! Cannot proceed.\"); return false;\n";
print "	}\n";
print "	document.file_edit.selected_items.value = \"delete\";\n";
print "	document.file_edit.action.value = \"file_edit1\";\n";
print "	document.file_edit.submit();\n";
print "}\n";
print "function cancel_update() {\n";
print "	document.file_edit.action.value = \"continue\";\n";
print "	document.file_edit.selected_items.value = \"\";\n";
print "	document.file_edit.submit();\n";
print "}\n";
print "function change_organism() {\n";
print "	document.file_edit.action.value = \"file_edit\";\n";
print "	document.file_edit.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print "<FORM NAME=file_edit ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST onSubmit=\"return upload_onsubmit();\">\n";
print "<p><FONT SIZE=6><b>Edit file attributes</b></FONT>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;\n";
print "<INPUT TYPE=button VALUE=\"Cancel & return to main menu\" onClick=cancel_update(); style=width:260px;><p>\n";
if($filename=~/^public-/){
	print "<b>Note:</b> Public files cannot be edited, but data can be copied to your own files.<p>\n";
}
print "<b>File type:</b> $file_type<p>\n";
print "<b>File name and description</b><p>\n";
print "&nbsp; &nbsp; &nbsp; <INPUT NAME=new_filename VALUE=\"$filename\" style=width:360px;> File name<br>\n";
print "&nbsp; &nbsp; &nbsp; <INPUT NAME=new_description VALUE=\"$description\" style=width:360px;> Description<br>\n";
print "&nbsp; &nbsp; &nbsp; <SELECT NAME=organismID style=width:360px; onChange=change_organism();>\n";
foreach my $ref (@organisms){
	print "<option value=$ref->[0]";
	if($ref->[0] == $organismID){ print " selected"; }
	print ">$ref->[3] ($ref->[2])\n";
}
print "</select> Organism (page is reloaded after change)<br>\n";
if($file_type eq "matrix"){
	my $platform = $hashFile{"series_platform_id"};
	print "&nbsp; &nbsp; &nbsp; <SELECT NAME=platform_new style=width:360px;>\n";
	for(my $i=0; $i<@annotation_list1; ++$i){ 
		my $name1 = $annotation_list1[$i]->[0];
		my $selected = "";
		if($name1 eq $platform || $name1 eq "public-".$platform){ $selected = " selected"; }
		if($name1 =~ /^public-GPL/){ $name1 =~ s/^public-//; }
		if($name1 =~ /^Gene symbols/){
			if($platform eq "symbol_$organismID" || $platform==$organismID){ $selected = " selected"; }
			print "<option value=\"symbol_$organismID\" $selected> $name1\n";
		}elsif($name1 =~ /^GenBank$/){
			if($platform eq "genbank_$organismID"){ $selected = " selected"; }
			print "<option value=\"genbank_$organismID\"$selected> $name1\n";
		}elsif($name1 =~ /^Ensembl$/){
			if($platform eq "ensembl_$organismID"){ $selected = " selected"; }
			print "<option value=\"ensembl_$organismID\"$selected> $name1\n";
		}elsif($name1 =~ /^Entrez$/){
			if($platform eq "entrez_$organismID"){ $selected = " selected"; }
			print "<option value=\"ensembl_$organismID\"$selected> $name1\n";
		}else{
			print "<option value=\"$name1\"$selected> $name1\n";
		}
	}
	print "</select> Gene expression platform<br>\n";
}
print "<p>";
my $name_lc = lc($item_name);
if(@item_names){
	if($filename !~ /^public-/){
		print "You can edit item titles/descriptions or delete/copy selected items.<br>";
	}
	if($file_type eq "matrix" && $filename !~ /^public-/){
		print "Note that samples with identical descriptions are treated as replications.<br>\n"
	}
	print "<p><INPUT TYPE=button VALUE=\"Check all items\" LANGUAGE=\"javascript\" onClick=check_all(); style=width:140px;>\n";
	print "<INPUT TYPE=button VALUE=\"Uncheck all items\" LANGUAGE=\"javascript\" onClick=uncheck_all(); style=width:140px;>\n";
	print "<h3>Edit $name_lc titles</h3>\n";
	print "<TABLE BORDER=0>\n";
	print "<TR><TD>No.<TD>Select<TD>$item_name title";
	if(@item_descriptions){ print "<TD>$item_name description"; }
	for(my $i=0; $i<@item_names; $i++){
		my $i1=$i+1;
		print "\n<TR><TD>$i1.";
		print "<TD><INPUT NAME=item_select$i TYPE=CHECKBOX>";
		print "<TD><INPUT NAME=item_name$i SIZE=40 VALUE=\"$item_names[$i]\">";
		if(@item_descriptions){
			print "<TD><INPUT NAME=item_description$i SIZE=40 VALUE=\"$item_descriptions[$i]\">";
		}
	}
	print "\n</TABLE>";
}
if(@geneset_attributes){
	print "<h3>Geneset attributes</h3>\n";
	if($filename !~ /^public-/){
		print "You can change attribute limits to filter genesets (original or copies)<br>";
	}
	print "<TABLE BORDER=0>\n";
	print "<TR><TD><b>No.<TD><b>Attribute<TD ALIGN=center><b>Minimum<TD ALIGN=center><b>Maximum\n";
	for(my $i=0; $i<@geneset_attributes; $i++){
		my $i1=$i+1;
		print "<TR><TD>$i1.<TD>$geneset_attributes[$i]\n";
		print "<TD><INPUT NAME=geneset_attribute_min$i VALUE=$geneset_attrib_min[$i] SIZE=10>\n";
		print "<TD><INPUT NAME=geneset_attribute_max$i VALUE=$geneset_attrib_max[$i] SIZE=10>\n";
	}
	print "\n</TABLE>";
}
print "<p><TABLE BORDER=0>";
if($filename !~ /^public-/){
	print "<TR><TD><INPUT TYPE=button value=\"Save changes\" onClick=save_changes(); style=width:250px;>\n";
}
if(@item_names){
	print "<TR><TD><INPUT TYPE=button VALUE=\"Copy selected items\" onClick=copy_items(); style=width:250px;><TD> to file:\n";
	print "<select name=copy_file onChange=update_description(); style=width:220px;>\n";
	print "<option value=0> $copy_file_list[0]->[0]\n";
	for(my $i=1; $i<@copy_file_list; ++$i){ 
		print "<option value=\"$copy_file_list[$i]->[0]\"> $copy_file_list[$i]->[0]\n";
	}
	print "</select>\n";
	print "&nbsp; &nbsp; &nbsp;Description: <INPUT NAME=description_copy_file style=width:220px;>\n";
	if($file_type eq "matrix"){
		print "<TR><TD><TD>Copy option: <select name=copy_option style=width:220px;>\n";
		print "<option value=all> Use all genes/probes\n";
		print "<option value=common> Use common genes/probes\n";
		print "<option value=original> Keep destination genes\n";
		print "</select>\n";
	}
	if($filename !~ /^public-/){
		print "<TR><TD><INPUT TYPE=button VALUE=\"Delete selected items\" onClick=delete_items(); style=width:250px;>\n";
	}
}
print "<TR><TD><INPUT TYPE=button VALUE=\"Cancel & return to main menu\" onClick=cancel_update(); style=width:250px;>\n";
print "</TABLE><br>\n";
print "If you select \"New file\" you will be prompted for file name, which shound be one-word without special characters (underscore allowed)<p>\n";

print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=action TYPE=hidden>\n";
print "<INPUT NAME=selected_items TYPE=hidden>\n";
print "<INPUT NAME=filename TYPE=hidden VALUE=\"$filename\">\n";
print "<INPUT NAME=file_type TYPE=hidden VALUE=\"$file_type\">\n";
print "<INPUT NAME=organismID_old TYPE=hidden VALUE=$organismID_old\">\n";
print "</FORM>\n";
print "</BODY></HTML>\n";
return;
}

#**************************************
sub file_edit1
#**************************************
{
my $organismID = $hashInput{"organismID"};
my $organismID_old = $hashInput{"organismID_old"};
my $filename = $hashInput{"filename"};
my $file_type = $hashInput{"file_type"};
my $new_description = $hashInput{"new_description"};
my $new_filename = $hashInput{"new_filename"};
my $selected_items = $hashInput{"selected_items"};

my $filename_full = "$PATH_DATA/$loginname-$filename";
if($filename=~ /^public-/){ $filename_full = "$PATH_DATA/$filename"; }
my @item_names=();
my @geneset_attributes;
my @geneset_attrib_min;
my @geneset_attrib_max;
my $geneset_filtered=0;
my %hashFile;
if($file_type eq "matrix"){
	parse_matrix_file($filename_full, \%hashFile,"continue");
	my $ref = $hashFile{"sample_title"};
	if($ref && ref($ref) eq 'ARRAY'){
		@item_names = @$ref;
	}
}elsif($file_type eq "geneset"){
	parse_geneset_file($filename_full, \%hashFile, "continue");
	my $ref = $hashFile{"geneset_names"};
	if($ref && ref($ref) eq 'ARRAY'){
		@item_names = @$ref;
	}
	$ref = $hashFile{"geneset_attributes"};
	if($ref && ref($ref) eq 'ARRAY'){
		@geneset_attributes = @$ref;
	}
	$ref = $hashFile{"geneset_attribute_min"};
	if($ref && ref($ref) eq 'ARRAY'){
		@geneset_attrib_min = @$ref;
	}
	$ref = $hashFile{"geneset_attribute_max"};
	if($ref && ref($ref) eq 'ARRAY'){
		@geneset_attrib_max = @$ref;
	}
	for(my $i=0; $i<@geneset_attributes; $i++){
		my $x = $hashInput{"geneset_attribute_min$i"};
		my $y = $hashInput{"geneset_attribute_max$i"};
		if($x != $geneset_attrib_min[$i] || $y != $geneset_attrib_max[$i]){
			$geneset_attrib_min[$i]=$x;
			$geneset_attrib_max[$i]=$y;
			$geneset_filtered=1;
		}
	}
}
if($selected_items =~ /^copy-/){
	my $file_copy = $selected_items;
	$file_copy =~ s/^copy-//;
	if($file_copy =~ /^public/){ error_message("Filename starts with -public-","continue"); }
	my $file_copy_full = "$PATH_DATA/$loginname-$file_copy";
	my $file_copy_description = $hashInput{"description_copy_file"};
	my @item_select=();
	my $count=0;
	for(my $i=0; $i<@item_names; $i++){
		if($hashInput{"item_select$i"} eq "on"){
			$item_select[$i] = 1;
			$count++;
		}
	}
	my @data_geneset;
	my $Ncopied;
	if($hashInput{"copy_file"} =~ /-\s*New\s*file\s*-/){
		if($count==@item_names && !$geneset_filtered){
			copy "$filename_full", "$file_copy_full";
			$Ncopied = $count;
		}else{
			open(INFO,'<',$filename_full) or error_message("Cannot open file","continue");
			open(OUT, ">$file_copy_full") or error_message("Cannot write to file","continue");
			if($file_type eq "matrix"){
				$Ncopied = $count;
				while(my $line = <INFO>){
					chop $line;
					my ($title,@data) = split(/\t/,$line);
					if(@data < 3){
						print OUT "$line\n";
					}else{
						my @list = apply_mask(\@data,\@item_select);
						if($list[0]<=$MISSING && median(\@list)==$MISSING){ next; }
						print OUT "$title\t".join("\t",@list)."\n";
					}
				}
				my $file_anova = "$loginname-anova-$file_copy";
				if(file_exist("$PATH_DATA/$file_anova")){
					unlink "$PATH_DATA/$file_anova";
				}
			}elsif($file_type eq "geneset"){
				my $count=-1;
				while(my $line = <INFO>){
					chop $line;
					if($line =~ /^\!/ || !$line){
						print OUT $line."\n";
						next;
					}
					my ($title,$description,@items)=split(/\t/,$line);
					if($title){
						$count++;
						if(@data_geneset){
							$Ncopied+=filter_geneset(\@data_geneset,\@geneset_attributes,\@geneset_attrib_min,\@geneset_attrib_max);
							@data_geneset=();
						}
					}
					if($item_select[$count]){
						if($geneset_filtered){
							push(@data_geneset,$line);
						}else{
							if($title){ $Ncopied++; }
							print OUT $line."\n";
						}
					}
				}
			}else{
				error_message("Unexpected file type $file_type","continue");				
			}
			close INFO;
			if(@data_geneset){
				$Ncopied+=filter_geneset(\@data_geneset,\@geneset_attributes,\@geneset_attrib_min,\@geneset_attrib_max);
				@data_geneset=();
			}
			close OUT;
		}
		if(!$Ncopied){
			unlink "$file_copy_full";
			error_message("No data fits filtering criteria","continue");
		}
		#Update configuration file
		open(OUT, ">$PATH_INFO/$loginname-config1.txt");
		open(INFO,"<$PATH_INFO/$loginname-config.txt");
		my $nLines;
		while(my $line = <INFO>){
			$line =~ s/\n$//;
			if(!$line){ next; }
			$nLines++;
			if($line =~ /^type_$file_type=$file_copy\s/){
				if($file_type eq "matrix" && file_exist("$PATH_DATA/$loginname-anova-$file_copy")){
					unlink "$PATH_DATA/$loginname-anova-$file_copy";
				}
			}else{
				print OUT $line."\n";
			}
		}
		close INFO;
		print OUT "type_$file_type=$file_copy\torganismID=$organismID";
		if($file_copy_description){ print OUT "\tdescription=$file_copy_description"; }
		print OUT "\tdate=$date_record\n";
		close OUT;
		my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
		if($nLines1 && $nLines1 > $nLines*0.9){
			copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
		}else{
			error_message("Failed to update configuration file!");
		}
		unlink("$PATH_INFO/$loginname-config1.txt");
	}
	else{   #Combine files
		if($file_type eq "geneset"){
			my %exist;
			open(INFO,'<',$file_copy_full) or error_message("Cannot read file","continue");
			while(my $line = <INFO>){
				chop $line;
				if($line =~ /^\!/ || !$line){ next; }
				my ($title,$description)=split(/\t/,$line);
				$exist{$title}=1;
			}
			close INFO;
			open(INFO,'<',$filename_full) or error_message("Cannot open file","continue");
			open(OUT, ">>$file_copy_full") or error_message("Cannot write to file","continue");
			my $count=-1;
			my $title1;
			while(my $line = <INFO>){
				chop $line;
				if($line =~ /^\!/ || !$line){
					next;
				}
				my ($title,$description,@items)=split(/\t/,$line);
				if($title){
					$title1=$title;
					$count++;
					if(@data_geneset){
						$Ncopied+=filter_geneset(\@data_geneset,\@geneset_attributes,\@geneset_attrib_min,\@geneset_attrib_max);
						@data_geneset=();
					}
				}
				if($item_select[$count] && !$exist{$title1}){
					if($geneset_filtered){
						push(@data_geneset,$line);
					}else{
						if($title){ $Ncopied++; }
						print OUT $line."\n";
					}
				}
			}
			close INFO;
			if($geneset_filtered && @data_geneset){
				$Ncopied+=filter_geneset(\@data_geneset,\@geneset_attributes,\@geneset_attrib_min,\@geneset_attrib_max);
				@data_geneset=();
			}
			close OUT;
			if(!$Ncopied){
				error_message("No data fits filtering criteria, copying cancelled","continue");
			}
			return;
		}
		#Combine matrix files
		my %hashMatrix1=();
		parse_matrix_file($file_copy_full,\%hashMatrix1);
		my $platform1 = $hashMatrix1{"series_platform_id"};
		my $platform2 = $hashFile{"series_platform_id"};
		my $platform;
		if($platform1 ne $platform2){ $platform = "symbol_$organismID"; }
		else{ $platform = $platform1; }
		my $ref = $hashMatrix1{"sample_title"};
		if(!$ref || ref($ref) ne 'ARRAY'){ error_message("No column headres","continue"); }
		my @headers1 = @$ref;
		my @headers2 = apply_mask(\@item_names,\@item_select);
		my %hash_annotation;
		if($platform1 ne $platform){
			if(!read_platform_annotation($platform1,\%hash_annotation)){
				error_message("Cannot read annotation.","continue");
			}
		}
		open (INFO,'<',$file_copy_full) or error_message("File1 not open","continue");
		my $header_text;
		while(my $line = <INFO>){
			if($line=~/^!Series_platform_id/){
				$line = "!Series_platform_id\t\"$platform\"\n";
			}elsif($line=~/^!Sample_title/){
				$line = "!Sample_title\t\"".join("\"\t\"",@headers1,@headers2)."\"\n";
			}elsif($line=~/^!Series_sample_id|^!Sample_data_row_count/){
				next;
			}
			$header_text .= $line;
			if($line =~ /^\!series_matrix_table_begin/i){ last; }
		}
		my $line = <INFO>; # skip_headers
		my %hashProbeID;
		my @probeID;
		my @dataAll1;
		my $irow=0;
		while(my $line = <INFO>){
			if($line =~ /^!/){ last; }
			chop $line;
			$line =~ s/\s+$//;
			my ($probe_id,@data) = split(/\t/, $line);
			$probe_id =~ s/\"//g;
			if($platform1 ne $platform){ $probe_id = $hash_annotation{$probe_id}; }
			if(!$probe_id || $probe_id=~/^none$/i){ next; }
			my $ref = $hashProbeID{$probe_id};
			if(!$ref){
				$hashProbeID{$probe_id} = [$irow++,median(\@data),0];
				push(@probeID,$probe_id);
				push(@dataAll1,\@data);
			}else{
				my $median = median(\@data);
				if($median > $ref->[1]){
					@{$dataAll1[$ref->[0]]} = @data;
					$ref->[1] = $median;
				}
			}
		}
		close INFO;
		if($platform2 ne $platform){
			if(!read_platform_annotation($platform2,\%hash_annotation)){
				error_message("Cannot read annotation","continue");
			}
		}
		my $copy_option = $hashInput{"copy_option"};
		open (INFO,'<',$filename_full) or error_message("File2 not open","continue");
		while(my $line = <INFO>){
			if($line =~ /^\!series_matrix_table_begin/i){ last; }
		}
		$line = <INFO>; # skip_headers
		my @blank1;
		foreach my $x (@headers1){ push(@blank1,$MISSING); }
		my @dataAll2;
		while(my $line = <INFO>){
			if($line =~ /^!/){ last; }
			chop $line;
			$line =~ s/\s+$//;
			my ($probe_id,@data1) = split(/\t/, $line);
			$probe_id =~ s/\"//g;
			my @data = apply_mask(\@data1,\@item_select);
			if($data[0]<=$MISSING && median(\@data)==$MISSING){ next; }
			if($platform2 ne $platform){ $probe_id = $hash_annotation{$probe_id}; }
			if(!$probe_id || $probe_id=~/^none$/i){ next; }
			my $ref = $hashProbeID{$probe_id};
			if(!$ref){
				if($copy_option eq "all"){
					$dataAll2[$irow] = \@data;
					$dataAll1[$irow] = \@blank1;
					$hashProbeID{$probe_id} = [$irow++,0,median(\@data)];
					push(@probeID,$probe_id);
				}
			}else{
				my $median = median(\@data);
				if($median > $ref->[2]){
					my $ii = $ref->[0];
					if(!$dataAll2[$ii]){
						$dataAll2[$ii] = \@data;
					}else{
						@{$dataAll2[$ii]} = @data;
					}
					$ref->[2] = $median;
				}
			}
		}
		close INFO;
		my @blank2;
		if($copy_option eq "all" || $copy_option eq "original"){
			foreach my $x (@headers2){ push(@blank2,$MISSING); }
			for(my $i=0; $i<$irow; $i++){
				if(!$dataAll2[$i]){
					$dataAll2[$i] = \@blank2;
				}
			}
		}
		my $fileID = get_outputID(1);
		open(OUT,">$PATH_OUTPUT/$fileID.txt") or error_message("Cannot open file","continue");
		print OUT $header_text;
		print OUT "Probe_id\t".join("\t",@headers1,@headers2)."\n";
		for(my $i=0; $i<$irow; $i++){
			if($dataAll1[$i] && $dataAll2[$i]){
				print OUT "$probeID[$i]\t".join("\t",@{$dataAll1[$i]},@{$dataAll2[$i]})."\n"
			}
		}
		print OUT "!series_matrix_table_end\n";
		close OUT;
		print "<HTML><HEAD><TITLE>ExAtlas</TITLE>\n";
		print_header();
		print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
		print "<!--\n";
		print "function count_samples(x){\n";
		print "	var nchecked = 0;\n";
		print "	for(i=0; i<document.file_edit.elements.length; ++i){\n";
		print "		var box = document.file_edit.elements[i]\n";
		print "		if(box.type==\"checkbox\"){\n";
		print "			if(x==1 && box.name.search(/^select-1/)>=0 && box.checked){ nchecked++; };\n";
		print "			if(x==2 && box.name.search(/^select-2/)>=0 && box.checked){ nchecked++; };\n";
		print "		}\n";
		print "	}\n";
		print "	return nchecked;\n";
		print "}\n";
		print "function finish_normalize() {\n";
		print "	if(document.file_edit.option_norm.selectedIndex==2 && (count_samples(1)==0 || count_samples(2)==0)){\n";
		print "		alert(\"You need to select at least one samples in each file (#1 and #2)! Cannot proceed.\"); return false;\n";
		print "	}\n";
		print "	if(document.file_edit.option_norm.selectedIndex!=2 && count_samples(1)+count_samples(2)>0){\n";
		print "		if(!confirm(\"You have checked samples that are not needed for selected normalization method.\\nDo you want to proceed and ignore checked samples?\")){\n";
		print "			return false;\n";
		print "		}\n";
		print "	}\n";
		print "	document.file_edit.action.value = \"file_edit2\";\n";
		print "	document.file_edit.submit();\n";
		print "}\n";
		print "function cancel_update() {\n";
		print "	document.file_edit.action.value = \"continue\";\n";
		print "	document.file_edit.submit();\n";
		print "}\n";
		print "// -->\n";
		print "</SCRIPT>\n";
		print "<FORM NAME=file_edit ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST onSubmit=\"return upload_onsubmit();\">\n";
		print "<p><FONT SIZE=5><b>Copying Samples from File #2 to File #1</b></FONT>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;\n";
		print "<INPUT TYPE=button VALUE=\"Cancel & return to main menu\" onClick=cancel_update(); style=width:260px;><p>\n";
		print "<TABLE BORDER=0>\n";
		print "<TR><TD><b>File #1:<TD>$file_copy<TD WIDTH=15><TD><b>Platform:<TD>$platform1<TD WIDTH=15><TD><b>Description:<TD>$file_copy_description\n";
		print "<TR><TD><b>File #2:<TD>$filename<TD><TD><b>Platform:<TD>$platform2<TD><TD><b>Description:<TD>$new_description\n";
		print "</TABLE><p>";
		print "<b>Data Normalization Option:</b>\n";
		print "<SELECT NAME=option_norm style=width:220px;>\n";
		print "<option value=quantile> Quantile\n";
		print "<option value=median> Median\n";
		print "<option value=columns> Selected columns\n";
		print "<option value=none> No normalization</SELECT>\n";
		print "<INPUT TYPE=button VALUE=\"Finish copying samples\" onClick=finish_normalize(); style=width:260px;><br>\n";
		print "<p><b>Notes:</b> (1) Quantile normalization is applied to each sample independently.";
		print "It is not recommended if combined files use different array plaforms<br>\n";
		print "(2) Median method equalizes median values in two data sets for each probe (or gene).<br>\n";
		print "(3) Option \"Selected columns\" equalizes median values in selected sets of columns for each probe (or gene).\n";
		print "If you use this option, then you have to select at least one column in each data set below.<p>\n";
		print "<TABLE BORDER=0><TR><TD><TD><b>Columns in File #1<TD><TD><b>Columns in File #2\n";
		my $ncol = @headers1;
		if($ncol < @headers2){ $ncol = @headers2; }
		for(my $i=0; $i<$ncol; $i++){
			print "<TR>";
			if($i<@headers1){ print "<TD><INPUT NAME=select-1-$i TYPE=CHECKBOX><TD>$headers1[$i]\n"; }
			else{ print "<TD><TD>\n"; }
			if($i<@headers2){ print "<TD WIDTH=20><TD><INPUT NAME=select-2-$i TYPE=CHECKBOX><TD>$headers2[$i]\n"; }
		}
		my $nCol1 = @headers1;
		print "</TABLE><p>";
		print "<INPUT TYPE=button VALUE=\"Finish copying samples\" onClick=finish_normalize(); style=width:260px;><br>\n";
		print "<INPUT TYPE=button VALUE=\"Cancel & return to main menu\" onClick=cancel_update(); style=width:260px;>\n";
		print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
		print "<INPUT NAME=action TYPE=hidden VALUE=file_edit2>\n";
		print "<INPUT NAME=fileID TYPE=hidden VALUE=$fileID>\n";
		print "<INPUT NAME=nCol1 TYPE=hidden VALUE=$nCol1>\n";
		print "<INPUT NAME=organismID TYPE=hidden VALUE=\"$organismID\">\n";
		print "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$file_copy\">\n";
		print "</FORM>\n";
		print "</BODY></HTML>\n";
		exit(0);
	}
	return;
}
if($filename=~/^public-/ && $loginname ne "public"){ error_message("Cannot edit public files"); }
my @item_descriptions;
my @item_delete=();
my $changed=0;
my $platform_old;
my $platform_new;
if($organismID != $organismID_old){
	$changed = 1;
}
if($file_type eq "matrix"){
	$platform_old = $hashFile{"series_platform_id"};
	$platform_new = $hashInput{"platform_new"};
	if($platform_new && $platform_new ne $platform_old){
		$changed = 1;
	}else{ $platform_new=$platform_old; }
}
elsif($file_type eq "geneset"){
	my $ref = $hashFile{"geneset_descriptions"};
	if($ref && ref($ref) eq 'ARRAY'){
		@item_descriptions = @$ref;
	}
}
for(my $i=0; $i<@item_names; $i++){
	my $x = $hashInput{"item_name$i"};
	my $y = $hashInput{"item_description$i"};
	my $z = $hashInput{"item_select$i"};
	if($x && $x ne $item_names[$i]){
		$item_names[$i] = $x;
		$changed = 1;
	}
	if($file_type eq "geneset" && $y && $y ne $item_descriptions[$i]){
		$item_descriptions[$i] = $y;
		$changed = 1;
	}
	if($z eq "on" && $selected_items =~ /^delete/){
		$item_delete[$i] = 1;
		$changed = 1;
	}
}
if($file_type =~ /^geneset|^matrix/ && $changed || $geneset_filtered){
	my $fileID = get_outputID(1);
	my $Ngenesets=0;
	my $count=-1;
	my $deleted=0;
	open(INFO,"<$PATH_DATA/$loginname-$filename") or error_message("Cannot open file $filename","continue");
	open(OUT, ">$PATH_OUTPUT/$fileID.txt") or error_message("Cannot open file $fileID.txt","continue");
	if($file_type eq "matrix"){
		while(my $line = <INFO>){
			if($line =~ /^\!Series_platform_id\t/i){
				print OUT "!Series_platform_id\t\"$platform_new\"\n";
			}elsif($line =~ /^\!Sample_title\t/i){
				if(!@item_delete){
					print OUT "!Sample_title\t\"".join("\"\t\"",@item_names)."\"\n";
				}else{
					my @list = apply_mask(\@item_names,\@item_delete,1);
					print OUT "!Sample_title\t\"".join("\"\t\"",@list)."\"\n";
				}
			}elsif($line =~ /^\!Series_sample_taxid\t/i){
				print OUT "!Series_sample_taxid\t\"$organismID\"\n";
			}elsif($line =~ /^\!Series_platform_taxid\t/i){
				print OUT "!Series_platform_taxid\t\"$organismID\"\n";
			}else{
				if(!@item_delete){
					print OUT $line;
				}else{
					chop $line;
					my ($title,@data) = split(/\t/,$line);
					if(@data < @item_names){
						print OUT "$line\n";
					}else{
						my @list = apply_mask(\@data,\@item_delete,1);
						my $missing = count_missing(\@list);
						if($missing<@list){
							print OUT "$title\t".join("\t",@list)."\n";
						}
					}
				}
			}
		}
		my $file_anova = "$loginname-anova-$filename";
		if(file_exist("$PATH_DATA/$file_anova")){
			unlink "$PATH_DATA/$file_anova";
		}
	}elsif($file_type eq "geneset"){
		my @data_geneset;
		while(my $line = <INFO>){
			chop $line;
			if($line =~ /^\!/ || !$line){
				if($line =~ /^\!Geneset_taxid\t/i){
					print OUT "!Geneset_taxid\t$organismID\n";
				}else{
					print OUT $line."\n";
				}
				next;
			}
			my ($title,$description,@items)=split(/\t/,$line);
			if($title){
				$count++;
				if(@data_geneset){
					$Ngenesets+=filter_geneset(\@data_geneset,\@geneset_attributes,\@geneset_attrib_min,\@geneset_attrib_max);
					@data_geneset=();
				}
			}
			if(!$item_delete[$count]){
				if($geneset_filtered){
					push(@data_geneset,$line);
				}else{
					if($title){ $Ngenesets++; }
					print OUT $line."\n";
				}
			}else{ $deleted++; }
		}
		if(@data_geneset){
			$Ngenesets+=filter_geneset(\@data_geneset,\@geneset_attributes,\@geneset_attrib_min,\@geneset_attrib_max);
			@data_geneset=();
		}
		if($geneset_filtered && !$Ngenesets){
			error_message("No data fits filtering criteria, filtering cancelled","continue");
		}
	}
	close INFO;
	close OUT;
	if($file_type ne "geneset" || $count+1-$deleted>0){
		copy "$PATH_OUTPUT/$fileID.txt", "$PATH_DATA/$loginname-$filename";
		unlink "$PATH_OUTPUT/$fileID.txt";
	}
}
my $text;
my $description;
my $nLines;
my $organismID_old;
my $record;
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!","continue");
while(my $line = <INFO>){
	chop $line;
	if(!$line){ next; }
	$nLines++;
	my %hash;
	read_config_line($line,\%hash);
	if($hash{"type_$file_type"} eq $filename){
		$description = $hash{"description"};
		$organismID_old = $hash{"organismID"};
		$record = $line;
	}else{
		$text .= $line."\n";
	}
}
close INFO;
if($changed || $new_filename && $new_filename ne $filename 
 || $new_description && $new_description ne $description
 || $organismID ne $organismID_old){
	my @items = split(/\t/,$record);
	my $descrFound=0;
	for(my $i=0; $i<@items; $i++){
		if($items[$i] =~ /^type_$file_type/ && $new_filename && $new_filename ne $filename){
			$items[$i] = "type_$file_type=$new_filename";
			system("mv","$PATH_DATA/$loginname-$filename","$PATH_DATA/$loginname-$new_filename");
			my $file_anova = "$loginname-anova-$filename";
			if(file_exist("$PATH_DATA/$file_anova")){
				system("mv","$PATH_DATA/$file_anova","$PATH_DATA/$loginname-anova-$new_filename");
			}
			$filename = $new_filename;
		}
		if($items[$i] =~ /^description/ && $new_description && $new_description ne $description){
			$items[$i] = "description=$new_description";
			$descrFound=1;
		}
		if($items[$i] =~ /^organismID/ && $organismID ne $organismID_old){
			$items[$i] = "organismID=$organismID";
		}
	}
	if(!$descrFound && $new_description){
		if(@items>=3){ splice(@items,2,0,"\tdescription=$new_description"); }
		else{ push(@items,"\tdescription=$new_description"); }
	}
	$record = join("\t",@items);
	if(open(OUT,">$PATH_INFO/$loginname-config1.txt")){
		print OUT $text.$record."\n";
		close OUT;
		my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
		if($nLines1 && $nLines1 > $nLines*0.9){
			copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
		}else{
			error_message("Failed to update configuration file!");
		}
		unlink "$PATH_INFO/$loginname-config1.txt";
	}
}
return;
}

#**************************************
sub file_edit2
#**************************************
{
my $organismID = $hashInput{"organismID"};
my $file_matrix = $hashInput{"file_matrix"};
my $fileID = $hashInput{"fileID"};
my $option_norm = $hashInput{"option_norm"};
my $nCol1 = $hashInput{"nCol1"};

if($file_matrix=~/^public-/ && $loginname ne "public"){ error_message("Cannot edit public files"); }
my $anova_file = "$PATH_DATA/$loginname-anova-$file_matrix";
my $file_matrix_full = "$PATH_DATA/$loginname-$file_matrix";
if(file_exist($anova_file)){
	unlink "$anova_file";
}
if($option_norm eq "none"){
	copy "$PATH_OUTPUT/$fileID.txt", "$PATH_DATA/$file_matrix_full";
	return;
}
if($option_norm eq "quantile"){
	my %hashMatrix;
	parse_matrix_file("$PATH_OUTPUT/$fileID.txt",\%hashMatrix);
	my $ncol = @{$hashMatrix{"sample_title"}};
	my $logFileID = get_outputID(1);
	$hashInput{"logFileID"} = $logFileID;
	$hashInput{"runID"} = $RUN_NORMALIZE;
	interrupt_program($ncol/200);
	return;
}
open(INFO,"<$PATH_OUTPUT/$fileID.txt") or error_message("Cannot open $fileID.txt");
open(OUT,">$PATH_DATA/$file_matrix_full") or error_message("Cannot open $file_matrix_full");
my $line;
my $done=0;
while($line=<INFO>){
	if($line =~ /^!series_normalized\t\"true\"/i){ next; }
	print OUT $line;
	if($line =~ /^!series_matrix_table_begin/i){ last; }
}
$line = <INFO>;
$line =~ s/\"//g;
print OUT $line;
my($junk,@headers) = split(/\t/,$line);
my $nCol = @headers;
my $nCol2 = $nCol-$nCol1;
my @med1;
my @med2;
my @Data1;
my @Data2;
my @probes;
my @select1;
my @select2;
if($option_norm eq "columns"){
	for(my $i=0; $i<$nCol1; $i++){
		if($hashInput{"select-1-$i"} eq "on"){
			$select1[$i] = 1;
		}
	}
	for(my $i=0; $i<$nCol2; $i++){
		if($hashInput{"select-2-$i"} eq "on"){
			$select2[$i] = 1;
		}
	}
}
while($line=<INFO>){
	if($line =~ /^!/){ last; }
	chop $line;
	$line =~ s/\s+$//;
	my($probe,@data1) = split(/\t/,$line);
	$probe =~ s/^\"//;
	$probe =~ s/\"$//;
	my @data2 = splice(@data1,$nCol1);
	my ($median1,$median2);
	if($option_norm eq "median"){
		$median1 = median(\@data1);
		$median2 = median(\@data2);
	}elsif($option_norm eq "columns"){
		my @columns1 = apply_mask(\@data1,\@select1);
		my @columns2 = apply_mask(\@data2,\@select2);
		$median1 = median(\@columns1);
		$median2 = median(\@columns2);
	}else{
		error_message("Unknown option $option_norm");
	}
	if($median1>0 && $median2>0){
		push(@med1,log($median1));
		push(@med2,log($median2));
		for(my $i=0; $i<@data2; $i++){
			my $x = $data2[$i];
			if($x > $MISSING){
				$x *= $median1/$median2;
				$data2[$i] = int(100000*$x+0.5)/100000;
			}
		}
	}
	if($median1>0){
		print OUT "$probe\t".join("\t",@data1,@data2)."\n";
	}else{
		push(@probes,$probe);
		push(@Data1,\@data1);
		push(@Data2,\@data2);
	}
}
close INFO;
my ($r,$n,$bb,$aa) = pearson_correlation(\@med2,\@med1);
for(my $i=0; $i<@probes; $i++){
	my $ref = $Data2[$i];
	for(my $j=0; $j<@$ref; $j++){
		my $x = $ref->[$j];
		if($x > $MISSING){
			my $y = exp(log($x)*$bb+$aa);
			$ref->[$j] = int(100000*$y+0.5)/100000;
		}elsif($x!=$MISSING){
			$ref->[$j] = $MISSING;
		}
	}
	print OUT "$probes[$i]\t".join("\t",@{$Data1[$i]},@$ref)."\n";
}
print OUT $line;
close OUT;
return;
}

#**************************************
sub file_edit3
#**************************************
{
my $logFileID = shift;
my $file_matrix = $hashInput{"file_matrix"};
my $fileID = $hashInput{"fileID"};
my $nCol1 = $hashInput{"nCol1"};

my $fileID1 = get_outputID(1);
my $response = system("$PATH_BIN/norm_new","$PATH_OUTPUT/$fileID.txt","$PATH_OUTPUT/$fileID1.txt","-n","$nCol1");
if($response){ error_message("Normalization failed!",$logFileID); }
copy "$PATH_OUTPUT/$fileID1.txt", "$PATH_DATA/$loginname-$file_matrix";
return;
}

#**************************************
sub  apply_mask
#**************************************
{
my $input = shift;
my $mask = shift;
my $reverse = shift;
if(!$input || !$mask || ref($input) ne 'ARRAY' || ref($mask) ne 'ARRAY'){ return; }
my $n = @$input;
my @output;
for(my $i=0; $i<$n; $i++){
	if(!$reverse && $mask->[$i] || $reverse && !$mask->[$i]){
		push(@output,$input->[$i]);
	}
}
return @output;
}

#**************************************
sub count_missing
#**************************************
{
my $ref = shift;
if(!$ref || ref($ref) ne 'ARRAY'){ return 0; }
my $count = 0;
foreach my $x (@$ref){
	if($x <= $MISSING){ $count++; }
}
return $count;
}

#**************************************
sub filter_geneset
#**************************************
{
my $data_geneset=shift;
my $attrib_names=shift;
my $attrib_min=shift;
my $attrib_max=shift;
if(!$data_geneset || ref($data_geneset) ne 'ARRAY'){ error_message("data_geneset in filter_geneset"); }
if(!$attrib_names || ref($attrib_names) ne 'ARRAY'){ error_message("attrib_names in filter_geneset"); }
if(!$attrib_min || ref($attrib_min) ne 'ARRAY'){ error_message("attrib_min in filter_geneset"); }
if(!$attrib_max || ref($attrib_max) ne 'ARRAY'){ error_message("attrib_max in filter_geneset"); }
my %hash;
for(my $i=0; $i<@$attrib_names; $i++){
	$hash{$attrib_names->[$i]} = $i+1;
}
my @removed = ();
for(my $i=1; $i<@$data_geneset; $i++){
	my ($blank,$name,@values) = split(/\t/,$data_geneset->[$i]);
	my $j = $hash{$name};
	if(!$j){ next; }
	my $min = $attrib_min->[$j-1];
	my $max = $attrib_max->[$j-1];
	for(my $i1=0; $i1<@values; $i1++){
		my $x = $values[$i1];
		if($x<$min || $x>$max){
			$removed[$i1]=1;
		}
	}
}
for(my $i=0; $i<@$data_geneset; $i++){
	my ($name,$descr,@values) = split(/\t/,$data_geneset->[$i]);
	for(my $j=@removed-1; $j>=0; $j--){
		if($removed[$j]){ splice(@values,$j,1); }
	}
	if(@values>0){
		print OUT join("\t",$name,$descr,@values)."\n";
	}else{
		return 0;
	}
}
return 1;
}


#**************************************
sub add_geneset 
#**************************************
{
my $fileID = $hashInput{"fileID"};
my $file_description = $hashInput{"description_copy_file"};
my $filename = $hashInput{"copy_to_geneset"};
my $geneset_name = $hashInput{"geneset_name_new"};
my $geneset_description = $hashInput{"geneset_description_new"};
my $organismID = $hashInput{"organismID"};

if($filename =~ /^public/){ error_message("Cannot add geneset to public-"); }
my $filename_full = "$PATH_DATA/$loginname-$filename";
my %hashFile;
my %hash_symbols;
open (INFO,"<$PATH_OUTPUT/$fileID.txt");
my $line;
my $direction;
while($line=<INFO>){
	chop $line;
	if(!$line){ next; }
	if($line !~/^!/){ last; }
	if($line =~/^!Number_of_over-expressed_genes/){ $direction=1; }
	elsif($line =~/^!Number_of_under-expressed_genes/){ $direction=-1; }
}
my @headers = split(/\t/,$line);
my ($isymb,$ilogr,$ifdr,$ifold)=(-1,-1,-1,-1);
for(my $i=0; $i<@headers; $i++){
	if($headers[$i]=~/symbol/i){ $isymb=$i; }
	elsif($headers[$i]=~/logratio/i){ $ilogr=$i; }
	elsif($headers[$i]=~/fold change/i){ $ifold=$i; }
	elsif($headers[$i]=~/FDR/i){ $ifdr=$i; }
}
if($isymb<0){ error_message("Column with gene symbols not found"); }
my @attrib_names=();
if($ilogr>=0 || $ifold>=0 && $direction){ push(@attrib_names,"logratio"); }
if($ifdr>=0){ push(@attrib_names,"FDR"); }
my @geneList;
my $log10 = log(10.0);
while(my $line=<INFO>){
	chop $line;
	if(!$line || $line=~/^!/){ next; }
	my @items = split(/\t/,$line);
	my $symbol = $items[$isymb];
	if($symbol && !$hash_symbols{$symbol}){
		my @attrib;
		if($ilogr>=0){ push(@attrib,$items[$ilogr]); }
		elsif($ifold>=0 && $direction){
			if($items[$ifold]<=1.0E-12){ push(@attrib,0); }
			else{ push(@attrib,$direction*int(10000*log($items[$ifold])/$log10)/10000); }
		}
		if($ifdr>=0){ push(@attrib,$items[$ifdr]); }
		push(@geneList,$symbol);
		if(@attrib_names){ $hash_symbols{$symbol}=\@attrib; }
		else{ $hash_symbols{$symbol}=1; }
	}
}
close INFO;
my %hash_attributes;
my ($use_logratio,$use_fdr)=(1,1);
if(file_exist($filename_full)){
	parse_geneset_file($filename_full, \%hashFile);
	my $ref = $hashFile{"geneset_names"};
	if($ref && ref($ref) eq 'ARRAY'){
		foreach my $name (@$ref){
			if($name eq $geneset_name){ error_message("Geneset named '$geneset_name' already exists in $filename"); }
		}
	}
	$ref = $hashFile{"geneset_attributes"};
	if($ref && ref($ref) eq 'ARRAY'){
		foreach my $x (@$ref){ $hash_attributes{$x}=1; }
	}
	if(!$hash_attributes{"logratio"}){ $use_logratio=0; }
	if(!$hash_attributes{"FDR"}){ $use_fdr=0; }
}else{
	#Update configuration file
	open(OUT, ">$PATH_INFO/$loginname-config1.txt");
	open(INFO,"<$PATH_INFO/$loginname-config.txt");
	my $nLines;
	while(my $line = <INFO>){
		$line =~ s/\n$//;
		if(!$line){ next; }
		$nLines++;
		if($line !~ /^type_geneset=$filename\s/){
			print OUT $line."\n";
		}
	}
	close INFO;
	print OUT "type_geneset=$filename\torganismID=$organismID";
	if($file_description){ print OUT "\tdescription=$file_description"; }
	print OUT "\tdate=$date_record\n";
	close OUT;
	my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
	if($nLines1 && $nLines1 > $nLines*0.9){
		copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
	}else{
		error_message("Failed to update configuration file!");
	}
	unlink("$PATH_INFO/$loginname-config1.txt");
}
open (OUT,">>$filename_full") or error_message("Cannot open filename_full");
print OUT "$geneset_name\t$geneset_description\t".join("\t",@geneList)."\n";
for(my $i=0; $i<@attrib_names; $i++){
	my $attr=$attrib_names[$i];
	if($attr eq "logratio" && !$use_logratio || $attr eq "FDR" && !$use_fdr){ next; }
	my @attr_list;
	foreach my $symbol (@geneList){
		my $x=0;
		my $ref = $hash_symbols{$symbol};
		if(ref($ref) eq 'ARRAY'){ $x=$ref->[$i]; }
		push(@attr_list,$x);
	}
	print OUT "\t$attr\t".join("\t",@attr_list)."\n";
}
close OUT;
terminal_window("<H3>Geneset '$geneset_name' is saved in file '$filename'</H3>");
exit(0);
}

#**************************************
sub anova_parameters
#**************************************
{
my $organismID = $hashInput{"organismID"};
my $file_matrix = $hashInput{"file_matrix"};
my $description_matrix = $hashInput{"description_matrix"};
my $file_matrix_full = "$loginname-$file_matrix";
my %hashMatrix;
my $file_platform = get_array_platform($file_matrix_full,\%hashMatrix);
if($file_platform=~/ /){ error_message("Wrong file platform $file_platform"); }
print "<HTML><HEAD><TITLE>ExAtlas</TITLE>\n";
print_header();
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "function run_anova() {\n";
print "	if(document.anova.prop_var.value > 0.5 || document.anova.prop_var.value < 0){\n";
print "		alert(\"Proportion of ignored error variance should be from 0 to 0.5. Cannot proceed.\"); return false;\n";
print "	}\n";
print "	if(document.anova.window_width.value > 1000 || document.anova.window_width.value < 3){\n";
print "		alert(\"Window width for averaging error variance should be from 3 to 1000. Cannot proceed.\"); return false;\n";
print "	}\n";
print "	var x = document.anova.error_model.value;\n";
print "	if(x!=1 && x!=2 && x!=3 && x!=4 && x!=5){\n";
print "		alert(\"Error model should be integer from 1 to 5. Cannot proceed.\"); return false;\n";
print "	}\n";
print "	var x = document.anova.df_bayesian.value;\n";
print "	if(x<4 || x>1000){\n";
print "		alert(\"Bayesian degrees of freedom should be integer from 4 to 1000. Cannot proceed.\"); return false;\n";
print "	}\n";
print "	document.anova.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print "<FORM NAME=anova ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST onSubmit=\"return upload_onsubmit();\">\n";
print "<p><FONT SIZE=5><b>Run ANOVA again</b></FONT>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;\n";
print "<INPUT TYPE=button VALUE=\" Cancel (close window) \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "<b>File name:</b> $file_matrix<br>\n";
my $description1 = add_hyperlinks($description_matrix);
print "<b>Description:</b> $description1<br>\n";
print "<b>Platform:</b> $file_platform<p>\n";
print "<INPUT NAME=cutoff VALUE=0 style=width:160px;> Cutoff (probes with maximum value below cutoff are ignored)<br>\n";
print "<INPUT NAME=z_outliers VALUE=8 style=width:160px;> Threshold z-value used to remove outliers (0 - 1000)<br>\n";
print "<INPUT NAME=prop_var VALUE=0.01 style=width:160px;> Proportion of probes with high error variances to ignore in error models (0 - 0.5)<br>\n";
print "<INPUT NAME=window_width VALUE=500 style=width:160px;> Number of probes in a sliding window to average error variance (4 - 1000)<br>\n";
print "<INPUT NAME=error_model VALUE=4 style=width:160px;> Error model number (see the list below)<br>\n";
print "<INPUT NAME=df_bayesian VALUE=10 style=width:160px;> Bayesian degrees of freedom (used with error models 3 or 5) (4 - 1000)<br>\n";
print "<select name=use_probeID style=width:160px;>\n";
print "<option value=0> Never\n";
print "<option value=1> If gene symbol is missing\n";
print "<option value=2> Always\n";
print "</select> use probe ID as gene symbol<br>\n";
print "<p><b>List of error models</b> <br>\n";
print "1 = Actual error variance for each probe<br>\n";
print "2 = Average error variance for probes with similar expression level<br>\n";
print "3 = Bayesian correction of error variance (Baldi & Long 2001)<br>\n";
print "4 = Maximum between actual and expected average error variances<br>\n";
print "5 = Maximum between actual and Bayesian error variances<p>\n";
print "<INPUT TYPE=button VALUE=\"Run ANOVA\" onClick=run_anova(); style=width:260px;><br>\n";
print "<INPUT TYPE=button VALUE=\"Cancel (close window)\" LANGUAGE=\"javascript\" onClick=\"window.close();\" style=width:260px;><p>\n";
print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=action TYPE=hidden VALUE=matrix_explore>\n";
print "<INPUT NAME=analysis TYPE=hidden VALUE=run_anova>\n";
print "<INPUT NAME=file_matrix TYPE=hidden VALUE=\"$file_matrix\">\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=\"$organismID\">\n";
print "</FORM>\n";
print "</BODY></HTML>\n";
exit(0);
}

#**************************************
sub  get_array_platform
#**************************************
{
my $file_matrix = shift;
my $hashMatrix = shift;

my $organismID = $hashInput{"organismID"};
if(!$hashMatrix){
	my %hashMatrix=();
	$hashMatrix = \%hashMatrix;
}
parse_matrix_file("$PATH_DATA/$file_matrix", $hashMatrix);
my $platform = $hashMatrix->{"series_platform_id"};
if($platform == $organismID){ $platform = "symbol_$organismID"; }
if($platform =~ /^gpl\d+$/){ $platform = uc($platform); }
if(!$platform){
	return "Platform not specified";
}
my $organismID1 = $hashMatrix->{"series_sample_taxid"};
if(!$organismID1){ return "No organism ID in $file_matrix"; }
if($organismID && $organismID1 ne $organismID){ return "Organism ID does not match"; }
if($platform eq "None"){
	return $platform;
}
my $file_platform = "$loginname-$platform"."_annot.txt";
if(file_exist("$PATH_DATA/$file_platform")){
	return $file_platform;
}
$file_platform = $platform."_annot.txt";
if(file_exist("$PATH_DATA/$file_platform")){
	return $file_platform;
}
$file_platform = "public-$platform"."_annot.txt";
if(file_exist("$PATH_DATA/$file_platform")){
	return $file_platform;
}
# Get array platform from GEO
if($platform =~ /^GPL\d+$/){
	$file_platform = "public-$platform"."_annot.txt";
	my $error_message = "";
	my $warning_message = "";
	my $return = download_platform_GEO($platform,$organismID,\$error_message,\$warning_message);
	if($return==0){
		return "$error_message $warning_message";
	}
	return $file_platform;
}
return "No such platform $platform";
}

#******************************
sub  parse_matrix_file
#******************************
{
my $file_matrix = shift;
my $hashSeries = shift;
my $comment = shift;

my @contributors;
#print "<p>Sub = parse_matrix_file\n";
open(INFO,'<',$file_matrix) or error_message("Cannot open $file_matrix",$comment);
while(my $line=<INFO>){
	chop $line;
	if(!$line){ next; }
	$line =~ s/[^[:ascii:]]//g;
	if($line =~ /^!series_matrix_table_begin/i){ last; }
	$line =~ s/^!//;
	$line =~ s/\"//g;
	$line =~ s/\s+$//g;
	my @items = split(/\t/,$line);
	if($items[0] =~ /^series/i){
		if($items[1] =~ /^keywords: /i){
			$items[1] =~ s/^keywords: //i;
			$hashSeries->{"keywords"} = $items[1];
		}elsif($items[0] =~ /^series_contributor/i){
			my($first,$scnd,$name) = split(/,/,$items[1]);
			$first = substr($first,0,1);
			push(@contributors,$name." $first");
		}else{
			my $key = lc($items[0]);
			my $value = $hashSeries->{$key};
			if($value){
				$hashSeries->{$key} = "$value,$items[1]";
			}else{
				$hashSeries->{$key} = $items[1];
			}
		}
	}elsif($items[0] =~ /^sample/i){
		my $key = shift(@items);
		if($items[0] =~ /\: / && $key !~ /sample_title/i){
			my($key1,$value) = split(/\: /,$items[0]);
			my $nn=0;
			while($key1 =~ / /g){ $nn++; }
			if(length($key1) <= 30 && $nn < 4 && $key1 !~ /[\(\)\{\}\@\#\$\^\&\*\%\!\?\+\.]/){
				$key .= ",$key1";
				$key =~ s/ /_/g;
				$key1 .= ": ";
				for(my $i=0; $i<@items; ++$i){
					$items[$i] =~ s/$key1//;
				}
			}
		}
		$hashSeries->{lc($key)} = \@items;
	}
}
if($hashSeries->{"series_platform_id"} =~ /,/ && $hashSeries->{"sample_platform_id"}){
	$hashSeries->{"series_platform_id"} = $hashSeries->{"sample_platform_id"}->[0];
}
if($hashSeries->{"series_platform_taxid"} =~ /,/ && $hashSeries->{"sample_platform_taxid"}){
	$hashSeries->{"series_platform_taxid"} = $hashSeries->{"sample_platform_taxid"}->[0];
}
if($hashSeries->{"series_sample_taxid"} =~ /,/ && $hashSeries->{"series_platform_taxid"} !~ /,/){
	$hashSeries->{"series_sample_taxid"} = $hashSeries->{"series_platform_taxid"};
}
$hashSeries->{"series_contributor"} = \@contributors;
my @DATA;
my $line=<INFO>;
my $count_rows=0;
my $Ncol=0;
while(my $line=<INFO>){
	if($line =~ /^!/){ last; }
	$line =~ s/[^[:ascii:]]//g;
	$count_rows++;
	if(@DATA < 200000){
		chop $line;
		my($id,@data) = split(/\t/,$line);
		if(!$Ncol){ $Ncol=@data; }
		foreach my $x (@data){
			if($x > 0.0000001){ push(@DATA,$x); }
		}
	}
}
@DATA = sort {$a<=>$b} @DATA;
my $xmin  = $DATA[int(0.005*@DATA)];
my $xmax  = $DATA[int(0.995*@DATA)];
$hashSeries->{"series_xmin"} = $xmin;
$hashSeries->{"series_xmax"} = $xmax;
my @row_counts;
for(my $i=0; $i<$Ncol; $i++){ push(@row_counts,$count_rows); }
$hashSeries->{"sample_data_row_count"} = \@row_counts;
#print "Xmax = $xmax<br>\n";
#print "Xmin = $xmin<br>\n";
close INFO;
return;
}

#**************************************
sub  parse_file_headers
#**************************************
{
my $filename = shift;
my $hash = shift;
my $comment = shift;

my $filename_short = $filename;
$filename_short =~ s/.+\///g;
open(INFO,'<',$filename) or error_message("Cannot open file '$filename_short' in parse_file_headers",$comment);
my $line;
while($line = <INFO>){
	if($line !~ /^\!/){ last; }
	chop $line;
	$line =~ s/[^[:ascii:]]//g;
	my ($key,@items)=split(/\t/,$line);
	$key =~ s/^\!//;
	$key = lc($key);
	if(@items==1){
		$hash->{$key} = $items[0];
	}elsif(@items>1){
		$hash->{$key} = \@items;
	}
}
chop $line;
my ($junk,@headers)=split(/\t/,$line);
if(@headers){
	$hash->{"column_titles"} = \@headers;
}
my @rows;
while(my $line = <INFO>){
	if($line =~ /^\!/){ last; }
	my ($row,$value)=split(/\t/,$line);
	push(@rows,$row);
}
$hash->{"row_titles"} = \@rows;
close INFO;
if($line !~ /^\!/){ return 0; }
return 1;
}

#**************************************
sub  parse_geneset_file
#**************************************
{
my $filename = shift;
my $hash = shift;
my $comment = shift;

open(INFO,'<',$filename) or error_message("Cannot open file '$filename' in parse_geneset_file",$comment);
my @names;
my @descriptions;
my %hashAttrib;
my @attribNames;
while(my $line = <INFO>){
	chop $line;
	if(!$line){ next; }
	if($line =~ /^\!/){
		my ($key,@items)=split(/\t/,$line);
		if(!@items){ next; }
		$key =~ s/^\!//;
		$key = lc($key);
		if(@items==1){
			$hash->{$key} = $items[0];
		}elsif(@items>1){
			$hash->{$key} = \@items;
		}
		next;
	}
	my ($name,$descr,@items)=split(/\t/,$line);
	if($name){
		push(@names,$name);
		push(@descriptions,$descr);
	}elsif($descr){
		my ($min,$max) = get_bounds(\@items);
		if($min == $MISSING){ next; }
		my $ref = $hashAttrib{$descr};
		if(!$ref){
			$hashAttrib{$descr}=[$min,$max];
			push(@attribNames,$descr);
		}else{
			my($min1,$max1) = @$ref;
			if($min<$min1){ $ref->[0] = $min; }
			if($max>$max1){ $ref->[1] = $max; }
		}
	}
}
close INFO;
my @attribMin;
my @attribMax;
foreach my $attr (@attribNames){
	my $ref = $hashAttrib{$attr};
	push(@attribMin,$ref->[0]);
	push(@attribMax,$ref->[1]);
}
$hash->{"geneset_names"} = \@names;
$hash->{"geneset_descriptions"} = \@descriptions;
if(@attribNames){
	$hash->{"geneset_attributes"} = \@attribNames;
	$hash->{"geneset_attribute_min"} = \@attribMin;
	$hash->{"geneset_attribute_max"} = \@attribMax;
}
return;
}


#***********************************
sub  get_gene_homolog
#***********************************
{
my $organismID1 = shift; #input
my $organismID = shift;  #output
my $geneHomolog = shift;
my $logfileID = shift;

if(!$geneHomolog){ error_message("Hash pointer homologene missing",$logfileID); }
if(!$organismID || !$organismID1){ error_message("In get_gene_homolog: organismID missing",$logfileID); }
%$geneHomolog=();
my @homologene;
open (INFO_HOM,"<$PATH_DATA/homologene.data") or error_message("Homologene file not opened",$logfileID);
while(my $line = <INFO_HOM>){
	chop $line;
	if(!$line){ next; }
	my($homologeneID,$speciesID,$geneID,$symbol,$proteinID,$protein) = split(/\t/,$line);
	if($speciesID==$organismID1 || $speciesID==$organismID){
		my $ii=0;
		if($speciesID==$organismID){ $ii=1; }
		$homologene[$homologeneID]->[$ii] = $symbol;
	}
}
close INFO_HOM;
foreach my $ref (@homologene){
	if($ref->[0] && $ref->[1]){
		$geneHomolog->{$ref->[0]} = $ref->[1];
	}
}
return;
}

#***********************************
sub   meta_analysis_combine
#***********************************
{
my $fileID = $hashInput{"fileID"};
my $file_matrix = $hashInput{"file_matrix"};
my $column1 = $hashInput{"select_column"};
my $column2 = $hashInput{"compare_column"};
my $organismID = $hashInput{"organismID"};
my $organismID1 = $hashInput{"organismID1"};
if(!$organismID1){ $organismID1=$organismID; }
my $command = $hashInput{"command"};
my $file_metaanalysis = $hashInput{"file_metaanalysis"};
my $description_metaanalysis = $hashInput{"description_metaanalysis"};
my $select_metaanalysis = $hashInput{"select_metaanalysis"};
my @metaanalysis_list = (["--- New file ---",""]);
my @file_list;
my %hashMatrixOrganism;
my @matrix_list;
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!");
while(my $line = <INFO>){
	chop $line;
	my %hash=();
	read_config_line($line,\%hash);
	my @items = split(/[=\t]/,$line);
	if($line=~/^type_metaanalysis/){
		my $file = $hash{"type_metaanalysis"};
		my $description = $hash{"description"};
		push(@metaanalysis_list,[$file,$description]);
	}elsif($line =~ /^type_/){
		if($line=~/^type_matrix/){
			my $file_matrix1 = $hash{"type_matrix"};
			my $org = $hash{"organismID"};
			$org =~ s/;.+$//;
			push(@matrix_list,[$file_matrix1,$hash{"description"},$org]);
			$hashMatrixOrganism{$file_matrix1}=$org;
		}
		push(@file_list,$items[1]);
	}
}
close INFO;
my $comment;
if($command eq "start_analysis"){
	meta_compute();
}elsif($command =~ /^print_table/){
	meta_analysis_table();
}elsif($command =~ /^save_geneset/){
	my $file_geneset = $hashInput{"file_geneset_new"};
	my $description_geneset = $hashInput{"description_geneset_new"};
	my $outputFileID = $hashInput{"outputFileID"};
	my $FDR_thresh = $hashInput{"FDR"};
	my $rowTable = $command;
	$rowTable =~ s/^save_geneset,//;
	my @geneset;
	open(INFO,"<$PATH_OUTPUT/$outputFileID.txt") or error_message("Cannot open table");
	my $line = <INFO>;
	while($line = <INFO>){
		chop $line;
		my ($id,$symbol,$title,$sumEff,$fold_combined,@data) = split(/\t/,$line);
		if($rowTable==4 || $data[$rowTable*2+1] > $FDR_thresh){ next; }
		if($sumEff>0){ $geneset[0]->{$symbol}=1; }
		else{ $geneset[1]->{$symbol}=1; }
	}
	close INFO;
	
	open (OUT, ">$PATH_DATA/$loginname-$file_geneset") or error_message("Cannot write to geneset");
	for(my $j=0; $j<2; $j++){
		my $refHash = $geneset[$j];
		if(!$refHash){ next; }
		my @symbols = sort keys %$refHash;
		my $setName = "metaanalysis_$method_name[$rowTable]$dirNames[$j]";
		$setName =~ s/ /_/g;
		print OUT "$setName\t$setName\t".join("\t",@symbols)."\n";
	}
	close OUT;
	open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!");
	open(OUT, ">$PATH_INFO/$loginname-config1.txt") or error_message("Cannot update configuration file!");
	my $nLines;
	while(my $line = <INFO>){
		chop $line;
		if(length($line)<3){ next; }
		$nLines++;
		my @items = split(/[=\t]/,$line);
		if($items[0] ne "type_geneset" || $items[1] ne $file_geneset){
			print OUT $line."\n";
		}
	}
	close INFO;
	print OUT "type_geneset=$file_geneset\torganismID=$organismID";
	if($description_geneset){ print OUT "\tdescription=$description_geneset"; }
	print OUT "\tdate=$date_record\n";
	close OUT;
	my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
	if($nLines1 && $nLines1 > $nLines*0.9){
		copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
	}else{
		error_message("Failed to update configuration file!");
	}
	unlink "$PATH_INFO/$loginname-config1.txt";
	$hashInput{"file_geneset"} = $file_geneset;
	$hashInput{"description_geneset"} = $description_geneset;
	geneset_explore();
	exit(0);
}elsif($command =~ /^save_metaanalysis/){
	$file_metaanalysis = $hashInput{"select_metaanalysis"};
	if($file_metaanalysis =~ /^New_file/){ $file_metaanalysis=$hashInput{"file_metaanalysis_new"}; }
	if(!$file_metaanalysis){ error_message("Metaanalysis file is blank"); }
	$description_metaanalysis = $hashInput{"description_metaanalysis_new"};
	copy "$PATH_OUTPUT/$fileID.txt", "$PATH_DATA/$loginname-$file_metaanalysis";
	open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!");
	open(OUT, ">$PATH_INFO/$loginname-config1.txt") or error_message("Cannot update configuration file!");
	my $nLines;
	while(my $line = <INFO>){
		chop $line;
		if(length($line)<3){ next; }
		$nLines++;
		my @items = split(/[=\t]/,$line);
		if($items[0] ne "type_metaanalysis" || $items[1] ne $file_metaanalysis){
			print OUT $line."\n";
		}
	}
	close INFO;
	print OUT "type_metaanalysis=$file_metaanalysis\torganismID=$organismID";
	if($description_metaanalysis){ print OUT "\tdescription=$description_metaanalysis"; }
	print OUT "\tdate=$date_record\n";
	close OUT;
	my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
	if($nLines1 && $nLines1 > $nLines*0.9){
		copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
	}else{
		error_message("Failed to update configuration file!");
	}
	unlink "$PATH_INFO/$loginname-config1.txt";
	push(@metaanalysis_list,[$file_metaanalysis,$description_metaanalysis]);
	$comment = "Meta-analysis \'$file_metaanalysis\' is saved";
}elsif($command =~ /^delete_metaanalysis/){
	if(!$select_metaanalysis){ error_message("Metaanalysis file is blank"); }
	my $found=0;
	open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!");
	open(OUT, ">$PATH_INFO/$loginname-config1.txt") or error_message("Cannot update configuration file!");
	my $nLines;
	while(my $line = <INFO>){
		chop $line;
		if(length($line)<3){ next; }
		$nLines++;
		my @items = split(/[=\t]/,$line);
		if($items[0] eq "type_metaanalysis" && $items[1] eq $select_metaanalysis){
			$found = 1;
		}else{
			print OUT $line."\n";
		}
	}
	close INFO;
	close OUT;
	$comment .= "Meta-analysis \'$select_metaanalysis\' is deleted";
	if(!$found){ $comment .= ". File not found"; }
	else{
		my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
		if($nLines1 && $nLines1 > $nLines*0.9){
			copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
		}else{
			error_message("Failed to update configuration file!");
		}
	}
	unlink "$PATH_INFO/$loginname-config1.txt";
	if(file_exist("$PATH_DATA/$loginname-$select_metaanalysis")){
		unlink "$PATH_DATA/$loginname-$select_metaanalysis";
	}
	for(my $i=0; $i<@metaanalysis_list; $i++){
		my $ref = $metaanalysis_list[$i];
		if($ref->[0] eq $select_metaanalysis){
			splice(@metaanalysis_list,$i,1);
			last;
		}
	}
	if($select_metaanalysis eq $file_metaanalysis){
		terminal_window("<H3>$comment<H3>");
	}
}elsif($command =~ /^load_metaanalysis/){
	foreach my $ref (@metaanalysis_list){
		if($ref->[0] eq $select_metaanalysis){
			$description_metaanalysis = $ref->[1];
			last;
		}
	}
	$file_metaanalysis = $select_metaanalysis;
	copy "$PATH_DATA/$loginname-$file_metaanalysis", "$PATH_OUTPUT/$fileID.txt";
	$comment = "Meta-analysis \'$file_metaanalysis\' is loaded";
}
@matrix_list = sort {lc($a->[0]) cmp lc($b->[0])} @matrix_list;
if($loginname ne "public" && open(INFO,"<$PATH_INFO/public-config.txt")){
	my @matrix_list1;
	while(my $line = <INFO>){
		chop $line;
		my %hash=();
		if(length($line) < 3) { next; }
		read_config_line($line,\%hash);
		my $file_matrix1 = $hash{"type_matrix"};
		if($file_matrix1){
			my $org = $hash{"organismID"};
			push(@matrix_list1,["public-$file_matrix1",$hash{"description"},$org]);
			$hashMatrixOrganism{"public-$file_matrix1"}=$org;
		}
	}
	push(@matrix_list, sort {lc($a->[0]) cmp lc($b->[0])} @matrix_list1);
}
close INFO;
my %hashOrganismID;
foreach my $ref (@matrix_list){
	$hashOrganismID{$ref->[2]}=1;
}
filter_list_by_organism(\@matrix_list, $organismID1);
my @anova_headers;
my @header_col;
my @meta_analysis;
my $warning_message;
my %hashData;
if($fileID && open(INFO,"<$PATH_OUTPUT/$fileID.txt")){
	while(my $line=<INFO>){
		chop $line;
		$hashData{$line}=1;
		my @items = split(/\t/,$line);
		push(@meta_analysis,\@items);
	}
	close INFO;
	if($command eq "delete_data"){
		my $fileID1 = get_outputID(1);
		open(OUT,">$PATH_OUTPUT/$fileID1.txt");
		my $i1=0;
		for(my $i=0; $i<@meta_analysis; ++$i){
			if($hashInput{"useSample_$i1"} eq "on"){
				splice(@meta_analysis,$i,1);
				$i--; $i1++;
				next;
			}
			print OUT join("\t",@{$meta_analysis[$i]})."\n";
			$i1++;
		}
		close OUT;
		if(!@meta_analysis){ error_message("All data deleted"); }
		copy "$PATH_OUTPUT/$fileID1.txt", "$PATH_OUTPUT/$fileID.txt";
	}
}
for(my $i=0; $i<@matrix_list; ++$i){
	my $file = $matrix_list[$i]->[0];
	my $file_anova1 = $file;
	$file_anova1 =~ s/^public-/public-anova-/;
	if($file_anova1 !~ /^public-/){
		$file_anova1 = "$loginname-anova-$file";
	}
	my @headers = get_anova_headers($file_anova1,1);
	if(!@headers){
		splice(@matrix_list,$i,1);
		$i--;
		next;
	}
	my $headers_list="\"Median profile\"";
	for(my $j=0; $j<@headers; ++$j){
		$headers_list .= ",\"".$headers[$j]."\""; 
	}
	push(@anova_headers,$headers_list);
	#print "A1 $file eq $file_matrix<br>\n";
	if($file eq $file_matrix){
		if($command eq "add_data" || !$fileID){
			my $header1 = $headers[$column1-1];
			my $header2 = "Median profile";
			if($column2>0){ $header2 = $headers[$column2-1]; }
			my $file_matrix_full = $file;
			if($file !~ /^public-/){
				$file_matrix_full = "$loginname-$file";
			}
			my %hashMatrix=();
			parse_matrix_file("$PATH_DATA/$file_matrix_full",\%hashMatrix);
			my $platform = $hashMatrix{"series_platform_id"};
			my @items = ($file_matrix,$platform,$organismID1,$column1,$column2,$header1,$header2);
			my $data_line = join("\t",@items);
			if($hashData{$data_line}){
				$warning_message = "Duplicated entry! Command ignored";
			}elsif(!$fileID){
				$fileID = get_outputID(1);
				file_append("$data_line","$PATH_OUTPUT/$fileID.txt",1);
				push(@meta_analysis,\@items);
			}else{
				file_append("$data_line","$PATH_OUTPUT/$fileID.txt");
				push(@meta_analysis,\@items);
			}
		}
	}
}
my $matrix_metaanalysis = $meta_analysis[0]->[0];
my $organism_metaanalysis = $hashMatrixOrganism{$matrix_metaanalysis};
if(!$fileID){ error_message("No fileID in meta-analysis"); }
my $file_list="";
foreach my $name (@file_list){
	if(!$file_list){ $file_list = "\"".$name."\""; }
	else{ $file_list .= ",\"".$name."\""; }
}
my ($items,$descriptions) = get_array_lists(\@matrix_list);
my ($ma_items,$ma_descriptions) = get_array_lists(\@metaanalysis_list);
my $n_pairs = @meta_analysis;

# Print page header
print "<HTML><HEAD><TITLE>ExAtlas - correlation</TITLE>\n";
print_header("update_description();");
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "var file_list = new Array($file_list);\n";
print "var matrix_description = new Array($descriptions);\n";
print "var ma_list = new Array($ma_items);\n";
print "var ma_description = new Array($ma_descriptions);\n";
print "var anova_headers = new Array();\n";
for(my $i=0; $i<@anova_headers; ++$i){
	print "anova_headers[$i] = new Array($anova_headers[$i]);\n";
}
print "function update_description() {\n";
print "	var index = document.metaanalysis.file_matrix.selectedIndex;\n";
print "	document.metaanalysis.target = \"\";\n";
print "	document.metaanalysis.select_column.options.length=0;\n";
print "	for(i=1; i<anova_headers[index].length; ++i){\n";
print "		var x = new Option(anova_headers[index][i],i);\n";
print "		document.metaanalysis.select_column.options[i-1] = x;\n";
print "	}\n";
print "	document.metaanalysis.compare_column.options.length=0;\n";
print "	for(i=0; i<anova_headers[index].length; ++i){\n";
print "		var x = new Option(anova_headers[index][i],i);\n";
print "		document.metaanalysis.compare_column.options[i] = x;\n";
print "	}\n";
print "	document.metaanalysis.select_column.selectedIndex=0;\n";
print "	document.metaanalysis.compare_column.selectedIndex=0;\n";
print "	document.metaanalysis.description_matrix.value = matrix_description[index];\n";
print "	index = document.metaanalysis.select_metaanalysis.selectedIndex;\n";
print "	document.metaanalysis.description_metaanalysis_new.value = ma_description[index];\n";
print "}\n";
print "function add_data(){\n";
print "	document.metaanalysis.target = \"\";\n";
print "	if(document.metaanalysis.select_column.selectedIndex==document.metaanalysis.compare_column.selectedIndex-1){\n";
print "		alert(\"You cannot compare expression data to itself\");\n";
print "		return(false);\n";
print "	}\n";
print "	document.metaanalysis.action.value =\"meta-analysis\";\n";
print "	document.metaanalysis.command.value =\"add_data\";\n";
print "	document.metaanalysis.submit();\n";
print "}\n";
print "function change_organism() {\n";
print "	document.metaanalysis.target = \"\";\n";
print "	document.metaanalysis.action.value =\"meta-analysis\";\n";
print "	document.metaanalysis.command.value =\"change_organism\";\n";
print "	document.metaanalysis.submit();\n";
print "}\n";
print "function count_samples(){\n";
print "	var nchecked = 0;\n";
print "	for(i=0; i<document.metaanalysis.elements.length; ++i){\n";
print "		var box = document.metaanalysis.elements[i]\n";
print "		if(box.type==\"checkbox\" && box.checked){ nchecked++; };\n";
print "	}\n";
print "	return nchecked;\n";
print "}\n";
print "function delete_data() {\n";
print "	if(count_samples()==0){\n";
print "		alert(\"No samples selected. Cannot proceed.\"); return false;\n";
print "	}\n";
print "	if(!confirm(\"Do you want to delete selected data?\")){\n";
print "		return false;\n";
print "	}\n";
print "	document.metaanalysis.target = \"\";\n";
print "	document.metaanalysis.action.value =\"meta-analysis\";\n";
print "	document.metaanalysis.command.value =\"delete_data\";\n";
print "	document.metaanalysis.submit();\n";
print "}\n";
print "function start_analysis() {\n";
print "	if($n_pairs<2){\n";
print "		alert(\"You need at least 2 pairs of data for meta-analysis\");\n";
print "		return(false);\n";
print "	}\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.metaanalysis.target = \"_BLANK\"+x;\n";
print "	document.metaanalysis.action.value =\"meta-analysis\";\n";
print "	document.metaanalysis.command.value =\"start_analysis\";\n";
print "	document.metaanalysis.submit();\n";
print "}\n";
print "function save_metaanalysis() {\n";
print " var index1 = document.metaanalysis.select_metaanalysis.selectedIndex;\n";
print " var file = document.metaanalysis.select_metaanalysis.options[index1].value;\n";
print "	var descrip = document.metaanalysis.description_metaanalysis_new.value;\n";
print " if(index1==0){\n";
print "		file = prompt(\"Provide name of a new file where to copy metaanalysis\");\n";
print "		if(!file){ return(false); }\n";
print "		var file1=file;\n";
print "		if(file.search(/\\.txt\$/)>=0){\n";
print "			if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "			file1=file.substring(0,file.length-4);\n";
print "		}\n";
print "		if(file1.search(/^[-\\w]+\$/)<0){\n";
print "			alert(\"File name should have neither spaces nor special characters\");\n";
print "			return(false);\n";
print "		}\n";
print "		if(file.search(/^public/i) >= 0){\n";
print "			alert(\"File name cannot start with 'public'\");\n";
print "			return(false);\n";
print "		}\n";
print "		for(i=0; i<file_list.length; ++i){\n";
print "			if(file == file_list[i]){\n";
print "				alert(\"File with this name already exists\"); return false;\n";
print "			}\n";
print "		}\n";
print "		if(descrip.search(/\\=|\\&/) >= 0){\n";
print "			alert(\"Description should not include character \'=\' or \'&\'\");\n";
print "			return false;\n";
print "		}\n";
print "		document.metaanalysis.file_metaanalysis_new.value = file;\n";
print "	}else if(!confirm(\"Do you want to overwrite/update existing metaanalysis file?\")){\n";
print "		return false;\n";
print "	}\n";
print "	document.metaanalysis.target = \"\";\n";
print "	document.metaanalysis.action.value =\"meta-analysis\";\n";
print "	document.metaanalysis.command.value =\"save_metaanalysis\";\n";
print "	document.metaanalysis.submit();\n";
print "}\n";
print "function load_metaanalysis() {\n";
print " if(document.metaanalysis.select_metaanalysis.selectedIndex==0){\n";
print "		alert(\"Use pull-down menu to select metaanalysis file that you saved before\");\n";
print "		return(false);\n";
print "	}\n";
if(@metaanalysis_list){
	print "	if(!confirm(\"Operation will overwrite your current data\\nDo you want to proceed?\")){\n";
	print "		return false;\n";
	print "	}\n";
}
print "	document.metaanalysis.target = \"\";\n";
print "	document.metaanalysis.action.value =\"meta-analysis\";\n";
print "	document.metaanalysis.command.value =\"load_metaanalysis\";\n";
print "	document.metaanalysis.submit();\n";
print "}\n";
print "function delete_metaanalysis() {\n";
print " if(document.metaanalysis.select_metaanalysis.selectedIndex==0){\n";
print "		alert(\"No file selected. Nothing to delete\");\n";
print "		return(false);\n";
print "	}\n";
	print "	if(!confirm(\"Do you want to delete selected file?\")){\n";
	print "		return false;\n";
	print "	}\n";
print "	document.metaanalysis.target = \"\";\n";
print "	document.metaanalysis.action.value =\"meta-analysis\";\n";
print "	document.metaanalysis.command.value =\"delete_metaanalysis\";\n";
print "	document.metaanalysis.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
my @FDR_list = (1,0.5,0.2,0.1,0.05,0.01,0.001,0.0001);
my @fold_change = (1,1.1,1.2,1.3,1.5,2,3,4,5,10);
print "<TABLE BORDER=0><TR><TD WIDTH=558>\n";
if($file_metaanalysis){
	print "<p style=font-size:24px><b>Meta-analysis file: $file_metaanalysis</b></p>\n";
	if($description_metaanalysis){ print "<b>Description:</b> $description_metaanalysis\n"; }
}
if($comment){ print "<br><font color=red><b>Comment:</b></font> $comment\n"; }
if($warning_message){ print "<br><font color=red><b>Warning:</b></font> $warning_message\n"; }
print "<p style=font-size:24px>Input data for meta-analysis</p>\n";
print "<TD VALIGN=TOP><INPUT TYPE=button VALUE=\"Cancel (close window)\" style=\"width: 200px;\" LANGUAGE=javascript onClick=window.close();>\n";
print "</TABLE>\n";
print "<FORM NAME=metaanalysis ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<TABLE BORDER=0><TR><TD WIDTH=558><p style=font-size:20px><b>1. Pairs of expression profiles</b>\n";
print "<TD><INPUT TYPE=button VALUE=\"Delete checked data\" style=\"width: 200px;\" onClick=delete_data();>\n";
print "</TABLE>\n";
print "<TABLE BORDER=0>\n";
for(my $i=0; $i<@meta_analysis; $i++){
	my $num = $i+1;
	my ($file_matrix1,$platform,$taxid,$col1,$col2,$header1,$header2) = @{$meta_analysis[$i]};
	print "<TR><TD WIDTH=55><center>($num) <INPUT NAME=useSample_$i TYPE=checkbox><TD COLSPAN=5><b>$file_matrix1</b> platform: $platform; organism: $hashOrganism{$taxid}\n";
	print "<TR><TD><TD WIDTH=120><i>Sample name:<TD WIDTH=218><font color=MediumBlue>$header1\n";
	print "<TD WIDTH=30><TD WIDTH=120><i>Compare with:<TD><font color=MediumBlue>$header2\n";
}
print "</TABLE><p>\n";
print "<TABLE BORDER=0><TR><TD WIDTH=558><p style=font-size:20px><b>2. Parameters of meta-analysis</b>\n";
print "<TD><INPUT TYPE=button VALUE=\"Start analysis\" style=\"width: 200px;\" onClick=start_analysis();>\n";
print "</TABLE>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD WIDTH=100><a href=../exatlas-help.html#fdr>FDR</a> threshold\n";
print "<TD WIDTH=330><select name=FDR style=\"width: 200px;\">\n";
for(my $i=0; $i<@FDR_list; ++$i){ 
	print "<option value=$FDR_list[$i]"; if($FDR_list[$i]==0.05){ print " selected"; } print "> $FDR_list[$i]\n";
}
print "</select>\n";
print "<TD WIDTH=124>Fold change:\n";
print "<TD><select name=fold_change style=\"width: 200px;\">\n";
for(my $i=0; $i<@fold_change; ++$i){ 
	print "<option value=$fold_change[$i]"; if($fold_change[$i]==2){ print " selected"; } print "> $fold_change[$i]\n";
}
print "</select>\n";
print "</TABLE>\n";
print "<p style=font-size:20px><b>3. Save or load meta-analysis</b></p>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD WIDTH=100>Select file:\n";
print "<TD WIDTH=330><SELECT NAME=select_metaanalysis onChange=update_description(); style=width:200px;>\n";
for(my $i=0; $i<@metaanalysis_list; ++$i){
	my $select;
	if($metaanalysis_list[$i]->[0] eq $file_metaanalysis){ $select=" selected"; }
	my $file = $metaanalysis_list[$i]->[0];
	if($file =~ /- new file -/i){ print "<option value=\"New_file\"$select> $file\n";  }
	else{ print "<option value=\"$file\"$select> $file\n"; }
}
print "</SELECT><TD WIDTH=124>Description:";
print "<TD><INPUT NAME=description_metaanalysis_new style=width:200px;>\n";
print "<TR><TD><TD COLSPAN=2><INPUT TYPE=button VALUE=\"Save metaanalysis\" onClick=save_metaanalysis(); style=width:200px;>";
print "<INPUT TYPE=button VALUE=\"Load metaanalysis\" onClick=load_metaanalysis(); style=width:200px;>\n";
print "<TD><INPUT TYPE=button VALUE=\"Delete metaanalysis\" onClick=delete_metaanalysis(); style=width:200px;>";
print "</TABLE><br>\n";
print "If you select \"New file\" you will be prompted for file name, which shound be one-word without special characters (underscore allowed)<p>\n";

print "<TABLE BORDER=0><TR><TD WIDTH=562><p style=font-size:20px><b>4. Add a pair of expression profiles</b>\n";
print "<TD><INPUT TYPE=button VALUE=\"Add data\" style=\"width: 200px;\" onClick=add_data();>\n";
print "</TABLE>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD WIDTH=100><b>Organism:</b><TD WIDTH=330>$hashOrganism{$organismID1}\n";
if(keys %hashOrganismID > 1){
	print "<TD>Select organism:\n";
	print "<TD><select name=organismID1 style=width:200px; onChange=change_organism();>\n";
	foreach my $ref (@organisms){
		if(!$hashOrganismID{$ref->[0]}){ next; }
		print "<option value=$ref->[0]";
		if($ref->[0] == $organismID1){ print " selected"; }
		print ">$ref->[3] ($ref->[2])\n";
	}
	print "</select> (page is reloaded after change)\n";
}else{
	print "<INPUT NAME=organismID1 TYPE=hidden VALUE=$organismID>\n";
}
print "<TR><TD>Select file:<TD><select name=file_matrix style=width:200px; onChange=update_description();>\n";
for(my $i=0; $i<@matrix_list; ++$i){
	print "<option value=\"$matrix_list[$i]->[0]\"";
	if($matrix_list[$i]->[0] eq $file_matrix){ print " selected"; }
	print "> $matrix_list[$i]->[0]\n";
}
print "</select><TD WIDTH=124>Description:<TD><INPUT NAME=description_matrix style=width:200px;>\n";
print "<TR><TD>Select data:<TD><select name=select_column style=width:200px;></select>\n";
print "<TD>Compare with:<TD><select name=compare_column style=width:200px;></select>\n";
print "</TABLE><p>\n";
print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=action TYPE=hidden VALUE=meta-analysis>\n";
print "<INPUT NAME=command TYPE=hidden>\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=\"$organismID\">\n";
print "<INPUT NAME=fileID TYPE=hidden VALUE=\"$fileID\">\n";
print "<INPUT NAME=file_metaanalysis_new TYPE=hidden>\n";
print "<INPUT NAME=organism_metaanalysis TYPE=hidden VALUE=$organism_metaanalysis>\n";
print "<INPUT NAME=matrix_metaanalysis TYPE=hidden VALUE=$matrix_metaanalysis>\n";
if($file_metaanalysis){ print "<INPUT NAME=file_metaanalysis TYPE=hidden VALUE=$file_metaanalysis>\n"; }
if($description_metaanalysis){ print "<INPUT NAME=description_metaanalysis TYPE=hidden VALUE=\"$description_metaanalysis\">\n"; }
print "</FORM><p>\n";
print "<HR NOSHADE></HR>\n";
print "<p><i>Please report any problems to <a href=mailto:sharoval\@mail.nih.gov>webmaster</a><p>\n";
print "<INPUT TYPE=button VALUE=\"Cancel (close window)\" style=\"width: 200px;\" LANGUAGE=javascript onClick=window.close();><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#***********************************
sub   meta_compute
#***********************************
{
my $fileID = $hashInput{"fileID"};
my $organismID = $hashInput{"organism_metaanalysis"};
my $file_matrix = $hashInput{"matrix_metaanalysis"};
my $FDR_thresh = $hashInput{"FDR"};
my $fold_thresh = $hashInput{"fold_change"};

my $count=0;
my $same_platform =-1;
open(INFO,"<$PATH_OUTPUT/$fileID.txt") or error_message("FileID not open");
while(my $line=<INFO>){
	chop $line;
	my ($file_matrix,$platform,@items) = split(/\t/,$line);
	$count++;
	if($same_platform==-1){ $same_platform = $platform; }
	elsif($same_platform ne $platform){ $same_platform =0; }
}
my $logFileID = get_outputID(3);
$hashInput{"logFileID"} = $logFileID;
$hashInput{"runID"} = $RUN_METAANALYSIS;
interrupt_program($count/10);
exit(0);
}

#***********************************
sub   meta_compute1
#***********************************
{
my $logFileID = shift;
my $webPageID = $hashInput{"logFileID"}+1;
my $outputFileID = $webPageID+1;

my $fileID = $hashInput{"fileID"};
my $organismID = $hashInput{"organism_metaanalysis"};
my $file_matrix = $hashInput{"matrix_metaanalysis"};
my $FDR_thresh = $hashInput{"FDR"};
my $fold_thresh = $hashInput{"fold_change"};
my $log10 = log(10);
my $log_thresh = log($fold_thresh)/$log10;

my @annotation;
my %lookup_row;
my @DATA;
my @meta_analysis;
my @Nreplications;
my $same_platform =-1;

open(INFO,"<$PATH_OUTPUT/$fileID.txt") or error_message("FileID not open",$logFileID);
while(my $line=<INFO>){
	chop $line;
	my ($file_matrix,$platform,@items) = split(/\t/,$line);
	push(@meta_analysis,[$file_matrix,$platform,@items]);
	if($same_platform==-1){ $same_platform = $platform; }
	elsif($same_platform ne $platform){ $same_platform =0; }
}
my %hashGenes;
my %hashAlias;
get_official_symbols($organismID,\%hashGenes,\%hashAlias);
my $metaDelta = 20;
if(@meta_analysis < 10){ $metaDelta = 1; }
elsif(@meta_analysis < 100){ $metaDelta = 5; }
for(my $imatrix=0; $imatrix<@meta_analysis; $imatrix++){
	if($logFileID && $imatrix%$metaDelta==0){
		file_append("Data pair #$imatrix","$PATH_OUTPUT/$logFileID.txt");
	}
	my($file_matrix,$platform,$organismID1,$column1,$column2,$header1,$header2) = @{$meta_analysis[$imatrix]};
	my $file_matrix_full = $file_matrix;
	my $file_anova = $file_matrix;
	$file_anova =~ s/^public-/public-anova-/;
	if($file_matrix !~ /^public-/){
		$file_matrix_full = "$loginname-$file_matrix";
		$file_anova = "$loginname-anova-$file_matrix";
	}
	my %geneHomolog=();
	my %hashGenes1;
	my %hashAlias1;
	if($organismID1 != $organismID){
		get_official_symbols($organismID1,\%hashGenes1,\%hashAlias1);
		get_gene_homolog($organismID1,$organismID,\%geneHomolog);
	}
	open(INFO,"<$PATH_DATA/$file_anova") or error_message("Cannot open ANOVA file",$logFileID);
	my $line = <INFO>;
	chop $line;
	my ($junk,$junk1,@headers) =split(/\t/,$line);
	if(!@headers){ error_message("Headers not found in anova file",$logFileID); }
	my @n_repl;
	my $nCol=0;
	my $nRep=0;
	while($nCol<@headers && $headers[$nCol] !~ /^Var\(/){ ++$nCol; }
	for(my $i=0; $i<$nCol; ++$i){
		$n_repl[$i] = 1;
		my @items = split(/ \(/,$headers[$i]);
		my $n = pop(@items);
		my $name = join(" (",@items);
		$n =~ s/\)$//;
		$nRep += $n;
		$n_repl[$i] = $n;
		if($n_repl[$i] < 1){ $n_repl[$i]=1; }
		$headers[$i] = $name;
	}
	my $icol = $column1-1;
	my $coefficient = 1.0/$n_repl[$icol];
	$Nreplications[$imatrix] = $n_repl[$icol];
	if($column2>0){ $coefficient += 1.0/$n_repl[$column2-1]; }
	else{ $coefficient += 1.0/$nCol; }
	my @quality;
	my $irow=0;
	while(my $line = <INFO>){
		chop $line;
		my ($probe,$expr,@data1) =split(/\t/, $line);
		my $F = $data1[$nCol+4];
		my $MSE = $data1[$nCol+3];
		my $FDR = $data1[$nCol+6];
		my $symbol = $data1[$nCol+8];
		my $gene_name = $data1[$nCol+9];
		if($organismID1==$organismID){
			if(!$hashGenes{$symbol}){
				my $symbol1 = $hashAlias{$symbol};
				if($symbol1){ $symbol=$symbol1; }
				else{ next; }
			}
		}else{
			if(!$hashGenes1{$symbol}){
				my $symbol1 = $hashAlias1{$symbol};
				if($symbol1){ $symbol=$symbol1; }
			}
			my $symbol1 = $geneHomolog{$symbol};
			if($symbol1){
				$symbol = $symbol1;
			}elsif($hashAlias{$symbol}){
				$symbol = $hashAlias{$symbol};
			}elsif(!$hashGenes{$symbol}){
				next;
			}
		}
		#if($nRep<=$nCol){ $F=$expr; }
		my $y = $data1[$icol];
		my $x;
		if($column2==0){
			splice(@data1,$nCol);
			$x = median(\@data1);
		}
		else{ $x = $data1[$column2-1]; }
		if($x==$MISSING || $y==$MISSING){ next; }
		$x = floor(10000*$x+0.5)/10000;
		$y = floor(10000*$y+0.5)/10000;
		my $z = 0;
		my $var = $MSE*$coefficient;
		if($var<=0){ $var=$EPS; }
		if($MSE*$coefficient > 0){
			$z = ($y-$x)/sqrt($var);
			$z = int(10000*$z+0.5)/10000;
		}
		my $id = $symbol;
		if($same_platform){ $id = $probe; }
		my $j = $lookup_row{$id};
		if(!defined($j)){
			$j = $irow;
			push(@annotation,[$id,$symbol,$gene_name]);
			$lookup_row{$id} = $irow++;
		}
		if(!$quality[$j]){
			$DATA[$j]->[$imatrix] = [$z,$y-$x,$var];
			$quality[$j] = $F;
		}elsif($F > $quality[$j]){
			@{$DATA[$j]->[$imatrix]} = ($z,$y-$x,$var);
			$quality[$j] = $F;
		}
	}
	close INFO;
}
my @results;
my @w;
for(my $j=0; $j<@meta_analysis; $j++){
	$w[$j] = sqrt($Nreplications[$j]);
}
for(my $i1=0; $i1<@DATA; $i1++){
	my ($chisq,$n,$sumz,$sumw,$sumww,$sumEff,$sumEff1,$sumwwSq)=(0,0,0,0,0,0,0,0,0,0,0);
	for(my $j=0; $j<@meta_analysis; $j++){
		my $ref = $DATA[$i1]->[$j];
		if(!$ref || ref($ref) ne 'ARRAY'){ next; }
		my ($z,$effect,$var) = @$ref;
		my $x;
		if(abs($z) > 7){
			$x = exp(log(abs($z))*1.83)*0.76;
		}else{
			$x = -log(2*(1-normal_distribution(abs($z))));
		}
		$chisq += $x;
		$n++;
		$sumz += $z*$w[$j];
		$sumw += $w[$j]*$w[$j];
		my $ww = 1/$var;
		$sumEff += $effect*$ww;
		$sumww += $ww;
		$sumwwSq += $ww*$ww;
	}
	$sumEff /= $sumww;
	my $Zfixed = $sumEff/sqrt(1.0/$sumww);
	if(abs($sumEff) < $log_thresh){ next; }
	$chisq *= 2;
	my $Zglobal = 0;
	if($sumw>0){ $Zglobal = $sumz/sqrt($sumw); }
	my $p1 = 2*(1-normal_distribution(abs($Zglobal)));	#Z-value method
	my $df = 2*$n;
	my $p2 = gammq(0.5*$df, 0.5*$chisq);		#Fisher's method
	my $p3 = 2*(1-normal_distribution(abs($Zfixed)));
	my $Q = 0;
	my $df1 = -1;
	my $C = $sumww - $sumwwSq/$sumww;
	my @effects;
	for(my $j=0; $j<@meta_analysis; $j++){
		my $ref = $DATA[$i1]->[$j];
		if(!$ref || ref($ref) ne 'ARRAY'){
			push(@effects,$MISSING);
			next;
		}
		my ($z,$effect,$var) = @$ref;
		$Q += ($effect-$sumEff)*($effect-$sumEff)/$var;
		push(@effects,$effect);
		$df1++;
	}
	my ($Zrand,$p4)=(0,1);
	if($df1 >0){
		my $Tsq = ($Q-$df1)/$C;
		if($Tsq<0){ $Tsq=0; }
		$sumww = 0;
		for(my $j=0; $j<@meta_analysis; $j++){
			my $ref = $DATA[$i1]->[$j];
			if(!$ref || ref($ref) ne 'ARRAY'){ next; }
			my ($z,$effect,$var) = @$ref;
			my $ww = 1/($Tsq+$var);
			$sumEff1 += $effect*$ww;
			$sumww += $ww;
		}
		$sumEff1 /= $sumww;
		$Zrand = $sumEff1/sqrt(1.0/$sumww);
		$p4 = 2*(1-normal_distribution(abs($Zrand)));
	}
	push(@results,[$i1,$n,$sumEff,$sumEff1,$p4,$p3,$p2,$p1,1,1,1,1,@effects]);
}
my $nresults=@results;
my $nGenes = @DATA;
for(my $imethod=0; $imethod<4; $imethod++){
	my $k = 4+$imethod;
	my @sorted = sort {$results[$a]->[$k]<=>$results[$b]->[$k]} 0..($nresults-1);
	my $FDR1 = 1;
	for(my $i=$nresults-1; $i>=0; $i--){
		my $j = $sorted[$i];
		my $p =  $results[$j]->[$k];
		my $rank = $i+1;
		my $FDR = $p*$nGenes/$rank;
		if($FDR > $FDR1){ $FDR = $FDR1; }
		else{ $FDR1 = $FDR; }
		$results[$j]->[$k+4] = $FDR;
	}
}
if(!$outputFileID || !open(OUT,">$PATH_OUTPUT/$outputFileID.txt")){
	error_message("cannot open output file1",$logFileID);
}
if($logFileID){
	file_append("Writing output file...","$PATH_OUTPUT/$logFileID.txt");
}
print OUT "Probe ID\tGene symbol\tGene name\tLogratio combined\tFold change combined\tp (Random)\tFDR (Random)\tp (Fixed)\tFDR (Fixed)\tp (Fisher)\tFDR (Fisher)\tp (Zvalue)\tFDR (Zvalue)";
for(my $j1=0; $j1<@meta_analysis; $j1++){
	my $num =$j1+1;
	print OUT "\tlogratio (Data#$num)\tfold change (Data#$num)";
}
print OUT "\n";
my $text;
my @Nsignif;
my @sorted = sort {abs($results[$a]->[2])<=>abs($results[$b]->[2])} 0..($nresults-1);
for(my $i=$nresults-1; $i>=0; $i--){
	my $j = $sorted[$i];
	my ($irow,$n,$sumEff,$sumEff1,$p4,$p3,$p2,$p1,$FDR4,$FDR3,$FDR2,$FDR1,@effects) = @{$results[$j]};
	my ($id,$symbol,$title) = @{$annotation[$irow]};
	my $fold_combined = exp(abs($sumEff)*$log10);
	if($fold_combined > 1){ $fold_combined = floor(1000*$fold_combined+0.5)/1000; }
	else{ $fold_combined = format_probability($fold_combined); }
	$sumEff = floor(10000*$sumEff+0.5)/10000;
	$p1 = format_probability($p1);
	$p2 = format_probability($p2);
	$p3 = format_probability($p3);
	$p4 = format_probability($p4);
	$FDR1 = format_probability($FDR1);
	$FDR2 = format_probability($FDR2);
	$FDR3 = format_probability($FDR3);
	$FDR4 = format_probability($FDR4);
	my $line = "$id\t$symbol\t$title\t$sumEff\t$fold_combined\t$p4\t$FDR4\t$p3\t$FDR3\t$p2\t$FDR2\t$p1\t$FDR1";
	for(my $j1=0; $j1<@effects; $j1++){
		my $x = $effects[$j1];
		if($x <= $MISSING){
			$line .= "\tN/A\tN/A";
		}else{
			my $fold = exp($x*$log10);
			if($fold > 1){ $fold = floor(1000*$fold+0.5)/1000; }
			else{ $fold = format_probability($fold); }
			$x = floor(10000*$x+0.5)/10000;
			$line .= "\t$x\t$fold";
		}
	}
	$line .= "\n";
	my $signif=0;
	if($sumEff > 0){
		if($FDR4 <= $FDR_thresh){ $Nsignif[0]++; $signif++; }
		if($FDR3 <= $FDR_thresh){ $Nsignif[1]++; $signif++; }
		if($FDR2 <= $FDR_thresh){ $Nsignif[2]++; $signif++; }
		if($FDR1 <= $FDR_thresh){ $Nsignif[3]++; $signif++; }
		if($signif){
			print OUT $line;
			$Nsignif[4]++;
		}
	}else{
		if($FDR4 <= $FDR_thresh){ $Nsignif[5]++; $signif++; }
		if($FDR3 <= $FDR_thresh){ $Nsignif[6]++; $signif++; }
		if($FDR2 <= $FDR_thresh){ $Nsignif[7]++; $signif++; }
		if($FDR1 <= $FDR_thresh){ $Nsignif[8]++; $signif++; }
		if($signif){
			$text .= $line;
			$Nsignif[9]++;
		}
	}
}
print OUT $text;
close OUT;

# Print page header
if(!$webPageID || !open(OUT,">$PATH_OUTPUT/$webPageID.txt")){
	error_message("cannot open web page ID",$logFileID);
}
print OUT "<HTML><HEAD><TITLE>ExAtlas - meta-analysis</TITLE>\n";
print OUT get_header();
print OUT "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print OUT "<!--\n";
print OUT "function show_table(i,j){\n";
print OUT "	document.metaresults.command.value =\"print_table,\"+i+\",\"+j;\n";
print OUT "	var x = Math.round(Math.random()*10000);\n";
print OUT "	document.metaresults.target = \"_BLANK\"+x;\n";
print OUT "	document.metaresults.submit();\n";
print OUT "}\n";
print OUT "// -->\n";
print OUT "</SCRIPT>\n";
print OUT "<TABLE BORDER=0>\n";
print OUT "<TR><TD WIDTH=450><p style=font-size:27px><b>Results of meta-analysis (summary)\n";
print OUT "<TD><INPUT TYPE=button VALUE=\"Close window\" style=\"width: 200px;\" LANGUAGE=javascript onClick=window.close();>\n";
print OUT "</TABLE>\n";
print OUT "<FORM NAME=metaresults ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print OUT "<p style=font-size:20px><b>1. Number of significant genes</b> (FDR <= $FDR_thresh and change >= $fold_thresh fold)<p>\n";
print OUT "Click on the number of genes in the table below to see the list of genes.<p>\n";
print OUT "Select effect format: <SELECT name=effect_type style=width:170px;>\n";
print OUT "<option value=0> Log-ratio (log10)<option value=1> Fold change\n";
print OUT "</select><p>\n";
print OUT "<TABLE BORDER=1>\n";
print OUT "<TR><TD WIDTH=160><b>Method of analysis<TD WIDTH=200 ALIGN=CENTER><b>Up-regulated genes<TD WIDTH=200 ALIGN=CENTER><b>Down-regulated genes";
my $Nmethod = @method_name;
for(my $i=0; $i<@method_name; $i++){
	print OUT "<TR><TD>$method_name[$i]";
	for(my $j=0; $j<2; $j++){
		print OUT "<TD ALIGN=center><a href=\"JavaScript:show_table($j,$i);\">$Nsignif[$j*$Nmethod+$i]</a>\n";
	}
}
print OUT "</TABLE><p style=font-size:20px><b>List of data used for meta-analysis<p>\n";
print OUT "<TABLE BORDER=0>\n";
for(my $i=0; $i<@meta_analysis; $i++){
	my $num = $i+1;
	my ($file_matrix1,$platform,$taxid,$col1,$col2,$header1,$header2) = @{$meta_analysis[$i]};
	print OUT "<TR><TD WIDTH=35><center>($num)<TD COLSPAN=5><b>$file_matrix1</b> platform: $platform; organism: $hashOrganism{$taxid}\n";
	print OUT "<TR><TD><TD WIDTH=120><i>Sample name:<TD><font color=MediumBlue>$header1\n";
	print OUT "<TD WIDTH=30><TD WIDTH=120><i>Compare with:<TD><font color=MediumBlue>$header2\n";
}
print OUT "</TABLE><p>\n";
print OUT "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
print OUT "<INPUT NAME=action TYPE=hidden VALUE=meta-analysis>\n";
print OUT "<INPUT NAME=command TYPE=hidden>\n";
print OUT "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
print OUT "<INPUT NAME=file_matrix TYPE=hidden VALUE=$file_matrix>\n";
print OUT "<INPUT NAME=FDR TYPE=hidden VALUE=$FDR_thresh>\n";
print OUT "<INPUT NAME=fileID TYPE=hidden VALUE=$fileID>\n";
print OUT "<INPUT NAME=outputFileID TYPE=hidden VALUE=$outputFileID>\n";
print OUT "</FORM><p>\n";
print OUT "<HR NOSHADE></HR><p>\n";
print OUT "<INPUT TYPE=button VALUE=\"Cancel (close window)\" style=\"width: 200px;\" LANGUAGE=javascript onClick=window.close();><p>\n";
print OUT "</BODY>\n";
print OUT "</HTML>\n";
close OUT;
return;
}

#***********************************
sub   meta_analysis_table
#***********************************
{
my $organismID = $hashInput{"organismID"};
my $outputFileID = $hashInput{"outputFileID"};
my $fileID = $hashInput{"fileID"};
my ($junk,$colTable,$rowTable) = split(/,/,$hashInput{"command"});
my $FDR_thresh = $hashInput{"FDR"};
my $effect_type = $hashInput{"effect_type"};
my $file_matrix = $hashInput{"file_matrix"};

my @tableText;
open(INFO,"<$PATH_OUTPUT/$outputFileID.txt") or error_message("Cannot open table");
my $line = <INFO>;
my ($id,$symbol,$title,$sumEff,$fold_combined,@data) = split(/\t/,$line);
$line = join("\t",$id,$symbol,$title,$sumEff,$fold_combined);
if($rowTable<4){
	$line .= "\tp-value\tFDR";
}else{
	for(my $i=1; $i<8; $i+=2){	
		$line .= "\t$data[$i]";
	}
}
my $N_data = int((@data-8)/2);
for(my $j=1; $j<=$N_data; $j++){
	$line .= "\tEffect (data#$j)";
}
push(@tableText,$line);
my $nGenes=0;
while($line = <INFO>){
	chop $line;
	my ($id,$symbol,$title,$sumEff,$fold_combined,@data) = split(/\t/,$line);
	#Filter by the direction of effect:
	if($colTable==0 && $sumEff<0 || $colTable==1 && $sumEff>0){ next; }
	my @effects = splice(@data,8);
	for(my $j=@effects-1; $j>=0; $j-=2){
		splice(@effects,$j-$effect_type,1);
	}
	if($rowTable<4){	#Specific method
		my @statistics = splice(@data,$rowTable*2,2);
		if($statistics[1] > $FDR_thresh){ next; }
		my $line1 = join("\t",$id,$symbol,$title,$sumEff,$fold_combined,@statistics,@effects);
		push(@tableText,$line1);
		$nGenes++;
	}else{			#Show all
		for(my $j=@data-2; $j>=0; $j-=2){
			splice(@data,$j,1);
		}
		my $line1 = join("\t",$id,$symbol,$title,$sumEff,$fold_combined,@data,@effects);
		push(@tableText,$line1);
		$nGenes++;
	}
}
close INFO;
my $fileTable = get_outputID(1);
open(OUT, ">$PATH_OUTPUT/$fileTable.txt") or error_message("Cannot open output");
foreach my $line (@tableText){ print OUT $line."\n"; }
close OUT;
my @geneset_list;
my $file_list="";
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Config file not found!");
while(my $line = <INFO>){
	chop $line;
	if(length($line) < 3) { next; }
	if($line =~ /^type_geneset=/){
		my %hash=();
		read_config_line($line,\%hash);
		push(@geneset_list,[$hash{"type_geneset"},$hash{"description"},$hash{"organismID"}]);
	}
	my @items = split(/[=\t]/,$line);
	if($items[0] =~ /^type_/){
		if(!$file_list){ $file_list = "\"".$items[1]."\""; }
		else{ $file_list .= ",\"".$items[1]."\""; }
	}
}
close INFO;
@geneset_list = sort {lc($a->[0]) cmp lc($b->[0])} @geneset_list;
if($loginname ne "public" && open(INFO,"<$PATH_INFO/public-config.txt")){
	my @geneset_list1=();
	while(my $line = <INFO>){
		chop $line;
		my %hash=();
		read_config_line($line,\%hash);
		my $file_geneset = $hash{"type_geneset"};
		if($file_geneset){
			push(@geneset_list1,["public-".$file_geneset,$hash{"description"},$hash{"organismID"}]);
		}
	}
	close INFO;
	push(@geneset_list, sort {lc($a->[0]) cmp lc($b->[0])} @geneset_list1);
}
filter_list_by_organism(\@geneset_list, $organismID);
my ($items,$descriptions) = get_array_lists(\@geneset_list);
my @direction_name=("over","under");

# Print page header
print "<HTML><HEAD><TITLE>ExAtlas: Meta-analysis - list of $direction_name[$colTable]-expressed genes</TITLE>\n";
if(@geneset_list && $nGenes>=5){
	print_header("update_description();");
}else{
	print_header();
}
print "<SCRIPT language=JavaScript>\n";
print "<!--\n";
print "geneset_list = new Array($items);\n";
print "geneset_description = new Array($descriptions);\n";
print "file_list = new Array($file_list);\n";
if(@geneset_list && $nGenes>=5){
	print "function geneset_overlap(file){\n";
	print "	document.form_metaanalysis.upload_geneset.value=file;\n";
	print "	document.form_metaanalysis.action.value=\"geneset_overlap\";\n";
	print "	var x = Math.round(Math.random()*10000);\n";
	print "	document.form_metaanalysis.target = \"_BLANK\"+x;\n";
	print "	document.form_metaanalysis.submit();\n";
	print "}\n";
	print "function update_description() {\n";
	print "	var index;\n";
	print "	index = document.form_metaanalysis.file_geneset1.selectedIndex;\n";
	print "	document.form_metaanalysis.description_geneset1.value = geneset_description[index];\n";
	print "}\n";
}
print "function save_genesets () {\n";
print "	var file = document.form_metaanalysis.file_geneset_new.value;\n";
print "	if(!file){\n";
print "		alert(\"Enter file name\"); return(false);\n";
print "	}\n";
print "	var file1=file;\n";
print "	if(file.search(/\\.txt\$/)>=0){\n";
print "		if(file==\".txt\"){ alert(\"File name '.txt' is not valid\"); return(false); }\n";
print "		file1=file.substring(0,file.length-4);\n";
print "	}\n";
print "	if(file1.search(/^[-\\w]+\$/)<0){\n";
print "		alert(\"File name should have neither spaces nor special characters\");\n";
print "		return(false);\n";
print "	}\n";
print "	if(file.search(/^public-/i) >= 0){\n";
print "		alert(\"File name cannot start with 'public-'\"); return(false);\n";
print "	}\n";
print "	for(i=0; i<file_list.length; ++i){\n";
print "		if(file == file_list[i]){\n";
print "			alert(\"A file with this name already exists\"); return false;\n";
print "		}\n";
print "	}\n";
print "	var descrip = document.form_metaanalysis.description_geneset_new.value;\n";
print "	if(descrip.search(/\\=|\\&/) >= 0){\n";
print "		alert(\"Description should not include character \'=\' or \'&\'\");\n";
print "		return false;\n";
print "	}\n";
print "	document.form_metaanalysis.action.value=\"meta-analysis\";\n";
print "	document.form_metaanalysis.command.value=\"save_geneset\"+\",$rowTable\";\n";
print "	document.form_metaanalysis.target = \"\";\n";
print "	document.form_metaanalysis.submit();\n";
print "}\n";
print "<!-- end script --></SCRIPT></HEAD>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD WIDTH=450><p style=font-size:27px><b>Table of $direction_name[$colTable]-expressed genes\n";
print "<TD><INPUT TYPE=button VALUE=\"Close window\" style=\"width: 170px;\" LANGUAGE=javascript onClick=window.close();>\n";
print "</TABLE><p>\n";
print "<b>Number of genes:</b> $nGenes<br>\n";
print "<b>Method of analysis:</b> $method_name[$rowTable]<br>";
print "<b>Note:</b> Effects in the table below are measured as log10(fold change).<br>";
my $x = int(10000*rand());
print "<b>Table of genes</b> as tab-delimited text: <a href=$HOME_ADDRESS/output/$fileTable.txt target=_BLANK$x>Table</a><p>\n";

print "<FORM NAME=form_metaanalysis ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
print "<b>Save lists of genes (both up- and down-regulated) as a geneset file:</b><br>\n";
print "File name: <INPUT NAME=file_geneset_new style=\"width:200px;\"> &nbsp; &nbsp; &nbsp; \n";
print "Description: <INPUT NAME=description_geneset_new style=width:200px;>\n";
print "<INPUT TYPE=button VALUE=\"Save genes\" onClick=save_genesets(); style=\"width:120px;\"><p>\n";
if(@geneset_list && $nGenes>=5){
	print  "<TABLE BORDER=0>\n";
	my $menu_text = menu_geneset_overlap(\@geneset_list,$fileTable);
	print $menu_text;
	print "</TABLE><p>\n";
}
print "<TABLE BORDER=0><TR>\n";
$line = shift(@tableText);
my @items = split(/\t/,$line);
for(my $i=0; $i<@items; $i++){
	print "<TD";
	if($i>2){ print " WIDTH=355"; }
	print "><b>$items[$i]";
}
print "\n";
foreach my $line (@tableText){
	print "<TR>";
	my @items = split(/\t/,$line);
	for(my $i=0; $i<@items; $i++){
		print "<TD>";
		if($i>2){ print "<p style=font-size:12px>"; }
		if($i==1){
			my $x = int(10000*rand());
			print "<a href=\"exatlas.cgi?category=1&search_term=$items[1]&action=matrix_explore1&analysis=search&file_matrix=$file_matrix&sessionID=$sessionID&organismID=$organismID\" target=_BLANK$x>$items[1]<a>";
		}else{
			print "$items[$i]";
		}
	}
	print "\n";
}
print "</TABLE>\n";
print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
print "<INPUT NAME=action TYPE=hidden>\n";
print "<INPUT NAME=command TYPE=hidden>\n";
print "<INPUT NAME=outputFileID TYPE=hidden VALUE=$outputFileID>\n";
print "<INPUT NAME=FDR TYPE=hidden VALUE=$FDR_thresh>\n";
print "<INPUT NAME=organismID TYPE=hidden VALUE=\"$organismID\">\n";
print "<INPUT NAME=upload_geneset TYPE=hidden>\n";
print "</FORM><p>\n";
print "<HR NOSHADE></HR><p>\n";
print "<INPUT TYPE=button VALUE=\"Close window\" style=\"width: 170px;\" LANGUAGE=javascript onClick=window.close();><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#***********************************
sub plot_histogram_horiz
#***********************************
{
my $plot_type = shift;
my $n = shift;
my $m = shift;
my $data = shift;
my $headers = shift;
my $plotID = shift;
my $yheader = shift;
my $forceMin = shift;
my $SE = shift;
my $n_repl = shift;

if($n<=0 || $m<=0){ error_message("Output file is empty!"); }
my $sort_histogram = $hashInput{"sort_histogram"};
if($sort_histogram eq "on" && $m==1){
	my @index = sort {$data->[$b]<=>$data->[$a]} 0..($n-1);
	my @data1; my @headers1;
	for(my $i=0; $i<$n; $i++){
		$data1[$i] = $data->[$index[$i]];
		$headers1[$i] = $headers->[$index[$i]];
	}
	$data = \@data1;
	$headers = \@headers1;
	if($SE){
		my @SE1;
		for(my $i=0; $i<$n; $i++){
			$SE1[$i] = $SE->[$index[$i]];
		}
		$SE = \@SE1;
	}
	if($n_repl){
		my @n_repl1;
		for(my $i=0; $i<$n; $i++){
			$n_repl1[$i] = $n_repl->[$index[$i]];
		}
		$n_repl = \@n_repl1;
	}
}
my $n1 = $n;
if($n > 300){ $n1=300; }

# Find minimum and maximum
my $miny = 1.0e10;
my $maxy = -1.0e10;
my $median=median($data);
my $nn=0;
for(my $i=0; $i<$n1; ++$i){
	my $y = $data->[$i];
	if($m==1){
		if($y <= $MISSING){ next; }
		$nn++;
		my $y1 = $y;
		if($SE){
			$y += $SE->[$i];
			$y1 -= $SE->[$i];
		}
		if($miny > $y1) { $miny=$y1; }
		if($maxy < $y) { $maxy=$y; }
	}else{
		if(ref($y) ne 'ARRAY'){ next; }
		if($plot_type eq "stacked"){
			if($y->[0]<0 || $y->[1]<0){ error_message("Stacked: negative"); }
			if($maxy < $y->[0]+$y->[1]) { $maxy=$y->[0]+$y->[1]; }
			next;
		}
		for(my $j=0; $j<$m; $j++){
			my $y1 = $y->[$j];
			if($y1 <= $MISSING){ next; }
			if($maxy < $y1) { $maxy=$y1; }
			if($miny > $y1) { $miny=$y1; }
		}
	}
}
if($miny==$maxy){
	if($maxy>0){ $miny=0; }
	else{ $maxy=0; }
}
if($forceMin){
	if($forceMin eq "zero" && $miny>0){ $miny=0; }
	elsif($miny>$forceMin){ $miny=$forceMin; }
}
if($median==$MISSING){ $median=$miny; }

my $spany = $maxy-$miny;
if($spany <= 0) { $spany=0.1; }
my $deltay;
if($spany < 1){ $deltay = 0.1; }
elsif($spany < 2) { $deltay = 0.2; }
elsif($spany < 5) { $deltay = 0.5; }
elsif($spany < 10) { $deltay = 1; }
elsif($spany < 20) { $deltay = 2; }
elsif($spany < 50) { $deltay = 5; }
elsif($spany < 100) { $deltay = 10; }
elsif($spany < 200) { $deltay = 20; }
elsif($spany < 500) { $deltay = 50; }
elsif($spany < 1000) { $deltay = 100; }
elsif($spany < 2000) { $deltay = 200; }
elsif($spany < 5000) { $deltay = 500; }
else { $deltay = 1000; }

if($miny >= 0) { $miny = int($miny/$deltay)*$deltay; }
else { $miny = int($miny/$deltay-1)*$deltay; }
$spany = (int(($maxy-$miny)/$deltay)+1)*$deltay;
if($spany <= 0) { $spany=0.1; }

my $drawID = $plotID+1;
my $max_length = 0;
my @headers1 = @$headers;
for(my $i=0; $i<$n1; ++$i){
	my $l = length($headers1[$i]);
	if($l > $maxHeaderLength){
		$headers1[$i] = substr($headers1[$i],0,$maxHeaderLength);
		$l = $maxHeaderLength;
	}
	if($max_length < $l){
		$max_length = $l;
	}
}
if(!open (OUT1, ">$PATH_OUTPUT/$drawID.txt")){
	print "Cannot write file $drawID.txt<p>\n";
}
my $SPACE;
if($n <= 10){ $SPACE = 30; }
elsif($n < 30){ $SPACE = 20; }
elsif($n < 70){ $SPACE = 15; }
elsif($n < 150){ $SPACE = 10; }
else{ $SPACE = int(1800/$n1); }

my $HGT = $max_length*6 + 250 + 35;
my $WID = $n1*$SPACE+70;
my $scale = 250/$spany;
my $y_center = $HGT/2-$max_length*6-20;
my $x_center = $WID/2-60;
my $nobjects = $n*10+10000;
print OUT1 "1\n$x_center\n$y_center\n$nobjects\n";

# Axis vertical
my $border = 6;
my $y1 = 0;
my $y2 = int($spany*$scale);
my $x1 = 0;
my $x2 = $n1*$SPACE;
print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x2 $y1\n";
print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x1 $y2\n";

# Make ticks on the vertical axis
my $x3 = $x1-3;
my $x4 = $x1-20;
for(my $i = 0; $i <= $spany/$deltay; ++$i) {
	my $y = $i*$deltay+$miny;
	my $y3 = int(($y-$miny)*$scale);
	if($i){ print OUT1 "$LINE $gray $thin $solid  $x1 $y3 $x2 $y3\n"; }
	print OUT1 "$LINE $black $thin $solid  $x3 $y3 $x1 $y3\n";
	print OUT1 "$TEXT $black $smallFont $x4 $y3 $y\n";
}
if($miny < 0){
	my $y3 = int(-$miny*$scale);
	print OUT1 "$LINE $black $thin $solid  $x1 $y3 $x2 $y3\n";
}
my $y5 = $y2 + 12;
if($yheader){
	print OUT1 "$TEXT $black $vertSmallFont -48 50 $yheader\n";
}
my @table2=();
my $y0 = int(-$miny*$scale);
if($y0 < 0){ $y0=0; }
if($plot_type eq "deviation"){
	$y0 = int(($median-$miny)*$scale);
	print OUT1 "$LINE $black $thin $solid  $x1 $y0 $x2 $y0\n";
}
for(my $i=0; $i<$n; ++$i){
	my $x1 = int(($i+0.2)*$SPACE);
	my $x2 = int(($i+0.8)*$SPACE);
	my $x4 = int(($x1+$x2)/2);
	my $x5 = $x4-1;
	my $x6 = $x4+1;
	if($i<$n1){
		my $y2 = -(length($headers1[$i])*6+5);
		my $x3 = $x1+4;
		print OUT1 "$TEXT $black $vertSmallFont $x3 $y2 $headers1[$i]\n";
	}
	my $y = $data->[$i];
	if($m>1 && $plot_type ne "stacked"){
		if(ref($y) ne 'ARRAY'){ next; }
		my $line = "<tr><td>$headers->[$i]";
		for(my $j=0; $j<$m; $j++){
			my $y1 = $y->[$j];
			if($y1 <= $MISSING){ $line .= "<td><center>n/a"; }
			$line .= "<td><center>$y1";
			if($i>=$n1 || $y1 <= $MISSING){ next; }
			my $y3 = int(($y1-$miny)*$scale);
			my $ic = $i%@colors;
			print OUT1 "$CIRCLE $black $colors[$ic] 3 $x4 $y3\n";
		}
		push(@table2,$line);
		next;
	}
	my ($yy1,$yy2)=(0,0);
	my $replic = "";
	if($n_repl){ $replic = "<td><center>$n_repl->[$i]"; }
	if($plot_type eq "stacked"){
		if(ref($y) eq 'ARRAY'){
			($yy1,$yy2) = ($y->[0],$y->[1]);
			$y = $yy1+$yy2;
		}else{
			$y = 0;
		}
	}
	if($y <= $MISSING){
		push(@table2,"<tr><td>$headers->[$i]<td><center>none$replic");
		next;
	}
	$y = int(10000*$y)/10000;
	push(@table2,"<tr><td>$headers->[$i]<td><center>$y$replic");
	if($i >= $n1){ next; }
	my $y3 = int(($y-$miny)*$scale);
	if($plot_type eq "bars"){
		print OUT1 "$BOX $black $ltgreen $x1 $y0 $x2 $y3\n";
	}elsif($plot_type eq "points"){
		print OUT1 "$CIRCLE $black $blue 3 $x4 $y3\n";
	}elsif($plot_type eq "deviation"){
		print OUT1 "$BOX $black $ltgreen $x1 $y0 $x2 $y3\n";
	}elsif($plot_type eq "stacked"){
		my $y4 = int(($yy1-$miny)*$scale);
		print OUT1 "$BOX $black $orange $x1 $y0 $x2 $y4\n";
		print OUT1 "$BOX $black $blue $x1 $y4 $x2 $y3\n";
	}
	if($SE && $SPACE>3){
		my $y1 = int(($data->[$i] + $SE->[$i] - $miny)*$scale);
		print OUT1 "$LINE $black $thin $solid  $x4 $y3 $x4 $y1\n";
		print OUT1 "$LINE $black $thin $solid  $x5 $y1 $x6 $y1\n";
		if($plot_type eq "deviation"){
			$y1 = int(($data->[$i] - $SE->[$i] - $miny)*$scale);
			print OUT1 "$LINE $black $thin $solid  $x4 $y3 $x4 $y1\n";
			print OUT1 "$LINE $black $thin $solid  $x5 $y1 $x6 $y1\n";
		}
	}
}
close (OUT1);
my $outputFile = "$PATH_OUTPUT/$plotID.gif";
system("$PATH_BIN/togif","$PATH_OUTPUT/$drawID.txt","$outputFile","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");
return(@table2);
}

#**************************************
sub  plot_matrix
#**************************************
{
my $fileID = shift;
my $matrix_name = shift;
my $logFileID = shift;

my $webPageID;
if($hashInput{"logFileID"} && $fileID){
	$webPageID = $hashInput{"logFileID"}+1;
}
if(!$fileID){ $fileID = $hashInput{"fileID"}; }
if(!$matrix_name){ $matrix_name = $hashInput{"matrix_name"}; }
my $display = $hashInput{"display"};
if($display){
	plot_profile($fileID,$display);
}
my $legendID = $hashInput{"legendID"};
my $scriptID = get_outputID(4);
my $imageID = $scriptID+1;
my $scaleScriptID = $scriptID+2;
my $scaleImageID = $scriptID+3;
my $autocorrelation=0;
my $file_output = $hashInput{"file_output"};
my $organismID = $hashInput{"organismID"};
my $imatrix = $hashInput{"matrix_plot"};
my $file_output_title = $hashInput{"file_output_title"};
my $description = $hashInput{"description"};
if($file_output){
	my %hashOutput;
	my $file_full = "$PATH_DATA/$file_output";
	if($file_output=~/^\d+\.txt/){ $file_full = "$PATH_OUTPUT/$file_output"; }
	elsif($file_output !~ /^public-/){ $file_full = "$PATH_DATA/$loginname-$file_output"; }
	parse_file_headers($file_full,\%hashOutput);
	my @data_name = ($hashOutput{"output_data_set1"},$hashOutput{"output_data_set2"});
	$autocorrelation = $hashOutput{"output_autocorrelation"};
	if(!$file_output_title){ $file_output_title=$file_output; }
	if(!$description){ $description=$hashOutput{"output_description"}; }
}else{
	$file_output = $hashInput{"file_matrix"};
	if(!$description){ $description = $hashInput{"description_matrix"}; }
}
open (INFO,"<$PATH_OUTPUT/$fileID.txt") or error_message("In plot matrix - cannot open $fileID.txt");
my $line = <INFO>;
chop $line;
my ($headerTop,@headers) = split(/\t/,$line);
my $nCol = @headers;
if(!$nCol){ terminal_window("<H3>Matrix has no data</H3>"); }
my $long_col_header=0;
for(my $i=0; $i<$nCol; ++$i){
	my $len = length($headers[$i]);
	if($len>$maxHeaderLength){ $headers[$i]=substr($headers[$i],0,$maxHeaderLength); $len=$maxHeaderLength; }
	if($long_col_header<$len){ $long_col_header=$len; }
}
my $nRow = 0;
my @data;
my $long_row_header = 0;
my @rows=();
my @abs_x;
while(my $line = <INFO>){
	if($line =~ /^!/){ last; }
	chop $line;
	my ($row_header,@data1)=split(/\t/, $line);
	my $len = length($row_header);
	if($len>$maxHeaderLength){ $row_header=substr($row_header,0,$maxHeaderLength); $len=$maxHeaderLength; }
	if($long_row_header<$len){ $long_row_header=$len; }
	for(my $i=0; $i<@data1; ++$i){
		my $x=$data1[$i];
		if($x==$MISSING || $autocorrelation && $x==1000){ next; }
		push(@abs_x, $x);
	}
	push(@data,\@data1);
	push(@rows,$row_header);
	++$nRow;
}
close INFO;
@abs_x = sort {$a<=>$b} @abs_x;
my $xmax = int(100*$abs_x[int(0.98*@abs_x)])/100;
if($xmax<0.5){ $xmax=0.5; }
my $maximum_value = $hashInput{"maximum_value"};
my $plot_scale = $hashInput{"plot_scale"};
my $no_black_lines=0;
if($hashInput{"no_black_lines"} eq "on"){ $no_black_lines=1; }
if($maximum_value > 0){ $xmax=$maximum_value; }
my @THRESH;
my $Nlevels=11;
my $logN = log($Nlevels);
for(my $i=0; $i<$Nlevels; ++$i){
	if($plot_scale eq "log"){
		$THRESH[$i] = $xmax*exp(0.5*($i+1-$Nlevels)/$logN);
	}else{
		$THRESH[$i] = $xmax*($i+1)/$Nlevels;
	}
}
#DRAWING VALUE SCALE
my $WID = 250;
my $HGT = 65;
my $x_center = 125;
my $y_center = 37;
my $nobjects = 1000;
my $scale=1;
open (OUT1, ">$PATH_OUTPUT/$scaleScriptID.txt");
print OUT1 "$scale\n$x_center\n$y_center\n$nobjects\n";
for(my $i=0; $i<$Nlevels; ++$i){
	my $x = $THRESH[$i];
	if($x>1000){ $x=int($x+0.5); }
	elsif($x>10){ $x=int(100*$x+0.5)/100; }
	elsif($x>0.1){ $x=int(10000*$x+0.5)/10000; }
	elsif($x>0.001){ $x=int(1000000*$x+0.5)/1000000; }
	my ($x1,$y1) = (131+$i*10, 28);
	print OUT1 "$TEXT $black $vertTinyFont $x1 $y1 $x\n";
	$x1 = 121-$i*10;
	$y1 = 23;
	print OUT1 "$TEXT $black $vertTinyFont $x1 $y1 -$x\n";
	$x1 = 130+$i*10;
	my $x2 = $x1+10;
	$y1 = 10;
	my $color = 57+$i;
	for(my $i1=0; $i1<10; $i1++){
		print OUT1 "$LINE $color $thin $solid  $x1 $y1 $x2 $y1\n";
		$y1++;
	}
	$x1 = 110-$i*10;
	my $x2 = $x1+10;
	$y1 = 10;
	$color = 56-$i;
	for(my $i1=0; $i1<10; $i1++){
		print OUT1 "$LINE $color $thin $solid  $x1 $y1 $x2 $y1\n";
		$y1++;
	}
}
close OUT1;
system("$PATH_BIN/togif","$PATH_OUTPUT/$scaleScriptID.txt","$PATH_OUTPUT/$scaleImageID.gif","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");

#REORDERING COLUMNS and ROWS
my $matrix_modified = 0;
my $select_column = $hashInput{"select_column_heatmap"};
my $move_column = $hashInput{"move_column"};
if($select_column){
	$matrix_modified = 1;
	$select_column--;
	if($move_column eq "end"){
		foreach my $ref (@data){
			my $x = splice(@$ref,$select_column,1);
			push(@$ref,$x);
		}
		my $x = splice(@headers,$select_column,1);
		push(@headers,$x);
	}elsif($move_column eq "delete"){
		foreach my $ref (@data){
			splice(@$ref,$select_column,1);
		}
		splice(@headers,$select_column,1);
		$nCol--;
	}else{
		$move_column--;
		if($move_column > $select_column){ $move_column--; }
		foreach my $ref (@data){
			my $x = splice(@$ref,$select_column,1);
			splice(@$ref,$move_column,0,$x);
		}
		my $x = splice(@headers,$select_column,1);
		splice(@headers,$move_column,0,$x);
	}
}
$select_column = $hashInput{"select_column_heatmap"};
my $select_row = $hashInput{"select_row_heatmap"};
my $move_row = $hashInput{"move_row"};
if($select_row && !$autocorrelation || $select_column && $autocorrelation){
	$matrix_modified = 1;
	if($autocorrelation && !$select_row){
		$select_row = $select_column;
		$move_row = $hashInput{"move_column"};
	}
	$select_row--;
	if($move_row eq "end"){
		my $x = splice(@data,$select_row,1);
		push(@data,$x);
		$x = splice(@rows,$select_row,1);
		push(@rows,$x);
	}elsif($move_row eq "delete"){
		splice(@data,$select_row,1);
		splice(@rows,$select_row,1);
		$nRow--;
	}else{
		$move_row--;
		if($move_row > $select_row){ $move_row--; }
		my $x = splice(@data,$select_row,1);
		splice(@data,$move_row,0,$x);
		$x = splice(@rows,$select_row,1);
		splice(@rows,$move_row,0,$x);
	}
}
my $start_column = int($hashInput{"start_column"});
my $end_column = int($hashInput{"end_column"});
my $n_range = $end_column-$start_column+1;
if($n_range < 1){ error_message("Wrong range of columns $start_column-$end_column"); }
if($start_column){
	$matrix_modified = 1;
	$start_column--;
	my $move_column = $hashInput{"move_column"};
	if($move_column eq "end"){
		foreach my $ref (@data){
			my @x = splice(@$ref,$start_column,$n_range);
			push(@$ref,@x);
		}
		my @x = splice(@headers,$start_column,$n_range);
		push(@headers,@x);
	}elsif($move_column eq "delete"){
		foreach my $ref (@data){
			splice(@$ref,$start_column,$n_range);
		}
		my @x = splice(@headers,$start_column,$n_range);
		$nCol -= $n_range;
	}else{
		$move_column--;
		if($move_column > $end_column){ $move_column -= $n_range; }
		foreach my $ref (@data){
			my @x = splice(@$ref,$start_column,$n_range);
			splice(@$ref,$move_column,0,@x);
		}
		my @x = splice(@headers,$start_column,$n_range);
		splice(@headers,$move_column,0,@x);
	}
}
$start_column = int($hashInput{"start_column"});
my $start_row = int($hashInput{"start_row"});
my $end_row = int($hashInput{"end_row"});
$n_range = $end_row-$start_row+1;
if($n_range < 1){ error_message("Wrong range of rows $start_row-$end_row"); }
$move_row = $hashInput{"move_row"};
if($start_row && !$autocorrelation || $start_column && $autocorrelation){
	$matrix_modified = 1;
	if($autocorrelation && !$start_row){
		$start_row = $start_column;
		$end_row = $end_column;
		$n_range = $end_row-$start_row+1;
		$move_row = $hashInput{"move_column"};
	}
	$start_row--;
	if($move_row eq "end"){
		my @x = splice(@data,$start_row,$n_range);
		push(@data,@x);
		@x = splice(@rows,$start_row,$n_range);
		push(@rows,@x);
	}elsif($move_row eq "delete"){
		splice(@data,$start_row,$n_range);
		splice(@rows,$start_row,$n_range);
		$nRow -= $n_range;
	}else{
		$move_row--;
		if($move_row > $start_row){ $move_row -= $n_range; }
		my @x = splice(@data,$start_row,$n_range);
		splice(@data,$move_row,0,@x);
		@x = splice(@rows,$start_row,$n_range);
		splice(@rows,$move_row,0,@x);
	}
}
if($matrix_modified){
	# Update modified matrix
	open (OUT, ">$PATH_OUTPUT/$fileID.txt");
	print OUT join("\t",$headerTop,@headers)."\n";
	for(my $i=0; $i<$nRow; ++$i){
		print OUT join("\t",$rows[$i],@{$data[$i]})."\n";
	}
	close OUT;
}
my $span_col = 12;
if($nRow > 150 && $nCol < 50){
	$span_col = int(500/$nCol);
	if($span_col > 50){ $span_col = 50; }
}
my $span_row = 10;
my $combine = 1;
my $nRowCombined = $nRow;
if($nCol > 150){
	$span_col = int(1500/$nCol);
	if($span_col < 1){ $span_col=1; }
}
$WID = $nCol*$span_col+$long_row_header*6+50;
$x_center = int($nCol*$span_col/2)-$long_row_header*3-5;
if($nRow > 150){
	$span_row = int(1500/$nRow);
	if($span_row<1){ $span_row=1; }
	if($span_row*$nRow > 1500){
		$combine = int($nRow/1200)+1;
		$nRowCombined = int($nRow/$combine)+1;
	}
	$WID = $nCol*$span_col+50;
	$x_center = int($nCol*$span_col/2)-10;
}
$HGT = $nRowCombined*$span_row+$long_col_header*6+50;
$y_center = int($nRowCombined*$span_row/2)+$long_col_header*3+5;
$nobjects = 3*$nCol*$nRowCombined*$span_row+10000;
$scale=1;

open (OUT1, ">$PATH_OUTPUT/$scriptID.txt");
print OUT1 "$scale\n$x_center\n$y_center\n$nobjects\n";
my $y0 = $nRowCombined*$span_row+5;
for (my $i = 0; $i < $nCol; ++$i){
	my $x1 = int($span_col*($i+0.5)-2);
	print OUT1 "$TEXT $black $vertSmallFont $x1 $y0 $headers[$i]\n";
}
my $max_color=0;
for(my $i=0; $i<$nCol; ++$i){
	my $x1 = $i*$span_col;
	my $x2 = $x1+$span_col;
	for(my $j=0; $j<$nRowCombined; ++$j){
		my $y1 = int(($nRowCombined-$j-1)*$span_row);
		my $x = $data[$j]->[$i];
		if($combine>1){
			my ($sum,$n) = (0,0);
			for(my $k=0; $k<$combine; ++$k){
				my $k1 = $j*$combine+$k;
				if($k1 >= $nRow){ last; }
				my $x1 = $data[$k1]->[$i];
				if($x1==$MISSING){ next; }
				$sum += $x1;
				$n++;
			}
			if(!$n){ next; }
			$x = $sum/$n;
		}
		if($x==$MISSING){ next; }
		my $color=0;
		my $icolor=0;
		if($x<0){
			while($icolor<$Nlevels && $x<-$THRESH[$icolor]){ $icolor++; }
			if($icolor){ $color = 57-$icolor; }
		}else{
			if($max_color < $THRESH[$icolor]){ $max_color = $THRESH[$icolor]; }
			while($icolor<$Nlevels && $x>$THRESH[$icolor]){ $icolor++; }
			if($icolor){ $color = 56+$icolor; }
		}
		for(my $i=0; $i<$span_row; $i++){
			print OUT1 "$LINE $color $thin $solid  $x1 $y1 $x2 $y1\n";
			$y1++;
		}
	}
}
if($nRow <= 150){
	my $xxx = $nCol*$span_col;
	if($no_black_lines==0){
		for (my $i = 0; $i <= $nRow; ++$i){
			my $y1 = $span_row*$i;
			print OUT1 "$LINE $black $thin $solid  0 $y1 $xxx $y1\n";
		}
		my $yyy = $nRow*$span_row;
		for (my $i = 0; $i <= $nCol; ++$i){
			my $x1 = $span_col*$i;
			my $y1 = $yyy+5;
			print OUT1 "$LINE $black $thin $solid  $x1 0 $x1 $yyy\n";
		}
	}
	for (my $i = 0; $i < $nRow; ++$i){
		my $y1 = int(($nRow-$i-0.5)*$span_row);
		my $len = length($rows[$i]);
		my $x1 = -$len*3-6;
		print OUT1 "$TEXT $black $smallFont $x1 $y1 $rows[$i]\n";
	}
}else{
	my $x1 = -20;
	my $x2 = -3;
	for(my $j=0; $j<$nRowCombined; ++$j){
		my $y1 = int(($nRowCombined-$j-0.5)*$span_row);
		print OUT1 "$LINE $gray $thin $solid  $x1 $y1 $x2 $y1\n";
	}
}
close (OUT1);
system("$PATH_BIN/togif","$PATH_OUTPUT/$scriptID.txt","$PATH_OUTPUT/$imageID.gif","$PATH_BIN/raster.par","$PATH_BIN/palette.par","$WID","$HGT");

my $text = "<HTML><HEAD><TITLE>ExAtlas matrix</TITLE>\n";
$text .= get_header();
$text .= "<SCRIPT LANGUAGE=\"JavaScript\">\n";
$text .= "<!--\n";
$text .= "function check_form(){\n";
$text .= "	document.modify_matrix.display.value=\"\";\n";
$text .= "	document.modify_matrix.target = \"\";\n";
$text .= "	if(document.modify_matrix.select_column_heatmap.selectedIndex == document.modify_matrix.move_column.selectedIndex+1){\n";
$text .= "		alert(\"You cannot move the column to the same place.\"); return false;\n";
$text .= "	}\n";
$text .= "	if(document.modify_matrix.select_row_heatmap.selectedIndex == document.modify_matrix.move_row.selectedIndex+1){\n";
$text .= "		alert(\"You cannot move the row to the same place.\"); return false;\n";
$text .= "	}\n";
$text .= "	if(document.modify_matrix.start_column.value !=\"\" && document.modify_matrix.select_column_heatmap.selectedIndex!=0){\n";
$text .= "		alert(\"You cannot specify columns by both pull-down menu and range. Use one of these methods.\"); return false;\n";
$text .= "	}\n";
$text .= "	if(document.modify_matrix.start_row.value !=\"\" && document.modify_matrix.select_row_heatmap.selectedIndex!=0){\n";
$text .= "		alert(\"You cannot specify row by both pull-down menu and range. Use one of these methods.\"); return false;\n";
$text .= "	}\n";
$text .= "	var startcol = new Number(document.modify_matrix.start_column.value);\n";
$text .= "	var endcol = new Number(document.modify_matrix.end_column.value);\n";
$text .= "	var startrow = new Number(document.modify_matrix.start_row.value);\n";
$text .= "	var endrow = new Number(document.modify_matrix.end_row.value);\n";
$text .= "	var maxrow = new Number($nRow);\n";
$text .= "	var maxcol = new Number($nCol);\n";
$text .= "	if(isNaN(startcol) || document.modify_matrix.start_column.value && startcol<1 || endcol > maxcol || startcol>endcol){\n";
$text .= "		alert(\"Error in the column range\"); return false;\n";
$text .= "	}\n";
$text .= "	if(isNaN(startrow) || document.modify_matrix.start_row.value && startrow<1 || endrow > maxrow || startrow>endrow){\n";
$text .= "		alert(\"Error in the row range\"); return false;\n";
$text .= "	}\n";
$text .= "	if(document.modify_matrix.move_column.selectedIndex+1 >= startcol && document.modify_matrix.move_column.selectedIndex+1 <= endcol){\n";
$text .= "		alert(\"Error: The target column is within the range of columts to be moved\"); return false;\n";
$text .= "	}\n";
$text .= "	if(document.modify_matrix.move_row.selectedIndex+1 >= startrow && document.modify_matrix.move_row.selectedIndex+1 <= endrow){\n";
$text .= "		alert(\"Error: The target row is within the range of rows to be moved\"); return false;\n";
$text .= "	}\n";
$text .= "	document.modify_matrix.submit();\n";
$text .= "}\n";
$text .= "function column(text){\n";
$text .= "	document.modify_matrix.display.value=\"column \"+text;\n";
$text .= "	var x = Math.round(Math.random()*10000);\n";
$text .= "	document.modify_matrix.target = \"_BLANK\"+x;\n";
$text .= "	document.modify_matrix.submit();\n";
$text .= "}\n";
$text .= "function row(text){\n";
$text .= "	document.modify_matrix.display.value=\"row \"+text;\n";
$text .= "	var x = Math.round(Math.random()*10000);\n";
$text .= "	document.modify_matrix.target = \"_BLANK\"+x;\n";
$text .= "	document.modify_matrix.submit();\n";
$text .= "}\n";
$text .= "// -->\n";
$text .= "</SCRIPT>\n";
$text .= "<FORM NAME=modify_matrix ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
my $ending1 = "for file '$file_output_title'";
if(!$file_output_title){ $ending1 = ""; }
$text .= "<h2>Heatmap $ending1</h2>\n";
if($matrix_name){
	$text .= "<b>Matrix name:</b> $matrix_name<br>\n";
}
$text .= "<b>Description:</b> $description &nbsp; &nbsp; <b>N columns:</b> $nCol &nbsp; &nbsp; <b>N rows:</b> $nRow<br>\n";
if($legendID){
	my $x = int(10000*rand());
	$text .= "<b>Legend:</b> <a href=$HOME_ADDRESS/output/$legendID.txt target=_blank$x>Legend file</a><br>\n";
}
my $x = int(10000*rand());
$text .= "Link to download output matrix as text: <a href=$HOME_ADDRESS/output/$fileID.txt target=_blank$x>Matrix_file</a><p>\n";
$text .= "If row or column headers are not seen, move the mouse over them and read at the bottom of the browser.<br>\n";
$text .= "<INPUT NAME=sort_histogram TYPE=checkbox> Sort profiles (if you click on headers of columns or rows)<p>\n";
$text .= "<IMG SRC=../output/$imageID.gif BORDER=0 useMap=#mapMatrix><br>\n";
$text .= "<IMG SRC=../output/$scaleImageID.gif BORDER=0><p>\n";
$text .= "<map NAME='mapMatrix'>\n";
my ($y1,$y2) = (10,$long_col_header*6+20);
for (my $i = 0; $i<$nCol; ++$i){
	my $x1 = $span_col*$i+30;
	if($nRow <= 150){
		$x1 += $long_row_header*6-5;
	}
	my $x2 = $x1+$span_col;
	my $num = $i+1;
	$text .= "<area shape=rect coords=$x1,$y1,$x2,$y2 href=\"javascript:column('$num $headers[$i]');\">\n";
}
if($nRow <=150){
	my ($x1,$x2) = (10,$long_row_header*6+20);
	for(my $i=0; $i<$nRow; ++$i){
		my $y1 = $i*$span_row + $long_col_header*6+30;
		my $y2 = $y1+$span_row;
		my $num = $i+1;
		$text .= "<area shape=rect coords=$x1,$y1,$x2,$y2 href=\"javascript:row('$num $rows[$i]');\">\n";
	}
}else{
	my ($x1,$x2) = (10,30);
	for(my $i=0; $i<$nRowCombined; ++$i){
		my $y1 = $i*$span_row + $long_col_header*6+30;
		my $y2 = $y1+$span_row;
		my $num = $i*$combine+1;
		$text .= "<area shape=rect coords=$x1,$y1,$x2,$y2 href=\"javascript:row('$num $rows[$i*$combine]');\">\n";
	}
}
$text .= "</map>\n";

$text .= "<TABLE><TR><TD COLSPAN=2 ALIGN=center><b>Move or delete a column</b><TD WIDTH=30><TD COLSPAN=2 ALIGN=center><b>Move or delete a row</b>\n";
$text .= "<TR><TD WIDTH=180>Select Column<TD><SELECT NAME=select_column_heatmap style=width:250px;>\n";
$text .= "<OPTION VALUE=0>-------Select-------\n";
for(my $i=1; $i<=@headers; ++$i){
	$text .= "<OPTION VALUE=$i>$i. $headers[$i-1]\n";
}
$text .= "</SELECT>\n";
$text .= "<TD><TD WIDTH=180>Select Row<TD><SELECT NAME=select_row_heatmap style=width:250px;>\n";
$text .= "<OPTION VALUE=0>-------Select-------\n";
for(my $i=1; $i<=@rows; ++$i){
	$text .= "<OPTION VALUE=$i>$i. $rows[$i-1]\n";
}
$text .= "</SELECT>\n";
$text .= "<TR><TD>Or enter column range:<TD><INPUT NAME=start_column SIZE=3> - <INPUT NAME=end_column SIZE=3>\n";
$text .= "<TD><TD>Or enter row range:<TD><INPUT NAME=start_row SIZE=3> - <INPUT NAME=end_row SIZE=3>\n";
$text .= "<TR><TD>Move before column<TD><SELECT NAME=move_column style=width:250px;>\n";
for(my $i=1; $i<=@headers; ++$i){
	$text .= "<OPTION VALUE=$i>$i. $headers[$i-1]\n";
}
$text .= "<OPTION VALUE=end>Move to the end";
$text .= "<OPTION VALUE=delete>Delete column";
$text .= "</SELECT>\n";

$text .= "<TD><TD>Move before row<TD><SELECT NAME=move_row style=width:250px;>\n";
for(my $i=1; $i<=@rows; ++$i){
	$text .= "<OPTION VALUE=$i>$i. $rows[$i-1]\n";
}
$text .= "<OPTION VALUE=end>Move to the end";
$text .= "<OPTION VALUE=delete>Delete row";
$text .= "</SELECT>\n";
$text .= "<TR><TD>Remove black lines<TD><INPUT TYPE=checkbox NAME=no_black_lines>\n";
$text .= "<TR><TD>Maximum value<TD><INPUT NAME=maximum_value value=$xmax style=width:250px;>\n";
$text .= "<TD><TD ALIGN=RIGHT>Color scale<TD><SELECT NAME=plot_scale style=width:250px;>";
$text .= "<OPTION VALUE=uniform>Uniform scale";
$text .= "<OPTION VALUE=log"; if($plot_scale eq "log"){ $text .= " selected"; } $text .= ">Log scale";
$text .= "</SELECT>\n";
$text .= "<TR><TD>&nbsp;<TR><TD><INPUT TYPE=BUTTON NAME=modify_matrix value=\" Re-plot the matrix \" onClick=check_form();>\n";
$text .= "<TD><TD><TD><INPUT TYPE=button VALUE=\" Cancel (close window) \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
$text .= "</TABLE>\n";
$text .= "<INPUT NAME=sessionID TYPE=hidden VALUE=$sessionID>\n";
$text .= "<INPUT NAME=action TYPE=hidden VALUE=plot_matrix>\n";
$text .= "<INPUT NAME=display TYPE=hidden>\n";
$text .= "<INPUT NAME=file_output TYPE=hidden VALUE=$file_output>\n";
$text .= "<INPUT NAME=fileID TYPE=hidden VALUE=$fileID>\n";
$text .= "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
$text .= "<INPUT NAME=description TYPE=hidden VALUE=\"$description\">\n";
$text .= "<INPUT NAME=matrix_name TYPE=hidden VALUE=\"$matrix_name\">\n";
if($legendID){ $text .= "<INPUT NAME=legendID TYPE=hidden VALUE=$legendID>\n"; }
$text .= "</FORM>\n";
if(!$webPageID){ print $text; exit(0); }
if(!open(OUT,">$PATH_OUTPUT/$webPageID.txt")){
	error_message("cannot open web page ID",$logFileID);
}
print OUT $text;
close OUT;
return;
}

#**************************************
sub  plot_profile
#**************************************
{
my $fileID = shift;
my $display_text = shift;
my ($type,$num,@items) = split(/ /,$display_text);
$num =~ s/\.$//;
my $name = join(' ',@items);

my $matrix_name = $hashInput{"matrix_name"};
open(INFO,"<$PATH_OUTPUT/$fileID.txt") or error_message("In plot_profile - cannot open file");
my @extracted;
my @headers;
if($type eq "column"){
	while(my $line = <INFO>){
		chop $line;
		if($line =~ /^!/ || !$line){ next; }
		my ($row_header,@data)=split(/\t/, $line);
		push(@extracted,$data[$num-1]);
		push(@headers,$row_header);
	}
}elsif($type eq "row"){
	my $count = -1;
	while(my $line = <INFO>){
		chop $line;
		if($line =~ /^!/ || !$line){ next; }
		++$count;
		if($count == 0){
			@headers=split(/\t/, $line);
		}elsif($count == $num){
			@extracted=split(/\t/, $line);
			last;
		}
	}
}else{ error_message("Display type = $type"); }
close INFO;
if(!@extracted){ error_message("No data found"); }
if($extracted[0] ne $name){ error_message("Names don't match $extracted[0] ne $name"); }
shift(@extracted);
shift(@headers);
my $N = @extracted;
my $plotID = get_outputID(2);
my @table1 = plot_histogram_horiz("bars",$N,1,\@extracted,\@headers,$plotID);
print "<HTML><HEAD><TITLE>ExAtlas: Output profile</TITLE>\n";
print_header();
print "<p><B><font size=+2>Profile of $display_text</font></B> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;\n";
print "<INPUT TYPE=button VALUE=\" Close window \" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
my $file_output = $hashInput{"file_output"};
my $description = $hashInput{"description"};
my $organismID = $hashInput{"organismID"};
if($file_output !~ /^\d\d\d\d\d\d\d\d/){
	print "<b><u>File name:</u></b> $file_output<br>\n";
}
if($matrix_name){
	print "<b><u>Matrix name:</u></b> $matrix_name<br>\n";
}
print "<b><u>Description:</u></b> $description<p>\n";
print "<IMG SRC=../output/$plotID.gif BORDER=0><p>\n";
print "<TABLE>\n";
print "<TR><TD><b>Header</b><TD><b>Value</b>\n";
for(my $i=0; $i<@table1; ++$i){
	print $table1[$i]."\n";
}
print "</TABLE>\n";
print "<HR NOSHADE></HR><p>\n";
print "<INPUT TYPE=button VALUE=\" Close window \" LANGUAGE=\"javascript\" onClick=\"window.close();\">\n";
print "</BODY></HTML>\n";
exit(0);
}

#**************************************
sub  add_hyperlinks
#**************************************
{
my $word = shift;
if($word=~/PMID:/){
	my @items = split(/PMID:/,$word);
	$word = $items[0];
	$word=~s/[\s_]+$//;
	for(my $i=1; $i<@items; $i++){
		$items[$i]=~s/[\s_]+$//;
		$items[$i]=~s/^[\s_]+//;
		my @items1 = split(/\D/,$items[$i]);
		$items[$i] =~ s/^\d+//;
		$word .= " PMID:<a href=https://www.ncbi.nlm.nih.gov/pubmed/?term=$items1[0]%5Buid%5D>$items1[0]</a> $items[$i]";
	}
}
if($word=~/[\s_]*GSE\d/){
	my @items = split(/GSE/,$word);
	$word = $items[0];
	$word=~s/[\s_]+$//;
	for(my $i=1; $i<@items; $i++){
		$items[$i]=~s/[\s_]+$//;
		my @items1 = split(/\D/,$items[$i]);
		$items[$i] =~ s/^\d+//;
		$word .= " <a href=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE$items1[0]>GSE$items1[0]</a> $items[$i]";
	}
}
return $word;
}

#**************************************
sub   detect_up_down
#**************************************
{
my $ref=shift;
if(!$ref || ref($ref) ne 'ARRAY' || @$ref==0){ return 0; }
my @updown=(0,0);
foreach my $name (@$ref){
	if($name =~ /_up$/i){ $updown[0]++; }
	elsif($name =~ /_down$/i){ $updown[1]++; }
	else{ return 0; }
}
if($updown[0]*$updown[1]){ return(1); }
return 0;
}

#**************************************
sub   filter_list_by_organism
#**************************************
{
my $ref=shift;
my $taxid=shift;
my $n = @$ref;
for(my $i=$n-1; $i>=0; $i--){
	my $taxid1 = $ref->[$i]->[2];
	my $found=0;
	foreach my $id (split(/;/,$taxid1)){
		if($id==$taxid){ $found=1; }
	}
	if(!$found){
		splice(@$ref,$i,1);
	}
}
return;
}

#**************************************
sub  check_configuration
#**************************************
{
my $organismID = $hashInput{"organismID"};
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!");
open(OUT, ">$PATH_INFO/$loginname-config1.txt") or error_message("Cannot update configuration file!");
my %hashDefault=();
my $changed=0;
my $nLines;
while(my $line = <INFO>){
	if(length($line)<3){ next; }
	$nLines++;
	if($line =~ /^mainPage=$organismID\t/){
		chop $line;
		read_config_line($line,\%hashDefault);
	}else{
		print OUT $line;
	}
}
close INFO;
foreach my $key (keys %hashInput){
	if($key =~/^file_(matrix|geneset|samples|output)$/){
		my $value = $hashDefault{$key};
		if($value ne $hashInput{$key}){
			$hashDefault{$key} = $hashInput{$key};
			$changed=1;
		}
	}
}
if($changed){
	print OUT "mainPage=$organismID";
	foreach my $key (keys %hashDefault){
		if($key && $key ne "mainPage"){
			print OUT "\t$key=$hashDefault{$key}";
		}
	}
	print OUT "\n";
}
close OUT;
if($changed){
	my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
	if($nLines1 && $nLines1 > $nLines*0.9){
		copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
	}else{
		error_message("Failed to update configuration file!");
	}
}
unlink("$PATH_INFO/$loginname-config1.txt");
return;
}

#**************************************
sub  file_delete
#**************************************
{
my $filename = $hashInput{"file_delete"};
my $file_type = $hashInput{"file_type"};
my $warning;
if($filename=~/^public-/ && $loginname ne "public"){ error_message("Cannot delete public files"); }
if(file_exist("$PATH_DATA/$loginname-$filename")){
	unlink "$PATH_DATA/$loginname-$filename";
	if(file_exist("$PATH_DATA/$loginname-anova-$filename")){
		unlink "$PATH_DATA/$loginname-anova-$filename";
	}
}elsif(file_exist("$PATH_DATA/$loginname-$filename"."_annot.txt")){
	unlink("$PATH_DATA/$loginname-$filename"."_annot.txt");
}else{
	$warning = "File $filename not found";
}
open(INFO,"<$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!");
open(OUT, ">$PATH_INFO/$loginname-config1.txt") or error_message("Cannot update configuration file!");
my $nLines;
while(my $line = <INFO>){
	chop $line;
	my %hash=();
	if(length($line)<3){ next; }
	$nLines++;
	read_config_line($line,\%hash);
	if($file_type eq "matrix" && $filename eq $hash{"type_matrix"} ||
	 $file_type eq "geneset" && $filename eq $hash{"type_geneset"} ||
	 $file_type eq "output" && $filename eq $hash{"type_output"} ||
	 $file_type eq "samples" && $filename eq $hash{"type_samples"} ||
	 $file_type eq "search" && $filename eq $hash{"type_search"}  ||
	 $file_type eq "annotation" && $filename eq $hash{"type_annotation"}){ next; }
	print OUT $line."\n";
}
close INFO;
close OUT;
my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
if($nLines1 && $nLines1 > ($nLines-1)*0.9){
	copy "$PATH_INFO/$loginname-config1.txt", "$PATH_INFO/$loginname-config.txt";
}else{
	error_message("Failed to update configuration file!");
}
unlink("$PATH_INFO/$loginname-config1.txt");
my @filelist = glob("$PATH_INFO/guest*");
foreach my $filename (@filelist){
	if($filename =~ /(guest|default)-config.txt/){ next; }
	my ($junk,$date_user) = split(/\/guest/,$filename);
	$date_user = substr($date_user,0,6);
	if($date_user>0 && $date > $date_user+2){
		unlink "$filename";
	}
}
@filelist = glob("$PATH_DATA/guest*");
foreach my $filename (@filelist){
	if($filename =~ /(guest|default)-/){ next; }
	my ($junk,$date_user) = split(/\/guest/,$filename);
	$date_user = substr($date_user,0,6);
	if($date_user>0 && $date > $date_user+2){
		unlink "$filename";
	}
}
if($warning){ print "Warning: $warning<br>\n"; }
return;
}

#**************************************
sub  file_download
#**************************************
{
my $filename = $hashInput{"file_download"};
my $add_symbols = 0;
my %annotation;
my $file_anova = $filename;
$file_anova =~ s/^public-/public-anova-/;
if($hashInput{"select_column"} && $filename !~ /anova-/){ $add_symbols = 1; }
if($filename !~ /^public-/){
	$file_anova = "$loginname-anova-$filename";
	$filename = "$loginname-$filename";
}
my $nCol=0;
if($add_symbols){
	open(INFO1,"<$PATH_DATA/$file_anova");
	my $line = <INFO1>;
	chop $line;
	my ($junk,$junk1,@headers) =split(/\t/,$line);
	if(!@headers){ error_message("Headers not found in $file_anova"); }
	while($nCol<@headers && $headers[$nCol] !~ /^Var\(/){ ++$nCol; }
	close INFO1;
	$nCol += 8;
	open(INFO1,"<$PATH_DATA/$file_anova");
}
if(!open(INFO,"<$PATH_DATA/$filename")){
	if(!open(INFO,"<$PATH_DATA/$filename"."annot.txt")){
		error_message("File $filename not found!");
	}
}
print "<pre>";
while(my $line=<INFO>){
	chop $line;
	print $line;
	if($line=~/^!/){ print "\n"; next; }
	if($add_symbols){
		my $line = <INFO1>;
		chop $line;
		my ($junk,$junk1,@data)=split(/\t/,$line);
		splice(@data,0,$nCol);
		print "\t".join("\t",@data);
	}
	print "\n";
}
close INFO;
if($add_symbols){ close INFO1; }
exit(0);
}

#***********************************
sub file_read
#***********************************
{
if(!$_[0]){ return ""; }
if(open (INFO_TEMP,"<",$_[1])){
	my $output=<OUT_TEMP>; close OUT_TEMP; return $output; 
}
return "";
}

#***********************************
sub file_exist
#***********************************
{
if(open (INFO_TEMP,"<$_[0]")){ close INFO_TEMP; return 1; }
return 0;
}

#***********************************
sub file_append
#***********************************
{
if(!$_[0] || !$_[1]){ return; }
my $xx = ">>";
if($_[2]){ $xx = ">";  }
if(open (OUT_TEMP,$xx,$_[1])){ print OUT_TEMP $_[0]."\n"; close OUT_TEMP; return; }
}

#***********************
sub validate
#***********************
{
my $loginname = shift;
my $passwd = shift;

my $line;
my @abcd;
if($SECURITY<1){ return(0); }
if ($loginname =~ /^guest/){ return(0); }
if (!$loginname){
	error_message("You need to put your name in the form!","register");
}
if(!open(INFO,"<$PATH_INFO/login.txt")){
	error_message("Validation failed!","register");
}
my $metaTerm = quotemeta($loginname);
while($line = <INFO>){
	if($line =~ /^$metaTerm\t/){ last; };
}
close INFO;
@abcd = split (/\t/, $line);
my $encryptedPsw = $passwd;
for(my $i=0; $i<47; $i++){
	$encryptedPsw = crypt $encryptedPsw, $abcd[1];
}
if($loginname ne $abcd[0] || $encryptedPsw ne $abcd[1]){
	if($hashInput{"action"} eq "update_password"){
		error_message("You did not enter correct current password. Operation failed");
	}
	my $text = "Incorrect login name or password\n";
	if($line =~ /^$metaTerm\t/){
		$text .= "<FORM METHOD=POST ACTION=$CGI_ADDRESS/exatlas.cgi>";
		$text .= "To reset your password, enter your login name and email\n";
		$text .= "<TABLE BORDER=0><TR><TD>Login name:<TD><INPUT NAME=loginname SIZE=20>";
		$text .= "<TR><TD>Email:<TD><INPUT NAME=email_reset SIZE=20>";
		$text .= "<INPUT TYPE=submit NAME=reset_password VALUE=\"Reset password\">";
		$text .= "</TABLE></FORM>\n";
	}
	error_message($text,"register");
}
return(0);
}

#***********************
sub reset_password
#***********************
{
my $loginname = shift;

my $register="register";
my $email_reset = $hashInput{"email_reset"};
my $line;
my @abcd;
if(!$loginname || !$email_reset){
	error_message("Invalid loginname or email address",$register);
}
if(!open(INFO,"<$PATH_INFO/login.txt")){
	error_message("Operation faled1",$register);
}
if(!open(OUT, ">$PATH_INFO/login-temp.txt")){
	close INFO;
	error_message("Operation faled2",$register);
}
my @letters = ('A' .. 'Z', 'a' .. 'z', '0' .. '9', '_');
my $salt = $letters[rand@letters] . $letters[rand@letters];
my $content="";
my $email_address;
while($line = <INFO>){
	chop $line;
	if(!$line){ next; }
	@abcd = split (/\t/, $line);
	if(!$abcd[0]){ next; }
	if($abcd[0] eq $loginname){
		if($abcd[4] !~ /^\S+@\S+$/ || $abcd[4] ne $email_reset){
			error_message("Invalid loginname or email address; reset cancelled",$register);
		}
		$email_address = $abcd[4];
		my $new_passwd;
		for(my $i=0; $i<8; $i++){
			$new_passwd .= $letters[rand@letters];
		}
		my $encryptedPsw = $new_passwd;
		for(my $i=0; $i<47; $i++){
			$encryptedPsw = crypt $encryptedPsw, $salt;
		}
		$abcd[1]=$encryptedPsw;
		$content .= "Your password was reset:\nloginname=$loginname passwrd=$new_passwd\n\nAutomated mail. Do not reply!\n";
		$line = join("\t",@abcd);
	}
	print OUT $line."\n";
}
close INFO;
close OUT;
my $text;
if($content){
	`echo "$content" | mailx -s "ExAtlas reset" $email_address`;
	$text = "Password was reset and sent to you by email!";
	copy "$PATH_INFO/login-temp.txt", "$PATH_INFO/login.txt";
}else{
	$text .= "Invalid loginname or email address; reset cancelled\n";
}
unlink "$PATH_INFO/login-temp.txt";
terminal_window("<H3>$text</H3>",$register);
}

#***********************************
sub read_configuration
#***********************************
{
my $filename = shift;
my $hashIni_ref = shift;

if(!open (INFO_TEMP,'<',$filename)){
	print "content-type: text/html","\n\n";
	error_message("Configuration file not found.","register");
}
while(my $line=<INFO_TEMP>){
	$line =~ s/\n$//;
	my ($keyword,$value) = split(/=/,$line);
	$hashIni_ref->{$keyword} = $value;
}
close INFO_TEMP;
return;
}

#**************************************
sub substitute_chars
#**************************************
{
my $word1 = shift;
if(!defined($word1)){ return ""; }

$word1 =~ s/^\s+//; $word1 =~ s/\s+$//;
$word1 =~ s/\+/ /g;
if($word1 =~ /%/){
	if($word1 =~ /%09/){ s/%09/\t/g; }
	elsif($word1 =~ /%[01][0-9A-F]/){ $word1 =~ s/%[01][0-9A-F]//g; }
	$word1 =~ s/%C[EF]%[89AB][0-9A-F]/_/g;
	$word1 =~ s/%2B/\+/g;
	$word1 =~ s/%2F/\//g;
	$word1 =~ s/%3A/:/g;
	$word1 =~ s/%7E/\~/g;
	$word1 =~ s/%60/\`/g;
	$word1 =~ s/%21/\!/g;
	$word1 =~ s/%23/\#/g;
	$word1 =~ s/%24/\$/g;
	$word1 =~ s/%5E/\^/g;
	$word1 =~ s/%26/\&/g;
	$word1 =~ s/%3D/=/g;
	$word1 =~ s/%22/"/g;
	$word1 =~ s/%27/'/g;
	$word1 =~ s/%3B/;/g;
	$word1 =~ s/%3F/\?/g;
	$word1 =~ s/%5C/\\/g;
	$word1 =~ s/%7C/|/g;
	$word1 =~ s/%3C/</g;
	$word1 =~ s/%3E/>/g;
	$word1 =~ s/%2D/-/g;
	$word1 =~ s/%2C/,/g;
	$word1 =~ s/%40/@/g;
	$word1 =~ s/%25/\%/g;
	$word1 =~ s/%28/\(/g;
	$word1 =~ s/%29/\)/g;
	$word1 =~ s/%AD/-/g;
	$word1 =~ s/%BO/o/g;
	$word1 =~ s/%B1/+-/g;
	$word1 =~ s/%B2/\^2/g;
	$word1 =~ s/%B3/\^3/g;
	$word1 =~ s/%B4/\'/g;
	$word1 =~ s/%B5/mu/g;
	$word1 =~ s/%B7/\*/g;
	$word1 =~ s/%B9/\^1/g;
	if($word1 =~ /%[0-9A-F][0-9A-F]/){
		$word1 =~ s/%[0-9A-F][0-9A-F]/_/g;
	}
}
return $word1;
}

#**************************************
sub register_new_user
#**************************************
{
my $firstname = uc $hashInput{"firstname"};
my $lastname = uc $hashInput{"lastname"};
my $new_name = $hashInput{"loginname"};
my $new_passwd = $hashInput{"passwd"};

if($new_name =~ /^guest|^public/i){
	error_message("This is not a valid login name (should not start with \'guest\' or \'public\'.\nSelect another name or login as \"guest\" without password.","register");
}
my $email = $hashInput{"email"};
if($SECURITY>1 && $new_name ne "administrator"){
	open(INFO,"<$PATH_INFO/login.txt") or error_message("Administrator should register first!","register");
	my $line;
	while($line = <INFO>){
		if($line =~ /^administrator\t/){ last; };
	}
	close INFO;
	if ($line !~ /^administrator\t/) {
		error_message("Administrator should register first!","register");
	}
}
if(!$firstname || !$lastname || !$new_name || !$new_passwd || !$email){
	error_message("Form not fully filled. Registration failed!","register");
}
if(length($firstname)>20 || length($lastname)>20 || length($new_name)>20 || length($new_passwd)>20 || length($email)>30){
	error_message("Too long strings!\nNames should be <20 characters, email <30 characters.\nRegistration failed!","register");
}
if(length($new_name)<4 || length($new_passwd)<4){
	error_message("Too short loginname or password!\nMinimum length is 5 characters.\nRegistration failed!","register");
}
if($new_name eq $new_passwd){
	error_message("Password matches loginname!\nPassword not acceptable","register");
}
my $passed=1;
my @items=split(/\@/,$email); if (@items!=2){ $passed=0; }
foreach my $word (split(/\./,$items[0].".".$items[1])){
	if($word !~ /^[-\w]/ || $word=~/^-|-$/){ $passed=0; }
}
if(!$passed){ error_message("Email not valid","register"); }
if(open(INFO,"<$PATH_INFO/login.txt")){
	my $line;
	my $metaTerm = quotemeta($new_name);
	while(my $line = <INFO>){
		my @items = split(/\t/,$line);
		if ($new_name eq $items[0]) {
			close INFO;
			error_message("Login name is already in use! Select another.","register");
		}
	}
	close INFO;
}
my @letters = ('A' .. 'Z', 'a' .. 'z', '0' .. '9', '_', '.');
my $salt = $letters[rand@letters] . $letters[rand@letters];
my $encryptedPsw = $new_passwd;
for(my $i=0; $i<47; $i++){
	$encryptedPsw = crypt $encryptedPsw, $salt;
}
file_append("$new_name\t$encryptedPsw\t$lastname\t$firstname\t$email\t$date","$PATH_INFO/login.txt");
$sessionID = start_session($new_name);
$loginname = $new_name;
$passwd = $new_passwd;
copy "$PATH_DATA/default-config.txt", "$PATH_INFO/$loginname-config.txt";
terminal_window("<H3>Your registration is successful</H3>","continue");
}

#**************************************
sub  read_config_line
#**************************************
{
my $line = shift;
my $hash_ref = shift;

$line =~ s/\n$//;
%$hash_ref=();
my @items = split(/\t/,$line);
foreach my $item (@items){
	my($key,$value) = split(/=/,$item);
	if($key && $value){
		$hash_ref->{$key}=$value;
	}
}
return;
}

#***********************************
sub plot_line
#***********************************
# variable: 0=x1 1=x2 2=color 3=width 4=y_shift 5=scale 6=direction 7=number1 8=number2
{
	my $y1;
	my $x1 = int($_[0]/$_[5]);
	my $x2 = int($_[1]/$_[5]);
	if($_[3] > 1) { $x2 -= 5*$_[6]; }
	my $ystart = int($_[3]/2);
	for(my $i=0; $i<$_[3]; ++$i){
		my $y2 = $_[4]+$ystart-$i;
		print INFO6 "$LINE $_[2] 1 0 $x1 $y2 $x2 $y2\n";
	}
	if($printNumbers){
		$y1 = $_[4] + 7;
		if($x2 > $x1 && $x2 < $x1+18) { $x2 = $x1+18; }
		elsif($x1 >= $x2 && $x1 < $x2+18) { $x1 = $x2+18; }
		if($_[7] != -1){
			print INFO6 "$TEXT $black $tinyFont $x1 $y1 $_[7]\n";
		}
		if($_[8] != -1){
			print INFO6 "$TEXT $black $tinyFont $x2 $y1 $_[8]\n";
		}
	}
}

#***********************************
sub plot_line1
#***********************************
# variable: 0=x1 1=x2 2=color 3=width 4=y_shift 5=scale 6=direction 7=number1 8=number2 9=script
{
	my $y1;
	my $x1 = int($_[0]/$_[5]);
	my $x2 = int($_[1]/$_[5]);
	if($_[3] > 1) { $x2 -= 5*$_[6]; }
	my $ystart = int($_[3]/2);
	my $ref = $_[9];
	for(my $i=0; $i<$_[3]; ++$i){
		my $y2 = $_[4]+$ystart-$i;
		$$ref .= "$LINE $_[2] 1 0 $x1 $y2 $x2 $y2\n";
	}
	if($printNumbers){
		$y1 = $_[4] + 7;
		if($x2 > $x1 && $x2 < $x1+18) { $x2 = $x1+18; }
		elsif($x1 >= $x2 && $x1 < $x2+18) { $x1 = $x2+18; }
		if($_[7] != -1){
			$$ref .= "$TEXT $black $tinyFont $x1 $y1 $_[7]\n";
		}
		if($_[8] != -1){
			$$ref .= "$TEXT $black $tinyFont $x2 $y1 $_[8]\n";
		}
	}
}

#**************************
sub   pearson_correlation
#**************************
{
my $xr = shift;
my $yr = shift;
my $restrict = shift;

if(!$xr || !$yr){ return (0,0); }
my $sx = 0;
my $sxx = 0;
my $sxy = 0;
my $sy = 0;
my $syy = 0;
my $n1 = @$xr;
my $n = 0;
for(my $i=0; $i<$n1; ++$i){
	my $x = $xr->[$i];
	my $y = $yr->[$i];
	if($x==$MISSING || $y==$MISSING){ next; }
	if($restrict==1 && $x<=0 || $restrict==2 && ($x<=0 || $y<=0) || 
	 $restrict==3 && $x>=0 || $restrict==4 && ($x>=0 || $y>=0)){
		next;
	} 
	$sx += $x;
	$sxx += $x*$x;
	$sxy += $x*$y;
	$sy += $y;
	$syy += $y*$y;
	++$n;
}
if($n<3){ return(0,0); }
my $vx = $sxx-$sx*$sx/$n;
my $vy = $syy-$sy*$sy/$n;
my $vxy = $vx*$vy;
my ($r,$bb,$aa)=(0,0,$sy/$n);
if($vxy>0){
	$r = ($sxy-$sx*$sy/$n)/sqrt($vxy);
	$bb = $r*sqrt($vy/$vx);
	$aa -= $bb*$sx/$n;
}
return ($r,$n,$bb,$aa);
}

#**************************
sub   median
#**************************
{
my $xr = shift;
if(!$xr){ return 0; }
my $n = @$xr;
if(!$n){ return 0; }
my $median = 0;
my @sorted = sort {$a<=>$b} @$xr;
while($sorted[0]==$MISSING){
	shift(@sorted);
	$n--;
}
if(!$n){ return $MISSING; }
my $i = $n/2;
if($i > int($i)){
	$i = int($i);
	$median = $sorted[$i];
}else{
	$median = ($sorted[$i] + $sorted[$i-1])/2;
}
return $median;
}

#**************************
sub   average
#**************************
{
my $xr = shift;
if(!$xr){ return 0; }
my $n = @$xr;
if(!$n){ return 0; }
my $avr = 0;
foreach my $x (@$xr){
	if($x==$MISSING){ $n--; }
	else{
		$avr += $x;
	}
}
if($n>0){
	$avr /= $n;
}
return $avr;
}

#*****************
sub find_in_list
#*****************
{
my $start = shift;
my $ref = shift;
my $N = shift;
my $pos = shift;

my $n = $N;
my $j = -1;
my $i = int($n/2);
while(1){
	my $start1;
	if(defined($pos)){ $start1 = $ref->[$i]->[$pos]; }
	else{ $start1 = $ref->[$i]; }
	if($start1==$start){
		last;
	}
	elsif($start > $start1){
		if($n-$i<=1){ last; }
		$j = $i;
		$i = int(($n+$i)/2);
	}
	else{
		if($i-$j<=1){ --$i; last; }
		$n = $i;
		$i = int(($j+$n)/2);
	}
}
return $i;
}

#**********************************************
sub  normal_distribution
#**********************************************
{
my $x = shift;
my $a1=-1.26551223;
my $a2= 1.00002368;
my $a3= 0.37409196;
my $a4= 0.09678418;
my $a5=-0.18628806;
my $a6= 0.27886807;
my $a7=-1.13520398;
my $a8= 1.48851587;
my $a9=-0.82215223;
my $a10=0.17087277;
my $z = abs($x/sqrt(2));
my $t = 1.0 / (1.0 + 0.5 * $z);
my $y = $t*exp(-$z * $z + $a1 + $t * ($a2 + $t * ($a3 + $t * ($a4 + $t * ($a5 + $t *
     ($a6 + $t * ($a7 + $t * ($a8 + $t * ($a9 + $t * $a10)))))))));
if($x < 0.0){ $y = 2.0 - $y; }
$y = 1.0 - 0.5 * $y;
return $y;
}

#*****************************************
sub  parse_data
#*****************************************
{
my $data_ref = shift;
my $hashForm_ref = shift;

my @lines = split(/\r?\n\r?/,$$data_ref);
if(@lines==1){
	@lines = split(/\n?\r\n?/,$$data_ref);
}
my $iline = 0;
while($iline<@lines && $lines[$iline] =~ /^----/){
	my($name,$attribute,$filename) = ("","","");
	++$iline;
	while($iline<@lines && $lines[$iline]){
		my @items = split(/; /,$lines[$iline]);
		foreach my $item (@items){
			if($item !~ /=/){ next; }
			my($type,$value) = split(/=/,$item);
			if($type eq "name"){
				$value =~ s/"//g;
				$name = $value;
			}
			if($type eq "filename"){
				$value =~ s/"//g;
				$filename = $value;
			}
		}
		++$iline;
	}
	while($iline<@lines && !$lines[$iline]){ ++$iline; }
	while($iline<@lines && $lines[$iline] !~ /^----/){
		$attribute .= "$lines[$iline]\n";
		$iline++;
	}
	$attribute =~ s/\n+$//;
	if($name && $attribute){
		$name =~ s/[^[:ascii:]]//g;
		$attribute =~ s/[^[:ascii:]]//g;
		if($filename){
			$filename =~ s/[^[:ascii:]]//g;
			$hashForm_ref->{$name} = [$filename,$attribute];
		}else{
			$hashForm_ref->{$name} = $attribute;
		}
	}
}
return;
}

#************************************
sub  terminate_task
#************************************
{
my $pid = $hashInput{"process_id"};
my $continue;
my $runID = $hashInput{"runID"};
if($runID == $RUN_MATRIX_UPLOAD || $runID == $RUN_NORMALIZE){
	$continue = "continue";
}
my $pid1 = pid_retrieve($pid);
if(!$pid1 || $pid1 !~ /^\d+$/){
	error_message("Process $pid - $pid1 not found",$continue);
}
my $response = `ps h --ppid $pid1 -o pid`;
my $response2 = `ps $pid1`;
if($response2 =~ /exatlas/){
	`kill $pid`;
	foreach my $ch (split(/\n/,$response)){
		$ch =~ s/^\s+//;
		if($ch =~ /\s/){ $ch =~ s/\s\.+$//; }
		my $response1 = `ps $ch`;
		if($response1 =~ /exatlas|pairwise|anova_|norm_new/){ `kill $ch`; }
	}
	terminal_window("<h3>Your task was terminated!</h3>",$continue);
}
error_message("Process apparently finished",$continue);
}

#**************************************
sub  get_outputID
#**************************************
{
my $block = shift;

if(!$block){ $block=1; }
my @letters = ('1'..'9');
my $outputID=$letters[rand@letters];
push(@letters,'0');
for(my $i=0; $i<14; $i++){
	$outputID .= $letters[rand@letters];
}
for(my $id=$outputID; $id<$outputID+$block; ++$id){
	my @list = glob("$PATH_OUTPUT/$id.*");
	if(@list){
		unlink "$PATH_OUTPUT/$id.*";
	}
}
return $outputID;
}

#*****************************************
sub  get_array_lists
#*****************************************
{
my $ref = shift;

my $items="";
my $descriptions="";
foreach my $ref1 (@$ref){
	my($item,$descr,$orgID,$date) = @$ref1;
	if(!$items){ $items = "\"".$item."\""; }
	else{ $items .= ",\"".$item."\""; }
	if(!$descriptions){ $descriptions = "\"".$descr."\""; }
	else{ $descriptions .= ",\"".$descr."\""; }
}
return($items,$descriptions);
}

#**************************************
sub   print_web_page
#**************************************
{
my $fileID = shift;
open(INFO_TEMP,"<$PATH_OUTPUT/$fileID.txt");
while(my $line=<INFO_TEMP>){ print $line; }
close INFO_TEMP;
exit(0);
}

#**************************************
sub  error_message
#**************************************
{
my $message = shift;
my $comment = shift;

my ($header,@lines) = split(/\n/,$message);
open(OUTERR, ">>$PATH_INFO/error_log.txt");
print OUTERR "$date_record\t$loginname\t".join(" - ",$header,@lines)."\n";
close OUTERR;

my $text = "<H3>ERROR: $header</H3>\n";
foreach my $line (@lines){
	$text .= "$line<br>\n";
}
if($comment =~ /^\d+$/){
	file_append("ERROR: $message","$PATH_OUTPUT/$comment.txt");
	exit(0);
}elsif($comment eq "register"){
	print "<HTML><HEAD><TITLE>ExAtlas error message</TITLE>\n";
	print_header();
	print $text."<p>\n";
	print "&nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; <A HREF=../index.html><IMG SRC=../images/button_home.gif BORDER=0></A><p>";
}else{
	terminal_window($text,$comment);
}
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  terminal_window
#**************************************
{
my $message = shift;
my $comment = shift;

my $organismID = $hashInput{"organismID"};
print "<HTML><HEAD><TITLE>ExAtlas response</TITLE>\n";
print_header();
print $message."\n";
print "<p><HR NOSHADE><p>\n";
print "<i>Please report any problems to <a href=mailto:sharoval\@mail.nih.gov>webmaster</a><p>\n";
if($comment eq "continue"){
	print "<FORM NAME=exatlas ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
	print "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
	print "<INPUT NAME=\"action\" TYPE=hidden VALUE=\"continue\">\n";
	print "<INPUT NAME=\"continue\" TYPE=submit VALUE=\"   Continue   \">\n";
	if($organismID){
		print "<INPUT NAME=organismID TYPE=hidden VALUE=$organismID>\n";
	}
	print "</FORM><p>\n";
}elsif($comment eq "register"){
	print "&nbsp; &nbsp; <A HREF=../index.html><IMG SRC=../images/button_home.gif BORDER=0></A><p>";
}else{
	print "<INPUT TYPE=button VALUE=\"Close window\" style=width:160px; LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
}
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  menu_geneset_overlap
#**************************************
{
my $geneset_list = shift;
my $fileID = shift;
my $text;

if(!$geneset_list || ref($geneset_list) ne "ARRAY"){ return(""); }
my @FDR_list = (1,0.5,0.2,0.1,0.05,0.01,0.001,0.0001);
my @fold_enrichment = (0.0001,1,1.5,2,3,4,5);
my @minN = (1,2,3,4,5,7,10);
my ($FDR,$fold_enrichment,$minN) = (0.05,2,5);
$text .="<TR><TD><b>Use geneset:<TD><select name=\"file_geneset1\" style=\"width:250px;\" onChange=update_description();>\n";
for(my $i=0; $i<@$geneset_list; ++$i){ 
	$text .="<option value=\"$geneset_list->[$i]->[0]\"";
	if($geneset_list->[$i]->[0] =~ /^public-GO/){ $text .=" selected"; }
	$text .="> $geneset_list->[$i]->[0]\n";
}
$text .="</select><td><INPUT NAME=\"description_geneset1\" style=width:250px;><TD>\n";
$text .="<INPUT TYPE=button VALUE=\" Overlap analysis \" onClick=\"geneset_overlap('$fileID.txt');\">\n";
$text .="<TR><TD><b>Parameters:<TD><select name=FDR>\n";
for(my $i=0; $i<@FDR_list; ++$i){ 
	$text .="<option value=$FDR_list[$i]"; if($FDR==$FDR_list[$i]){ $text .=" selected"; } $text .="> $FDR_list[$i]\n";
}
$text .="</select> FDR threshold\n";
$text .="<TD><select name=fold_enrichment>\n";
for(my $i=0; $i<@fold_enrichment; ++$i){ 
	$text .="<option value=$fold_enrichment[$i]"; if($fold_enrichment==$fold_enrichment[$i]){ $text .=" selected"; } $text .="> $fold_enrichment[$i]\n";
}
$text .="</select> Fold enrichment threshold\n";
$text .="<TD><select name=minimum_genes>\n";
for(my $i=0; $i<@minN; ++$i){ 
	$text .="<option value=$minN[$i]"; if($minN==$minN[$i]){ $text .=" selected"; } $text .="> $minN[$i]\n";
}
$text .="</select> N genes (min)\n";
$text .="<TR><TD><TD COLSPAN=3><INPUT TYPE=CHECKBOX NAME=use_attribute> Use gene attributes (if available)\n";
return $text;
}

#**************************************
sub  get_header
#**************************************
{
my $onload = shift;

my $text = "";
$text .= "</HEAD><BODY BGCOLOR=white";
if($onload){
	$text .= " onLoad=\"$onload\"";
}
$text .= ">\n";
$text .= "<TABLE BGCOLOR=black>\n";
$text .= "<TR><TD WIDTH=168 background=../images/head3.gif><IMG SRC=../images/head1.jpg BORDER=0></TD>\n";
$text .= "<TD WIDTH=678 background=../images/head4.gif ALIGN=CENTER><IMG SRC=../images/head2.gif BORDER=0 ALT=\"ExAtlas: DNA matrix finder\" useMap=#mapHead2></TD></TR>\n";
$text .= "</TABLE>\n";
$text .= "<map NAME='mapHead2'>\n";
$text .= "<area shape=rect href=../index.html coords=0,0,677,80>\n";
$text .= "</map>\n";
return $text;
}

#**************************************
sub  print_header
#**************************************
{
my $onload = shift;
my $text = get_header($onload);
print $text;
return;
}

#***********************************
sub  clean_up
#***********************************
{
my $date_cleaned;
if(open(INFO,"<$PATH_INFO/cleaned_date.txt")){
	$date_cleaned=<INFO>;
	chomp $date_cleaned;
	close INFO;
}
if($date_cleaned==$date){ return; }
my @ls_info = glob("$PATH_INFO/guest*");
foreach my $filename (@ls_info){
	$filename =~ s/^.+guest/guest/;
	if($filename =~ /guest-config\.txt/){ next; }
	my $date_user = substr($filename,5,6);
	if($filename =~ /^guest/ && $date_user && $date > $date_user+1){
		unlink "$PATH_INFO/$filename";
	}
}
my @ls_data = glob("$PATH_DATA/guest*");
foreach my $filename (@ls_data){
	$filename =~ s/^.+guest/guest/;
	if($filename =~ /guest-config\.txt/){ next; }
	my $date_user = substr($filename,5,6);
	if($filename =~ /^guest/ && $date_user && $date > $date_user+1){
		unlink "$PATH_DATA/$filename";
	}
}
my $command = "ls -l --time-style=+\"%y%m%d\" $PATH_OUTPUT";
my @ls_output = `$command`;
foreach my $line (@ls_output){
	my @items = split(/\s+/,$line);
	my $n = @items;
	my $filename = $items[$n-1];
	my $filedate = $items[$n-2];
	if($filename =~ /^\d\d\d\d/ && $filedate < $date-1){
		unlink "$PATH_OUTPUT/$filename";
	}
}
file_append("$date","$PATH_INFO/cleaned_date.txt",1);
return;
}

#**************************************
sub  get_bounds
#**************************************
{
my $ref=shift;
if(!$ref || ref($ref) ne "ARRAY" || !@$ref){ return $MISSING; }
my $n = @$ref;
my ($min,$max,$aver)=(1.0E15,-1.0E15,0);
foreach my $x (@$ref){
	if($x==$MISSING){ $n--; next; }
	if($min>$x){ $min=$x; }
	elsif($max<$x){ $max=$x; }
	$aver += $x;
}
if($n){ $aver /= $n; }
return($min,$max,$aver);
}

#**********************************************/
sub  gammln
#**********************************************
{
my $xx = shift;
my @cof=(76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5);
my $y = $xx;
my $x = $xx;
my $tmp=$x+5.5;
$tmp -= ($x+0.5)*log($tmp);
my $ser=1.000000000190015;
for (my $j=0;$j<=5;$j++){
	$ser += $cof[$j]/++$y;
}
return -$tmp+log(2.5066282746310005*$ser/$x);
}

#***********************************************/
sub  gammq
#***********************************************/
{
my $aaa = shift;
my $x = shift;

my ($gamser,$gammcf,$gln)=(0,0,0);
if ($x < 0.0 || $aaa <= 0.0){ print "Invalid arguments in routine gammq\n"; }
if ($x < ($aaa+1.0)) {
	gser(\$gamser,$aaa,$x,\$gln);
	return(1.0-$gamser);
}
gcf(\$gammcf,$aaa,$x,\$gln);
return $gammcf;
}

#***********************************************/
sub  gser
#***********************************************/
{
my $gamser=shift;
my $aaa=shift;
my $x=shift;
my $gln=shift;

$$gln=gammln($aaa);
if ($x <= 0.0) {
	if ($x < 0.0){ print "x < 0 in gser\n"; }
	$$gamser=0.0;
	return;
}
my $ap=$aaa;
my $del=1.0/$aaa;
my $sum=$del;
for(my $n=1;$n<=$ITMAX;$n++) {
	++$ap;
	$del *= $x/$ap;
	$sum += $del;
	if(abs($del) < abs($sum)*$EPS) {
		$$gamser=$sum*exp(-$x+$aaa*log($x)-($$gln));
		return;
	}
}
print "a too large, ITMAX too small in gser\n";
return;
}

#***********************************************/
sub   gcf
#***********************************************/
{
my $gammcf=shift;
my $aaa=shift;
my $x=shift;
my $gln=shift;

$$gln=gammln($aaa);
my $bbb=$x+1.0-$aaa;
my $c=1.0/$FPMIN;
my $d=1.0/$bbb;
my $h=$d;
my $i;
for ($i=1;$i<=$ITMAX;$i++) {
	my $an = -$i*($i-$aaa);
	$bbb += 2.0;
	$d=$an*$d+$bbb;
	if (abs($d) < $FPMIN){ $d=$FPMIN; }
	$c=$bbb+$an/$c;
	if (abs($c) < $FPMIN){ $c=$FPMIN; }
	$d=1.0/$d;
	my $del=$d*$c;
	$h *= $del;
	if (abs($del-1.0) < $EPS){ last; }
}
if ($i > $ITMAX){ print "aaa too large, ITMAX too small in gcf\n"; }
$$gammcf=exp(-$x+$aaa*log($x)-($$gln))*$h;
return;
}

#***********************************************/
sub   format_probability
#***********************************************/
{
my $x=shift;
if($x>0.01){ $x=floor(10000*$x+0.5)/10000; }
elsif($x>0.0001){ $x=floor(1000000*$x+0.5)/1000000; }
elsif($x>1.0e-99){ $x = sprintf("%.2e",$x); }
else{ $x=0; }
return $x;
}

#**********************************************************
sub acos { atan2( sqrt(1 - $_[0] * $_[0]), $_[0] ) }
#**********************************************************

#**********************************************************
sub  exp_transform
#**********************************************************
{
my $fileID = shift;

my $fileID1 = get_outputID(1);
open(INFO_1,"<$PATH_OUTPUT/$fileID.txt");
open(OUT_1,">$PATH_OUTPUT/$fileID1.txt");
my $line =<INFO_1>;
print OUT_1 $line;
while(my $line =<INFO_1>){
	chop $line;
	my ($row,@data) = split(/\t/,$line);
	for(my $i; $i<@data; $i++){
		my $x = exp($data[$i]);
		if($x>1000){ $x = int($x-1); }
		elsif($x>10){ $x = int(100*($x-1))/100; }
		else{ $x = int(1000*($x-1))/1000; }
		$data[$i] = $x;
	}
	print OUT_1 "$row\t".join("\t",@data)."\n";
}
close INFO_1;
close OUT_1;
copy "$PATH_OUTPUT/$fileID1.txt", "$PATH_OUTPUT/$fileID.txt";
return;
}

#**********************************************************
sub  get_line_counts
#**********************************************************
{
my $filename = shift;
my $response = `wc -l $filename`;
$response =~ s/\D.+$//;
return $response;
}

#**********************************************************
sub  parse_logfile
#**********************************************************
{
my $filename = shift;
my $response;
if($filename && open(INFO,"<$filename")){
	while(my $line=<INFO>){
		if($line =~ /task_stopped|task completed|error/i){ $response .= $line; }
	}
}
return $response;
}

#**********************************************************
sub  admin_job
#**********************************************************
{
my $logFileID = get_outputID(1);

file_append("Task started","$PATH_OUTPUT/$logFileID.txt",1);
my $admin_param = $hashInput{"admin_param"};   #Use: 0 0 1
my @items = split(/\s+/,$admin_param);
my $pid = fork();
if($pid){
	$|++;
	$pid = pid_record($pid);
	print "<HTML><HEAD><TITLE>ExAtlas - processing</TITLE>\n";
	print_header();
	print "<h3>Running ExAtlas update</h3>\n";
	my $x = int(10000*rand());
	print "The status of your task can be checked here: <a href=$HOME_ADDRESS/output/$logFileID.txt target=_blank$x>Log file</a><p>\n";
	print "<FORM ACTION=$CGI_ADDRESS/exatlas.cgi METHOD=POST>\n";
	print "<INPUT NAME=terminate_task TYPE=submit VALUE=\"  Cancel the task  \">\n";
	print "<INPUT NAME=sessionID TYPE=hidden VALUE=\"$sessionID\">\n";
	print "<INPUT NAME=process_id TYPE=hidden VALUE=\"$pid\">\n";
	print "<INPUT NAME=action TYPE=hidden VALUE=interrupt_program>\n";
	print "</FORM><p><HR NOSHADE>\n";
	print "</BODY>\n";
	print "</HTML>\n";
	exit(0);
}
# Child process
if($items[0]>=1){
	my @command = ("$PATH_BIN/update_platforms.pl","array_platforms.txt","-data","../../exatlasInfo/data","-bin","$PATH_BIN","-log","$PATH_OUTPUT/$logFileID.txt");
	if($items[0]>9){ push(@command,"-ex"); }
	system(@command);
	my $response = parse_logfile("$PATH_OUTPUT/$logFileID.txt");
	if($response =~ /task_stopped/){ print "Task stopped!<br>$response<br>\n"; exit(0); }
	print "Finished1<br>\n";
}
if($items[1]>=1){
	system("$PATH_BIN/GEO_update.pl","GEO_series_summary.txt","-data","../../exatlasInfo/data","-log","$PATH_OUTPUT/$logFileID.txt");
	my $response = parse_logfile("$PATH_OUTPUT/$logFileID.txt");
	if($response =~ /task_stopped/){ print "Task stopped!<br>$response<br>\n"; exit(0); }
	print "Finished2<br>\n";
}
if($items[2]>=1){
	my @command = ("$PATH_BIN/update_exatlas.pl","-info","../../exatlasInfo/info","-data", "../../exatlasInfo/data","-log","$PATH_OUTPUT/$logFileID.txt","-skip","1,2");
	if($items[3]>0){ push(@command,"-update",$items[3]); }
	system(@command);
	my $response = parse_logfile("$PATH_OUTPUT/$logFileID.txt");
	if($response =~ /task_stopped/){ print "Task stopped!<br>$response<br>\n"; exit(0); }
	print "Finished3<br>\n";
}
exit(0);

#my %hashFiles=("soma","Atsumi","AAA","Kanamori","juho","Holzenspies","JJ","Junghyun","Yann","Tapponnier","arch","AroonC","tlee","Tyrone","cuic","ChangyiC","Thea","VatsveenT","hsaw","SawakiH","BW","WeidenbuschB","amit","amitveeru","rlai","rosalindpmlai","lili","LGreger","Anbu","PalanisamyA","devi","LeimarembiN","ruba","MajzoubR","FRM","RenaultF");
#foreach my $key (keys %hashFiles){
#	my $change = $hashFiles{$key};
#	my @ls_data = glob("$PATH_INFO/$key-*");
#	my $nn = @ls_data;
#	foreach my $filename (@ls_data){
#		my $filename1 = $filename;
#		$filename1 =~ s/\/$key-/\/$change-/;
#		`mv $filename $filename1`;
#	}
#}
#exit(0);
}

#**************************************
sub change_password_form
#**************************************
{
my $loginname = shift;
my $sessionID = shift;

print "<HEAD><title>Password Change Form</title>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\"><!--\n";
print "function check_onsubmit() {\n";
print "	if(!document.passwd_change.passwd.value){\n";
print "		alert(\"You need to enter existing password\");\n";
print "		return false;\n";
print "	}\n";
print "	var x = new String(document.passwd_change.new_passwd.value)\n";
print "	if(!x || x.length<5 || x.length>20){\n";
print "		alert(\"You need to enter new password (from 5 to 20 characters)\"+x);\n";
print "		return false;\n";
print "	}\n";
print "	if(x.search(/ |\\/|\\\|\\(|\\)|\"|\'|\`/) >= 0){\n";
print "		alert(\"Password should not include spaces, slashes, parentheses, quotes\");\n";
print "		return(false);\n";
print "	}\n";
print "	if(document.passwd_change.new_passwd.value != document.passwd_change.new_passwd1.value){\n";
print "		alert(\"New passwords do not match. Please retype.\");\n";
print "		return false;\n";
print "	}\n";
print "	document.passwd_change.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT></HEAD>\n";
print "<BODY BGCOLOR=WHITE>\n";
print "<H2>Change password</H2>\n";
print "<FORM NAME=passwd_change METHOD=POST ACTION=$CGI_ADDRESS/exatlas.cgi>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD>Current password:<TD><input  type=\"password\" name=\"passwd\" style=width:200px;>\n";
print "<TR><TD>New password:<TD><input type=\"password\" name=\"new_passwd\" value=\"\" style=width:200px;> Length 5-20 characters, no spaces, slashes, parentheses\n";
print "<TR><TD>Retype new password:<TD><input type=\"password\" name=\"new_passwd1\" value=\"\" style=width:200px;>\n";
print "<TR><TD><TD><INPUT TYPE=BUTTON VALUE = \"Change password\" onClick=\"check_onsubmit();\" style=width:200px;><p>\n";
print "<TR><TD><TD><INPUT TYPE=BUTTON VALUE = \"Cancel (close)\" onClick=\"window.close();\" style=width:200px;><p>\n";
print "</TABLE>\n";
print "<INPUT TYPE=hidden NAME=loginname VALUE=\"$loginname\">\n";
print "<INPUT TYPE=hidden NAME=sessionID VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=update_password>\n";
print "</FORM><p>\n";
exit(0);
}

#**************************************
sub update_password
#**************************************
{
my $loginname = shift;
my $new_passwd = $hashInput{"new_passwd"};

my @letters = ('A' .. 'Z', 'a' .. 'z', '0' .. '9', '_', '.');
open (INFO,"<$PATH_INFO/login.txt") or die $!;
open (OUT, ">$PATH_INFO/login_temp.txt") or die $!;
while(my $line = <INFO>){
	chomp $line;
	my($name,$passwd,$lastname,$firstname,$email)=split(/\t/, $line);
	if($name eq $loginname){
		my $salt = $letters[rand@letters] . $letters[rand@letters];
 		my $encryptedPsw = $new_passwd;
		for(my $i=0; $i<47; $i++){
			$encryptedPsw = crypt $encryptedPsw, $salt;
		}
		$line = join("\t",$name,$encryptedPsw,$lastname,$firstname,$email);
	}
	print OUT "$line\n";
}
close INFO;
close OUT;
copy "$PATH_INFO/login_temp.txt", "$PATH_INFO/login.txt";
unlink "$PATH_INFO/login_temp.txt";
terminal_window("<H2>Your password has been changed</H2>");
}

#**************************************
sub start_session
#**************************************
{
my $loginname = shift;
my $sessionID;
my @letters = ('A' .. 'Z', 'a' .. 'z', '0' .. '9', '_', '.');
my $n=@letters;
for(my $i=0; $i<12; $i++){
	$sessionID .= $letters[rand $n];
}

#Update session_list file
my $date0 = `date \'+%s\'`;
open(OUT, ">$PATH_INFO/temp_#session.txt");
open(INFO,"<$PATH_INFO/session_#list.txt");
while(my $line = <INFO>){
	chomp $line;
	if(!$line){ next; }
	my ($name,$session,$remote,$date1) = split(/\t/,$line);
	if($name ne $loginname){
		if(!$date1){
			print OUT "$name\t$session\t$remote\t$date0\n";
		}elsif($date0-$date1<43200){
			print OUT $line."\n";
		}
	}
}
close INFO;
my $salt = $letters[rand@letters] . $letters[rand@letters];
my $encryptedID = $sessionID;
for(my $i=0; $i<47; $i++){
	$encryptedID = crypt $encryptedID, $salt;
}
print OUT "$loginname\t$encryptedID\t$remoteID\t$date0\n";
close OUT;
copy "$PATH_INFO/temp_#session.txt", "$PATH_INFO/session_#list.txt";
unlink "$PATH_INFO/temp_#session.txt";
return $sessionID;
}

#**************************************
sub check_sessionID
#**************************************
{
my $sessionID = shift;

my $date0 = `date \'+%s\'`;
open(INFO,"<$PATH_INFO/session_#list.txt") or error_message("Session list file not found!");
while(my $line = <INFO>){
	chomp $line;
	my @items = split(/\t/,$line);
	if(length($items[1])<2){ next; }
	my $encryptedID = $sessionID;
	for(my $i=0; $i<47; $i++){
		$encryptedID = crypt $encryptedID, $items[1];
	}
	if($items[1] eq $encryptedID && $remoteID eq $items[2] && $date0-$items[3]<43200){
		close INFO;
		return($items[0]);
	}
}
close INFO;
error_message("Your session ID has expired. Please login again");
}

#**************************************
sub pid_record
#**************************************
{
my $pid = shift;
my $pid_replacement = get_outputID(1);
my $date0 = `date \'+%s\'`;
#print "D1 $pid $pid_replacement<br>\n";
file_append("$loginname\t$pid_replacement\t$pid\t$date0\n","$PATH_INFO/session_#list.txt");
return $pid_replacement;
}

#**************************************
sub pid_retrieve
#**************************************
{
my $pid_replacement = shift;
my $pid;

if($pid_replacement<1000){ error_message("process id lost"); }
my $date0 = `date \'+%s\'`;
open(OUT, ">$PATH_INFO/temp_#session.txt");
open(INFO,"<$PATH_INFO/session_#list.txt");
while(my $line = <INFO>){
	chomp $line;
	my ($name,$pid1,$pid2,$date1) = split(/\t/,$line);
	if($name eq $loginname && $pid1==$pid_replacement){
		$pid = $pid2;
	}else{
		print OUT $line."\n";
	}
}
close INFO;
close OUT;
if(!$pid){ error_message("record not found"); }
copy "$PATH_INFO/temp_#session.txt", "$PATH_INFO/session_#list.txt";
unlink "$PATH_INFO/temp_#session.txt";
#print "E1 $pid_replacement $pid<br>\n";
return $pid;
}

#**************************************
sub  check_form
#**************************************
{
my $key=shift;
my $value=shift;

#print "A1 key=$key $value<br>\n";
unless($key =~ /^[-_.a-zA-Z0-9]+$/ && length($key)<=80){ return 1; }
if($key =~ /sendmail/i){ return 1; }
if(!$value){ return 0; }
my $changed=0;
if($key eq "upload_filename"){
	if($hashInput{"file_type"} =~ /^(genelist|expression)$/){ return 0; }
	unless(ref($value) eq 'ARRAY'){ return 1; }
	if($value->[0] =~ /\\/){ $value->[0] =~ s/^.+\\//; }
	if($value->[0] =~ /\//){ $value->[0] =~ s/^.+\///; }
	if($value->[0] =~ /\s/){ $value->[0] =~ s/\s/_/g; }
	if($value->[0] =~ /\./){
		if($value->[0] =~ /\.xls/){ error_message("Excel files cannot be uploaded\nSave your file as a tab-delimited text"); }
		elsif($value->[0] =~ /\.(doc|ppt)/){ error_message("Only tab-delimited text files can be uploaded"); }
		$value->[0] =~ s/\.txt$/\@txt/; $value->[0] =~ s/\./_/g; $value->[0] =~ s/\@/./;
	}
	if(!$value->[0]){ error_message("Invalid file name"); }
	my $email1 = "sharoval\@mail.nih.gov";
	my $filename = $value->[0];
	unless($filename =~ /^[-.\w]+$/ && length($filename)<=150){
		print "File name should not include special characters or spaces (length<150). Rename file<br>\n";
		return 1;
	}
	if($value->[1] =~ /\{/){ $value->[1] =~ s/\{/(/g; $value->[1] =~ s/\}/)/g; }
	if($value->[1] =~ /\[/){ $value->[1] =~ s/\[/(/g; $value->[1] =~ s/\]/)/g; }
	if($value->[1] =~ /<[^>]+>/){ $value->[1] =~ s/<[^>]+>//g; }
	if($value->[1] =~ /['<\?\r\^\%\~\&\`\%\$\|]/){ $value->[1] =~ s/['<\?\r\^\%\~\&\`\%\$\|]//g; }
	if($value->[1] =~ /\@/){ $value->[1] =~ s/\@/.at./g; }
	if($value->[1] =~ /sendmail/i){ $value->[1] =~ s/sendmail/send mail/gi; }
 	unless($value->[1] =~ /^["\-\+\s>,.:;\!\#\~\=\*\(\)\/\\\w]+$/){
		my @items = split(/["\-\+\s>,.:;\!\#\~\=\*\(\)\/\\\w]+/,$value->[1]);
		while(@items && !$items[0]){ shift(@items); }
		print "<b>Unexpected string \`$items[0]\` in the file. Remove it before uploading.</b><br>\n";
		`echo "File(content) - $items[0]" | mailx -s "ExAtlas failed" $email1`;
		return 1;
	}
	$value->[1] =~ s/[\s_]+$//;
	for(my $i=0; $i<@unacceptable; $i++){ 
		if($value->[1] =~ /$unacceptable[$i]/){
			$value->[1] =~ s/$unacceptable[$i]//g;
		}
	}
	return 0;
}
if($key =~ /^(new_)*passwd1*$/){
	if($hashInput{"loginname"}=~/^guest/i){ return 0; }
	if(length($value)>20 || length($value)<3){ return 1; }
	return 0;
}
elsif($key =~ /^(last|first)name$/){
	if($value =~ /^\s/){ $value =~ s/^\s+//; $changed=1; }
	if($value =~ /\s$/){ $value =~ s/\s+$//; $changed=1; }
	if($value =~ /\s\s/){ $value =~ s/\s\s*/ /g; $changed=1; }
	if($changed){ $hashInput{$key}=$value; }
	#if($value =~ /\s/){ $value =~ s/\s+/_/g; }
	unless($value =~ /^[- \w\.]+$/){ print "Remove special characters from your $key<br>\n"; return 1; }
	return 0;
}
elsif($key eq "loginname"){
	unless($value =~ /^[+\w\@.\-]+$/ && length($value)<=30){ print "Login name should be less than 30 characters, do not use special characters and spaces<br>\n"; return 1; }
	return 0;
}
elsif($key eq "sessionID"){
	unless($value =~ /^[._a-zA-Z0-9]+$/ && length($value)<=20 && length($value) > 4){ return 1; }
	return 0;
}
elsif($key eq "selected_items"){
	if($value eq "copy-",$hashInput{"copy_file"}){ return 0; }
	if($value =~ /[ \+,;:=\@\!\#\(\)]/){ $value =~ s/[ \+,;:\=\@\!\#\(\)]/_/g; $changed=1; }
	unless($value =~ /^[-.\w]+$/){ return 1; }
	if($value =~ /\./){
		$value =~ s/\.txt$/\@txt/; 
		$value =~ s/\./_/g; $value =~ s/\@/./g; $changed=1;
	}
	return 0;
}
if($value =~ /[\s_]+$|^[\s_]+/){ $value =~ s/[\s_]+$|^[\s_]+//; $changed=1; }
for(my $i=0; $i<@unacceptable; $i++){
	if($value =~ /$unacceptable[$i]/){
		$value =~ s/$unacceptable[$i]//g; $changed=1;
	}
}
if($key eq "pasted_text"){
	if($value =~ /[ ]+\n/){ $value =~ s/[ \t]+\n/\n/g; $changed=1; }
	if($value =~ /\{/){ $value =~ s/\{/(/g; $value =~ s/\}/)/g; $changed=1; }
	if($value =~ /\[/){ $value =~ s/\[/(/g; $value =~ s/\]/)/g; $changed=1; }
	if($value =~ /<[^>]+>/){ $value =~ s/<[^>]+>//g; $changed=1; }
	if($value =~ /['<\?\r\^\%\~\&\`\$\|]/){ $value =~ s/['<\?\r\^\%\~\&\`\$\|]//g; $changed=1; }
	if($value =~ /\@/){ $value =~ s/\@/.at./g; $changed=1; }
	if($value =~ /sendmail/i){ $value =~ s/sendmail/send mail/gi; $changed=1; }
	unless($value =~ /^[-\s">+,.:;!#*=~\(\)\/\\\w]+$/){ return 1; }
	if($changed){ $hashInput{$key}=$value; }
	return 0;
}
if($value =~ /sendmail|\/etc|foo;|\/passwd/i){ return 1; }
if($value =~ /^\d+\.\s+\w+/){ $value =~ s/^\d+\.\s+//; $changed=1; }
if($value =~ /\{/){ $value =~ s/\{/(/g; $value =~ s/\}/)/g; $changed=1; }
if($value =~ /\[/){ $value =~ s/\[/(/g; $value =~ s/\]/)/g; $changed=1; }
if($value =~ /[\"\'\n\r]/){ $value =~ s/[\"\'\n\r]//g; $changed=1; }
if($value =~ /<[^>]+>/){ $value =~ s/<[^>]+>//g; $changed=1; }
if($value =~ /[<>'"\?\%\$\`]/){ $value =~ s/[<>'"\?\%\$\`]//g; $changed=1; }
if($value =~ /^\s+/){ $value =~ s/^\s+//; $changed=1; }
if($value =~ /\s+$/){ $value =~ s/\s+$//; $changed=1; }
if(length($value) > 250){
	$value = substr($value,0,250);
	$changed=1;
}
unless($value =~ /^[ \-\+=_.,;:@\!\%\#\*\w\(\)\/\\]+$/){ return 1; }
#remove special characters from beginning and end
if($value =~ /^[ _=.,;:\@\!\#\(\)\/\\]+/){ $value =~ s/^[ +_=.,;:\@\!\#\*\(\)\/\\]+//; $changed=1; }
if($value =~ /[ _=,;:\@\#\(\/\\]+$/){ $value =~ s/[ _=,;:\@\#\*\(\/\\]+$//; $changed=1; }

if($key eq "command" && $value=~/^print_table,\d+,\d+$|^save_geneset,\d+$/){ return 0; }
if($value =~/ / || $key =~ /description|_header|_term|[_\-]title|^item_name|^terminate_t|^reset_pa|^continue$|^display$|^search_term$|^geneset_name|^matrix_name|^admin_param/){
	if($value =~ /^-+ *new file *-+$/i){ return 0; }
	if($key !~/description|_header|_term|[_\-]title|^item_name|^terminate_t|^reset_pa|^continue$|^display$|^search_term$|^geneset_name|^matrix_name|^admin_param/){ return 1; }
	if($value =~ /^-+ *new file *-+$/i && $key=~/filename|^file_|upload_geneset|_file$|select_metaanalysis|^copysamples|^copy_to_geneset|^matrix_metaanalysis/){
		return 0;
	}
	if($value =~ /\@/){ $value =~ s/\@/.at./g; $changed=1; }
	if($value =~ /;/){ $value =~ s/;/,/g; $changed=1; }
}elsif($key=~/email/){
	my $passed=1;
	my @items=split(/\@/,$value); if (@items!=2){ $passed=0; }
	foreach my $word (split(/\./,$items[0].".".$items[1])){
		if($word !~ /^[-\w]/ || $word=~/^-|-$/){ $passed=0; }
	}
	if(!$passed){ error_message("Email not valid"); }
}elsif($key=~/filename|^rename|^new_file|^source_matrix|^file_|upload_geneset|_file$|select_metaanalysis|^copysamples|^copy_to_geneset|^matrix_metaanalysis/){
	if($value =~ /^TEMPORARY!-\d+\.txt$/){ return 0; }
	if($value =~ /[ +,;:=@!#\(\)]/){ $value =~ s/[ +,;:=@!#\(\)]/_/g; $changed=1; }
	if($value =~ /\.+$/){ $value =~ s/\.+$//; $changed=1; }
	unless($value =~ /^[\-.\w]+$/){ return 1; }
	if($value =~ /\./){
		$value =~ s/\.txt$/\@txt/; 
		$value =~ s/\./_/g; $value =~ s/\@/./g; $changed=1;
	}
}else{
	if(length($value) > 100){
		$value = substr($value,0,100);
		$changed=1;
	}
	unless($value =~ /^[-\+.\w]+$/){ return 1; }
	if($value =~ /\./){
		my $count=0; $count++ while($value =~ /\./g);
		if($count!=1){ return 1; }
		if($value !~ /e/i){
			if($value !~/^[-+]*[\.\d]+$/ || $value eq "."){ return 1; }
		}elsif($value !~ /^[\.\d][\.\d]+e[-\+]*\d+$/i){
			return 1;
		}		
	}
}
#print "A3 $key $value<br>\n";
if($changed){ $hashInput{$key}=$value; }
return 0;
}


