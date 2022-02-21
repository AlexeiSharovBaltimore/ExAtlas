#!/usr/local/bin/perl

use strict;
use Net::FTP;
my $usage = "$0 output logFile\n";
my $arg=0;
my $output_file = $ARGV[$arg++] or die $usage;
my $logfile;
my $start = 0;
my $PATH_DATA=".";
my $RNAseq_option=0;
while(my $option = $ARGV[$arg++]){
	if(uc($option) eq "-START"){ $start = $ARGV[$arg++] or die "start number missing"; }
	elsif(uc($option) eq "-DATA"){ $PATH_DATA = $ARGV[$arg++] or die "path data missing"; }
	elsif(uc($option) eq "-LOG"){ $logfile = $ARGV[$arg++] or die "logfile missing"; }
	elsif(uc($option) eq "-RNASEQ"){ $RNAseq_option=1; }
	else{ die "ERROR: Wrong option $option\n"; }
}

my @exclude = ("GSE39655","GSE48099","GSE48190","GSE48377","GSE48386","GSE48482","GSE48483","GSE48484","GSE49045","GSE49082","GSE49123","GSE49215","GSE49276","GSE49280","GSE49317","GSE49357","GSE49396","GSE49452","GSE49602","GSE49614","GSE49615","GSE50024","GSE50025","GSE5013","GSE50252","GSE50253","GSE50441","GSE50495","GSE50496","GSE50558","GSE50666","GSE50667","GSE50668","GSE50727","GSE5082","GSE50973","GSE51158","GSE51265","GSE51538","GSE51589","GSE51640","GSE51641","GSE51711","GSE5173","GSE51734","GSE51973","GSE51976","GSE52107","GSE52147","GSE52155","GSE52221","GSE52442","GSE53021","GSE53096","GSE53141","GSE53261","GSE53425","GSE53426","GSE53436","GSE53445","GSE5347","GSE53488","GSE53577","GSE53626","GSE54445","GSE54446","GSE54509","GSE54948","GSE5511","GSE55134","GSE55175","GSE55230","GSE55232","GSE55586","GSE55638","GSE55640","GSE56043","GSE56316","GSE56518","GSE56685","GSE56781","GSE57100","GSE57117","GSE57118","GSE57229","GSE57277","GSE57612","GSE57981","GSE58356","GSE58951","GSE59032","GSE59150","GSE59244","GSE59546","GSE59678","GSE60720","GSE61122","GSE61129","GSE61660",
	"GSE62875","GSE64509","GSE66881");
my %hashExclude;
foreach my $x (@exclude){ $hashExclude{$x}=1; }

my $host = "ftp.ncbi.nlm.nih.gov";
my $user = "anonymous";
my $password = "user\@comcast.net";
my $f = Net::FTP->new($host, Port => 21, Passive => 0) or die "Can't open $host\n";
$f->login($user, $password) or die "Can't log $user in\n";
$f->binary();

if(!$start && open(INFO, "$PATH_DATA/$output_file")){
	my $lastID=0;
	while(my $line=<INFO>){
		if($line =~ /^>/){
			my ($ID,$junk) = split(/\t/,$line);
			$ID =~ s/^>GSE//;
			if($lastID < $ID){ $lastID=$ID; }
		}
	}
	close INFO;
	$start = $lastID+1;
}
open(OUT, ">>$PATH_DATA/$output_file") or die $!;
if($logfile){
	open(OUTLOG, ">>$logfile") or die $!;
}
open(FTPLOG, ">>$PATH_DATA/FTP_log.txt") or die $!;
if($start){
	if($logfile){ print OUTLOG "Starting ID = GSE$start\n"; }
	else{ print "Starting ID = GSE$start\n"; }
}
my @organisms = ([9606,15000,"Homo sapiens","human"],
[10090,12000,"Mus musculus","mouse","mice"],
[10116,12000,"Rattus norvegigus","rat"],
[9544,12000,"Macaca mulatta","rhesus monkey"],
[9541,12000,"Macaca fascicularis","macaque","monkey"],
[9598,12000,"Pan troglodytes","chimpanzee","monkey","primate"],
[9615,12000,"Canis familiaris","dog"],
[9940,12000,"Ovis aries","sheep"],
[9823,12000,"Sus scrofa","pig"],
[9913,12000,"Bos taurus","cow"],
[9796,12000,"Equus caballus","horse"],
[9986,12000,"Oryctolagus cuniculus","rabbit"],
[9031,12000,"Gallus gallus","chiken"],
[9103,10000,"Meleagris gallopavo","turkey"],
[8355,10000,"Xenopus laevis","xenopus"],
[7955,10000,"Danio rerio","zebrafish"],
[8022,8000,"Oncorhynchus mykiss","rainbow trout","trout"],
[8030,8000,"Salmo salar","salmon","fish"],
[7227,10000,"Drosophila melanogaster","fruit fly","insect","fly"],
[6239,8000,"Caenorhabditis elegans","nematode"],
[3702,8000,"Arabidopsis thaliana","thale cress"],
[4530,10000,"Oryza sativa","rice"],
[3847,10000,"Glycine max","soybean"],
[4081,10000,"Lycopersicon esculentum","tomato"],
[4577,10000,"Zea mays","maize"],
[5476,4000,"Candida albicans","yeast","fungi"],
[4932,4000,"Saccharomyces cerevisiae","yeast","fungi"],
[90371,3000,"Salmonella typhimurium","salmonella","protista"],
[562,3000,"Escherichia coli","bacterium"],
[300852,2000,"Thermus thermophilus","bacterium"],
[1773,2000,"Mycobacterium tuberculosis","bacterium"],
[1280,2000,"Staphylococcus aureus","bacterium"],
[1423,2000,"Bacillus subtilis","bacterium"]);

my %hashPlatforms;
if($RNAseq_option){
	print OUTLOG "Reading RNAseq platforms\n";
	if(open(INFO,"$PATH_DATA/sequencing_platforms.txt")){
		while(my $line=<INFO>){
			chop $line;
			my ($platformID,$title,$species,$nProbes) = split(/\t/,$line);
			$hashPlatforms{$platformID} = [$title,$title,$species,$nProbes];
		}
		close INFO;
	}
}else{
	print OUTLOG "Reading array platforms\n";
	if(open(INFO,"$PATH_DATA/array_platforms.txt")){
		while(my $line=<INFO>){
			chop $line;
			my ($platformID,$title,$technology,$species,$nProbes,$taxid,$type) = split(/\t/,$line);
			$hashPlatforms{$platformID} = [$title,$technology,$species,$nProbes];
		}
		close INFO;
	}
}
my %hashOrganism;
foreach my $ref (@organisms){
	my @items = @$ref;
	my $key = shift(@items);
	$hashOrganism{$key} = \@items;
}

#Download data series
my $dir = "geo/series";
$f->cwd($dir) or die "Can't cwd to $dir\n";
my $count=0;
my $MAX_COUNT=15000000;
my @parts = $f->ls;
my @partNo;
foreach my $part (@parts){
	my $x = $part;
	$x =~ s/nnn/999/;
	$x =~ s/^GSE//;
	push(@partNo,$x);
}
my @sorted = sort {$partNo[$a]<=>$partNo[$b]} 0..(@partNo-1);
foreach my $ipart (@sorted){
	my $part = $parts[$ipart];
	my $x = $partNo[$ipart];
	if($x < $start){ next; }
	if(!$f->cwd("/geo/series/$part")){ print "Can't cwd to $part\n"; next; }
	my @series = $f->ls;
	foreach my $series (sort @series){
		my $partNo = $series;
		$partNo =~ s/^GSE//;
		if($partNo < $start){ next; }
		if($hashExclude{$series}){ next; }  # not expression and too big!
		$count++;
		#if($count > $MAX_COUNT){ last; }
		if(!$f->cwd("/geo/series/$part/$series")){
			print "Can't cwd to $series\n";
			next;
		}
		my @seriesContent = $f->ls;
		my $found=0;
		foreach my $x (@seriesContent){ if($x eq "matrix"){ $found=1; }}
		if(!$found){
			if($logfile){ print OUTLOG "ERROR: No matrix dir in $series\n"; }
			else{ print "ERROR: No matrix dir in $series\n"; }
			next;
		}
		if(!$f->cwd("matrix")){ print "ERROR: Can't cwd to $series\n"; next; }
		my @matrix = $f->ls;
		if(!@matrix){
			if($logfile){ print OUTLOG "ERROR: Matrix in $series not found\n"; }
			else{ print "ERROR: Matrix in $series not found\n"; }
			next;
		}
		foreach my $matrix_file (sort @matrix){
			my %hashSeries=();
			my $file_size = $f->size($matrix_file);
			if($file_size > 500000000){
				if($logfile){ print OUTLOG "$matrix_file too big $file_size > 500M\n"; }
				else{ print "$matrix_file too big $file_size > 500M\n"; }
				next;
			}
			if(!$f->get($matrix_file,"$PATH_DATA/$matrix_file")){
				if($logfile){ print OUTLOG "Can't get $matrix_file\n"; }
				else{ print "Can't get $matrix_file\n"; }
				next;
			}
			my $count=0;
			while(!$f && $count < 3){
				$f = Net::FTP->new($host);
				$f->login($user, $password);
				$f->binary();
				$f->cwd("/geo/series/$part/$series/matrix");
				$f->get($matrix_file,"$PATH_DATA/$matrix_file");
				$count++;
			}
			if(!$f){
				if($logfile){ print OUTLOG "Error reading $series\ntask_stopped\n"; }
				else{ print "Error reading $series\ntask_stopped\n"; }
				exit(0);
			}
			if(!file_exist("$PATH_DATA/$matrix_file")){
				if($logfile){ print OUTLOG "Error: $series missing\ntask_stopped\n"; }
				else{ print "Error: $series missing\ntask_stopped\n"; }
				exit(0);
			}
			parse_matrix_file("$PATH_DATA/$matrix_file", \%hashSeries);

			# FILTERING DATA
			my $submission_date = $hashSeries{"series_submission_date"};
			my ($month,$day,$year) = split(/\s+/,$submission_date);
			my $platform = $hashSeries{"series_platform_id"};
			if(!$platform){ if($logfile){ print OUTLOG "$year\t$series\tno_platform\n"; } next; }
			my $series_type = $hashSeries{"series_type"};
			if($series_type !~ /expression profiling|RNA profiling by array/i || $series_type =~ /sequencing/){
				if($logfile){ print OUTLOG "$year\t$series\t$platform\tType = $series_type\n"; } 
				next; 
			}
			my $taxid = $hashSeries{"series_sample_taxid"};
			if(!$taxid){
				if($logfile){ print OUTLOG "$year\t$series\t$platform\tno taxid\n"; }
				next;
			}
			my $ref = $hashSeries{"sample_title"};
			if(!$ref || ref($ref) ne 'ARRAY'){
				if($logfile){ print OUTLOG "$year\t$series\t$platform\tNo sample info\n"; } 
				next; 
			}
			$platform = uc($platform);
			my $nrows = $hashSeries{"sample_data_row_count"};
			if($nrows && ref($nrows) eq 'ARRAY'){ $nrows = $nrows->[0]; }
			if(!$nrows){
				if($logfile){ print OUTLOG "$year\t$series\t$platform\tNo rows\n"; }
				print FTPLOG "$year\t$series\t$platform\tNo rows\n"; 
				next;
			}
			my $refOrg;
			foreach my $id (split(/;/,$taxid)){
				if($hashOrganism{$taxid}){
					$refOrg = $hashOrganism{$taxid}; last;
				}
			}
			if(!$refOrg){ if($logfile){ print OUTLOG "$year\t$series\t$platform\tother_taxid = $taxid\n"; } next; }
			my $ref = $hashSeries{"sample_type"};
			if(!$ref || ref($ref) ne 'ARRAY' || $ref->[0] !~ /RNA/i){
				if($logfile){ print OUTLOG "$year\t$series\t$platform\tsample not RNA but $ref->[0]\n"; }
				print FTPLOG "$year\t$series\t$platform\tsample not RNA but $ref->[0]\n"; 
				next;
			}
			my $ref = $hashSeries{"sample_geo_accession"};
			if(!$ref || ref($ref) ne 'ARRAY' || @$ref==0){
				if($logfile){ print OUTLOG "$year\t$series\t$platform\tNo sample information\n"; }
				print FTPLOG "$year\t$series\t$platform\tNo sample information\n";
				next;
			}
			my @samples = @$ref;
			my $Nsamples = @samples;
			if(!$Nsamples){ 
				if($logfile){ print OUTLOG "$year\t$series\t$platform\tNo samples\n"; }
				print FTPLOG "$year\t$series\t$platform\tNo samples\n";
				next;
			}

			#PROCESSING DATA
			#print "$count - $series\n";
			if($logfile){ print OUTLOG "$year\t$series\t$platform\t$hashPlatforms{$platform}->[0]\tOK\n"; }
			else{ print "$year\t$series\t$platform\t$hashPlatforms{$platform}->[0]\tOK\n"; }
			my $key = "series_title";
			print OUT ">$series\t$platform\t$hashSeries{$key}\n";
			my $ref = $hashPlatforms{$platform};
			if($ref){ print OUT "platform: $ref->[0]\n"; }
			my $keywords = $hashSeries{"keywords"};
			print OUT "keywords:$keywords\n";
			print OUT "taxid:$taxid\n";
			$key = "series_contributor";
			my $ref = $hashSeries{$key};
			if($ref && ref($ref) eq 'ARRAY'){
				print OUT "contributor:$ref->[0]";
				my $nn = @$ref;
				if($nn>1){ print OUT ",$ref->[$nn-1]"; }
				print OUT "\n";
			}
			my $inst = $hashSeries{"series_contact_institute"};
			if($inst){ print OUT "institute: $inst\n"; }
			print OUT "samples:$Nsamples\n";
			if($year){ print OUT "year:$year\n"; }
			for(my $i=0; $i<$Nsamples; ++$i){
				#PROCESSING SAMPLES
				$key = "sample_channel_count";
				my $nChannels = $hashSeries{$key}->[$i];
				$samples[$i] = uc($samples[$i]);
				my ($startCh,$endCh)=(1,$nChannels);
				if($nChannels==2){
					$key = "sample_source_name_ch2";
					if($hashSeries{$key}->[0] =~ /reference/){
						$endCh = 1;
					}
					$key = "sample_source_name_ch1";
					if($hashSeries{$key}->[0] =~ /reference/){
						$startCh = 2;
					}
				}
				my $ref = $hashSeries{"sample_title"};
				for(my $j=$startCh; $j<=$endCh; ++$j){
					print OUT "$samples[$i]\t$ref->[$i]";
					if($startCh < $endCh){
						$key = "sample_source_name_ch$j";
						if($hashSeries{$key}){ print OUT "\tsource:$hashSeries{$key}->[$i]"; }
					}
					my $ref1 = $hashSeries{"sample_taxid_ch$j"};
					if(!$ref1 && ref($ref) eq 'ARRAY'){
						my $taxid1 = $ref1->[$i];
						if($taxid1 && $taxid1 ne $taxid){ print OUT "\ttaxid:$taxid1"; }
					}
					$ref1 = $hashSeries{"sample_characteristics_ch$j,cell_line"};
					if(!$ref1 && ref($ref) eq 'ARRAY'){
						print OUT "\tcell_line$ref->[$i]";
					}
					$ref1 = $hashSeries{"sample_characteristics_ch$j,pooled_sample"};
					if(!$ref1 && ref($ref) eq 'ARRAY'){
						print OUT "\tpooled:$ref->[$i]";
					}
					print OUT "\n";
				}
			}
		}
	}
}
close OUT;
close FTPLOG;
if($logfile){ close OUTLOG; }
print $f->message.".<br>\n";
$f->quit();
exit(0);

#******************************
sub  parse_matrix_file
#******************************
{
my $matrix_file = shift;
my $hashSeries = shift;

my $compressed = $matrix_file;
`gzip -d $matrix_file`;
my @contributors;
my $uncompressed = $matrix_file;
$uncompressed =~ s/\.gz$//;
open(INFO, $uncompressed) or die "File $uncompressed not found";
while(my $line=<INFO>){
	chop $line;
	if(!$line){ next; }
	if($line =~ /^\!series_matrix_table_begin/i){ last; }
	$line =~ s/^\!//;
	$line =~ s/\"//g;
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
		for(my $i=0; $i<@items; ++$i){
			if($items[$i] =~ /\: /){
				$items[$i] =~ s/[\[\]\(\)\{\}\@\#\$\^\&\*\%\!\?\+\.]//g;
				my($key1,$value) = split(/\: /,$items[$i]);
				$items[$i] = $value;
			}
		}
		$hashSeries->{lc($key)} = \@items;
	}
}
my $platid = $hashSeries->{"series_platform_id"};
my $platid1 = $hashSeries->{"sample_platform_id"};
if($platid1){ $platid1 = $platid1->[0]; }
my $taxid = $hashSeries->{"series_sample_taxid"};
my $taxid1 = $hashSeries->{"sample_taxid_ch1"};
if($taxid1){ $taxid1 = $taxid1->[0]; }

if($platid =~ /,/ && $platid1){
	$hashSeries->{"series_platform_id"} = $platid1;
}
if($taxid =~ /,/ && $taxid1){
	$hashSeries->{"series_sample_taxid"} = $taxid1;
}
$hashSeries->{"series_contributor"} = \@contributors;
close INFO;
`rm $uncompressed`;
return;
}

#***********************************
sub file_exist
#***********************************
{
if(open (INFO_TEMP,"<$_[0]")){ close INFO_TEMP; return 1; }
return 0;
}


