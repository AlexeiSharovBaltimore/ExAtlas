#!/usr/local/bin/perl

#bin\update_platforms.pl array_platforms.txt -data data -bin bin
use strict;
my $usage = "$0 platformFile logFile\n";
my $arg=0;
my $platforms_file = $ARGV[$arg++] or die $usage;
my $start = 0;
my $logfile;
my $PATH_DATA=".";
my $PATH_BIN=".";
my $check_existing=0;
while(my $option = $ARGV[$arg++]){
	if(uc($option) eq "-START"){ $start = $ARGV[$arg++] or die "start number missing"; }
	elsif(uc($option) eq "-DATA"){ $PATH_DATA = $ARGV[$arg++] or die "path data missing"; }
	elsif(uc($option) eq "-BIN"){ $PATH_BIN = $ARGV[$arg++] or die "path bin missing"; }
	elsif(uc($option) eq "-LOG"){ $logfile = $ARGV[$arg++] or die "logfile missing"; }
	elsif(uc($option) eq "-EX"){ $check_existing=1; }
	else{ die "ERROR: Wrong option $option\n"; }
}

my @organisms = ([9606,9000,"Homo sapiens","Human"],
[10090,9000,"Mus musculus","Mouse"],
[10116,9000,"Rattus norvegicus","Rat"],
[9544,9000,"Macaca mulatta","Rhesus monkey"],
[9541,9000,"Macaca fascicularis","Macaque"],
[9598,9000,"Pan troglodytes","Chimpanzee"],
[9615,9000,"Canis familiaris","Dog"],
[9940,9000,"Ovis aries","Sheep"],
[9823,9000,"Sus scrofa","Pig"],
[9913,9000,"Bos taurus","Cow"],
[9796,9000,"Equus caballus","Horse"],
[9986,9000,"Oryctolagus cuniculus","Rabbit"],
[9031,9000,"Gallus gallus","Chicken"],
[9103,8000,"Meleagris gallopavo","Turkey"],
[8355,8000,"Xenopus laevis","Xenopus frog"],
[7955,8000,"Danio rerio","Zebrafish"],
[8022,8000,"Oncorhynchus mykiss","Rainbow trout"],
[8030,8000,"Salmo salar","Salmon"],
[7227,8000,"Drosophila melanogaster","Fruit fly"],
[6239,7000,"Caenorhabditis elegans","Nematode"],
[3702,7000,"Arabidopsis thaliana","Thale cress"],
[4530,7000,"Oryza sativa","Rice"],
[3847,7000,"Glycine max","Soybean"],
[4081,7000,"Lycopersicon esculentum","Tomato"],
[4577,7000,"Zea mays","Maize"],
[4932,4000,"Saccharomyces cerevisiae","Yeast"],
[562,4000,"Escherichia coli","Bacterium"],
[300852,2000,"Thermus thermophilus","Bacterium"],
[1280,2000,"Staphylococcus aureus", "Bacterium"]);
my %hashOrganismNum;
my %hashOrganismID;
for(my $i=0; $i<@organisms; $i++){
	my $ref = $organisms[$i];
	$hashOrganismNum{$ref->[2]}=$i+1;
	$hashOrganismID{$ref->[0]}=$i+1;
}

if($logfile){
	open(OUT, ">>$logfile");
	if($check_existing){
		check_existing_platforms();
	}
}

my @output;
my $max_id = 0;
if($start){
	$max_id=$start;
}elsif(open(INFO, "$PATH_DATA/$platforms_file")){
	my $line = <INFO>;
	while(my $line = <INFO>){
		chop $line;
		my($ID,$title,$technology,$taxonomy,$rows,$taxid,$type) = split(/\t/,$line);
		my($found)=0;
		foreach my $id (split(/;/,$taxid)){
			if($hashOrganismID{$id}){ $found=1; last; }
		}
		if(!$found){ print OUT "ERR upfate_platforms.pl File $ID: organism $taxid not found\n"; next; }
		my $num = $ID;
		$num =~ s/^GPL//;
		if($max_id < $num){ $max_id=$num; }
	}
	close INFO;
}
print OUT "\nStarting from GPL$max_id\n\n";
print "Platforms starting from GPL$max_id<br>\n";
open(OUT1, ">>$PATH_DATA/$platforms_file");
my $missing=0;
my $lastProcessed;
use LWP::Simple;
for(my $idNum = $max_id+1; $missing < 5; $idNum++){
	my $platform = "GPL".$idNum;
	my $len = length($platform);
	if($len==8){ $platform = "GPL".$idNum; }
	elsif($len==7){ $platform = "GPL0".$idNum; }
	elsif($len==6){ $platform = "GPL00".$idNum; }
	elsif($len==5){ $platform = "GPL000".$idNum; }
	elsif($len==4){ $platform = "GPL0000".$idNum; }
	my $content = get("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$platform") or error_message("Unable to get page");
	if($content =~ /Could not find a public or private accession|Accession "$platform" was deleted by the GEO staff/){
		if($content =~ /Could not find a public or private accession/){
			$missing++;
		}
		next;
	}
	$lastProcessed = $platform;
	$missing=0;
	my @lines = split(/\n/,$content);
	my ($title,$technology,$taxonomy,$rows,$taxid,$type,$manufacturer,$description,$Nsamples);
	my $iline = 200;
	while($iline<@lines){
		if($lines[$iline++] =~ /^\<\/tr\>\<\/table\>\<\/form\>/i){ last; }
	}
	my $found = 0;
	my $min_rows = 0;
	while($iline<@lines){
		my $line = $lines[$iline++];
		if($line =~ /^<tr valign="top"><td nowrap>Title/i){
			$title = $lines[$iline++];
			$title =~ s/<[^>]+>//g;
			$title =~ s/^\[[^\]]*\]\s*//;
		}
		elsif($line =~ /^<tr valign="top"><td nowrap>Technology type/i){
			$technology = $lines[$iline++];
			$technology =~ s/<[^>]+>//g;
		}
		elsif($line =~ /^<tr valign="top"><td nowrap>Manufacturer/i){
			$manufacturer = $lines[$iline++];
			$manufacturer =~ s/<[^>]+>//g;
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
				my $iii = $hashOrganismNum{$species[$i]};
				if($iii){
					$found=1;
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
	if(!$found){ next; }
	if($technology && $technology !~ /^in situ oligonucleotide|^oligonucleotide beads/ ||
	 !$taxonomy ||
	 $title =~ /CGH|SNP|HELP|ROMA|CNV|HiSeq/ || $title =~ /promoter|tiling|chip-on-chip|chr[123456789XY _]/i ||
	 $description =~ /CGH|SNP|HELP|CNV|MIP/ || $description =~ /promoter|tiling|chip-on-chip|chr[123456789XY _]/i){
		#print "SKIPPED:\t$platform\t$rows\t$type\t$taxid\t$technology\n";
		next;
	}
	$type = "mRNA";
	if($title =~ /miR|microR|piRNA|PIWI|piwi/){ $type = "miRNA"; }
	elsif($rows > 200000 || $rows<$min_rows){
		print OUT "SKIPPED:\t$platform\t$rows\t$type\t$taxid\t$technology\n";
		next;
	}
	my @items = ($platform,$title,$technology,$taxonomy,$rows,$taxid,$type,$Nsamples,$manufacturer);
	push(@output,\@items);
	print OUT1 join("\t",@items)."\n";
	my $file1 = "public-$platform"."_annot.txt";
	if($PATH_DATA && $PATH_BIN){
		my $taxid1 = $taxid;
		$taxid1 =~ s/;/,/g;
		if($type eq "mRNA"){
			my $command = "$PATH_BIN/download_platform.pl $platform $file1 $taxid1 -data $PATH_DATA";
			#my $command = "$PATH_BIN\\download_platform.pl $platform $file1 $taxid1 -data $PATH_DATA";
			#print OUT "$command\n";
			print OUT `$command`;
		}else{
			my $command = "$PATH_BIN/download_miRNA_platform.pl $platform $file1 $taxid1 -data $PATH_DATA";
			#my $command = "$PATH_BIN\\download_miRNA_platform.pl $platform $file1 $taxid1 -data $PATH_DATA";
			#print OUT "$command\n";
			print OUT `$command`;
		}
		`chmod 666 $PATH_DATA/$file1`;
	}
	my ($countAll,$countSymbols)=(0,0);
	if(open(INFO1, "$PATH_DATA/$file1")){
		my $line = <INFO1>;
		while($line = <INFO1>){
			$countAll++;
			my($probeID,$symbol,$title1) = split(/\t/,$line);
			if($symbol){ $countSymbols++; }
		}
		close INFO1;
	}
	if($countAll>=10 && $logfile && ($type eq "miRNA" || $countSymbols>=10)){
		print OUT "$platform\t$rows\t$type\t$taxid\t$technology\n";
	}
}
close OUT1;
if($logfile){ close OUT; }
exit(0);

#***********************************
sub  check_existing_platforms
#***********************************
{
open(INFO, "$PATH_DATA/$platforms_file");
my $line = <INFO>;
while(my $line = <INFO>){
	chop $line;
	my($ID,$title,$technology,$taxonomy,$rows,$taxid,$type) = split(/\t/,$line);
	my $orgMin;
	foreach my $org (split(/;/,$taxid)){
		my $ii = $hashOrganismID{$taxid};
		if($ii || $orgMin>$ii){ $orgMin=$ii; }
	}
	if(!$orgMin || $orgMin>19){ next; }
	my ($countAll,$countSymbols)=(0,0);
	my $file = "public-$ID"."_annot.txt";
	my $loaded=0;
	if(!file_exist($file)){
		my $taxid1 = $taxid;
		$taxid1 =~ s/;/,/g;
		my $command = "$PATH_BIN/download_platform.pl $ID $file $taxid1 -data $PATH_DATA";
		if($type eq "miRNA"){
			$command = "$PATH_BIN/download_miRNA_platform.pl $ID $file $taxid1 -data $PATH_DATA";
		}
		`$command`;
		$loaded=1;
	}
	if(open(INFO1, "$PATH_DATA/$file")){
		my $line = <INFO1>;
		while($line = <INFO1>){
			$countAll++;
			my($probeID,$symbol,$title1) = split(/\t/,$line);
			if($symbol){ $countSymbols++; }
		}
		close INFO1;
	}
	print OUT join("\t",$ID,$title,$technology,$rows,$taxid,$type,$countAll,$countSymbols)."\n";
}
close INFO;
return;
}

#***********************************
sub file_exist
#***********************************
{
	if(open (INFO_TEMP, $_[0])){ close INFO_TEMP; 1; }
	else { 0; }
}


