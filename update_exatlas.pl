#!/usr/local/bin/perl

#MSIGDB database is manually downloaded from http://thebiogrid.org/download.php
#File is unzipped and renamed as msigdb.entrez.gmt
#BIOGRID database (by species, tab2) is manually downloaded from http://thebiogrid.org/download.php

use strict;
my $arg=0;
my $start = 0;
my $update_option=0;
my $skip_option=0;
my $logfile;
my ($PATH_INFO,$PATH_DATA);
while(my $option = $ARGV[$arg++]){
	if(uc($option) eq "-DATA"){ $PATH_DATA = $ARGV[$arg++] or error_message("path data missing",1); }
	elsif(uc($option) eq "-INFO"){ $PATH_INFO = $ARGV[$arg++] or error_message("path info missing",1); }
	elsif(uc($option) eq "-UPDATE"){ $update_option = $ARGV[$arg++] or error_message("update option missing",1); }
	elsif(uc($option) eq "-SKIP"){ $skip_option = $ARGV[$arg++] or error_message("skip option missing",1); }
	elsif(uc($option) eq "-LOG"){ $logfile = $ARGV[$arg++] or error_message("logfile missing",1); }
	else{ error_message("ERROR: Wrong option $option",1); }
}
if(!$PATH_DATA || !$PATH_INFO){ error_message("update_exatlas.pl -data path -info path [-log, -update, -skip]",1); }
my ($UPDATE_ALL,$UPDATE_GEO,$UPDATE_NCBI,$UPDATE_PLATFORM_LIST,
   $UPDATE_BIOGRID,$UPDATE_MSIGDB,$UPDATE_PUBLIC_DATA,$UPDATE_PLATFORM_ANNOTATIONS)=(0..7);
my ($SKIP_BIOGRID,$SKIP_MSIGDB,$SKIP_PUBLIC_DATA)=(1..3);

my $debug = 0;
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
my %hashOrganism;
for(my $i=0; $i<@organisms; ++$i){
	$hashOrganism{$organisms[$i]->[0]} = $i+1;
}
if($logfile){ open(OUTLOG, ">>$logfile") or error_message("Can't open $!",1); }
open(FTPLOG, ">>$PATH_DATA/FTP_log.txt") or error_message("Can't open $!",1);

my @geneAnnot;
my %symbolsUC;
my %hashAcc;
if(!$update_option || $update_option==$UPDATE_NCBI){
	update_NCBI_genes();
	update_ensembl_annotations();
	parse_geneinfo(1);
	parse_gene2accession(1);
	update_symbol_annotations();
	update_GO_genesets();
	`rm $PATH_DATA/temp6*`;
	if($update_option){ exit(0); }
}elsif($update_option==$UPDATE_BIOGRID || $update_option==$UPDATE_MSIGDB){
	parse_geneinfo(0);
	parse_gene2accession(0);
}
my @homologene=();
my %gene2homologene=();
if((!$update_option || $update_option==$UPDATE_BIOGRID) && $skip_option !~ /$SKIP_BIOGRID/){
	#First upload manually the new version!
	load_homologene(\@homologene,%gene2homologene);
	update_BIOGRID_protein_interaction();
	`rm $PATH_DATA/temp6*`;
	if($update_option){ exit(0); }
}
if((!$update_option || $update_option==$UPDATE_MSIGDB) && $skip_option !~ /$SKIP_MSIGDB/){
	#First upload manually the new version!
	load_homologene(\@homologene,%gene2homologene);
	update_MSIGDB_genesets();
	`rm $PATH_DATA/temp6*`;
	if($update_option){ exit(0); }
}
my %hashPlatform;
if((!$update_option || $update_option>=$UPDATE_PUBLIC_DATA) && $skip_option !~ /$SKIP_PUBLIC_DATA/){
	load_platform_annotations(\%hashPlatform);
	if($update_option==$UPDATE_PLATFORM_ANNOTATIONS){ exit(0); }
	update_public_data_files(\%hashPlatform,0);
	`rm $PATH_DATA/temp6*`;
}
if($logfile){ print OUTLOG "Task completed\n"; close OUTLOG; }
close FTPLOG;
exit(0);

#***********************************
sub file_exist
#***********************************
{
if(open (INFO_TEMP, $_[0])){ close INFO_TEMP; 1; }
else { 0; }
}

#***********************************
sub error_message
#***********************************
{
if($logfile){ print OUTLOG $_[0]."\n"; }
else{ print $_[0]."\n"; }
if($_[1]){
	if($logfile){ print OUTLOG "task_stopped\n"; }
	print "Task stopped\n"; exit(0);
}
return;
}

#***********************************
sub  load_homologene
#***********************************
{
my $homologene = shift;
my $gene2homologene = shift;

if(!open(INFO, "$PATH_DATA/homologene.data")){ error_message("Not found homologene",1); }
while(my $line = <INFO>){
	chop $line;
	if(!$line){ next; }
	my($homologeneID,$speciesID,$geneID,$symbol,$proteinID,$protein) = split(/\t/,$line);
	my $ii = $hashOrganism{$speciesID};
	if(!$ii){ next; }
	$gene2homologene{"$speciesID,$geneID"} = $homologeneID;
	$homologene[$homologeneID]->[$ii-1] = $geneID;
}
close INFO;
if($logfile){ print OUTLOG "Loaded Homologene\n"; }
return;
}

{ #Start a static block
my %hashSymbol;
my %hashAlies;
my $orgID_old;
#***********************************
sub  update_symbols_in_geneset
#***********************************
{
my $file_geneset=shift;
my $orgID=shift;

if(!$orgID){ error_message("Cannot update $file_geneset"); return; }
if($orgID ne $orgID_old){
	$orgID_old = $orgID;
	if(!open (INFO, "$PATH_DATA/gene_info_$orgID.txt")){ error_message("Cannot update $file_geneset: gene_info_$orgID is missing"); return; }
	while(my $line = <INFO>){
		my ($orgID1,$geneID,$symbol,$junk,$alies_list,$junk1,$chr,$band,$title)=split(/\t/, $line);
		$hashSymbol{uc($symbol)} = $symbol;
		if($alies_list !~ /^-+$/){
			my @items=split(/\|/, $alies_list);
			foreach my $alies (@items){
				if(!$hashAlies{$alies}){
					$hashAlies{$alies} = $symbol;
				}
			}
		}
	}
	close INFO;
}
if(!%hashSymbol){ error_message("Cannot update $file_geneset: No symbols"); return; }
if(!open (INFO, "$PATH_DATA/$file_geneset")){ error_message("Cannot update $file_geneset: file not found",1); }
if(!open (OUT, ">$PATH_DATA/temp633.txt")){ error_message("Cannot update $file_geneset: cannot open output file",1); }
my $changed=0;
my ($Ntotal,$Nchanged,$Nremoved)=(0,0,0);
while(my $line = <INFO>){
	if($line=~/^[#!]/){
		print OUT $line;
		next;
	}
	chop $line;
	$line =~ s/\s$//;
	my @items = split(/\t/, $line);
	if(!$items[0]){
		print OUT "$line\n";
		next;
	}
	$Ntotal += @items-2;
	for(my $i=@items-1; $i>=2; $i--){
		my $symbol = $items[$i];
		my $symNew = $hashSymbol{uc($symbol)};
		if($symNew && $symbol ne $symNew){
			$items[$i] = $symNew; $changed=1;
			$Nchanged++;
		}
	}
	if(@items>2){
		print OUT join("\t",@items)."\n";
	}
}
close INFO;
close OUT;
if($logfile){ print OUTLOG "$file_geneset\tgeneset\t$orgID\t$Ntotal\t$Nchanged\n"; }
if($changed){
	#`copy $PATH_DATA\\temp634.txt $PATH_DATA\\$file_geneset`;
	`cp $PATH_DATA/temp634.txt $PATH_DATA/$file_geneset`;
	`rm $PATH_DATA/temp634.txt`;
}
return;
}
} #End of static block

{ #Start a static block
my %hashSymbol;
my %hashAlies;
my %hashEntrez;
my %hashGenbank;
my $orgID_old;
#***********************************
sub  update_symbols_in_arrays
#***********************************
{
my $file_anova=shift;
my $orgID=shift;

if(!$orgID){ error_message("Cannot update $file_anova"); return; }
if($orgID ne $orgID_old){
	$orgID_old = $orgID;
	foreach my $organismID (split(/[;,]/,$orgID)){
		if(open(INFO, "$PATH_DATA/gene_info_$organismID.txt")){
			while(my $line = <INFO>){
				my ($orgID1,$geneID,$symbol,$junk,$alies_list,$junk1,$chr,$band,$title)=split(/\t/, $line);
				if(!$hashSymbol{$symbol}){
					$hashSymbol{$symbol} = [$title,$geneID];
				}
				if(!$hashEntrez{$geneID}){
					$hashEntrez{$geneID} = $symbol;
				}
				if($alies_list !~ /^-+$/){
					my @items=split(/\|/, $alies_list);
					foreach my $alies (@items){
						if(!$hashAlies{$alies}){
							$hashAlies{$alies} = $symbol;
						}
					}
				}
			}
			close INFO;
		}
		if(open (INFO1, "$PATH_DATA/symbol2refseq_$organismID.txt")){
			while(my $line = <INFO1>){
				chop $line;
				$line =~ s/\s$//;
				my ($symbol,$title,$genbank)=split(/\t/,$line);
				foreach my $acc (split(/,/,$genbank)){
					$hashGenbank{$acc} = $symbol;
				}
			}
			close INFO1;
		}
	}
}
if(!%hashSymbol){ error_message("Cannot update $file_anova: No symbols"); return; }
if(!open (INFO, "$PATH_DATA/$file_anova")){ error_message("Cannot update $file_anova: file not found",1); return; }
if(!open (OUT, ">$PATH_DATA/temp633.txt")){ error_message("Cannot open output temp633",1); close INFO; return; }
my ($isymbol,$igenbank,$ientrez,$iensembl,$ientrezNew,$Ncol)=(0,0,0,0,0);
my ($Nadded,$Nchanged,$Ntotal) = (0,0,0,0);
my $update=0;
while(my $line = <INFO>){
	chop $line;
	$line =~ s/\s$//;
	if(!$line){ next; }
	if($line=~/^[\#\!]/){
		print OUT $line."\n";
		next;
	}
	my @items = split(/\t/, $line);
	my $changed=0;
	my $symbolChanged=0;
	for(my $i=0; $i<@items; $i++){
		$items[$i] =~ s/^[\s-\/]+//;
		$items[$i] =~ s/^(null|#*n\/*d)$//i;
		if($items[$i] =~ /\/\/\//){
			my @items1 = split(/\s*\/\/\/\s*/,$items[$i]); $items[$i]=$items1[0];
		}
		if($items[$i] =~ /\/\//){
			my @items1 = split(/\s*\/\/\s*/,$items[$i]); $items[$i]=$items1[0];
		}
		if($items[$i] =~ /\|/){
			my @items1 = split(/\s*\|\s*/,$items[$i]); $items[$i]=$items1[0];
		}
	}
	if(!$isymbol){
		$Ncol = @items;
		while($isymbol<@items && $items[$isymbol] !~ /symbol/i){ $isymbol++; }
		if($isymbol>=@items){ error_message("Symbols missing in $file_anova"); close INFO; close OUT; return; }
		my $i=$isymbol+1;
		while($i<@items && $items[$i] !~ /genbank|^acc|gb_acc|refseq/i){ $i++; }
		if($i<@items){ $igenbank=$i; }
		$i=$isymbol+1;
		while($i<@items && $items[$i] !~ /^entrez|^gene[_ ]id$/i){ $i++; }
		if($i<@items){ $ientrez=$i; }
		$i=$isymbol+1;
		while($i<@items && $items[$i] !~ /^ensembl/i){ $i++; }
		if($i<@items){ $iensembl=$i; }
		print OUT $line;
		if(!$ientrez && ($isymbol||$igenbank)){ print OUT "\tEntrez"; $ientrezNew=$Ncol++; $update=1; }
		print OUT "\n";
		next;
	}
	my $symbol = $items[$isymbol];
	if($symbol =~ /^n\/*a$|^NULL$/i){ $symbol =""; $items[$isymbol] =""; $changed=1; if($debug){print OUTLOG "A1 $symbol\n";} }
	if($items[$isymbol+1] =~ /^n\/*a$|^NULL$/i){ $items[$isymbol+1] =""; $changed=1; if($debug){print OUTLOG "A2 $symbol\n";} }
	my $symbolGenbank;
	my $symbolEntrez;
	#print "A2 $ientrez $igenbank $items[$ientrez] $items[$igenbank]\n";
	if($ientrez && $items[$ientrez]){
		my @entrez = split(/[,\|\/]+/,$items[$ientrez]);
		if(remove_redundant(\@entrez)){
			$items[$ientrez] = join(',',@entrez);
			$changed = 1; if($debug){print OUTLOG "A3 $symbol\n";} 
		}
		foreach my $id (@entrez){
			my $symbol1 = $hashEntrez{$id};
			if($symbol1 && ($symbol1 eq $symbol || !$symbolEntrez)){ $symbolEntrez=$symbol1; }
		}
	}
	if($igenbank && $items[$igenbank]){
		if($items[$ientrez]=~ /^gb[:\|\/]/){ $items[$ientrez]=~ s/^gb[:\|\/]//; $changed = 1; if($debug){print OUTLOG "A4 $symbol\n";} }
		my @genbank = split(/,/,$items[$igenbank]);
		if(remove_redundant(\@genbank)){
			$items[$igenbank] = join(',',@genbank);
			$changed = 1; if($debug){print OUTLOG "A5 $symbol\n";} 
		}
		foreach my $acc (@genbank){
			$acc =~ s/\.\d+$//;
			my $symbol1 = $hashGenbank{$acc};
			if($symbol1 && ($symbol1 eq $symbol || !$symbolGenbank)){ $symbolGenbank=$symbol1; }
		}
	}
	if($iensembl && $items[$iensembl]){
		if($items[$iensembl]=~ /^\d+.+ENS/){ $items[$iensembl]=~ s/^\d+.+ENS/ENS/; $changed = 1; if($debug){print OUTLOG "A6 $symbol\n";} }
	}
	#if($symbol && ($symbolGenbank || $symbolEntrez)){
	#	print "A1 $symbol - $symbolGenbank - $symbolEntrez\n";
	#}
	if($symbolGenbank =~ /^n\/*a$/i){ $symbolGenbank=""; }
	if($symbolEntrez =~ /^n\/*a$/i){ $symbolEntrez=""; }
	if(!$symbol){
		my $sym1 = $items[0];
		$sym1 =~ s/_at$//;
		$sym1 =~ s/_[sx]$//;
		if($hashSymbol{$sym1}){
			$symbol=$sym1;
			$items[$isymbol] = $sym1;
		}elsif($hashGenbank{$sym1}){
			$symbol=$hashGenbank{$sym1};
			if($igenbank && !$items[$igenbank]){ $items[$igenbank]=$sym1; }
		}
		if($symbol){
			$symbolChanged=1; $Nadded++;
			$items[$isymbol+1] = "";
		}
	}
	if($symbol){
		if($symbol =~ /^["']/){ $symbol =~ s/^["']\s*//; $symbol =~ s/\s*["']$//; $symbolChanged=1; }
		if($symbol =~ /[:|;]/){ $symbol =~ s/[:|;].+$//; $symbolChanged=1; }
		if(!$hashSymbol{$symbol}){
			if($symbol =~ /^\d+\-mar$/i){
				$symbol =~ s/-mar$//i;
				if($items[$isymbol+1] =~ /matrix/i && $symbol<3){ $symbol = "Mar".$symbol; }
				else{ $symbol = "March".$symbol; }
				$symbolChanged=1;
			}
			if($symbol =~ /^\d+\-sep$/i){
				$symbol =~ s/-sep$//i;
				if($symbol==15){ $symbol = "Sep".$symbol; }
				else{ $symbol = "Sept".$symbol; }
				$symbolChanged=1;
			}
			my $symbol1 = $symbol; $symbol1 =~ s/orf/ORF/;
			if($symbol =~ /\-mir\-\d|^mir\d/i){ }
			elsif(($orgID>9000 && $orgID<=9999) && $symbol ne uc($symbol) && $symbol1 ne uc($symbol)){
				$symbol = uc($symbol); $symbolChanged=1; if($debug){print OUTLOG "A11 $symbol\n";} 
			}elsif($orgID==10090 || $orgID==10116){
				if($symbol =~ /\d\dRIK$/){
					$symbol =~ s/RIK$/Rik/; $symbolChanged=1; if($debug){print OUTLOG "A7 $symbol\n";}
				}
				elsif($symbol =~ /\d\d\d\d+|^(LOC|MGM)\d+$|^COX\d|^ND\d|^CYT|^nd$|\(ROSA\)|\d\d[rR]ik$|^ATP\d|^l\d+\w/){ }
				elsif($symbol =~ /^D[\dxXyY]/ && $symbol =~ /\w\d|ertd|zems|chr|xiu[mx]|xpas|xbay|xsmh|xyhgu/i){ }
				elsif($symbol =~ /^\w+\(.+\)/){ }
				elsif($symbol =~ /^B*\d/ && $symbol =~ /\w\d|\||\d\./ || length($symbol)<4){ }
				elsif($symbol =~ /^D(el|p)\([hxy\d]+/i){ }
				elsif($symbol =~ /^H2\./){ }
				elsif($symbol =~ /^(H2|Ig[ghkl]|L1Md|Mj|Mt|RatNPn|rp\d+|rt\d+|gb_rt\d|rsa|tcr\w|trav.+)\-/i){ }
				elsif($symbol =~ /^L1Md\-|^tcl(tuw|lub|pa)|^V\-\d\d|^V\d\d+|^n\-r\d+/i){ }
				else{
					my $second = substr($symbol,1,length($symbol)-1);
					my $first = uc(substr($symbol,0,1));
					if($second ne lc($second)){
						$second = lc($second); $symbol=$first.$second; 
						$symbolChanged=1; if($debug){print OUTLOG "A8 $symbol\n";}
					}
				}
			}
		}
		$Ntotal++;
		if(!$hashSymbol{$symbol}){
			if($symbolGenbank){
				if($symbol ne $symbolGenbank){
					$symbol = $symbolGenbank; $symbolChanged = 1;
					$items[$isymbol+1] = "";
					$symbolChanged = 1; if($debug){print OUTLOG "A9 $symbol\n";} 
				}
			}elsif($symbolEntrez){
				if($symbol ne $symbolEntrez){
					$symbol = $symbolEntrez; $symbolChanged = 1;
					$items[$isymbol+1] = "";
					$symbolChanged = 1; if($debug){print OUTLOG "A10 $symbol\n";} 
				}
			}else{
				my $symbol1 = $hashAlies{$symbol};
				if($symbol1){
					$symbol = $symbol1;
					$items[$isymbol+1] = "";
					$symbolChanged = 1; if($debug){print OUTLOG "A12 $symbol\n";} 
				}
			}
		}
		if(!$items[$isymbol+1] && $hashSymbol{$symbol} && $hashSymbol{$symbol}->[0]){
			$items[$isymbol+1] = $hashSymbol{$symbol}->[0];
		}
	}elsif($symbolGenbank){
		$symbol = $symbolGenbank; $symbolChanged = 1; if($debug){print OUTLOG "A13 $symbol\n";} 
		$Nadded++;
		$Ntotal++;
	}elsif($symbolEntrez){
		$symbol = $symbolEntrez; $symbolChanged = 1; if($debug){print OUTLOG "A14 $symbol\n";} 
		$Nadded++;
		$Ntotal++;
	}
	if($symbol && $hashSymbol{$symbol} && $ientrezNew){
		$changed = 1; if($debug){print OUTLOG "A15 $symbol\n";} 
		$items[$ientrezNew] = $hashSymbol{$symbol}->[1];
	}
	if($symbolChanged && $symbol){
		$Nchanged++;
		$changed = 1; if($debug){print OUTLOG "A16 $symbol\n";} 
		$items[$isymbol] = $symbol;
		if($hashSymbol{$symbol}){
			$items[$isymbol+1] = $hashSymbol{$symbol}->[0];
		}
	}
	if($changed){
		print OUT join("\t",@items)."\n";
		$update=1;
	}else{
		print OUT $line."\n";
	}
}
close INFO;
close OUT;
if($update){
	`cp $PATH_DATA/temp633.txt $PATH_DATA/$file_anova`;
	#`copy $PATH_DATA\\temp633.txt $PATH_DATA\\$file_anova`;
	`rm $PATH_DATA/temp633.txt`;
}
my $file_type = "anova";
if($file_anova =~ /_annot.txt$/){ $file_type = "platform"; }
if($logfile){ print OUTLOG "$file_anova\t$file_type\t$orgID\t$Ntotal\t$Nchanged\t$Nadded\n"; }
return;
}
} #Close static block

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

#***********************************
sub update_NCBI_genes
#***********************************
{
use Net::FTP;
my @NCBIfiles = ("gene2ensembl","gene2unigene","gene2go","gene_info","gene2accession");
my $host = "ftp.ncbi.nlm.nih.gov";
my $user = "anonymous";
my $password = "user\@comcast.net";
my $f;

$f = Net::FTP->new($host, Port => 21, Passive => 0) or error_message("Can't open $host",1);
$f->login($user, $password) or error_message("Can't log $user in");
$f->binary();
$f->cwd("/gene/DATA") or error_message("Can't cwd to /gene/DATA");
my $date0 = `date "+%s"`;
my @filesLoaded;
foreach my $file (@NCBIfiles){
	if(file_exist("$PATH_DATA/$file.txt")){
		my $fileDate = `find $PATH_DATA/$file.txt -maxdepth 0 -printf "%Ts"`;
		if($date0 < $fileDate+864000){ next; } #10 days
	}
	push(@filesLoaded, $file);
	if($file =~ /^gene2unigene/){
		$f->get($file, "$PATH_DATA/$file") or error_message("cannot load $file",1);
	}else{
		$f->get("$file.gz", "$PATH_DATA/$file.gz") or error_message("cannot load $file",1);
	}
	if($logfile){ print OUTLOG "File $file loaded\n"; }
}
if(file_exist("$PATH_DATA/homologene.data")){
	my $fileDate = `find $PATH_DATA/homologene.data -maxdepth 0 -printf "%Ts"`;
	if($date0 >= $fileDate+864000){
		$f->cwd("/pub/HomoloGene") or error_message("Can't cwd to HomoloGene",1);
		my @versions = sort $f->ls;
		while($versions[@versions-1] !~ /^build/){ pop(@versions); }
		$f->cwd("/pub/HomoloGene/$versions[@versions-1]") or error_message("Can't cwd to homologene",1);
		$f->get("homologene.data", "$PATH_DATA/homologene.data") or error_message("$versions[@versions-1] cannot load homologene",1);
		if($logfile){ print OUTLOG "File homologene.data loaded\n"; }
	}
}
$f->quit();
foreach my $file (@filesLoaded){
	if($file !~ /^gene2unigene/){
		`gzip -df $PATH_DATA/$file.gz`;
	}
	`mv $PATH_DATA/$file $PATH_DATA/$file.txt`;
	`chmod 666 $PATH_DATA/$file.txt`;
}
if(@filesLoaded){
	if($logfile){ print OUTLOG "NCBI info updated\n"; }
	print "NCBI info updated<br>\n";
}

$host = "ftp.ensembl.org";
$f = Net::FTP->new($host, Port => 21, Passive => 0) or error_message("Can't open $host",1);
$f->login($user, $password) or error_message("Can't log $user in",1);
$f->binary();
my @ENSfiles;
foreach my $ref (@organisms){
	my $organismID = $ref->[0];
	my $species_name = lc($ref->[2]);
	$species_name =~ s/ /_/g;
	if(file_exist("$PATH_DATA/ensembl_gtf_$organismID.txt")){
		my $fileDate = `find $PATH_DATA/ensembl_gtf_$organismID.txt -maxdepth 0 -printf "%Ts"`;
		if($date0 < $fileDate+864000){ next; } #10 days
	}
	if($f->cwd("/pub/current_gtf/$species_name/")){
		my @files = $f->ls;
		while(@files && $files[0] !~ /\d\.gtf\.gz$/){ shift(@files); }
		if(!@files){ error_message("Cannot find GTF for $species_name"); next; }
		$f->get("$files[0]", "$PATH_DATA/ensembl_gtf_$organismID.txt.gz");
		push(@ENSfiles,"ensembl_gtf_$organismID.txt");
		if($logfile){ print OUTLOG "Ensembl for $species_name loaded\n"; }
	}
}
$f->quit();
foreach my $file (@ENSfiles){
	`gzip -df $PATH_DATA/$file.gz`;
	`chmod 666 $PATH_DATA/$file`;
}
if(@ENSfiles){
	if($logfile){ print OUTLOG "ENSEMBL info updated\n"; }
	print "ENSEMBL info updated<br>\n";
}
return;
}

#***********************************
sub parse_geneinfo
#***********************************
{
my $option=shift;
my $id_old;
if($logfile){ print OUTLOG "Parsing gene_info...\n"; }
if(!open (INFO, "$PATH_DATA/gene_info.txt")){ error_message("Not found gene_info",1); }
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^[#!]/){ next; }
	my ($organismID,$geneID,$symbol,$junk,$alias,$alternative,$junk1,$junk2,$title,$geneType)=split(/\t/, $line);
	my $ii = $hashOrganism{$organismID};
	if(!$ii){ next; }
	$symbolsUC{uc($symbol)}->[$ii-1] = [$symbol,$geneID];
	$geneAnnot[$ii-1]->{$geneID} = [$symbol,$alias,$alternative,$title];
	if($organismID != $id_old && $option){
		if($id_old){
			close OUT;
			`cp $PATH_DATA/temp632.txt $PATH_DATA/gene_info_$id_old.txt`;
		}
		if(!open (OUT, ">$PATH_DATA/temp632.txt")){ error_message("Not open temp632"); return; }
	}
	if($option){ print OUT "$line\n"; }
	$id_old = $organismID;
}
close INFO;
if($option){ close OUT; }
if($id_old && $option){
	`cp $PATH_DATA/temp632.txt $PATH_DATA/gene_info_$id_old.txt`;
}
if($logfile){ print OUTLOG "Gene_info parsed\n"; }
return;
}

#***********************************
sub update_BIOGRID_protein_interaction
#***********************************
{
if($logfile){ print OUTLOG "Processing BIOGRID\n"; }
my %hashBIOGRID;
my %hashOrgBIOGRID;
if(!@homologene){ error_message("Missing homologene in BIOGRID",1); }
if(!%gene2homologene){ error_message("Missing gene2homologene in BIOGRID",1); }
if(!@geneAnnot){ error_message("Missing geneAnnot in BIOGRID",1); return; }
if(!open (INFO, "$PATH_DATA/BIOGRID-ALL.tab2.txt")){ error_message("Not found BIOGRID",1); }
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^[#!]/){ next; }
	my ($BioGRID,$GeneIDa,$GeneIDb,$BGidInta,$BGidIntb,$NameInta,$NameIntb,$SymbolA,$SymbolB,$AliasA,$AliasB,$ExperSys,$ExperType,$Author,$Pubmed,$orgIDa,$orgIDb)=split(/\t/, $line);
	my $iia = $hashOrganism{$orgIDa};
	my $iib = $hashOrganism{$orgIDb};
	if(!$iia || !$iib){ next; }
	for(my $ii=0; $ii<@organisms; $ii++){
		my $organismID = $organisms[$ii]->[0];
		my $ii = $hashOrganism{$organismID};
		my ($symbolA, $symbolB, $geneIDa1, $geneIDb1);
		my $homologeneID = $gene2homologene{"$orgIDa,$GeneIDa"};
		if($homologeneID && $homologene[$homologeneID]->[$ii-1]){
			$geneIDa1 = $homologene[$homologeneID]->[$ii-1];
			my $ref = $geneAnnot[$ii-1];
			if(!$ref || ref($ref) ne 'HASH'){ error_message("Bad geneAnnot1 BIOGRID",1); }
			my $ref1 = $ref->{$geneIDa1};
			if($ref1 && ref($ref1) eq 'ARRAY'){
				$symbolA = $ref1->[0];
			}
		}
		if(!$symbolA){
			my $ref = $geneAnnot[$iia-1];
			if(!$ref || ref($ref) ne 'HASH'){ error_message("Bad geneAnnot2 BIOGRID",1); }
			my $ref1 = $ref->{$geneIDa1};
			if($ref1 && ref($ref1) eq 'ARRAY'){
				my $symbol1 = uc($ref1->[0]);
				my $ref2 = $symbolsUC{$symbol1};
				if($ref2 && $ref2->[$ii-1]){
					($symbolA,$geneIDa1) = @{$ref2->[$ii-1]};
				}
			}
		}
		if(!$symbolA){ next; }
		$homologeneID = $gene2homologene{"$orgIDb,$GeneIDb"};
		if($homologeneID && $homologene[$homologeneID]->[$ii-1]){
			$geneIDb1 = $homologene[$homologeneID]->[$ii-1];
			my $ref = $geneAnnot[$ii-1]->{$geneIDb1};
			if($ref && ref($ref) eq 'ARRAY'){ $symbolB = $ref->[0]; }
		}
		if(!$symbolB){
			my $ref = $geneAnnot[$iib-1];
			if(!$ref || ref($ref) ne 'HASH'){ error_message("Bad geneAnnot2 BIOGRID",1); }
			my $ref1 = $ref->{$geneIDb1};
			if($ref1 && ref($ref1) eq 'ARRAY'){
				my $symbol1 = uc($ref1->[0]);
				my $ref2 = $symbolsUC{$symbol1};
				if($ref2 && $ref2->[$ii-1]){
					($symbolB,$geneIDb1) = @{$ref2->[$ii-1]};
				}
			}
		}
		if(!$symbolB){ next; }
		if($ii<=12 && ($iia>12 || $iib>12)){ next; }
		if($ii>12 && $iia != $ii && $iib != $ii){ next; }
		if(!$hashOrgBIOGRID{$organismID}){ $hashOrgBIOGRID{$organismID}=1; }
		my $ref1 = $hashBIOGRID{"$organismID,$geneIDa1"};
		if(!$ref1){
			$hashBIOGRID{"$organismID,$geneIDa1"}->{$symbolB}=1;
		}elsif(!$ref1->{$symbolB}){
			$ref1->{$symbolB}=1;
		}
		my $ref2 = $hashBIOGRID{"$organismID,$geneIDb1"};
		if(!$ref2){
			$hashBIOGRID{"$organismID,$geneIDb1"}->{$symbolA}=1;
		}elsif(!$ref2->{$symbolA}){
			$ref2->{$symbolA}=1;
		}
	}
}
close INFO;
my $found=0;
for(my $ii=0; $ii<@organisms; $ii++){
	my $organismID = $organisms[$ii]->[0];
	if(!open (OUT, ">$PATH_DATA/temp631.txt")){ error_message("Not open temp631"); return; }
	print OUT "!Protein Interaction (BioGRID) gene sets for $organisms[$ii]->[2]\n";
	my $Nsets=0;
	foreach my $geneID (sort keys %{$geneAnnot[$ii]}){
		my $ref = $hashBIOGRID{"$organismID,$geneID"};
		if(!$ref){ next; }
		my @symbol_list = keys %$ref;
		if(@symbol_list < 5 || @symbol_list > 3000){ next; }
		print OUT "$geneAnnot[$ii]->{$geneID}->[0]-interaction\t\t".join("\t",sort @symbol_list)."\n";
		$Nsets++;
	}
	close OUT;
	if($Nsets > 10){
		my $file_name = "$PATH_DATA/public-Protein_interaction_$organismID";
		$found=1;
		`cp -f $PATH_DATA/temp631.txt $file_name`;
	}
}
if($found){
	`chmod 666 $PATH_DATA/public-Protein_interaction*`;
}
return;
}

#***********************************
sub   update_ensembl_annotations
#***********************************
{
my $id_old;
my $Nspecies=0;
foreach my $organismID (keys %hashOrganism){
	if(!open(INFO, "$PATH_DATA/ensembl_gtf_$organismID.txt")){ next; }
	if(!open (OUT, ">$PATH_DATA/temp542.txt")){ error_message("Not open temp542",1); }
	print OUT "ProbeID\tsymbol\tgene_title";
	my $ii = $hashOrganism{$organismID};
	my %hashAlias;
	my $ref = $geneAnnot[$ii-1];
	if($ref && ref($ref) eq "HASH"){
		print OUT "\tEntrez";
		foreach my $geneID (keys %$ref){
			my ($symbol,$alias_list,$alternative,$title) = @{$ref->{$geneID}};
			if($alias_list eq "-"){ next; }
			foreach my $alias (split(/\|/,$alias_list)){
				if(!$hashAlias{$alias}){
					$hashAlias{$alias}=$symbol;
				}
			}
		}
	}
	print OUT "\n";
	my $count=0;
	my %hashDone;
	while(my $line = <INFO>){
		chop $line;
		my @items=split(/\t/, $line);
		if($items[2] !~ /^gene|^transcript/){ next; }
		my ($ensGene,$ensTrans,$symbol);
		foreach my $descrip (split(/; /, $items[8])){
			my($term,$value) = split(/ /,$descrip);
			if($term eq "gene_id"){ $ensGene = $value; }
			elsif($term eq "transcript_id"){ $ensTrans = $value; }
			elsif($term eq "gene_name"){ $symbol = $value; last; }
		}

		$ensGene =~ s/\"//g; $ensTrans =~ s/\"//g; $symbol =~ s/\"//g;
		$ensGene =~ s/\.\d+$//; $ensTrans =~ s/\.\d+$//;
		#print "A2 $ensGene\t$ensTrans\t$symbol\n";
		if(!$symbol || $symbol eq "-"){ next; }
		my $refSymbol = $symbolsUC{uc($symbol)};
		if(!$refSymbol){
			my $symbol1 = $hashAlias{$symbol};
			if($symbol1){ $symbol = $symbol1; }
			$refSymbol = $symbolsUC{uc($symbol)};
		}
		my ($symbol1,$symbol2,$alias,$alternative,$title,$geneID);
		if($refSymbol && $refSymbol->[$ii-1] && $ref && ref($ref) eq "HASH"){
			($symbol1,$geneID) = @{$refSymbol->[$ii-1]};
			if($symbol1 && $symbol1 ne $symbol){ $symbol = $symbol1; }
			($symbol2,$alias,$alternative,$title) = @{$ref->{$geneID}};
		}
		if($items[2] eq "gene" && !$hashDone{$ensGene} && $symbol){
			print OUT "$ensGene\t$symbol\t$title\t$geneID\n";
			$hashDone{$ensGene} = 1;
			$count++;
		}
		elsif($items[2] eq "transcript" && !$hashDone{$ensTrans} && $symbol){
			print OUT "$ensTrans\t$symbol\t$title\t$geneID\n";
			$hashDone{$ensTrans} = 1;
			$count++;
		}
	}
	close OUT;
	close INFO;
	if($count > 10){
		my $annotFile = "ensembl_$organismID"."_annot.txt";
		`cp $PATH_DATA/temp542.txt $PATH_DATA/$annotFile`;
		$Nspecies++;
	}
}
close INFO;
if($logfile){ print OUTLOG "Ensembl annotations done for $Nspecies species\n"; }
return;
}

#***********************************
sub   update_GO_genesets
#***********************************
{
if(!@geneAnnot){ error_message("Missing geneAnnot in gene2accession",1);  }
if($logfile){ print OUTLOG "Processing GO annotations\n"; }
my %hashGO;
my %hashGOterm;
my %hashOrgGO;
if(!open (INFO, "$PATH_DATA/gene2go.txt")){ error_message("Missing geneAnnot in gene2go",1); }
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^[#!]/){ next; }
	my ($organismID,$geneID,$GO,$evidence,$qual,$GOterm)=split(/\t/, $line);
	if(!$hashOrganism{$organismID}){ next; }
	if(!$hashOrgGO{$organismID}){ $hashOrgGO{$organismID}=1; }
	if(!$hashGOterm{$GO}){ $hashGOterm{$GO}=$GOterm; }
	my $ref = $hashGO{"$organismID,$GO"};
	if(!$ref){
		$hashGO{"$organismID,$GO"}->{$geneID}=1;
		$ref = $hashGO{"$organismID,$GO"};
	}
	if(!$ref->{$geneID}){
		$ref->{$geneID}=1;
	}
}
close INFO;
my $found = 0;
foreach my $organismID (keys %hashOrgGO){
	my $ii = $hashOrganism{$organismID};
	open (OUT, ">$PATH_DATA/temp631.txt") or error_message("Can't open $!",1);
	print OUT "!Gene Ontology gene sets for $organisms[$ii-1]->[2]\n";
	foreach my $GO (sort keys %hashGOterm){
		my $ref = $hashGO{"$organismID,$GO"};
		if(!$ref || ref($ref) ne 'HASH'){ next; }
		my @id_list = keys %$ref;
		if(@id_list < 5 || @id_list > 3000){ next; }
		my @symbols;
		foreach my $geneID (@id_list){
			my $ref = $geneAnnot[$ii-1]->{$geneID};
			if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot GO genesets",1); }
			push(@symbols,$ref->[0]);
		}
		print OUT "$GO\t$hashGOterm{$GO}\t".join("\t",sort @symbols)."\n";
	}
	close OUT;
	my $file_name = "$PATH_DATA/public-GO_geneset_$organismID";
	$found = 1;
	`cp -f $PATH_DATA/temp631.txt $file_name`;
}
if($found){
	`chmod 666 $PATH_DATA/public-GO_geneset*`;
}
return;
}

#***********************************
sub   update_MSIGDB_genesets
#***********************************
{
if(!@homologene){ error_message("Missing homologene in MSIGDB",1); }
if(!%gene2homologene){ error_message("Missing gene2homologene in MSIGDB",1); }
if(!@geneAnnot){ error_message("Missing geneAnnot in MSIGDB",1); }
if($logfile){ print OUTLOG "Processing MSIGDB\n"; }
my %hashMSIGDB;
my %hashMIR;
my %hashGNF2;
my %hashKEGG;
my %hashREACTOME;
my %hashPID;
my $section=0;
if(!open (INFO, "$PATH_DATA/msigdb.entrez.gmt")){ error_message("Missing msigdb",1); }
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^[#!]/){ next; }
	my ($ID,$html,@entrez)=split(/\t/, $line);
	if($section==0 && $ID =~/^KEGG_/){ $section=1; }
	elsif($section==1 && $ID =~/^MODULE_/){ $section=2; }
	elsif($section==2 && $ID !~/^(REACTOME|PID|MODULE)_/){ $section=3; }
	elsif($section==3 && $ID =~/^GSE\d\d/){ $section=4; }
	if($section==0){
		if($ID=~/,(MIR|LET)-\d/i){ $hashMIR{"$ID\t$html"} = \@entrez; }
		elsif($ID=~/^GNF2/){ $hashGNF2{"$ID\t$html"} = \@entrez; }
	}elsif($section==1){
		if($ID=~/^KEGG_/){ $ID=~s/^KEGG_//; $ID=lc($ID); $hashKEGG{"$ID\t$html"} = \@entrez; }
	}elsif($section==2){
		if($ID=~/^REACTOME_/){ $ID=~s/^REACTOME_//; $ID=lc($ID); $hashREACTOME{"$ID\t$html"} = \@entrez; }
		elsif($ID=~/^PID_/i){ $ID=~s/^PID_//; $ID=lc($ID); $hashPID{"$ID\t$html"} = \@entrez; }
	}elsif($section==3){
		$ID=lc($ID); $hashMSIGDB{"$ID\t$html"} = \@entrez;
	}elsif($section==4){
		if($ID !~/^GSE\d\d/){ $ID=lc($ID); $hashMSIGDB{"$ID\t$html"} = \@entrez; }
	}
}
close INFO;
my $found=0;
for(my $i=0; $i<12; ++$i){
	my @countGenes;
	my $organismID = $organisms[$i]->[0];
	open (OUT, ">$PATH_DATA/temp630.txt") or error_message("Can't open $!",1);
	print OUT "!KEGG pathway database (part of MSIGDB) for $organisms[$i]->[2]\n";
	foreach my $ID (sort keys %hashKEGG){
		my @id_list = @{$hashKEGG{$ID}};
		if(@id_list < 5 || @id_list > 3000){ next; }
		my @symbols;
		foreach my $id_human (@id_list){
			my $homologeneID = $gene2homologene{"9606,$id_human"};
			my $symbol;
			if($homologeneID && $homologene[$homologeneID]->[$i]){
				my $ref = $geneAnnot[$i]->{$homologene[$homologeneID]->[$i]};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				$symbol = $ref->[0];
				$countGenes[0]++;
			}else{
				my $ref = $geneAnnot[0]->{$id_human};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				my $symbol1 = uc($ref->[0]);
				if($symbolsUC{$symbol1}->[$i]){
					$symbol = $symbolsUC{$symbol1}->[$i]->[0];
					$countGenes[1]++;
				}
			}
			if($symbol){
				push(@symbols,$symbol);
			}
		}
		if(@symbols >=4){
			print OUT "$ID\t".join("\t",sort @symbols)."\n";
		}
	}
	close OUT;
	my $file_name = "$PATH_DATA/public-KEGG_pathways_$organismID";
	$found=1;
	`cp -f $PATH_DATA/temp630.txt $file_name`;
}
if($found){
	`chmod 666 $PATH_DATA/public-KEGG_pathways*`;
}
$found=0;
for(my $i=0; $i<12; ++$i){
	my @countGenes;
	my $organismID = $organisms[$i]->[0];
	open (OUT, ">$PATH_DATA/temp630.txt") or error_message("Can't open $!",1);
	print OUT "!PID pathway database (part of MSIGDB) for $organisms[$i]->[2]\n";
	foreach my $ID (sort keys %hashPID){
		my @id_list = @{$hashPID{$ID}};
		if(@id_list < 5 || @id_list > 3000){ next; }
		my @symbols;
		foreach my $id_human (@id_list){
			my $homologeneID = $gene2homologene{"9606,$id_human"};
			my $symbol;
			if($homologeneID && $homologene[$homologeneID]->[$i]){
				my $ref = $geneAnnot[$i]->{$homologene[$homologeneID]->[$i]};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				$symbol = $ref->[0];
				$countGenes[0]++;
			}else{
				my $ref = $geneAnnot[0]->{$id_human};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				my $symbol1 = uc($ref->[0]);
				if($symbolsUC{$symbol1}->[$i]){
					$symbol = $symbolsUC{$symbol1}->[$i]->[0];
					$countGenes[1]++;
				}
			}
			if($symbol){
				push(@symbols,$symbol);
			}
		}
		if(@symbols >=4){
			print OUT "$ID\t".join("\t",sort @symbols)."\n";
		}
	}
	close OUT;
	my $file_name = "$PATH_DATA/public-PID_pathways_$organismID";
	$found=1;
	`cp -f $PATH_DATA/temp630.txt $file_name`;
}
if($found){
	`chmod 666 $PATH_DATA/public-PID_pathways*`;
}
$found=0;
for(my $i=0; $i<12; ++$i){
	my @countGenes;
	my $organismID = $organisms[$i]->[0];
	open (OUT, ">$PATH_DATA/temp630.txt") or error_message("Can't open $!",1);
	print OUT "!MIR pathway database (part of MSIGDB) for $organisms[$i]->[2]\n";
	foreach my $ID (sort keys %hashMIR){
		my @id_list = @{$hashMIR{$ID}};
		if(@id_list < 5 || @id_list > 3000){ next; }
		my @symbols;
		foreach my $id_human (@id_list){
			my $homologeneID = $gene2homologene{"9606,$id_human"};
			my $symbol;
			if($homologeneID && $homologene[$homologeneID]->[$i]){
				my $ref = $geneAnnot[$i]->{$homologene[$homologeneID]->[$i]};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				$symbol = $ref->[0];
				$countGenes[0]++;
			}else{
				my $ref = $geneAnnot[0]->{$id_human};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				my $symbol1 = uc($ref->[0]);
				if($symbolsUC{$symbol1}->[$i]){
					$symbol = $symbolsUC{$symbol1}->[$i]->[0];
					$countGenes[1]++;
				}
			}
			if($symbol){
				push(@symbols,$symbol);
			}
		}
		if(@symbols >=4){
			print OUT "$ID\t".join("\t",sort @symbols)."\n";
		}
	}
	close OUT;
	my $file_name = "$PATH_DATA/public-MIR_targets_$organismID";
	$found=1;
	`cp -f $PATH_DATA/temp630.txt $file_name`;
}
if($found){
	`chmod 666 $PATH_DATA/public-MIR_targets*`;
}
$found=0;
for(my $i=0; $i<12; ++$i){
	my @countGenes;
	my $organismID = $organisms[$i]->[0];
	open (OUT, ">$PATH_DATA/temp630.txt") or error_message("Can't open $!",1);
	print OUT "!REACTOME pathway database (part of MSIGDB) for $organisms[$i]->[2]\n";
	foreach my $ID (sort keys %hashREACTOME){
		my @id_list = @{$hashREACTOME{$ID}};
		if(@id_list < 5 || @id_list > 3000){ next; }
		my @symbols;
		foreach my $id_human (@id_list){
			my $homologeneID = $gene2homologene{"9606,$id_human"};
			my $symbol;
			if($homologeneID && $homologene[$homologeneID]->[$i]){
				my $ref = $geneAnnot[$i]->{$homologene[$homologeneID]->[$i]};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				$symbol = $ref->[0];
				$countGenes[0]++;
			}else{
				my $ref = $geneAnnot[0]->{$id_human};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				my $symbol1 = uc($ref->[0]);
				if($symbolsUC{$symbol1}->[$i]){
					$symbol = $symbolsUC{$symbol1}->[$i]->[0];
					$countGenes[1]++;
				}
			}
			if($symbol){
				push(@symbols,$symbol);
			}
		}
		if(@symbols >=4){
			print OUT "$ID\t".join("\t",sort @symbols)."\n";
		}
	}
	close OUT;
	my $file_name = "$PATH_DATA/public-REACTOME_pathways_$organismID";
	$found=1;
	`cp -f $PATH_DATA/temp630.txt $file_name`;
}
if($found){
	`chmod 666 $PATH_DATA/public-REACTOME_pathways*`;
}
$found=0;
for(my $i=0; $i<12; ++$i){
	my @countGenes;
	my $organismID = $organisms[$i]->[0];
	open (OUT, ">$PATH_DATA/temp630.txt") or error_message("Can't open $!",1);
	print OUT "!MSIGDB pathway database (part of MSIGDB) for $organisms[$i]->[2]\n";
	foreach my $ID (sort keys %hashMSIGDB){
		my @id_list = @{$hashMSIGDB{$ID}};
		if(@id_list < 5 || @id_list > 3000){ next; }
		my @symbols;
		foreach my $id_human (@id_list){
			my $homologeneID = $gene2homologene{"9606,$id_human"};
			my $symbol;
			if($homologeneID && $homologene[$homologeneID]->[$i]){
				my $ref = $geneAnnot[$i]->{$homologene[$homologeneID]->[$i]};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				$symbol = $ref->[0];
				$countGenes[0]++;
			}else{
				my $ref = $geneAnnot[0]->{$id_human};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				my $symbol1 = uc($ref->[0]);
				if($symbolsUC{$symbol1}->[$i]){
					$symbol = $symbolsUC{$symbol1}->[$i]->[0];
					$countGenes[1]++;
				}
			}
			if($symbol){
				push(@symbols,$symbol);
			}
		}
		if(@symbols >=4){
			print OUT "$ID\t".join("\t",sort @symbols)."\n";
		}
	}
	close OUT;
	my $file_name = "$PATH_DATA/public-MSIGDB_genesets_$organismID";
	$found=1;
	`cp -f $PATH_DATA/temp630.txt $file_name`;
}
if($found){
	`chmod 666 $PATH_DATA/public-MSIGDB_genesets*`;
}
$found=0;
for(my $i=0; $i<12; ++$i){
	my @countGenes;
	my $organismID = $organisms[$i]->[0];
	open (OUT, ">$PATH_DATA/temp630.txt") or error_message("Can't open $!",1);
	print OUT "!GNF2 pathway database (part of MSIGDB) for $organisms[$i]->[2]\n";
	foreach my $ID (sort keys %hashGNF2){
		my @id_list = @{$hashGNF2{$ID}};
		if(@id_list < 5 || @id_list > 3000){ next; }
		my @symbols;
		foreach my $id_human (@id_list){
			my $homologeneID = $gene2homologene{"9606,$id_human"};
			my $symbol;
			if($homologeneID && $homologene[$homologeneID]->[$i]){
				my $ref = $geneAnnot[$i]->{$homologene[$homologeneID]->[$i]};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				$symbol = $ref->[0];
				$countGenes[0]++;
			}else{
				my $ref = $geneAnnot[0]->{$id_human};
				if(!$ref || ref($ref) ne 'ARRAY'){ error_message("Bad geneAnnot1 MSIGDB genesets",1); }
				my $symbol1 = uc($ref->[0]);
				if($symbolsUC{$symbol1}->[$i]){
					$symbol = $symbolsUC{$symbol1}->[$i]->[0];
					$countGenes[1]++;
				}
			}
			if($symbol){
				push(@symbols,$symbol);
			}
		}
		if(@symbols >=4){
			print OUT "$ID\t".join("\t",sort @symbols)."\n";
		}
	}
	close OUT;
	my $file_name = "$PATH_DATA/public-GNF2_coregulated_$organismID";
	$found=1;
	`cp -f $PATH_DATA/temp630.txt $file_name`;
}
if($found){
	`chmod 666 $PATH_DATA/public-GNF2_coregulated*`;
}
return;
}

#***********************************
sub   parse_gene2accession
#***********************************
{
my $option=shift;
if($option && !@geneAnnot){ error_message("Missing geneAnnot in gene2accession",1); }
if($logfile){ print OUTLOG "Parsing gene2accession $option\n"; }
my %hash;
if(!open(INFO, "$PATH_DATA/gene2accession.txt")){ error_message("Missing gene2accession",1); }
my $geneID_old;
my $orgID_old;
my @status_list=("VALIDATED","REVIEWED","INFERRED","PROVISIONAL","PREDICTED","MODEL","SUPPRESSED","-","NA");
my %hash_status;
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^[#!]/){ next; }
	my ($organismID,$geneID,$status,$RNAacc)=split(/\t/, $line);
	if($RNAacc eq "-"){ next; }
	if(!$hashOrganism{$organismID}){ next; }
	if(!$orgID_old && $option){
		open(OUT, ">$PATH_DATA/temp621.txt");
		open(OUT1, ">$PATH_DATA/temp622.txt");
		print OUT "Symbol\tGene title\tGenBank\n";
		print OUT1 "ProbeID\tSymbol\tGene title\tEntrez\n";
	}
	if($geneID ne $geneID_old && $geneID_old){
		my @accList;
		my @best;
		foreach my $st (@status_list){
			my $ref = $hash{$st};
			if(!$ref){ next; }
			my @acc = keys %$ref;
			push(@accList,@acc);
			if(!@best){ @best = @acc; }
		}
		if(@best){
			$hashAcc{"$orgID_old,$geneID_old"} = join(',',@best);
		}
		%hash = ();
		if($option){
			my $ref = $geneAnnot[$hashOrganism{$orgID_old}-1];
			if(!$ref || ref($ref) ne 'HASH'){ error_message("Bad geneAnnot1",1); }
			my $ref1 = $ref->{$geneID_old};
			if(!$ref1 || ref($ref1) ne 'ARRAY'){ error_message("Bad geneAnnot2 $geneID_old"); }
			else{
				my ($symbol,$alias,$alternative,$title) = @$ref1;
				print OUT "$symbol\t$title\t".join(',',@accList)."\n";
				foreach my $acc (@accList){
					print OUT1 "$acc\t$symbol\t$title\t$geneID_old\n";
				}
			}
		}
		if($organismID != $orgID_old && $option){
			close OUT;
			close OUT1;
			my $file1 = "$PATH_DATA/genbank_$orgID_old"."_annot.txt";
			`cp $PATH_DATA/temp621.txt $PATH_DATA/symbol2refseq_$orgID_old.txt`;
			`cp $PATH_DATA/temp622.txt $file1`;
			open(OUT, ">$PATH_DATA/temp621.txt");
			open(OUT1, ">$PATH_DATA/temp622.txt");
		}		
	}
	$RNAacc =~ s/\.\d+$//;
	$hash{$status}->{$RNAacc}=1;
	$geneID_old = $geneID;
	$orgID_old = $organismID;
}
close INFO;
if($geneID_old){
	my @accList;
	my @best;
	foreach my $st (@status_list){
		my $ref = $hash{$st};
		if(!$ref){ next; }
		my @acc = keys %$ref;
		push(@accList,@acc);
		if(!@best){ @best = @acc; }
	}
	if(@best){
		$hashAcc{"$orgID_old,$geneID_old"} = join(',',@best);
	}
	if($option){
		my $ref = $geneAnnot[$hashOrganism{$orgID_old}-1];
		if(!$ref || ref($ref) ne 'HASH'){ error_message("Bad geneAnnot1",1); }
		my $ref1 = $ref->{$geneID_old};
		if(!$ref1 || ref($ref1) ne 'ARRAY'){ error_message("Bad geneAnnot2 $geneID_old"); }
		else{
			my ($symbol,$alias,$alternative,$title) = @$ref1;
			print OUT "$symbol\t$title\t".join(',',@accList)."\n";
			foreach my $acc (@accList){
				print OUT1 "$acc\t$symbol\t$title\t$geneID_old\n";
			}
		}
	}
	if($orgID_old && $option){
		close OUT;
		close OUT1;
		my $file1 = "$PATH_DATA/genbank_$orgID_old"."_annot.txt";
		`cp $PATH_DATA/temp621.txt $PATH_DATA/symbol2refseq_$orgID_old.txt`;
		`cp $PATH_DATA/temp622.txt $file1`;
	}		
}
if($logfile){ print OUTLOG "gene2accession parsed\n"; }
return;
}

#***********************************
sub   update_symbol_annotations
#***********************************
{
if(!%hashAcc){ error_message("Missing hashAcc",1); }
if($logfile){ print OUTLOG "Updating symbol annotations\n"; }
for(my $i=0; $i<@organisms; ++$i){
	my $organismID = $organisms[$i]->[0];
	my @aliasSymbols;
	my %hashSymbols;
	my %hashAlias;
	open (OUT, ">$PATH_DATA/temp621.txt") or error_message("Can't open $!",1);
	print OUT "ProbeID\tSymbol\tGene title\tGenBank\tEntrez\n";
	foreach my $geneID (keys %{$geneAnnot[$i]}){
		my ($symbol,$alias_list,$alternative,$title) = @{$geneAnnot[$i]->{$geneID}};
		my $acc = $hashAcc{"$organismID,$geneID"};
		my $text = "\t$symbol\t$title\t$acc\t$geneID\n";
		print OUT $symbol.$text;
		$hashSymbols{$symbol}=1;
		if($alias_list eq "-"){ next; }
		foreach my $alias (split(/\|/,$alias_list)){
			push(@aliasSymbols,$alias.$text);
			push(@{$hashAlias{$alias}},$symbol);
		}
	}
	foreach my $alias (keys %hashAlias){
		my @symbols = @{$hashAlias{$alias}};
		if(@symbols==1){ $hashAlias{$alias} = $symbols[0]; next; }
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
		#print "$alias @symbols\t@similarity\n";
		my @sorted = sort {$similarity[$b]<=>$similarity[$a]} 0..(@symbols-1);
		if($similarity[$sorted[0]]>1){
			#print "$alias $similarity[$sorted[0]] - $symbols[$sorted[0]]\n";
			$hashAlias{$alias} = $symbols[$sorted[0]];
			next;
		}
		while(@symbols>1 && uc(substr($symbols[0],0,2)) ne uc(substr($symbols[1],0,2))){
			shift(@symbols);
		}
		if(@symbols==1){ undef($hashAlias{$alias}); next; }
		my $count=2;
		for(my $i=2; $i<@symbols; $i++){ 
			if(uc(substr($symbols[$i],0,2)) eq uc(substr($symbols[0],0,2))){ $count++; }
			else{ last; }
		}
		#print "DD $alias @symbols $count\n";
		if($count >= @symbols/2){ 
			$hashAlias{$alias} = $symbols[0];
		}else{
			undef($hashAlias{$alias});
		}
	}
	foreach my $text (@aliasSymbols){
		my ($alias,$symbol,$junk) = split(/\t/,$text);
		if($hashSymbols{$alias}){ next; }
		if($hashAlias{$alias} eq $symbol){ print OUT $text; }
#		elsif(!$hashAlias{$alias} || ref($hashAlias{$alias}) eq 'ARRAY'){
#			print $text;
#		}
	}
	close OUT;
	my $file_name = "$PATH_DATA/symbol_$organismID"."_annot.txt";
	`cp -f $PATH_DATA/temp621.txt $file_name`;
}
return;
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
sub load_platform_annotations
#***********************************
{
my $hashPlatform = shift;

my %hash_platform_config;
if(open(INFO, "$PATH_INFO/public-config.txt")){
	while(my $line = <INFO>){
		$line =~ s/\n$//;
		my %hash=();
		read_config_line($line,\%hash);
		my $org = $hash{"organismID"};
		my $array_type = $hash{"array_type"};
		my $platform = $hash{"type_annotation"};
		if($platform){
			$hash_platform_config{$platform} = [$array_type,$org];
		}
	}
	close INFO;
}
my @platforms_missing_config;
if(!open (INFO, "$PATH_DATA/array_platforms.txt")){ error_message("Not found: array_platforms!",1); }
while(my $line = <INFO>){
	chop $line;
	if(!$line || $line=~/^acc/i){ next; }
	my($platformID,$title,$technology,$taxonomy,$rows,$taxid,$type) = split(/\t/,$line);
	$hashPlatform->{$platformID}=[$type,$taxid];
	my $ref1 = $hash_platform_config{$platformID};
	if($ref1){
		if($ref1->[0] ne $type){
			error_message("Wrong array platform type: $ref1->[0] ne $type\t$platformID");
		}
		if($ref1->[1] ne $taxid){
			error_message("Wrong array platform taxid: $ref1->[1] ne $taxid\t$platformID");
		}
	}else{
		push(@platforms_missing_config,[$platformID,$title,$taxid,$type]);
	}
}
close INFO;
if(@platforms_missing_config && file_exist("$PATH_INFO/public-config.txt")){
	open(OUT, ">>$PATH_INFO/public-config.txt");
	foreach my $ref (@platforms_missing_config){
		my ($platformID,$title,$taxid,$array_type) = @$ref;
		if(file_exist("$PATH_DATA/public-$platformID"."_annot.txt")){
			$hash_platform_config{$platformID} = [$array_type,$taxid];
			print OUT "type_annotation=$platformID\tdescription=$title\torganismID=$taxid\tarray_type=$array_type\n";
		}
	}
	close OUT;
}
if($logfile){ print OUTLOG "Platforms loaded\n"; }
return;
}

#***********************************
sub   update_public_data_files
#***********************************
{
my $hashPlatform=shift;
my $option=shift;

if(!$hashPlatform || ref($hashPlatform) ne 'HASH'){ error_message("Missing hashPlatform in update_public_files",1); }
if($logfile){ print OUTLOG "Updating public files\n"; }
my @matrix_list;
my @geneset_list;
my @annotation_list;
if(!open(INFO, "$PATH_INFO/public-config.txt")){ error_message("Missing public-config",1);  }
while(my $line = <INFO>){
	$line =~ s/\n$//;
	my %hash=();
	if(length($line) < 3) { next; }
	read_config_line($line,\%hash);
	my $descr = $hash{"description"};
	my $org = $hash{"organismID"};
	my $type = $hash{"array_type"};
	#if($hash{"type_annotation"} !~ /GPL15190/){ next; }
	if($line=~/^type_matrix/){
		if($option && $option !=2){ next; }
		my $file_matrix = $hash{"type_matrix"};
		push(@matrix_list,["public-anova-".$file_matrix,$descr,$org,$type]);
	}elsif($line=~/^type_geneset/){
		if($option && $option !=3){ next; }
		my $file_geneset = $hash{"type_geneset"};
		#if($file_geneset =~ /^(GO|MSIGDB)_geneset|^Protein_interaction|^(KEGG|PID|REACTOME|BIOCARTA)_pathway|^MIR_targets|^GNF2_coregulated/){ next; }
		push(@geneset_list,["public-".$file_geneset,$descr,$org]);
	}elsif($line=~/^type_annotation/){
		if($option && $option !=1){ next; }
		my $file_annotation = $hash{"type_annotation"};
		$file_annotation =~ s/_annot.txt$//;
		if($hashPlatform->{$file_annotation} && $hashPlatform->{$file_annotation}->[1] && $hashPlatform->{$file_annotation}->[1] ne $org){
			error_message("Wrong organism ID $org for $file_annotation");
		}
		#Do not process array for miRNA:
		if($org>=9000 && $org<=12000 && $type eq "mRNA"){
			push(@annotation_list,["public-".$file_annotation."_annot.txt",$descr,$org]);
		}
	}
}
close INFO;
if($logfile){ print OUTLOG "File name\tFile_type\tOrganimsID\tNsymbols\tNchanged\tNadded\n"; }
foreach my $ref (sort {$a->[2] cmp $b->[2]} @annotation_list,@matrix_list){
	update_symbols_in_arrays($ref->[0],$ref->[2]);
}
foreach my $ref (sort {$a->[2] cmp $b->[2]} @geneset_list){
	update_symbols_in_geneset($ref->[0],$ref->[2]);
}
return;
}


