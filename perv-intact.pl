#!/usr/bin/perl -w 

use strict;

my $file = shift;

####
# variables
####
my $mdelg = 30; # tolerant size of deletion in gag;
my $mdelp = 30; # tolerant size of deletion in pol;
my $mdele = 30; # tolerant size of deletion in env;

my $stpg = 1; # tolerant shift of stop codon in gag;
my $stpp = 2; # tolerant shift of stop codon in pol;
my $stpe = 7; # tolerant shift of stop codon in env;
####



my %list;
$list{gag}= "585-2159";
#$list{pol}= "2307-5744";
$list{pol}= "2307-5747";
$list{env}= "5620-7536";

my $list;
foreach my $i (keys %list){
  my @tmp = split /\-/, $list{$i};
  $list->{$i}->{s} = $tmp[0];
  $list->{$i}->{l} = $tmp[1];
}

my $seq;
my $n = 0;

if($file =~ /.gz$/){
  open FILE, "gunzip -cd $file";
}else{
  open FILE, "$file";
}
while(<FILE>){
  s/\A\s+//;
  s/\s+\Z//;
  my $line = $_;

  if(/^\>([\W\w]+)$/){
    my $name = $1;
    $n++;
    $name =~ s/\,//g;
    $seq->{$n}->{n} = $name;
    $seq->{$n}->{s} = "";
  }else{
    $seq->{$n}->{s} = $seq->{$n}->{s} . uc($line); 
  }
}
close FILE;

printf "#sequence ID,,";
printf "length(gag),length(pol),length(env),";

printf "no-start(gag),no-start(pol),no-start(env),";
printf "frameshift(gag),frameshift(pol),frameshift(env),";
printf "psc(gag),psc(pol),psc(env),";
printf "no-stop(gag),no-stop(pol),no-stop(env),";
printf "large-del(gag),large-del(pol),large-del(env),";

printf "#ins(gag),#ins(pol),#ins(env),";
printf "#del(gag),#del(pol),#del(env),";

printf "AA sequence(Gag),";
printf "AA sequence(Pol),";
printf "AA sequence(Env),";
print "\n";

my @ref1 = split //, $seq->{1}->{s};
my @ref2 = split //, $seq->{2}->{s};
my @ref3 = split //, $seq->{3}->{s};
for(my $i = 1; $i <= $n; $i++){
  my %len;
  my %delC;
  my %insC;
  my %diffC;
  my %delB;
  my %insB;
  my %diffB;
  my %delA;
  my %insA;
  my %diffA;
  my %ss1;
  my %ss2;
  my @tmp = split //, $seq->{$i}->{s};
  my $j = 0;
  for(my $i2 = 0; $i2 <= $#ref1; $i2++){
    if($ref1[$i2] =~ /[A-Z]/){
      $j++;
    }

    foreach my $id ("gag", "pol", "env"){
      if(! exists $len{$id}){ $len{$id} = 0; }
      if(! exists $ss1{$id}){ $ss1{$id} = ""; }
      if(! exists $ss2{$id}){ $ss2{$id} = ""; }
      if(! exists $insC{$id}){ $insC{$id} = 0; }
      if(! exists $insB{$id}){ $insB{$id} = 0; }
      if(! exists $insA{$id}){ $insA{$id} = 0; }
      if(! exists $delC{$id}){ $delC{$id} = 0; }
      if(! exists $delB{$id}){ $delB{$id} = 0; }
      if(! exists $delA{$id}){ $delA{$id} = 0; }
      if(! exists $diffC{$id}){ $diffC{$id} = 0; }
      if(! exists $diffB{$id}){ $diffB{$id} = 0; }
      if(! exists $diffA{$id}){ $diffA{$id} = 0; }

      if($j >= $list->{$id}->{s} && $j <= $list->{$id}->{l}){
        if($tmp[$i2] =~ /[A-Z]/){
          $len{$id} = $len{$id} + 1;
          $ss1{$id} = $ss1{$id} . $tmp[$i2];

          if(!($ref1[$i2] =~ /[A-Z]/)){ $insC{$id} = $insC{$id} + 1; }
          if(!($ref2[$i2] =~ /[A-Z]/)){ $insB{$id} = $insB{$id} + 1; }
          if(!($ref3[$i2] =~ /[A-Z]/)){ $insA{$id} = $insA{$id} + 1; }
        }else{
          if($ref1[$i2] =~ /[A-Z]/){ $delC{$id} = $delC{$id} + 1; }
          if($ref2[$i2] =~ /[A-Z]/){ $delB{$id} = $delB{$id} + 1; }
          if($ref3[$i2] =~ /[A-Z]/){ $delA{$id} = $delA{$id} + 1; }
        }
        $ss2{$id} = $ss2{$id} . $tmp[$i2];

        if($ref1[$i2] ne $tmp[$i2]){ $diffC{$id} = $diffC{$id} + 1; }
        if($ref2[$i2] ne $tmp[$i2]){ $diffB{$id} = $diffB{$id} + 1; }
        if($ref3[$i2] ne $tmp[$i2]){ $diffA{$id} = $diffA{$id} + 1; }

      }
    }
  }

  my %aa;
  my %psc;
  my %nss;
  my %nsl;
  foreach my $id ("gag", "pol", "env"){
    $aa{$id} = &na2aa($ss1{$id});
    my $m = 0;
    if($id eq "gag"){ $m = $stpg; }
    if($id eq "pol"){ $m = $stpp; }
    if($id eq "env"){ $m = $stpe; }

    if($aa{$id} =~ /^([\W\w]+)([\W\w]{$m})$/){
      my $s1 = $1;
      my $s2 = $2;
      if($s1 =~ /\*/){ $psc{$id} = 1; }
      else{ $psc{$id} = 0; }

      if(!($s2 =~ /\*/)){ $nsl{$id} = 1; }
      else{ $nsl{$id} = 0; }
    }
  
    if($ss2{$id} =~ /^ATG/){ $nss{$id} = 0; }
    else{ $nss{$id} = 1; }
  }

  printf "%s,", $seq->{$i}->{n};

  my $nnC = $diffC{gag} + $diffC{pol} + $diffC{env};
  my $nnB = $diffB{gag} + $diffB{pol} + $diffB{env};
  my $nnA = $diffA{gag} + $diffA{pol} + $diffA{env};

  my %ins;
  my %del;
  if($nnC > $nnB){ 
    if($nnA > $nnB){ 
      print "B-like,";
      $ins{gag} = $insB{gag};
      $ins{pol} = $insB{pol};
      $ins{env} = $insB{env};
      $del{gag} = $delB{gag};
      $del{pol} = $delB{pol};
      $del{env} = $delB{env};
    #}elsif($nnB >= $nnA){ 
    }else{
      print "A(or C)-like,";
      $ins{gag} = $insA{gag};
      $ins{pol} = $insA{pol};
      $ins{env} = $insA{env};
      $del{gag} = $delA{gag};
      $del{pol} = $delA{pol};
      $del{env} = $delA{env};
    }
#  }elsif($nnB >= $nnC){ 
  }else{
    if($nnA > $nnC){
      print "C(or A)-like,";
      $ins{gag} = $insC{gag};
      $ins{pol} = $insC{pol};
      $ins{env} = $insC{env};
      $del{gag} = $delC{gag};
      $del{pol} = $delC{pol};
      $del{env} = $delC{env};
    #}elsif($nnC >= $nnA){ 
    }else{
      print "A(or C)-like,";
      $ins{gag} = $insA{gag};
      $ins{pol} = $insA{pol};
      $ins{env} = $insA{env};
      $del{gag} = $delA{gag};
      $del{pol} = $delA{pol};
      $del{env} = $delA{env};
    }
  }

  printf "%d,%d,%d,", $len{gag}, $len{pol}, $len{env};

  if($nss{gag} == 1){ print "+,"; }
  else{ print "-,"; }
  if($nss{pol} == 1){ print "+,"; }
  else{ print "-,"; }
  if($nss{env} == 1){ print "+,"; }
  else{ print "-,"; }

  my $fs = abs($ins{gag} - $del{gag});
  if(($fs % 3) == 0){ print "-,"; }
  else{               print "+,"; }
  $fs = abs($ins{pol} - $del{pol});
  if(($fs % 3) == 0){ print "-,"; }
  else{               print "+,"; }
  $fs = abs($ins{env} - $del{env});
  if(($fs % 3) == 0){ print "-,"; }
  else{               print "+,"; }
  
  if($psc{gag} == 1){ print "+,"; }
  else{ print "-,"; }
  if($psc{pol} == 1){ print "+,"; }
  else{ print "-,"; }
  if($psc{env} == 1){ print "+,"; }
  else{ print "-,"; }

  if($nsl{gag} == 1){ print "+,"; }
  else{ print "-,"; }
  if($nsl{pol} == 1){ print "+,"; }
  else{ print "-,"; }
  if($nsl{env} == 1){ print "+,"; }
  else{ print "-,"; }

  if($del{gag} > $mdelg){ print "+,"; }
  else{ print "-,"; }
  if($del{pol} > $mdelp){ print "+,"; }
  else{ print "-,"; }
  if($del{env} > $mdele){ print "+,"; }
  else{ print "-,"; }

  printf "%d,", $ins{gag};
  printf "%d,", $ins{pol};
  printf "%d,", $ins{env};

  printf "%d,", $del{gag};
  printf "%d,", $del{pol};
  printf "%d,", $del{env};

  printf "%s,", &na2aa($ss1{gag});
  printf "%s,", &na2aa($ss1{pol});
  printf "%s,", &na2aa($ss1{env});
  print "\n";
}



sub na2aa
{
  my $s = shift;

  my @s = split //, lc ($s);

  my $aa = "";

  for(my $i = 0; $i <= $#s; $i = $i + 3){
    my $na = $s[$i];
    if(($i+1) <= $#s){ 
       $na = $na . $s[$i+1];
    }
    if(($i+2) <= $#s){ 
       $na = $na . $s[$i+2];
    }

    my $flag = 0;
    if($na eq "ttt"){ $aa = $aa . "F"; $flag++; }
    if($na eq "ttc"){ $aa = $aa . "F"; $flag++; }
    if($na eq "tta"){ $aa = $aa . "L"; $flag++; }
    if($na eq "ttg"){ $aa = $aa . "L"; $flag++; }

    if($na eq "tct"){ $aa = $aa . "S"; $flag++; }
    if($na eq "tcc"){ $aa = $aa . "S"; $flag++; }
    if($na eq "tca"){ $aa = $aa . "S"; $flag++; }
    if($na eq "tcg"){ $aa = $aa . "S"; $flag++; }

    if($na eq "tat"){ $aa = $aa . "Y"; $flag++; }
    if($na eq "tac"){ $aa = $aa . "Y"; $flag++; }
    if($na eq "taa"){ $aa = $aa . "*"; $flag++; }
    if($na eq "tag"){ $aa = $aa . "*"; $flag++; }

    if($na eq "tgt"){ $aa = $aa . "C"; $flag++; }
    if($na eq "tgc"){ $aa = $aa . "C"; $flag++; }
    if($na eq "tga"){ $aa = $aa . "*"; $flag++; }
    if($na eq "tgg"){ $aa = $aa . "W"; $flag++; }


    if($na eq "ctt"){ $aa = $aa . "L"; $flag++; }
    if($na eq "ctc"){ $aa = $aa . "L"; $flag++; }
    if($na eq "cta"){ $aa = $aa . "L"; $flag++; }
    if($na eq "ctg"){ $aa = $aa . "L"; $flag++; }

    if($na eq "cct"){ $aa = $aa . "P"; $flag++; }
    if($na eq "ccc"){ $aa = $aa . "P"; $flag++; }
    if($na eq "cca"){ $aa = $aa . "P"; $flag++; }
    if($na eq "ccg"){ $aa = $aa . "P"; $flag++; }

    if($na eq "cat"){ $aa = $aa . "H"; $flag++; }
    if($na eq "cac"){ $aa = $aa . "H"; $flag++; }
    if($na eq "caa"){ $aa = $aa . "Q"; $flag++; }
    if($na eq "cag"){ $aa = $aa . "Q"; $flag++; }

    if($na eq "cgt"){ $aa = $aa . "R"; $flag++; }
    if($na eq "cgc"){ $aa = $aa . "R"; $flag++; }
    if($na eq "cga"){ $aa = $aa . "R"; $flag++; }
    if($na eq "cgg"){ $aa = $aa . "R"; $flag++; }


    if($na eq "att"){ $aa = $aa . "I"; $flag++; }
    if($na eq "atc"){ $aa = $aa . "I"; $flag++; }
    if($na eq "ata"){ $aa = $aa . "I"; $flag++; }
    if($na eq "atg"){ $aa = $aa . "M"; $flag++; }

    if($na eq "act"){ $aa = $aa . "T"; $flag++; }
    if($na eq "acc"){ $aa = $aa . "T"; $flag++; }
    if($na eq "aca"){ $aa = $aa . "T"; $flag++; }
    if($na eq "acg"){ $aa = $aa . "T"; $flag++; }

    if($na eq "aat"){ $aa = $aa . "N"; $flag++; }
    if($na eq "aac"){ $aa = $aa . "N"; $flag++; }
    if($na eq "aaa"){ $aa = $aa . "K"; $flag++; }
    if($na eq "aag"){ $aa = $aa . "K"; $flag++; }

    if($na eq "agt"){ $aa = $aa . "S"; $flag++; }
    if($na eq "agc"){ $aa = $aa . "S"; $flag++; }
    if($na eq "aga"){ $aa = $aa . "R"; $flag++; }
    if($na eq "agg"){ $aa = $aa . "R"; $flag++; }


    if($na eq "gtt"){ $aa = $aa . "V"; $flag++; }
    if($na eq "gtc"){ $aa = $aa . "V"; $flag++; }
    if($na eq "gta"){ $aa = $aa . "V"; $flag++; }
    if($na eq "gtg"){ $aa = $aa . "V"; $flag++; }

    if($na eq "gct"){ $aa = $aa . "A"; $flag++; }
    if($na eq "gcc"){ $aa = $aa . "A"; $flag++; }
    if($na eq "gca"){ $aa = $aa . "A"; $flag++; }
    if($na eq "gcg"){ $aa = $aa . "A"; $flag++; }

    if($na eq "gat"){ $aa = $aa . "D"; $flag++; }
    if($na eq "gac"){ $aa = $aa . "D"; $flag++; }
    if($na eq "gaa"){ $aa = $aa . "E"; $flag++; }
    if($na eq "gag"){ $aa = $aa . "E"; $flag++; }

    if($na eq "ggt"){ $aa = $aa . "G"; $flag++; }
    if($na eq "ggc"){ $aa = $aa . "G"; $flag++; }
    if($na eq "gga"){ $aa = $aa . "G"; $flag++; }
    if($na eq "ggg"){ $aa = $aa . "G"; $flag++; }

    if($na eq "---"){ $aa = $aa . "-"; $flag++; }

    if($flag == 0){ $aa = $aa . "X"; }
  }

  return ($aa);
}

