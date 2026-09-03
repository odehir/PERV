#!/usr/bin/perl -w

use strict;

####
# variables
####
my $threthold = 50; # threthold of probability (%);

my $dataM;
my $dataI;
my %lenI;
my %name;
my %omit;
my $id = 0;
my %list;
my %list2;
my $ll;


my $file = shift;

if($file =~ /.gz$/){
  open FILE, "gunzip -cd $file |";
}else{
  open FILE, "$file";
}
while(<FILE>){
  s/^\s+//;
  s/\s+$//;
  my $line = $_;

  if(/^\@/){ next; }
  if(/^$/){ next; }


  my @tmp = split /\t+/, $line;

  if($tmp[9] =~ /[Nn]/){ next; }
  if($tmp[9] =~ / /){ next; }

  my $flag  = $tmp[1];
  if($flag != 0 && $flag != 16){ next; }

  my $start = $tmp[3];

  my $cigar = $tmp[5];
  my $cigar1 = $tmp[5];
     $cigar1 =~ s/[A-Z]/ /g;
     $cigar1 =~ s/\s+$//;
  my $cigar2 = $tmp[5];
     $cigar2 =~ s/\d+/ /g;
     $cigar2 =~ s/^\s+//;
  my @cigar1 = split /\s+/, $cigar1;
  my @cigar2 = split /\s+/, $cigar2;
  

  $id++;
  if(!exists $list2{$tmp[0]}){ $list2{$tmp[0]} = 1; }
  else{ next; }

  my $n = $start - 1;
  my $n2 = 0;
  my $seq = $tmp[9];
     $seq = lc ($seq);
  my @seq = split //, $seq;
  my $tmp = "";
  for(my $i = 0; $i <= $#cigar2; $i++){
    if($cigar2[$i] eq "S"){
#      $n  = $n  + $cigar1[$i];
#      $n2 = $n2 + $cigar1[$i];
      for(my $j = 0; $j < $cigar1[$i]; $j++){
        if($n >= (1 - 1)){
          $tmp = $tmp . uc($seq[$n2]);
        }
        $n2++;
      }
    }
    if($cigar2[$i] eq "H"){
      ;
    }
    if($cigar2[$i] eq "M"){
      for(my $j = 0; $j < $cigar1[$i]; $j++){
        if($n >= (1 - 1)){
          $tmp = $tmp . $seq[$n2];
        }
        $n++;
        $n2++;
      }
    }
    if($cigar2[$i] eq "I"){
      for(my $j = 0; $j < $cigar1[$i]; $j++){
        if($n >= (1)){
          $tmp = $tmp . uc ($seq[$n2]);
        }
        $n2++;
      }
    }
    if($cigar2[$i] eq "D"){
      for(my $j = 0; $j < $cigar1[$i]; $j++){
        if($n >= (1 - 1)){
          if($seq[$n2]){
            $tmp = $tmp . "-";
          }
        }
        $n++;
      }
    }
  }

  if($tmp =~ /[Nn]/){ next; }

  my $mm = "";
  my $ml = "";
  for(my $j = 11; $j <= $#tmp; $j++){
    if($tmp[$j] =~ /^MM\:Z\:([\W\w]+)\;$/){
      $mm = $1;
    }
    if($tmp[$j] =~ /^ML\:B\:C\,([\W\w]+)$/){
      $ml = $1;
    }
  }
  my $mmp;
  my %mmpid;
  my $mmm;
  my %mmmid;
  my @mm = split /\;/, $mm;
  my $mmaxp = 0;
  my $mmaxm = 0;
  for(my $j = 0; $j <= $#mm; $j++){
    my @tmp = split /\,/, $mm[$j];
    my $id = shift @tmp;
    my @ml  = split /\,/, $ml;

    if($id =~ /^C\+/){
      $mmaxp++;
      $mmpid{$mmaxp} = $id;
      my $nn = 0;
      for(my $k = 0; $k <= $#tmp; $k++){
        $nn = $nn + $tmp[$k] + 1;
        $mmp->{$k}->{pos}->{$mmaxp} = $nn;
        $mmp->{$k}->{prob}->{$mmaxp} = $ml[($#tmp + 1)*$j + $k];
      }
    }
    if($id =~ /^C\-/){
      $mmaxm++;
      $mmmid{$mmaxm} = $id;
      my $nn = 0;
      for(my $k = 0; $k <= $#tmp; $k++){
        $nn = $nn + $tmp[$k] + 1;
        $mmm->{$k}->{pos}->{$mmaxm} = $nn;
        $mmm->{$k}->{prob}->{$mmaxm} = $ml[($#tmp + 1)*$j + $k];
      }
    }
  }
  
  if($flag == 0){
    my @tmp2 = split //, $tmp;
    my @tmp3p;
    my @tmp3m;
    for(my $j = 0; $j <= $#tmp2; $j++){
      if($tmp2[$j] =~ /[Cc]/){ push @tmp3p, $j; }
      if($tmp2[$j] =~ /[Gg]/){ push @tmp3m, $j; }
    }
  
    if($mmaxp > 0){
      foreach my $j (sort {$main::a <=> $main::b} keys %{$mmp}){
        my $tpos = $mmp->{$j}->{pos}->{1};
        my $pos  = $tmp3p[$tpos - 1];
        my $aa = $tmp2[$pos];

        my $pp = 0;
        foreach my $k (sort {$main::a <=> $main::b} keys %{$mmp->{$j}->{pos}}){
          my $id = "";
          if($mmpid{$k} =~ /^C\+([a-z])[\.\?]?$/){
            $id = $1;
          }
          my $tpp = $mmp->{$j}->{prob}->{$k};
          if(($tpp + 1)/256*100 >= $pp && ($tpp + 1)/256*100 >= $threthold){
            $pp = ($tpp + 1)/256*100;
            if($aa =~ /[a-z]/){ $aa = lc($id); }
            if($aa =~ /[A-Z]/){ $aa = uc($id); }
          }
        }

        $tmp2[$pos] = $aa;
      }
    }

    if($mmaxm > 0){
      foreach my $j (sort {$main::a <=> $main::b} keys %{$mmm}){
        my $tpos = $mmm->{$j}->{pos}->{1};
        my $pos  = $tmp3m[$tpos - 1];
        my $aa = $tmp2[$pos];

        my $pp = 0;
        foreach my $k (sort {$main::a <=> $main::b} keys %{$mmm->{$j}->{pos}}){
          my $id = "";
          if($mmmid{$k} =~ /^C\+([a-z])[\.\?]?$/){
            $id = $1;
          }
          my $tpp = $mmm->{$j}->{prob}->{$k};
          if(($tpp + 1)/256*100 >= $pp && ($tpp + 1)/256*100 >= $threthold){
            $pp = ($tpp + 1)/256*100;
            if($aa =~ /[a-z]/){ $aa = lc($id); }
            if($aa =~ /[A-Z]/){ $aa = uc($id); }
          }
        }

        $tmp2[$pos] = $aa;
      }
    }

    $tmp = "";
    for(my $j = 0; $j <= $#tmp2; $j++){
      $tmp = $tmp . $tmp2[$j];
    }
  }
  elsif($flag == 16){
    my @tmp2 = split //, $tmp;
    my @tmp3p;
    my @tmp3m;
    for(my $j = $#tmp2; $j >= 0; $j--){
      if($tmp2[$j] =~ /[Gg]/){ push @tmp3p, $j; }
      if($tmp2[$j] =~ /[Cc]/){ push @tmp3m, $j; }
    }
 
    if($mmaxp > 0){
      foreach my $j (sort {$main::a <=> $main::b} keys %{$mmp}){
        my $tpos = $mmp->{$j}->{pos}->{1};
        my $pos  = $tmp3p[$tpos - 1];
        my $aa = $tmp2[$pos];

        my $pp = 0;
        foreach my $k (sort {$main::a <=> $main::b} keys %{$mmp->{$j}->{pos}}){
          my $id = "";
          if($mmpid{$k} =~ /^C\+([a-z])[\.\?]?$/){
            $id = $1;
          }
          my $tpp = $mmp->{$j}->{prob}->{$k};
          if(($tpp + 1)/256*100 >= $pp && ($tpp + 1)/256*100 >= $threthold){
            $pp = ($tpp + 1)/256*100;
            if($aa =~ /[a-z]/){ $aa = lc($id); }
            if($aa =~ /[A-Z]/){ $aa = uc($id); }
          }
        }

        $tmp2[$pos] = $aa;
      }
    }

    if($mmaxm > 0){
      foreach my $j (sort {$main::a <=> $main::b} keys %{$mmm}){
        my $tpos = $mmm->{$j}->{pos}->{1};
        my $pos  = $tmp3m[$tpos - 1];
        my $aa = $tmp2[$pos];
  
        my $pp = 0;
        foreach my $k (sort {$main::a <=> $main::b} keys %{$mmm->{$j}->{pos}}){
          my $id = "";
          if($mmmid{$k} =~ /^C\+([a-z])[\.\?]?$/){
            $id = $1;
          }
          my $tpp = $mmm->{$j}->{prob}->{$k};
          if(($tpp + 1)/256*100 >= $pp && ($tpp + 1)/256*100 >= $threthold){
            $pp = ($tpp + 1)/256*100;
            if($aa =~ /[a-z]/){ $aa = lc($id); }
            if($aa =~ /[A-Z]/){ $aa = uc($id); }
          }
        }

        $tmp2[$pos] = $aa;
      }
    }

    $tmp = "";
    for(my $j = 0; $j <= $#tmp2; $j++){
      $tmp = $tmp . $tmp2[$j];
    }
  }

  if($tmp =~ /^[A-Z]+([a-z][\W\w]+)$/){ $tmp = $1; }
  if($tmp =~ /^([\W\w]+[a-z])[A-Z]+$/){ $tmp = $1; }
#  $tmp =~ s/[A-Z]//g;

  if($tmp =~ /^$/ || $tmp =~ /^\-+$/){ ; }
  else{
    $name{$id} = $tmp[0];
    printf ">%s\n", $tmp[0] . ";" . $tmp[2];
    printf "%s", $tmp;
    print "\n";
  }
}
close FILE;
