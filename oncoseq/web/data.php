<?php
  header("content-type: text/plain");
  header("Pragma: no-cache");
  header("Expires: 0");
  
  echo "{ \"aaData\": [\n";
  
  $IFP = fopen("/Library/WebServer/Documents/flamingskull/file","r") or die("Error: file not found.");

  fgets($IFP); //header
  
  $line = fgets($IFP); //line1 - need to in order to do commas
  $line = str_replace("\"","",$line);
  $tokens = explode("\t",chop($line)); //chomp and tokenize
  echo "[\"<img src=./style/images/details_open.png>\",\"" . implode("\",\"",$tokens) . "\"]";
  
  //table body
  while($line = fgets($IFP))
  {
    $line = htmlentities($line);
    $tokens = explode("\t",chop($line)); //chomp and tokenize
    echo ",\n[\"<img src=./style/images/details_open.png>\",\"" . implode("\",\"",$tokens) . "\"]";
  }
  feof($IFP) ? fclose($IFP) : die("Error: fgets failed.\n"); //if error, close or die
  echo "\n]}";
  
  return;
?>