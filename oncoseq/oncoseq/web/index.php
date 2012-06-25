<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
  <meta http-equiv="content-type" content="text/html; charset=utf-8" />
  <link rel="shortcut icon" type="image/ico" href="favicon.ico" />
  <title>vtable</title>
  <style type="text/css" title="currentStyle">
    @import "style/style.css";
  </style>
  <script src="jquery/js/jquery.js"></script>
  <script src="jquery/js/jquery.dataTables.js"></script>
  <script src="javascript/row_expand.js"></script>
  <script src="javascript/footer_filter.js"></script>
  <script type="text/javascript">
  $(document).ready(function() {
    oTable = $('#vtable').dataTable({
      "bJQueryUI": true,
      "bScrollInfinite": true,
      "bScrollCollapse": false,
      "iDisplayLength": 40,
      "sScrollY": "200px",
      "sScrollX": "200px",
      "sAjaxSource": "data.php",
      "bDeferRender": true,
      "fnDrawCallback": fnDrawCallback,
      "aoColumnDefs": [{"bSortable": false, "aTargets": [0]}],
      "aaSorting": [[1, 'asc']]
    });
    
    listenToDetails(); // listener for details expansion 
    listenToFooter();  // listener for footer clicking 
  });
  </script>
</head>
<body id="page_body" class="ex_highlight_row">
  <div id="page_body_container" style="width:100%">    
    <div id="table_container">
      <table cellpadding="0" cellspacing="0" border="0" class="display" id="vtable">
      <?php
        $IFP = fopen("/Library/WebServer/Documents/flamingskull/file","r") or die("Error: file not found.");
        
        //header line
        $line = fgets($IFP);
        $htokens = explode("\t",chop($line));
        echo "<thead><tr><th></th>";
        array_map(create_function('$n', 'echo "<th style=\"word-wrap:break-word\">$n</th>";'),$htokens);
        echo "</tr></thead><tfoot><th></th>\n";
        
        $i = 2;
        //footer
        foreach($htokens as $n)
        {
          echo "<th><input type=\"text\" name=\"";
          echo $n;
          echo "\" value=\"Search ";
          echo $n;
          echo "\" tabindex=$i class=\"search_init\"/></th>"; //search_init_numeric
          $i++;
        }
        echo "</tr></tfoot>";
		?>
      </table>
    </div>
    </div>
    <!--vtable built by Brendan Veeneman-->
    <!--DataTables designed and created by <a href="http://www.sprymedia.co.uk">Allan Jardine</a> &copy; 2007-2011<br> -->
</body>
</html>