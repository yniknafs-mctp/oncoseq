//Support functions to provide a little bit of 'user friendliness' to the textboxes in 
//the footer
function listenToFooter()
{ 
  var asInitVals = new Array();
  var numericCols = new Array();

  $("tfoot input").keyup( function () { // perform filtering    
    if(numericCols[$("tfoot input").index(this)] == 1 &&
       (this.value.substring(0,1) == ">" ||
        this.value.substring(0,1) == "<")) //numeric filter
    {
      oTable.fnDraw();
    }
    else //normal filter
    {
      oTable.fnFilter( this.value, $("tfoot input").index(this) + 1 );
    }
  });

  $("tfoot input").each( function (i) { //saves initial values for reloading
    asInitVals[i] = this.value;
    if ( this.className == "search_init_numeric" )
    {
      numericCols[i] = 1;
    }
    else
    {
      numericCols[i] = 0;
    }
  });
   
  $("tfoot input").focus( function () { //blanks the box for entry
    if ( this.className == "search_init" || this.className == "search_init_numeric")
    {
      this.className = "";
      this.value = "";
    }
  });
  
  $("tfoot input").blur( function (i) { //restores gray text when exited
    if ( this.value == "" )
    {
      this.value = asInitVals[$("tfoot input").index(this)];

      if(numericCols[$("tfoot input").index(this)] == 1)
      {
        this.className = "search_init_numeric";
      }
      else
      {
        this.className = "search_init";
      }
    }
  });
  
  //Fancy filtering by brendan veeneman XXX - this could be more efficient
  // nb - this sits on TOP of normal filtering
  $.fn.dataTableExt.afnFiltering.push(
    function( oSettings, aData, iDataIndex ) {
      
      for(i in numericCols)
      {
        if(numericCols[i] == 0){ continue; }
        
        var input_val = $("tfoot input")[i].value;
        var table_val = aData[i*1 + 1] == "-" ? 0 : aData[i*1 + 1];
                
        if ( input_val.substring(0,1) == "<" && table_val >= (input_val.substring(1) * 1)) /* lt */
        { return false; }
        else if ( input_val.substring(0,1) == ">" && table_val <= (input_val.substring(1) * 1)) /* gt */
        { return false; }
      }
      
      return true;
    }
  );
}