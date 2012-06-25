// Formating function for row details
function fnFormatDetails ( oTable, nTr )
{
  var aData = oTable.fnGetData( nTr );
  var sOut = '<table cellpadding="3" cellspacing="0" style="padding-left:10px;">';
  sOut += '<tr><td>Rendering engine:</td><td>'+aData[1]+' '+aData[5]+'</td></tr>';
  sOut += '<tr><td>Link to source:</td><td>Could provide a link here</td></tr>';
  sOut += '<tr><td>Extra info:</td><td>And any further details here (images etc)</td></tr>';
  sOut += '</table>';
  return sOut;
}

function fnDrawCallback ()
{
  var tableWrapper = $('#vtable_wrapper');
  tableWrapper.find('.dataTables_scrollBody').css('height', $(document).height() - 160); //303px
  tableWrapper.css('width', '');
  var tableDataContent = tableWrapper.find('.dataTables_scrollBody')[0];
  var tableHasScrollBar = (tableDataContent.scrollHeight > tableDataContent.clientHeight);
  //fixes header aligment issues in all major browsers
  if (tableHasScrollBar) {
    tableWrapper.find('.dataTables_scrollHeadInner').css('margin-right', '15px');
    tableWrapper.find('.dataTables_scrollFootInner').css('margin-right', '15px');
  }
  //fixes table data stretchyness in <IE8
  if (navigator.appVersion.indexOf("MSIE 6") != -1 || navigator.appVersion.indexOf("MSIE 7") != -1) {
    $('#' + tableId).css('width', '');
    tableWrapper.find('.dataTables_scrollBody').css('overflow-x', 'hidden');
  }
}

// Event listener for opening and closing details
function listenToDetails ()
{
  $('#vtable tbody td img').live('click', function () {
    var nTr = $(this).parents('tr')[0];
    if ( oTable.fnIsOpen(nTr) ) // This row is already open - close it
    {
      this.src = "./style/images/details_open.png";
      oTable.fnClose( nTr );
    }
    else // Open this row
    {
      this.src = "./style/images/details_close.png";
      oTable.fnOpen( nTr, fnFormatDetails(oTable, nTr), 'details' );
    }
  });
}