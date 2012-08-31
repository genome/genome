(function() {
  var prepareDataTables;
  prepareDataTables = function(data) {
    var col, columns, val, _i, _len;
    columns = data.aoColumns;
    for (_i = 0, _len = columns.length; _i < _len; _i++) {
      col = columns[_i];
      if (col.mDataProp === "command class") {
        col.fnRender = function(obj) {
          var cellVal;
          cellVal = obj.aData[obj.iDataColumn];
          return "<a href='/view/genome/task/status.html?id=" + data.wutgiTaskIds[obj.iDataRow] + "'>" + cellVal + "</a>";
        };
      }
      if (col.mDataProp === "status") {
        col.fnRender = function(obj) {
          var cellVal;
          cellVal = obj.aData[obj.iDataColumn];
          return "<span class='task-status " + cellVal + "'>" + cellVal + "</span>";
        };
      }
      val = col.mDataProp;
      $("#table-header").append($("<th>" + val + "</th>"));
    }
    $("#loading-task-info").hide();
    return $("#task-list").dataTable({
      "oLanguage": { "sUrl": "/res/tgi.txt" },
      aaData: data.aaData,
      aoColumns: data.aoColumns,
      bJQueryUI: true,
      bFilter: false
    });
  };
  $(function() {
    return $.ajax({
      url: '/view/genome/task/set/data-table.json',
      dataType: 'json',
      success: function(data) {
        prepareDataTables(data);
        return $('#task-list_length').parent().prepend($('<h3>Upload Tasks</h3>'));
      }
    });
  });
}).call(this);
