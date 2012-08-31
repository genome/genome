prepareDataTables = (data) ->
    columns = data.aoColumns
    for col in columns
        if col.mDataProp == "command class"
            col.fnRender = (obj) ->
                cellVal = obj.aData[obj.iDataColumn]
                "<a href='/view/genome/task/status.html?id=#{data.wutgiTaskIds[obj.iDataRow]}'>#{cellVal}</a>"
        if col.mDataProp == "status"
            col.fnRender = (obj)->
                cellVal = obj.aData[obj.iDataColumn]
                "<span class='task-status #{cellVal}'>#{cellVal}</span>"
        val = col.mDataProp
        $("#table-header").append($("<th>#{val}</th>"))
    $("#loading-task-info").hide()
    $("#task-list").dataTable({aaData: data.aaData, aoColumns: data.aoColumns, bJQueryUI: true, bFilter:false})

$ ()->
    $.ajax
        url: '/view/genome/task/set/data-table.json'
        dataType: 'json'
        success: (data) ->
            prepareDataTables(data)
            $('#task-list_length').parent().prepend($('<h3>Upload Tasks</h3>'))
