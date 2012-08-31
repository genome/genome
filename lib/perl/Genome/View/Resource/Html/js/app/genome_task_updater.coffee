# updater

$ ->
    # adapted from this helpful answer
    # http://stackoverflow.com/questions/901115/get-query-string-values-in-javascript/2480180#2480180
    get_url_param = (name) ->
        results = new RegExp('[\\?&]' + name + '=([^&#]*)').exec(window.location.href)
        if (!results) 
            return 0
        else
            return results[1] || 0

    json_url = "/view/genome/task/output.json?id=#{get_url_param('id')}"

    window.do_update = ()->
        $.ajax
           url: json_url
           type: 'GET'
           success: (data, textStatus, jqXHR) ->
                $('#task-status').text(data.status)

                if (data.status != 'submitted')
                    $('#stdout-content').text(data.stdout)
                    $('#stderr-content').text(data.stderr)
                    $('#time-started').text(data.time_started)
                    $('#time-finished').text(data.time_started)

                if (data.status == 'succeeded' || data.status == 'failed')
                    $('.spinner').hide()
                else
                    $('.spinner').show()
                    setTimeout("do_update()", 10000);

    $('.spinner').show()
    do_update()
