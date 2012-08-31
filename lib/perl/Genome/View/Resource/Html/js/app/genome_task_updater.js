(function() {
  $(function() {
    var get_url_param, json_url;
    get_url_param = function(name) {
      var results;
      results = new RegExp('[\\?&]' + name + '=([^&#]*)').exec(window.location.href);
      if (!results) {
        return 0;
      } else {
        return results[1] || 0;
      }
    };
    json_url = "/view/genome/task/output.json?id=" + (get_url_param('id'));
    window.do_update = function() {
      return $.ajax({
        url: json_url,
        type: 'GET',
        success: function(data, textStatus, jqXHR) {
          $('#task-status').text(data.status);
          if (data.status !== 'submitted') {
            $('#stdout-content').text(data.stdout);
            $('#stderr-content').text(data.stderr);
            $('#time-started').text(data.time_started);
            $('#time-finished').text(data.time_started);
          }
          if (data.status === 'succeeded' || data.status === 'failed') {
            return $('.spinner').hide();
          } else {
            $('.spinner').show();
            return setTimeout("do_update()", 10000);
          }
        }
      });
    };
    $('.spinner').show();
    return do_update();
  });
}).call(this);
