(function() {
  var __hasProp = Object.prototype.hasOwnProperty, __extends = function(child, parent) {
    for (var key in parent) { if (__hasProp.call(parent, key)) child[key] = parent[key]; }
    function ctor() { this.constructor = child; }
    ctor.prototype = parent.prototype;
    child.prototype = new ctor;
    child.__super__ = parent.prototype;
    return child;
  };
  window.action_mode = 'PUT';
  $(function() {
    var Column, ColumnView, List, ListView, construct_from_json, fieldCount, listView, m, templateContents, url;
    templateContents = '<script type="text/template" id="nomenclature-column-template">\n  <label>Title:</label>\n  <input type="text" class="nomenclature-column-name text" value="<%= model.get(\'name\') %>">\n  <label>Type:</label>\n  <select class="column-type-select">\n    <option value="integer">integer</option>\n    <option value="real">real</option>\n    <option value="numeric">numeric (real or integer)</option>\n    <option value="date">date</option>\n    <option value="string">string</option>\n    <option value="enumerated">enumerated</option>\n  </select>\n  <button class="delete">Remove Column</button>\n</script>\n\n<script type="text/template" id="nomenclature-column-enumerated-template">\n  <div class="enumerated-choices">\n    <strong>Enumerated Column Definition:</strong><br/>\n    <% _.each(model.get("enumerated_values"), function(e) { %>\n\n    <div class="enum-row"><input type="text" class="enumerated-value-entry text" value="<%= e %>">\n\n    <button class="remove-enum">remove</button>\n    <button class="add-new-enum">add</button>\n  </div>\n\n  <% }); %>\n</div>\n</script>';
    $("#script-container").append($(templateContents));
    Backbone.sync = function(method, model, callback_bundle) {
      return callback_bundle.success();
    };
    fieldCount = 0;
    Column = (function() {
      __extends(Column, Backbone.Model);
      function Column() {
        Column.__super__.constructor.apply(this, arguments);
      }
      Column.prototype.defaults = {
        name: 'column',
        type: 'string',
        id: '',
        use_count: 0,
        enumerated_values: [],
        enumerated_value_ids: [],
        enumerated_value_use_counts: []
      };
      return Column;
    })();
    List = (function() {
      __extends(List, Backbone.Collection);
      function List() {
        List.__super__.constructor.apply(this, arguments);
      }
      List.prototype.model = Column;
      return List;
    })();
    ColumnView = (function() {
      __extends(ColumnView, Backbone.View);
      function ColumnView() {
        ColumnView.__super__.constructor.apply(this, arguments);
      }
      ColumnView.prototype.tagName = 'li';
      ColumnView.prototype.events = {
        'click .delete': 'remove',
        'blur input.nomenclature-column-name': 'nameChanged',
        'change select.column-type-select': 'typeChanged',
        'click .add-new-enum': 'addNewEnum',
        'click .remove-enum': 'removeEnum',
        'blur input.enumerated-value-entry': 'enumChanged'
      };
      ColumnView.prototype.nameChanged = function() {
        var new_name;
        new_name = $(this.el).children('.nomenclature-column-name')[0].value;
        return this.model.set({
          "name": new_name
        });
      };
      ColumnView.prototype.typeChanged = function() {
        var type;
        type = $(this.el).children('.column-type-select')[0].value;
        if (type === 'enumerated') {
          this.model.set({
            "enumerated_values": ['Untitled']
          }, {
            "silent": true
          });
        }
        return this.model.set({
          "type": type
        });
      };
      ColumnView.prototype.addNewEnum = function() {
        var types;
        types = this.model.get("enumerated_values");
        types.push("untitled");
        this.model.set({
          "enumerated_values": void 0
        }, {
          "silent": true
        });
        return this.model.set({
          "enumerated_values": types
        });
      };
      ColumnView.prototype.getEnumIndexForInitiator = function(initiator) {
        var entire_set, index, parent_row, _ref;
        parent_row = initiator.parent(".enum-row")[0];
        entire_set = $($(parent_row).parent('.enumerated-choices')).children('.enum-row');
        for (index = 0, _ref = entire_set.length; 0 <= _ref ? index < _ref : index > _ref; 0 <= _ref ? index++ : index--) {
          if (entire_set.get(index) === parent_row) {
            return index;
          }
        }
      };
      ColumnView.prototype.enumChanged = function(k) {
        var index, initiator, types, value;
        initiator = $(k.target);
        index = this.getEnumIndexForInitiator(initiator);
        value = initiator.val();
        types = this.model.get("enumerated_values");
        types[index] = value;
        this.model.set({
          "enumerated_values": void 0
        }, {
          silent: true
        });
        return this.model.set({
          "enumerated_values": types
        });
      };
      ColumnView.prototype.removeEnum = function(k) {
        var ids, index, initiator, types;
        initiator = $(k.target);
        index = this.getEnumIndexForInitiator(initiator);
        types = this.model.get("enumerated_values");
        ids = this.model.get("enumerated_value_ids");
        types.splice(index, 1);
        ids.splice(index, 1);
        this.model.set({
          "enumerated_values": void 0
        }, {
          silent: true
        });
        return this.model.set({
          "enumerated_values": types
        });
      };
      ColumnView.prototype.initialize = function() {
        _.bindAll(this, 'render', 'unrender', 'remove');
        this.model.bind('change', this.render);
        return this.model.bind('remove', this.unrender);
      };
      ColumnView.prototype.render = function() {
        var enumerated_template, i, item_template, kids, _ref;
        item_template = _.template($('#nomenclature-column-template').html(), {
          model: this.model
        });
        i = $(item_template);
        $(this.el).html(i);
        $(this.el).children(".column-type-select")[0].value = this.model.get('type');
        if (this.model.get('type') === 'enumerated') {
          enumerated_template = _.template($('#nomenclature-column-enumerated-template').html(), {
            model: this.model
          });
          $(this.el).append(enumerated_template);
          kids = $(this.el).find(".enum-row");
          for (i = 0, _ref = kids.size() - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
            if (parseInt(this.model.get('enumerated_value_use_counts')[i]) > 0) {
              $(kids[i]).children('.remove-enum').attr('disabled', 'true');
              $(kids[i]).children('.enumerated-value-entry').attr('disabled', 'true');
              $('#why-cant-i-delete').show();
            }
          }
        }
        if (parseInt(this.model.get('use_count')) > 0) {
          $(this.el).children(".delete").attr('disabled', 'true');
          $(this.el).children(".column-type-select").attr('disabled', 'true');
          $(this.el).children(".nomenclature-column-name").attr('disabled', 'true');
          $('#why-cant-i-delete').show();
        }
        return this;
      };
      ColumnView.prototype.unrender = function() {
        return $(this.el).remove();
      };
      ColumnView.prototype.remove = function() {
        return this.model.destroy();
      };
      return ColumnView;
    })();
    ListView = (function() {
      __extends(ListView, Backbone.View);
      function ListView() {
        ListView.__super__.constructor.apply(this, arguments);
      }
      ListView.prototype.el = $('#nomenclature-columns');
      ListView.prototype.events = {
        'click button#add': 'addColumn'
      };
      ListView.prototype.initialize = function() {
        _.bindAll(this, 'render', 'addColumn', 'appendColumn');
        this.collection = new List;
        this.collection.bind('add', this.appendColumn);
        this.collection.bind('reset', this.render);
        return this.render();
      };
      ListView.prototype.appendColumn = function(column) {
        var columnView;
        columnView = new ColumnView({
          model: column
        });
        return $('ul#nomenclature-list', this.el).append(columnView.render().el);
      };
      ListView.prototype.render = function() {
        $(this.el).html('');
        $(this.el).append("<ul id='nomenclature-list'></ul>");
        $(this.el).append("<div id='nomenclature-add'><button id='add'>Add column</button></div>");
        return _(this.collection.models).each(function(column) {
          return appendColumn(column, this);
        });
      };
      ListView.prototype.addColumn = function() {
        var column, name;
        fieldCount++;
        column = new Column;
        name = column.get('name');
        name = name + fieldCount;
        column.set({
          "name": "name",
          name: name
        });
        return this.collection.add(column);
      };
      return ListView;
    })();
    listView = new ListView;
    if (window.location.hash) {
      m = window.location.hash.split('=');
      if (m[0] === '#id') {
        url = '/view/Genome/Nomenclature/detail.json?id=' + m[1];
        $.ajax({
          url: url,
          type: 'GET',
          success: function(data, textStatus, jqXHR) {
            return construct_from_json(data);
          }
        });
      }
    }
    construct_from_json = function(object) {
      var name;
      window.action_mode = 'POST';
      name = object.name;
      $("#nomenclature-name-input").val(name);
      $("#nomenclature-empty-equivalent-input").val(object.empty_equivalent);
      window.nomenclature_id = object.id;
      $(".title h1").html("Edit Nomenclature: " + name);
      $("#directions").html("Use the form below to edit the nomenclature " + name);
      document.title = "Edit Nomenclature: " + name;
      _(object.fields).each(function(i) {
        return listView.collection.add(i);
      });
      return listView.render;
    };
    $('.save-nomenclature').bind('click', function() {
      var ajax_data, empty_equiv, jsonToPost, k, name;
      name = $("#nomenclature-name-input").val();
      if (name === "") {
        alert("You can't save a nomenclature without a name!");
        return;
      }
      empty_equiv = $("#nomenclature-empty-equivalent-input").val();
      if (empty_equiv === "") {
        alert("Please specify an empty-equivalent string!");
        return;
      }
      if (listView.collection.length === 0) {
        alert("You need at least one column to create a nomenclature!");
        return;
      }
      m = {
        "name": name,
        "empty_equivalent": empty_equiv,
        "fields": listView.collection
      };
      jsonToPost = JSON.stringify(m);
      k = JSON.parse(jsonToPost);
      ajax_data = {
        json: JSON.stringify(k)
      };
      if (window.nomenclature_id !== void 0) {
        ajax_data.id = window.nomenclature_id;
      }
      $('#save-spinner').show();
      return $.ajax({
        url: '/view/genome/nomenclature',
        type: window.action_mode,
        dataType: 'json',
        data: ajax_data,
        success: function(response) {
          return window.location = "/view/genome/nomenclature/set/status.html";
        },
        error: function(response) {
          return alert("Sorry, an error occurred trying to save this nomenclature.");
        },
        complete: function() {
          return $('#save-spinner').hide();
        }
      });
    });
    return $('.load-nomenclature').bind('click', function() {
      return listView.collection.reset(load_json);
    });
  });
}).call(this);
