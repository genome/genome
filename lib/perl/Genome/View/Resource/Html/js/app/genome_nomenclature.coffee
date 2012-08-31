window.action_mode='PUT'

$ ->

      templateContents = '''
        <script type="text/template" id="nomenclature-column-template">
          <label>Title:</label>
          <input type="text" class="nomenclature-column-name text" value="<%= model.get('name') %>">
          <label>Type:</label>
          <select class="column-type-select">
            <option value="integer">integer</option>
            <option value="real">real</option>
            <option value="numeric">numeric (real or integer)</option>
            <option value="date">date</option>
            <option value="string">string</option>
            <option value="enumerated">enumerated</option>
          </select>
          <button class="delete">Remove Column</button>
        </script>

        <script type="text/template" id="nomenclature-column-enumerated-template">
          <div class="enumerated-choices">
            <strong>Enumerated Column Definition:</strong><br/>
            <% _.each(model.get("enumerated_values"), function(e) { %>

            <div class="enum-row"><input type="text" class="enumerated-value-entry text" value="<%= e %>">

            <button class="remove-enum">remove</button>
            <button class="add-new-enum">add</button>
          </div>

          <% }); %>
        </div>
        </script>
      '''

      $("#script-container").append($(templateContents))
      Backbone.sync = (method, model, callback_bundle) ->
        callback_bundle.success()

      fieldCount = 0

      class Column extends Backbone.Model
        defaults:
          name: 'column'
          type: 'string'
          id: ''
          use_count: 0
          enumerated_values: []
          enumerated_value_ids: []
          enumerated_value_use_counts: []

      class List extends Backbone.Collection
        model: Column

      class ColumnView extends Backbone.View
        tagName: 'li'

        events:
          'click .delete': 'remove'
          'blur input.nomenclature-column-name': 'nameChanged'
          'change select.column-type-select': 'typeChanged'
          'click .add-new-enum': 'addNewEnum'
          'click .remove-enum': 'removeEnum'
          'blur input.enumerated-value-entry': 'enumChanged'

        nameChanged: ->
          new_name = $(@el).children('.nomenclature-column-name')[0].value
          @model.set({"name": new_name})

        typeChanged: ->
          type = $(@el).children('.column-type-select')[0].value
          if (type == 'enumerated')
            @model.set({"enumerated_values":['Untitled']}, {"silent": true})
          
          @model.set({"type": type})

        addNewEnum: ->
            types = @model.get "enumerated_values"
            types.push("untitled")
            @model.set({"enumerated_values": undefined}, {"silent" : true})
            @model.set({"enumerated_values": types})

        getEnumIndexForInitiator: (initiator) ->
            parent_row = initiator.parent(".enum-row")[0]
            entire_set = $($(parent_row).parent('.enumerated-choices')).children('.enum-row')
            
            for index in [0...entire_set.length]
                if entire_set.get(index) == parent_row
                    return index

        enumChanged: (k) ->
            initiator = $(k.target)
            index = @getEnumIndexForInitiator(initiator) 
            value = initiator.val()

            types = @model.get "enumerated_values"
            types[index] = value
            @model.set({"enumerated_values": undefined}, {silent: true})
            @model.set({"enumerated_values": types})
            
        removeEnum: (k) ->
            initiator = $(k.target)
           
            index = @getEnumIndexForInitiator(initiator) 
            types = @model.get "enumerated_values"
            ids = @model.get "enumerated_value_ids"
            types.splice(index, 1)
            ids.splice(index, 1)
            @model.set({"enumerated_values": undefined}, {silent: true})
            @model.set({"enumerated_values": types})


        initialize: ->
          _.bindAll @, 'render', 'unrender', 'remove'

          @model.bind 'change', @render
          @model.bind 'remove', @unrender

        render: ->
          item_template = _.template($('#nomenclature-column-template').html(), {model:@model})
          i = $(item_template)
        
          $(@el).html(i)
          $(@el).children(".column-type-select")[0].value = @model.get('type')
          if @model.get('type') == 'enumerated'
            enumerated_template =  _.template($('#nomenclature-column-enumerated-template').html(), {model:@model}) 
            $(@el).append(enumerated_template)
            kids = $(@el).find(".enum-row")
            for i in [0..kids.size()-1]
                if parseInt(@model.get('enumerated_value_use_counts')[i]) > 0
                    $(kids[i]).children('.remove-enum').attr('disabled','true')
                    $(kids[i]).children('.enumerated-value-entry').attr('disabled','true')
                    $('#why-cant-i-delete').show()
        
          if parseInt(@model.get('use_count')) > 0
            $(@el).children(".delete").attr('disabled','true')
            $(@el).children(".column-type-select").attr('disabled','true')
            $(@el).children(".nomenclature-column-name").attr('disabled','true')
            $('#why-cant-i-delete').show()
          this

        unrender: ->
          $(@el).remove()

        remove: ->
          this.model.destroy()

      class ListView extends Backbone.View
        el: $('#nomenclature-columns')

        events:
          'click button#add': 'addColumn'

        initialize: ->
          _.bindAll @, 'render', 'addColumn', 'appendColumn'

          @collection = new List
          @collection.bind 'add', @appendColumn
          @collection.bind 'reset', @render

          @render()
        
        appendColumn: (column) ->
          columnView = new ColumnView model: column
          $('ul#nomenclature-list', @el).append columnView.render().el

        render: ->
          $(@el).html('')
          $(@el).append "<ul id='nomenclature-list'></ul>"
          $(@el).append "<div id='nomenclature-add'><button id='add'>Add column</button></div>"
          _(@collection.models).each (column) -> appendColumn column, @

        addColumn: ->
          fieldCount++
          column = new Column
          name = column.get('name')
          name = name + fieldCount
          column.set({"name", name})
          @collection.add column


      listView = new ListView
    
      if (window.location.hash) 
        m = window.location.hash.split('=')
        if m[0] == '#id'
          url = '/view/Genome/Nomenclature/detail.json?id=' + m[1]
          $.ajax
           url: url
           type: 'GET'
           success: (data, textStatus, jqXHR) ->
            construct_from_json(data)
        
      construct_from_json = (object) -> 
        window.action_mode = 'POST'
        name = object.name
        $("#nomenclature-name-input").val(name)
        $("#nomenclature-empty-equivalent-input").val(object.empty_equivalent)
        window.nomenclature_id = object.id
        $(".title h1").html("Edit Nomenclature: #{name}")
        $("#directions").html("Use the form below to edit the nomenclature #{name}")
        document.title = "Edit Nomenclature: #{name}"
        _(object.fields).each (i) ->
            listView.collection.add(i)
        listView.render

      $('.save-nomenclature').bind 'click', ->
        name = $("#nomenclature-name-input").val()
        if name == ""
            alert "You can't save a nomenclature without a name!"
            return

        empty_equiv = $("#nomenclature-empty-equivalent-input").val()
        if empty_equiv == ""
            alert "Please specify an empty-equivalent string!"
            return

        if listView.collection.length == 0
            alert "You need at least one column to create a nomenclature!"
            return

        m = {"name": name, "empty_equivalent" : empty_equiv, "fields" : listView.collection}
        jsonToPost = JSON.stringify(m)
        k = JSON.parse(jsonToPost)
        ajax_data = {json:JSON.stringify(k)}
        if window.nomenclature_id != undefined
            ajax_data.id = window.nomenclature_id
        $('#save-spinner').show()
        $.ajax
           url: '/view/genome/nomenclature'
           type: window.action_mode
           dataType: 'json'
           data: ajax_data
           success: (response ) ->
            window.location="/view/genome/nomenclature/set/status.html"
           error: (response) ->
            alert("Sorry, an error occurred trying to save this nomenclature.")
           complete: ->
            $('#save-spinner').hide();

      $('.load-nomenclature').bind 'click', ->
        #load_json  = "[{'name':  'hello', 'type': 'string', 'enumerated_values': []},{'name':  'there', 'type': 'enumerated', 'enumerated_values': ['woo', 'there']} ]"
        listView.collection.reset(load_json)


          
      
