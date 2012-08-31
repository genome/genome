(function(global){     
    var api = ["log","debug","info","warn","error","assert","dir","dirxml",
"trace","group","groupCollapsed","groupEnd","time","timeEnd",
"profile","profileEnd","count","exception","table"],
      log, i, len
   ;

   if (typeof global.console == "undefined" || !global.console) {
      try { global.console = {}; } catch (err) { }
   }

   log = (typeof global.console.log != "undefined") ?
      global.console.log :
      function(){}
   ;

   for (i=0, len=api.length; i<len; i++) {
      if (typeof global.console[api[i]] == "undefined" ||
         !global.console[api[i]])
      {
         try { global.console[api[i]] = log; } catch (err) { }
      }
   }
})(window);