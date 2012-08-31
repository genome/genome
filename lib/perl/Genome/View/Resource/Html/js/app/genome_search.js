

$('document').ready(function() {

    $('#searchBox').focus();

    $('#searchForm').submit(function() { 
        if ($('#searchBox').val() == '') { return false; }
        return true;
    });

    $('#authUser').change(function() {
        initSearchPage($(this).html());
    });
});


function initSearchPage(userName) {
    $("#project_link").attr("href","/view/genome/sys/user/status.html?id=" + userName + "@genome.wustl.edu");
}




