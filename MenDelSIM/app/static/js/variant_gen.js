// $("#gene").click(function () {
//     console.log("hello")
//     $.getJSON('/get_transcript',
//         function (data) {
//         });
//     return false;
// });
$(function () {
    $('button#gene').bind('click', function () {
        $.getJSON($SCRIPT_ROOT + '/_get_transcript', {
            gene_name: $('input[name="gene_name"]').val()
        }, function (data) {
            console.log(data)
        });
        return false;
    });
});