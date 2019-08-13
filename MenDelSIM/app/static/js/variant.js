var gud_host = "http://127.0.0.1:5000"
var mendelsim_host = "http://127.0.0.1:5001"

$(function () {
    $("#getTranscripts").click(function () {
        var genome = $('#genomeSelect').val()
        var name = $("#genes").val()
        var url = mendelsim_host + '/get_transcripts/' + genome + '/' + name
        $.getJSON(url, function (data) {
            for (i in data) {

            }
        });
    });
});