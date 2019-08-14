var gud_host = "http://127.0.0.1:5000"
var mendelsim_host = "http://127.0.0.1:5001"

$(function () {
    $("#getTranscripts").click(function () {
        $("#transcriptSelect").empty();
        $("#noTranscript").hide();
        var genome = $('#genomeSelect').val()
        var name = $("#genes").val()
        var url = mendelsim_host + '/get_transcripts/' + genome + '/' + name
        $.getJSON(url, function (data) {
            if (!data) {
                $("#noTranscript").show();
            } else {
                var new_option;
                data.forEach(function (i) {
                    exonStart = i.qualifiers.exonStarts.split(",")[0]
                    exonEnd = i.qualifiers.exonEnds.split(",")
                    chrom = i.chrom
                    exonEnd = exonEnd[exonEnd.length - 2]
                    new_option = document.createElement("option")
                    new_option.id = i.qualifiers.uid
                    new_option.setAttribute('chrom', chrom);
                    new_option.innerHTML = "Accession: " + i.qualifiers.name + "&emsp;&emsp;Chrom: " + chrom + "&emsp;&emsp;Exon Start: " + exonStart + "&emsp;&emsp;Exon End: " + exonEnd
                    $("#transcriptSelect").append(new_option);
                })
            }
        })
    });
});