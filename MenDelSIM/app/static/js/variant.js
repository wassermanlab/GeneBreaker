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
// on genome change
$(function () {
    $("#genomeSelect").change(function () {
        $("#transcriptSelect").empty();
    });
});
// show the general info -> var1
$(function () {
    $("#generalInfoNext").click(function () {
        // remove previous warnings
        $(".g-warning").remove()
        // check that all fields are full and valid 
        var sex = $("#sexSelect").val()
        var trancript_uid = $("#transcriptSelect option:Selected").attr("id")
        var chrom = $("#transcriptSelect option:Selected").attr("chrom")
        if (trancript_uid === undefined) {
            $("#generalInfoNext").before("<div class='alert alert-warning g-warning' role='alert'>No transcript selected.</div>");
            return;
        } else if (chrom === "chrY" && sex != "XY") {
            $("#generalInfoNext").before("<div class='alert alert-warning g-warning' role='alert'>Cannot select a transcript on the Y chromosome if the proband's sex is not XY.</div>");
            return;
        }
        // switch screens
        $(".variant_general_info").hide()
        $(".var_1").show()
    });
});
// on transcript select reset options for var1 and var2
$(function () {
    $("#transcriptSelect").change(function () {
        // fill var zygosity   
        $('.var_1').find(".zygositySelect").children().remove()
        if ((chrom === "chrY" || chrom === "chrX") && sex === "XY") {
            $('.var_1').find(".zygositySelect").append('<option class="zygosityHemi">Hemizygous</option>')
        } else {
            $('.var_1').find(".zygositySelect").append('<option class="zygosityHomo">Homozygous</option>')
            $('.var_1').find(".zygositySelect").append('<option class="zygosityHet">Heterozygous</option>')
        }
        // fill var1 options
        $('.var_1').find(".typeSelect").prop("selectedIndex", 0);
        $('.var_1').find(".regionSelect").prop("selectedIndex", 0);
    })
})
// var1 -> general info
$(function () {
    $("#var_1Back").click(function () {
        $(".var_1").hide()
        $(".variant_general_info").show()
    });
});
// var1 -> var2
$(function () {
    $("#var_1Next").click(function () {
        $(".var_1").hide()
        $(".var_2").show()
    });
});
// var2 -> var1
$(function () {
    $("#var_2Back").click(function () {
        $(".var_2").hide()
        $(".var_1").show()
    });
});