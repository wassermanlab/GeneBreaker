$(function () {
    $('button#geneTranscriptB').bind('click', function () {
        $.getJSON($SCRIPT_ROOT + '/_get_transcript', {
            gene_name: $("#geneName").val()
        }, function (data) {
            $("#geneTranscriptL").empty();
            for (i = 0; i < data.length; i++) {
                $('<option/>', {
                    'data-chrom': data[i].chrom,
                    'data-uid': data[i].uid,
                    'data-name': data[i].name,
                    class: 'geneTranscript',
                    text: data[i].accession
                }).appendTo('#geneTranscriptL');
            }
        });
        return false;
    });
});
//--onClick of "Next on Global Information"--//
$(function () {
    $('a#toVariant').click(function () {
        $(".invalid-feedback").hide();
        error = global_error_checks()
        if (error) {
            return true;
        }
        url = $("a#toVariant").attr("data-href");
        gene_uid    = $('#geneTranscriptL').children("option:selected").attr('data-uid'),
        chrom       = $('#geneTranscriptL').children("option:selected").attr('data-chrom'),
        sex         = $('#probandSex').children("option:selected").val()
        window.location.href = url + "?gene_uid=" + gene_uid + "&chrom=" + chrom + "&sex=" + sex;
    });
});

//--server side global errors--//
function global_error_checks() {
    var error = false;
    var chrom = $('#geneTranscriptL').children("option:selected").attr('data-chrom');
    var sex = $('#probandSex').children("option:selected").val();
    if (chrom === undefined) {
        $(".geneTranscriptE").show();
        error = true
    }
    if (sex === undefined || sex === "Select") {
        $(".sexE").show();
        error = true
    }
    if (sex === "XX-Female" && chrom === "chrY") {
        error = true
        $(".sexYchromE").show();
    }
    if (error === true) {
        $(".missing-error").show();
        return true;
    }
    return false;
}
//--server side global errors--//
function var_basic_error_checks() {
    var error = false;
    var type = $('#var1_type').children("option:selected").val();
    var region = $('#var1_region').children("option:selected").val();
    var cregion = $('#var1_region_custom').val();
    var zygosity = $('#var1_zygosity').children("option:selected").val();

    if (type === "Select") {
        $(".missing-var-type").show();
        error = true
    }
    if (region === "Select") {
        $(".missing-var-region").show();
        error = true
    }
    if (zygosity === "Select") {
        $(".missing-var-zygosity").show();
        error = true
    }
    var reg = /^chr([1-9]|[XY]|1[0-9]|2[0-2])\d+-\d+$/;
    if (region === "chrZ:int-int" && !reg.test(cregion)) {
        $(".missing-var-cregion").show();
        error = true
    }
    if (error === true) {
        $(".missing-var-missing").show();
        return true;
    }
    return false;
}