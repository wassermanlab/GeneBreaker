$(function () {
    $('button#gene').bind('click', function () {
        $.getJSON($SCRIPT_ROOT + '/_get_transcript', {
            gene_name: $('input[name="gene_name"]').val()
        }, function (data) {
            $("#gene_transcript").empty();
            for (i = 0; i < data.length; i++) {
                $('<option/>', {
                    'data-chrom': data[i].chrom,
                    'data-uid': data[i].uid,
                    'data-name': data[i].name,
                    class: 'gene_transcript',
                    text: data[i].accession
                }).appendTo('#gene_transcript');
            }
        });
        return false;
    });
});
//--onClick of "Next on Global Information"--//
$(function () {
    $('a#next-global').click(function () {
        $(".invalid-feedback").hide();
        error = global_error_checks()
        if (error) {
            return true;
        }
        var here = $(this).parent();
        var next = $(this).parent().next();
        here.hide();
        var chrom = $('#gene_transcript').children("option:selected").attr('data-chrom');
        var sex = $('#sex').children("option:selected").val();
        if (sex === "XY-Male" && (chrom === "chrX" || chrom === "chrX")) {
            $('#homo').hide()
            $('#het').hide()
        }
        $('#hemi').hide()
        next.show();
    });
});

//--server side global errors--//
function global_error_checks() {
    var error = false;
    var chrom = $('#gene_transcript').children("option:selected").attr('data-chrom');
    var sex = $('#sex').children("option:selected").val();
    if (chrom === undefined) {
        $(".gene-transcript-error").show();
        error = true
    }
    if (sex === undefined || sex === "Select") {
        $(".sex-error").show();
        error = true
    }
    if (sex === "XX-Female" && chrom === "chrY") {
        error = true
        $(".sex-y-gene-error").show();
    }
    if (error === true) {
        $(".missing-error").show();
        return true;
    }
    return false;
}
//--existing not existing radios--//
$(function () {
    $('input[type=radio][name=ne_radios]').change(function () {
        if (this.value == 'e_radio') {
            $('#var1_type').prop('selectedIndex', 0);
            $('.new').hide();
            $('.existing').show();
        }
        else if (this.value == 'n_radio') {
            $('#var1_type').prop('selectedIndex', 0);
            $('.existing').hide();
            $('.new').show();
        }
    });
});
// custom stuff
$(function () {
    $("#var1_region").change(function () {
        var value = $("#var1_region option:selected").val();
        $('#var1_region_custom').hide();
        if (value === "chrZ:int-int") {
            $('#var1_region_custom').show();
            return true;
        }
    });
});
//--onClick of "fill in details "--//
$(function () {
    $('a#var1-details').click(function () {
        $(".invalid-feedback").hide();
        // check if region, zygosity, variant type is filled out
        // lock region, zygosity, variant type, existing
        // unhide variant type + ne  
                
        
        
        
        error = var_error_checks()
        if (error) {
            return true;
        }
    });
});