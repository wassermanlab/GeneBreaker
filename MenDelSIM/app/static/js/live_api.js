
var resources;
$(function () {
    $.getJSON("http://127.0.0.1:5000/json_docs", function (json) {
        resources = json;
    });
})

// create table row 
function create_row(param, params) {
    var newRow = document.createElement("tr");
    newRow.id = param;
    $("#paramRows").append(newRow);
    // create check  box
    var check = document.createElement("td");
    check.id = param + "Check";
    var form_div = document.createElement("div");
    form_div.id = param + "DivCheck";
    form_div.classList.add("form-check");
    var check_input = document.createElement("input");
    check_input.classList.add("form-check-input");
    check_input.id = param + "RowCheck";
    check_input.setAttribute("type", "checkbox");
    // create key box
    var key = document.createElement("td");
    key.id = param + "Key";
    key.innerHTML = param;
    // create value box
    var value = document.createElement("td");
    value.id = param + "Val";
    var input = document.createElement("input");
    input.type = "text";
    input.classList.add("form-control");
    input.id = param + "Input";
    // create description box
    var desc = document.createElement("td");
    desc.id = param + "Desc";
    desc.innerHTML = params[param]['DESCRIPTION'];
    // create required  box
    var req = document.createElement("td");
    req.innerHTML = params[param]['REQUIRED'];
    req.id = param + "Req";
    // add elements to dom
    $("#" + newRow.id).append(check);
    $("#" + check.id).append(form_div);
    $("#" + form_div.id).append(check_input);
    if (params[param]['REQUIRED']) {
        $("#" + check_input.id).prop("checked", true);
        $("#" + check_input.id).prop("disabled", true);
    }
    $("#" + newRow.id).append(key);
    $("#" + newRow.id).append(value);
    $("#" + value.id).append(input);
    $("#" + newRow.id).append(desc);
    $("#" + newRow.id).append(req);
}

// on change of select resource 
// build url
// populate param list
$(function () {
    $("#resourceSelect").change(function () {
        $("#paramRows").empty();
        var resource = $("#resourceSelect").val()
        if (resource === "select"){
            build_url()
            return
        }
        var params = resources[resource]['PARAMS']
        var keys = Object.keys(params)
        if (keys.includes("genome")) {
            create_row("genome", params)
        }
        if (keys.includes("chrom")) {
            create_row("chrom", params)
        }
        if (keys.includes("start")) {
            create_row("start", params)
        }
        if (keys.includes("end")) {
            create_row("end", params)
        }
        if (keys.includes("location")) {
            create_row("location", params)
        }

        for (param in params) {
            if (!["chrom", "start", "end", "location", "genome"].includes(param)) {
                create_row(param, params)
            }
        }
        build_url()
    });
})

// builds url based on selections on page
function build_url() {
    url = ""
    // get selected resource
    var resource = $("#resourceSelect").val()
    if (resource == "select") {
        $("#url").val(url);
        return;
    }
    // for each selected row get key value pair
    var parameters = {};
    var params;
    var check;
    var key;
    var val;
    $("#paramRows").children().each(function (index) {
        params = this.id;
        check = $("#" + params + "RowCheck").prop("checked")
        key = $("#" + params + "Key").text()
        val = $("#" + params + "Input").val()
        if (check) {
            parameters[key] = val;
        }
    });
    
    // get genome only and remove from rest of parameters 
    genome = parameters["genome"]
    delete parameters.genome
    
    // combine and set get request row
    url = resource
    parameters = Object.entries(parameters);
    if (parameters.length != 0) {
        url = url + "?" + parameters[0][0] + "=" + parameters[0][1];
        for (var i = 1; i < parameters.length; i++) {
            url = url + "&" + parameters[i][0] + "=" + parameters[i][1];
        }
    }

    url = url.replace("{genome}", genome)

    $("#url").val(url);
}

// on selection/deselction of parameter 
// build url

$(document).on('change', '.form-check-input', function () {
    build_url()
});

// on change of parameter value
// build url
$(document).on('change', '.form-control', function () {
    build_url()
});

// on send 
$(function () {
    $("#sendButton").click(function () {
        url = $("#url").val();
        url = "http://127.0.0.1:5000" + url;
        $(".responseCode").html("Loading ...")
        $.ajax({
            url: url,
            dataType: 'json',
            success: function( data ) {
                $(".responseCode").html(JSON.stringify(data, null, 2));
            },
            error: function( data ) {
                error = data['responseText'].match(/>(.*?)<\//g);
                json = {status: error[0].substr(1, error[0].length-3),
                     statusText: error[1].substr(1, error[1].length-3),
                     error: error[2].substr(1, error[2].length-3)}
                $(".responseCode").html(JSON.stringify(json, null, 2));
            }
          });



        // $.getJSON(url, function (json) {
        //     $(".responseCode").html(JSON.stringify(json, null, 2))
        // });
    });
});

// on copy
$(function () {
    $("#copyButton").click(function () {
        var copyText = document.querySelector("#url");
        copyText.select();
        document.execCommand("copy");
    });
});

