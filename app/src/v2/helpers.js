import { host } from '../host'
function check_var(type, region, zygosity, clinvar_id, start, end, length, element, snv_type, str_id, variant) {
  let errors = [];
  if (variant === "false") {
    return errors;
  }
  // if region, type or zygosity === ""
  if (region === "" || type === "" || zygosity === "") {
    errors.push("Ensure all fields are filled.")
  }
  switch (type) {
    case "CLINVAR":
      if (clinvar_id === "") {
        errors.push("Ensure a clinvar variant is selected.")
      }
      break;
    case "CLINGEN":
    case "CNV":
      if (start === "" || end === "" || length === "") {
        errors.push("Ensure all fields for CNVs are filled out.")
      }
      break;
    case "INDEL":
      if (start === "" || length === "") {
        errors.push("Ensure all fields for Indels are filled out.")
      }
      break;
    case "MEI":
      if (start === "" || element === "") {
        errors.push("Ensure all MEI fields are filled out.")
      }
      break;
    case "SNV":
      if (start === "" || snv_type === "") {
        errors.push("Ensure all fields for SNVs are filled out.")
      }
      break;
    case "STR":
      if (str_id === "" || length === "") {
        errors.push("Ensure all STR fields are filled out.")
      }
      break;
    default:
      break;
  }
  return errors;
}

function check_page_1(chrom, gene_uid, sex) {
  let errors = []
  if (gene_uid === "") {
    errors.push("Must select transcript.")
  }
  if (chrom === "Y" && sex === "XX") {
    errors.push("cannot select Y chromosome gene with XX proband.")
  }
  return errors;
}

export function check_errors(props) {
  let variant = props.page - 1;
  switch (props.page) {
    case 1:
      return check_page_1(props.chrom, props.gene_uid, props.sex);
    case 2:
      return check_var(
        props["type_" + variant],
        props["region_" + variant],
        props["zygosity_" + variant],
        props["clinvar_id_" + variant],
        props["start_" + variant],
        props["end_" + variant],
        props["length_" + variant],
        props["element_" + variant],
        props["snv_type_" + variant],
        props["str_id_" + variant],
        "true");
    case 3:
      return check_var(
        props["type_" + variant],
        props["region_" + variant],
        props["zygosity_" + variant],
        props["clinvar_id_" + variant],
        props["start_" + variant],
        props["end_" + variant],
        props["length_" + variant],
        props["element_" + variant],
        props["snv_type_" + variant],
        props["str_id_" + variant],
        props.var2);
    case 4:
      break
    default:
      break
  }
}

function create_variant(props, variant) {
  const region = props[("region_" + variant)];
  const zygosity = props[("zygosity_" + variant)];
  let type = props[("type_" + variant)];
  let impact;
  switch (type) {
    case "CLINVAR":
      impact = { CLINVAR_ID: parseInt(props[("clinvar_id_" + variant)]) }
      break;
    case "CLINGEN":
    case "CNV":
      type = "CNV"
      impact = {
        START: parseInt(props[("start_" + variant)]),
        END: parseInt(props[("end_" + variant)]),
        COPY_CHANGE: parseInt(props[("length_" + variant)])
      }
      break;
    case "INDEL":
      impact = {
        START: (props[("start_" + variant)] === "ANY" ? "ANY" : parseInt(props[("start_" + variant)])),
        INDEL_AMOUNT: parseInt(props[("length_" + variant)])
      }
      break;

    case "MEI":
      impact = {
        START: (props[("start_" + variant)] === "ANY" ? "ANY" : parseInt(props[("start_" + variant)])),
        ELEMENT: props[("element_" + variant)]
      }
      break;
    case "SNV":
      impact = {
        START: (props[("start_" + variant)] === "ANY" ? "ANY" : parseInt(props[("start_" + variant)])),
        SNV_TYPE: props[("snv_type_" + variant)]
      }
      break;
    case "STR":
      impact = {
        STR_ID: parseInt(props[("str_id_" + variant)]),
        LENGTH: parseInt(props[("length_" + variant)])
      }
      break;
    default:
      break;
  }
  return {
    REGION: region,
    ZYGOSITY: zygosity,
    TYPE: type,
    IMPACT: impact
  }
}

export async function get_variants(props) {
  let json = {
    GENE_UID: parseInt(props.gene_uid),
    GENOME: props.genome,
    SEX: props.sex,
  }
  json["VAR1"] = create_variant(props, 1)
  if (props.var2 === "true") {
    json["VAR2"] = create_variant(props, 2)
  } else {
    json["VAR2"] = "None"
  }
  // console.log(JSON.stringify(json))
  const rawResponse = await fetch(host + 'design_variants', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json'
    },
    body: JSON.stringify(json),
  });
  const vcf = await rawResponse.json();
  return vcf;
}
