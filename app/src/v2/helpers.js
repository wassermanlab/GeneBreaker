function check_var(type, region, zygosity, clinvar_id, start, end, length, element, snv_type, str_id, variant=true) {
  let errors = [];
  if (variant === false) {
    return errors;
  } 
  // if region, type or zygosity === ""
  if (region === "" || type === "" || zygosity === "") {
    errors.push("Insure all fields are filled.")
  }
  switch (type) {
    case "CLINVAR":
      if (clinvar_id === "") {
        errors.push("Insure a clinvar variant is selected.")
      }
      break;
    case "CLINGEN":
    case "CNV":
      if (start === "" || end === "" || length === "") {
        errors.push("Insure all fields for CNVs are filled out.")
      }
      break;
    case "INDEL":
      if (start === "" || length === "") {
        errors.push("Insure all fields for Indels are filled out.")
      }
      break;
    case "MEI":
      if (start === "" || element === "") {
        errors.push("Insure all MEI fields are filled out.")
      }
      break;
    case "SNV":
      if (start === "" || snv_type === "" ) {
        errors.push("Insure all fields for CNVs are filled out.")
      }
      break;
    case "STR":
      if (str_id === "" || length === "") {
        errors.push("Insure all STR fields are filled out.")
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
  if (chrom === "chrY" && sex === "XX") {
    errors.push("cannot select Y chromosome gene with XX proband.")
  }
  return errors;
}

function check_errors(props) {
  switch (props.page) {
    case 1:
      return check_page_1(props.chrom, props.gene_uid, props.sex)
    case 2:
    case 3:
      const variant = props.page - 1;
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
export default check_errors;
