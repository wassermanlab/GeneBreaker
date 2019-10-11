function check_general(chrom, sex, transcript) {
  let errors = [];
  if (transcript === "") {
    errors.push("No transcript is selected.")
  }
  if (sex === "XX" && chrom === "chrY") {
    errors.push("Cannot select Y chromosome gene with XX proband.")
  }
  return errors;
}
function check_variant(variant, type, region, zygosity, customStart, customEnd,
  clinvar_id, start, end, clingen_id, length, element, snv_type, str_id) {
  if (!variant) {
    return []
  }
  let errors = [];
  if (type === "") {
    errors.push("No type selected")
  }
  if (zygosity === "") {
    errors.push("No zygosity selected")
  }
  if (region === "") {
    errors.push("No region selected")
  } else if (region === "CUSTOM" && parseInt(customStart) >= parseInt(customEnd)) {
    errors.push("Invalid custom region.")
  }
  let start_error = false;
  if (!((parseInt(start) >= 1) || (start === "ANY"))) {
    start_error = true;
  }
  switch (type) { // check types specific stuff 
    case "cnv":
      if (parseInt(start) >= end) {
        errors.push("start of cnv must be less than end.")
      }
      if (parseInt(length) < -1 || parseInt(length) === 0) {
        errors.push("copy change cannot be less than -1 or equal 0.")
      }
      break;
    case "clingen":
      if (parseInt(start) >= end) {
        errors.push("start of cnv must be less than end.")
      }
      if (parseInt(length) < -1 || parseInt(length) === 0) {
        errors.push("copy change cannot be less than -1 or equal 0.")
      }
      break;
    case "clinvar":
      if (clinvar_id === "") {
        errors.push("clinvar variant must be selected.")
      }
      break;
    case "indel":
      if (start_error) {
        errors.push("must pick valid start of indel.")
      }
      if (parseInt(length) > 200 || parseInt(length) < -200 || parseInt(length) === 0) {
        errors.push("indel length  must be between -200 and 200 and cannot equal to 0.")
      }
      break;
    case "mei":
      if (start_error) {
        errors["mei_start"] = "must pick valid start for mei."
      }
      break;
    case "snv":
      if (start_error) {
        errors.push("must pick valid start for snv.")
      }
      if (snv_type === "") {
        errors.push("must select valid snv type.")
      }
      break;
    case "str":
      if (start_error) {
        errors.push("must pick valid start for str.")
      }
      if (length === 0) {
        errors.push("length of str cannot equal 0")
      }
      break;
    default:
      break;
  }
  return errors;
}
function check_family(var1, var2, family){
  let errors = []
  const chrom = var1.chrom;
  let single = false;
  if (var2 === "" && var1.proband === "0/1") {
    single = true;
  }
  for (let m in family) {
    const member = family[m];
    if (member.sex === "XX" && chrom === "chrY" && member.var1) {
      errors.push(m+" cannot have a Y chromosome variant being a XX female.")
    } 
    if (member.sex === "XY" && member.var2 && member.var1 && chrom === "chrY") {  
      errors.push(m+" cannot have two Y chromosome variants being a XY male.")
    }
    if (member.sex === "XY" && member.var2 && member.var1 && chrom === "chrX") {  
      errors.push(m+" cannot have two X chromosome variants being a XY male.")
    }
    if (member.var1 && member.var2 && single) {
      errors.push(m+" cannot have variants because there is only 1 variant to choose from.")
    }
  }
  return errors;
}
function check_errors(props) {
  let errors;
  if (props.page === 1) {
    errors = check_general(props.chrom, props.sex, props.gene_uid)
  } else if (props.page === 2) {
    errors = check_variant(true, props.type_1, props.region_1, props.zygosity_1,
      props.customStart_1, props.customEnd_1, props.clinvar_id_1,
      props.start_1, props.end_1, props.clingen_id_1, props.length_1, props.element_1,
      props.snv_type_1, props.str_id_1)
  } else if (props.page === 3) {
    errors = check_variant(props.var2, props.type_2, props.region_2, props.zygosity_2,
      props.customStart_2, props.customEnd_2, props.clinvar_id_2,
      props.start_2, props.end_2, props.clingen_id_2, props.length_2, props.element_2,
      props.snv_type_2, props.str_id_2)
  } else if (props.page === 4) {
    errors = check_family(props.vars.var1, props.vars.var2, props.family)
  }
  return errors;
}
export default check_errors;
