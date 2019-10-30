function check_var(variant, region, zygosity, clinvar_id, start, end, clingen_id, length, element, snv_type, str_id){
  let errors = []
  return errors;
}

function check_page_1(chrom, gene_uid, sex){
  let errors = []
  if (gene_uid === "") {
    errors.push("Must select transcript.")
  }
  if (chrom === "chrY" && sex === "XX") {
    errors.push("cannot select Y chromosome gene with XX proband.")
  }
  return errors;
}

function check_errors(props){
  switch(props.page){
    case 1:
      return check_page_1(props.chrom, props.gene_uid, props.sex)
    case 2:
    case 3:
      const variant = props.page -1 ;  
      return check_var(
        variant, 
        props["region_" + variant],
        props["zygosity_" + variant],
        props["clinvar_id_" + variant],
        props["start_" + variant],
        props["end_" + variant],
        props["clingen_id_" + variant],
        props["length_" + variant],
        props["element_" + variant],
        props["snv_type_" + variant],
        props["str_id_" + variant]);
    case 4: 
      break
    default:
      break
  }
}
export default check_errors;
