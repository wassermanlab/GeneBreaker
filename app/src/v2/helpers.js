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
      break
    case 3:
      break
    case 4: 
      break
    default:
      break
  }
}
export default check_errors;
