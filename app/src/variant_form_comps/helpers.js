function check_general(chrom, sex, transcript){
    let errors=[];
    if (transcript === ""){
        errors.push("No transcript is selected.")
    } 
    if (sex === "XX" && chrom === "chrY") {
        errors.push("Cannot select Y chromosome gene with XX proband.")
    }
    return errors;
}

function check_errors(props) {
    let errors;
    if(props.page === 1) {
        errors = check_general(props.chrom, props.sex, props.gene_uid)
    }
    return errors;
}
export default check_errors;
