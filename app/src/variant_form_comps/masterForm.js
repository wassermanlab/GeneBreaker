import React from 'react';
import GlobalInfo from './globalInfo'
import VariantInfo from './variantInfo'
import NavButtons from './navButtons'
import check_errors from './helpers.js'
import SNV from './snv'
import CNV from './cnv'
import Indel from './indel'
import MEI from './mei'
// import VariantInfo from './variantInfo'
// import FInfo from './familyInfo'

class MasterForm extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      errors: [],
      page: 1,
      // general state
      gene_uid: "",
      genome: "hg38",
      gene_name: "",
      chrom: "",
      sex: "XX",
      // var 1 state
      type_1: "",
      region_1: "",
      customStart_1: "",
      customEnd_1: "",
      zygosity_1: "",
      clinvar_id_1: "",
      start_1: "",
      end_1: "",
      clingen_id_1: "",
      length_1: "",
      element_1: "",
      snv_type_1: "",
      str_id_1: "",
      //var2 state
      var2: false,
      type_2: "",
      region_2: "",
      customStart_2: "",
      customEnd_2: "",
      zygosity_2: "",
      clinvar_id_2: "",
      start_2: "",
      end_2: "",
      clingen_id_2: "",
      length_2: "",
      element_2: "",
      snv_type_2: "",
      str_id_2: "",
      vars: {
        header: "##fileDate=2019-09-13↵##source=variant_simulator↵#…M	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	PROBAND↵",
        var1: {
          alt: "T",
          chrom: "chr17",
          filter: ".",
          format: "GT",
          id: "450264",
          info: ".",
          pos: "72121406",
          proband: "0/1",
          qual: ".",
          ref: "C",
        },
        var2: {
          alt: "A",
          chrom: "chr17",
          filter: ".",
          format: "GT",
          id: "450264",
          info: ".",
          pos: "72121406",
          proband: "0/1",
          qual: ".",
          ref: "C",
        }
      },
      family: [{ relationship: "", sex: "", var1: "", var2: "" }]
    };

    this.handleInputChange = this.handleInputChange.bind(this);
    this.handleInputChangeGeneral = this.handleInputChangeGeneral.bind(this);
    this.next = this.next.bind(this)
    this.back = this.back.bind(this)
  }

  // sets page to page+1
  next() {
    const errors = check_errors(this.state) 
    if (errors.length !== 0 ) {
      this.setState({ errors: errors });
      return null;
    }
    let currentPage = this.state.page;
    currentPage = currentPage + 1;
    this.setState({ page: currentPage, errors: [] });
  }
  // sets page to page-1
  back() {
    let currentPage = this.state.page;

    currentPage = currentPage - 1;

    this.setState({ page: currentPage });
  }
  handleInputChange(event) {
    const target = event.target;
    const value = target.value;
    const name = target.name;

    // state change for everything else 
    this.setState({
      [name]: value
    },
      () => console.log(this.state));
  }
  handleInputChangeGeneral(event) {
    const target = event.target;
    const value = target.value;
    const name = target.name;

    // state change for everything else 
    if (name === "gene_name" || name === "genome") {
      this.setState({ [name]: value, gene_uid: "", chrom: "" },
        () => console.log(this.state));
    } else if (name === "gene_uid" && value !== "") {
      const index = target.selectedIndex;
      const option = target.childNodes[index];
      const chrom = option.getAttribute('chrom');
      this.setState({
        gene_uid: value,
        chrom: chrom
      },
        () => console.log(this.state));
    } else { this.setState({ [name]: value }, () => console.log(this.state)); }
  }

  render() {
    return (
      <form>
        {/* global */}
        <GlobalInfo page={this.state.page} gene_uid={this.state.gene_uid} genome={this.state.genome}
          gene_name={this.state.gene_name} sex={this.state.sex} onChange={this.handleInputChangeGeneral} />
        {/* var1 */}
         <VariantInfo page={this.state.page} var={1} sex={this.state.sex} chrom={this.state.chrom}
          type={this.state.type_1} region={this.state.region_1}
          customStart={this.state.customStart_1} customEnd={this.state.customEnd_1} zygosity={this.state.zygosity_1} onChange={this.handleInputChange} />
          {/* var2 */}
        {/*<ClinVar type={this.state.type_1} page={this.state.page} var={1}
          clinvar_id={this.state.clinvar_id_1} onChange={this.handleInputChange} />
        <ClinGen type={this.state.type_1} page={this.state.page} var={1}
          clingen_id={this.state.clingen_id_1} onChange={this.handleInputChange}/>
           <STR type={this.state.type_1} page={this.state.page} var={1} str_id={this.state.str_id_1}  onChange={this.handleInputChange} />
           */}
        <CNV type={this.state.type_1} page={this.state.page} var={1}
          start={this.state.start_1} end={this.state.end_1} length={this.state.length_1}  onChange={this.handleInputChange} />
        <Indel type={this.state.type_1} page={this.state.page} var={1}
          start={this.state.start_1} length={this.state.length_1}  onChange={this.handleInputChange} />
        <MEI type={this.state.type_1} page={this.state.page} var={1}
          element={this.state.element_1} onChange={this.handleInputChange} />
        <SNV type={this.state.type_1} page={this.state.page} var={1}
          start={this.state.start_1} snv_type={this.state.snv_type_1} onChange={this.handleInputChange} />
        {/* var2 */}
        {/* <VariantInfo page={this.state.page} var={2} sex={this.state.sex} chrom={this.state.chrom}
          type={this.state.type_2} region={this.state.region_2}
          customStart={this.state.customStart_2} customEnd={this.state.customEnd_2} zygosity={this.state.zygosity_2} />
        <ClinVar type={this.state.type_2} page={this.state.page} var={2}
          clinvar_id={this.state.clinvar_id_2} />
        <ClinGen type={this.state.type_2} page={this.state.page} var={2}
          clingen_id={this.state.clingen_id_2} />
        <CNV type={this.state.type_2} page={this.state.page} var={2}
          start={this.state.start_2} end={this.state.end_2} length={this.state.length_2} />
        <Indel type={this.state.type_2} page={this.state.page} var={2}
          start={this.state.start_2} length={this.state.length_2} />
        <MEI type={this.state.type_2} page={this.state.page} var={2}
          element={this.state.element_2} />
        <SNV type={this.state.type_2} page={this.state.page} var={2}
          start={this.state.start_2} snv_type={this.state.snv_type_2} />
        <STR type={this.state.type_2} page={this.state.page} var={2} str_id={this.state.str_id_2} /> */}
        {/* family */}
        {/* <FamilyInfo /> */}
        {/* buttons */}
        {this.state.errors.map((item, index) => (
            <div key={"error_"+index} className="alert alert-danger" role="alert">{item}</div>
          ))}
        <NavButtons page={this.state.page} next={this.next} />
      </form>
    );
  }
}

export default MasterForm;
