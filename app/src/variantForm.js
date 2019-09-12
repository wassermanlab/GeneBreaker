import React from 'react';
import GeneralInfo from './generalInfo'
import VariantInfo from './variantInfo'

class VFrom extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      page: 1,
      // general state
      gene_uid: "",
      genome: "hg38",
      gene_name: "",
      chrom: "",
      sex: "XX",
      // var 1 state
      var1_type: "",
      var1_region: "",
      var1_customStart: "",
      var1_customEnd: "",
      var1_zygosity: "",
      var1_clinvar_id: "",
      var1_start: "",
      var1_end: "",
      var1_clingen_uid: "",
      var1_length: "",
      var1_element: "",
      var1_snv_type: "",
      var1_str_id: "",
      //var2 state
      var2: false,
      var2_type: "",
      var2_region: "",
      var2_customStart: "",
      var2_customEnd: "",
      var2_zygosity: "",
      var2_clinvar_id: "",
      var2_start: "",
      var2_end: "",
      var2_clingen_uid: "",
      var2_length: "",
      var2_element: "",
      var2_snv_type: "",
      var2_str_id: "",
    };

    this.handleInputChange = this.handleInputChange.bind(this);
    this.submit = this.submit.bind(this)
    this.next = this.next.bind(this)
    this.back = this.back.bind(this)
  }
  get_region(variant) {
    if (this.state["var" + variant + "_region"] === "CUSTOM") {
      return this.state.chrom + ":" + this.state["var" + variant + "_customStart"] + "-" + this.state["var" + variant + "_customEnd"];
    }
    return this.state["var" + variant + "_region"];

  }
  get_impact(variant) {
    let impact = {}
    switch (this.state["var" + variant + "_type"]) {
      case "clinvar":
        impact.CLINVAR_ID = parseInt(this.state["var" + variant + "_clinvar_id"])
        break
      case "clingen":
        impact.START = parseInt(this.state["var" + variant + "_start"])
        impact.END = parseInt(this.state["var" + variant + "_end"])
        impact.COPY_CHANGE = parseInt(this.state["var" + variant + "_length"])
        break
      case "cnv":
        impact.START = parseInt(this.state["var" + variant + "_start"])
        impact.END = parseInt(this.state["var" + variant + "_end"])
        impact.COPY_CHANGE = parseInt(this.state["var" + variant + "_length"])
        break
      case "indel":
        impact.START = (this.state["var" + variant + "_start"] === "ANY") ? "ANY" : parseInt(this.state["var" + variant + "_start"])
        impact.INDEL_AMOUNT = parseInt(this.state["var" + variant + "_length"])
        break
      case "mei":
        impact.START = (this.state["var" + variant + "_start"] === "ANY") ? "ANY" : parseInt(this.state["var" + variant + "_start"])
        impact.ELEMENT = this.state["var" + variant + "_element"]
        break
      case "snv":
        impact.START = (this.state["var" + variant + "_start"] === "ANY") ? "ANY" : parseInt(this.state["var" + variant + "_start"])
        impact.SNV_TYPE = this.state["var" + variant + "_snv_type"]
        break
      case "str":
        impact.STR_ID = parseInt(this.state["var" + variant + "_str_id"])
        impact.LENGTH = parseInt(this.state["var" + variant + "_length"])
        break
      default:
        break
    }
    return impact;
  }
  async submit() {

    let config = {
      GENE_UID: parseInt(this.state.gene_uid),
      GENOME: this.state.genome,
      SEX: this.state.sex,
      VAR1: {
        TYPE: this.state.var1_type.toUpperCase(),
        REGION: this.get_region(1),
        IMPACT: this.get_impact(1),
        ZYGOSITY: this.state.var1_zygosity.toUpperCase()
      },
      VAR2: {
        TYPE: this.state.var2_type.toUpperCase(),
        REGION: this.get_region(2),
        IMPACT: this.get_impact(2),
        ZYGOSITY: this.state.var2_zygosity.toUpperCase()
      }
    }
    if (!this.state.var2) {
      config["VAR2"] = "None"
    }
    
    console.log(JSON.stringify(config))
    const rawResponse = await fetch('http://127.0.0.1:5001/design_variants', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify(config),
    });
    const vcf = await rawResponse.json();
    // console.log(vcf)
    this.props.history.push('/family', {vcf: vcf})




  }
  // sets page to page+1
  next() {
    let currentPage = this.state.page;
    currentPage = currentPage + 1;
    if (currentPage === 3) {
      this.setState({
        page: currentPage,
        var2: true
      });
    }
    this.setState({ page: currentPage });
  }
  // sets page to page-1
  back() {
    let currentPage = this.state.page;

    currentPage = currentPage - 1;
    if (currentPage === 2) {
      this.setState({
        page: currentPage,
        var2: false
      });
    }
    this.setState({ page: currentPage });
  }
  handleInputChange(event) {
    const target = event.target;
    const value = target.value;
    const name = target.name;

    //custom state change to deal with chrom
    if (name === "var1_clingen_uid" && value !== "") {
      const index = target.selectedIndex;
      const option = target.childNodes[index];
      const start = option.getAttribute('start');
      const end = option.getAttribute('end');
      const copy_change = option.getAttribute('copy_change');
      this.setState({
        var1_start: start,
        var1_end: end,
        var1_length: copy_change,
        var1_clingen_uid: value
      },
        () => console.log(this.state));
      return null;
    }
    if (name === "var2_clingen_uid" && value !== "") {
      const index = target.selectedIndex;
      const option = target.childNodes[index];
      const start = option.getAttribute('start');
      const end = option.getAttribute('end');
      const copy_change = option.getAttribute('copy_change');
      this.setState({
        var2_start: start,
        var2_end: end,
        var2_length: copy_change,
        var2_clingen_uid: value
      },
        () => console.log(this.state));
      return null;
    }

    if (name === "gene_uid" && value !== "") {
      const index = target.selectedIndex;
      const option = target.childNodes[index];
      const chrom = option.getAttribute('chrom');
      this.setState({
        gene_uid: value,
        chrom: chrom
      },
        () => console.log(this.state));
      return null;
    }
    // state change for everything else 
    this.setState({
      [name]: value
    },
      () => console.log(this.state));
  }
  render() {
    return (
      <form>
        {/* general: genome, sex, gene_uid, chr, transcript */}
        <GeneralInfo
          page={this.state.page}
          chrom={this.state.chrom}
          genome={this.state.genome}
          sex={this.state.sex}
          gene_name={this.state.gene_name}
          gene_uid={this.state.gene_uid}
          handleInputChange={this.handleInputChange}
          next={this.next}
          back={this.back}
        />
        {/* var1 */}
        {/* //page, chrom, sex, gene_uid, handleInputChange, next, back
        //var, var2, type, region, zygosity, clinvar_id, start, end, length, element, snv_type, motif */}
        <VariantInfo
          page={this.state.page}
          chrom={this.state.chrom} sex={this.state.sex} gene_uid={this.state.gene_uid} genome={this.state.genome}
          handleInputChange={this.handleInputChange} next={this.next} back={this.back} submit={this.submit}
          var={1} type={this.state.var1_type} region={this.state.var1_region} zygosity={this.state.var1_zygosity}
          customStart={this.state.var1_customStart} customEnd={this.state.var1_customEnd}
          start={this.state.var1_start} end={this.state.var1_end} length={this.state.var1_length}
          clingen_uid={this.state.var1_clingen_uid} clinvar_id={this.state.var1_clinvar_id} 
          element={this.state.var1_element} snv_type={this.state.var1_snv_type} str_id={this.state.var1_str_id}
        />
        {/* var2 */}
        <VariantInfo
          page={this.state.page}
          chrom={this.state.chrom} sex={this.state.sex} gene_uid={this.state.gene_uid} genome={this.state.genome}
          handleInputChange={this.handleInputChange} next={this.next} back={this.back} submit={this.submit}
          var={2} type={this.state.var2_type} region={this.state.var2_region} zygosity={this.state.var2_zygosity}
          customStart={this.state.var2_customStart} customEnd={this.state.var2_customEnd}
          start={this.state.var2_start} end={this.state.var2_end} length={this.state.var2_length}
          clingen_uid={this.state.var2_clingen_uid} clinvar_id={this.state.var2_clinvar_id}
          element={this.state.var2_element} snv_type={this.state.var2_snv_type} str_id={this.state.var2_str_id}
        />
      </form>
    );
  }
}

export default VFrom;
