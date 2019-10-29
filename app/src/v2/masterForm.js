import React from 'react';
// import { saveAs } from 'file-saver';
// import { host } from '../host'
import Nav from '../nav';
import GeneralInfo from './generalInfo';
import VariantInfo from './variantInfo';
import NavButtons from './navButtons';

class MasterForm2 extends React.Component {
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
      var2: "false",
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
      vars: {},
      family: {}
    };

    this.handleInputChange = this.handleInputChange.bind(this);
    this.next = this.next.bind(this)
    this.back = this.back.bind(this)
  }

  // sets page to page+1
  next() {
    // const errors = []
    // if (errors.length !== 0) {
    //   this.setState({ errors: errors });
    //   return null;
    // }
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
  
  resetState(name) {
    switch (name) {
      case "genome":
      case "gene_name": // resets everything past var 1 and gene_uid/chrom
        this.setState({
          gene_uid: "",
          chrom: ""
        });
        /* falls through */
      case "gene_uid": // resets everything past var 1 
        this.setState({
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
          str_id_1: ""
        });
        /* falls through */
      case "zygosity_1": // resets everything past var 2
        this.setState({
          var2: "false",
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
          str_id_2: ""
        });
        break
      case "type_1": // resets everything past var 1 
        console.log("type1")
        this.setState({
          clinvar_id_1: "",
          start_1: "",
          end_1: "",
          clingen_id_1: "",
          length_1: "",
          element_1: "",
          snv_type_1: "",
          str_id_1: ""
        });
        break
      case "type_2": // resets everything past var 1 
        this.setState({
          clinvar_id_2: "",
          start_2: "",
          end_2: "",
          clingen_id_2: "",
          length_2: "",
          element_2: "",
          snv_type_2: "",
          str_id_2: ""
        });
        break
      default:
        break
    }
  }

  handleInputChange(event) {
    const target = event.target;
    const value = target.value;
    const name = target.name;
    this.resetState(name)
    this.setState({
      [name]: value
    }, () => { console.log(this.state); });
  }

  render() {
    return (
      <div className="master-background">
        <Nav />
        <div className="container">
          <div className="formDiv">
            <form>
              {/* generalInfo */}
              <GeneralInfo page={this.state.page} gene_uid={this.state.gene_uid} genome={this.state.genome}
                gene_name={this.state.gene_name} sex={this.state.sex} onChange={this.handleInputChange} />
              {/* var1Info */}
              <VariantInfo
                var={1} page={this.state.page} type={this.state.type_1} region={this.state.region_1}
                customStart={this.state.customStart_1} customEnd={this.state.customEnd_1} zygosity={this.state.zygosity_1}
                clinvar_id={this.state.clinvar_id_1} start={this.state.start_1} end={this.state.end_1}
                clingen_id={this.state.clingen_id_1} length={this.state.length_1} element={this.state.element_1}
                snv_type={this.state.snv_type_1} str_id={this.state.str_id_1} genome={this.state.genome}
                chrom={this.state.chrom} gene_uid={this.state.gene_uid} sex={this.state.sex} onChange={this.handleInputChange} />
              {/* var2Info */}
              <VariantInfo
                zygosity_1={this.state.zygosity_1} var={2} page={this.state.page} type={this.state.type_2} region={this.state.region_2}
                customStart={this.state.customStart_2} customEnd={this.state.customEnd_2} zygosity={this.state.zygosity_2}
                clinvar_id={this.state.clinvar_id_2} start={this.state.start_2} end={this.state.end_2}
                clingen_id={this.state.clingen_id_2} length={this.state.length_2} element={this.state.element_2}
                snv_type={this.state.snv_type_2} str_id={this.state.str_id_2} genome={this.state.genome}
                chrom={this.state.chrom} gene_uid={this.state.gene_uid} sex={this.state.sex} onChange={this.handleInputChange} />
              {/* familyInfo */}
              {/* buttons */}
              <NavButtons page={this.state.page} next={this.next} back={this.back} />
            </form>
          </div >
        </div >
      </div>
    );
  }
}

export default MasterForm2;
