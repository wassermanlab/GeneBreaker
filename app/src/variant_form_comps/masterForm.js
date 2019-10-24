import React from 'react';
import GlobalInfo from './globalInfo'
import VariantInfo from './variantInfo'
import NavButtons from './navButtons'
import check_errors from './helpers.js'
import SNV from './snv'
import CNV from './cnv'
import Indel from './indel'
import MEI from './mei'
import ClinGen from './clingen'
import STR from './str'
import Clinvar from './clinvar'
import FInfo from './familyInfo'
import { saveAs } from 'file-saver';
import Nav from '../nav';
import './masterForm.css';
import host from '../config'

function Errors(props) {
  return (

    props.errors.map((item, index) => (
      <div key={"error_" + index} className="alert alert-danger" role="alert">{item}</div>
    ))

  )
}

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
    this.handleInputChangeGeneral = this.handleInputChangeGeneral.bind(this);
    this.handelInputChangeClearFields = this.handelInputChangeClearFields.bind(this);
    this.next = this.next.bind(this)
    this.back = this.back.bind(this)
    this.get_vars = this.get_vars.bind(this)
    this.addFamily = this.addFamily.bind(this)
    this.removeFamily = this.removeFamily.bind(this)
    this.handleFamilyCheckChange = this.handleFamilyCheckChange.bind(this)
    this.downloadFile = this.downloadFile.bind(this)
  }

  async downloadFile(event) {
    const errors = check_errors(this.state)
    if (errors.length !== 0) {
      this.setState({ errors: errors });
      return null;
    }

    let fam = this.state.family
    // two different variants
    // two the same
    // 1 variant
    if (this.state.vars.var2 === "" && this.state.vars.var1.Proband === "0/1") {
      fam['proband'] = { relationship: "sibling", sex: this.state.sex, var1: true, var2: false, affected: true }
    } else {
      fam['proband'] = { relationship: "sibling", sex: this.state.sex, var1: true, var2: true, affected: true }
    }

    const file_type = event.target.value;
    const rawResponse = await fetch(host + 'get_file?filetype=' + file_type, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({ var1: this.state.vars.var1, var2: this.state.vars.var2, family: fam }),
    });
    const blob = await rawResponse.blob();
    saveAs(blob, "test." + file_type)
  }

  handleFamilyCheckChange(event) {
    const checked = event.target.checked;
    const name = event.target.name.split("_");
    const family_key = name[1]
    const key = name[0]
    let fam = this.state.family;
    fam[family_key][key] = checked;
    this.setState({
      family: fam
    }, () => console.log(this.state))
  }

  removeFamily(event) {
    const member = event.target.value;
    let fam = this.state.family;
    delete fam[member];
    this.setState({
      family: fam
    }, () => console.log(this.state))
  }

  addFamily(event) {
    const member = event.target.value;
    const keys = Object.keys(this.state.family);
    const brothers = keys.filter((m) => { return m.startsWith("brother") });
    const sisters = keys.filter((m) => { return m.startsWith("sister") });
    let fam = this.state.family;
    switch (member) {
      case "m":
        if (!keys.includes('mother')) { fam.mother = { relationship: "mother", sex: "XX", var1: false, var2: false, affected: false } }
        break
      case "f":
        if (!keys.includes('father')) { fam.father = { relationship: "father", sex: "XY", var1: false, var2: false, affected: false } }
        break
      case "b":
        fam["brother" + (brothers.length + 1)] = { relationship: "sibling", sex: "XY", var1: false, var2: false, affected: false }
        break
      case "s":
        fam["sister" + (sisters.length + 1)] = { relationship: "sibling", sex: "XX", var1: false, var2: false, affected: false }
        break
      default:
        break
    }
    this.setState({
      family: fam
    }, () => console.log(this.state))
  }

  get_region(variant) {
    if (this.state["region_" + variant] === "CUSTOM") {
      return this.state.chrom + ":" + this.state["customStart_" + variant] + "-" + this.state["customEnd_" + variant];
    }
    return this.state["region_" + variant];
  }

  get_impact(variant) {
    let impact = {}
    switch (this.state["type_" + variant]) {
      case "clinvar":
        impact.CLINVAR_ID = parseInt(this.state["clinvar_id_" + variant])
        break
      case "clingen":
        impact.START = parseInt(this.state["start_" + variant])
        impact.END = parseInt(this.state["end_" + variant])
        impact.COPY_CHANGE = parseInt(this.state["length_" + variant])
        break
      case "cnv":
        impact.START = parseInt(this.state["start_" + variant])
        impact.END = parseInt(this.state["end_" + variant])
        impact.COPY_CHANGE = parseInt(this.state["length_" + variant])
        break
      case "indel":
        impact.START = (this.state["start_" + variant] === "ANY") ? "ANY" : parseInt(this.state["start_" + variant])
        impact.INDEL_AMOUNT = parseInt(this.state["length_" + variant])
        break
      case "mei":
        impact.START = (this.state["start_" + variant] === "ANY") ? "ANY" : parseInt(this.state["start_" + variant])
        impact.ELEMENT = this.state["element_" + variant]
        break
      case "snv":
        impact.START = (this.state["start_" + variant] === "ANY") ? "ANY" : parseInt(this.state["start_" + variant])
        impact.SNV_TYPE = this.state["snv_type_" + variant]
        break
      case "str":
        impact.STR_ID = parseInt(this.state["str_id_" + variant])
        impact.LENGTH = parseInt(this.state["length_" + variant])
        break
      default:
        break
    }
    return impact;
  }

  async get_vars() {

    let config = {
      GENE_UID: parseInt(this.state.gene_uid),
      GENOME: this.state.genome,
      SEX: this.state.sex,
      VAR1: {
        TYPE: this.state.type_1.toUpperCase(),
        REGION: this.get_region(1),
        IMPACT: this.get_impact(1),
        ZYGOSITY: this.state.zygosity_1.toUpperCase()
      },
      VAR2: {
        TYPE: this.state.type_2.toUpperCase(),
        REGION: this.get_region(2),
        IMPACT: this.get_impact(2),
        ZYGOSITY: this.state.zygosity_2.toUpperCase()
      }
    }
    if (this.state.var2 === "false") {
      config["VAR2"] = "None"
    }
    const rawResponse = await fetch(host + 'design_variants', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify(config),
    });
    const vcf = await rawResponse.json();
    if ("error" in vcf) {
      this.setState({ errors: [vcf["error"]] })
    } else {
      this.setState({ vars: vcf, page: 4 },
        () => console.log(this.state));
    }
  }
  // sets page to page+1
  next() {
    const errors = check_errors(this.state)
    if (errors.length !== 0) {
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
  handleInputChangeClinGen(event) {
    const target = event.target;
    const name = target.name;
    const index = target.selectedIndex;
    const option = target.childNodes[index];
    const start = option.getAttribute('start');
    const end = option.getAttribute('end');
    const copy_change = option.getAttribute('copy_change');

    // state change for everything else 
    if (name === "clingen_id_1") {
      this.setState({
        start_1: start,
        end_1: end,
        copy_change_1: copy_change,
      },
        () => console.log(this.state));
    } else {
      this.setState({
        start_2: start,
        end_2: end,
        copy_change_2: copy_change,
      },
        () => console.log(this.state));
    }

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
  handelInputChangeClearFields(event) {
    const target = event.target;
    const value = target.value;
    const name = target.name;
    const num = name[name.length - 1];
    // state change for everything else 
    this.setState({
      [name]: value,
      ["clinvar_id_" + num]: "",
      ["start_" + num]: "",
      ["end_" + num]: "",
      ["clingen_id_" + num]: "",
      ["length_" + num]: "",
      ["element_" + num]: "",
      ["snv_type_" + num]: "",
      ["str_id_" + num]: "",
    },
      () => console.log(this.state));
  }

  render() {
    return (
      <div className="master-background">
        <Nav />
        <div className="container">
          <div className="formDiv">
            <form>
              <div className="row">
                <div className="col-sm">
                  {/**************** global ****************/}
                  <GlobalInfo page={this.state.page} gene_uid={this.state.gene_uid} genome={this.state.genome}
                    gene_name={this.state.gene_name} sex={this.state.sex} onChange={this.handleInputChangeGeneral} />
                  {/**************** var1 ****************/}
                  <VariantInfo page={this.state.page} var={1} sex={this.state.sex} chrom={this.state.chrom}
                    type={this.state.type_1} region={this.state.region_1}
                    customStart={this.state.customStart_1} customEnd={this.state.customEnd_1} zygosity={this.state.zygosity_1}
                    onChange={this.handleInputChange} onChangeClear={this.handelInputChangeClearFields} />
                  <CNV type={this.state.type_1} page={this.state.page} var={1}
                    start={this.state.start_1} end={this.state.end_1} length={this.state.length_1} onChange={this.handleInputChange} />
                  <Indel type={this.state.type_1} page={this.state.page} var={1}
                    start={this.state.start_1} length={this.state.length_1} onChange={this.handleInputChange} />
                  <MEI type={this.state.type_1} page={this.state.page} var={1}
                    element={this.state.element_1} onChange={this.handleInputChange} />
                  <SNV type={this.state.type_1} page={this.state.page} var={1}
                    start={this.state.start_1} snv_type={this.state.snv_type_1} onChange={this.handleInputChange} />
                  <ClinGen type={this.state.type_1} page={this.state.page} var={1} region={this.state.region_1}
                    chrom={this.state.chrom} customStart={this.state.customStart_1} customEnd={this.state.customEnd_1}
                    clingen_id={this.state.clingen_id_1} onChange={this.handleInputChangeClinGen}
                    genome={this.state.genome} gene_uid={this.state.gene_uid} />
                  <STR type={this.state.type_1} page={this.state.page} var={1} str_id={this.state.str_id_1} onChange={this.handleInputChange}
                    genome={this.state.genome} gene_uid={this.state.gene_uid} region={this.state.region_1}
                    chrom={this.state.chrom} customStart={this.state.customStart_1} customEnd={this.state.customEnd_1} />
                  <Clinvar type={this.state.type_1} page={this.state.page} var={1}
                    clinvar_id={this.state.clinvar_id_1} onChange={this.handleInputChange}
                    genome={this.state.genome} gene_uid={this.state.gene_uid} region={this.state.region_1}
                    chrom={this.state.chrom} customStart={this.state.customStart_1} customEnd={this.state.customEnd_1} />
                  {/**************** var2 ****************/}
                  <VariantInfo var2={(this.state.zygosity_1 !== "heterozygous") ? false : this.state.var2} page={this.state.page} var={2} sex={this.state.sex} chrom={this.state.chrom}
                    type={this.state.type_2} region={this.state.region_2}
                    customStart={this.state.customStart_2} customEnd={this.state.customEnd_2} zygosity={this.state.zygosity_2}
                    onChange={this.handleInputChange} onChangeClear={this.handelInputChangeClearFields} />
                  <CNV type={this.state.type_2} page={this.state.page} var={2}
                    start={this.state.start_2} end={this.state.end_2} length={this.state.length_2} onChange={this.handleInputChange} />
                  <Indel type={this.state.type_2} page={this.state.page} var={2}
                    start={this.state.start_2} length={this.state.length_2} onChange={this.handleInputChange} />
                  <MEI type={this.state.type_2} page={this.state.page} var={2}
                    element={this.state.element_2} onChange={this.handleInputChange} />
                  <SNV type={this.state.type_2} page={this.state.page} var={2}
                    start={this.state.start_2} snv_type={this.state.snv_type_2} onChange={this.handleInputChange} />
                  <ClinGen type={this.state.type_2} page={this.state.page} var={2} region={this.state.region_2}
                    chrom={this.state.chrom} customStart={this.state.customStart_2} customEnd={this.state.customEnd_2}
                    clingen_id={this.state.clingen_id_2} onChange={this.handleInputChangeClinGen}
                    genome={this.state.genome} gene_uid={this.state.gene_uid} />
                  <STR type={this.state.type_2} page={this.state.page} var={2} str_id={this.state.str_id_2} onChange={this.handleInputChange}
                    genome={this.state.genome} gene_uid={this.state.gene_uid} region={this.state.region_2}
                    chrom={this.state.chrom} customStart={this.state.customStart_2} customEnd={this.state.customEnd_2} />
                  <Clinvar type={this.state.type_2} page={this.state.page} var={2}
                    clinvar_id={this.state.clinvar_id_2} onChange={this.handleInputChange}
                    genome={this.state.genome} gene_uid={this.state.gene_uid} region={this.state.region_2}
                    chrom={this.state.chrom} customStart={this.state.customStart_2} customEnd={this.state.customEnd_1} />
                  {/**************** family ****************/}
                  <FInfo page={this.state.page} family={this.state.family} vars={this.state.vars} sex={this.state.sex}
                    onAdd={this.addFamily} onRemove={this.removeFamily} onChange={this.handleFamilyCheckChange} downloadFile={this.downloadFile} />

                </div>
              </div>
              {/**************** buttons ****************/}
              <div className="row">
                <div className="col-sm">
                  <NavButtons page={this.state.page} next={this.next} back={this.back} get_vars={this.get_vars} />
                </div>
              </div>
              <div className="row errors">
                <div className="col-sm">
                  <Errors errors={this.state.errors} />
                </div>
              </div>
            </form>
          </div >
        </div >
      </div>
    );
  }
}

export default MasterForm;
