import React from 'react';
import NavButtons from './navButtons'
import CNV from './variant_components/cnv'
import Clinvar from './variant_components/clinvar'
import Clingen from './variant_components/clingen'
import Indel from './variant_components/indel'
import Mei from './variant_components/mei'
import Snv from './variant_components/snv'
import Str from './variant_components/str'
import Zygosity from './variant_components/zygosity'

// todo: move to separate component 
function VariantErrors(props) {
  const error = [];
  for (var key in props.errors) {
    error.push(<div key={key} className="alert alert-danger" role="alert">{props.errors[key]}</div>)
  }
  return error;
}

class VariantInfo extends React.Component {
  //page, chrom, sex, gene_uid, handleInputChange, next, back, genome
  //var, var2, type, region, zygosity, clinvar_id, start, end, length, element, snv_type, motif

  constructor(props) {
    super(props);
    this.state = {
      errors: {},
    };

    this.getOptions = this.getOptions.bind(this);
    this.customNext = this.customNext.bind(this);

  }

  // error check for the whole page!
  customNext() {
    errors = {}
    for (var val of ["region", "zygosity", "type"]) { // checks that region zygosity and type are filled
      if (this.props[val] === "") {
        errors["no_" + val] = "no " + val + " selected."
      }
    }
    if (this.props.region === "CUSTOM") { // check custom region
      if (parseInt(this.props.customStart) >= parseInt(this.props.customEnd)) {
        errors["custom_region"] = "start of custom region must be less than end."
      }
    }
    switch (this.props.type) { // check types specific stuff 
      case cnv:
        if (parseInt(this.props.start) >= this.props.end) {
          errors["cnv_start_end"] = "start of cnv must be less than end."
        }
        if (parseInt(this.props.length) < -1 || parseInt(this.props.length) === 0) {
          errors["cnv_copy_change"] = "copy change cannot be less than -1 or equal 0."
        }
        break;
      case clingen:
        if (parseInt(this.props.start) >= this.props.end) {
          errors["cnv_start_end"] = "start of cnv must be less than end."
        }
        if (parseInt(this.props.length) < -1 || parseInt(this.props.length) === 0) {
          errors["cnv_copy_change"] = "copy change cannot be less than -1 or equal 0."
        }
        break;
      case clinvar:
        if (this.props.clinvar_id === "") {
          errors["no_clinvar_id"] = "clinvar variant must be selected."
        }
        break;
      case indel:
        if (parseInt(this.props.start) < 1) {
          errors["indel_start"] = "must pick valid start of indel."
        }
        if (parseInt(this.props.length) > 200 || parseInt(this.props.length) < -200 || parseInt(this.props.length) == 0) {
          errors["indel_length"] = "indel length  must be between -200 and 200 and cannot equal to 0."
        }
        break;
      case mei:
        if (parseInt(this.props.start) < 1) {
          errors["mei_start"] = "must pick valid start for mei."
        }
        break;
      case snv:
        if (parseInt(this.props.start) < 1) {
          errors["snv_start"] = "must pick valid start for snv."
        }
        if (this.props.snv_type === "") {
          errors["snv_type"] = "must select valid snv type."
        }
        break;
      case str:
        if (parseInt(this.props.start) < 1) {
          errors["str_start"] = "must pick valid start for str."
        }
        if (this.props.length === 0) {
          errors["str_length"] = "length of str cannot equal 0"
        }
        break;
      default:
        break;
    }
    this.setState({ errors: errors })
    if (errors == {}) {
      this.props.next()
    }
  }

  render() {
    // render only if page is 1 or 2 
    if (this.props.page !== 2 && this.props.page !== 3) {
      return null;
    }

    return (

      // {/* general: genome, sex, gene_uid, chr, transcript */}
      <React.Fragment>
        {/* region */}
        <div className="form-group">
          <label>Region</label>
          <select className="form-control"
            name={"var" + this.props.var + "_region"}
            value={this.props.region}
            onChange={this.props.handleInputChange}>
            <option value="">Select</option>
            <option value="CODING">Coding</option>
            <option value="GENIC">Genic</option>
            <option value="UTR">UTR</option>
            <option value="INTRONIC">Intronic</option>
            <option value="CUSTOM">Custom</option>
          </select>
        </div>
        {/* customRegion */}
        {this.props.region === "CUSTOM" ?
          <div>
            <label>Custom region</label>
            <div className="input-group" >
              <input type="text" className="form-control" value={this.props.chrom} disabled />
              <input type="text" className="form-control"
                name={"var" + this.props.var + "_customStart"} value={this.props.customStart} onChange={this.props.handleInputChange} />
              <input type="text" className="form-control"
                name={"var" + this.props.var + "_customEnd"} value={this.props.customEnd} onChange={this.props.handleInputChange} />
            </div>
            <small className="form-text text-muted">
              chrZ:start-end
          </small>
          </div>
          : null}
        {/* type */}
        <div className="form-group">
          <label>Type</label>
          <select className="form-control"
            name={"var" + this.props.var + "_type"}
            value={this.props.type}
            onChange={this.props.handleInputChange}>
            <option value="">Select</option>
            <option value="clinvar">ClinVar Variant</option>
            <option value="clingen">Clingen Copy Number Variant</option>
            <option value="cnv">Copy Number Variant</option>
            <option value="indel">Indel</option>
            <option value="mei">Mobile Element Insertion</option>
            <option value="snv">Single Nucleotide Variant</option>
            <option value="str">Short Tandem Repeat</option>
          </select>
        </div>
        {/* impact */}
        <Clinvar
          type={this.props.type}
          handleInputChange={this.props.handleInputChange}
          clinvar_id={this.props.clinvar_id}
          var={this.props.var}
          gene_uid={this.props.gene_uid}
          genome={this.props.genome}
          region={this.props.region === "CUSTOM" ?
            (this.props.chrom + ":" + this.props.customStart + "-" + this.props.customEnd) :
            this.props.region} />
        <Clingen
          type={this.props.type}
          handleInputChange={this.props.handleInputChange}
          clingen_uid={this.props.clingen_uid}
          var={this.props.var}
          gene_uid={this.props.gene_uid}
          genome={this.props.genome}
          region={this.props.region} />
        <Str
          type={this.props.type}
          handleInputChange={this.props.handleInputChange}
          str_id={this.props.str_id}
          var={this.props.var}
          gene_uid={this.props.gene_uid}
          genome={this.props.genome}
          region={this.props.region === "CUSTOM" ?
            (this.props.chrom + ":" + this.props.customStart + "-" + this.props.customEnd) :
            this.props.region}
          length={this.props.length} />
        <CNV
          type={this.props.type}
          handleInputChange={this.props.handleInputChange}
          start={this.props.start}
          end={this.props.end}
          copy_change={this.props.length}
          var={this.props.var} />
        <Indel
          type={this.props.type}
          handleInputChange={this.props.handleInputChange}
          start={this.props.start}
          length={this.props.length}
          var={this.props.var}
        />
        <Mei
          type={this.props.type}
          handleInputChange={this.props.handleInputChange}
          start={this.props.start}
          element={this.props.element}
          var={this.props.var} />
        <Snv
          type={this.props.type}
          handleInputChange={this.props.handleInputChange}
          start={this.props.start}
          end={this.props.end}
          var={this.props.var} />
        {/* zygosity */}
        <div className="form-group">
          <label>Zygosity</label>
          <select className="form-control"
            name={"var" + this.props.var + "_zygosity"}
            value={this.props.zygosity}
            onChange={this.props.handleInputChange}>
            <option value="">Select</option>
            <Zygosity chrom={this.props.chrom} sex={this.props.sex} var={this.props.var} />
          </select>
        </div>
        <VariantErrors
          errors={this.state.errors}
        />
        <NavButtons
          next={this.props.back}
          back={this.props.back}
          page={this.props.page} />
      </React.Fragment>
    );
  }


}

export default VariantInfo;
