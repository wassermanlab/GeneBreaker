import React, { useState, useEffect } from 'react';
import NavButtons from './navButtons'
import CNV from './variant_components/cnv'
import Indel from './variant_components/indel'
import Mei from './variant_components/mei'
import Snv from './variant_components/snv'
import Str from './variant_components/str'
import Zygosity from './variant_components/zygosity'


// todo: move to separate component 
function Clinvar(props) {
  if (props.type !== "clinvar") {
    return null;
  }
  return (
    <React.Fragment>
    </React.Fragment>
  )
}
function Clingen(props) {
  if (props.type !== "clingen") {
    return null;
  }
  return (
    <React.Fragment>
    </React.Fragment>
  )
}

class VariantInfo extends React.Component {
  //page, chrom, sex, gene_uid, handleInputChange, next, back, genome
  //var, var2, type, region, zygosity, clinvar_id, start, end, copy_change, length, element, snv_type, motif

  constructor(props) {
    super(props);

    this.getOptions = this.getOptions.bind(this);
    this.customNext = this.customNext.bind(this);

  }

  getOptions(event) {
    return null;
  }

  // error check for the whole page!
  customNext() {
    return null;
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
          region={this.props.region} />
        <Clingen
          type={this.props.type}
          handleInputChange={this.props.handleInputChange}
          start={this.props.start}
          end={this.props.end}
          copy_change={this.props.copy_change}
          var={this.props.var}
          gene_uid={this.props.var}
          genome={this.props.genome}
          region={this.props.region} />
        <Str
          type={this.props.type}
          handleInputChange={this.props.handleInputChange}
          str_id={this.props.str_id}
          var={this.props.var}
          gene_uid={this.props.gene_uid}
          genome={this.props.genome}
          region={this.props.region==="CUSTOM"? 
          (this.props.chrom+":"+this.props.customStart+"-"+this.props.customEnd): 
          this.props.region}
          length={this.props.length} />
        <CNV
          type={this.props.type}
          handleInputChange={this.props.handleInputChange}
          start={this.props.start}
          end={this.props.end}
          copy_change={this.props.copy_change}
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
        <NavButtons
          next={this.props.back}
          back={this.props.back}
          page={this.props.page} />
      </React.Fragment>
    );
  }


}

export default VariantInfo;
