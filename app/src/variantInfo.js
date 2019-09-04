import React from 'react';
import NavButtons from './navButtons'

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
function Str(props) {
  if (props.type !== "str") {
    return null;
  }
  return (
    <React.Fragment>

    </React.Fragment>
  )
}
function CNV(props) {
  //start, end, copy_change
  if (props.type !== "cnv") {
    return null;
  }
  return (
    <React.Fragment>
      <label>Start position</label>
      <input type="number" className="form-control" min="1"
        name={"var" + props.var + "_start"}
        value={props.start}
        onChange={props.handleInputChange} />
      <label>End position</label>
      <input type="number" className="form-control" min="2"
        name={"var" + props.var + "_end"}
        value={props.end}
        onChange={props.handleInputChange} />
      <label>Copy Change</label>
      <input type="number" className="form-control" min="-1"
        name={"var" + props.var + "_copy_change"}
        value={props.copy_change}
        onChange={props.handleInputChange} />
      <small className="form-text text-muted">
        input an integer copy number change where -1 is a deletion and positive numbers greater than 1 indicated multiplications.
      </small>
    </React.Fragment>
  )
}
function Indel(props) {
  // start, length
  if (props.type !== "indel") {
    return null;
  }
  return (<React.Fragment>
    <label>Start position</label>
    <input type="text" className="form-control"
      name={"var" + props.var + "_start"}
      value={props.start}
      onChange={props.handleInputChange} />
    <small className="form-text text-muted">input the start position of your variant or "ANY" if you want it anywhere within the region.
    </small>
    <label>Length</label>
    <input type="number" className="form-control" min="-200" max="200"
      name={"var" + props.var + "_length"}
      value={props.length}
      onChange={props.handleInputChange} />
    <small className="form-text text-muted">
      input an integer length for your indel where negative numbers are deletions and positive numbers insertions, indels must be between -200 and 200.
      </small>
  </React.Fragment>)
}
function Mei(props) {
  // element, start
  if (props.type !== "mei") {
    return null;
  }
  return (<React.Fragment>
    <label>Start position</label>
    <input type="text" className="form-control"
      name={"var" + props.var + "_start"}
      value={props.start}
      onChange={props.handleInputChange} />
    <small className="form-text text-muted">input the start position of your variant or "ANY" if you want it anywhere within the region.
    </small>
    <div className="form-group">
      <label>Element</label>
      <select className="form-control"
        name={"var" + props.var + "_element"}
        value={props.element}
        onChange={props.handleInputChange}>
        <option value="">Select</option>
        <option value="alu">ALU</option>
        <option value="line">LINE</option>
        <option value="sva">SVA</option>
      </select>
    </div>
  </React.Fragment>)
}
function Snv(props) {
  //start, type
  if (props.type !== "snv") {
    return null;
  }
  return (<React.Fragment>
    <label>Start position</label>
    <input type="text" className="form-control"
      name={"var" + props.var + "_start"}
      value={props.start}
      onChange={props.handleInputChange} />
    <small className="form-text text-muted">input the start position of your variant or "ANY" if you want it anywhere within the region.
    </small>
    <div className="form-group">
      <label>SNV impact</label>
      <select className="form-control"
        name={"var" + props.var + "_snv_type"}
        value={props.snv_type}
        onChange={props.handleInputChange}>
        <option value="">Select</option>
        <option value="stoploss">Stop Loss</option>
        <option value="missense">Missense</option>
        <option value="nonsense">Nonsense</option>
        <option value="synonymous">Synonymous</option>
        <option value="A">A</option>
        <option value="T">T</option>
        <option value="G">G</option>
        <option value="C">C</option>
        <option value="ANY">ANY</option>
      </select>
    </div>
  </React.Fragment>)
}


function Zygosity(props) {
  // XY sex and X or Y gene 
  if (props.sex === "XY" && (props.chrom === "chrY" || props.chrom === "chrX")) {
    return <option value="hemizygous">Hemizygous</option>
  } else if (props.var === 2) { // variant 2
    return <option value="heterozygous">Heterozygous</option>
  } else { // variant 1
    return (
      <React.Fragment>
        <option value="homozygous">Homozygous</option>
        <option value="heterozygous">Heterozygous</option>
      </React.Fragment>)
  }
}


class VariantInfo extends React.Component {
  //page, chrom, sex, gene_uid, handleInputChange, next, back
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
            <option value="coding">Coding</option>
            <option value="genic">Genic</option>
            <option value="utr">UTR</option>
            <option value="intronic">Intronic</option>
            <option value="custom">Custom</option>
          </select>
        </div>
        {/* customRegion */}
        {this.props.region === "custom" ?
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
          handleInputChange={this.props.handleInputChange} />
        <Clingen
          type={this.props.type}
          handleInputChange={this.props.handleInputChange} />
        <Str
          type={this.props.type}
          handleInputChange={this.props.handleInputChange} />
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
