import React from 'react';
import SelectComp from './selectComp'
import Zygosity from './zygosity'

function VariantInfo(props) {
  if (props.page !== 2 && props.page !== 3) {
    return null;
  } else if (props.page === 2 && props.var !== 1) {
    return null;
  } else if (props.page === 3 && props.var !== 2) {
    return null;
  }
  const second_variant = <SelectComp
    title={"Do you want a second variant?"}
    name={"var2"}
    value={props.var2}
    onChange={props.onChange}
    options={[{ value: "false", text: "No" },
    { value: "true", text: "Yes" }]} />
  const region = <SelectComp
    title={"Region"}
    name={"region_" + props.var}
    value={props.region}
    onChange={props.onChangeClear}
    options={[{ value: "", text: "Select" },
    { value: "CODING", text: "Coding" },
    { value: "GENIC", text: "Genic" },
    { value: "UTR", text: "Untranslated Region" },
    { value: "INTRONIC", text: "Intronic" },
    { value: "CUSTOM", text: "Custom" }]} />
  const custom_pos = (<div>
    <label>Custom region</label>
    <div className="input-group" >
      <input type="text" className="form-control" value={props.chrom} disabled />
      <input type="text" className="form-control"
        name={"customStart_" + props.var} value={props.customStart} onChange={props.onChange} />
      <input type="text" className="form-control"
        name={"customEnd_" + props.var} value={props.customEnd} onChange={props.onChange} />
    </div>
    <small className="form-text text-muted">chrZ:start-end</small>
  </div>)
  const type = <SelectComp
    title={"Type"}
    name={"type_" + props.var}
    value={props.type}
    onChange={props.onChangeClear}
    options={[{ value: "", text: "Select" },
    { value: "clinvar", text: "ClinVar Variant" },
    { value: "clingen", text: "Clingen Copy Number Variant" },
    { value: "cnv", text: "Copy Number Variant" },
    { value: "indel", text: "Indel" },
    { value: "mei", text: "Mobile Element Insertion" },
    { value: "snv", text: "Single Nucleotide Variant" },
    { value: "str", text: "Short Tandem Repeat" }]} />
  const zygosity = (<div className="form-group">
    <label>Zygosity</label>
    <select className="form-control"
      name={"zygosity_" + props.var}
      value={props.zygosity}
      onChange={props.onChange}>
      <option value="">Select</option>
      <Zygosity chrom={props.chrom} sex={props.sex} var={props.var} />
    </select>
  </div>)
  if (props.page === 2) {
    return (
      <React.Fragment>
        {region}
        {props.region === "CUSTOM" ? custom_pos : null}
        {type}
        {zygosity}
      </React.Fragment >
    )
  } else if (!props.var2) {
    return <React.Fragment>No second variant can be made with homozygous/hemizygous first variant.</React.Fragment>
  }
  else {
    return (
      <React.Fragment>
        {second_variant}
        {props.var2==="true"? region: null}
        {(props.region === "CUSTOM" &&  props.var2==="true")? custom_pos : null}
        {props.var2==="true"? type: null}
        {props.var2==="true"? zygosity: null}
      </React.Fragment >
    )
  }
}


export default VariantInfo;