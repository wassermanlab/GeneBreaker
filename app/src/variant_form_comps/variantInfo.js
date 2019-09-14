import React from 'react';
import SelectComp from './selectComp'
import Zygosity from './zygosity'

function VariantInfo(props) {
  if (props.page !== 2 && props.page !== 3) {
    return null;
  }
  return (
    <React.Fragment>
      {/* Region */}
      <SelectComp
        title={"Region"}
        name={"region_" + props.var}
        value={props.region}
        onChange={props.onChangeClear}
        options={[{ value: "", text: "Select" },
        { value: "CODING", text: "Coding" },
        { value: "GENIC", text: "Genic" },
        { value: "UTR", text: "Untranslated Region" },
        { value: "INTRONIC", text: "Intronic" },
        { value: "CUSTOM", text: "Custom" }]}
      />
      {/* customRegion */}
      {props.region === "CUSTOM" ?
        <div>
          <label>Custom region</label>
          <div className="input-group" >
            <input type="text" className="form-control" value={props.chrom} disabled />
            <input type="text" className="form-control"
              name={"customStart_" + props.var} value={props.customStart} onChange={props.onChange} />
            <input type="text" className="form-control"
              name={"customEnd_" + props.var} value={props.customEnd} onChange={props.onChange} />
          </div>
          <small className="form-text text-muted">
            chrZ:start-end
          </small>
        </div>
        : null}
      {/* type */}
      <SelectComp
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
        { value: "str", text: "Short Tandem Repeat" }]}
      />
      {/* zygosity */}
      <div className="form-group">
        <label>Zygosity</label>
        <select className="form-control"
          name={"zygosity_" + props.var}
          value={props.zygosity}
          onChange={props.onChange}>
          <option value="">Select</option>
          <Zygosity chrom={props.chrom} sex={props.sex} var={props.var} />
        </select>
      </div>

    </React.Fragment >
  )
}


export default VariantInfo;