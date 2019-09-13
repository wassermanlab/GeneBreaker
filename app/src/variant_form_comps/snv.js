import React from 'react';
import SelectComp from './selectComp'

function Snv(props) {
    //start, type
    if (props.type !== "snv" || props.page -1 !== props.var) {
      return null;
    }
    return (<React.Fragment>
      <label>Start position</label>
      <input type="text" className="form-control"
        name={"var" + props.var + "_start"}
        value={props.start}
        onChange={props.onChange} />
      <small className="form-text text-muted">input the start position of your variant or "ANY" if you want it anywhere within the region.
      </small>
      <SelectComp
          title={"SNV Impact"}
          name={"var" + props.var + "_snv_type"}
          value={props.snv_type}
          onChange={props.onChange}
          options={[{value: "", text: "Select"},
          {value: "STOPLOSS", text: "Stoploss"},
          {value: "MISSENSE", text: "Missense"},
          {value: "NONSENSE", text: "Nonsense"},
          {value: "SYNONYMOUS", text: "Synonymous"},
          {value: "A", text: "A"},
          {value: "T", text: "T"},
          {value: "G", text: "G"},
          {value: "C", text: "C"},
          {value: "ANY", text: "Any"}]}
          />
    </React.Fragment>)
  }

export default Snv;