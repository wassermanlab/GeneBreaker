import React from 'react';
import SelectComp from '../selectComp'

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
      <SelectComp
          title={"SNV Impact"}
          name={"var" + props.var + "_snv_type"}
          value={props.snv_type}
          onChange={props.handleInputChange}
          options={[{value: "", text: "Select"},
          {value: "stoploss", text: "Stoploss"},
          {value: "missense", text: "Missense"},
          {value: "nonsense", text: "Nonsense"},
          {value: "synonymous", text: "Synonymous"},
          {value: "A", text: "A"},
          {value: "T", text: "T"},
          {value: "G", text: "G"},
          {value: "C", text: "C"},
          {value: "any", text: "Any"}]}
          />
    </React.Fragment>)
  }

export default Snv;