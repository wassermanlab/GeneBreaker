import React from 'react';

function Mei(props) {
    // element, start
    if (props.type !== "mei" || props.page -1 !== props.var) {
      return null;
    }
    return (<React.Fragment>
      <label>Start position</label>
      <input type="text" className="form-control"
        name={"start_"  + props.var}
        value={props.start}
        onChange={props.onChange} />
      <small className="form-text text-muted">input the start position of your variant or "ANY" if you want it anywhere within the region.
      </small>
      <div className="form-group">
        <label>Element</label>
        <select className="form-control"
          name={"element_"  + props.var}
          value={props.element}
          onChange={props.onChange}>
          <option value="">Select</option>
          <option value="ALU">ALU</option>
          <option value="LINE">LINE</option>
          <option value="SVA">SVA</option>
        </select>
      </div>
    </React.Fragment>)
  }
export default Mei;