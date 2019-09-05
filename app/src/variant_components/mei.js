import React from 'react';

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
export default Mei;