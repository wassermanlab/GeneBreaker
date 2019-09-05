import React from 'react';

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

export default Snv;