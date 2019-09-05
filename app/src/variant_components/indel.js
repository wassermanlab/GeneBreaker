import React from 'react';

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

export default Indel;