import React from 'react';

function Indel(props) {
    // start, length
    if (props.type !== "indel" || props.page -1 !== props.var) {
      return null;
    }
    return (<React.Fragment>
      <label>Start position</label>
      <input type="text" className="form-control"
        name={"start_" + props.var }
        value={props.start}
        onChange={props.onChange} />
      <small className="form-text text-muted">input the start position of your variant or "ANY" if you want it anywhere within the region.
      </small>
      <label>Length</label>
      <input type="number" className="form-control" min="-200" max="200"
        name={"length_" + props.var }
        value={props.length}
        onChange={props.onChange} />
      <small className="form-text text-muted">
        input an integer length for your indel where negative numbers are deletions and positive numbers insertions, indels must be between -200 and 200.
        </small>
    </React.Fragment>)
  }

export default Indel;