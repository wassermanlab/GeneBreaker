import React from 'react';

function CNV(props) {
    //start, end, copy_change
    if (props.type !== "cnv" || props.page -1 !== props.var) {
      return null;
    }
    return (
      <React.Fragment>
        <label>Start position</label>
        <input type="number" className="form-control" min="1"
          name={"start_" + props.var}
          value={props.start}
          onChange={props.onChange} />
        <label>End position</label>
        <input type="number" className="form-control" min="2"
          name={"end_" + props.var}
          value={props.end}
          onChange={props.onChange} />
        <label>Copy Change</label>
        <input type="number" className="form-control" min="-1"
          name={"length_" + props.var}
          value={props.copy_change}
          onChange={props.onChange} />
        <small className="form-text text-muted">
          input an integer copy number change where -1 is a deletion and positive numbers are additions.
        </small>
      </React.Fragment>
    )
  }


export default CNV;